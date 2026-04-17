
#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "reaction.hpp"
#include "syncfile.hpp"
#include <unordered_set>
#include "TLorentzVector.h"

// Helper function to check if a string contains a substring (case-sensitive)
bool contains(const std::string& str, const std::string& substr) {
  return str.find(substr) != std::string::npos;
}

template <class CutType>
size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<SyncFile>& _sync, int thread_id, const std::string& output_filename) {
  // Get the number of events in this thread
  size_t num_of_events = (int)_chain->GetEntries();

  // Determine the beam energy from the environment variable
  float beam_energy = 10.6;
  if (getenv("BEAM_E") != NULL) beam_energy = atof(getenv("BEAM_E"));

  // Determine the data processing type (gen, rec, or exp) based on output filename
  bool is_gen_data = contains(output_filename, "gen");
  bool is_rec_data = contains(output_filename, "rec");
  bool is_exp_data = contains(output_filename, "exp");

  // Ensure only one data type is selected
  if ((is_gen_data + is_rec_data + is_exp_data) > 1) {
    throw std::invalid_argument("Output filename must specify exactly one data type: gen, rec, or exp.");
  }

  // Determine the topology based on output filename
  bool is_topology_excl = contains(output_filename, "excl");
  bool is_topology_mProt = contains(output_filename, "mProt");
  bool is_topology_mPip = contains(output_filename, "mPip");
  bool is_topology_mPim = contains(output_filename, "mPim");

  // Ensure only one topology is selected
  if ((is_topology_excl + is_topology_mProt + is_topology_mPip + is_topology_mPim) > 1) {
    throw std::invalid_argument("Output filename must specify exactly one topology: excl, mProt, mPip, or mPim.");
  }

  // Print some information for each thread
  std::cout << "=============== " << MAGENTA << "Thread " << thread_id << DEF << " =============== " << BLUE
            << num_of_events << " Events " << DEF << "===============\n";

  // Make a data object which all the branches can be accessed from
  // auto data = is_gen_data || is_rec_data ? std::make_shared<Branches12>(_chain, true) : std::make_shared<Branches12>(_chain);
  auto data = (is_gen_data || is_rec_data) 
               ? std::make_shared<Branches12>(_chain, true)  // For gen and rec
               : std::make_shared<Branches12>(_chain);       // For exp


  // Total number of events "Processed"
  size_t total = 0;
  size_t total_twopion_events = 0;

  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    // Get current event
    _chain->GetEntry(current_event);
    
    // If we are the 0th thread print the progress of the thread every 1000 events
    if (thread_id == 0 && current_event % 1000 == 0)
      std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;
      
      // ----- Process Generated Data -----
      if (is_gen_data) {
        // ----- Generated reaction class -----
        auto mc_event = std::make_shared<MCReaction>(data, beam_energy, "gen");
        for (int part = 1; part < data->mc_npart(); part++) {
          if (data->mc_pid(part) == PIP) {
            mc_event->SetMCPip(part);
          } else if (data->mc_pid(part) == PROTON) {
            mc_event->SetMCProton(part);
          } else if (data->mc_pid(part) == PIM) {
            mc_event->SetMCPim(part);
          }
        }

          // ----- Generated data output -----
        csv_data output;
        output.event = current_event;
        output.w_mc = mc_event->W_mc();
        output.q2_mc = mc_event->Q2_mc();
        output.weight_gen = mc_event->weight();

        _sync->write(output);
        total++;  // Increment for all events when processing gen data
      }

      // ----- Process Reconstructed Data -----
      else if (is_rec_data || is_exp_data) {
        int statusPim = -9999;
        int statusPip = -9999;
        int statusProt = -9999;
        float vertex_hadron[3][3];

        // Make cuts
        if (is_rec_data && data->mc_npart() < 1) continue;
        auto dt = std::make_shared<Delta_T>(data);
        auto cuts = std::make_shared<Pass2_Cuts>(data);
        if (!cuts->ElectronCuts()) continue;
        
        total++;  // Increment only if the event is processed with rec cuts

        // ----- Reconstructed reaction class -----
        auto event = std::make_shared<Reaction>(data, beam_energy, is_rec_data ? "rec" : "exp");
        for (int part = 1; part < data->gpart(); part++) {
          dt->dt_calc(part);
          if (cuts->IsProton(part)) {
            event->SetProton(part);
            statusProt = abs(data->status(part));
          } else if (cuts->IsPip(part)) {
            event->SetPip(part);
            statusPip = abs(data->status(part));
          } else if (cuts->IsPim(part)) {
            event->SetPim(part);
            statusPim = abs(data->status(part));
          } else {
            event->SetOther(part);
          }
        }
        
        double q2_min_analysis = -1.0, q2_max_analysis = 30.0;
        double w_min_analysis = 1.0, w_max_analysis = 2.5;

        // Dynamically set W and Q2 limits based on BEAM_E
        if (getenv("BEAM_E") != NULL) {
          double beam_energy = atof(getenv("BEAM_E"));
          if (beam_energy < 3) {
            q2_max_analysis = 1.0;
            w_max_analysis = 3.5;
            w_min_analysis = 0.9;
          } else if (beam_energy < 11) {
            q2_min_analysis = 0.0;
            q2_max_analysis = 12.0;
            w_min_analysis = 1.0;
            w_max_analysis = 2.5;
          } else if (beam_energy < 24) {
            q2_min_analysis = 2.0;
            q2_max_analysis = 30.0;
            w_min_analysis = 1.0;
            w_max_analysis = 2.5;
          }
        }

        // Update the condition to use the dynamically set limits
        if ((is_topology_excl && event->TwoPion_exclusive()) ||
            (is_topology_mProt && event->TwoPion_missingProt()) ||
            (is_topology_mPip && event->TwoPion_missingPip()) ||
            (is_topology_mPim && event->TwoPion_missingPim())
            ) {
          if (event->W() > w_min_analysis && event->W() < w_max_analysis && 
              event->Q2() > q2_min_analysis && event->Q2() < q2_max_analysis && 
              event->weight() > 0.0) {
            total_twopion_events++;
          // }
        // }
      // }

        // if ((is_topology_excl && event->TwoPion_exclusive()) ||
        //     (is_topology_mProt && event->TwoPion_missingProt()) ||
        //     (is_topology_mPip && event->TwoPion_missingPip()) ||
        //     (is_topology_mPim && event->TwoPion_missingPim())) {
        //   if (event->W() > 1.0 && event->W() < 2.5 && event->Q2() > 2.0 && event->Q2() < 30.0 && event->weight() > 0.0) {
            
        //     total_twopion_events++;

          

            // ----- Reconstructed data output -----
            csv_data output;
            output.event = current_event;
            output.w = event->W();
            output.q2 = event->Q2();
            output.weight_rec = event->weight() * 1e4;

            output.mm2_mPim = event->MM2_mPim();
            output.mm2_mPip = event->MM2_mPip();
            output.mm2_mProt = event->MM2_mProt();
            output.mm2_exclusive = event->MM2_exclusive();

            output.pim_mom_miss = event->pim_momentum();
            output.pim_mom_meas = event->pim_momentum_measured();
            output.pip_mom_miss = event->pip_momentum();
            output.pip_mom_meas = event->pip_momentum_measured();
            output.prot_mom_miss = event->prot_momentum();
            output.prot_mom_meas = event->prot_momentum_measured();
            output.excl_mom = event->Mom_excl();

            output.pim_theta_miss = event->pim_theta_lab();
            output.pim_theta_meas = event->pim_theta_lab_measured();
            output.pip_theta_miss = event->pip_theta_lab();
            output.pip_theta_meas = event->pip_theta_lab_measured();
            output.prot_theta_miss = event->prot_theta_lab();
            output.prot_theta_meas = event->prot_theta_lab_measured();

            output.prot_theta_angle_btwn_P = event->prot_theta_angle_btwn_P();
            output.pip_theta_angle_btwn_P = event->pip_theta_angle_btwn_P();
            output.pim_theta_angle_btwn_P = event->pim_theta_angle_btwn_P();

            _sync->write(output);
          }
        }
      }
    }

    std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
    std::cout << " total no of events = " << total << std::endl;
    std::cout << " total no of twopion events = " << total_twopion_events << std::endl;

    return num_of_events;
  }

#endif
