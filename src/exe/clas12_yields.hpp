
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
#include <algorithm>

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

  std::cout << "Branches12 constructed for thread "
            << thread_id << std::endl;
  // Total number of events "Processed"
  // size_t total = 0;
  // std::atomic<size_t> total{0};
  // size_t total_twopion_events = 0;

  static std::atomic<size_t> total_events{0};
  static std::atomic<size_t> total_twopion_events{0};

  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
  // size_t current_event;

    // Get current event
    _chain->GetEntry(current_event);

    if (current_event < 10) {
        std::cout
            << "event=" << current_event
            << " gpart=" << data->gpart()
            << " mc_npart=" << data->mc_npart()
            << " run=" << data->run()
            << " event=" << data->event()
            << std::endl;
    }
    
    // If we are the 0th thread print the progress of the thread every 1000 events
    if (thread_id == 0 && current_event % 1000 == 0)
      std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;
      
    // if (current_event % 1000 == 0)
    //   std::cout << "Thread " << thread_id
    //             << " " << (100 * current_event / num_of_events)
    //             << " %\n";

      // ----- Process Generated Data -----
      if (is_gen_data) {
        csv_data output;

        output.pid_prot_mc = 0;
        output.pid_pip_mc = 0;
        output.pid_pim_mc = 0;
        // ----- Generated reaction class -----
        auto mc_event = std::make_shared<MCReaction>(data, beam_energy, "gen");
        for (int part = 1; part < data->mc_npart(); part++) {
          int mc_pid = data->mc_pid(part);

          if (data->mc_pid(part) == PIP) {
            mc_event->SetMCPip(part);
            output.pid_pip_mc = mc_pid;
          } else if (data->mc_pid(part) == PROTON) {
            mc_event->SetMCProton(part);
            output.pid_prot_mc = mc_pid;
          } else if (data->mc_pid(part) == PIM) {
            mc_event->SetMCPim(part);
            output.pid_pim_mc = mc_pid;
          }
        }

          // ----- Generated data output -----
        // csv_data output;
        // output.event = current_event;
        output.mc_run = data->mc_run();
        output.mc_event = data->mc_event();
        output.w_mc = mc_event->W_mc();
        output.q2_mc = mc_event->Q2_mc();
        output.weight_gen = mc_event->weight();

        _sync->write(output);
        // total++;  // Increment for all events when processing gen data
        // total.fetch_add(1);
        total_events.fetch_add(1, std::memory_order_relaxed);
      }

      // ----- Process Reconstructed Data -----
      else if (is_rec_data || is_exp_data) {
        int statusPim = -9999;
        int statusPip = -9999;
        int statusProt = -9999;
        float vertex_hadron[3][3];

        csv_data output;

        output.pid_prot_rec = 0;
        output.pid_pip_rec = 0;
        output.pid_pim_rec = 0;

        output.pid_prot_mc = 0;
        output.pid_pip_mc = 0;
        output.pid_pim_mc = 0;

        // int prot_idx = -1;
        // int pip_idx  = -1;
        // int pim_idx  = -1;

        // Make cuts
        if (is_rec_data && data->mc_npart() < 1) continue;
        auto dt = std::make_shared<Delta_T>(data);
        auto cuts = std::make_shared<Pass2_Cuts>(data);

        if (current_event == 0) {
            std::cout << "ElectronCuts first event = "
                      << cuts->ElectronCuts()
                      << std::endl;
        }
        
        if (!cuts->ElectronCuts()) continue;
        
        // total++;  // Increment only if the event is processed with rec cuts
        // total.fetch_add(1);
        total_events.fetch_add(1, std::memory_order_relaxed);

        // ----- Reconstructed reaction class -----
        auto event = std::make_shared<Reaction>(data, beam_energy, is_rec_data ? "rec" : "exp");

        // std::vector<int> proton_candidates;
        // std::vector<int> pip_candidates;
        // std::vector<int> pim_candidates;

        std::vector<std::pair<int,double>> proton_score;
        std::vector<std::pair<int,double>> pip_score;
        std::vector<std::pair<int,double>> pim_score;

        for (int part = 1; part < data->gpart(); part++) {
            dt->dt_calc(part);

            double p2 =
                data->px(part)*data->px(part) +
                data->py(part)*data->py(part) +
                data->pz(part)*data->pz(part);

            if (cuts->IsProton(part)) {
                proton_score.emplace_back(part, p2);
            }

            else if (cuts->IsPip(part)) {
                pip_score.emplace_back(part, p2);
            }

            else if (cuts->IsPim(part)) {
                pim_score.emplace_back(part, p2);
            }
        }

        int prot_idx = -1;
        int pip_idx  = -1;
        int pim_idx  = -1;

        auto best = [&](auto &vec) {
            if (vec.empty()) return -1;

            return std::max_element(
                vec.begin(),
                vec.end(),
                [](auto &a, auto &b) {
                    return a.second < b.second;
                }
            )->first;
        };

        prot_idx = best(proton_score);
        pip_idx  = best(pip_score);
        pim_idx  = best(pim_score);

        if (prot_idx >= 0) event->SetProton(prot_idx);
        if (pip_idx  >= 0) event->SetPip(pip_idx);
        if (pim_idx  >= 0) event->SetPim(pim_idx);

        auto get_mc_pid = [&](int rec_index) {
            if (rec_index < 0) return -999;
            int mc_idx = data->rectoGen_mcindex(rec_index);
            
            if (mc_idx >= 0 && mc_idx < data->mc_npart()) {
                return data->mc_pid(mc_idx);
            }

            return -999;
        };

        if (prot_idx >= 0) {
            output.pid_prot_rec = data->pid(prot_idx);
            output.pid_prot_mc  = get_mc_pid(prot_idx);
        }

        if (pip_idx >= 0) {
            output.pid_pip_rec = data->pid(pip_idx);
            output.pid_pip_mc  = get_mc_pid(pip_idx);
        }

        if (pim_idx >= 0) {
            output.pid_pim_rec = data->pid(pim_idx);
            output.pid_pim_mc  = get_mc_pid(pim_idx);
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
            // total_twopion_events++;
            total_twopion_events.fetch_add(1, std::memory_order_relaxed);

            // ----- Reconstructed data output -----
            // output.event = current_event;
            // output.run = data->run();
            output.event = data->event();
            output.w = event->W();
            output.q2 = event->Q2();
            output.weight_rec = event->weight(); // * 1e4;

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

    // std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
    // std::cout << " total no of events = " << total << std::endl;
    // std::cout << " total no of twopion events = " << total_twopion_events << std::endl;
    std::cout << "Percent = " << 100.0 * total_events.load() / num_of_events << std::endl;
    std::cout << " total no of processed events = " << total_events.load() << std::endl;
    std::cout << " total no of twopion events = " << total_twopion_events.load() << std::endl;

    return num_of_events;
  }

#endif





// #ifndef MAIN_H_GUARD
// #define MAIN_H_GUARD

// #include <iostream>
// #include <string>
// #include "TFile.h"
// #include "TH1.h"
// #include "branches.hpp"
// #include "colors.hpp"
// #include "cuts.hpp"
// #include "reaction.hpp"
// #include "syncfile.hpp"
// #include <unordered_set>
// #include "TLorentzVector.h"

// // Helper function to check if a string contains a substring (case-sensitive)
// bool contains(const std::string& str, const std::string& substr) {
//   return str.find(substr) != std::string::npos;
// }

// template <class CutType>
// size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<SyncFile>& _sync, int thread_id, const std::string& output_filename) {
//   // Get the number of events in this thread
//   size_t num_of_events = (int)_chain->GetEntries();

//   // Determine the beam energy from the environment variable
//   float beam_energy = 10.6;
//   if (getenv("BEAM_E") != NULL) beam_energy = atof(getenv("BEAM_E"));

//   // Determine the data processing type (gen, rec, or exp) based on output filename
//   bool is_gen_data = contains(output_filename, "gen");
//   bool is_rec_data = contains(output_filename, "rec");
//   bool is_exp_data = contains(output_filename, "exp");

//   // Ensure only one data type is selected
//   if ((is_gen_data + is_rec_data + is_exp_data) > 1) {
//     throw std::invalid_argument("Output filename must specify exactly one data type: gen, rec, or exp.");
//   }

//   // Determine the topology based on output filename
//   bool is_topology_excl = contains(output_filename, "excl");
//   bool is_topology_mProt = contains(output_filename, "mProt");
//   bool is_topology_mPip = contains(output_filename, "mPip");
//   bool is_topology_mPim = contains(output_filename, "mPim");

//   // Ensure only one topology is selected
//   if ((is_topology_excl + is_topology_mProt + is_topology_mPip + is_topology_mPim) > 1) {
//     throw std::invalid_argument("Output filename must specify exactly one topology: excl, mProt, mPip, or mPim.");
//   }

//   // Print some information for each thread
//   std::cout << "=============== " << MAGENTA << "Thread " << thread_id << DEF << " =============== " << BLUE
//             << num_of_events << " Events " << DEF << "===============\n";

//   // std::cout << "PID MODE = " << (csv_data::use_thrown_pid ? "THROWN" : "REC") << std::endl;

//   // std::cout << "Thread " << thread_id << " processing " << num_of_events << " events on " << std::this_thread::get_id() << std::endl;

//   // Make a data object which all the branches can be accessed from
//   // auto data = is_gen_data || is_rec_data ? std::make_shared<Branches12>(_chain, true) : std::make_shared<Branches12>(_chain);
//   auto data = (is_gen_data || is_rec_data) 
//                ? std::make_shared<Branches12>(_chain, true)  // For gen and rec
//                : std::make_shared<Branches12>(_chain);       // For exp

//   data->init();
//   if (is_gen_data || is_rec_data) data->initMC();

//   // Total number of events "Processed"
//   size_t total = 0;
//   size_t total_twopion_events = 0;

//   // For each event
//   for (size_t current_event = 0; current_event < num_of_events; current_event++) {
//     // Get current event
//     // _chain->GetEntry(current_event);
//     data->GetEntry(current_event);

//     // if (current_event == 0) {
//     //     std::cout << "rectoGen_pindex(0)=" << data->rectoGen_pindex(0) << std::endl;
//     //     std::cout << "rectoGen_mcindex(0)=" << data->rectoGen_mcindex(0) << std::endl;
//     // }

//     if (current_event == 0) {
//         std::cout << "px[0] = " << data->px(0) << std::endl;
//     }
    
//     // If we are the 0th thread print the progress of the thread every 1000 events
//     if (thread_id == 0 && current_event % 1000 == 0)
//       std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;
      
//     // if (current_event % 1000 == 0)
//     //   std::cout << "Thread " << thread_id
//     //             << " " << (100 * current_event / num_of_events)
//     //             << " %\n";

//       // csv_data output;

//       // ----- Process Generated Data -----
//       if (is_gen_data) {
//         csv_data output;

//         output.pid_prot_mc = 0;
//         output.pid_pip_mc = 0;
//         output.pid_pim_mc = 0;
//         // ----- Generated reaction class -----
//         auto mc_event = std::make_shared<MCReaction>(data, beam_energy, "gen");
//         for (int part = 1; part < data->mc_npart(); part++) {
//           int mc_pid = data->mc_pid(part);

//           if (current_event == 0 && part == 1) {
//               std::cout << "px = " << data->px(part) << std::endl;
//               std::cout << "branch test done" << std::endl;
//           }

//           if (data->mc_pid(part) == PIP) {
//             mc_event->SetMCPip(part);
//             output.pid_pip_mc = mc_pid;
//           } else if (data->mc_pid(part) == PROTON) {
//             mc_event->SetMCProton(part);
//             output.pid_prot_mc = mc_pid;
//           } else if (data->mc_pid(part) == PIM) {
//             mc_event->SetMCPim(part);
//             output.pid_pim_mc = mc_pid;
//           }
//         }

//           // ----- Generated data output -----
//         // csv_data output;
//         // output.event = current_event;
//         output.run = data->mc_run();
//         output.event = data->mc_event();
//         output.w_mc = mc_event->W_mc();
//         output.q2_mc = mc_event->Q2_mc();
//         output.weight_gen = mc_event->weight();

//         _sync->write(output);
//         total++;  // Increment for all events when processing gen data
//       }

//       // ----- Process Reconstructed Data -----
//       else if (is_rec_data || is_exp_data) {
//         int statusPim = -9999;
//         int statusPip = -9999;
//         int statusProt = -9999;
//         float vertex_hadron[3][3];

//         csv_data output;

//         // Make cuts
//         if (is_rec_data && data->mc_npart() < 1) continue;
//         auto dt = std::make_shared<Delta_T>(data);
//         auto cuts = std::make_shared<Pass2_Cuts>(data);
//         if (!cuts->ElectronCuts()) continue;
        
//         total++;  // Increment only if the event is processed with rec cuts

//         // std::unordered_map<int,int> rec_to_mc_map;

//         // int nmatch = data->rectoGen_n();

//         // for (int i = 0; i < data->rectoGen_n(); i++) {

//         //     int rec_i = data->rectoGen_pindex(i);
//         //     int mc_i  = data->rectoGen_mcindex(i);

//         //     if (part == rec_i) {
//         //         mc_pid = data->mc_pid(mc_i);
//         //     }
//         // }

//         std::unordered_map<int,int> rec_to_mc_map;

//         for (int i = 0; i < data->rectoGen_n(); i++) {
//             int rec_i = data->rectoGen_pindex(i);
//             int mc_i  = data->rectoGen_mcindex(i);

//             rec_to_mc_map[rec_i] = mc_i;
//         }

//         // ----- Reconstructed reaction class -----
//         auto event = std::make_shared<Reaction>(data, beam_energy, is_rec_data ? "rec" : "exp");

//         struct Candidate {
//           int part;
//           int rec_pid;
//           int mc_pid;
//           int pid_type; // NEW: PROTON, PIP, PIM
//           double px, py, pz;
//           int charge;
//           double dp;
//         };

//         std::vector<Candidate> proton_cands;
//         std::vector<Candidate> pip_cands;
//         std::vector<Candidate> pim_cands;

//         for (int part = 1; part < data->gpart(); part++) {
//           dt->dt_calc(part);

//           // int rec_pid = data->pid(part);
//           // int mc_pid  = data->mc_pid(part);

//           int rec_pid = data->pid(part);
//           int mc_pid = -9999;

//           auto it = rec_to_mc_map.find(part);
//           if (it != rec_to_mc_map.end()) {
//               int mc_index = it->second;

//               if (mc_index >= 0 && mc_index < data->mc_npart()) {
//                   mc_pid = data->mc_pid(mc_index);
//               }
//           }

//           Candidate c = {};
//           c.part = part;
//           c.rec_pid = rec_pid;
//           c.mc_pid = mc_pid;
//           c.px = data->px(part);
//           c.py = data->py(part);
//           c.pz = data->pz(part);
//           c.charge = data->charge(part);
//           c.dp = 0.0; // optional later

//           if (cuts->IsProton(part)) {
//               c.pid_type = PROTON;
//               proton_cands.push_back(c);
//               event->SetProton(part);
//           }
//           else if (cuts->IsPip(part)) {
//               c.pid_type = PIP;
//               pip_cands.push_back(c);
//               event->SetPip(part);
//           }
//           else if (cuts->IsPim(part)) {
//               c.pid_type = PIM;
//               pim_cands.push_back(c);
//               event->SetPim(part);
//           }
//         }
        
//         double q2_min_analysis = -1.0, q2_max_analysis = 30.0;
//         double w_min_analysis = 1.0, w_max_analysis = 2.5;

//         // Dynamically set W and Q2 limits based on BEAM_E
//         if (getenv("BEAM_E") != NULL) {
//           double beam_energy = atof(getenv("BEAM_E"));
//           if (beam_energy < 3) {
//             q2_max_analysis = 1.0;
//             w_max_analysis = 3.5;
//             w_min_analysis = 0.9;
//           } else if (beam_energy < 11) {
//             q2_min_analysis = 0.0;
//             q2_max_analysis = 12.0;
//             w_min_analysis = 1.0;
//             w_max_analysis = 2.5;
//           } else if (beam_energy < 24) {
//             q2_min_analysis = 2.0;
//             q2_max_analysis = 30.0;
//             w_min_analysis = 1.0;
//             w_max_analysis = 2.5;
//           }
//         }

//         // Update the condition to use the dynamically set limits // doing this later rather than earlier lessens the possibility of lost good events. **** change loop
//         if ((is_topology_excl && event->TwoPion_exclusive()) ||
//             (is_topology_mProt && event->TwoPion_missingProt()) ||
//             (is_topology_mPip && event->TwoPion_missingPip()) ||
//             (is_topology_mPim && event->TwoPion_missingPim())
//             ) {
//           if (event->W() > w_min_analysis && event->W() < w_max_analysis && 
//               event->Q2() > q2_min_analysis && event->Q2() < q2_max_analysis && 
//               event->weight() > 0.0) {
//             total_twopion_events++;
//           // }
//         // }
//       // }

//         // if ((is_topology_excl && event->TwoPion_exclusive()) ||
//         //     (is_topology_mProt && event->TwoPion_missingProt()) ||
//         //     (is_topology_mPip && event->TwoPion_missingPip()) ||
//         //     (is_topology_mPim && event->TwoPion_missingPim())) {
//         //   if (event->W() > 1.0 && event->W() < 2.5 && event->Q2() > 2.0 && event->Q2() < 30.0 && event->weight() > 0.0) {
            
//         //     total_twopion_events++;

//             csv_data row = {};

//             // if (cuts->IsProton(part)) {
//             //     c.pid_type = PROTON;
//             //     proton_cands.push_back(c);
//             // }
//             // else if (cuts->IsPip(part)) {
//             //     c.pid_type = PIP;
//             //     pip_cands.push_back(c);
//             // }
//             // else if (cuts->IsPim(part)) {
//             //     c.pid_type = PIM;
//             //     pim_cands.push_back(c);
//             // }

//             // if (csv_data::WRITE_ALL_CANDIDATES) {
//               for (auto &p : proton_cands) {
//                   csv_data row;
//                   row.event = data->event();
//                   row.run = data->run();

//                   row.part = p.part;
//                   row.rec_pid = p.rec_pid;
//                   row.mc_pid = p.mc_pid;

//                   row.pid_type = p.pid_type;
//                   row.is_misid = (p.mc_pid != -9999 && p.rec_pid != p.mc_pid);

//                   row.px = p.px;
//                   row.py = p.py;
//                   row.pz = p.pz;

//                   _sync->write(row);
//               }

//               for (auto &p : pip_cands) {
//                   csv_data row;
//                   row.event = data->event();
//                   row.run = data->run();

//                   row.part = p.part;
//                   row.rec_pid = p.rec_pid;
//                   row.mc_pid = p.mc_pid;

//                   row.pid_type = p.pid_type;
//                   row.is_misid = (p.mc_pid != -9999 && p.rec_pid != p.mc_pid);

//                   row.px = p.px;
//                   row.py = p.py;
//                   row.pz = p.pz;

//                   _sync->write(row);
//               }

//               for (auto &p : pim_cands) {
//                   csv_data row;
//                   row.event = data->event();
//                   row.run = data->run();

//                   row.part = p.part;
//                   row.rec_pid = p.rec_pid;
//                   row.mc_pid = p.mc_pid;

//                   row.pid_type = p.pid_type;
//                   row.is_misid = (p.mc_pid != -9999 && p.rec_pid != p.mc_pid);

//                   row.px = p.px;
//                   row.py = p.py;
//                   row.pz = p.pz;

//                   // if (!std::isfinite(c.px) || !std::isfinite(c.py) || !std::isfinite(c.pz)) continue;

//                   _sync->write(row);
//               }
//             // }
            
//             // if (WRITE_BEST_PAIR_ONLY) {

//             //   double best_sum = 1e9;
//             //   int best_p = -1;
//             //   int best_pi = -1;

//             //   for (auto &p : proton_cands) {
//             //       for (auto &pi : pip_cands) {

//             //           double dp_sum =
//             //               pow(p.px - pi.px, 2) +
//             //               pow(p.py - pi.py, 2) +
//             //               pow(p.pz - pi.pz, 2);

//             //           if (dp_sum < best_sum) {
//             //               best_sum = dp_sum;
//             //               best_p = p.part;
//             //               best_pi = pi.part;
//             //           }
//             //       }
//             //   }

//             //   if (best_p >= 0 && best_pi >= 0) {
//             //       output.pid_prot_rec = data->pid(best_p);
//             //       output.pid_pip_rec  = data->pid(best_pi);
//             //   }
//             // }


//             // ----- Reconstructed data output -----
//             // csv_data output; // moved this line up
//             // output.event = current_event;
//             output.run = data->run();
//             output.event = data->event();
//             output.w = event->W();
//             output.q2 = event->Q2();
//             output.weight_rec = event->weight(); // * 1e4;

//             output.mm2_mPim = event->MM2_mPim();
//             output.mm2_mPip = event->MM2_mPip();
//             output.mm2_mProt = event->MM2_mProt();
//             output.mm2_exclusive = event->MM2_exclusive();

//             output.pim_mom_miss = event->pim_momentum();
//             output.pim_mom_meas = event->pim_momentum_measured();
//             output.pip_mom_miss = event->pip_momentum();
//             output.pip_mom_meas = event->pip_momentum_measured();
//             output.prot_mom_miss = event->prot_momentum();
//             output.prot_mom_meas = event->prot_momentum_measured();
//             output.excl_mom = event->Mom_excl();

//             output.pim_theta_miss = event->pim_theta_lab();
//             output.pim_theta_meas = event->pim_theta_lab_measured();
//             output.pip_theta_miss = event->pip_theta_lab();
//             output.pip_theta_meas = event->pip_theta_lab_measured();
//             output.prot_theta_miss = event->prot_theta_lab();
//             output.prot_theta_meas = event->prot_theta_lab_measured();

//             output.prot_theta_angle_btwn_P = event->prot_theta_angle_btwn_P();
//             output.pip_theta_angle_btwn_P = event->pip_theta_angle_btwn_P();
//             output.pim_theta_angle_btwn_P = event->pim_theta_angle_btwn_P();
//           }
//         }
//       }
//     }

//     std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
//     std::cout << " total no of events = " << total << std::endl;
//     std::cout << " total no of twopion events = " << total_twopion_events << std::endl;

//     return num_of_events;
//   }

// #endif
