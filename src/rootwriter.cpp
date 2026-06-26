#include "rootwriter.hpp"
#include <iostream>

RootWriter::RootWriter(const std::string& filename, bool isGenerated)
{
    generated = isGenerated;
    std::cout << "Constructing RootWriter " << this << std::endl;
    file = new TFile(filename.c_str(),"RECREATE");
    file->cd();
    tree = new TTree("events","Analysis tree");
    tree->SetDirectory(file);
    createBranches();
    // if (generated)
    //     createGenBranches();
    // else   
    //     createRecBranches();
}

RootWriter::~RootWriter()
{
    std::cout << "Destroying RootWriter " << this << std::endl;
    std::cout << "Entries = " << tree->GetEntries() << std::endl;

    file->cd();
    tree->Write("", TObject::kOverwrite);
    // tree->Write();
    // file->Write();
    file->Close();
    delete file;
}

void RootWriter::Fill(const csv_data& data)
{
    buffer = csv_data{};  // 🔥 reset everything
    buffer = data;
    tree->Fill();
}

// void RootWriter::createGenBranches() {

//     // Generated output branches
//     tree->Branch("run",        &buffer.mc_run);
//     tree->Branch("event",      &buffer.mc_event);
//     tree->Branch("w_mc",       &buffer.w_mc);
//     tree->Branch("q2_mc",      &buffer.q2_mc);
//     tree->Branch("weight_mc",     &buffer.weight_gen);

//     tree->Branch("pid_prot_mc", &buffer.pid_prot_mc);
//     tree->Branch("pid_pip_mc",  &buffer.pid_pip_mc);
//     tree->Branch("pid_pim_mc",  &buffer.pid_pim_mc);
// }

// void RootWriter::createRecBranches() {

//     // Reconstructed output branches
//     tree->Branch("event", &buffer.event);

//     tree->Branch("w", &buffer.w);
//     tree->Branch("q2", &buffer.q2);
//     tree->Branch("weight", &buffer.weight_rec);

//     tree->Branch("pid_prot_rec", &buffer.pid_prot_rec);
//     tree->Branch("pid_pip_rec",  &buffer.pid_pip_rec);
//     tree->Branch("pid_pim_rec",  &buffer.pid_pim_rec);

//     tree->Branch("pid_prot_mc", &buffer.pid_prot_mc);
//     tree->Branch("pid_pip_mc",  &buffer.pid_pip_mc);
//     tree->Branch("pid_pim_mc",  &buffer.pid_pim_mc);

//     tree->Branch("mm2_mPim", &buffer.mm2_mPim);
//     tree->Branch("mm2_mPip", &buffer.mm2_mPip);
//     tree->Branch("mm2_mProt", &buffer.mm2_mProt);
//     tree->Branch("mm2_excl", &buffer.mm2_exclusive);

//     tree->Branch("pim_mom_miss", &buffer.pim_mom_miss);
//     tree->Branch("pim_mom_meas", &buffer.pim_mom_meas);
//     tree->Branch("pip_mom_miss", &buffer.pip_mom_miss);
//     tree->Branch("pip_mom_meas", &buffer.pip_mom_meas);
//     tree->Branch("prot_mom_miss", &buffer.prot_mom_miss);
//     tree->Branch("prot_mom_meas", &buffer.prot_mom_meas);
//     tree->Branch("excl_mom", &buffer.excl_mom);

//     tree->Branch("pim_theta_miss", &buffer.pim_theta_miss);
//     tree->Branch("pim_theta_meas", &buffer.pim_theta_meas);
//     tree->Branch("pip_theta_miss", &buffer.pip_theta_miss);
//     tree->Branch("pip_theta_meas", &buffer.pip_theta_meas);
//     tree->Branch("prot_theta_miss", &buffer.prot_theta_miss);
//     tree->Branch("prot_theta_meas", &buffer.prot_theta_meas);

//     tree->Branch("pim_theta_angle_btwn_P", &buffer.pim_theta_angle_btwn_P);
//     tree->Branch("pip_theta_angle_btwn_P", &buffer.pip_theta_angle_btwn_P);
//     tree->Branch("prot_theta_angle_btwn_P", &buffer.prot_theta_angle_btwn_P);
// }

void RootWriter::createBranches() {

    if (generated) {

        // Generated output branches
        tree->Branch("run",        &buffer.mc_run);
        tree->Branch("event_mc",      &buffer.mc_event);
        tree->Branch("w_mc",       &buffer.w_mc);
        tree->Branch("q2_mc",      &buffer.q2_mc);
        tree->Branch("weight_mc",     &buffer.weight_gen);

        tree->Branch("pid_prot_mc", &buffer.pid_prot_mc);
        tree->Branch("pid_pip_mc",  &buffer.pid_pip_mc);
        tree->Branch("pid_pim_mc",  &buffer.pid_pim_mc);

    } else {

        // Reconstructed output branches
        tree->Branch("event", &buffer.event);

        tree->Branch("w", &buffer.w);
        tree->Branch("q2", &buffer.q2);
        tree->Branch("weight", &buffer.weight_rec);

        tree->Branch("pid_prot_rec", &buffer.pid_prot_rec);
        tree->Branch("pid_pip_rec",  &buffer.pid_pip_rec);
        tree->Branch("pid_pim_rec",  &buffer.pid_pim_rec);

        tree->Branch("pid_prot_mc", &buffer.pid_prot_mc);
        tree->Branch("pid_pip_mc",  &buffer.pid_pip_mc);
        tree->Branch("pid_pim_mc",  &buffer.pid_pim_mc);

        tree->Branch("mm2_mPim", &buffer.mm2_mPim);
        tree->Branch("mm2_mPip", &buffer.mm2_mPip);
        tree->Branch("mm2_mProt", &buffer.mm2_mProt);
        tree->Branch("mm2_excl", &buffer.mm2_exclusive);

        tree->Branch("pim_mom_miss", &buffer.pim_mom_miss);
        tree->Branch("pim_mom_meas", &buffer.pim_mom_meas);
        tree->Branch("pip_mom_miss", &buffer.pip_mom_miss);
        tree->Branch("pip_mom_meas", &buffer.pip_mom_meas);
        tree->Branch("prot_mom_miss", &buffer.prot_mom_miss);
        tree->Branch("prot_mom_meas", &buffer.prot_mom_meas);
        tree->Branch("excl_mom", &buffer.excl_mom);

        tree->Branch("pim_theta_miss", &buffer.pim_theta_miss);
        tree->Branch("pim_theta_meas", &buffer.pim_theta_meas);
        tree->Branch("pip_theta_miss", &buffer.pip_theta_miss);
        tree->Branch("pip_theta_meas", &buffer.pip_theta_meas);
        tree->Branch("prot_theta_miss", &buffer.prot_theta_miss);
        tree->Branch("prot_theta_meas", &buffer.prot_theta_meas);

        tree->Branch("pim_theta_angle_btwn_P", &buffer.pim_theta_angle_btwn_P);
        tree->Branch("pip_theta_angle_btwn_P", &buffer.pip_theta_angle_btwn_P);
        tree->Branch("prot_theta_angle_btwn_P", &buffer.prot_theta_angle_btwn_P);
    }
}