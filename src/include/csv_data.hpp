#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>
#include <iomanip> // For std::setprecision

struct csv_data {
  // int event;  // Added for event number 
  short electron_sector;
  float sf;
  float w;
  float w_mc;
  double elec_prime_m2;
  double elec_m2;
  double elec_energy_rec;
  float q2_mc;
  float w_had;
  float w_diff;
  float elec_mom_rec;
  float elec_theta_rec;
  float elec_phi_rec;

  float energy_x_mu;
  float mom_x_mu;

  float q2;

  float gen_elec_E;
  float gen_elec_mom;
  float gen_elec_theta;
  float gen_elec_phi;

  float gen_pim_mom;
  float gen_pim_theta;
  float gen_pim_phi;

  float gen_pip_mom;
  float gen_pip_theta;
  float gen_pip_phi;

  float gen_prot_mom;
  float gen_prot_theta;
  float gen_prot_phi;

  float vertex_x;
  float vertex_y;
  float vertex_z;

  float vertex_had[3][3];

  float weight_gen;
  float weight_rec;

  float mm2_mPim;
  float mm2_mPip;
  float mm2_mProt;
  float mm2_exclusive;

  float pim_mom_miss;
  float pim_mom_meas;
  float pip_mom_miss;
  float pip_mom_meas;
  float prot_mom_miss;
  float prot_mom_meas;
  float excl_mom;

  float pim_theta_miss;
  float pim_theta_meas;
  float pip_theta_miss;
  float pip_theta_meas;
  float prot_theta_miss;
  float prot_theta_meas;

  float prot_theta_angle_btwn_P;
  float pip_theta_angle_btwn_P;
  float pim_theta_angle_btwn_P;

  int status_Elec;
  int status_Pim;
  int status_Pip;
  int status_Prot;

  int pid_prot_rec;
  int pid_pip_rec;
  int pid_pim_rec;

  int pid_prot_mc;
  int pid_pip_mc;
  int pid_pim_mc;

  int run = 0;
  int event = 0;

  int mc_run = 0;
  int mc_event = 0;

  int n_prot_rec = 0;
  int n_pip_rec = 0;
  int n_pim_rec = 0;

  int n_prot_mc = 0;
  int n_pip_mc = 0;
  int n_pim_mc = 0;

  // inline static bool use_thrown_pid = true;

  // Declare the static flag to choose between generated and reconstructed data
  static bool isGenerated;

  // Static function to set the flag based on the filename
  static void setOutputType(const std::string &filename) {
    if (filename.find("gen") != std::string::npos) {
      isGenerated = true;  // Generated data
    } else if (filename.find("rec") != std::string::npos || filename.find("exp") != std::string::npos) {
      isGenerated = false; // Reconstructed or experimental data
    } else {
      // Default case or error handling
      isGenerated = true;  // Default to generated if no keyword is found
    }
  }

  // Static functions can be called without making a new struct
  static std::string header() {
    if (isGenerated) {
      // ----- Generated -----
      return "run,event,w_mc,q2_mc,weight,"
                "pid_prot_mc,pid_pip_mc,pid_pim_mc";
    } else {
      // ----- Reconstructed -----
      return "run,event,w,q2,weight,"
                "pid_prot_rec,pid_pip_rec,pid_pim_rec,"
                "pid_prot_mc,pid_pip_mc,pid_pim_mc,"
                "mm2_mPim,mm2_mPip,mm2_mProt,mm2_excl,"
                "pim_mom_miss,pim_mom_meas,pip_mom_miss,pip_mom_meas,prot_mom_miss,prot_mom_meas,excl_mom,"
                "pim_theta_miss,pim_theta_meas,pip_theta_miss,pip_theta_meas,prot_theta_miss,prot_theta_meas,"
                "pim_theta_angle_btwn_P,pip_theta_angle_btwn_P,prot_theta_angle_btwn_P";
    }
  }

  friend std::ostream& operator<<(std::ostream& os, const csv_data& data) {
    os << std::setprecision(7);

    if (isGenerated) {
      // ----- Generated -----
      os << data.mc_run << ",";
      os << data.mc_event << ",";
      os << data.w_mc << ",";
      os << data.q2_mc << ",";
      os << data.weight_gen << ",";

      os << data.pid_prot_mc << ",";
      os << data.pid_pip_mc << ",";
      os << data.pid_pim_mc << ",";
    } else {
      // ----- Reconstructed -----
      os << data.run << ",";
      os << data.event << ",";
      os << data.w << ",";
      os << data.q2 << ",";
      os << data.weight_rec << ",";

      os << data.pid_prot_rec << ",";
      os << data.pid_pip_rec << ",";
      os << data.pid_pim_rec << ",";

      os << data.pid_prot_mc << ",";
      os << data.pid_pip_mc << ",";
      os << data.pid_pim_mc << ",";

      os << data.mm2_mPim << ",";
      os << data.mm2_mPip << ",";
      os << data.mm2_mProt << ",";
      os << data.mm2_exclusive << ",";

      os << data.pim_mom_miss << ",";
      os << data.pim_mom_meas << ",";
      os << data.pip_mom_miss << ",";
      os << data.pip_mom_meas << ",";
      os << data.prot_mom_miss << ",";
      os << data.prot_mom_meas << ",";
      os << data.excl_mom << ",";

      os << data.pim_theta_miss << ",";
      os << data.pim_theta_meas << ",";
      os << data.pip_theta_miss << ",";
      os << data.pip_theta_meas << ",";
      os << data.prot_theta_miss << ",";
      os << data.prot_theta_meas << ",";

      os << data.pim_theta_angle_btwn_P << ",";
      os << data.pip_theta_angle_btwn_P << ",";
      os << data.prot_theta_angle_btwn_P << ",";
    }

    return os;
  };
};

#endif
