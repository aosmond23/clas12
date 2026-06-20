
#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <iostream>
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "branches.hpp"
#include "constants.hpp"
#include "physics.hpp"
#include <cstdlib>
#include <TRandom3.h>

class Reaction {
 protected:
  std::shared_ptr<Branches12> _data;

  double _beam_energy = 10.6041; // default beam energy

  Reaction() {
      const char* env = std::getenv("BEAM_E");
      if (env) {
          _beam_energy = std::atof(env);
      }
  }

  std::string _data_type;
  // double _beam_energy = beam_energy;
  std::unique_ptr<TLorentzVector> _beam;
  std::unique_ptr<TLorentzVector> _elec;
  std::unique_ptr<TLorentzVector> _gamma;
  std::unique_ptr<TLorentzVector> _target;
  std::unique_ptr<TLorentzVector> _prot;
  std::unique_ptr<TLorentzVector> _pip;
  std::unique_ptr<TLorentzVector> _pim;
  std::unique_ptr<TLorentzVector> _other;
  std::unique_ptr<TLorentzVector> _neutron;

  // std::vector<std::unique_ptr<TLorentzVector>> _prot;
  // std::vector<std::unique_ptr<TLorentzVector>> _pip;
  // std::vector<std::unique_ptr<TLorentzVector>> _pim;

  std::unique_ptr<TLorentzVector> _elecUnSmear;

  std::vector<int> _prot_indices;
  std::vector<int> _pip_indices;
  std::vector<int> _pim_indices;

  bool _mc = false;
  bool is_gen_data = false;
  bool is_rec_data = false;
  bool is_exp_data = false;

  bool _hasE = false;
  bool _hasP = false;
  bool _hasPip = false;
  bool _hasPim = false;
  bool _hasOther = false;
  bool _hasNeutron = false;

  short _numProt = 0;
  short _numPip = 0;
  short _numPim = 0;
  short _numPos = 0;
  short _numNeg = 0;
  short _numNeutral = 0;
  short _numOther = 0;

  short _sector = -1;

  float _MM = NAN;
  float _MM2_mPim = NAN;
  float _MM2_exclusive = NAN;
  float _excl_Energy = NAN;
  float _excl_Mom = NAN;;
  float _MM2_mPip = NAN;
  float _MM2_mProt = NAN;

  float _W = NAN;
  float _Q2 = NAN;
  double _emu_prime_mag2 = NAN;;
  double _emu_mag2 = NAN;;
  double _elec_energy = NAN;
  double _elec_mom = NAN;
  float _theta_e = NAN;
  float _elec_phi = NAN;

  float _inv_Ppip = NAN;
  float _inv_Ppim = NAN;
  float _inv_pip_pim = NAN;

  float _prot_status = NAN;
  float _pip_status = NAN; 
  float _pim_status = NAN;
  float _elec_status = NAN;

  int _sectorElec = -1;



  void SetElec();


 public:
  // Reaction(){};
  Reaction(const std::shared_ptr<Branches12> &data, float beam_energy, const std::string &data_type);
  ~Reaction();
  
  inline bool mc() { return _mc; }
    
  inline float weight() {
    if (_data_type == "rec" || _data_type == "gen") {
      return _data->mc_weight(); // / 1e4; // ***** search mc_weight and edit that / 1e4 too   
    }
    if (_data_type == "exp") {
      return 1.0;
    }
    return 0.0;
  }

  // ******************** add smearing function here **************************

  /// smearing fx's function
  void SmearingFunc(int part_id, int status_part, double p, double theta, double phi, double &pNew, double &thetaNew,
                    double &phiNew)
  {
          float syst_multi_fact = 1.0;

          // Constants
          const double pS1 = 0.0184291 - 0.0110083 * theta + 0.00227667 * pow(theta, 2) - 0.000140152 * pow(theta, 3) +
                              3.07424e-6 * pow(theta, 4);
          const double pR = 0.02 * sqrt(pow(pS1 * p, 2) + pow(0.02 * theta, 2));
          const double thetaR = 2.5 * sqrt(pow((0.004 * theta + 0.1) * (pow(p, 2) + 0.13957 * 0.13957) / pow(p, 2), 2));
          const double phiS1 = 0.85 - 0.015 * theta;
          const double phiS2 = 0.17 - 0.003 * theta;
          const double phiR = 3.5 * sqrt(pow(phiS1 * sqrt(pow(p, 2) + 0.13957 * 0.13957) / pow(p, 2), 2) + pow(phiS2, 2));

          // Generate new values
          if (part_id == ELECTRON)
          {
                  phiNew = phi + phiR * gRandom->Gaus(0, 1) * 0.4 * syst_multi_fact;       // * ((-0.4632) * w_val + (2.0038) * w_val + (-0.9035) * w_val); /// 0.4 was the origonal from pass1
                  thetaNew = theta + thetaR * gRandom->Gaus(0, 1) * 0.4 * syst_multi_fact; // * ((-0.4632) * w_val + (2.0038) * w_val + (-0.9035) * w_val);
                  pNew = p + pR * gRandom->Gaus(0, 1) * p * 0.4 * syst_multi_fact;         // * ((-0.4632) * w_val + (2.0038) * w_val + (-0.9035) * w_val);
          }
          else if (part_id == PROTON)
          {
                  double fact_cd = 0;
                  double fact_fd = 0;
                  double fact_cd1 = 0;
                  double fact_fd1 = 0;
                  if (status_part > 4000)
                  {
                          // fact_cd = (0.001948) * pow(p, 3) + (-0.038387) * pow(p, 2) + (0.181457) * p + (1.427647); ////// new pass2

                          fact_cd = (0.000821) * pow(p, 3) + (-0.016500) * pow(p, 2) + (0.103611) * p + (1.393237); ////// old
                          fact_cd1 = syst_multi_fact;
                          //(0.001536) * pow(p, 3) + (-0.024778) * pow(p, 2) + (0.119853) * p + (0.832939);

                          phiNew = phi + 1 / (fact_cd * fact_cd1) * phiR * gRandom->Gaus(0, 1);
                          thetaNew = theta + 1 / (fact_cd * fact_cd1) * thetaR * gRandom->Gaus(0, 1);
                          pNew = p + 1 / (fact_cd * fact_cd1) * pR * gRandom->Gaus(0, 1) * p;
                          // std::cout << "mom " << p << "prot fact_cd : " << 1 / fact_cd << std::endl;
                  }
                  else if (status_part <= 4000)
                  {
                          // fact_fd = (-0.000448) * pow(p, 3) + (0.009443) * pow(p, 2) + (-0.082798) * p + (1.947335); ///// new pass2
                          fact_fd = (0.000264) * pow(p, 3) + (-0.006454) * pow(p, 2) + (0.032683) * p + (1.658142); ////// old
                          fact_fd1 = syst_multi_fact;
                          //(0.000051) * pow(p, 3) + (-0.001569) * pow(p, 2) + (0.015891) * p + (0.966351);

                          phiNew = phi + 1 / (fact_fd * fact_fd1) * phiR * gRandom->Gaus(0, 1);
                          thetaNew = theta + 1 / (fact_fd * fact_fd1) * thetaR * gRandom->Gaus(0, 1);
                          pNew = p + 1 / (fact_fd * fact_fd1) * pR * gRandom->Gaus(0, 1) * p;

                          // std::cout << "mom " << p << "prot fact_fd : " << 1 / fact_fd << std::endl;
                  }
          }

          else if (part_id == PIP)
          {
                  double fact_cd = 0;
                  double fact_fd = 0;
                  double fact_cd1 = 0;
                  double fact_fd1 = 0;
                  if (status_part > 4000)
                  {
                          // fact_cd = (0.000891) * pow(p, 3) + (-0.011991) * pow(p, 2) + (0.006993) * p + (1.866801); //// new

                          fact_cd = (0.000981) * pow(p, 3) + (-0.016882) * pow(p, 2) + (0.046752) * p + (1.720426); /// old
                          fact_cd1 = syst_multi_fact;
                          //(-0.000104) * pow(p, 3) + (0.000998) * pow(p, 2) + (-0.008019) * p + (1.105314);

                          // std::cout << "mom " << p << "pip fact_cd : " << 1 / fact_cd << std::endl;

                          phiNew = phi + 1 / (fact_cd * fact_cd1) * phiR * gRandom->Gaus(0, 1);
                          thetaNew = theta + 1 / (fact_cd * fact_cd1) * thetaR * gRandom->Gaus(0, 1);
                          pNew = p + 1 / (fact_cd * fact_cd1) * pR * gRandom->Gaus(0, 1) * p;
                  }
                  else if (status_part <= 4000)
                  {
                          // fact_fd = (0.000369) * pow(p, 3) + (-0.012047) * pow(p, 2) + (0.080091) * p + (1.386573); // new

                          fact_fd = (0.000085) * pow(p, 3) + (-0.003096) * pow(p, 2) + (0.023553) * p + (1.509910); /// old
                          fact_fd1 = syst_multi_fact;
                          //(-0.000006) * pow(p, 3) + (-0.001310) * pow(p, 2) + (0.023171) * p + (0.890554);

                          // std::cout << "mom " << p << "pip fact_fd : " << 1 / fact_fd << std::endl;

                          phiNew = phi + 1 / (fact_fd * fact_fd1) * phiR * gRandom->Gaus(0, 1);
                          thetaNew = theta + 1 / (fact_fd * fact_fd1) * thetaR * gRandom->Gaus(0, 1);
                          pNew = p + 1 / (fact_fd * fact_fd1) * pR * gRandom->Gaus(0, 1) * p;
                  }
          }

          else if (part_id == PIM)
          {
                  double fact_cd = 0;
                  double fact_fd = 0;
                  double fact_cd1 = 0;
                  double fact_fd1 = 0;
                  if (status_part > 4000)
                  {
                          // fact_cd = (-0.001201) * pow(p, 3) + (0.014156) * pow(p, 2) + (-0.088521) * p + (2.054783); ////// new
                          fact_cd = (-0.001788) * pow(p, 3) + (0.025796) * pow(p, 2) + (-0.136577) * p + (2.007917); ////// old
                          fact_cd1 = syst_multi_fact;
                          //(-0.001327) * pow(p, 3) + (0.019826) * pow(p, 2) + (-0.097667) * p + (1.308904);
                          // std::cout << "mom " << p << "pim fact_cd : " << 1 / fact_cd << std::endl;

                          phiNew = phi + 1 / (fact_cd * fact_cd1) * phiR * gRandom->Gaus(0, 1);
                          thetaNew = theta + 1 / (fact_cd * fact_cd1) * thetaR * gRandom->Gaus(0, 1);
                          pNew = p + 1 / (fact_cd * fact_cd1) * pR * gRandom->Gaus(0, 1) * p;
                  }
                  else if (status_part <= 4000)
                  {
                          // fact_fd = (0.000062) * pow(p, 3) + (-0.005154) * pow(p, 2) + (0.072283) * p + (1.324110); //// new
                          fact_fd = (0.000760) * pow(p, 3) + (-0.021295) * pow(p, 2) + (0.171180) * p + (1.238299); ////// old
                          fact_fd1 = syst_multi_fact;
                          //(0.000249) * pow(p, 3) + (-0.007461) * pow(p, 2) + (0.067686) * p + (0.817653);

                          // std::cout << "mom " << p << "pim fact_cd : " << 1 / fact_fd << std::endl;

                          phiNew = phi + 1 / (fact_fd * fact_fd1) * phiR * gRandom->Gaus(0, 1);
                          thetaNew = theta + 1 / (fact_fd * fact_fd1) * thetaR * gRandom->Gaus(0, 1);
                          pNew = p + 1 / (fact_fd * fact_fd1) * pR * gRandom->Gaus(0, 1) * p;
                  }
          }
  }
        // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  // **********************************************

  // Check lists when you swich from mc to exp or vice-versa
  // 3. from if (data->mc_npart() < 1) to all particle set up im mc events.

  // momentum correction
  void SetMomCorrElec();
  double dpp(float px, float py, float pz, int sec_mom_corr, int ivec);
  double elec_mom();
  double elec_E();
  double elec_theta();
  double elec_phi();
  double elec_prime_mass2();
  double elec_mass2();

  void SetProton(int i);
  void SetPip(int i);
  void SetPim(int i);
  void SetOther(int i);
  void SetNeutron(int i);

  // missingPim
  float pim_momentum();
  float pim_theta_lab();
  float pim_Phi_lab();
  float pim_momentum_measured();
  float pim_theta_lab_measured();
  float pim_Phi_lab_measured();
  float pim_theta_angle_btwn_P();

  // missingPip
  float pip_momentum();
  float pip_theta_lab();
  float pip_Phi_lab();
  float pip_momentum_measured();
  float pip_theta_lab_measured();
  float pip_Phi_lab_measured();
  float pip_theta_angle_btwn_P();

  // missingProt
  float prot_momentum();
  float prot_theta_lab();
  float prot_Phi_lab();
  float prot_momentum_measured();
  float prot_theta_lab_measured();
  float prot_Phi_lab_measured();
  float prot_theta_angle_btwn_P();

  void CalcMissMass();
  float MM();
  float MM2_mPim();
  float MM2_exclusive();
  float Energy_excl();
  float Mom_excl();
  float MM2_mPip();
  float MM2_mProt();

  float w_hadron();
  float w_difference();

  void invMassPpip(const TLorentzVector &prot, const TLorentzVector &pip);
  void invMassPpim(const TLorentzVector &prot, const TLorentzVector &pip);
  void invMasspippim(const TLorentzVector &prot, const TLorentzVector &pip);

  float inv_Ppip();
  float inv_Ppim();
  float inv_pip_pim();

  inline float W() { return _W; }
  inline float Q2() { return _Q2; }

  inline short sec() { return _data->dc_sec(0); }
  inline int det() { return abs(_data->status(0) / 1000); }

  // // ----------------------- FOR EXCLUSIVE TOPO -----------------------
  // inline bool TwoPion_missingPim() {
  //   bool _channelTwoPi = true;
  //   // _channelTwoPi &= ((_numProt == 1 && _numPip == 1) && (_hasE && _hasP  && _hasPip) && !TwoPion_exclusive());
  //   _channelTwoPi &= ((_numProt >= 1 && _numPip >= 1) && (_hasE && _hasP && _hasPip)); // &&!_hasPim));
  //   return _channelTwoPi;
  // }
  // inline bool TwoPion_exclusive() {
  //   bool _channelTwoPi_excl = true;
  //   // _channelTwoPi_excl &= ((_numProt == 1 && _numPip == 1 && _numPim == 1) &&
  //   _channelTwoPi_excl &= ((_numProt >= 1 && _numPip >= 1 && _numPim >= 1) &&
  //                          (_hasE && _hasP && _hasPip && _hasPim /*&& !_hasNeutron && !_hasOther*/));
  //   return _channelTwoPi_excl;
  // }
  // inline bool TwoPion_missingPip() {
  //   bool _channelTwoPi_mpip = true;
  //   _channelTwoPi_mpip &=
  //       // ((_numProt == 1 && _numPim == 1) && (_hasE && _hasP && _hasPim /*&&!_hasPip && !_hasNeutron && !_hasOther*/));
  //       ((_numProt >= 1 && _numPim >= 1) && (_hasE && _hasP && _hasPim /*&&!_hasPip && !_hasNeutron && !_hasOther*/));
  //   return _channelTwoPi_mpip;
  // }
  // inline bool TwoPion_missingProt() {
  //   bool _channelTwoPi_mprot = true;
  //   _channelTwoPi_mprot &=
  //       // ((_numPip == 1 && _numPim == 1) && (_hasE && _hasPip && _hasPim /*&&!_hasP  && !_hasOther*/));
  //       ((_numPip >= 1 && _numPim >= 1) && (_hasE && _hasPip && _hasPim /*&&!_hasP  && !_hasOther*/));
  //   return _channelTwoPi_mprot;
  // }

// ----------------------- FOR MISSING TOPOS -----------------------

  inline bool TwoPion_missingPim() {
    bool _channelTwoPi = true;
    // _channelTwoPi &= ((_numProt == 1 && _numPip == 1) && (_hasE && _hasP  && _hasPip) && !TwoPion_exclusive());
    _channelTwoPi &= ((_numProt >= 1 && _numPip >= 1) && (_hasE && _hasP && _hasPip &&!_hasPim));
    return _channelTwoPi;
  }
  inline bool TwoPion_exclusive() {
    bool _channelTwoPi_excl = true;
    // _channelTwoPi_excl &= ((_numProt == 1 && _numPip == 1 && _numPim == 1) &&
    _channelTwoPi_excl &= ((_numProt >= 1 && _numPip >= 1 && _numPim >= 1) &&
                           (_hasE && _hasP && _hasPip && _hasPim /*&& !_hasNeutron && !_hasOther*/));
    return _channelTwoPi_excl;
  }
  inline bool TwoPion_missingPip() {
    bool _channelTwoPi_mpip = true;
    _channelTwoPi_mpip &=
        // ((_numProt == 1 && _numPim == 1) && (_hasE && _hasP && _hasPim /*&&!_hasPip && !_hasNeutron && !_hasOther*/));
        ((_numProt >= 1 && _numPim >= 1) && (_hasE && _hasP && _hasPim &&!_hasPip /*&& !_hasNeutron && !_hasOther*/));
    return _channelTwoPi_mpip;
  }
  inline bool TwoPion_missingProt() {
    bool _channelTwoPi_mprot = true;
    _channelTwoPi_mprot &=
        // ((_numPip == 1 && _numPim == 1) && (_hasE && _hasPip && _hasPim /*&&!_hasP  && !_hasOther*/));
        ((_numPip >= 1 && _numPim >= 1) && (_hasE && _hasPip && _hasPim &&!_hasP  /*&& !_hasOther*/));
    return _channelTwoPi_mprot;
  }

  // inline bool TwoPion_exclusive() {
  //   bool _channelTwoPi_excl = true;
  //   _channelTwoPi_excl &= (_hasE && _hasP && _hasPip && _hasPim &&
  //           _numProt >= 1 && _numPip  >= 1 && _numPim  >= 1);
  //   return _channelTwoPi_excl;
  // }

  // inline bool TwoPion_missingProt() {
  //   bool _channelTwoPi_mprot = true;
  //   _channelTwoPi_mprot &= (_hasE &&
  //           _numProt < 1 && _numPip  >= 1 && _numPim  >= 1);
  //   return _channelTwoPi_mprot;
  // }

  // inline bool TwoPion_missingPip() {
  //   bool _channelTwoPi_mpip = true;
  //   _channelTwoPi_mpip &= (_hasE && 
  //           _numProt >= 1 && _numPip  < 1 && _numPim  >= 1);
  //   return _channelTwoPi_mpip;
  // }

  // inline bool TwoPion_missingPim() {
  //   bool _channelTwoPi_mpim = true;
  //   _channelTwoPi_mpim &= (_hasE &&
  //           _numProt >= 1 && _numPip  >= 1 && _numPim  < 1);
  //   return _channelTwoPi_mpim;
  // }

  const TLorentzVector &e_mu() { return *_beam; }
  const TLorentzVector &e_mu_prime() { return *_elec; }
  const TLorentzVector &gamma() { return *_gamma; }
};

class MCReaction : public Reaction {
 private:
  float _weight_mc = NAN;
  float _W_mc = NAN;
  float _Q2_mc = NAN;

  std::unique_ptr<TLorentzVector> _elec_mc;
  std::unique_ptr<TLorentzVector> _gamma_mc;
  std::unique_ptr<TLorentzVector> _prot_mc;
  std::unique_ptr<TLorentzVector> _pip_mc;
  std::unique_ptr<TLorentzVector> _pim_mc;
  std::unique_ptr<TLorentzVector> _other_mc;

  float _MM2_exclusive_mc = NAN;
  float _excl_Energy_mc = NAN;

 public:
  void SetMCProton(int i);
  void SetMCPip(int i);
  void SetMCPim(int i);
  void SetMCOther(int i);

  MCReaction(const std::shared_ptr<Branches12> &data, float beam_energy, const std::string &data_type);
  void SetMCElec();
  void CalcMissMass_mc();

  inline float weight() { return _data->mc_weight(); } // / 1e4; }
  inline float W_mc() { return _W_mc; }
  inline float Q2_mc() { return _Q2_mc; }
  inline float MM2_exclusive_mc() const { return _MM2_exclusive_mc; }

  float elec_mom_mc_gen();
  float pim_mom_mc_gen();
  float pip_mom_mc_gen();
  float prot_mom_mc_gen();

  float elec_E_mc_gen();

  float MCinv_Ppip();
  float MCinv_Ppim();
  float MCinv_pip_pim();

  float elec_theta_mc_gen();
  float pim_theta_mc_gen();
  float pip_theta_mc_gen();
  float prot_theta_mc_gen();

  float elec_phi_mc_gen();
  float pim_phi_mc_gen();
  float pip_phi_mc_gen();
  float prot_phi_mc_gen();

  // void CalcMissMass_mc();

};

#endif
