
#ifndef BRANCHES_H
#define BRANCHES_H
#include <iostream>
#include <vector>
#include "TChain.h"
#include "constants.hpp"

using v_int = std::vector<int> *;
using v_float = std::vector<float> *;

class Branches12 {
 private:
  std::shared_ptr<TChain> _tree;

  // Gotten from root's MakeClass
  Int_t _NRUN;
  Int_t _NEVENT;
  Float_t _EVNTime;
  Int_t _TYPE;
  Int_t _TRG;
  Float_t _BCG;
  Float_t _STTime;
  Float_t _RFTime;
  Int_t _Helic;
  Int_t _EvCAT;
  Int_t _NPGP;
  Double_t _LT;
  Float_t _PTIME;

  Int_t _run;
  Int_t _event;
  Int_t _unixtime;
  Float_t _trigger;
  Float_t _timestamp;
  Int_t _type;
  Int_t _mode;
  Float_t _torus;
  Float_t _solenoid;
  Int_t _category;
  Int_t _topology;
  Float_t _beamCharge;
  Double_t _liveTime;
  Float_t _startTime;
  Int_t _helicity;
  Int_t _helicityRaw;
  Float_t _procTime;
  Int_t _hel_run;
  Int_t _hel_event;
  Float_t _hel_timestamp;
  Int_t _hel_helicity;
  Int_t _hel_helicityRaw;
  Int_t _hel_pair;
  Int_t _hel_pattern;
  Int_t _hel_status;

  bool _is_mc;
  int _mc_run;
  int _mc_npart = 0;
  int _mc_event;
  int _mc_type;
  float _mc_helicity;
  float _mc_weight;

  v_int _mc_pid;
  v_float _mc_px;
  v_float _mc_py;
  v_float _mc_pz;
  v_float _mc_vx;
  v_float _mc_vy;
  v_float _mc_vz;
  v_float _mc_vt;

  v_int _pid;
  v_float _p;
  v_float _p2;
  v_float _px;
  v_float _py;
  v_float _pz;
  v_float _vx;
  v_float _vy;
  v_float _vz;
  v_int _charge;
  v_float _beta;
  v_float _chi2pid;
  v_int _status;
  v_int _dc_sec;
  v_float _dc_r1_x;
  v_float _dc_r1_y;
  v_float _dc_r1_z;
  v_float _dc_r2_x;
  v_float _dc_r2_y;
  v_float _dc_r2_z;
  v_float _dc_r3_x;
  v_float _dc_r3_y;
  v_float _dc_r3_z;
  v_float _cvt_x;
  v_float _cvt_y;
  v_float _cvt_z;
  v_float _fmt_x;
  v_float _fmt_y;
  v_float _fmt_z;
  v_float _ec_tot_energy;
  v_float _ec_pcal_energy;
  v_int _ec_pcal_sec;
  v_float _ec_pcal_time;
  v_float _ec_pcal_path;
  v_float _ec_pcal_x;
  v_float _ec_pcal_y;
  v_float _ec_pcal_z;
  v_float _ec_pcal_hx;
  v_float _ec_pcal_hy;
  v_float _ec_pcal_hz;
  v_float _ec_pcal_lu;
  v_float _ec_pcal_lv;
  v_float _ec_pcal_lw;
  v_float _ec_pcal_du;
  v_float _ec_pcal_dv;
  v_float _ec_pcal_dw;
  v_float _ec_pcal_m2u;
  v_float _ec_pcal_m2v;
  v_float _ec_pcal_m2w;
  v_float _ec_pcal_m3u;
  v_float _ec_pcal_m3v;
  v_float _ec_pcal_m3w;
  v_float _ec_ecin_energy;
  v_int _ec_ecin_sec;
  v_float _ec_ecin_time;
  v_float _ec_ecin_path;
  v_float _ec_ecin_x;
  v_float _ec_ecin_y;
  v_float _ec_ecin_z;
  v_float _ec_ecin_hx;
  v_float _ec_ecin_hy;
  v_float _ec_ecin_hz;
  v_float _ec_ecin_lu;
  v_float _ec_ecin_lv;
  v_float _ec_ecin_lw;
  v_float _ec_ecin_du;
  v_float _ec_ecin_dv;
  v_float _ec_ecin_dw;
  v_float _ec_ecin_m2u;
  v_float _ec_ecin_m2v;
  v_float _ec_ecin_m2w;
  v_float _ec_ecin_m3u;
  v_float _ec_ecin_m3v;
  v_float _ec_ecin_m3w;
  v_float _ec_ecout_energy;
  v_int _ec_ecout_sec;
  v_float _ec_ecout_time;
  v_float _ec_ecout_path;
  v_float _ec_ecout_x;
  v_float _ec_ecout_y;
  v_float _ec_ecout_z;
  v_float _ec_ecout_hx;
  v_float _ec_ecout_hy;
  v_float _ec_ecout_hz;
  v_float _ec_ecout_lu;
  v_float _ec_ecout_lv;
  v_float _ec_ecout_lw;
  v_float _ec_ecout_du;
  v_float _ec_ecout_dv;
  v_float _ec_ecout_dw;
  v_float _ec_ecout_m2u;
  v_float _ec_ecout_m2v;
  v_float _ec_ecout_m2w;
  v_float _ec_ecout_m3u;
  v_float _ec_ecout_m3v;
  v_float _ec_ecout_m3w;
  v_float _cc_nphe_tot;
  v_int _cc_ltcc_sec;
  v_float _cc_ltcc_nphe;
  v_float _cc_ltcc_time;
  v_float _cc_ltcc_path;
  v_float _cc_ltcc_theta;
  v_float _cc_ltcc_phi;
  v_float _cc_ltcc_x;
  v_float _cc_ltcc_y;
  v_float _cc_ltcc_z;
  v_int _cc_htcc_sec;
  v_float _cc_htcc_nphe;
  v_float _cc_htcc_time;
  v_float _cc_htcc_path;
  v_float _cc_htcc_theta;
  v_float _cc_htcc_phi;
  v_float _cc_htcc_x;
  v_float _cc_htcc_y;
  v_float _cc_htcc_z;
  v_int _cc_rich_sec;
  v_float _cc_rich_nphe;
  v_float _cc_rich_time;
  v_float _cc_rich_path;
  v_float _cc_rich_theta;
  v_float _cc_rich_phi;
  v_float _cc_rich_x;
  v_float _cc_rich_y;
  v_float _cc_rich_z;
  v_int _sc_ftof_1a_sec;
  v_float _sc_ftof_1a_time;
  v_float _sc_ftof_1a_path;
  v_float _sc_ftof_1a_energy;
  v_int _sc_ftof_1a_component;
  v_float _sc_ftof_1a_x;
  v_float _sc_ftof_1a_y;
  v_float _sc_ftof_1a_z;
  v_float _sc_ftof_1a_hx;
  v_float _sc_ftof_1a_hy;
  v_float _sc_ftof_1a_hz;
  v_int _sc_ftof_1b_sec;
  v_float _sc_ftof_1b_time;
  v_float _sc_ftof_1b_path;
  v_float _sc_ftof_1b_energy;
  v_int _sc_ftof_1b_component;
  v_float _sc_ftof_1b_x;
  v_float _sc_ftof_1b_y;
  v_float _sc_ftof_1b_z;
  v_float _sc_ftof_1b_hx;
  v_float _sc_ftof_1b_hy;
  v_float _sc_ftof_1b_hz;
  v_int _sc_ftof_2_sec;
  v_float _sc_ftof_2_time;
  v_float _sc_ftof_2_path;
  v_float _sc_ftof_2_energy;
  v_int _sc_ftof_2_component;
  v_float _sc_ftof_2_x;
  v_float _sc_ftof_2_y;
  v_float _sc_ftof_2_z;
  v_float _sc_ftof_2_hx;
  v_float _sc_ftof_2_hy;
  v_float _sc_ftof_2_hz;
  v_float _sc_ctof_time;
  v_float _sc_ctof_path;
  v_float _sc_ctof_energy;
  v_int _sc_ctof_component;
  v_float _sc_ctof_x;
  v_float _sc_ctof_y;
  v_float _sc_ctof_z;
  v_float _sc_ctof_hx;
  v_float _sc_ctof_hy;
  v_float _sc_ctof_hz;
  v_float _sc_cnd_time;
  v_float _sc_cnd_path;
  v_float _sc_cnd_energy;
  v_int _sc_cnd_component;
  v_float _sc_cnd_x;
  v_float _sc_cnd_y;
  v_float _sc_cnd_z;
  v_float _sc_cnd_hx;
  v_float _sc_cnd_hy;
  v_float _sc_cnd_hz;
  v_float _ft_cal_energy;
  v_float _ft_cal_time;
  v_float _ft_cal_path;
  v_float _ft_cal_x;
  v_float _ft_cal_y;
  v_float _ft_cal_z;
  v_float _ft_cal_dx;
  v_float _ft_cal_dy;
  v_float _ft_cal_radius;
  v_float _ft_hodo_energy;
  v_float _ft_hodo_time;
  v_float _ft_hodo_path;
  v_float _ft_hodo_x;
  v_float _ft_hodo_y;
  v_float _ft_hodo_z;
  v_float _ft_hodo_dx;
  v_float _ft_hodo_dy;
  v_float _ft_hodo_radius;

  // List of branches
  TBranch *b_run;                   //!
  TBranch *b_event;                 //!
  TBranch *b_unixtime;              //!
  TBranch *b_trigger;               //!
  TBranch *b_timestamp;             //!
  TBranch *b_type;                  //!
  TBranch *b_mode;                  //!
  TBranch *b_torus;                 //!
  TBranch *b_solenoid;              //!
  TBranch *b_category;              //!
  TBranch *b_topology;              //!
  TBranch *b_beamCharge;            //!
  TBranch *b_liveTime;              //!
  TBranch *b_startTime;             //!
  TBranch *b_RFTime;                //!
  TBranch *b_helicity;              //!
  TBranch *b_helicityRaw;           //!
  TBranch *b_procTime;              //!
  TBranch *b_hel_run;               //!
  TBranch *b_hel_event;             //!
  TBranch *b_hel_timestamp;         //!
  TBranch *b_hel_helicity;          //!
  TBranch *b_hel_helicityRaw;       //!
  TBranch *b_hel_pair;              //!
  TBranch *b_hel_pattern;           //!
  TBranch *b_hel_status;            //!
  TBranch *b_pid;                   //!
  TBranch *b_p;                     //!
  TBranch *b_p2;                    //!
  TBranch *b_px;                    //!
  TBranch *b_py;                    //!
  TBranch *b_pz;                    //!
  TBranch *b_vx;                    //!
  TBranch *b_vy;                    //!
  TBranch *b_vz;                    //!
  TBranch *b_vt;                    //!
  TBranch *b_charge;                //!
  TBranch *b_beta;                  //!
  TBranch *b_chi2pid;               //!
  TBranch *b_status;                //!
  TBranch *b_dc_sec;                //!
  TBranch *b_dc_r1_x;               //!
  TBranch *b_dc_r1_y;               //!
  TBranch *b_dc_r1_z;               //!
  TBranch *b_dc_r2_x;               //!
  TBranch *b_dc_r2_y;               //!
  TBranch *b_dc_r2_z;               //!
  TBranch *b_dc_r3_x;               //!
  TBranch *b_dc_r3_y;               //!
  TBranch *b_dc_r3_z;               //!
  TBranch *b_cvt_x;                 //!
  TBranch *b_cvt_y;                 //!
  TBranch *b_cvt_z;                 //!
  TBranch *b_fmt_x;                 //!
  TBranch *b_fmt_y;                 //!
  TBranch *b_fmt_z;                 //!
  TBranch *b_ec_tot_energy;         //!
  TBranch *b_ec_pcal_energy;        //!
  TBranch *b_ec_pcal_sec;           //!
  TBranch *b_ec_pcal_time;          //!
  TBranch *b_ec_pcal_path;          //!
  TBranch *b_ec_pcal_x;             //!
  TBranch *b_ec_pcal_y;             //!
  TBranch *b_ec_pcal_z;             //!
  TBranch *b_ec_pcal_hx;            //!
  TBranch *b_ec_pcal_hy;            //!
  TBranch *b_ec_pcal_hz;            //!
  TBranch *b_ec_pcal_lu;            //!
  TBranch *b_ec_pcal_lv;            //!
  TBranch *b_ec_pcal_lw;            //!
  TBranch *b_ec_pcal_du;            //!
  TBranch *b_ec_pcal_dv;            //!
  TBranch *b_ec_pcal_dw;            //!
  TBranch *b_ec_pcal_m2u;           //!
  TBranch *b_ec_pcal_m2v;           //!
  TBranch *b_ec_pcal_m2w;           //!
  TBranch *b_ec_pcal_m3u;           //!
  TBranch *b_ec_pcal_m3v;           //!
  TBranch *b_ec_pcal_m3w;           //!
  TBranch *b_ec_ecin_energy;        //!
  TBranch *b_ec_ecin_sec;           //!
  TBranch *b_ec_ecin_time;          //!
  TBranch *b_ec_ecin_path;          //!
  TBranch *b_ec_ecin_x;             //!
  TBranch *b_ec_ecin_y;             //!
  TBranch *b_ec_ecin_z;             //!
  TBranch *b_ec_ecin_hx;            //!
  TBranch *b_ec_ecin_hy;            //!
  TBranch *b_ec_ecin_hz;            //!
  TBranch *b_ec_ecin_lu;            //!
  TBranch *b_ec_ecin_lv;            //!
  TBranch *b_ec_ecin_lw;            //!
  TBranch *b_ec_ecin_du;            //!
  TBranch *b_ec_ecin_dv;            //!
  TBranch *b_ec_ecin_dw;            //!
  TBranch *b_ec_ecin_m2u;           //!
  TBranch *b_ec_ecin_m2v;           //!
  TBranch *b_ec_ecin_m2w;           //!
  TBranch *b_ec_ecin_m3u;           //!
  TBranch *b_ec_ecin_m3v;           //!
  TBranch *b_ec_ecin_m3w;           //!
  TBranch *b_ec_ecout_energy;       //!
  TBranch *b_ec_ecout_sec;          //!
  TBranch *b_ec_ecout_time;         //!
  TBranch *b_ec_ecout_path;         //!
  TBranch *b_ec_ecout_x;            //!
  TBranch *b_ec_ecout_y;            //!
  TBranch *b_ec_ecout_z;            //!
  TBranch *b_ec_ecout_hx;           //!
  TBranch *b_ec_ecout_hy;           //!
  TBranch *b_ec_ecout_hz;           //!
  TBranch *b_ec_ecout_lu;           //!
  TBranch *b_ec_ecout_lv;           //!
  TBranch *b_ec_ecout_lw;           //!
  TBranch *b_ec_ecout_du;           //!
  TBranch *b_ec_ecout_dv;           //!
  TBranch *b_ec_ecout_dw;           //!
  TBranch *b_ec_ecout_m2u;          //!
  TBranch *b_ec_ecout_m2v;          //!
  TBranch *b_ec_ecout_m2w;          //!
  TBranch *b_ec_ecout_m3u;          //!
  TBranch *b_ec_ecout_m3v;          //!
  TBranch *b_ec_ecout_m3w;          //!
  TBranch *b_cc_nphe_tot;           //!
  TBranch *b_cc_ltcc_sec;           //!
  TBranch *b_cc_ltcc_nphe;          //!
  TBranch *b_cc_ltcc_time;          //!
  TBranch *b_cc_ltcc_path;          //!
  TBranch *b_cc_ltcc_theta;         //!
  TBranch *b_cc_ltcc_phi;           //!
  TBranch *b_cc_ltcc_x;             //!
  TBranch *b_cc_ltcc_y;             //!
  TBranch *b_cc_ltcc_z;             //!
  TBranch *b_cc_htcc_sec;           //!
  TBranch *b_cc_htcc_nphe;          //!
  TBranch *b_cc_htcc_time;          //!
  TBranch *b_cc_htcc_path;          //!
  TBranch *b_cc_htcc_theta;         //!
  TBranch *b_cc_htcc_phi;           //!
  TBranch *b_cc_htcc_x;             //!
  TBranch *b_cc_htcc_y;             //!
  TBranch *b_cc_htcc_z;             //!
  TBranch *b_cc_rich_sec;           //!
  TBranch *b_cc_rich_nphe;          //!
  TBranch *b_cc_rich_time;          //!
  TBranch *b_cc_rich_path;          //!
  TBranch *b_cc_rich_theta;         //!
  TBranch *b_cc_rich_phi;           //!
  TBranch *b_cc_rich_x;             //!
  TBranch *b_cc_rich_y;             //!
  TBranch *b_cc_rich_z;             //!
  TBranch *b_sc_ftof_1a_sec;        //!
  TBranch *b_sc_ftof_1a_time;       //!
  TBranch *b_sc_ftof_1a_path;       //!
  TBranch *b_sc_ftof_1a_energy;     //!
  TBranch *b_sc_ftof_1a_component;  //!
  TBranch *b_sc_ftof_1a_x;          //!
  TBranch *b_sc_ftof_1a_y;          //!
  TBranch *b_sc_ftof_1a_z;          //!
  TBranch *b_sc_ftof_1a_hx;         //!
  TBranch *b_sc_ftof_1a_hy;         //!
  TBranch *b_sc_ftof_1a_hz;         //!
  TBranch *b_sc_ftof_1b_sec;        //!
  TBranch *b_sc_ftof_1b_time;       //!
  TBranch *b_sc_ftof_1b_path;       //!
  TBranch *b_sc_ftof_1b_energy;     //!
  TBranch *b_sc_ftof_1b_component;  //!
  TBranch *b_sc_ftof_1b_x;          //!
  TBranch *b_sc_ftof_1b_y;          //!
  TBranch *b_sc_ftof_1b_z;          //!
  TBranch *b_sc_ftof_1b_hx;         //!
  TBranch *b_sc_ftof_1b_hy;         //!
  TBranch *b_sc_ftof_1b_hz;         //!
  TBranch *b_sc_ftof_2_sec;         //!
  TBranch *b_sc_ftof_2_time;        //!
  TBranch *b_sc_ftof_2_path;        //!
  TBranch *b_sc_ftof_2_energy;      //!
  TBranch *b_sc_ftof_2_component;   //!
  TBranch *b_sc_ftof_2_x;           //!
  TBranch *b_sc_ftof_2_y;           //!
  TBranch *b_sc_ftof_2_z;           //!
  TBranch *b_sc_ftof_2_hx;          //!
  TBranch *b_sc_ftof_2_hy;          //!
  TBranch *b_sc_ftof_2_hz;          //!
  TBranch *b_sc_ctof_time;          //!
  TBranch *b_sc_ctof_path;          //!
  TBranch *b_sc_ctof_energy;        //!
  TBranch *b_sc_ctof_component;     //!
  TBranch *b_sc_ctof_x;             //!
  TBranch *b_sc_ctof_y;             //!
  TBranch *b_sc_ctof_z;             //!
  TBranch *b_sc_ctof_hx;            //!
  TBranch *b_sc_ctof_hy;            //!
  TBranch *b_sc_ctof_hz;            //!
  TBranch *b_sc_cnd_time;           //!
  TBranch *b_sc_cnd_path;           //!
  TBranch *b_sc_cnd_energy;         //!
  TBranch *b_sc_cnd_component;      //!
  TBranch *b_sc_cnd_x;              //!
  TBranch *b_sc_cnd_y;              //!
  TBranch *b_sc_cnd_z;              //!
  TBranch *b_sc_cnd_hx;             //!
  TBranch *b_sc_cnd_hy;             //!
  TBranch *b_sc_cnd_hz;             //!
  TBranch *b_ft_cal_energy;         //!
  TBranch *b_ft_cal_time;           //!
  TBranch *b_ft_cal_path;           //!
  TBranch *b_ft_cal_x;              //!
  TBranch *b_ft_cal_y;              //!
  TBranch *b_ft_cal_z;              //!
  TBranch *b_ft_cal_dx;             //!
  TBranch *b_ft_cal_dy;             //!
  TBranch *b_ft_cal_radius;         //!
  TBranch *b_ft_hodo_energy;        //!
  TBranch *b_ft_hodo_time;          //!
  TBranch *b_ft_hodo_path;          //!
  TBranch *b_ft_hodo_x;             //!
  TBranch *b_ft_hodo_y;             //!
  TBranch *b_ft_hodo_z;             //!
  TBranch *b_ft_hodo_dx;            //!
  TBranch *b_ft_hodo_dy;            //!
  TBranch *b_ft_hodo_radius;        //!

 public:
  Branches12(){};
  Branches12(const std::shared_ptr<TChain> &tree);
  Branches12(const std::shared_ptr<TChain> &tree, bool mc);
  ~Branches12(){};
  int GetEntry(long evnt) { return _tree->GetEntry(evnt); };
  bool mc();
  void mc_branches();
  void init();
  void initMC();
  int gpart();
  int pid(int i);
  float p(int i);
  float p2(int i);
  float px(int i);
  float py(int i);
  float pz(int i);
  float vx(int i);
  float vy(int i);
  float vz(int i);
  int charge(int i);
  float beta(int i);
  float chi2pid(int i);
  int status(int i);

  int run() { return _run; }
  int event() { return _NEVENT; }

  // DC
  int dc_sec(int i);
  float dc_r1_x(int i);
  float dc_r1_y(int i);
  float dc_r1_z(int i);
  float dc_r2_x(int i);
  float dc_r2_y(int i);
  float dc_r2_z(int i);
  float dc_r3_x(int i);
  float dc_r3_y(int i);
  float dc_r3_z(int i);
  // SC
  int sc_ftof_1a_sec(int i);
  float sc_ftof_1a_time(int i);
  float sc_ftof_1a_path(int i);
  float sc_ftof_1a_energy(int i);
  int sc_ftof_1a_component(int i);
  float sc_ftof_1a_x(int i);
  float sc_ftof_1a_y(int i);
  float sc_ftof_1a_z(int i);
  float sc_ftof_1a_hx(int i);
  float sc_ftof_1a_hy(int i);
  float sc_ftof_1a_hz(int i);
  int sc_ftof_1b_sec(int i);
  float sc_ftof_1b_time(int i);
  float sc_ftof_1b_path(int i);
  float sc_ftof_1b_energy(int i);
  int sc_ftof_1b_component(int i);
  float sc_ftof_1b_x(int i);
  float sc_ftof_1b_y(int i);
  float sc_ftof_1b_z(int i);
  float sc_ftof_1b_hx(int i);
  float sc_ftof_1b_hy(int i);
  float sc_ftof_1b_hz(int i);
  int sc_ftof_2_sec(int i);
  float sc_ftof_2_time(int i);
  float sc_ftof_2_path(int i);
  float sc_ftof_2_energy(int i);
  int sc_ftof_2_component(int i);
  float sc_ftof_2_x(int i);
  float sc_ftof_2_y(int i);
  float sc_ftof_2_z(int i);
  float sc_ftof_2_hx(int i);
  float sc_ftof_2_hy(int i);
  float sc_ftof_2_hz(int i);
  float sc_ctof_time(int i);
  float sc_ctof_path(int i);
  float sc_ctof_energy(int i);
  int sc_ctof_component(int i);

  float ec_tot_energy(int i);
  float ec_pcal_energy(int i);
  int ec_pcal_sec(int i);
  float ec_pcal_time(int i);
  float ec_pcal_path(int i);
  float ec_pcal_x(int i);
  float ec_pcal_y(int i);
  float ec_pcal_z(int i);
  float ec_pcal_hx(int i);
  float ec_pcal_hy(int i);
  float ec_pcal_hz(int i);
  float ec_pcal_lu(int i);
  float ec_pcal_lv(int i);
  float ec_pcal_lw(int i);
  float ec_pcal_du(int i);
  float ec_pcal_dv(int i);
  float ec_pcal_dw(int i);
  float ec_pcal_m2u(int i);
  float ec_pcal_m2v(int i);
  float ec_pcal_m2w(int i);
  float ec_pcal_m3u(int i);
  float ec_pcal_m3v(int i);
  float ec_pcal_m3w(int i);

  float ec_ecin_energy(int i);
  int ec_ecin_sec(int i);
  float ec_ecin_time(int i);
  float ec_ecin_path(int i);
  float ec_ecin_x(int i);
  float ec_ecin_y(int i);
  float ec_ecin_z(int i);
  float ec_ecin_hx(int i);
  float ec_ecin_hy(int i);
  float ec_ecin_hz(int i);
  float ec_ecin_lu(int i);
  float ec_ecin_lv(int i);
  float ec_ecin_lw(int i);
  float ec_ecin_du(int i);
  float ec_ecin_dv(int i);
  float ec_ecin_dw(int i);
  float ec_ecin_m2u(int i);
  float ec_ecin_m2v(int i);
  float ec_ecin_m2w(int i);
  float ec_ecin_m3u(int i);
  float ec_ecin_m3v(int i);
  float ec_ecin_m3w(int i);
  float ec_ecout_energy(int i);
  int ec_ecout_sec(int i);
  float ec_ecout_time(int i);
  float ec_ecout_path(int i);
  float ec_ecout_x(int i);
  float ec_ecout_y(int i);
  float ec_ecout_z(int i);
  float ec_ecout_hx(int i);
  float ec_ecout_hy(int i);
  float ec_ecout_hz(int i);
  float ec_ecout_lu(int i);
  float ec_ecout_lv(int i);
  float ec_ecout_lw(int i);
  float ec_ecout_du(int i);
  float ec_ecout_dv(int i);
  float ec_ecout_dw(int i);
  float ec_ecout_m2u(int i);
  float ec_ecout_m2v(int i);
  float ec_ecout_m2w(int i);
  float ec_ecout_m3u(int i);
  float ec_ecout_m3v(int i);
  float ec_ecout_m3w(int i);

  int mc_run();
  int mc_event();
  int mc_type();
  int mc_helicity();
  float mc_weight();
  int mc_npart();
  float mc_ebeam();
  int mc_pid(int i);
  float mc_px(int i);
  float mc_py(int i);
  float mc_pz(int i);
  float mc_vx(int i);
  float mc_vy(int i);
  float mc_vz(int i);
  float mc_vt(int i);
  float cvt_x(int i);
  float cvt_y(int i);
  float cvt_z(int i);
  float fmt_x(int i);
  float fmt_y(int i);
  float fmt_z(int i);

  float cc_nphe_tot(int i);
  int cc_ltcc_sec(int i);
  float cc_ltcc_nphe(int i);
  float cc_ltcc_time(int i);
  float cc_ltcc_path(int i);
  float cc_ltcc_theta(int i);
  float cc_ltcc_phi(int i);
  float cc_ltcc_x(int i);
  float cc_ltcc_y(int i);
  float cc_ltcc_z(int i);
  int cc_htcc_sec(int i);
  float cc_htcc_nphe(int i);
  float cc_htcc_time(int i);
  float cc_htcc_path(int i);
  float cc_htcc_theta(int i);
  float cc_htcc_phi(int i);
  float cc_htcc_x(int i);
  float cc_htcc_y(int i);
  float cc_htcc_z(int i);
  int cc_rich_sec(int i);
  float cc_rich_nphe(int i);
  float cc_rich_time(int i);
  float cc_rich_path(int i);
  float cc_rich_theta(int i);
  float cc_rich_phi(int i);
  float cc_rich_x(int i);
  float cc_rich_y(int i);
  float cc_rich_z(int i);

  float sc_cnd_time(int i);
  float sc_cnd_path(int i);
  float sc_cnd_energy(int i);
  int sc_cnd_component(int i);
  float sc_cnd_x(int i);
  float sc_cnd_y(int i);
  float sc_cnd_z(int i);
  float sc_cnd_hx(int i);
  float sc_cnd_hy(int i);
  float sc_cnd_hz(int i);
  float ft_cal_energy(int i);
  float ft_cal_time(int i);
  float ft_cal_path(int i);
  float ft_cal_x(int i);
  float ft_cal_y(int i);
  float ft_cal_z(int i);
  float ft_cal_dx(int i);
  float ft_cal_dy(int i);
  float ft_cal_radius(int i);
  float ft_hodo_energy(int i);
  float ft_hodo_time(int i);
  float ft_hodo_path(int i);
  float ft_hodo_x(int i);
  float ft_hodo_y(int i);
  float ft_hodo_z(int i);
  float ft_hodo_dx(int i);
  float ft_hodo_dy(int i);
  float ft_hodo_radius(int i);
};

class MCBranches12 : public Branches12 {
 public:
  MCBranches12(const std::shared_ptr<TChain> &tree);
};
#endif
