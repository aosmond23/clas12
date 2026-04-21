#include "cuts.hpp"
#include <iostream>
#include "TFile.h"
#include "histogram.hpp"
#include "reaction.hpp"

Cuts::Cuts(const std::shared_ptr<Branches12> &data) : _data(data) { _dt = std::make_shared<Delta_T>(data); }
Cuts::Cuts(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt) : _data(data), _dt(dt) {}
Cuts::~Cuts() {}

// // /////////////////////////// both data dt cuts ////////////// {exp,sim}-> prot, pip, pim /////////////////

double dt_cut_fd_up[2][3][6] = {{{-0.0055, 0.10236, -0.7266, 2.447, -3.926, 2.893},
                                 {-0.003983, 0.0719, -0.4958, 1.634, -2.598, 2.059},
                                 {-0.002695, 0.0391, -0.2146, 0.56, -0.7153, 0.743}},
                                {{-0.003056, 0.05838, -0.429, 1.513, -2.592, 2.297},
                                 {-0.0013075, 0.02531, -0.1893, 0.6865, -1.227, 1.377},
                                 {-0.0009727, 0.0131, -0.0642, 0.1498, -0.2065, 0.6343}}};

double dt_cut_fd_down[2][3][6] = {{{0.005783, -0.1055, 0.7295, -2.377, 3.672, -2.697},
                                   {0.005386, -0.0987, 0.6816, -2.203, 3.326, -2.277},
                                   {0.00419, -0.06027, 0.3281, -0.8403, 1.025, -0.8545}},
                                  {{0.004017, -0.07574, 0.5454, -1.866, 3.047, -2.426},
                                   {0.002022, -0.0358, 0.2452, -0.823, 1.369, -1.379},
                                   {0.002928, -0.04208, 0.2229, -0.5396, 0.6313, -0.7734}}};

double dt_cut_cd_up[2][3][3] = {{{0.06, -0.3303, 0.7656},
                                 {0.01909, 0.00434, 0.428},
                                 {-0.005173, 0.09735, 0.3801}},
                                {{0.05585, -0.2876, 0.714},
                                 {0.0458, -0.1715, 0.5796},
                                 {0.014305, -0.04828, 0.5063}}};

double dt_cut_cd_down[2][3][3] = {{{-0.04974, 0.286, -0.73},
                                   {-0.03662, 0.1555, -0.4675},
                                   {-0.02998, 0.1083, -0.4429}},
                                  {{-0.02893, 0.1836, -0.649},
                                   {-0.0731, 0.303, -0.58745},
                                   {-0.03882, 0.1859, -0.519}}};


// bool Pass2_Cuts::ElectronCuts() {
//   bool _elec = true;
//   // selecting good electrons, this is constant for beam energies, other cuts may vary for beam energy
//   if (!_elec) return false;

//   _elec &= (_data->gpart() > 0);  // Number of good particles is greater than 0
//   _elec &= (_data->gpart() < 20); // but why tho; arbitrary :) ; look at root file for gpart

//   _elec &= (_data->pid(0) == ELECTRON); // setting ==ELECTRON is using their electron 
//   _elec &= (_data->charge(0) == NEGATIVE);

//   _elec &= (2000 <= abs(_data->status(0)) && abs(_data->status(0)) < 4000); // only using forward detector 
//   return _elec;
// }

bool Pass2_Cuts::IsPip(int i) 
{
  int is_mc = 0;
  if (is_gen_data)
  {
    is_mc = 1;
  }

  if (_data->gpart() <= i) return false;
  bool _pip = true;
  _pip &= (_data->charge(i) == POSITIVE);
  _pip &= (_data->pid(i) == PIP);

  // the following is copied from Krishna's code
        
  // _pip &= (abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.4);
  _pip &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 6000);

  // // // // // min/max mom cuts
  if (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000)
  {
    _pip &= (_data->p(i) > 0.5);
    // _pip &= (_data->p(i) < 4.6);
    _pip &= (_dt->dt_Pi(i) < (dt_cut_fd_up[is_mc][1][0] * pow(_data->p(i), 5) + dt_cut_fd_up[is_mc][1][1] * pow(_data->p(i), 4) +
                              dt_cut_fd_up[is_mc][1][2] * pow(_data->p(i), 3) + dt_cut_fd_up[is_mc][1][3] * pow(_data->p(i), 2) +
                              dt_cut_fd_up[is_mc][1][4] * pow(_data->p(i), 1) + dt_cut_fd_up[is_mc][1][5]));
    _pip &= (_dt->dt_Pi(i) > (dt_cut_fd_down[is_mc][1][0] * pow(_data->p(i), 5) + dt_cut_fd_down[is_mc][1][1] * pow(_data->p(i), 4) +
                              dt_cut_fd_down[is_mc][1][2] * pow(_data->p(i), 3) + dt_cut_fd_down[is_mc][1][3] * pow(_data->p(i), 2) +
                              dt_cut_fd_down[is_mc][1][4] * pow(_data->p(i), 1) + dt_cut_fd_down[is_mc][1][5]));

    _pip &= DC_fiducial_cut_XY(i, 2);
  }
  // }
  else if (abs(_data->status(i)) >= 4000)
  {
    _pip &= (_data->p(i) > 0.2);
    // _pip &= (_data->p(i) < 1.7);
    _pip &= (_dt->dt_Pi(i) < (dt_cut_cd_up[is_mc][1][0] * pow(_data->p(i), 2) + dt_cut_cd_up[is_mc][1][1] * _data->p(i) + dt_cut_cd_up[is_mc][1][2]));
    _pip &= (_dt->dt_Pi(i) > (dt_cut_cd_down[is_mc][1][0] * pow(_data->p(i), 2) + dt_cut_cd_down[is_mc][1][1] * _data->p(i) + dt_cut_cd_down[is_mc][1][2]));
    _pip &= CD_fiducial_had(i);
  }
  // _pip &= (_data->p(i) > 0.2);
  _pip &= Hadron_Delta_vz_cut(i);
  _pip &= Hadron_Chi2pid_cut(i);
  return _pip;
}

bool Pass2_Cuts::IsProton(int i) 
{
  int is_mc = 0;
  if (is_gen_data)
  {
    is_mc = 1;
  }
  if (_data->gpart() <= i) return false;
  bool _proton = true;
  _proton &= (_data->charge(i) == POSITIVE);
  _proton &= (_data->pid(i) == PROTON);
  // _proton &= (abs(_dt->dt_P(i)) < 0.5 || abs(_dt->dt_ctof_P(i)) < 0.4);
  // // // _proton &= !(abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.2);
  _proton &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 6000);
  if (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000)
  {
    _proton &= (_data->p(i) > 0.4);
    // _proton &= (_data->p(i) < 4.5);
    _proton &= (_dt->dt_P(i) < (dt_cut_fd_up[is_mc][0][0] * pow(_data->p(i), 5) + dt_cut_fd_up[is_mc][0][1] * pow(_data->p(i), 4) +
                                dt_cut_fd_up[is_mc][0][2] * pow(_data->p(i), 3) + dt_cut_fd_up[is_mc][0][3] * pow(_data->p(i), 2) +
                                dt_cut_fd_up[is_mc][0][4] * pow(_data->p(i), 1) + dt_cut_fd_up[is_mc][0][5]));

    _proton &= (_dt->dt_P(i) > (dt_cut_fd_down[is_mc][0][0] * pow(_data->p(i), 5) + dt_cut_fd_down[is_mc][0][1] * pow(_data->p(i), 4) +
                                dt_cut_fd_down[is_mc][0][2] * pow(_data->p(i), 3) + dt_cut_fd_down[is_mc][0][3] * pow(_data->p(i), 2) +
                                dt_cut_fd_down[is_mc][0][4] * pow(_data->p(i), 1) + dt_cut_fd_down[is_mc][0][5]));

    _proton &= DC_fiducial_cut_XY(i, 1);
  }
  // }
  else if (abs(_data->status(i)) >= 4000)
  // else if ((_dt->isCtof()))
  {
    _proton &= (_data->p(i) > 0.2); /// this 0.4 look harse when we do missing Pim channel
    // _proton &= (_data->p(i) < 2.0);
    _proton &= (_dt->dt_P(i) < (dt_cut_cd_up[is_mc][0][0] * pow(_data->p(i), 2) + dt_cut_cd_up[is_mc][0][1] * _data->p(i) + dt_cut_cd_up[is_mc][0][2]));
    _proton &= (_dt->dt_P(i) > (dt_cut_cd_down[is_mc][0][0] * pow(_data->p(i), 2) + dt_cut_cd_down[is_mc][0][1] * _data->p(i) + dt_cut_cd_down[is_mc][0][2]));
    _proton &= CD_fiducial_had(i);
  }

  // _proton &= (_data->p(i) > 0.2);
  // _proton &= (abs(_data->chi2d(i)) < 0.5);
  _proton &= Hadron_Delta_vz_cut(i);
  _proton &= Hadron_Chi2pid_cut(i);

  return _proton;
}

bool Pass2_Cuts::IsPim(int i) 
{
  int is_mc = 0;
  if (is_gen_data)
  {
    is_mc = 1;
  }
  if (_data->gpart() <= i) return false;
  bool _pim = true;
  _pim &= (_data->charge(i) == NEGATIVE);
  _pim &= (_data->pid(i) == PIM);

  _pim &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 6000);
  if (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000)
  {
    _pim &= (_data->p(i) > 0.4);
    // _pim &= (_data->p(i) < 4.5);
    _pim &= (_dt->dt_Pi(i) < (dt_cut_fd_up[is_mc][2][0] * pow(_data->p(i), 5) + dt_cut_fd_up[is_mc][2][1] * pow(_data->p(i), 4) +
                              dt_cut_fd_up[is_mc][2][2] * pow(_data->p(i), 3) + dt_cut_fd_up[is_mc][2][3] * pow(_data->p(i), 2) +
                              dt_cut_fd_up[is_mc][2][4] * pow(_data->p(i), 1) + dt_cut_fd_up[is_mc][2][5]));

    _pim &= (_dt->dt_Pi(i) > (dt_cut_fd_down[is_mc][2][0] * pow(_data->p(i), 5) + dt_cut_fd_down[is_mc][2][1] * pow(_data->p(i), 4) +
                              dt_cut_fd_down[is_mc][2][2] * pow(_data->p(i), 3) + dt_cut_fd_down[is_mc][2][3] * pow(_data->p(i), 2) +
                              dt_cut_fd_down[is_mc][2][4] * pow(_data->p(i), 1) + dt_cut_fd_down[is_mc][2][5]));
  }
  else if (abs(_data->status(i)) >= 4000)
  {

    _pim &= (_data->p(i) > 0.2);
    // _pim &= (_data->p(i) < 1.9);
    _pim &= (_dt->dt_Pi(i) < (dt_cut_cd_up[is_mc][2][0] * pow(_data->p(i), 2) + dt_cut_cd_up[is_mc][2][1] * _data->p(i) + dt_cut_cd_up[is_mc][2][2]));
    _pim &= (_dt->dt_Pi(i) > (dt_cut_cd_down[is_mc][2][0] * pow(_data->p(i), 2) + dt_cut_cd_down[is_mc][2][1] * _data->p(i) + dt_cut_cd_down[is_mc][2][2]));
  }
  // _pim &= (_data->p(i) > 0.2);

  _pim &= Hadron_Delta_vz_cut(i);
  _pim &= Hadron_Chi2pid_cut(i);

  return _pim;
}

bool Pass2_Cuts::ElectronCuts() {
  bool cut = true;
  if (!cut) return false;

  cut &= (_data->gpart() > 0);
  cut &= (_data->gpart() < 20);
  cut &= (_data->pid(0) == ELECTRON);
  cut &= (_data->charge(0) == NEGATIVE);
  cut &= (2000 <= abs(_data->status(0)) && abs(_data->status(0)) < 4000);

  cut &= (_data->p(0) > 1.50);
  // cut &= CC_nphe_cut();
  cut &= DC_fiducial_cut_XY(0, 0);
  cut &= EC_sampling_fraction_cut();
  cut &= PCAL_minimum_energy();
  cut &= PCAL_fiducial_cut_HX_HY();
  cut &= DC_z_vertex_cut();
  cut &= EC_hit_position_fiducial_cut_homogeneous();

  return cut;
}
// removed bc already in cooking
// bool Pass2_Cuts::CC_nphe_cut()
// {
//   float nphe_min = 2;
//   return (_data->cc_nphe_tot(0) > nphe_min);
// }

bool Pass2_Cuts::DC_fiducial_cut_XY(int i, int pid)
{
        bool _dc_fid_cut = true;
        /// supergaus
        // use minparams_in and maxparams_in for Y limits of a valid hit in given X; gives Y range for specific X
        // min_r for inner circular cut (min radius cut)
        double minparams_in[3][6][3][2] = {
            {{{-0.58933, 15.03000}, {-0.60189, 22.39868}, {-0.54181, 21.96126}},
             {{-0.58242, 14.04545}, {-0.59240, 20.85846}, {-0.55192, 22.29423}},
             {{-0.55927, 14.26636}, {-0.57437, 21.95077}, {-0.51302, 21.37775}},
             {{-0.56533, 14.35000}, {-0.57253, 21.57033}, {-0.52352, 22.35879}},
             {{-0.58582, 16.10909}, {-0.58703, 23.30440}, {-0.52440, 23.26868}},
             {{-0.57624, 13.84455}, {-0.59262, 21.42824}, {-0.53071, 20.85082}}},

            {{{-0.51812, 4.80882}, {-0.60660, 19.32861}, {-0.60392, 34.35024}},
             {{-0.52912, 5.39632}, {-0.61995, 20.57325}, {-0.62008, 37.14595}},
             {{-0.50012, 4.10882}, {-0.59236, 18.29823}, {-0.58771, 32.53500}},
             {{-0.52647, 5.16029}, {-0.60426, 18.88900}, {-0.60434, 34.83024}},
             {{-0.50600, 4.00000}, {-0.58717, 17.29909}, {-0.57177, 28.34119}},
             {{-0.52388, 4.91618}, {-0.60842, 19.32450}, {-0.61696, 36.42214}}},

            {{{-0.48371, 2.71544}, {-0.56539, 14.80074}, {-0.56644, 28.90231}},
             {{-0.49800, 3.45000}, {-0.56931, 14.84191}, {-0.56402, 27.89462}},
             {{-0.49165, 3.19853}, {-0.55265, 13.67279}, {-0.56789, 30.39264}},
             {{-0.48529, 2.72206}, {-0.56480, 14.85221}, {-0.57158, 30.11088}},
             {{-0.46906, 1.67941}, {-0.55137, 13.52941}, {-0.54323, 24.85703}},
             {{-0.51318, 4.11324}, {-0.56931, 14.84191}, {-0.56697, 28.84165}}}

        };

        double maxparams_in[3][6][3][2] = {
            {{{0.57079, -13.46727}, {0.57108, -19.36121}, {0.52978, -20.38791}},
             {{0.56255, -14.79273}, {0.58448, -23.36066}, {0.52522, -23.29478}},
             {{0.62473, -17.85364}, {0.59042, -23.65901}, {0.53396, -24.24066}},
             {{0.58242, -14.04545}, {0.59934, -21.59780}, {0.55005, -22.42995}},
             {{0.57624, -15.64455}, {0.59525, -23.87989}, {0.53313, -24.02225}},
             {{0.59309, -15.06545}, {0.60044, -22.17527}, {0.54753, -22.55247}}},

            {{{0.53865, -5.86103}, {0.61408, -19.00965}, {0.64226, -38.73976}},
             {{0.52876, -5.36985}, {0.61309, -19.36563}, {0.61519, -34.45357}},
             {{0.52388, -4.91618}, {0.60462, -18.49771}, {0.61353, -34.83833}},
             {{0.53606, -5.61691}, {0.61106, -18.60558}, {0.62413, -35.58310}},
             {{0.51347, -4.49779}, {0.60016, -18.02074}, {0.59738, -31.99500}},
             {{0.52918, -5.16324}, {0.60592, -18.18320}, {0.61468, -33.58214}}},

            {{{0.49447, -3.13529}, {0.55402, -12.47426}, {0.57062, -26.22209}},
             {{0.49800, -3.45000}, {0.56431, -13.96176}, {0.57756, -28.62341}},
             {{0.49235, -2.92647}, {0.55402, -12.61544}, {0.56486, -26.74143}},
             {{0.49818, -3.30074}, {0.56431, -13.45588}, {0.56947, -26.46780}},
             {{0.48971, -2.89044}, {0.55402, -12.61544}, {0.55198, -24.45110}},
             {{0.49818, -3.30074}, {0.56431, -13.45588}, {0.56565, -25.70758}}}};

        double min_r[3][6] = {{35.775, 35.765, 35.787, 36.065, 35.887, 36.104},
                              {53.531, 53.314, 52.660, 53.439, 53.473, 53.412},
                              {62.403, 63.254, 63.941, 63.396, 63.273, 62.715}};

        // identify sector where particle was detected
        // short dc_sector = (_data->dc_sec(i) - 1);

        int raw_sector = _data->dc_sec(i);

        if (raw_sector <= 0 || raw_sector > 6) {
            return false;
        }

        short dc_sector = raw_sector - 1;

        // region 1
        double X1 = _data->dc_r1_x(i);
        double Y1 = _data->dc_r1_y(i);
        float X1_new = X1 * cos(DEG2RAD * (-60 * (dc_sector))) - Y1 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y1 = X1 * sin(DEG2RAD * (-60 * (dc_sector))) + Y1 * cos(DEG2RAD * (-60 * (dc_sector)));

        X1 = X1_new;
        int region_1 = 1;

        double calc_min1 = minparams_in[pid][dc_sector][region_1 - 1][0] * X1 + minparams_in[pid][dc_sector][region_1 - 1][1];
        double calc_max1 = maxparams_in[pid][dc_sector][region_1 - 1][0] * X1 + maxparams_in[pid][dc_sector][region_1 - 1][1];
        double DC_r1 = sqrt(X1 * X1 + Y1 * Y1);
        double Min_r1 = min_r[region_1 - 1][dc_sector];
        // region 2
        double X2 = _data->dc_r2_x(i);
        double Y2 = _data->dc_r2_y(i);
        float X2_new = X2 * cos(DEG2RAD * (-60 * (dc_sector))) - Y2 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y2 = X2 * sin(DEG2RAD * (-60 * (dc_sector))) + Y2 * cos(DEG2RAD * (-60 * (dc_sector)));
        X2 = X2_new;
        int region_2 = 2;
        double calc_min2 = minparams_in[pid][dc_sector][region_2 - 1][0] * X2 + minparams_in[pid][dc_sector][region_2 - 1][1];
        double calc_max2 = maxparams_in[pid][dc_sector][region_2 - 1][0] * X2 + maxparams_in[pid][dc_sector][region_2 - 1][1];
        double DC_r2 = sqrt(X2 * X2 + Y2 * Y2);
        double Min_r2 = min_r[region_2 - 1][dc_sector];
        // region 3
        double X3 = _data->dc_r3_x(i);
        double Y3 = _data->dc_r3_y(i);
        float X3_new = X3 * cos(DEG2RAD * (-60 * (dc_sector))) - Y3 * sin(DEG2RAD * (-60 * (dc_sector)));
        Y3 = X3 * sin(DEG2RAD * (-60 * (dc_sector))) + Y3 * cos(DEG2RAD * (-60 * (dc_sector)));
        X3 = X3_new;
        int region_3 = 3;
        double calc_min3 = minparams_in[pid][dc_sector][region_3 - 1][0] * X3 + minparams_in[pid][dc_sector][region_3 - 1][1];
        double calc_max3 = maxparams_in[pid][dc_sector][region_3 - 1][0] * X3 + maxparams_in[pid][dc_sector][region_3 - 1][1];
        double DC_r3 = sqrt(X3 * X3 + Y3 * Y3);
        double Min_r3 = min_r[region_3 - 1][dc_sector];

        // particle must be within all fiducial bands in all 3 regions and pass inner circular cut
        return ((Y1 > calc_min1) && (Y1 < calc_max1) && (Y2 > calc_min2) && (Y2 < calc_max2) && (Y3 > calc_min3) &&
                (Y3 < calc_max3) && (DC_r1 > Min_r1) && (DC_r2 > Min_r2) && (DC_r3 > Min_r3));
}

bool Pass2_Cuts::EC_sampling_fraction_cut()
{
  int isec = (_data->ec_pcal_sec(0) - 1);
  double upper_lim_total = 0;
  double lower_lim_total = 0;
  int is_mc = 0; // manually change this to 1 for sim, 0 for exp until below part works
  if (is_gen_data)   // having issues with this part of the code but it works otherwise (i think); update: fixed!! 5/14/25
  {
    is_mc = 1;
  }

  ///// Both  exp data and sim data ////////
  double mean_minus_3_5_sigma[2][6][3] = {{{-0.0001186, 0.0001892, 0.1942}, {-0.000856, 0.01084, 0.1637}, {-0.001184, 0.014046, 0.1593}, {-0.001268, 0.01918, 0.1287}, {-0.0002744, 0.003532, 0.1844}, {-0.001039, 0.012505, 0.1593}},
                                          {{-0.00058, 0.00687, 0.19312}, {-0.00088, 0.01022, 0.18360}, {-0.00089, 0.00941, 0.18832}, {-0.00066, 0.00888, 0.18466}, {-0.00066, 0.00798, 0.18884}, {-0.00055, 0.00685, 0.19319}}};
  double mean_plus_3_5_sigma[2][6][3] = {{{-0.0004027, 0.001746, 0.2903}, {-9.36e-05, -0.000999, 0.2979}, {-0.0003238, 0.00101, 0.2957}, {-4.303e-05, -0.0004702, 0.2954}, {-0.0001818, 0.003223, 0.2742}, {-0.0002906, 0.0015335, 0.2883}},
                                         {{-0.00002, -0.00078, 0.29991}, {0.00023, -0.00396, 0.31026}, {0.00010, -0.00156, 0.30077}, {0.00017, -0.00400, 0.31052}, {0.00018, -0.00342, 0.30823}, {0.00012, -0.00297, 0.30706}}};

  for (Int_t k = 0; k < 6; k++)
  {
    if (isec == k)
    {
      upper_lim_total = mean_plus_3_5_sigma[is_mc][k][0] * pow(_data->p(0), 2) + (mean_plus_3_5_sigma[is_mc][k][1]) * _data->p(0) +
                        mean_plus_3_5_sigma[is_mc][k][2];

      lower_lim_total = mean_minus_3_5_sigma[is_mc][k][0] * pow(_data->p(0), 2) + (mean_minus_3_5_sigma[is_mc][k][1]) * _data->p(0) +
                        mean_minus_3_5_sigma[is_mc][k][2];
    }
  }
  bool pass_band = _data->ec_tot_energy(0) / _data->p(0) <= upper_lim_total &&
                   _data->ec_tot_energy(0) / _data->p(0) >= lower_lim_total;
  // bool pass_band = true;
  // bool pass_triangle = false;

  // if (_data->p(0) < 4.5)
  // {
  //   pass_triangle = true;
  // }
  // else
  // {
  //   // pass_triangle = (_data->ec_ecin_energy(0) / _data->p(0)) > (0.2 - _data->ec_pcal_energy(0) / _data->p(0));
  //   pass_triangle = true;
  // }

  if (pass_band) // && pass_triangle)
    return true;
  else
    return false;
}

// error with cc_nphe_tot; changed now to ec_pcal_energy
bool Pass2_Cuts::PCAL_minimum_energy()
{
  double edep_tight = 0.06, edep_medium = 0.07, edep_loose = 0.09;
  return (_data->ec_pcal_energy(0) > edep_medium);
}

bool Pass2_Cuts::PCAL_fiducial_cut_HX_HY()
{
  double minparams_pcal_in[6][2] = {{-0.52452, 20.33242}, {-0.51548, 18.38758}, {-0.49609, 19.04455}, {-0.51318, 22.13909}, {-0.50361, 20.48697}, {-0.51821, 19.48394}};
  double maxparams_pcal_in[6][2] = {{0.52494, -20.38030}, {0.50706, -22.01970}, {0.50900, -21.77000}, {0.51967, -19.31667}, {0.52082, -23.08091}, {0.52288, -20.65061}};

  double min_radius[6] = {71.245, 71.587, 72.142, 73.101, 72.025, 72.921};

  // short pcal_sector = (_data->ec_pcal_sec(0) - 1);

  int raw_sector = _data->ec_pcal_sec(0);

  if (raw_sector <= 0 || raw_sector > 6) {
      return false;
  }

  short pcal_sector = raw_sector - 1;

  double HX = _data->ec_pcal_hx(0);
  double HY = _data->ec_pcal_hy(0);

  float HX_new = HX * cos(DEG2RAD * (-60 * (pcal_sector))) - HY * sin(DEG2RAD * (-60 * (pcal_sector)));
  HY = HX * sin(DEG2RAD * (-60 * (pcal_sector))) + HY * cos(DEG2RAD * (-60 * (pcal_sector)));

  HX = HX_new;

  double H_radius = sqrt(HX * HX + HY * HY);
  double Min_radius = min_radius[pcal_sector];

  double calc_min = minparams_pcal_in[pcal_sector][0] * HX + minparams_pcal_in[pcal_sector][1];
  double calc_max = maxparams_pcal_in[pcal_sector][0] * HX + maxparams_pcal_in[pcal_sector][1];

  return ((HY > calc_min) && (HY < calc_max) && (H_radius > Min_radius));
}

bool Pass2_Cuts::DC_z_vertex_cut(){
  int pcal_sector = _data->ec_pcal_sec(0);
  float partvz = _data->vz(0);

  float vz_min_sect_inb[] = {-10, -10, -10, -10, -10, -10};
  float vz_max_sect_inb[] = {5, 5, 5, 5, 5, 5};

  float vz_min_sect[6];
  float vz_max_sect[6];

  for (int i = 0; i < 6; i++)
  {
    vz_min_sect[i] = vz_min_sect_inb[i];
    vz_max_sect[i] = vz_max_sect_inb[i];
  }

  int isec = pcal_sector - 1;
  float vz_min = vz_min_sect[isec];
  float vz_max = vz_max_sect[isec];

  return partvz > vz_min && partvz < vz_max;
}

bool Pass2_Cuts::EC_hit_position_fiducial_cut_homogeneous()
{
        // Cut using the natural directions of the scintillator bars/ fibers:
        double min_v_tight_inb[6] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
        double min_v_med_inb[6] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
        double min_v_loose_inb[6] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
        
        double max_v_tight_inb[6] = {400, 400, 400, 400, 400, 400};
        double max_v_med_inb[6] = {400, 400, 400, 400, 400, 400};
        double max_v_loose_inb[6] = {400, 400, 400, 400, 400, 400};
        
        double min_w_tight_inb[6] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
        double min_w_med_inb[6] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
        double min_w_loose_inb[6] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
        
        double max_w_tight_inb[6] = {400, 400, 400, 400, 400, 400};
        double max_w_med_inb[6] = {400, 400, 400, 400, 400, 400};
        double max_w_loose_inb[6] = {400, 400, 400, 400, 400, 400};
        
        int isec = (_data->ec_pcal_sec(0) - 1);
        double min_v = min_v_med_inb[isec];
        double max_v = max_v_med_inb[isec];
        double min_w = min_w_med_inb[isec];
        double max_w = max_w_med_inb[isec];
        return (_data->ec_pcal_lv(0) > min_v && _data->ec_pcal_lv(0) < max_v && _data->ec_pcal_lw(0) > min_w &&
                _data->ec_pcal_lw(0) < max_w);
}

bool Pass2_Cuts::CD_fiducial_had(int i)
{
  bool pass_fiducial = true;
  // int pid = _data->pid(i);
  //        if (pid == 2212)
  {
    double momT = sqrt(_data->px(i) * _data->px(i) + _data->py(i) * _data->py(i));
    double theta = atan2(momT, _data->pz(i)) * 180 / PI;
    double phi = atan2(_data->py(i), _data->px(i)) * 180 / PI;

    double fiducial_phi_width = 3; // 3 is used by andrew
    double fiducial_phi_shift = 0;
    double fiducial_momT_start = 0.15;
    double fiducial_phi_central = (-asin(fiducial_momT_start / momT) - (PI / 2)) * 180 / PI;

    if ((fabs(phi - fiducial_phi_central - fiducial_phi_shift) < fiducial_phi_width) ||
        (fabs(phi - fiducial_phi_central - fiducial_phi_shift - 120) < fiducial_phi_width) ||
        (fabs(phi - fiducial_phi_central - fiducial_phi_shift - 240) < fiducial_phi_width)) // || (theta < 40) || (theta > 125))
    {
      pass_fiducial = false;
    }
  }
  return pass_fiducial;
}

// /** Delta VZ cut for hadrons
//  * @param pid hadron PID code
//  * @param dvz difference between Vz of hadron candidate and electron
//  */
bool Pass2_Cuts::Hadron_Delta_vz_cut(int i)
{
  int pid = _data->pid(i);
  float dvz = (_data->vz(i) - _data->vz(0));
  return dvz > -20 && dvz < 20;
}

// chi2pid cut for hadrons

bool Pass2_Cuts::Hadron_Chi2pid_cut(int i)
{
  bool isstrict = false;
  float chi2pid = _data->chi2pid(i);
  float p = _data->p(i);
  int pid = _data->pid(i);
  int status = abs(_data->status(i));

  if (status < 4000)
    return abs(chi2pid) < 5.0; /// trying very loose cuts
  else
  {
    return abs(chi2pid) < 7.0;
  }
}