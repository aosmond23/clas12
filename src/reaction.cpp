
#include "reaction.hpp"

Reaction::Reaction(const std::shared_ptr<Branches12> &data, float beam_energy, const std::string &data_type)
    : _data(data), _beam_energy(beam_energy), _data_type(data_type) {
  
  if (data_type == "rec")
      _mc = true;
  else
      _mc = false;
      
  _data = data;
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_energy;
  _sector = data->dc_sec(0);

  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  _gamma = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elec = std::make_unique<TLorentzVector>();
  _elecUnSmear = std::make_unique<TLorentzVector>();
  this->SetElec();

  _prot = std::make_unique<TLorentzVector>();
  _pip = std::make_unique<TLorentzVector>();
  _pim = std::make_unique<TLorentzVector>();
  _other = std::make_unique<TLorentzVector>();
  _neutron = std::make_unique<TLorentzVector>();

}

Reaction::~Reaction() {}

void Reaction::SetElec() {
  _hasE = true;

  _sectorElec = _data->dc_sec(0);
  _elec_status = abs(_data->status(0));

  if (_mc)
  {

          _elecUnSmear->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);

          double W_unsmear = physics::W_calc(*_beam, *_elecUnSmear);

          double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear;

          pUnSmear = _elecUnSmear->P();

          thetaUnSmear = _elecUnSmear->Theta() * 180 / PI;

          if (_elecUnSmear->Phi() > 0)
                  phiUnSmear = _elecUnSmear->Phi() * 180 / PI;
          else if (_elecUnSmear->Phi() < 0)
                  phiUnSmear = (_elecUnSmear->Phi() + 2 * PI) * 180 / PI;

          ////////////////////////////////////////////////////////////////

          // Generate new values
          Reaction::SmearingFunc(ELECTRON, _elec_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

          _pxPrimeSmear = _elecUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                          sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
          _pyPrimeSmear = _elecUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                          sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
          _pzPrimeSmear =
              _elecUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

          // _elecSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_E);  // smeared
          _elec->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_E); // smeared
          // _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);  // unsmeared

          *_gamma += *_beam - *_elec; // be careful you are commenting this only to include the momentum correction

          // // // // // Can calculate W and Q2 here (useful for simulations as sim do not have elec mom corrections)
          _W = physics::W_calc(*_beam, *_elec);
          _Q2 = physics::Q2_calc(*_beam, *_elec);

          _elec_mom = _elec->P();
          _elec_energy = _elec->E();
          _theta_e = _elec->Theta() * 180 / PI;

          // if (_elec->Phi() > 0)
          //         _phi_elec = _elec->Phi() * 180 / PI;
          // else if (_elec->Phi() < 0)
          //         _phi_elec = (_elec->Phi() + 2 * PI) * 180 / PI;
  }
  else
  {
          // fe = objMomCorr->dppC(_data->px(0), _data->py(0), _data->pz(0), _data->dc_sec(0), 0) + 1;
          // _elec->SetXYZM(_data->px(0) * fe, _data->py(0) * fe, _data->pz(0) * fe,
                          // MASS_E); // this is new electron mom corrections aug 2022
          // _elec->SetXYZM(_data->px(0) * fe, _data->py(0) * fe, _data->pz(0) * fe,
          //                MASS_E); // elec and mom corr elec are SAME !!!!!
          _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E); // unsmeared

          *_gamma += *_beam - *_elec;
          _W = physics::W_calc(*_beam, *_elec);
          _Q2 = physics::Q2_calc(*_beam, *_elec);

          _elec_mom = _elec->P();
          _elec_energy = _elec->E();
          _theta_e = _elec->Theta() * 180 / PI;
  }
}

// // before adding electron smearing, this was after _hasE
//   _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);
//   *_gamma += *_beam - *_elec;
//   // // Can calculate W and Q2 here (useful for simulations as sim do not have elec mom corrections)

//   _W = physics::W_calc(*_beam, *_elec);
//   _Q2 = physics::Q2_calc(*_beam, *_elec);

//   _emu_prime_mag2 = _elec->M2();
//   _emu_mag2 = _beam->M2();
//   _elec_energy = _elec->E();

//   _elec_mom = _elec->P();
//   _theta_e = _elec->Theta() * 180.0 / PI;

//   if (_elec->Phi() > 0)
//     _elec_phi = _elec->Phi() * 180 / PI;
//   else if (_elec->Phi() < 0)
//     _elec_phi = (_elec->Phi() + 2 * PI) * 180 / PI;

// }


double Reaction::elec_prime_mass2() {
  if (_emu_prime_mag2 != _emu_prime_mag2) SetElec();

  return _emu_prime_mag2;
}double Reaction::elec_mass2() {
  if (_emu_mag2 != _emu_mag2) SetElec();

  return _emu_mag2;
}

double Reaction::elec_E() {
  if (_elec_energy != _elec_energy) SetElec();
  // std::cout << " emec mom " << _elec_mom << std::endl;

  return _elec_energy;
}
double Reaction::elec_mom() {
  if (_elec_mom != _elec_mom) SetElec();
  // std::cout << " emec mom " << _elec_mom << std::endl;

  return _elec_mom;
}
double Reaction::elec_theta() {
  if (_theta_e != _theta_e) SetElec();

  return _theta_e;
}
double Reaction::elec_phi() {
  if (_elec_phi != _elec_phi) SetElec();

  return _elec_phi;
}

void Reaction::SetProton(int i) {
  _numProt++;
  _numPos++;
  _hasP = true;

  _prot_status = abs(_data->status(i));

  // // auto proton = std::make_unique<TLorentzVector>();
  _prot = std::make_unique<TLorentzVector>();

  _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
  // // proton->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);

  // // _prot.push_back(std::move(proton));

  _prot_indices.clear();
  _prot_indices.push_back(i);

  // _prot = std::make_unique<TLorentzVector>();

  if (_mc)
  {
      // smear ONLY MC reco
      double pUnSmear = _prot->P();
      double thetaUnSmear = _prot->Theta() * 180 / PI;

      double phiUnSmear = (_prot->Phi() > 0)
          ? _prot->Phi() * 180 / PI
          : (_prot->Phi() + 2 * PI) * 180 / PI;

      double pSmear, thetaSmear, phiSmear;

      SmearingFunc(PROTON, _prot_status,
                  pUnSmear, thetaUnSmear, phiUnSmear,
                  pSmear, thetaSmear, phiSmear);

      double px = _prot->Px();
      double py = _prot->Py();
      double pz = _prot->Pz();

      double px_s = px * (pSmear / pUnSmear)
                      * sin(DEG2RAD * thetaSmear)
                      / sin(DEG2RAD * thetaUnSmear)
                      * cos(DEG2RAD * phiSmear)
                      / cos(DEG2RAD * phiUnSmear);

      double py_s = py * (pSmear / pUnSmear)
                      * sin(DEG2RAD * thetaSmear)
                      / sin(DEG2RAD * thetaUnSmear)
                      * sin(DEG2RAD * phiSmear)
                      / sin(DEG2RAD * phiUnSmear);

      double pz_s = pz * (pSmear / pUnSmear)
                      * cos(DEG2RAD * thetaSmear)
                      / cos(DEG2RAD * thetaUnSmear);

      _prot->SetXYZM(px_s, py_s, pz_s, MASS_P);
  }
}

void Reaction::SetPip(int i) {
  _numPip++;
  _numPos++;
  _hasPip = true;

  _pip_status = abs(_data->status(i));

  // auto pip = std::make_unique<TLorentzVector>();
  _pip = std::make_unique<TLorentzVector>();
  
  _pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  // pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);

  // _pip.push_back(std::move(pip));

  _pip_indices.clear();
  _pip_indices.push_back(i);

  if (_mc)
  {
      // smear ONLY MC reco
      double pUnSmear = _pip->P();
      double thetaUnSmear = _pip->Theta() * 180 / PI;

      double phiUnSmear = (_pip->Phi() > 0)
          ? _pip->Phi() * 180 / PI
          : (_pip->Phi() + 2 * PI) * 180 / PI;

      double pSmear, thetaSmear, phiSmear;

      SmearingFunc(PIP, _pip_status,
                  pUnSmear, thetaUnSmear, phiUnSmear,
                  pSmear, thetaSmear, phiSmear);

      double px = _pip->Px();
      double py = _pip->Py();
      double pz = _pip->Pz();

      double px_s = px * (pSmear / pUnSmear)
                      * sin(DEG2RAD * thetaSmear)
                      / sin(DEG2RAD * thetaUnSmear)
                      * cos(DEG2RAD * phiSmear)
                      / cos(DEG2RAD * phiUnSmear);

      double py_s = py * (pSmear / pUnSmear)
                      * sin(DEG2RAD * thetaSmear)
                      / sin(DEG2RAD * thetaUnSmear)
                      * sin(DEG2RAD * phiSmear)
                      / sin(DEG2RAD * phiUnSmear);

      double pz_s = pz * (pSmear / pUnSmear)
                      * cos(DEG2RAD * thetaSmear)
                      / cos(DEG2RAD * thetaUnSmear);

      _pip->SetXYZM(px_s, py_s, pz_s, MASS_PIP);
  }
}

void Reaction::SetPim(int i) {
  _numPim++;
  _numNeg++;
  _hasPim = true;

  _pim_status = abs(_data->status(i));

  // auto pim = std::make_unique<TLorentzVector>();
  _pim = std::make_unique<TLorentzVector>();

  _pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
  // pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);

  // _pim.push_back(std::move(pim));

  _pim_indices.clear();
  _pim_indices.push_back(i);

  if (_mc)
  {
      // smear ONLY MC reco
      double pUnSmear = _pim->P();
      double thetaUnSmear = _pim->Theta() * 180 / PI;

      double phiUnSmear = (_pim->Phi() > 0)
          ? _pim->Phi() * 180 / PI
          : (_pim->Phi() + 2 * PI) * 180 / PI;

      double pSmear, thetaSmear, phiSmear;

      SmearingFunc(PIM, _pim_status,
                  pUnSmear, thetaUnSmear, phiUnSmear,
                  pSmear, thetaSmear, phiSmear);

      double px = _pim->Px();
      double py = _pim->Py();
      double pz = _pim->Pz();

      double px_s = px * (pSmear / pUnSmear)
                      * sin(DEG2RAD * thetaSmear)
                      / sin(DEG2RAD * thetaUnSmear)
                      * cos(DEG2RAD * phiSmear)
                      / cos(DEG2RAD * phiUnSmear);

      double py_s = py * (pSmear / pUnSmear)
                      * sin(DEG2RAD * thetaSmear)
                      / sin(DEG2RAD * thetaUnSmear)
                      * sin(DEG2RAD * phiSmear)
                      / sin(DEG2RAD * phiUnSmear);

      double pz_s = pz * (pSmear / pUnSmear)
                      * cos(DEG2RAD * thetaSmear)
                      / cos(DEG2RAD * thetaUnSmear);

      _pim->SetXYZM(px_s, py_s, pz_s, MASS_PIM);
  }
}

void Reaction::SetNeutron(int i) {
  _numNeutral++;
  _hasNeutron = true;
  _neutron->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_N);
}

void Reaction::SetOther(int i) {
  if (_data->pid(i) == NEUTRON) {
    SetNeutron(i);
  } else {
    _numOther++;
    _hasOther = true;
    _other->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), mass[_data->pid(i)]);
  }
}

float Reaction::w_hadron() {
  if (TwoPion_exclusive())
    return ((*_prot) + (*_pip) + (*_pim)).Mag();
  else
    return NAN;
}
float Reaction::w_difference() {
  if (TwoPion_exclusive())
    return (physics::W_calc(*_beam, *_elec) - ((*_prot) + (*_pip) + (*_pim)).Mag());
  else
    return NAN;
}


void Reaction::CalcMissMass() {
  auto mm_mpim = std::make_unique<TLorentzVector>();
  auto mm_mpip = std::make_unique<TLorentzVector>();
  auto mm_mprot = std::make_unique<TLorentzVector>();
  auto mm_excl = std::make_unique<TLorentzVector>();

  // *mm += (*_gamma + *_target);
  //  _MM = mm->M();

  if (TwoPion_missingPim()) {
    *mm_mpim += (*_gamma + *_target);
    *mm_mpim -= *_prot;
    *mm_mpim -= *_pip;
    _MM2_mPim = mm_mpim->M2();
  
  

    // *mm -= *_pim;

    

  }
  if (TwoPion_exclusive()) {
    *mm_excl += (*_gamma + *_target);
    *mm_excl -= *_prot;
    *mm_excl -= *_pip;
    *mm_excl -= *_pim;


    _MM2_exclusive = mm_excl->M2();
    _excl_Energy = mm_excl->E();
    _excl_Mom = mm_excl->P();
  }

  if (TwoPion_missingPip()) {
  *mm_mpip += (*_gamma + *_target);
  *mm_mpip -= *_prot;
  *mm_mpip -= *_pim;
  _MM2_mPip = mm_mpip->M2();
  }
  if (TwoPion_missingProt()) {
  *mm_mprot += (*_gamma + *_target);
  *mm_mprot -= *_pip;
  *mm_mprot -= *_pim;
  _MM2_mProt = mm_mprot->M2(); //print here
  }
}
float Reaction::MM() {
  if (_MM != _MM) CalcMissMass();
  return _MM;
}
float Reaction::MM2_mPim() {
  if (_MM2_mPim != _MM2_mPim) CalcMissMass();
  return _MM2_mPim;
}
float Reaction::MM2_exclusive() {
  if (_MM2_exclusive != _MM2_exclusive) CalcMissMass();
  return _MM2_exclusive;
}
float Reaction::MM2_mPip() {
  if (_MM2_mPip != _MM2_mPip) CalcMissMass();
  return _MM2_mPip;
}
float Reaction::MM2_mProt() {
  if (_MM2_mProt != _MM2_mProt) CalcMissMass();
  return _MM2_mProt;
}
float Reaction::Energy_excl() {
  if (_excl_Energy != _excl_Energy) CalcMissMass();
  //  std::cout << "_x_mu_p  " << _x_mu->E() << '\n';
  //  if (_x_mu_E > 0)
  return _excl_Energy;
  // else
  // return NAN;
}
float Reaction::Mom_excl() {
  if (_excl_Mom != _excl_Mom) CalcMissMass();
  return _excl_Mom;
  // else
  // return NAN;
}
float Reaction::pim_momentum() {


  if (TwoPion_missingPim()) {
  // if (TwoPion_exclusive()) {
    auto missingpim_ = std::make_unique<TLorentzVector>();
    *missingpim_ += *_gamma + *_target - *_prot - *_pip;

    return missingpim_->P();

  } else
    return NAN;
}
float Reaction::pim_theta_lab() {
  // if (_rec_pim_theta != _rec_pim_theta) CalcMissMass();

  if (TwoPion_missingPim()) {
  // if (TwoPion_exclusive()) {
    auto missingpim_ = std::make_unique<TLorentzVector>();
    *missingpim_ += *_gamma + *_target - *_prot - *_pip;

    return missingpim_->Theta() * 180.0 / PI;
    // return _rec_pim_theta;
  } else
    return NAN;
}
float Reaction::pim_Phi_lab() {
  // if (_rec_pim_phi != _rec_pim_phi) CalcMissMass();

  if (TwoPion_missingPim()) {
  // if (TwoPion_exclusive()) {
    auto missingpim_ = std::make_unique<TLorentzVector>();
    *missingpim_ += *_gamma + *_target - *_prot - *_pip;

    if (missingpim_->Phi() > 0)
      return missingpim_->Phi() * 180 / PI;
    else if (missingpim_->Phi() < 0)
      return (missingpim_->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
    // return _rec_pim_phi;
  } else
    return NAN;
}
float Reaction::pim_momentum_measured() {
  if (TwoPion_exclusive())
    return _pim->P();
  else
    return NAN;
}

float Reaction::pim_theta_lab_measured() {
  if (TwoPion_exclusive())
    return _pim->Theta() * 180.0 / PI;
  else
    return NAN;
}

float Reaction::pim_Phi_lab_measured() {
  if (TwoPion_exclusive()) {
    if (_pim->Phi() > 0)
      return _pim->Phi() * 180 / PI;
    else if (_pim->Phi() < 0)
      return (_pim->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}

float Reaction::pim_theta_angle_btwn_P() {
  if (TwoPion_exclusive()) { // && TwoPion_missingPim()) {
    auto missingpim_ = std::make_unique<TLorentzVector>();
    *missingpim_ += *_gamma + *_target - *_prot - *_pip;

    TVector3 p_miss = missingpim_->Vect();
    TVector3 p_meas = _pim->Vect();

    return p_miss.Angle(p_meas) * 180.0 / PI;
  } else {
    return NAN;
  }
}

////////////////mPip
float Reaction::pip_momentum() {
  if (TwoPion_missingPip()) {
  // if (TwoPion_exclusive()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    *missingpip_ += *_gamma + *_target - *_prot - *_pim;

    return missingpip_->P();
  } else
    return NAN;
}
float Reaction::pip_theta_lab() {
  if (TwoPion_missingPip()) {
  // if (TwoPion_exclusive()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    *missingpip_ += *_gamma + *_target - *_prot - *_pim;
    return missingpip_->Theta() * 180.0 / PI;
  } else
    return NAN;
}
float Reaction::pip_Phi_lab() {
  if (TwoPion_missingPip()) {
  // if (TwoPion_exclusive()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    *missingpip_ += *_gamma + *_target - *_prot - *_pim;

    if (missingpip_->Phi() > 0)
      return missingpip_->Phi() * 180 / PI;
    else if (missingpip_->Phi() < 0)
      return (missingpip_->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float Reaction::pip_momentum_measured() {
  if (TwoPion_exclusive())
    return _pip->P();
  else
    return NAN;
}

float Reaction::pip_theta_lab_measured() {
  if (TwoPion_exclusive())
    return _pip->Theta() * 180.0 / PI;
  else
    return NAN;
}

float Reaction::pip_Phi_lab_measured() {
  if (TwoPion_exclusive()) {
    if (_pip->Phi() > 0)
      return _pip->Phi() * 180 / PI;
    else if (_pip->Phi() < 0)
      return (_pip->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}

float Reaction::pip_theta_angle_btwn_P() {
  if (TwoPion_exclusive()) { // && TwoPion_missingPip()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    *missingpip_ += *_gamma + *_target - *_prot - *_pim;

    TVector3 p_miss = missingpip_->Vect();
    TVector3 p_meas = _pip->Vect();

    return p_miss.Angle(p_meas) * 180.0 / PI;
  } else {
    return NAN;
  }
}


////////////////mProt
float Reaction::prot_momentum() {
  if (TwoPion_missingProt()) {
  // if (TwoPion_exclusive()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    *missingprot_ += *_gamma + *_target - *_pip - *_pim;

    return missingprot_->P();
  } else
    return NAN;
}
float Reaction::prot_theta_lab() {
  if (TwoPion_missingProt()) {
  // if (TwoPion_exclusive()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    *missingprot_ += *_gamma + *_target - *_pip - *_pim;

    return missingprot_->Theta() * 180.0 / PI;
  } else
    return NAN;
}
float Reaction::prot_Phi_lab() {
  if (TwoPion_missingProt()) {
  // if (TwoPion_exclusive()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    *missingprot_ += *_gamma + *_target - *_pip - *_pim;

    if (missingprot_->Phi() > 0)
      return missingprot_->Phi() * 180 / PI;
    else if (missingprot_->Phi() < 0)
      return (missingprot_->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float Reaction::prot_momentum_measured() {
  if (TwoPion_exclusive())
    return _prot->P();
  else
    return NAN;
}

float Reaction::prot_theta_lab_measured() {
  if (TwoPion_exclusive())
    return _prot->Theta() * 180.0 / PI;
  else
    return NAN;
}

float Reaction::prot_Phi_lab_measured() {
  if (TwoPion_exclusive()) {
    if (_prot->Phi() > 0)
      return _prot->Phi() * 180 / PI;
    else if (_prot->Phi() < 0)
      return (_prot->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}

float Reaction::prot_theta_angle_btwn_P() {
  if (TwoPion_exclusive()) { // removed this part to test theta cuts -->    && TwoPion_missingProt()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    *missingprot_ += *_gamma + *_target - *_pip - *_pim;

    TVector3 p_miss = missingprot_->Vect();
    TVector3 p_meas = _prot->Vect();

    return p_miss.Angle(p_meas) * 180.0 / PI;
  } else {
    return NAN;
  }
}

// calculate invariant masses
void Reaction::invMassPpim(const TLorentzVector &prot, const TLorentzVector &pip)
{
  auto delta0 = std::make_unique<TLorentzVector>();
  auto missingpim_ = std::make_unique<TLorentzVector>();
  *missingpim_ += *_gamma + *_target - prot - pip;
  *delta0 += prot + *missingpim_;
  _inv_Ppim = delta0->M();
}
void Reaction::invMasspippim(const TLorentzVector &prot, const TLorentzVector &pip)
{
  auto rho0 = std::make_unique<TLorentzVector>();
  auto missingpim_ = std::make_unique<TLorentzVector>();
  *missingpim_ += *_gamma + *_target - prot - pip;
  *rho0 += pip + *missingpim_;
  _inv_pip_pim = rho0->M();
}
void Reaction::invMassPpip(const TLorentzVector &prot, const TLorentzVector &pip)
{
  auto deltaPP = std::make_unique<TLorentzVector>();
  *deltaPP += prot + pip;
  _inv_Ppip = deltaPP->M();
}

float Reaction::inv_Ppip()
{
  return _inv_Ppip;
}

float Reaction::inv_Ppim()
{
  return _inv_Ppim;
}
float Reaction::inv_pip_pim()
{
  return _inv_pip_pim;
}

MCReaction::MCReaction(const std::shared_ptr<Branches12> &data, float beam_energy, const std::string &data_type)
    : Reaction(data, beam_energy, data_type) {

  _data = data;
  if (!_data->mc()) _data->mc_branches();

  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_energy;
  _weight_mc = _data->mc_weight(); // * 1e4; dont need this anymore!! twopeg update in late dec 2025 seems to work with increased precision!
  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  _gamma_mc = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elec_mc = std::make_unique<TLorentzVector>();
  this->SetMCElec();

  _prot_mc = std::make_unique<TLorentzVector>();
  _pip_mc = std::make_unique<TLorentzVector>();
  _pim_mc = std::make_unique<TLorentzVector>();
  _other_mc = std::make_unique<TLorentzVector>();
  //_neutron = std::make_unique<TLorentzVector>();
}
// Reaction::~Reaction() {} // why this is not here
void MCReaction::SetMCElec() 
{
  _elec_mc->SetXYZM(_data->mc_px(0), _data->mc_py(0), _data->mc_pz(0), MASS_E);
  *_gamma_mc += *_beam - *_elec_mc;

  // Can calculate W and Q2 here
  _W_mc = physics::W_calc(*_beam, *_elec_mc);
  _Q2_mc = physics::Q2_calc(*_beam, *_elec_mc);
}

void MCReaction::SetMCProton(int i) { _prot_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_P); }

void MCReaction::SetMCPip(int i) { _pip_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIP); }

void MCReaction::SetMCPim(int i) { _pim_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIM); }

void MCReaction::SetMCOther(int i) {
  _other_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i),
  mass[_data->pid(i)]);
}

float MCReaction::elec_mom_mc_gen() {
  return _elec_mc->P();
}

float MCReaction::elec_E_mc_gen() { return _elec_mc->E(); }
float MCReaction::elec_theta_mc_gen() { return _elec_mc->Theta() * 180 / PI; }
float MCReaction::elec_phi_mc_gen() {
  if (_elec_mc->Phi() >= 0)
    return (_elec_mc->Phi() * 180 / PI);
  else if (_elec_mc->Phi() < 0)
    return ((_elec_mc->Phi() + 2 * PI) * 180 / PI);
  else
    return NAN;
}

float MCReaction::pim_mom_mc_gen() {
  // if (Reaction::TwoPion_exclusive())
  return _pim_mc->P();
  // else
  //   return NAN;
}
float MCReaction::pip_mom_mc_gen() {
  // if (Reaction::TwoPion_exclusive())
  return _pip_mc->P();
  // else
  //   return NAN;
}
float MCReaction::prot_mom_mc_gen() {
  // if (Reaction::TwoPion_exclusive())
  return _prot_mc->P();
  // else
  //   return NAN;
}
float MCReaction::pim_theta_mc_gen() {
  // if (Reaction::TwoPion_exclusive())
    return _pim_mc->Theta() * 180 / PI;
  // else
  //   return NAN;
}
float MCReaction::pip_theta_mc_gen() {
  // if (Reaction::TwoPion_exclusive())
    return _pip_mc->Theta() * 180 / PI;
  // else
  //   return NAN;
}
float MCReaction::prot_theta_mc_gen() {
  // if (TwoPion_exclusive())
    return _prot_mc->Theta() * 180 / PI;
  // else
  //   return NAN;
}

float MCReaction::pim_phi_mc_gen() {
  if (_pim_mc->Phi() >= 0)
    return (_pim_mc->Phi() * 180 / PI);
  else if (_pim_mc->Phi() < 0)
    return ((_pim_mc->Phi() + 2 * PI) * 180 / PI);
  else
    return NAN;
}
float MCReaction::pip_phi_mc_gen() {
  if (_pip_mc->Phi() >= 0)
    return (_pip_mc->Phi() * 180 / PI);
  else if (_pip_mc->Phi() < 0)
    return ((_pip_mc->Phi() + 2 * PI) * 180 / PI);
  else
    return NAN;
}
float MCReaction::prot_phi_mc_gen() {
  if (_prot_mc->Phi() >= 0)
    return (_prot_mc->Phi() * 180 / PI);
  else if (_prot_mc->Phi() < 0)
    return ((_prot_mc->Phi() + 2 * PI) * 180 / PI);
  else
    return NAN;
}

// calculate invariant mass
float MCReaction::MCinv_Ppip()
{
  auto deltaPPMC = std::make_unique<TLorentzVector>();
  TLorentzVector *prot_mc = _prot_mc.get();
  TLorentzVector *pip_mc = _pip_mc.get();
  *deltaPPMC += *prot_mc + *pip_mc;
  return deltaPPMC->M();
}
float MCReaction::MCinv_Ppim()
{
  auto delta0MC = std::make_unique<TLorentzVector>();
  TLorentzVector *prot_mc = _prot_mc.get();
  TLorentzVector *pim_mc = _pim_mc.get();
  *delta0MC += *prot_mc;
  *delta0MC += *pim_mc;
  return delta0MC->M();
}
float MCReaction::MCinv_pip_pim()
{
  auto rho0MC = std::make_unique<TLorentzVector>();
  TLorentzVector *pip_mc = _pip_mc.get();
  TLorentzVector *pim_mc = _pim_mc.get();
  *rho0MC += *pip_mc + *pim_mc;
  return rho0MC->M();
}