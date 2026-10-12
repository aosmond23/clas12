#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD
#include "branches.hpp"
#include "constants.hpp"
#include "deltat.hpp"
#include "physics.hpp"

class Cuts {
protected:
    std::shared_ptr<Branches12> _data = nullptr;
    std::shared_ptr<Delta_T> _dt = nullptr;

public:
    Cuts(const std::shared_ptr<Branches12>& data);
    Cuts(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt);
    ~Cuts();
};

class Pass2_Cuts : public Cuts 
{
    private:
        bool is_gen_data;
        bool is_rec_data;
        bool is_exp_data;
        
    public:
    
    Pass2_Cuts(const std::shared_ptr<Branches12>& data) : Cuts(data) {}
    Pass2_Cuts(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt) : Cuts(data, dt){}

    bool ElectronCuts();
    // bool FiducialCuts();
    bool IsPip(int i);
    bool IsProton(int i);
    bool IsPim(int i);

    // bool CC_nphe_cut();
    bool DC_fiducial_cut_XY(int i, int pid);
    bool EC_sampling_fraction_cut();
    bool PCAL_minimum_energy();
    bool PCAL_fiducial_cut_HX_HY();
    bool DC_z_vertex_cut();
    bool EC_hit_position_fiducial_cut_homogeneous();

    bool CD_fiducial_had(int i);
    bool Hadron_Delta_vz_cut(int i);
    bool Hadron_Chi2pid_cut(int i);

};
#endif
