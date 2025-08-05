#pragma once

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <vector>
#include "G4ParticleDefinition.hh"
#include "unordered_map"

// INITIAL PARTICLES CAN BE: PHOTON, ELECTRON, POSITRON, MUONN, MUONP, NEUTRON, PROTON
// RECORD PARTICLES ARE: PHOTON, ELECTRON, POSITRON, MUONN, MUONP

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
struct geant4_initial_cosmic_ray {
    G4ThreeVector momentum_ini;
    G4ThreeVector position_ini;
    double time;
    double energy;
    G4ParticleDefinition *g4_particle;
};
typedef unsigned int uint;

class Settings {
private:

    Settings() = default; // Private so that it can not be called

    Settings(Settings const &) {}

    // copy constructor is private
    // assignment operator is private
    static Settings *instance;

public:

    static Settings *getInstance();

public:

    const double MAX_POSSIBLE_EVENT_RUNTIME = 30 * 60.0; // 30 minutes

    const bool MAKE_EVERYTHING_VACCUM = false; // if true, sets all density to very low value, to be able to check/debug cosmic ray sampling
    G4String MODE = "run";
    //G4String MODE="visu";

    const bool spread_initial_sampling = true;
    const bool pinhole_initial_sampling = false; // if true, will sample from a point with isotropic direction and record only within the "light beam" projection
    const double max_angle_pinhole_degree = 180.0;
    const double delta_angle_pinhole_degree = 75.0;

    const bool USE_MAX_TIME = true;

    const uint i_ph = 0;
    const uint i_el = 1;
    const uint i_po = 2;
    const uint i_mn = 3;
    const uint i_mp = 4;
    const uint i_ne = 5;
    const uint i_pr = 6;

    const G4String INITIAL_PARTICLE_NAMES[7]{"photon", "electron", "positron", "muonN", "muonP", "neutron", "proton"};

    std::vector<std::string> particules_names_record = {"photon", "electron", "positron"};

    enum PDG_nb : int {
        pdg_phot = 22,
        pdg_elec = 11,
        pdg_posi = -11,
        pdg_muN = 13,
        pdg_muP = -13,
        pdg_neut = 2112,
        pdg_prot = 2212
    };

    const std::vector<int> PDG_LIST_INITIAL = {pdg_phot, pdg_elec, pdg_posi, pdg_neut, pdg_prot};
    // std::vector<int> PDG_LIST_INITIAL = {pdg_phot, pdg_elec, pdg_posi};
     std::vector<int> PDG_LIST_ALL = {pdg_phot, pdg_elec, pdg_posi, pdg_muN, pdg_muP, pdg_neut, pdg_prot};

    std::unordered_map<std::string, int> particle_name_to_index_mapping{
            {INITIAL_PARTICLE_NAMES[0], i_ph},
            {INITIAL_PARTICLE_NAMES[1], i_el},
            {INITIAL_PARTICLE_NAMES[2], i_po},
            {INITIAL_PARTICLE_NAMES[3], i_mn},
            {INITIAL_PARTICLE_NAMES[4], i_mp},
            {INITIAL_PARTICLE_NAMES[5], i_ne},
            {INITIAL_PARTICLE_NAMES[6], i_pr},
    };

    std::unordered_map<int, std::string> particle_index_to_name_mapping{
            {i_ph, INITIAL_PARTICLE_NAMES[0]},
            {i_el, INITIAL_PARTICLE_NAMES[1]},
            {i_po, INITIAL_PARTICLE_NAMES[2]},
            {i_mn, INITIAL_PARTICLE_NAMES[3]},
            {i_mp, INITIAL_PARTICLE_NAMES[4]},
            {i_ne, INITIAL_PARTICLE_NAMES[5]},
            {i_pr, INITIAL_PARTICLE_NAMES[6]},
    };

    std::unordered_map<std::string, int> particle_name_to_PDG_mapping{
            {INITIAL_PARTICLE_NAMES[0], pdg_phot},
            {INITIAL_PARTICLE_NAMES[1], pdg_elec},
            {INITIAL_PARTICLE_NAMES[2], pdg_posi},
            {INITIAL_PARTICLE_NAMES[3], pdg_muN},
            {INITIAL_PARTICLE_NAMES[4], pdg_muP},
            {INITIAL_PARTICLE_NAMES[5], pdg_neut},
            {INITIAL_PARTICLE_NAMES[6], pdg_prot},
    };

    std::unordered_map<int, std::string> particle_PDG_to_name_mapping{
            {pdg_phot, INITIAL_PARTICLE_NAMES[0]},
            {pdg_elec, INITIAL_PARTICLE_NAMES[1]},
            {pdg_posi, INITIAL_PARTICLE_NAMES[2]},
            {pdg_muN,  INITIAL_PARTICLE_NAMES[3]},
            {pdg_muP,  INITIAL_PARTICLE_NAMES[4]},
            {pdg_neut, INITIAL_PARTICLE_NAMES[5]},
            {pdg_prot, INITIAL_PARTICLE_NAMES[6]},
    };

    std::unordered_map<int, int> particle_PDG_to_index_mapping{
            {pdg_phot, i_ph},
            {pdg_elec, i_el},
            {pdg_posi, i_po},
            {pdg_muN,  i_mn},
            {pdg_muP,  i_mp},
            {pdg_neut, i_ne},
            {pdg_prot, i_pr},
    };

    std::unordered_map<int, int> particle_index_to_PDG_mapping{
            {i_ph, pdg_phot},
            {i_el, pdg_elec},
            {i_po, pdg_posi},
            {i_mn, pdg_muN},
            {i_mp, pdg_muP},
            {i_ne, pdg_neut},
            {i_pr, pdg_prot},
    };

    std::unordered_map<int, int> PARMA_ID_to_PDG_mapping{
            {33, pdg_phot},
            {31, pdg_elec},
            {32, pdg_posi},
            {30, pdg_muN},
            {29, pdg_muP},
            {0, pdg_neut},
            {1, pdg_prot},
    };

    int INITIAL_SAMPLE_PDG_TYPE = pdg_phot; // just initialization, will use all values in the main loop

    double CURRENT_WEIGHT = 0.0;

    double POTENTIAL_VALUE = 0.0;   // MV, can be overwritten by input argument

    const double MAX_POSSIBLE_TIME = 1.0 * second;

    double RECORD_ALTITUDE = 10.0; // km
    int RECORD_POSITION = 0; // -1 = below e-field region, 0 = middle of e-field region, 1 = above e-field region
//    const double RECORD_ALTITUDE = 0.141; // km

    // Y axis is direction of altitude (positive towards increasing altitude)
    double EFIELD_REGION_Y_FULL_LENGTH = 1.0;   // km // just initialization, will be changed later
    double EFIELD_REGION_Y_CENTER = 10.0;   // km // just initialization, will be changed later
    double EFIELD_XZ_HALF_SIZE = 40.0; // km

    const double DELTA_ALT_INISAMPLED_RECORD = 1.5; // km

    double CR_SAMPLING_XZ_HALF_SIZE = 40.0; // km
    double CR_GENERATION_ALT_MIN = 50.0; // km, used only by PARMA, just initialization, will be changed later
    double CR_GENERATION_ALT_MAX = 50.1; // km
    double RECORD_XZ_HALF_SIZE = 10.0; // km , could be updated depending on the settings
    const double SOIL_ALT_MAX = 0.0001; // km

    const double drOverR = 0.2;

    enum Efield_State {
        efield_ON, efield_OFF
    };

    std::vector<double> WEIGHTS{0,0,0,0,0,0,0};

    // list of PDG number of particles that we want to generate from parma and be recorded

    const double WORLD_MAX_ALT = 40.0;     // km

    const double GLOBAL_MAX_STEP = 10.0 * meter;
    const bool USE_GLOBAL_MAX_STEP = false;
    const double STEP_MAX_VAL_RECORD_REGION = 5.0 * meter;
    const bool USE_RECORD_REGION_STEP_MAX = false;

    const bool ATMOS_LAYERS_OUTPUT_TO_FILE = false; // write info about atmosphere layers on a file (for test / debug)

    const bool USE_STACKING_ACTION = false; // will mimic time oriented simulation if set to true
    // can be a bad idea to set it on, because it can use a lot of memory compared to default G4 behaviour
    // always turned off if potential is 0

    const double ENERGY_MIN_RECORD = 25.0 * keV;
    const double ENERGY_MAX_RECORD = 1000. * MeV;

    ulong RANDOM_SEED = 250; // just for initialization, will be replace at beginning of main

    Efield_State current_efield_status = efield_ON;
    Efield_State initial_efield_status = efield_ON;

    //    const double longitude = 130.5; // ILDAS positron event coordinates
    //    (Australia)
    //    const double latitude = -13.5;  // ILDAS positron event coordinates
    //    (Australia)
    //    const double longitude = -103.5; // deg, FEGS glow coordinates (Colorado)
    //    const double latitude = 39.5;   // deg, FEGS glow coordinates (Colorado)

    const double longitude = 91.172; // Tibet
    const double latitude = 29.65;   // Tibet

//    const double longitude = 136.7; // Kanazawa
    //   const double latitude = 36.55;   // Kanazawa

    const double CR_GENERATION_ENER_MIN = 0.04;   // MeV
    const double CR_GENERATION_ENER_MAX = 8.95e5 * 0.99; // MeV
    const int year = 2016;
    const int month = 1;
    const int day = 20;

    ulong NB_EVENT = 0; // initialisation is important here
    double VARIABLE_TIME_LIMIT_GLOBAL = 0;

};
