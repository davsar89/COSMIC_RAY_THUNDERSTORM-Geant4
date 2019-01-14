#pragma once

#include "AnalysisManager.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "fortran.hh"
#include "globals.hh"
#include <locale.h>
#include <numeric>
#include <stdlib.h>
#include "Settings.hh"
#include <random>
#include "myUtils.hh"
#include "Randomize.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct cosmic_ray_parma_output {
    double cos_zenith_angle;
    double altitude;
    double energy;
    G4ParticleDefinition * g4_part_type;
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Cosmic_Ray_Generator_PARMA {
public:

    Cosmic_Ray_Generator_PARMA();

    ~Cosmic_Ray_Generator_PARMA() = default;

    geant4_initial_cosmic_ray generate_Cosmic_ray();

    std::vector<double> Calculate_Weights_from_PARMA();

private:

    double rand_double();

    const int i_phot = 0;
    const int i_elec = 1;
    const int i_posi = 2;
    const int i_muN = 3;
    const int i_muP = 4;
    const int i_neut = 5;
    const int i_prot = 6;

    Settings *settings = Settings::getInstance();

    enum Parma_ID : int {
        parma_phot = 33,
        parma_elec = 31,
        parma_posi = 32,
        parma_muN = 30,
        parma_muP = 29,
        parma_neut = 0,
        parma_prot = 1
    };

    std::vector<int> parmaID_list_ALL{parma_phot,
                                      parma_elec,
                                      parma_posi,
                                      parma_neut,
                                      parma_prot};

    ////////////////////////////////////////////////

    int nb_part_type_wanted = 5;

    double min_cr_ener = settings->CR_GENERATION_ENER_MIN;      // MeV
    double max_cr_ener = settings->CR_GENERATION_ENER_MAX;      // MeV
    double min_alt = settings->CR_GENERATION_ALT_MIN;       // km
    double max_alt = settings->CR_GENERATION_ALT_MAX;       // km
    double min_cosAng = -1.0;
    double max_cosAng = 1.0;

    int iyear = settings->year;
    int imonth = settings->month;
    int iday = settings->day;
    double glat = settings->latitude;
    double glong = settings->longitude;

    int nebin = 512; // size of energy mesh (log)
    int nabin = 128; // size of angle mesh (linear)
    int naltbin = 2;   // size of altitude mesh (linear) can be changed latter
    // the bigger the numbers, more precise will be the sampling, but the more memory it will take

    int counter = 0;      // indexing of the sampled cosmic rays

    const static int nb_int = 51000; // number of cosmic to generate. If all are used, it will generate again
    double output_energies[nb_int]{};
    double output_cosangles[nb_int]{};
    double output_altitudes[nb_int]{};
    int output_types[nb_int]{};

    // Functions

    G4ThreeVector CR_direction_rand_sample(const double &cos_Sampled);

    G4ThreeVector sample_CR_secondary_position(const double &altitude);

    cosmic_ray_parma_output sample_One_CR_from_PARMA();

    void Generate_CR_samples_list_from_PARMA();

    G4ThreeVector direction_rand_sample_isotropic(const double &max_angle_degree, const double& delta_angle_degree);
};
