// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#pragma once

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4Box.hh"
#include "G4ParticleTable.hh"
#include <numeric>
#include "AnalysisManager.hh"
#include "fortran.hh"
#include <stdlib.h>
#include <locale.h>

class G4Event;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


struct parameter_set
{
    int type;
    double cos_zenith_angle;
    double altitude;
    double energy;
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction: public G4VUserPrimaryGeneratorAction
{


    public:

        PrimaryGeneratorAction();
        ~PrimaryGeneratorAction();

    public:

        virtual void
        GeneratePrimaries(G4Event *);
        G4ParticleGun *
        GetParticleGun()
        {
            return fParticleGun;
        }

        G4double
        Sample_one_RREA_gammaray_energy(const G4double MinEner,
                                        const G4double MaxEner,
                                        const G4double E_cut);

        G4double
        Sample_one_PL_energy(const G4double MinEner,
                             const G4double MaxEner,
                             const G4double index);

    private:

        G4ParticleGun *fParticleGun; // pointer a to G4 service class
        G4ThreeVector isotropic_direction_rand_sample();


        G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

        G4ThreeVector sample_gaussian_beam_direction(G4ThreeVector position, G4double OpeningAngle);

        G4ThreeVector sample_isotropic_beam_direction(G4ThreeVector position, G4double Opening_Ang);

        G4double sample_CR_particle_energy(G4int particleType);

        G4double interpolate(const std::vector<G4double> &xData, const std::vector<G4double> &yData, G4double x, bool extrapolate);

        G4ThreeVector sample_CR_secondary_position(double altitude);

        G4double sum_vector(std::vector<G4double> vec);

        G4ParticleDefinition *Gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
        G4ParticleDefinition *Electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");
        G4ParticleDefinition *Positron = G4ParticleTable::GetParticleTable()->FindParticle("e+");
        G4ParticleDefinition *Neutron = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
        G4ParticleDefinition *Proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");

        AnalysisManager *analysis = AnalysisManager::getInstance();

        G4double min_cr_ener = 0.02 ; // MeV
        G4double max_cr_ener = 0.95e6; // MeV

        G4double min_cosAng = -1.0;
        G4double max_cosAng = 1.0;

        G4double min_alt = 19.0;
        G4double max_alt = 20.0;
        int nebin = 512; // size of energy mesh (will be log)
        int nabin = 128; // size of angle mesh (linear)
        int naltbin = 1; // size of altitude mesh (linear)
        // the bigger the numbers, more precise will be the sampling, but the more memory it will take
        //        int nebin = 768; // size of energy mesh (will be log)
        //        int nabin = 128; // size of angle mesh (linear)
        //        int naltbin = 256; // size of altitude mesh (linear)

        G4String asciiFileName1;
        std::ofstream asciiFile;

        template<typename T>
        std::vector<double> linspace(T start_in, T end_in, int num_in)
        {

            std::vector<double> linspaced;

            double start = static_cast<double>(start_in);
            double end = static_cast<double>(end_in);
            double num = static_cast<double>(num_in);

            if (num == 0)
                {
                    return linspaced;
                }

            if (num == 1)
                {
                    linspaced.push_back(start);
                    return linspaced;
                }

            double delta = (end - start) / (num - 1);

            for (int i = 0; i < num - 1; ++i)
                {
                    linspaced.push_back(start + delta * i);
                }

            linspaced.push_back(end); // I want to ensure that start and end
            // are exactly the same as the input
            return linspaced;
        }

        std::vector<double> logspace(double min, double max, int nb);

        G4ThreeVector CR_direction_rand_sample(double cos_Sampled);

        parameter_set sample_parameter_set_from_PARMA_distrutions_cumu();

        ///////

        int seed_cr_smpl = 1; // value just for initialisation
        const static int nb_int = 10000000;
        int counter = nb_int + 10; // initialisation, has to be >= nb_int
        double output_energies[nb_int];
        double output_cosangles[nb_int];
        double output_altitudes[nb_int];
        int output_types[nb_int];

        void generator_output_for_test();
        void generator_output_for_test_momentum(const G4ThreeVector &mom);
        G4String name_outFile_mom = "./tests/cr_sampl_test_mom.txt";
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
