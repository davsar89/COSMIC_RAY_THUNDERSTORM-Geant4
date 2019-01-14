#include "cosmic_ray_generator_PARMA.hh"

extern "C" {
void gen_parma_cr_(const int *,
                   double[],       // MeV
                   double[],       // cosine of zenith angle (e.g. ang=1.0 for vertical direction, ang=0.0 for holizontal direction)
                   double[],       // km
                   int[],
                   const int *,    // Particle ID (Particle ID, 31:e-, 32:e+, 33:photon)
                   const int *,    // size of energy mesh (will be log)
                   const int *,    // size of angle mesh (linear)
                   const int *,    // size of altitude mesh (linear)
                   const double *, // minimum altitude (km)
                   const double *, // maximum altitude (km)
                   const double *, // emin MeV
                   const double *, // emax MeV
                   const int *,    // iyear
                   const int *,    // imonth
                   const int *,    // iday
                   const double *, // glat deg -90 =< glat =< 90
                   const double *, // glong deg -180 =< glat =< 180
                   const int *,          // wanted Particle ID list (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
                   const int *);   // number of wanted Particle ID list

// See custom_subroutines.f90 for more info

void get_parma_particles_weights_(double[],
                                  int[],
                                  const int *,
                                  const double *,
                                  const double *,
                                  const double *,
                                  const double *,
                                  const int *,    // iyear
                                  const int *,    // imonth
                                  const int *,    // iday
                                  const double *, // glat deg -90 =< glat =< 90
                                  const double *);   // number of wanted Particle ID list
}

Cosmic_Ray_Generator_PARMA::Cosmic_Ray_Generator_PARMA() {
    setlocale(LC_ALL, "C"); // just in case

    if (settings->spread_initial_sampling) {
        naltbin = 30;
    } else {
        naltbin = 2;
    }

    settings->CURRENT_WEIGHT = 1.0;

    nb_part_type_wanted = parmaID_list_ALL.size();

//    std::vector<double> weights = Calculate_Weights_from_PARMA();
//
//    G4cout << "Weights:" << G4endl;
//    for (int ii = 0; ii < weights.size(); ++ii) {
//        G4cout << weights[ii] << G4endl;
//    }
//    G4cout << "  " << G4endl;

    // set up the list of indexes of particles we want

    // first call to PARMA to generate the list of cosmic rays
    Generate_CR_samples_list_from_PARMA();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

geant4_initial_cosmic_ray Cosmic_Ray_Generator_PARMA::generate_Cosmic_ray() {

    geant4_initial_cosmic_ray g4_cosmic = {};

    //            sampled_set = read_particles[index_sampling_part];
    cosmic_ray_parma_output sampled_set = sample_One_CR_from_PARMA();

    G4ThreeVector position_ini{};

    if (settings->pinhole_initial_sampling) {
        position_ini = {0.0,
                        sampled_set.altitude * CLHEP::km,
                        0.0};
    } else {
        position_ini = sample_CR_secondary_position(sampled_set.altitude);
    }

    // from PARMA OUTPUT:
    //   cos(theta)=1,  indicates  the  vertical  downward  direction,
    //   while  90  degree,  i.e.  cos(theta)=0, indicates the horizontal direction.
    G4ThreeVector momentum_ini{};

    if (settings->pinhole_initial_sampling) {
        momentum_ini = direction_rand_sample_isotropic(settings->max_angle_pinhole_degree, settings->delta_angle_pinhole_degree);
    } else {
        momentum_ini = CR_direction_rand_sample(-1.0 * sampled_set.cos_zenith_angle);
    }
    // multiplication by -1 is important, to make sure that when sampled_set.cos_zenith_angle is 1, the particle is sampled vertical downward

    const double time = 0.0;
    const double energy = sampled_set.energy * MeV;

    g4_cosmic.energy = energy;
    g4_cosmic.time = time;
    g4_cosmic.momentum_ini = momentum_ini;
    g4_cosmic.position_ini = position_ini;
    g4_cosmic.g4_particle = sampled_set.g4_part_type;

    return g4_cosmic;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Following angular distributions compuer from PARMA code
G4ThreeVector Cosmic_Ray_Generator_PARMA::CR_direction_rand_sample(const double &cos_Sampled) {

    G4ThreeVector momentum_ini{};

    // if cos_Sampled == 1 (zenith) then direction should be (0,1,0)

    double uu = cos_Sampled;
    double theta = rand_double() * 2.0 * CLHEP::pi;

    momentum_ini.setX(sqrt(1.0 - uu * uu) * cos(theta));
    momentum_ini.setY(uu);
    momentum_ini.setZ(sqrt(1.0 - uu * uu) * sin(theta));

    return momentum_ini;
}

G4ThreeVector Cosmic_Ray_Generator_PARMA::direction_rand_sample_isotropic(const double &max_angle_degree, const double &delta_angle_degree) {

    G4ThreeVector momentum_ini{};

    // if cos_Sampled == 1 (zenith) then direction should be (0,1,0)

    const double min_angle_degree = max_angle_degree - delta_angle_degree;

    double uu = std::cos(min_angle_degree * degree) + (rand_double() * (std::cos(max_angle_degree * degree) - std::cos(min_angle_degree * degree)));
    double theta = rand_double() * 2.0 * CLHEP::pi;

    momentum_ini.setX(sqrt(1.0 - uu * uu) * cos(theta));
    momentum_ini.setY(uu);
    momentum_ini.setZ(sqrt(1.0 - uu * uu) * sin(theta));

    return momentum_ini;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector Cosmic_Ray_Generator_PARMA::sample_CR_secondary_position(const double &altitude) {
    double r1 = rand_double();
    double r2 = rand_double();

    G4ThreeVector position = {(r1 - 0.5) * settings->CR_SAMPLING_XZ_HALF_SIZE * 2.0 * CLHEP::km,
                              altitude * CLHEP::km,
                              (r2 - 0.5) * settings->CR_SAMPLING_XZ_HALF_SIZE * 2.0 * CLHEP::km};

    return position;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// uses cumulative distribution sampling
cosmic_ray_parma_output Cosmic_Ray_Generator_PARMA::sample_One_CR_from_PARMA() {

    if (counter >= nb_int - 3) // if all the particles generated from Parma have been already used, generate new ones
    {
        Generate_CR_samples_list_from_PARMA();

        counter = 0;

        G4cout << "Generated " << nb_int << " new random cosmic ray particles." << G4endl;
    }

    const double eRand = output_energies[counter];
    const double cos_angRand = output_cosangles[counter];
    const double alt_Rand = output_altitudes[counter];
    const int parmaID_sampled = parmaID_list_ALL[output_types[counter]-1];
    G4ParticleDefinition *g4_part_type = G4ParticleTable::GetParticleTable()->FindParticle(settings->PARMA_ID_to_PDG_mapping[parmaID_sampled]);
    counter++;

    if ((alt_Rand < min_alt) || (alt_Rand > max_alt)) {
        G4cout << "Sampled altitude is not in the fixed altitude range. Aborting." << G4endl;
        std::abort();
    }

    if ((cos_angRand < min_cosAng) || (cos_angRand > max_cosAng)) {
        G4cout << "Sampled cosine of angle is not between -1 and 1. Aborting." << G4endl;
        std::abort();
    }

    if ((eRand < min_cr_ener) || (eRand > max_cr_ener)) {
        G4cout << "Energy is out of range. Aborting." << G4endl;
        std::abort();
    }

    if (g4_part_type == nullptr) {
        G4cout << "g4_part_type Pointer is null" << G4endl;
        std::abort();
    }

    cosmic_ray_parma_output spld_set{cos_angRand, alt_Rand, eRand, g4_part_type};

    return spld_set;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Cosmic_Ray_Generator_PARMA::Generate_CR_samples_list_from_PARMA() {

    int nb = nb_int;

    int *parmaID_list_wanted2 = &parmaID_list_ALL[0]; // trick to transform a C++ vector into a C array

    nb_part_type_wanted = parmaID_list_ALL.size();

    const int the_seed = myUtils::generate_a_unique_ID_int32();

    gen_parma_cr_(&the_seed,
                  output_energies,
                  output_cosangles,
                  output_altitudes,
                  output_types,
                  &nb,
                  &nebin,
                  &nabin,
                  &naltbin,
                  &min_alt,
                  &max_alt,
                  &min_cr_ener,
                  &max_cr_ener,
                  &iyear,
                  &imonth,
                  &iday,
                  &glat,
                  &glong,
                  parmaID_list_wanted2,
                  &nb_part_type_wanted);

    counter = 0;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<double> Cosmic_Ray_Generator_PARMA::Calculate_Weights_from_PARMA() {

    nb_part_type_wanted = parmaID_list_ALL.size();

    double output_w[nb_part_type_wanted];

    int *parmaID_list_wanted2 = &parmaID_list_ALL[0]; // trick to transform a C++ vector into a C array

    get_parma_particles_weights_(output_w,
                                 parmaID_list_wanted2,
                                 &nb_part_type_wanted,
                                 &min_alt,
                                 &max_alt,
                                 &min_cr_ener,
                                 &max_cr_ener,
                                 &iyear,
                                 &imonth,
                                 &iday,
                                 &glat,
                                 &glong);

    std::vector<double> output;
    output.clear();
    output.reserve(7);
    for (int ii = 0; ii < nb_part_type_wanted; ++ii) {
        output.push_back(output_w[ii]);
    }

    if (settings->PDG_LIST_INITIAL.size() == 3) {
        output.clear();
        double sum = output_w[0] + output_w[1] + output_w[2];
        output.push_back(output_w[0] / sum);
        output.push_back(output_w[1] / sum);
        output.push_back(output_w[2] / sum);
    }

    return output;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Cosmic_Ray_Generator_PARMA::rand_double() {
    return G4UniformRand();
}