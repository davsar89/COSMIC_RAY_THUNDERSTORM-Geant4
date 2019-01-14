// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Sphere.hh"
#include "G4Material.hh"
#include "G4Box.hh"

#include "Settings.hh"

// !!!!!!!!!!! WARNING : !!!!!!!!!!!!!!!!!!!!!!!
// external fucntion remove because of bug producing NaN values (that is not present when running the fortran code, or in c++ with the debugger)

// interface with fortran PARMA code from http://phits.jaea.go.jp/expacs/
//extern "C" {
//    void getparmaspec_(int *ip, // Particle ID (Particle ID, 31:e-, 32:e+, 33:photon)
//                       double *ener, // MeV
//                       double *cosang, // cosine of zenith angle (e.g. ang=1.0 for vertical direction, ang=0.0 for holizontal direction)
//                       double *Alti, // km
//                       double *glat, // Latitude (deg), -90 =< glat =< 90
//                       double *glong, // Longitude (deg), -180 =< glong =< 180
//                       int *iyear, // tear
//                       int *imonth, // month
//                       int *iday,  // day
//                       double *outputFlux); // Angular Differential Flux(/cm2/s/(MeV/n)/sr)
//    // See custom_subroutines.f90 for more info
//}

// energies,cosangles,altitudes,types
extern "C" {
    void gen_parma_cr_(const int *,
                       double [], // MeV
                       double [], // cosine of zenith angle (e.g. ang=1.0 for vertical direction, ang=0.0 for holizontal direction)
                       double [], // km
                       int [],
                       const int *, // Particle ID (Particle ID, 31:e-, 32:e+, 33:photon)
                       const int *,// size of energy mesh (will be log)
                       const int *,// size of angle mesh (linear)
                       const int *,// size of altitude mesh (linear)
                       const double *,// minimum altitude (km)
                       const double *);// maximum altitude (km)
    // See custom_subroutines.f90 for more info
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(), fParticleGun(0)
{
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);

    // reading input file
    //    G4cout << "Reading input file" << G4endl;
    //    LoadInputDataFile("./input/sampled_particles2.txt");
    //    G4cout << "Input file read" << G4endl;

    setlocale(LC_ALL, "C"); // just in case

    if (Settings::WRITE_MOM_OUTPUT_FOR_TEST)
        {
            std::ofstream asciiFile7;
            asciiFile7.open(name_outFile_mom, std::ios::trunc);
            asciiFile7.close();
        }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
    G4double energy = 0.;
    G4ThreeVector position_ini(0., 0., 0.);
    G4ThreeVector momentum_ini(0., 0., 0.);

    G4ParticleDefinition *particle = 0;

    G4double time = 0.;

    parameter_set sampled_set;
    sampled_set.type = -10; // dummy initialization to remove warning from compiler
    /////////////////////////////////////////////////////////////////////

    if (Settings::initial_state == Settings::parma) // from expacs tabulated values of photon, electron, positron spectra
        {
            //            sampled_set = read_particles[index_sampling_part];
            sampled_set = sample_parameter_set_from_PARMA_distrutions_cumu();

            position_ini = sample_CR_secondary_position(sampled_set.altitude);
            // from PARMA OUTPUT:
            //   cos(theta)=1,  indicates  the  vertical  downward  direction,
            //   while  90  degree,  i.e.  cos(theta)=0, indicates the horizontal direction.
            momentum_ini = CR_direction_rand_sample(-1.0 * sampled_set.cos_zenith_angle);
            // multiplication by -1 is important, to make sure that when sampled_set.cos_zenith_angle is 1, the particle is sampled vertical downward


            //            G4cout << momentum_ini[1] << G4endl;

            if (Settings::WRITE_MOM_OUTPUT_FOR_TEST) generator_output_for_test_momentum(momentum_ini);

            //            G4cout << " Sampled type : " << TYPEE << G4endl;

            if (sampled_set.type == 0)                   // RREA gamma rays from a given distance (approximative spectrum)
                {
                    energy = sampled_set.energy * MeV;
                    particle = Gamma;
                }
            else if (sampled_set.type == -1)
                {
                    energy = sampled_set.energy * MeV;
                    particle = Electron;
                }
            else if (sampled_set.type == 1)
                {
                    energy = sampled_set.energy * MeV;
                    particle = Positron;
                }
            else if (sampled_set.type == 20)
                {
                    energy = sampled_set.energy * MeV;
                    particle = Neutron;
                }
            else if (sampled_set.type == 21)
                {
                    energy = sampled_set.energy * MeV;
                    particle = Proton;
                }
            else
                {
                    G4cout << "ERROR : sampled type is not 0, -1, 1, 20 or 21. Aborting." << G4endl;
                    std::abort();
                }
        }
    else if (Settings::initial_state == Settings::CR_proton) // shooting CR proton
        {
            energy = 30.*GeV;
            position_ini = G4ThreeVector(0, Settings::ALTITUDE_CR_PROTON_SAMPLE, 0);
            momentum_ini = G4ThreeVector(0., -1., 0.);
            particle = Proton;
        }
    else
        {
            G4cout << "Settings::initial_state is not defined. Aborting" << '\n';
            std::abort();
        }
    //    G4cout << energy / keV << G4endl;

    // pos_dir is pair with 3 vector positon and 3 vector momentum direction
    //    position_ini = G4ThreeVector(0, Settings::ALTITUDE_PARAMETER, 0);
    //    momentum_ini = isotropic_direction_rand_sample();
    //    momentum_ini = sample_isotropic_beam_direction(position_ini, OpeningAngle);

    fParticleGun->SetParticleEnergy(energy);
    fParticleGun->SetParticlePosition(position_ini);
    fParticleGun->SetParticleMomentumDirection(momentum_ini);
    fParticleGun->SetParticleTime(time);
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->GeneratePrimaryVertex(anEvent);
} // PrimaryGeneratorAction::GeneratePrimaries

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector PrimaryGeneratorAction::sample_CR_secondary_position(double altitude)
{
    G4ThreeVector position = {(G4UniformRand() - 0.5) *Settings::CR_SAMPLING_XY_HALF_SIZE * 2.0 * CLHEP::km,
                              altitude *CLHEP::km,
                              (G4UniformRand() - 0.5) *Settings::CR_SAMPLING_XY_HALF_SIZE * 2.0 * CLHEP::km
                             };
    return position;
}


// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction::Sample_one_RREA_gammaray_energy(const G4double MinEner, const G4double MaxEner, const G4double cut_ener)
{
    // random samples the energy of one RREA gamma ray
    // RREA gamma spectrum is approximately = 1/E * exp (-E / 7.3MeV)
    // (rejection method)
    const G4double H = cut_ener;
    const G4double pMax = 1. / MinEner * exp(-MinEner / H);
    const G4double pMin = 1. / MaxEner * exp(-MaxEner / H);
    G4double pOfeRand = 0;
    G4double pRand    = 1;
    G4double eRand;

    while (pOfeRand < pRand)   // rejection
        {
            pRand = pMin + (pMax - pMin) * G4UniformRand();
            eRand = MinEner + (MaxEner - MinEner) * G4UniformRand();
            pOfeRand = 1. / eRand * exp(-eRand / H);
        }

    return eRand;

    // G4double energy=MinEner*(pow(((MaxEner)/(MinEner)),G4UniformRand())); // code for sampling 1/E
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction::Sample_one_PL_energy(const G4double MinEner, const G4double MaxEner, const G4double index)
{
    // random samples the energy of one RREA gamma ray
    // RREA gamma spectrum is approximately = 1/E * exp (-E / 7.3MeV)
    // (rejection method)

    const G4double pMax = 1. / std::pow(MinEner, index) ;
    const G4double pMin = 1. / std::pow(MaxEner, index) ;
    G4double pOfeRand = 0;
    G4double pRand    = 1;
    G4double eRand;

    while (pOfeRand < pRand)   // rejection
        {
            pRand = pMin + (pMax - pMin) * G4UniformRand();
            eRand = MinEner + (MaxEner - MinEner) * G4UniformRand();
            pOfeRand = 1. / std::pow(eRand, index) ;
        }

    return eRand;

    // G4double energy=MinEner*(pow(((MaxEner)/(MinEner)),G4UniformRand())); // code for sampling 1/E
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector PrimaryGeneratorAction::sample_gaussian_beam_direction(G4ThreeVector position, G4double Opening_Ang)
{
    G4double R_max          = 10000.;                               // -> maximum angle is atan(10000) = 89.9943 degrees
    G4double sigma_sample_R = std::tan(Opening_Ang);
    G4double R_try          = R_max + 10.; // just for initialization
    G4double X_try = 0., Z_try = 0.;// just for initialization

    while (R_try > R_max)
        {
            X_try = CLHEP::RandGauss::shoot(0., sigma_sample_R); // gaussian position sample
            Z_try = CLHEP::RandGauss::shoot(0., sigma_sample_R); // gaussian position sample

            R_try = sqrt(X_try * X_try + Z_try * Z_try);
        }

    G4ThreeVector position_virtual_particle = position + G4ThreeVector(0., 1., 0.) + G4ThreeVector(1., 0., 0.) * X_try + G4ThreeVector(0., 0., 1.) * Z_try;

    G4ThreeVector momentumDirection = (position_virtual_particle - position);
    momentumDirection = momentumDirection / momentumDirection.mag();

    return momentumDirection;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// isotropic whitin given angle
G4ThreeVector PrimaryGeneratorAction::sample_isotropic_beam_direction(G4ThreeVector position, G4double Opening_Ang)
{
    G4double sigma_sample_R = std::tan(Opening_Ang);
    G4double R_try          = sigma_sample_R * 5.; // just for initialization
    G4double X_try = 0., Z_try = 0.;// just for initialization

    while (R_try > sigma_sample_R)
        {
            X_try = sigma_sample_R * (G4UniformRand() * 2. - 1.); // gaussian position sample
            Z_try = sigma_sample_R * (G4UniformRand() * 2. - 1.); // gaussian position sample

            R_try = sqrt(X_try * X_try + Z_try * Z_try);
        }

    G4ThreeVector position_virtual_particle = position + G4ThreeVector(0., 1., 0.) + G4ThreeVector(1., 0., 0.) * X_try + G4ThreeVector(0., 0., 1.) * Z_try;

    G4ThreeVector momentumDirection = (position_virtual_particle - position);
    momentumDirection = momentumDirection / momentumDirection.mag();

    return momentumDirection;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// isotropic over 4 pi steradians
G4ThreeVector PrimaryGeneratorAction::isotropic_direction_rand_sample()
{
    G4ThreeVector momentum_ini;

    G4double uu    = G4UniformRand();
    G4double theta = G4UniformRand() * 2. * CLHEP::pi;

    momentum_ini.setX(sqrt(1. - uu * uu) * cos(theta));
    momentum_ini.setY(sqrt(1. - uu * uu) * sin(theta));
    momentum_ini.setZ(uu);

    return momentum_ini;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Following angular distributions compuer from PARMA code
G4ThreeVector PrimaryGeneratorAction::CR_direction_rand_sample(double cos_Sampled)
{
    G4ThreeVector momentum_ini;

    // if cos_Sampled == 1 (zenith) then direction should be (0,1,0)

    G4double uu    = cos_Sampled;
    G4double theta = G4UniformRand() * 2. * CLHEP::pi;

    momentum_ini.setX(sqrt(1. - uu * uu) * cos(theta));
    momentum_ini.setY(uu);
    momentum_ini.setZ(sqrt(1. - uu * uu) * sin(theta));

    return momentum_ini;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// uses cumulative distribution sampling

parameter_set PrimaryGeneratorAction::sample_parameter_set_from_PARMA_distrutions_cumu()
{
    G4double eRand = 0.0, cos_angRand = 0.0, alt_Rand = 0.0;
    int idx_rand = 0;

    //    int idx_angRand = 0;
    //    int idx_alt_Rand = 0;
    //energies,cosangles,altitudes,types

    if (counter >= nb_int - 3)
        {
            seed_cr_smpl = (int) std::round(G4UniformRand() * 10000.0);
            int seed_ = seed_cr_smpl;
            int nb = nb_int;
            gen_parma_cr_(&seed_, output_energies, output_cosangles, output_altitudes, output_types, &nb, &nebin, &nabin, &naltbin, &min_alt, &max_alt);
            seed_cr_smpl++;

            counter = 0;

            if (Settings::CR_GENRATOR_write_output_FOR_TEST) generator_output_for_test();

            G4cout << "Generated " << nb_int << " random cosmic ray particles." << G4endl;
        }

    eRand = output_energies[counter];
    cos_angRand = output_cosangles[counter];
    alt_Rand = output_altitudes[counter];
    idx_rand = output_types[counter];
    counter++;

    if (alt_Rand < min_alt || alt_Rand > max_alt)
        {
            G4cout << "Sampled altitude is not in the fixed altitude range. Aborting." << G4endl;
            std::abort();
        }

    if (cos_angRand < min_cosAng || cos_angRand > max_cosAng)
        {
            G4cout << "Sampled cosine of angle is not between -1 and 1. Aborting." << G4endl;
            std::abort();
        }

    // photon electron positron neutron proton
    if (idx_rand < 1 || idx_rand > 5)
        {
            G4cout << "Sampled index is not 1, 2, 3, 4 or 5 . Aborting." << G4endl;
            std::abort();
        }

    if (eRand < min_cr_ener || eRand > max_cr_ener)
        {
            G4cout << "Energy is out of range. Aborting." << G4endl;
            std::abort();
        }

    parameter_set spld_set;
    spld_set.energy = eRand;
    spld_set.cos_zenith_angle = cos_angRand;
    spld_set.altitude = alt_Rand;

    if (idx_rand == 1)
        {
            spld_set.type = 0; // photon
        }

    else if (idx_rand == 2)
        {
            spld_set.type = -1; // electron
        }

    else if (idx_rand == 3)
        {
            spld_set.type = 1; // positron
        }
    else if (idx_rand == 4)
        {
            spld_set.type = 20; // neutron
        }
    else if (idx_rand == 5)
        {
            spld_set.type = 21; // proton
        }
    else
        {
            G4cout << "ERROR in sample_parameter_set_from_PARMA_distrutions_cumu : idx_rand is not 1, 2 or 3. Aborting.";
            std::abort();
        }

    return spld_set;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::generator_output_for_test()
{
    std::ofstream asciiFile6;
    asciiFile6.open("./tests/cr_sampl_test_" + std::to_string(seed_cr_smpl) + ".txt", std::ios::trunc);

    for (int ii = 0; ii < nb_int; ++ii)
        {
            asciiFile6 << output_types[ii] << " " << output_energies[ii] << " " << output_altitudes[ii] << " " << output_cosangles[ii] << G4endl;
        }

    asciiFile6.close();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::generator_output_for_test_momentum(const G4ThreeVector &mom)
{
    std::ofstream asciiFile7;
    asciiFile7.open(name_outFile_mom, std::ios::app);
    asciiFile7 << mom[0] << " " << mom[1] << " " << mom[2] << G4endl;
    asciiFile7.close();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<double> PrimaryGeneratorAction::logspace(double xmin, double xmax, int nbb)
{
    std::vector<double> XX = linspace(log10(xmin), log10(xmax), nbb);

    for (double & X : XX)
        {
            X = std::pow(10.0, X);
        }

    return XX;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction::interpolate(const std::vector<G4double> &xData, const std::vector<G4double> &yData, G4double x, bool extrapolate)
{
    int size = xData.size();

    int i = 0;                                                                  // find left end of interval for interpolation

    if (x >= xData[size - 2])                                                   // special case: beyond right end
        {
            i = size - 2;
        }
    else
        {
            while (x > xData[i + 1]) i++;
        }

    double xL = xData[i], yL = yData[i], xR = xData[i + 1], yR = yData[i + 1];  // points on either side (unless beyond ends)

    if (!extrapolate)                                                           // if beyond ends of array and not extrapolating
        {
            if (x < xL) yR = yL;

            if (x > xR) yL = yR;
        }

    double dydx = (yR - yL) / (xR - xL);                                        // gradient

    return yL + dydx * (x - xL);                                                // linear interpolation
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction::sum_vector(std::vector<G4double> vec)
{
    G4double sum_of_elems = 0.;

    for (auto & n : vec)
        sum_of_elems += n;

    return sum_of_elems;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// exact values from PARMA
//G4double PrimaryGeneratorAction::get_PARMA_density(int idx_rand, double eRand, double cos_angRand, double alt_Rand)
//{

//    int ip = id_list[idx_rand];

//    int iyear = 2016;

//    int imonth = 3;
//    int iday = 1;
//    double glat = Settings::latitude;
//    double glong = Settings::longitude;

//    double outputFlux = 10.0;

//    getparmaspec_(&ip, &eRand, &cos_angRand, &alt_Rand, &glat, &glong, &iyear, &imonth, &iday, &outputFlux);

//    return outputFlux;
//}


//parameter_set PrimaryGeneratorAction::sample_parameter_set_from_PARMA_tabulated_particle_list()
//{
//    parameter_set spld_set;

//    return spld_set;
//}

//void PrimaryGeneratorAction::LoadInputDataFile(const G4String &filename)
//{
//    G4int type_tmp;
//    G4double ene_tmp, cos_ang, alt_tmp;

//    std::ifstream file(filename);

//    parameter_set set_tmp;

//    if (file)
//        {
//            while (!file.eof())
//                {

//                    file >> type_tmp >> ene_tmp >> cos_ang  >> alt_tmp;

//                    if (std::isnan(type_tmp) || std::isnan(ene_tmp) || std::isnan(cos_ang) || std::isnan(alt_tmp)
//                        || std::isinf(type_tmp) || std::isinf(ene_tmp) || std::isinf(cos_ang) || std::isinf(alt_tmp))
//                        {
//                            G4cout << "Error : one of the cosmic ray read parameter is NaN. Aborting" << G4endl;
//                            std::abort();
//                        }

//                    set_tmp.type = type_tmp;
//                    set_tmp.energy = ene_tmp;
//                    set_tmp.cos_zenith_angle = cos_ang;
//                    set_tmp.altitude = alt_tmp;

//                    read_particles.push_back(set_tmp);

//                }
//        }
//    else
//        {
//            G4cout << "ERROR : impossible to read the input file." << G4endl;
//        }

//    max_index_asdouble = (double(read_particles.size()) - 1.0);
//}

//val_index PrimaryGeneratorAction::getGeneration(int nbin, std::vector<double> high, std::vector<double> table)
//{

//    // implicit real * 8(a - h, o - z)
//    // dimension high(0: nbin)
//    // dimension table(0: nbin)

//    val_index v_i;

//    G4double randd = G4UniformRand(); // random number

//    for (int ii = 1; ii <= nbin - 1; ++ii)
//        {

//            if (randd >= table[ii])
//                {
//                    v_i.ind = ii; // bin ID
//                    break;
//                }
//        }

//    randd = G4UniformRand(); // random number
//    v_i.val = high[v_i.ind - 1] * randd + high[v_i.ind] * (1.0 - randd);

//    return v_i;
//}


//// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//// uses simple (unefficient) rejection sampling
//parameter_set PrimaryGeneratorAction::sample_parameter_set_from_PARMA_distrutions()
//{

//    G4double pOfeRand = 0.0;
//    G4double pRand    = 1.0;
//    G4double eRand = 0.0, cos_angRand = 0.0, alt_Rand = 0.0;
//    int idx_rand = 0;
//    //    int idx_angRand = 0;
//    //    int idx_alt_Rand = 0;

//    while ((pOfeRand * eRand < pRand) || std::isnan(pOfeRand))  // rejection sampling with log trick to get more efficient sampling
//        {
//            idx_rand = std::floor(G4UniformRand() * 3.0);

//            // energy is not sampled by index because the grid is not uniromfly spaced
//            eRand = std::log(min_cr_ener) + (std::log(max_cr_ener) - std::log(min_cr_ener)) * G4UniformRand();
//            eRand = exp(eRand);

//            //            idx_angRand = std::floor(G4UniformRand() * double(nb_cosAngle));
//            //            idx_alt_Rand = std::floor(G4UniformRand() * double(nb_altitudes));

//            cos_angRand = min_cosAng + (max_cosAng - min_cosAng) * G4UniformRand();
//            alt_Rand = min_alt + (max_alt - min_alt) * G4UniformRand();

//            pOfeRand = get_PARMA_density(idx_rand, eRand, cos_angRand, alt_Rand);
//            //pOfeRand = get_PARMA_density_approx(idx_rand, eRand, idx_angRand, idx_alt_Rand);

//            pRand = Parma_CR_data_min + (Parma_CR_data_max - Parma_CR_data_min) * G4UniformRand();
//        }

//    //    G4cout << "could sample one !" << G4endl;

//    parameter_set spld_set;

//    spld_set.energy = eRand;
//    spld_set.cos_zenith_angle = cos_angRand;
//    spld_set.altitude = alt_Rand;

//    if (idx_rand == 0)    spld_set.type = 0;

//    if (idx_rand == 1)    spld_set.type = -1;

//    if (idx_rand == 2)    spld_set.type = 1;


//    return spld_set;

//}
