// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************

#include "AnalysisManager.hh"
#include <fstream>
#include <iomanip>
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "Settings.hh"

// class following singleton pattern

AnalysisManager *AnalysisManager::instance = 0;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AnalysisManager *AnalysisManager::getInstance()  // singleton lazy initialization
{
    if (instance == 0) instance = new AnalysisManager;

    return instance;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
AnalysisManager::AnalysisManager() // constructor
{

    RECORDED_OUTPUT_STRINGS.clear();

    G4String first_part = "./output/Part_ener_mom_dists_" ;
    G4String potential_str = "_" + std::to_string(int(Settings::POTENTIAL_VALUE));
    G4String efield_alt_str = "_" + std::to_string(int(Settings::EFIELD_REGION_ALT_CENTER * 10.));
    G4String efield_length_str = "_" + std::to_string(int(Settings::EFIELD_REGION_LEN * 100.));
    G4String rec_alt_part = "_" + std::to_string(int(Settings::RECORD_ALTITUDES[0] * 10.)) + "km";

    asciiFileName1 = first_part + std::to_string(Settings::RANDOM_SEED) + rec_alt_part + potential_str + efield_alt_str + efield_length_str + ".out";

    G4cout << "Creating output file : " << asciiFileName1 << G4endl;

    // open / close with trunc to clean the file (if it does exist) or to create it

    asciiFile_analysis.open(asciiFileName1, std::ios::trunc);

    if (asciiFile_analysis.is_open())
        {
            asciiFile_analysis.close();
        }
    else
        {
            G4cout << G4endl << "ERROR : cannot open output file. Aborting" << G4endl;
            std::abort();
        }

    ENER_GRID = get_ener_grid();
    MOM_GRID = get_MOM_grid();

    for (uint jj = 0; jj < 3; ++jj) // loop over photon electrons and positrons
        {
            for (uint ii = 0; ii < ENER_GRID.size(); ++ii)
                {
                    SPEC_PART[jj].push_back(0);
                }

            for (uint ii = 0; ii < MOM_GRID.size(); ++ii)
                {
                    PART_MOM_X[jj].push_back(0);
                    PART_MOM_Y[jj].push_back(0);
                    PART_MOM_Z[jj].push_back(0);
                }
        }

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AnalysisManager::~AnalysisManager()
{

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int AnalysisManager::NB_OUTPUT() const
{
    return NB_OUTPUT_;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AnalysisManager::add_NB_OUTPUT()
{
    NB_OUTPUT_++;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void AnalysisManager::write_output_file_endOf_program()
{
    asciiFile_analysis.open(asciiFileName1, std::ios::out | std::ios::app);

    if (asciiFile_analysis.is_open())
        {
            asciiFile_analysis  << std::scientific << std::setprecision(6)
                                << Settings::RANDOM_SEED
                                << " " << Settings::NB_EVENT
                                << " " << Settings::EFIELD_REGION_ALT_CENTER
                                << " " << Settings::EFIELD_REGION_LEN
                                << " " << Settings::POTENTIAL_VALUE
                                << " " << Settings::TILT
                                << " " << Settings::RREA_PART_NB_LIMIT_HAS_BEEN_REACHED
                                << " " << photon_counter_up
                                << " " << electron_counter_up
                                << " " << positron_counter_up
                                << " " << photon_counter_down
                                << " " << electron_counter_down
                                << " " << positron_counter_down;

            for (uint jj = 0; jj < 3; ++jj)
                {
                    for (uint ii = 0; ii < ENER_GRID.size(); ++ii)
                        {
                            asciiFile_analysis << " " << SPEC_PART[jj][ii];
                        }

                    for (uint ii = 0; ii < MOM_GRID.size(); ++ii)
                        {
                            asciiFile_analysis << " " << PART_MOM_X[jj][ii];
                        }

                    for (uint ii = 0; ii < MOM_GRID.size(); ++ii)
                        {
                            asciiFile_analysis << " " << PART_MOM_Y[jj][ii];
                        }

                    for (uint ii = 0; ii < MOM_GRID.size(); ++ii)
                        {
                            asciiFile_analysis << " " << PART_MOM_Z[jj][ii];
                        }
                }

            asciiFile_analysis << G4endl;
        }
    else
        {
            G4cout << G4endl << "ERROR : cannot open output file. Aborting." << G4endl;
            std::abort();
        }

    asciiFile_analysis.close();
}


//======================================================================
// Returns interpolated value at x from parallel arrays ( xData, yData )
//   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
//   boolean argument extrapolate determines behaviour beyond ends of array (if needed)
G4double AnalysisManager::interpolate(std::vector<G4double> &xData, std::vector<G4double> &yData, G4double x, bool extrapolate)
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


G4double AnalysisManager::get_scale(G4double alt)
{
    // returns the atmospheric relative scale compared to sea level (>1)
    std::vector<G4double> alt_list =
    {
        0,
        0.500000000000000,
        1.000000000000000,
        1.500000000000000,
        2.000000000000000,
        2.500000000000000,
        3.000000000000000,
        3.500000000000000,
        4.000000000000000,
        4.500000000000000,
        5.000000000000000,
        5.500000000000000,
        6.000000000000000,
        6.500000000000000,
        7.000000000000000,
        7.500000000000000,
        8.000000000000000,
        8.500000000000000,
        9.000000000000000,
        9.500000000000000,
        10.000000000000000,
        10.500000000000000,
        11.000000000000000,
        11.500000000000000,
        12.000000000000000,
        12.500000000000000,
        13.000000000000000,
        13.500000000000000,
        14.000000000000000,
        14.500000000000000,
        15.000000000000000,
        15.500000000000000,
        16.000000000000000,
        16.500000000000000,
        17.000000000000000,
        17.500000000000000,
        18.000000000000000,
        18.500000000000000,
        19.000000000000000,
        19.500000000000000,
        20.000000000000000,
    };

    std::vector<G4double> scale_list =
    {
        1.000000000000000,
        1.059301380991064,
        1.121238177128117,
        1.184377838328792,
        1.249042145593870,
        1.317304778260430,
        1.388415672913118,
        1.463688404983724,
        1.543377914546100,
        1.628371628371629,
        1.719862833025587,
        1.818435364663227,
        1.925291598996014,
        2.041647095663066,
        2.168634624979211,
        2.307556184746062,
        2.460377358490566,
        2.628502318081032,
        2.813376483279396,
        3.018518518518519,
        3.245395719263315,
        3.496916063287745,
        3.774240231548481,
        4.080100125156446,
        4.414353419092755,
        4.781811514484781,
        5.180770758839889,
        5.613430908308223,
        6.082089552238807,
        6.589186457806973,
        7.133479212253830,
        7.720544701006513,
        8.348271446862997,
        9.024221453287199,
        9.753178758414361,
        10.541632983023444,
        11.398601398601400,
        12.325141776937620,
        13.334696799263730,
        14.431164231961047,
        15.627996164908916,
    };

    if ((alt > 20.) || (alt < 0.))
        {
            G4cout << "ERROR in get_scale : altitude is not between 0 and 20. Aborting.";
            std::abort();
        }

    return interpolate(alt_list, scale_list, alt, false);

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::fill_histogram(G4int idx_part, const std::vector<double> &GRID, std::array<std::vector<uint>, 3> &COUNTS, const G4double value)
{
    for (uint ii = 0; ii < GRID.size() - 1; ++ii)
        {
            if (value >= GRID[ii] && value < GRID[ii + 1])
                {
                    COUNTS[idx_part][ii]++;
                    return;
                }
        }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::fill_histogram_ener(G4int idx_part, const G4double value)
{
    fill_histogram(idx_part, ENER_GRID, SPEC_PART, value);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::fill_histogram_momX(G4int idx_part, const G4double value)
{
    fill_histogram(idx_part, MOM_GRID, PART_MOM_X, value);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::fill_histogram_momY(G4int idx_part, const G4double value)
{
    fill_histogram(idx_part, MOM_GRID, PART_MOM_Y, value);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::fill_histogram_momZ(G4int idx_part, const G4double value)
{
    fill_histogram(idx_part, MOM_GRID, PART_MOM_Z, value);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::check_if_should_use_stacking_action()
{
    // checks if E-field is above RREA, and if so, sets Settings::USE_STACKING_ACTION to true

    G4double EFIELD_VAL = (Settings::POTENTIAL_VALUE * megavolt) / (Settings::EFIELD_REGION_LEN * km) ;

    G4double margin_factor = 1.0;
    G4double EFIELD_VAL_KV_PER_M = std::abs(EFIELD_VAL / (kilovolt / meter)) * margin_factor;

    G4double atmos_scale = get_scale(Settings::EFIELD_REGION_ALT_CENTER); // should be > 1

    Settings::DELTA_T = 27.3 / (EFIELD_VAL_KV_PER_M * atmos_scale - 277.0) * atmos_scale; // in micro second for the stacking action that set up the fake time oriented

    Settings::DELTA_T = Settings::DELTA_T / 10.0 * microsecond;

    double DELTA_T_tmp = Settings::DELTA_T ; // variable for debug mode check

    if (atmos_scale < 1.0)
        {
            G4cout << "ERROR in create_thunderstorm_electric_fields : atmos_scale is <1. Aborting" << G4endl;
            std::abort();
        }

    // if below RREA threshold, no need to have the stacking action for fake time oriented
    // and the calculed time step (DELTA_T) should be negative
    if (Settings::DELTA_T < 0.0)
        {
            Settings::USE_STACKING_ACTION = false;
        }
    else
        {
            Settings::USE_STACKING_ACTION = true;
        }
}

std::vector<double> AnalysisManager::get_ener_grid()
{
    // energy grid, it is already in the Geant4 unit system (50 keV to 100 MeV with 256 values)
    std::vector<double> energy_grid
    {
        0.0500,
        0.0515,
        0.0531,
        0.0547,
        0.0563,
        0.0580,
        0.0598,
        0.0616,
        0.0635,
        0.0654,
        0.0674,
        0.0694,
        0.0715,
        0.0737,
        0.0759,
        0.0782,
        0.0806,
        0.0830,
        0.0855,
        0.0881,
        0.0908,
        0.0935,
        0.0963,
        0.0992,
        0.1022,
        0.1053,
        0.1085,
        0.1118,
        0.1152,
        0.1187,
        0.1223,
        0.1260,
        0.1298,
        0.1337,
        0.1378,
        0.1419,
        0.1462,
        0.1506,
        0.1552,
        0.1599,
        0.1647,
        0.1697,
        0.1749,
        0.1801,
        0.1856,
        0.1912,
        0.1970,
        0.2030,
        0.2091,
        0.2154,
        0.2219,
        0.2287,
        0.2356,
        0.2427,
        0.2500,
        0.2576,
        0.2654,
        0.2734,
        0.2817,
        0.2902,
        0.2990,
        0.3081,
        0.3174,
        0.3270,
        0.3369,
        0.3471,
        0.3576,
        0.3684,
        0.3795,
        0.3910,
        0.4028,
        0.4150,
        0.4276,
        0.4405,
        0.4539,
        0.4676,
        0.4817,
        0.4963,
        0.5113,
        0.5268,
        0.5427,
        0.5592,
        0.5761,
        0.5935,
        0.6115,
        0.6300,
        0.6490,
        0.6687,
        0.6889,
        0.7097,
        0.7312,
        0.7533,
        0.7761,
        0.7996,
        0.8238,
        0.8487,
        0.8744,
        0.9009,
        0.9281,
        0.9562,
        0.9851,
        1.0149,
        1.0456,
        1.0773,
        1.1099,
        1.1435,
        1.1780,
        1.2137,
        1.2504,
        1.2882,
        1.3272,
        1.3674,
        1.4088,
        1.4514,
        1.4953,
        1.5405,
        1.5871,
        1.6352,
        1.6846,
        1.7356,
        1.7881,
        1.8422,
        1.8980,
        1.9554,
        2.0145,
        2.0755,
        2.1383,
        2.2030,
        2.2696,
        2.3383,
        2.4091,
        2.4820,
        2.5570,
        2.6344,
        2.7141,
        2.7962,
        2.8808,
        2.9680,
        3.0578,
        3.1503,
        3.2456,
        3.3438,
        3.4450,
        3.5492,
        3.6566,
        3.7673,
        3.8812,
        3.9987,
        4.1197,
        4.2443,
        4.3727,
        4.5050,
        4.6413,
        4.7818,
        4.9264,
        5.0755,
        5.2291,
        5.3873,
        5.5503,
        5.7182,
        5.8912,
        6.0695,
        6.2531,
        6.4423,
        6.6372,
        6.8380,
        7.0449,
        7.2581,
        7.4777,
        7.7039,
        7.9370,
        8.1771,
        8.4246,
        8.6795,
        8.9421,
        9.2126,
        9.4913,
        9.7785,
        10.0744,
        10.3792,
        10.6932,
        11.0168,
        11.3501,
        11.6935,
        12.0473,
        12.4118,
        12.7873,
        13.1742,
        13.5728,
        13.9835,
        14.4066,
        14.8425,
        15.2915,
        15.7542,
        16.2309,
        16.7220,
        17.2279,
        17.7491,
        18.2862,
        18.8394,
        19.4094,
        19.9967,
        20.6017,
        21.2251,
        21.8672,
        22.5289,
        23.2105,
        23.9128,
        24.6363,
        25.3817,
        26.1496,
        26.9408,
        27.7559,
        28.5957,
        29.4609,
        30.3523,
        31.2706,
        32.2168,
        33.1915,
        34.1958,
        35.2304,
        36.2963,
        37.3945,
        38.5259,
        39.6916,
        40.8925,
        42.1297,
        43.4044,
        44.7177,
        46.0707,
        47.4646,
        48.9007,
        50.3802,
        51.9045,
        53.4750,
        55.0929,
        56.7598,
        58.4771,
        60.2464,
        62.0693,
        63.9472,
        65.8820,
        67.8754,
        69.9290,
        72.0448,
        74.2246,
        76.4703,
        78.7840,
        81.1677,
        83.6236,
        86.1537,
        88.7604,
        91.4459,
        94.2127,
        97.0632,
        100.0000
    };

    return energy_grid;
}


std::vector<double> AnalysisManager::get_MOM_grid()
{
    std::vector<double> mom_grid
    {
        -1.0000,
        -0.9843,
        -0.9685,
        -0.9528,
        -0.9370,
        -0.9213,
        -0.9055,
        -0.8898,
        -0.8740,
        -0.8583,
        -0.8425,
        -0.8268,
        -0.8110,
        -0.7953,
        -0.7795,
        -0.7638,
        -0.7480,
        -0.7323,
        -0.7165,
        -0.7008,
        -0.6850,
        -0.6693,
        -0.6535,
        -0.6378,
        -0.6220,
        -0.6063,
        -0.5906,
        -0.5748,
        -0.5591,
        -0.5433,
        -0.5276,
        -0.5118,
        -0.4961,
        -0.4803,
        -0.4646,
        -0.4488,
        -0.4331,
        -0.4173,
        -0.4016,
        -0.3858,
        -0.3701,
        -0.3543,
        -0.3386,
        -0.3228,
        -0.3071,
        -0.2913,
        -0.2756,
        -0.2598,
        -0.2441,
        -0.2283,
        -0.2126,
        -0.1969,
        -0.1811,
        -0.1654,
        -0.1496,
        -0.1339,
        -0.1181,
        -0.1024,
        -0.0866,
        -0.0709,
        -0.0551,
        -0.0394,
        -0.0236,
        -0.0079,
        0.0079,
        0.0236,
        0.0394,
        0.0551,
        0.0709,
        0.0866,
        0.1024,
        0.1181,
        0.1339,
        0.1496,
        0.1654,
        0.1811,
        0.1969,
        0.2126,
        0.2283,
        0.2441,
        0.2598,
        0.2756,
        0.2913,
        0.3071,
        0.3228,
        0.3386,
        0.3543,
        0.3701,
        0.3858,
        0.4016,
        0.4173,
        0.4331,
        0.4488,
        0.4646,
        0.4803,
        0.4961,
        0.5118,
        0.5276,
        0.5433,
        0.5591,
        0.5748,
        0.5906,
        0.6063,
        0.6220,
        0.6378,
        0.6535,
        0.6693,
        0.6850,
        0.7008,
        0.7165,
        0.7323,
        0.7480,
        0.7638,
        0.7795,
        0.7953,
        0.8110,
        0.8268,
        0.8425,
        0.8583,
        0.8740,
        0.8898,
        0.9055,
        0.9213,
        0.9370,
        0.9528,
        0.9685,
        0.9843,
        1.0000,
    };

    return mom_grid;
}
