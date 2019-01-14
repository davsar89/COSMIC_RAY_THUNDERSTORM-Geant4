//

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


#pragma once

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "g4csv.hh"
#include <array>
#include "Settings.hh"
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <limits.h>
#include <unistd.h>
#include <unordered_map>

// #include "g4root.hh"
// #include "g4xml.hh"

// following singleton pattern

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

struct record {
    uint ID;
    bool upward;
    bool downward;
};

class AnalysisManager {
private:

    AnalysisManager();

    static AnalysisManager *instance;

private:


public:

    ~AnalysisManager();

    static AnalysisManager *getInstance();

    void write_output_file();

    void add_NB_OUTPUT();

    void clean_output_file();

    void fill_histogram_E(const G4int idx_part_record,
                          const G4int idx_z_angle,
                          const G4double &value);

    static const uint nb_particules = 3;
    static const uint ngrid_ener = 256;
    static const uint ngrid_angle = 13;

    std::unordered_map< std::string, std::vector<int> > PART_SPEC = {};

    int counter_up[nb_particules];
    int counter_down[nb_particules];
    int counter_total[nb_particules];

    std::stringstream output_file_buffer;

    std::string fileName_output;

    uint NB_OUTPUT = 0;


private:

    Settings *settings = Settings::getInstance();

    std::ofstream asciiFile_output;

    G4double interpolate(std::vector<G4double> &xData,
                         std::vector<G4double> &yData,
                         G4double x,
                         bool extrapolate);

    ////////////////////////////////////////////////////
    ///
    ///
    ////////////////////////////////////////////////////

    const uint size_grid_ener = ngrid_ener;
    const std::vector<double>
            ENER_GRID{
            0.00300,
            0.00313,
            0.00326,
            0.00339,
            0.00353,
            0.00368,
            0.00384,
            0.00400,
            0.00416,
            0.00434,
            0.00452,
            0.00471,
            0.00491,
            0.00511,
            0.00533,
            0.00555,
            0.00578,
            0.00602,
            0.00628,
            0.00654,
            0.00681,
            0.00710,
            0.00739,
            0.00770,
            0.00803,
            0.00836,
            0.00871,
            0.00908,
            0.00946,
            0.00985,
            0.01026,
            0.01069,
            0.01114,
            0.01161,
            0.01209,
            0.01260,
            0.01313,
            0.01368,
            0.01425,
            0.01484,
            0.01547,
            0.01611,
            0.01679,
            0.01749,
            0.01822,
            0.01899,
            0.01978,
            0.02061,
            0.02147,
            0.02237,
            0.02331,
            0.02428,
            0.02530,
            0.02636,
            0.02746,
            0.02861,
            0.02981,
            0.03105,
            0.03235,
            0.03371,
            0.03512,
            0.03659,
            0.03812,
            0.03971,
            0.04138,
            0.04311,
            0.04491,
            0.04679,
            0.04875,
            0.05079,
            0.05292,
            0.05513,
            0.05744,
            0.05984,
            0.06235,
            0.06496,
            0.06767,
            0.07051,
            0.07346,
            0.07653,
            0.07974,
            0.08307,
            0.08655,
            0.09017,
            0.09395,
            0.09788,
            0.10197,
            0.10624,
            0.11069,
            0.11532,
            0.12015,
            0.12518,
            0.13042,
            0.13587,
            0.14156,
            0.14749,
            0.15366,
            0.16009,
            0.16679,
            0.17377,
            0.18104,
            0.18862,
            0.19651,
            0.20474,
            0.21331,
            0.22224,
            0.23154,
            0.24123,
            0.25132,
            0.26184,
            0.27280,
            0.28422,
            0.29612,
            0.30851,
            0.32142,
            0.33487,
            0.34889,
            0.36349,
            0.37871,
            0.39456,
            0.41107,
            0.42827,
            0.44620,
            0.46487,
            0.48433,
            0.50460,
            0.52572,
            0.54772,
            0.57065,
            0.59453,
            0.61941,
            0.64534,
            0.67235,
            0.70049,
            0.72981,
            0.76035,
            0.79217,
            0.82533,
            0.85987,
            0.89586,
            0.93335,
            0.97242,
            1.01312,
            1.05552,
            1.09970,
            1.14572,
            1.19368,
            1.24363,
            1.29568,
            1.34991,
            1.40641,
            1.46528,
            1.52660,
            1.59050,
            1.65706,
            1.72642,
            1.79867,
            1.87395,
            1.95239,
            2.03410,
            2.11923,
            2.20793,
            2.30034,
            2.39662,
            2.49692,
            2.60143,
            2.71031,
            2.82374,
            2.94193,
            3.06505,
            3.19334,
            3.32699,
            3.46624,
            3.61131,
            3.76245,
            3.91993,
            4.08399,
            4.25492,
            4.43300,
            4.61854,
            4.81184,
            5.01323,
            5.22305,
            5.44165,
            5.66940,
            5.90669,
            6.15390,
            6.41146,
            6.67981,
            6.95938,
            7.25065,
            7.55412,
            7.87028,
            8.19968,
            8.54286,
            8.90041,
            9.27292,
            9.66103,
            10.06537,
            10.48664,
            10.92554,
            11.38282,
            11.85923,
            12.35557,
            12.87270,
            13.41146,
            13.97278,
            14.55759,
            15.16687,
            15.80166,
            16.46301,
            17.15204,
            17.86991,
            18.61783,
            19.39705,
            20.20888,
            21.05469,
            21.93590,
            22.85399,
            23.81051,
            24.80706,
            25.84532,
            26.92703,
            28.05402,
            29.22818,
            30.45147,
            31.72597,
            33.05381,
            34.43723,
            35.87854,
            37.38018,
            38.94467,
            40.57463,
            42.27282,
            44.04208,
            45.88539,
            47.80585,
            49.80669,
            51.89126,
            54.06309,
            56.32581,
            58.68324,
            61.13933,
            63.69822,
            66.36421,
            69.14177,
            72.03559,
            75.05053,
            78.19164,
            81.46423,
            84.87378,
            88.42604,
            92.12697,
            95.98279,
            100.00000,
            1000000.00000
    };

    const uint size_grid_mom = ngrid_angle;
public:
    const std::vector<double>
            ZENITH_ANGLE_GRID{0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180};

    void initialize_energy_spectra_for_angles();

    void prepare_output_file_folders();

    void write_in_output_file_endOfRun();
};
