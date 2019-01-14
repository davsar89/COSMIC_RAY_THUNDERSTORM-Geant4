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

// class following singleton pattern

AnalysisManager *AnalysisManager::instance = nullptr;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AnalysisManager *AnalysisManager::getInstance() // singleton lazy initialization
{
    if (instance == nullptr) {
        instance = new AnalysisManager;
    }

    return instance;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
AnalysisManager::AnalysisManager() // constructor
{
    // initialization of record arrays

    prepare_output_file_folders();

    memset(counter_up, 0, sizeof(counter_up));
    memset(counter_down, 0, sizeof(counter_down));
    memset(counter_total, 0, sizeof(counter_total));

    initialize_energy_spectra_for_angles();

    // sanity checks

    if (ENER_GRID.size() != size_grid_ener) {
        G4cout << "Error in AnalysisManager.hh/cc : energy grid vector does not have 256 elements." << G4endl;
        std::abort();
    }

    if (ZENITH_ANGLE_GRID.size() != size_grid_mom) {
        G4cout << "Error in AnalysisManager.hh/cc : momentum-direction grid vector does not have 128 elements." << G4endl;
        std::abort();
    }

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AnalysisManager::~AnalysisManager() = default;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AnalysisManager::add_NB_OUTPUT() {

    NB_OUTPUT++;

#ifndef NDEBUG // debug mode
//    G4cout << "Number of outputs : " << NB_OUTPUT << G4endl;
#endif // ifndef NDEBUG
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void AnalysisManager::write_output_file() {

    prepare_output_file_folders();

    clean_output_file();

    for (uint i_record_types = 0; i_record_types < 3; ++i_record_types) {

        output_file_buffer << std::scientific << std::setprecision(5)
                           << settings->RANDOM_SEED
                           << " " << settings->NB_EVENT
                           << " " << 777
                           << " " << settings->particle_index_to_PDG_mapping[i_record_types]
                           << " " << settings->EFIELD_REGION_Y_CENTER
                           << " " << settings->EFIELD_REGION_Y_FULL_LENGTH
                           << " " << settings->POTENTIAL_VALUE
                           << " " << 1.0
                           << " " << settings->RECORD_POSITION
                           << " " << settings->CR_GENERATION_ALT_MIN;

        output_file_buffer << G4endl << G4endl;

        output_file_buffer << counter_total[i_record_types] << " ";

        output_file_buffer << G4endl << G4endl;

        output_file_buffer << std::fixed << std::setprecision(1);

        for (uint i_angle = 0; i_angle < ngrid_angle; ++i_angle) {

            output_file_buffer << ZENITH_ANGLE_GRID[i_angle] << " ";
        }

        output_file_buffer << std::scientific << std::setprecision(3);

        output_file_buffer << G4endl << G4endl;

        for (const auto &ener : ENER_GRID) {
            output_file_buffer << ener << " ";
        }

        output_file_buffer << G4endl << G4endl;

        for (uint i_angle = 0; i_angle < ngrid_angle; ++i_angle) {

            std::string str = std::to_string(i_record_types) + "_" + std::to_string(i_angle);

            for (uint ii = 0; ii < ngrid_ener; ++ii) {
                output_file_buffer << PART_SPEC[str][ii] << " ";
            }
            output_file_buffer << G4endl;
        }

        output_file_buffer << G4endl << G4endl;

    }

    write_in_output_file_endOfRun();

}

// ======================================================================
// Returns interpolated value at x from parallel arrays ( xData, yData )
//   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
//   boolean argument extrapolate determines behaviour beyond ends of array (if needed)
G4double AnalysisManager::interpolate(std::vector<G4double> &xData, std::vector<G4double> &yData, G4double x, bool extrapolate) {
    int size = static_cast<int>(xData.size());

    int i = 0;                // find left end of interval for interpolation

    if (x >= xData[size - 2]) // special case: beyond right end
    {
        i = size - 2;
    } else {
        while (x > xData[i + 1]) i++;
    }

    double xL = xData[i], yL = yData[i], xR = xData[i + 1], yR = yData[i + 1]; // points on either side (unless beyond ends)

    if (!extrapolate)                                                          // if beyond ends of array and not extrapolating
    {
        if (x < xL) {
            yR = yL;
        }

        if (x > xR) {
            yL = yR;
        }
    }

    double dydx = (yR - yL) / (xR - xL); // gradient

    return yL + dydx * (x - xL);         // linear interpolation
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::fill_histogram_E(const G4int idx_part_record, const G4int idx_z_angle, const G4double &value) {

    std::string str = std::to_string(idx_part_record) + "_" + std::to_string(idx_z_angle);

    for (uint ii = 0; ii < ngrid_ener - 1; ++ii) {
        if (value >= ENER_GRID[ii] && value < ENER_GRID[ii + 1]) {
            PART_SPEC[str][ii] += 1;
            return;
        }

        // overflow, i.e. > 100 MeV
        if (value >= 100.0) {
            PART_SPEC[str][ngrid_ener - 1] += 1;
            return;
        }
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// open / close with trunc to clean the files (if it does exist) or to create it
void AnalysisManager::clean_output_file() {

    asciiFile_output.open(fileName_output, std::ios::trunc);
    if (asciiFile_output.is_open()) {
        asciiFile_output.close();
    } else {
        G4cout << G4endl << "ERROR : cannot open output file " + fileName_output + ". Aborting" << G4endl;
        std::abort();
    }

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::initialize_energy_spectra_for_angles() {

    for (uint i_record_type = 0; i_record_type < settings->particules_names_record.size(); ++i_record_type) {
        for (uint i_zenith_angle = 0; i_zenith_angle < ngrid_angle; ++i_zenith_angle) {
            G4String str = std::to_string(i_record_type) + "_" + std::to_string(i_zenith_angle);
            std::vector<int> zero_vector(ngrid_ener, 0);
            PART_SPEC[str] = zero_vector;
        }
    }

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::prepare_output_file_folders() {

    char hostname[HOST_NAME_MAX];
    int result = gethostname(hostname, HOST_NAME_MAX);

    const G4String hostname_str = G4String(hostname);
    const G4String to_find = "iftrom";

    bool IS_FRAM = true;
    if (!(hostname_str.find(to_find) == std::string::npos)) {
        IS_FRAM = false;
    }

    G4cout << hostname_str << G4endl;

    G4String BASE_DIR = "";

    if (IS_FRAM) {
        if (settings->pinhole_initial_sampling) {
            BASE_DIR = "/cluster/work/users/dsarria/SIMULATION_DATAFILES/COSMIC_THUNDER_PH/";
        } else {
            BASE_DIR = "/cluster/work/users/dsarria/SIMULATION_DATAFILES/COSMIC_THUNDER/";
        }
        settings->MODE = "run";
    } else {

        if (settings->pinhole_initial_sampling) {
            BASE_DIR = "./output_PH/";
        } else {
            BASE_DIR = "./output/";
        }
    }

    G4cout << "Base dir: " << BASE_DIR << G4endl;


    G4String second_part;
    G4String first_part;

    G4String pot = std::to_string(int(settings->POTENTIAL_VALUE)) + "MV";

    G4String rec_alt = '/' + std::to_string(int(settings->RECORD_POSITION)) + "/";

    G4String EFIELD_CENTER_AND_FULLSIZE = std::to_string(int(settings->EFIELD_REGION_Y_CENTER)) + "km"
                                          + '_' + std::to_string(int(settings->EFIELD_REGION_Y_FULL_LENGTH)) + "km";

    G4String dir_seed = std::to_string(settings->RANDOM_SEED) + "/";

    first_part = BASE_DIR + pot + rec_alt + EFIELD_CENTER_AND_FULLSIZE + "/" + dir_seed;


    G4String mkdir_str = "mkdir -p " + first_part;
    system(mkdir_str);

    second_part = first_part + "ALL" + "_ener_mom_dists_";

    G4String potential_str = "_" + std::to_string(int(settings->POTENTIAL_VALUE));
    G4String efield_alt_str = "_" + std::to_string(int(settings->EFIELD_REGION_Y_CENTER * 10.));
    G4String efield_length_str = "_" + std::to_string(int(settings->EFIELD_REGION_Y_FULL_LENGTH * 100.));
    G4String rec_alt_part = "_" + std::to_string(int(settings->RECORD_POSITION));

    fileName_output = second_part + std::to_string(settings->RANDOM_SEED) + rec_alt_part + potential_str + efield_alt_str + efield_length_str + ".out";

    G4cout << "Reserving output file : " << fileName_output << G4endl;

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void AnalysisManager::write_in_output_file_endOfRun() {

    asciiFile_output.open(fileName_output, std::ios::app);

    if (asciiFile_output.is_open()) {
        asciiFile_output << output_file_buffer.str();
    }

    asciiFile_output.close();
}