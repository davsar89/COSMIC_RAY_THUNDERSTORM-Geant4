
#include "G4RunManager.hh"
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "G4UImanager.hh"
#include "PhysicsList.hh"
#include "myExceptionHandler.hh"

#define G4VIS_USE
#define G4UI_USE

#ifdef G4VIS_USE

# include "G4VisExecutive.hh"

#endif // ifdef G4VIS_USE

#ifdef G4UI_USE

# include "G4UIExecutive.hh"

#endif // ifdef G4UI_USE

#include <chrono>
#include "G4PhysListFactory.hh"
#include "myUtils.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void set_mode_to_run_if_cluster();

int main(int argc, char **argv) {

    const char *env_var[2] = {"CC", "CXX"};
    char *env_val[2];
    for (int i = 0; i < 2; ++i) {
        env_val[i] = getenv(env_var[i]);
        if (env_val[i] != NULL)
            G4cout << "Variable = " << env_var[i] << ", Value= " << env_val[i] << G4endl;
        else {
            G4cout << env_var[i] << " doesn't exist" << G4endl;
        }
    }

    const double WT1 = myUtils::get_wall_time();

    G4cout << "Total available RAM : " << myUtils::check_available_RAM() << " GB" << G4endl;

    Settings *settings = Settings::getInstance();

    if (settings->spread_initial_sampling == false && settings->pinhole_initial_sampling == false) {
        G4cout << "Error. Both spread_initial_sampling and pinhole_initial_sampling are false" << G4endl;
        return 0;
    }

    setlocale(LC_ALL, "C");  // just in case, to avoid potential bug for PARMA ("," <-> ".")

    G4String NB_PARTICLES_TO_SHOOT = "100000";

    if (argc >= 4) {
        settings->POTENTIAL_VALUE = std::stod(argv[1]);
        NB_PARTICLES_TO_SHOOT = argv[2];
        settings->RECORD_POSITION = std::stoi(argv[3]);
        settings->EFIELD_REGION_Y_CENTER = std::stod(argv[4]);
        settings->EFIELD_REGION_Y_FULL_LENGTH = std::stod(argv[5]);
    } else {
        G4cout << "WARNING: Executable should have 6 input arguments. Now running with default settings..." << G4endl;
    }

//    if (settings->POTENTIAL_VALUE == -20.0) {
//        settings->POTENTIAL_VALUE = -21.2;
//    }

    const double rec_position_d = double(settings->RECORD_POSITION);
    settings->RECORD_ALTITUDE = settings->EFIELD_REGION_Y_CENTER + rec_position_d * settings->EFIELD_REGION_Y_FULL_LENGTH / 1.999;
    if (settings->RECORD_ALTITUDE < -0.2 || settings->RECORD_ALTITUDE > 20.2) {
        G4cout << "Record altitude is below 0 or above 20 km, aborting." << G4endl;
        return 0;
    }


    if (settings->pinhole_initial_sampling) {
        if (settings->RECORD_ALTITUDE > settings->EFIELD_REGION_Y_CENTER) {
            settings->CR_GENERATION_ALT_MIN = settings->RECORD_ALTITUDE + settings->DELTA_ALT_INISAMPLED_RECORD;
            settings->CR_GENERATION_ALT_MAX = settings->CR_GENERATION_ALT_MIN + 0.05;
        } else {
            settings->CR_GENERATION_ALT_MIN =
                    settings->EFIELD_REGION_Y_CENTER + settings->EFIELD_REGION_Y_FULL_LENGTH / 1.999 + settings->DELTA_ALT_INISAMPLED_RECORD;
            settings->CR_GENERATION_ALT_MAX = settings->CR_GENERATION_ALT_MIN + 0.05;
        }
        const double alt_diff_between_sample_and_record = std::abs(settings->CR_GENERATION_ALT_MIN - settings->RECORD_ALTITUDE);
        const double fact = std::tan(settings->delta_angle_pinhole_degree * CLHEP::degree);
        settings->RECORD_XZ_HALF_SIZE = alt_diff_between_sample_and_record * fact;
        settings->EFIELD_XZ_HALF_SIZE = settings->RECORD_XZ_HALF_SIZE * 1.5; // km
        settings->CR_SAMPLING_XZ_HALF_SIZE = 0.0; // km
        G4cout << alt_diff_between_sample_and_record << " " << settings->RECORD_XZ_HALF_SIZE << " " << settings->EFIELD_XZ_HALF_SIZE << G4endl;

    } else if (settings->spread_initial_sampling) {
        if (settings->RECORD_ALTITUDE > settings->EFIELD_REGION_Y_CENTER) {
            settings->CR_GENERATION_ALT_MIN =
                    settings->EFIELD_REGION_Y_CENTER - settings->EFIELD_REGION_Y_FULL_LENGTH / 2.0 - settings->DELTA_ALT_INISAMPLED_RECORD;
            settings->CR_GENERATION_ALT_MAX = settings->RECORD_ALTITUDE + settings->DELTA_ALT_INISAMPLED_RECORD;
        } else if (settings->RECORD_ALTITUDE <= settings->EFIELD_REGION_Y_CENTER) {
            settings->CR_GENERATION_ALT_MIN = settings->RECORD_ALTITUDE - settings->DELTA_ALT_INISAMPLED_RECORD;
            settings->CR_GENERATION_ALT_MAX =
                    settings->EFIELD_REGION_Y_CENTER + settings->EFIELD_REGION_Y_FULL_LENGTH / 2.0 + settings->DELTA_ALT_INISAMPLED_RECORD;
        }
    } else {
        if (settings->RECORD_ALTITUDE > settings->EFIELD_REGION_Y_CENTER) {
            settings->CR_GENERATION_ALT_MIN =
                    settings->RECORD_ALTITUDE + settings->DELTA_ALT_INISAMPLED_RECORD;
            settings->CR_GENERATION_ALT_MAX = settings->CR_GENERATION_ALT_MIN + 0.05;
        } else {
            settings->CR_GENERATION_ALT_MIN =
                    settings->EFIELD_REGION_Y_CENTER + settings->EFIELD_REGION_Y_FULL_LENGTH / 2.0 + settings->DELTA_ALT_INISAMPLED_RECORD;
            settings->CR_GENERATION_ALT_MAX = settings->CR_GENERATION_ALT_MIN + 0.05;
        }
    }

    if (settings->CR_GENERATION_ALT_MIN < 0.05) {
        settings->CR_GENERATION_ALT_MIN = 0.05;
    }

    if (settings->POTENTIAL_VALUE != 0.) {
        settings->initial_efield_status = settings->efield_ON;
        settings->current_efield_status = settings->efield_ON;
    } else if (settings->POTENTIAL_VALUE == 0.) {
        settings->initial_efield_status = settings->efield_OFF;
        settings->current_efield_status = settings->efield_OFF;
    }

    if (settings->DELTA_ALT_INISAMPLED_RECORD < settings->EFIELD_REGION_Y_FULL_LENGTH / 2.0) {
        G4cout << "DELTA_ALT_INISAMPLED_RECORD is too small compared to EFIELD_REGION_Y_FULL_LENGTH" << G4endl;
        return 0;
    }

    ////////////////////////////////////////////////////////////////////////////////

    settings->RANDOM_SEED = myUtils::generate_a_unique_ID();

    set_mode_to_run_if_cluster();

    G4cout << "First random seed for GEANT4: " << settings->RANDOM_SEED << G4endl;

    if (settings->EFIELD_REGION_Y_CENTER > 20.0) {
        G4cout << "ERROR. Efield should be set with a center below 20 km" << G4endl;
        return 0;
    }

    // choose the Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::MixMaxRng);
    long seeds[2];
    seeds[0] = settings->RANDOM_SEED;
    seeds[1] = myUtils::generate_a_unique_ID();
    CLHEP::HepRandom::setTheSeeds(seeds, 2);

    AnalysisManager *analysis = AnalysisManager::getInstance();

    auto *runManager = new G4RunManager;

    runManager->SetVerboseLevel(0);

    // set mandatory initialization classes
    auto *det = new DetectorConstruction;
    runManager->SetUserInitialization(det);

//    myExceptionHandler* myExceptionHand = new myExceptionHandler(0);



    runManager->SetVerboseLevel(0);

//    G4bool EM_ONLY = false;
//    G4VUserPhysicsList *phys = new PhysicsList(EM_ONLY);
//    phys->SetVerboseLevel(0);
//    runManager->SetUserInitialization(phys);

//    G4VUserPhysicsList *physicsList = physListFactory->GetReferencePhysList("QGSP_BERT_HP");

    if (settings->PDG_LIST_INITIAL.size() == 3) {
        const G4bool EM_ONLY = true;
        G4VUserPhysicsList *phys = new PhysicsList(EM_ONLY);
        phys->SetVerboseLevel(0);
        runManager->SetUserInitialization(phys);
    } else {
        auto *physListFactory = new G4PhysListFactory();
        G4VUserPhysicsList *physicsList = physListFactory->GetReferencePhysList("QGSP_BERT");
        physicsList->SetDefaultCutValue(10 * m);
        physicsList->SetCutsWithDefault();
        runManager->SetUserInitialization(physicsList);
    }

    runManager->SetUserInitialization(new ActionInitialization(det));

    // get the pointer to the User Interface manager
    G4UImanager *UI = G4UImanager::GetUIpointer();
    UI->SetVerboseLevel(0);

    if (settings->MODE == "run") {
        // Initialize G4 kernel
        runManager->Initialize();

        UI->ApplyCommand(G4String("/tracking/verbose 0"));
        UI->ApplyCommand(G4String("/stepping/verbose 0"));

        UI->ApplyCommand("/run/beamOn " + NB_PARTICLES_TO_SHOOT);

        analysis->write_output_file();

    } else if (settings->MODE == "visu") // define visualization and UI terminal for interactive mode
    {
        runManager->Initialize();
#ifdef G4VIS_USE
        G4VisManager *visManager = new G4VisExecutive;
        visManager->SetVerboseLevel(0);
        visManager->Initialize();
#endif // ifdef G4VIS_USE
#ifdef G4UI_USE
        G4UIExecutive *ui_ex = new G4UIExecutive(argc, argv);
        UI->ApplyCommand("/control/execute vis.mac");
        UI->SetVerboseLevel(0);
        ui_ex->SessionStart();
        delete ui_ex;
#endif // ifdef G4UI_USE
#ifdef G4VIS_USE
        delete visManager;
#endif // ifdef G4VIS_USE
    } else {
        G4cout << G4endl << "ERROR : Mode should be set to 'visu' or 'run'. Aborting. " << G4endl;
        std::abort();
    }

    // job termination

    delete runManager;

    G4double WT2 = myUtils::get_wall_time();

    G4cout << "GEANT4 run finished" << G4endl;
    G4cout << "Potential: " << settings->POTENTIAL_VALUE << G4endl;
    G4cout << "Time taken: " << (WT2 - WT1) / 1.0e6 << " seconds" << G4endl;

    return 0;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void set_mode_to_run_if_cluster() {

    Settings *settings = Settings::getInstance();

    char hostname[HOST_NAME_MAX];
    int result = gethostname(hostname, HOST_NAME_MAX);

    const G4String hostname_str = G4String(hostname);
    const G4String to_find = "iftrom";

    bool IS_FRAM = true;
    if (!(hostname_str.find(to_find) == std::string::npos)) {
        IS_FRAM = false;
    }

    G4cout << hostname_str << G4endl;

    if (IS_FRAM) {
        settings->MODE = "run";
    }
}