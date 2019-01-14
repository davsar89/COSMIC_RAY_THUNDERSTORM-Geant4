#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "AnalysisManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

#ifdef G4VIS_USE
# include "G4VisExecutive.hh"
#endif // ifdef G4VIS_USE

#ifdef G4UI_USE
# include "G4UIExecutive.hh"
#endif // ifdef G4UI_USE

#include <fstream>
#include "time.h"
#include "G4UnitsTable.hh"
#include "Settings.hh"
#include <locale.h>
#include "myExceptionHandler.hh"
#include "G4PhysListFactory.hh"

#include "G4VUserPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"

#include <fstream>
#include <chrono>
// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// ------------------------------------------------------------------------

double get_wall_time_main()
// returns time in seconds
{
    struct timeval tv;

    gettimeofday(&tv, NULL);
    return tv.tv_sec + (tv.tv_usec / 1000000.0);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char **argv)
{
    setlocale(LC_ALL, "C"); // just in case, to avoid potential bug ("," <-> ".")

    G4String Mode                  = "run"; // visualization (visu) or run
    G4String NB_PARTICLES_TO_SHOOT = "20";
    G4int    NB_PARTICLES_TO_GET   = 1000;

    G4double WT1 = get_wall_time_main();

    std::chrono::high_resolution_clock m_clock;
    long start = std::chrono::duration_cast<std::chrono::nanoseconds>(m_clock.now().time_since_epoch()).count();;
    G4cout << start << " ns" << G4endl;

    Settings::RANDOM_SEED        = start;

    if (argc >= 4)
        {
            NB_PARTICLES_TO_GET          = std::stoi(argv[1]);
            Settings::EFIELD_REGION_ALT_CENTER = std::stod(argv[2]);
            Settings::EFIELD_REGION_LEN = std::stod(argv[3]);
            Settings::POTENTIAL_VALUE    = std::stod(argv[4]);
            Settings::TILT    = std::stod(argv[5]);

            Settings::RECORD_ALTITUDES.clear();

            for (int kk = 6; kk < argc; ++kk) // loop over all requested record altitudes
                {
                    Settings::RECORD_ALTITUDES.push_back(std::stod(argv[kk]));
                    //                    G4cout << argv[kk] << G4endl;
                }

            Mode = "run";

            if (Settings::EFIELD_REGION_ALT_CENTER > 20.0)
                {
                    G4cout << "ERROR. For this version of the code, Efield should be set with a center below 20 km" << G4endl;
                    std::abort();
                }
            else Settings::current_efield_status = Settings::BELOW_PLANE;

            if (Settings::POTENTIAL_VALUE == 0.)
                {
                    Settings::current_efield_status = Settings::OFF;
                }

            Settings::initial_efield_status = Settings::current_efield_status;

        }
    else
        {
            Mode = "run";
            Settings::RECORD_ALTITUDES.push_back(13.0); // km
            //            Settings::RECORD_ALTITUDES.push_back(12.0);
            //            Settings::RECORD_ALTITUDES.push_back(12.5);
        }

    if (Settings::EFIELD_REGION_ALT_CENTER > 20.0)
        {
            G4cout << "ERROR. For this version of the code, Efield should be set with a center below 20 km" << G4endl;
            std::abort();
        }

    else Settings::current_efield_status = Settings::BELOW_PLANE;

    // choose the Random engine
    //  G4Random::setTheEngine(new CLHEP::MTwistEngine);
    // CLHEP::HepRandom::setTheEngine(new CLHEP::MTwistEngine);
    // random seed with system time
    CLHEP::HepRandom::setTheSeed(start);
    // CLHEP::HepRandom::showEngineStatus();

    AnalysisManager *analysis = AnalysisManager::getInstance();

    analysis->check_if_should_use_stacking_action();

#ifdef G4MULTITHREADED
    G4MTRunManager *runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());
#else
    G4RunManager *runManager = new G4RunManager;
#endif

    //    myExceptionHandler *ExceptionHandler = new myExceptionHandler(0); // G4ExceptionHandler class with verbosity option (0 -> silent)

    runManager->SetVerboseLevel(0);
    // set mandatory initialization classes
    DetectorConstruction *det = new DetectorConstruction;
    runManager->SetUserInitialization(det);

    G4bool EM_ONLY = false;
    PhysicsList *phys = new PhysicsList(EM_ONLY);
    phys->SetVerboseLevel(0);
    runManager->SetUserInitialization(phys);

    //    G4PhysListFactory *physListFactory = new G4PhysListFactory();
    //    G4VUserPhysicsList *physicsList =   physListFactory->GetReferencePhysList("Shielding");
    //    runManager ->SetUserInitialization(physicsList);

    //    G4PhysListFactory *physListFactory = new G4PhysListFactory();
    //    G4VUserPhysicsList *physicsList =   physListFactory->GetReferencePhysList("QGSP_BERT");
    //    runManager ->SetUserInitialization(physicsList);

    runManager->SetUserInitialization(new ActionInitialization(det));

    // get the pointer to the User Interface manager
    G4UImanager *UI = G4UImanager::GetUIpointer();
    UI->SetVerboseLevel(0);

    if (Mode == "run")
        {
            // Initialize G4 kernel
            runManager->Initialize();
            //     UI->ApplyCommand("/run/printProgress 1000000");
            //            AnalysisManager *analysis = AnalysisManager::getInstance();

            while (analysis->NB_OUTPUT() < NB_PARTICLES_TO_GET)
                {
                    UI->ApplyCommand("/run/beamOn " + NB_PARTICLES_TO_SHOOT);
                }

            // to make sure we recorded at least 5 positron for each run  to avoid potential bias (put 5 to have a bit of margin)
            //            while (analysis->NB_OUTPUT() < 5)
            //                {
            //                    UI->ApplyCommand("/run/beamOn 5");
            //                }

            analysis->write_output_file_endOf_program();

            //                }
        }
    else if (Mode == "visu")     // define visualization and UI terminal for interactive mode
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
        }
    else
        {
            G4cout << G4endl << "ERROR : Mode should be set to 'visu' or 'run'. Aborting. " << G4endl;
            std::abort();
        }

    // job termination
    //
    delete runManager;

    G4double WT2 = get_wall_time_main();

    G4cout << "GEANT4 run finished" << G4endl;
    G4cout << "Time taken: " << WT2 - WT1 << G4endl;
    return 0;
} // main




