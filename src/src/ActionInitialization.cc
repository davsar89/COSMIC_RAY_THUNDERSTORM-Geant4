#include "ActionInitialization.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization(DetectorConstruction *detector)
    : G4VUserActionInitialization(), fDetector(detector) {}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{
    PrimaryGeneratorAction *primary = new PrimaryGeneratorAction();
    SetUserAction(primary);
    RunAction *runAction = new RunAction();
    SetUserAction(runAction);
    EventAction *event = new EventAction();
    SetUserAction(event);

    //    if (Settings::efield_status == Settings::ON)
    //        {
    //    G4bool use_stckg_act_tmp = Settings::USE_STACKING_ACTION; // variable for debug

    if (Settings::USE_STACKING_ACTION)
        {
            G4UserStackingAction *stackingAction = new BaseStackingAction();
            SetUserAction(stackingAction);
        }

    //        }

    //  TrackingAction *trackingAction = new TrackingAction(fDetector);
    //  SetUserAction(trackingAction);
    SteppingAction *steppingAction = new SteppingAction(fDetector, event);
    SetUserAction(steppingAction);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
