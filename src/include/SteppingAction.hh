#pragma once

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
#include "AnalysisManager.hh"
#include "RegionInformation.hh"

#include <CLHEP/Units/SystemOfUnits.h>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

#include "Settings.hh"

#include <time.h>
#include <sys/time.h>

class DetectorConstruction;
class EventAction;

struct detected_part
{
    int direction;
    int ID;
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SteppingAction: public G4UserSteppingAction
{
    public:
        SteppingAction(DetectorConstruction *, EventAction *);
        ~SteppingAction();

        virtual void
        UserSteppingAction(const G4Step *aStep);

    private:
        AnalysisManager *analysis = AnalysisManager::getInstance();
        DetectorConstruction *fDetector;
        EventAction *fEventAction;
        G4StepPoint *thePrePoint = 0;

        const G4int PDG_phot = 22;
        const G4int PDG_elec = 11;
        const G4int PDG_posi = -11;

        const G4int NB_alt = Settings::RECORD_ALTITUDES.size();
        double get_wall_time();

        double computation_length_for_event_limit = 600.; // 10 minutes

        G4int part_ID;
        G4int previous_part_ID;
        G4bool is_not_recorded_ID(const detected_part &det_part);
        G4bool not_contains(const detected_part &x, const std::vector<detected_part> &v);

        std::vector<detected_part> RECORDED_LIST;


        bool is_inside_eField_region(const G4double &alt, const G4double &xx, const G4double &zz);
        G4double alt_min = Settings::EFIELD_REGION_ALT_CENTER - Settings::EFIELD_REGION_LEN / 2.0; // km
        G4double alt_max = Settings::EFIELD_REGION_ALT_CENTER + Settings::EFIELD_REGION_LEN / 2.0; // km
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
