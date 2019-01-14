#pragma once

#include "AnalysisManager.hh"
#include "G4ThreeVector.hh"
#include "G4UserSteppingAction.hh"
#include "RegionInformation.hh"
#include "globals.hh"
#include <vector>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <CLHEP/Units/SystemOfUnits.h>

#include "Settings.hh"

#include <sys/time.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include "myUtils.hh"

class DetectorConstruction;

class EventAction;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct index_found {
    uint index;
    bool found;
};

class SteppingAction : public G4UserSteppingAction {
public:

    SteppingAction(DetectorConstruction *, EventAction *);

    ~SteppingAction() override;

    void UserSteppingAction(const G4Step *aStep) override;

private:

    Settings *settings = Settings::getInstance();

    AnalysisManager *analysis = AnalysisManager::getInstance();

    DetectorConstruction *fDetector = nullptr;
    EventAction *fEventAction = nullptr;

    G4StepPoint *thePrePoint = nullptr;
    G4StepPoint *thePostPoint = nullptr;
    G4Track *theTrack = nullptr;


    bool is_inside_eField_region(const G4double &alt,
                                 const G4double &xx,
                                 const G4double &zz);

    const double EFIELD_alt_min = settings->EFIELD_REGION_Y_CENTER - settings->EFIELD_REGION_Y_FULL_LENGTH / 2.0; // km
    const double EFIELD_alt_max = settings->EFIELD_REGION_Y_CENTER + settings->EFIELD_REGION_Y_FULL_LENGTH / 2.0; // km
    int find_zenith_angle_index(const double za, const std::vector<double> &vector);

    double MAX_TIME=-100; // us
    double calculate_MAX_TIME(const double efield_center,
                              const double efield_fullsize,
                              const double initial_CR_altitude,
                              const double record_altitude);
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
