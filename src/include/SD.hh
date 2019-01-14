#pragma once
#include "G4Track.hh"
#include "G4VSensitiveDetector.hh"
#include <fstream>
#include "Settings.hh"
#include "AnalysisManager.hh" // singleton
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

// IMPORTANT : THERE IS ONE SENSITIVE DETECTOR PER RECORD ALTITUDE (LAYER)

using namespace std;

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

//struct detected_part
//{
//    int direction;
//    int ID;
//};

class SensitiveDet : public G4VSensitiveDetector
{

    public:
        SensitiveDet(G4String name, const G4int ID, const G4double alti_in_km);
        virtual ~SensitiveDet();

        virtual void Initialize(G4HCofThisEvent *HCE);
        virtual G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist);
        virtual void EndOfEvent(G4HCofThisEvent *HCE);

        G4int get_ID_SD() const;

        G4double get_RECORD_ALT_IN_KM() const;

    private:

        G4int evtctr;
        G4StepPoint *thePrePoint = nullptr;

        AnalysisManager *analysis = AnalysisManager::getInstance();

        void given_altitude_particle_record(const G4Step *aStep);

        const G4int PDG_phot = 22;
        const G4int PDG_elec = 11;
        const G4int PDG_posi = -11;

        const G4int NB_alt = Settings::RECORD_ALTITUDES.size();

        G4int ID_SD = 0; // just initialisation
        G4double RECORD_ALT_IN_KM = 0; // jsut initialisation

//        std::vector<detected_part> RECORDED_LIST;

//        G4bool is_not_recorded_ID(const detected_part &ID_par);
//        G4bool not_contains(const detected_part &x, const std::vector<detected_part> &v);


};

