//
// ********************************************************************
#include "StackingAction.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

// Added to mimic time dependent simulation, and to count the number of energetic electron every 5 microseconds
// and turn OFF the electric field if there is more than 100000 (-> RREA)

// can be a bad idea to set it on, because it can use a lot of memory compared to default G4 behaviour

BaseStackingAction::BaseStackingAction()
{
    TIME_STEP = Settings::DELTA_T ;
}

BaseStackingAction::~BaseStackingAction()
{  }

void BaseStackingAction::PrepareNewEvent()
{
    LIST_ENERGETIC_PART_IDS.clear();

    VARIABLE_TIME_LIMIT = 11.0 * millisecond; // propagation time of light from 50 km to about 15 km
    EVENT_NB++;

    Settings::current_efield_status = Settings::initial_efield_status ;

    Settings::VARIABLE_TIME_LIMIT = VARIABLE_TIME_LIMIT;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ClassificationOfNewTrack BaseStackingAction::ClassifyNewTrack(const G4Track *aTrack)
{
    if (Settings::current_efield_status == Settings::OFF) return fUrgent;

    if (!is_inside_eField_region(aTrack->GetPosition().y() / km,
                                 aTrack->GetPosition().x() / km,
                                 aTrack->GetPosition().z() / km))
        {
            return fUrgent;    // ignore if we are not in the E-field region
        }

    // Warning : dot not use "aTrack->GetStep()->GetPreStepPoint()", this is not the stepping action, step may not exist

    const G4int ID = aTrack->GetTrackID();

    // if energy > 500 keV and if the ID is not already saved
    if (does_not_contain(ID, LIST_ENERGETIC_PART_IDS))
        {
            if (aTrack->GetKineticEnergy() > ENER_THRES)
                {
                    LIST_ENERGETIC_PART_IDS.push_back(ID); // save the ID

                    check_PART_NB_LIMIT();
                }
        }

    if (aTrack->GetGlobalTime() > VARIABLE_TIME_LIMIT)
        {
            return fWaiting;
        }

    // default classification
    return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BaseStackingAction::print_status()
{
    G4cout << "current max time : " << VARIABLE_TIME_LIMIT / microsecond << " microsecond" << G4endl;
    G4cout << ">0.1 MeV electron count : " << LIST_ENERGETIC_PART_IDS.size() << G4endl;
    G4cout << "Event nb : " << EVENT_NB << G4endl;
    //    G4cout << part_name << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BaseStackingAction::NewStage() // called when the "Urgent stack" is empty
// and the particles in the "waiting stack" are transfered to the "Urgent stack"
{
    //    check_PART_NB_LIMIT();

    if (Settings::current_efield_status == Settings::OFF) return;

    LIST_ENERGETIC_PART_IDS.clear();

    VARIABLE_TIME_LIMIT += TIME_STEP;
    //    LIST_ENERGETIC_PART_IDS.clear();
    //    print_status();
    Settings::VARIABLE_TIME_LIMIT = VARIABLE_TIME_LIMIT;
    return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BaseStackingAction::check_PART_NB_LIMIT()
{
    if (LIST_ENERGETIC_PART_IDS.size() > ENERGETIC_PART_NB_LIMIT)
        {
            Settings::current_efield_status = Settings::OFF;
            LIST_ENERGETIC_PART_IDS.clear();
            Settings::RREA_PART_NB_LIMIT_HAS_BEEN_REACHED = 1;
        }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

bool BaseStackingAction::does_not_contain(const G4int ID, const std::vector<G4int> &LIST_IDS)
{
    if (LIST_IDS.empty()) return true;

    if (std::find(LIST_IDS.begin(), LIST_IDS.end(), ID) != LIST_IDS.end()) return false;
    else return true;
}


bool BaseStackingAction::is_inside_eField_region(const G4double &alt, const G4double &xx, const G4double &zz)
// alt assumed in km
{
    if (alt > alt_min && alt < alt_max
        && abs(xx) < Settings::EFIELD_XY_HALF_SIZE
        && abs(zz) < Settings::EFIELD_XY_HALF_SIZE)
        {
            return true;
        }
    else
        {
            return false;
        }
}

