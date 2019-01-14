#include "SD.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4StepPoint.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "CLHEP/Units/SystemOfUnits.h"

//

// IMPORTANT : THERE IS ONE SENSITIVE DETECTOR PER RECORD ALTITUDE (LAYER)

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SensitiveDet::SensitiveDet(G4String name, const G4int ID, const G4double alti_in_km) : G4VSensitiveDetector(name)
{
    ID_SD = ID;
    RECORD_ALT_IN_KM = alti_in_km;

    //    RECORDED_LIST.clear();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SensitiveDet::~SensitiveDet()
{
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SensitiveDet::Initialize(G4HCofThisEvent *) // executed at begin of each event
{
    //    RECORDED_LIST.clear();
    Settings::current_efield_status = Settings::initial_efield_status;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SensitiveDet::ProcessHits(G4Step *aStep, G4TouchableHistory * /*ROhist*/)
{
    given_altitude_particle_record(aStep);

    return true;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SensitiveDet::EndOfEvent(G4HCofThisEvent * /*HCE*/)
{
    // RK : see EndOfEventAction method in EventAction.cc

    //    RECORDED_LIST.clear(); // redundant, but won't hurt
    Settings::current_efield_status = Settings::initial_efield_status;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SensitiveDet::given_altitude_particle_record(const G4Step *aStep)
{
    //    const G4int PDG_nb = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();

    //    if (!(PDG_nb == PDG_phot || PDG_nb == PDG_posi || PDG_nb == PDG_elec)) return;

    //    thePrePoint = aStep->GetPreStepPoint();

    //    if (std::abs(thePrePoint->GetPosition().x()) > Settings::EFIELD_XY_HALF_SIZE * km)
    //        {
    //            return;
    //        }

    //    if (std::abs(thePrePoint->GetPosition().z()) > Settings::EFIELD_XY_HALF_SIZE * km)
    //        {
    //            return;
    //        }

    //    //    if (thePrePoint->GetStepStatus() == fGeomBoundary) // if the particle has just entered the volume ; should not be necessary, but won't hurt
    //    //        {
    //    const G4double ener = thePrePoint->GetKineticEnergy();
    //    const G4int ID_part = aStep->GetTrack()->GetTrackID();
    //    const G4double momy = thePrePoint->GetMomentumDirection().y();

    //    // WARNING : PARTICLES ARE ALLOWED TO BE RECORED TWO TIMES IF COMING FROM DIFFERENT DIRECTIONS

    //    G4int sign_momy = -1;

    //    if (momy > 0.0)
    //        {
    //            sign_momy = 1;
    //        }
    //    else
    //        {
    //            sign_momy = -1;
    //        }

    //    detected_part det_part;
    //    det_part.direction = sign_momy;
    //    det_part.ID = ID_part;

    //    if (ener > Settings::ENERGY_MIN && ener < Settings::ENERGY_MAX && is_not_recorded_ID(det_part))
    //        //            if (aStep->GetTrack()->GetKineticEnergy() > Settings::ENERGY_MIN)
    //        {

    //            if (momy > 0.0)
    //                {
    //                    if (PDG_nb == 22)
    //                        {
    //                            analysis->photon_counter_up++;
    //                        }
    //                    else if (PDG_nb == 11)
    //                        {
    //                            analysis->electron_counter_up++;
    //                        }
    //                    else if (PDG_nb == -11)
    //                        {
    //                            analysis->positron_counter_up++;
    //                        }
    //                }

    //            if (momy < 0.0)
    //                {
    //                    if (PDG_nb == 22)
    //                        {
    //                            analysis->photon_counter_down++;
    //                        }
    //                    else if (PDG_nb == 11)
    //                        {
    //                            analysis->electron_counter_down++;
    //                        }
    //                    else if (PDG_nb == -11)
    //                        {
    //                            analysis->positron_counter_down++;
    //                        }
    //                }

    //            analysis->add_NB_OUTPUT();



    //            //            analysis->register_particle(aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding(),
    //            //                                        ID_part,
    //            //                                        thePrePoint->GetPosition().x() / meter,
    //            //                                        thePrePoint->GetPosition().y() / meter,
    //            //                                        thePrePoint->GetPosition().z() / meter,
    //            //                                        thePrePoint->GetMomentumDirection().x(),
    //            //                                        momy,
    //            //                                        thePrePoint->GetMomentumDirection().z(),
    //            //                                        ener / keV,
    //            //                                        thePrePoint->GetGlobalTime() / microsecond,
    //            //                                        RECORD_ALT_IN_KM);

    //            RECORDED_LIST.push_back(det_part);
    //        }

    //        }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SensitiveDet::get_RECORD_ALT_IN_KM() const
{
    return RECORD_ALT_IN_KM;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int SensitiveDet::get_ID_SD() const
{
    return ID_SD;
}

//// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//G4bool SensitiveDet::is_not_recorded_ID(const detected_part &det_part)
//{
//    return not_contains(det_part, RECORDED_LIST);
//}

//// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//G4bool SensitiveDet::not_contains(const detected_part &x, const std::vector<detected_part> &v)
//{
//    if (v.empty()) return true;

//    for (auto & part : v)
//        {
//            if (part.direction == x.direction && part.ID == x.ID)
//                {
//                    return false;
//                }
//        }

//    return true;
//}
