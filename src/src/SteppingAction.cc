#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4RunManager.hh"
#include "Run.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction *det, EventAction *event) : G4UserSteppingAction(), fDetector(det), fEventAction(event) {

    if (settings->USE_MAX_TIME) {
        MAX_TIME = calculate_MAX_TIME(settings->EFIELD_REGION_Y_CENTER,
                                      settings->EFIELD_REGION_Y_FULL_LENGTH,
                                      settings->CR_GENERATION_ALT_MAX,
                                      settings->RECORD_ALTITUDE);
    }

    G4cout << G4endl;



}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction() = default;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step *aStep) {
    //////// AVOID PUTTING THIS LINE, it will produce incorrect energy record
    //        if (aStep->GetTrack()->GetTrackStatus() != fAlive) return;
    ////////////////////////////////////////////

    theTrack = aStep->GetTrack();
    thePrePoint = aStep->GetPreStepPoint();
    thePostPoint = aStep->GetPostStepPoint();

    ///////////////////////
    const int PDG_num = theTrack->GetParticleDefinition()->GetPDGEncoding();
    ///////////////////////

//    if (settings->CPU_TIME_LIMIT_PER_EVENT_HAS_BEEN_REACHED) {
//        theTrack->SetTrackStatus(fStopAndKill);
//        return;
//    }

    if (thePrePoint->GetGlobalTime() > settings->MAX_POSSIBLE_TIME) {
#ifndef NDEBUG
        G4cout << G4endl << "PARTICLE IS TERMINATED BECAUSE IT EXCEEDED THE TIME LIMIT. " << thePrePoint->GetGlobalTime() / us << G4endl << G4endl;
#endif // ifndef NDEBUG
        theTrack->SetTrackStatus(fStopAndKill);
        return;
    }


    if (thePrePoint->GetPosition().y() < double(settings->SOIL_ALT_MAX * km - 0.1 * km)) {
        //G4cout << " KILLED BECAUSE REACHED MIN ALTITUDE " << G4endl;
        theTrack->SetTrackStatus(fStopAndKill);
        return;
    }

    if (settings->USE_MAX_TIME) {
        if (thePrePoint->GetGlobalTime() > MAX_TIME * CLHEP::microsecond && PDG_num!=22) { // keeping photons anyways
            //G4cout << " KILLED BECAUSE REACHED MAX TIME " << G4endl;
            theTrack->SetTrackStatus(fStopAndKill);
            return;
        }
    }

    if (thePrePoint->GetPosition().y() > (settings->WORLD_MAX_ALT - 2) * km) {
        theTrack->SetTrackStatus(fStopAndKill);
        //G4cout << " KILLED BECAUSE REACHED MAX ALTITUDE " << G4endl;
        return;
    }

    const double x_tmp = thePrePoint->GetPosition().x();
    const double z_tmp = thePrePoint->GetPosition().z();
    const double radial_dist = sqrt(x_tmp * x_tmp + z_tmp * z_tmp);

    if (radial_dist > settings->RECORD_XZ_HALF_SIZE * km) {
        //G4cout << " SKIPED BECAUSE OUTSIDE OF HORIZONTAL RECORD DISK " << G4endl;
        return;
    }


    if (PDG_num == -12 || PDG_num == 12) { // killing if neutrino
        //G4cout << " KILLED BECAUSE IT IS NEUTRINO " << G4endl;
        theTrack->SetTrackStatus(fStopAndKill);
        return;
    }

    if (PDG_num == 2212 && thePrePoint->GetKineticEnergy() < 50.0 * MeV) {
        //G4cout << " KILLED BECAUSE IT LOW ENERGY PROTON " << G4endl;
        theTrack->SetTrackStatus(fStopAndKill);
        return;
    }

    // cleaning particles below 8 keV to improve performance
    if (thePrePoint) {
        // avoid killing positrons below low energy threshold to make sure we still have annihilation

        if (thePrePoint->GetKineticEnergy() < settings->ENERGY_MIN_RECORD) {
            if (PDG_num != settings->pdg_posi) {
                //G4cout << " KILLED BECAUSE BELOW  LOW_ENERGY_THRES" << G4endl;
                theTrack->SetTrackStatus(fStopAndKill);
                return;
            }
            // using fStopAndAlive for positrons makes a bug that saturates the RAM...
        }
    }

    // skipping the rest if current particle is not photon electron or positron
    if (!(PDG_num == 22 || PDG_num == 11 || PDG_num == -11)) {
        //G4cout << " SKIPPED RECORD BECAUSE NOT PHOTON ELECTRON OR POSITRON" << G4endl;
        return;
    }

    const int part_index_record = settings->particle_PDG_to_index_mapping[PDG_num];

    /// stacking action check, if pseudo time oriented simulation is activated (slow)
    if (settings->USE_STACKING_ACTION) {
        if (thePrePoint) {
            if (settings->current_efield_status == settings->efield_ON) {
                if (theTrack->GetGlobalTime() > settings->VARIABLE_TIME_LIMIT_GLOBAL) {
                    theTrack->SetTrackStatus(fSuspend);
                    return;
                }
            }
        }
    }

    // check if particle should be recorded

    if (thePrePoint->GetStepStatus() == fGeomBoundary || thePostPoint->GetStepStatus() == fGeomBoundary) {

        const G4String &vol_name = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();

        if (vol_name.find("record_phys") != std::string::npos && thePrePoint->GetPosition().y() == settings->RECORD_ALTITUDE * km) {

            const double ener = thePrePoint->GetKineticEnergy() * CLHEP::MeV;

            const double momy = thePrePoint->GetMomentumDirection().y();

            // WARNING : PARTICLES ARE ALLOWED TO BE RECORED 2 TIMES if they cross the detection altitude with upwards and downwards momentum

            if ((ener >= settings->ENERGY_MIN_RECORD) && (ener <= settings->ENERGY_MAX_RECORD)) {

//                if (analysis->check_record(theTrack->GetTrackID(), momy)) {

                analysis->add_NB_OUTPUT();

                const double za = std::acos(momy) * 180.0 / CLHEP::pi;

                const int i_za = find_zenith_angle_index(za, analysis->ZENITH_ANGLE_GRID);

                analysis->fill_histogram_E(part_index_record, i_za, ener / CLHEP::MeV);

                analysis->counter_total[part_index_record] += 1;

                if (momy > 0.0) {
                    analysis->counter_up[part_index_record] += 1;
                } else {
                    analysis->counter_down[part_index_record] += 1;
                }
//                }
            }
        }
    }
}

// ------------------------------------------------------------------------

bool SteppingAction::is_inside_eField_region(const double &alt, const double &xx, const double &zz)
// alt, xx, zz assumed in km
{
    return alt > EFIELD_alt_min && alt < EFIELD_alt_max && std::abs(xx) < settings->EFIELD_XZ_HALF_SIZE && std::abs(zz) < settings->EFIELD_XZ_HALF_SIZE;
}

int SteppingAction::find_zenith_angle_index(const double za, const std::vector<double> &vector) {
    for (uint ii = 0; ii < vector.size() - 1; ++ii) {
        if (vector[ii] <= za && vector[ii + 1] >= za) {
            return ii;
        }
    }
    G4cout << "ERROR: zenith angle out of possible range. Aborting." << G4endl;
    std::abort();
}

// ------------------------------------------------------------------------
// takes into account twice speed-of-light propagation from emission to record or emission to efield-edge
// and 5 RREA times or given E-field center density
// output un MICROSECONDS
double SteppingAction::calculate_MAX_TIME(const double efield_center,
                                          const double efield_fullsize,
                                          const double initial_CR_altitude,
                                          const double record_altitude) {


    // max distance for propagation in atmosphere
    const double efield_top = efield_center + efield_fullsize / 2.0;
    const double efield_bottom = efield_center + efield_fullsize / 2.0;

    std::vector<double> distances{settings->CR_GENERATION_ALT_MAX,
                                  efield_center,
                                  efield_top,
                                  efield_bottom,
                                  record_altitude,
                                  settings->CR_GENERATION_ALT_MIN,
    };

    std::vector<double> all_possible_diff{};

    for (const double &dist1 : distances) {
        for (const double &dist2 : distances) {
            all_possible_diff.push_back(std::abs(dist1 - dist2));
        }
    }

    double max_Dist = 0;

    for (const double &dist : all_possible_diff) {
        if (dist > max_Dist) {
            max_Dist = dist;
        }
    }

    // max distance for 4 back-and-forth in e-field region

    const double efield_region_lengths = efield_fullsize * 8.0;

    // accumulation of all distances with some margin

    const double max_distance_accumulated = 5.5 * max_Dist + 6.0 * efield_region_lengths; //

    const double c_km_per_us = 0.299792458; // speed of light in km per microsecond

    return max_distance_accumulated / c_km_per_us;
}
