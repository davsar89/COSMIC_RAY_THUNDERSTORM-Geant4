// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "Run.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Settings.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction() : G4UserEventAction()
{
    if (std::abs(Settings::POTENTIAL_VALUE / Settings::EFIELD_REGION_LEN) < 51.0)
        {
            print_nb = 100;
        }
    else
        {
            print_nb = 100;
        }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() {}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *)
{
    Settings::NB_EVENT++;

    Settings::current_efield_status = Settings::initial_efield_status ; // redundant with SD.CC, just to make sure

    if (Settings::USE_WALL_TIME_LIMIT_FOR_EVENT) Settings::wall_T_begin_event = get_wall_time();

#ifndef NDEBUG // debug mode

    if (Settings::NB_EVENT % print_nb == 0)
        {
            G4cout << "Begin of event : " << Settings::NB_EVENT << G4endl;
        }

#endif

    Settings::RREA_PART_NB_LIMIT_HAS_BEEN_REACHED = 0;


    if (Settings::TIME_EVENT_DURATIONS)
        {
            time_begin_event = get_wall_time();
        }

} // EventAction::BeginOfEventAction

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *)
{
    //    analysis->WriteOutputFile_endOfEvent();

    if (Settings::TIME_EVENT_DURATIONS)
        {
            time_end_event = get_wall_time();

            double duration = time_end_event - time_begin_event;

            if (duration > max_event_duration)
                {
                    max_event_duration = duration;
                }

            G4cout << "Max event duration: " << max_event_duration << G4endl;
        }

} // EventAction::EndOfEventAction


// ------------------------------------------------------------------------

double EventAction::get_wall_time() const
// returns time in seconds
{
    struct timeval tv;

    gettimeofday(&tv, NULL);
    return tv.tv_sec + (tv.tv_usec / 1000000.0);
}
