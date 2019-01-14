// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include <chrono>
#include "EventAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "Run.hh"
#include "myUtils.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction() : G4UserEventAction() {

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() {}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *) {
    settings->NB_EVENT++;
//#ifndef NDEBUG // debug mode

    if (settings->NB_EVENT % print_nb == 0) {

        auto end = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        G4cout << "Begin of event : " << settings->NB_EVENT << "; current time: " << std::ctime(&end_time) << G4endl;
    }

//#endif // ifndef NDEBUG

    settings->current_efield_status = settings->initial_efield_status;

} // EventAction::BeginOfEventAction

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *) {


} // EventAction::EndOfEventAction