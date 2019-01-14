
#pragma once

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "Settings.hh"
#include "AnalysisManager.hh"

#include <time.h>
#include <sys/time.h>

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction: public G4UserEventAction
{
    public:

        EventAction();
        ~EventAction();

    public:

        virtual void
        BeginOfEventAction(const G4Event *);
        virtual void
        EndOfEventAction(const G4Event *);

    private:

        G4int nb_detectors = 3;

        G4int print_nb = 1; // just initialisation

        AnalysisManager *analysis = AnalysisManager::getInstance();
        double get_wall_time() const;

        double time_begin_event = -5.;
        double time_end_event = -5.;

        double max_event_duration = -10.;
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
