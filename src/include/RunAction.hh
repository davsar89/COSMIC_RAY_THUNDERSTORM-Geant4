#pragma once

#include "G4UserRunAction.hh"
#include "globals.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction;
class Run;
class PrimaryGeneratorAction;
class HistoManager;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction: public G4UserRunAction
{
    public:

        RunAction();
        ~RunAction();

    public:

        virtual void
        BeginOfRunAction(const G4Run *);
        virtual void
        EndOfRunAction(const G4Run *);

    private:

};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
