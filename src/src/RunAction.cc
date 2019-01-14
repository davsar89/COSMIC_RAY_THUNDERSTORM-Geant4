#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"


#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <iomanip>
#include "AnalysisManager.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
    : G4UserRunAction()
{
    // Book predefined histograms
    //    fHistoManager = new HistoManager();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
    //    delete fHistoManager;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run *) {}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run *)
{
    //

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
