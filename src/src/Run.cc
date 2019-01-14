#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"


#include "G4ProcessTable.hh"
#include "G4RadioactiveDecay.hh"
#include "G4TwoVector.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction *det)
    : G4Run(), fDetector(det), fParticle(0)
{
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run() {}


// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{

} // Run::EndOfRun

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
