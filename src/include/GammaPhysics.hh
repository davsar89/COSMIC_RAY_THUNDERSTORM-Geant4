// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#pragma once

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class GammaPhysics: public G4VPhysicsConstructor
{
    public:

        GammaPhysics(const G4String &name = "gamma");
        ~GammaPhysics();

    public:

        virtual void
        ConstructParticle() {}

        virtual void
        ConstructProcess();
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
