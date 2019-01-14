#pragma once

#include "G4VUserActionInitialization.hh"
#include "StackingAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

class DetectorConstruction;
class G4VSteppingVerbose;

/// Action initialization class.
///

class ActionInitialization: public G4VUserActionInitialization
{
    public:

        ActionInitialization(DetectorConstruction *detector);
        virtual ~ActionInitialization();

        virtual void BuildForMaster() const;
        virtual void Build() const;

    private:

        DetectorConstruction *fDetector;
};
