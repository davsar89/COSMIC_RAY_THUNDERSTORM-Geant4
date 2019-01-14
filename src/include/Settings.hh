#pragma once
#include <vector>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Settings
{
    enum Initial_State { parma, CR_proton };
    enum Efield_State { ABOVE_PLANE, BELOW_PLANE, OFF };
    // All these variables are put here to be shared amounts source files
    // (not very c++ but easier to implement)

    // simulation records
    extern G4long NB_EVENT;

    extern const double ENERGY_MIN;
    extern const double ENERGY_MAX;

    // variables storing simulation input parameters
    extern G4long   RANDOM_SEED;
    extern const G4double ALTITUDE_CR_PROTON_SAMPLE; // meters

    extern G4double POTENTIAL_VALUE;

    extern std::vector<G4double> RECORD_ALTITUDES;

    extern G4double EFIELD_REGION_LEN;
    extern G4double EFIELD_REGION_ALT_CENTER;

    extern const G4double CLOUD_TOP_ALTITUDE;

    extern const Initial_State initial_state;
    extern Efield_State current_efield_status;
    extern Efield_State initial_efield_status;

    extern const G4double longitude;
    extern const G4double latitude;
    extern const G4double drOverR;

    extern G4double TILT;

    extern const bool USE_MAX_STEP_FOR_EFIELD;
    extern const G4double MAX_STEP_INSIDE_EFIELD;

    extern const bool USE_MAX_STEP_FOR_RECORD;

    extern const bool CR_GENRATOR_write_output_FOR_TEST;
    extern const bool ATMOS_LAYERS_OUTPUT_TO_FILE;
    extern const bool WRITE_MOM_OUTPUT_FOR_TEST;

    extern bool USE_STACKING_ACTION;
    extern G4double DELTA_T;

    extern const bool USE_WALL_TIME_LIMIT_FOR_EVENT;

    extern double wall_T_begin_event;

    extern double VARIABLE_TIME_LIMIT;

    extern G4double SIZE_RECORD_LAYER;

    extern const G4double CR_SAMPLING_XY_HALF_SIZE;
    extern const G4double EFIELD_XY_HALF_SIZE ;

    extern int RREA_PART_NB_LIMIT_HAS_BEEN_REACHED;

    extern const bool TIME_EVENT_DURATIONS;
}
