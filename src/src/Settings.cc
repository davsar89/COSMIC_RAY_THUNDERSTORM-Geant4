#include "G4UnitsTable.hh"
#include "Settings.hh"
#include "G4SystemOfUnits.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Settings
{
    // All these variables are put here to be shared amounts source file
    // (not very c++ but easier to implement)
    //
    G4long NB_EVENT = 0; // initialisation is important here

    const double ENERGY_MIN = 300. * keV;
    const double ENERGY_MAX = 100. * MeV;

    // variables storing simulation input parameters
    G4long   RANDOM_SEED        = 250; // 22698
    const G4double ALTITUDE_CR_PROTON_SAMPLE = 199.*km; // meters, this value is just for initialization

    G4double POTENTIAL_VALUE = 250.0 ; // MV, may be override by input argument

    G4double EFIELD_REGION_LEN = 2.0 ; // km
    G4double EFIELD_REGION_ALT_CENTER = 12.0 ; // km

    std::vector<G4double> RECORD_ALTITUDES; // km

    G4double TILT = 0.0 * degree;

    const G4double CLOUD_TOP_ALTITUDE = 13.; //km, only used for ABOVE scenario

    const Initial_State initial_state = parma;
    Efield_State current_efield_status = BELOW_PLANE;
    Efield_State initial_efield_status = BELOW_PLANE;

    //    const G4double longitude = 130.5; // ILDAS positron event coordinates (Australia)
    //    const G4double latitude = -13.5;  // ILDAS positron event coordinates (Australia)
    const G4double drOverR = 0.001;

    const G4double longitude = -103.5; // FEGS glow coordinates (Colorado)
    const G4double latitude = 39.5;    // FEGS glow coordinates (Colorado)

    const bool USE_MAX_STEP_FOR_EFIELD = false;
    const G4double MAX_STEP_INSIDE_EFIELD = 0.2; // cm, only for electron and positrons

    const bool USE_MAX_STEP_FOR_RECORD = false; // if turned on, may cause a bug that blocks the code at some event

    const bool CR_GENRATOR_write_output_FOR_TEST = false; // test or not the CR particle generator (i.e. output the list of particles)
    const bool ATMOS_LAYERS_OUTPUT_TO_FILE = false;
    const bool WRITE_MOM_OUTPUT_FOR_TEST = false;

    bool USE_STACKING_ACTION = false;  // will mimic time oriented simulation if set to true
    // can be a bad idea to set it on, because it can use a lot of memory compared to default G4 behaviour

    G4double DELTA_T = -77.88; // microsecond , value just for initialization

    const bool USE_WALL_TIME_LIMIT_FOR_EVENT = false;

    double wall_T_begin_event = 0;

    double VARIABLE_TIME_LIMIT = 0;

    G4double SIZE_RECORD_LAYER = 10 * km;
    //
    int RREA_PART_NB_LIMIT_HAS_BEEN_REACHED = 0; // a flag indication if the limit of particles for RREA has been reached
    // meaning that real multiplciation factor is higher than the one obtained

    const G4double CR_SAMPLING_XY_HALF_SIZE = 50.0; // km
    const G4double EFIELD_XY_HALF_SIZE = 8.0; // km

    const bool TIME_EVENT_DURATIONS = false;

}
