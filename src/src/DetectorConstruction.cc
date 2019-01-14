// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include <iostream>
#include <fstream>

#include <G4GDMLParser.hh>
#include <string>
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4UnionSolid.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4SubtractionSolid.hh"

#include "G4EllipticalTube.hh"

// #include <CADMesh.hh>

#include <G4Sphere.hh>
#include "RegionInformation.hh"

#include "G4RegionStore.hh"
#include <cmath>

using namespace std;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//C     INPUT VARIABLES:
//C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
//C              (Year ignored in current model)
//C        SEC - UT(SEC)
//C        ALT - ALTITUDE(KM)
//C        GLAT - GEODETIC LATITUDE(DEG)
//C        GLONG - GEODETIC LONGITUDE(DEG)
//C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
//C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
//C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
//C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
//C           - ARRAY CONTAINING:
//C             (1) DAILY AP
//C             (2) 3 HR AP INDEX FOR CURRENT TIME
//C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
//C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
//C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
//C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
//C                    TO CURRENT TIME
//C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
//C                    TO CURRENT TIME
//C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
//C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
//C                 MASS 17 IS Anomalous O ONLY.)
//C
//C     NOTES ON INPUT VARIABLES:
//C        UT, Local Time, and Longitude are used independently in the
//C        model and are not of equal importance for every situation.
//C        For the most physically realistic calculation these three
//C        variables should be consistent (STL=SEC/3600+GLONG/15).
//C        The Equation of Time departures from the above formula
//C        for apparent local time can be included if available but
//C        are of minor importance.
//c
//C        F107 and F107A values used to generate the model correspond
//C        to the 10.7 cm radio flux at the actual distance of the Earth
//C        from the Sun rather than the radio flux at 1 AU. The following
//C        site provides both classes of values:
//C        ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
//C
//C        F107, F107A, and AP effects are neither large nor well
//C        established below 80 km and these parameters should be set to
//C        150., 150., and 4. respectively.
//C
//C     OUTPUT VARIABLES:
//C        D(1) - HE NUMBER DENSITY(CM-3)
//C        D(2) - O NUMBER DENSITY(CM-3)
//C        D(3) - N2 NUMBER DENSITY(CM-3)
//C        D(4) - O2 NUMBER DENSITY(CM-3)
//C        D(5) - AR NUMBER DENSITY(CM-3)
//C        D(6) - TOTAL MASS DENSITY(GM/CM3)
//C        D(7) - H NUMBER DENSITY(CM-3)
//C        D(8) - N NUMBER DENSITY(CM-3)
//C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
//C        T(1) - EXOSPHERIC TEMPERATURE
//C        T(2) - TEMPERATURE AT ALT

// IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T
#include <fortran.hh>

// extrernal fortran subroutine to get MSIS atmospheric densities
extern "C" {
    void gtd7_(INTEGER &IYD, // YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
               REAL &SEC, // UT(SEC)
               REAL &ALT, // ALTITUDE(KM)
               REAL &GLAT, // GEODETIC LATITUDE(DEG)
               REAL &GLONG, // GEODETIC LONGITUDE(DEG)
               REAL &STL, // LOCAL APPARENT SOLAR TIME
               REAL &F107A, // 81 day AVERAGE OF F10.7 FLUX (centered on day DDD
               REAL &F107, // DAILY F10.7 FLUX FOR PREVIOUS DAY
               REAL &AP,  // MAGNETIC INDEX(DAILY)
               INTEGER &MASS, // MASS NUMBER
               REAL *D, REAL *T); // OUTPUT VARIABLES temperatures
}

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction()
{
    sWorld = 0;
    lWorld = 0;
    pWorld = 0;

    if (Settings::ATMOS_LAYERS_OUTPUT_TO_FILE)
        {
            asciiFile.open("./tests/alt_dens.txt", std::ios::trunc);

            asciiFile << "ID // altitude (km) // range (km) // density (g/cm2) // has_efield // is_record" << G4endl;
        }

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    if (Settings::ATMOS_LAYERS_OUTPUT_TO_FILE) asciiFile.close();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{


    // Cleanup old geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    //
    // World Material : Vaccum
    G4NistManager *man = G4NistManager::Instance();
    G4Material *vac = man->FindOrBuildMaterial("G4_Galactic");

    ////////////////////////////////////////
    /// // assigning default region information
    RegionInformation *Reginfo_def = new RegionInformation();
    Reginfo_def->Set_World();
    G4Region *defaultRegion = (*(G4RegionStore::GetInstance()))[0]; // the default (world) region is index 0 in the region store
    defaultRegion->SetUserInformation(Reginfo_def);

    ////////////////////////////////////////
    // World Definition
    //
    sWorld = new G4Box("sWorld", World_XZ_size * 1.01, World_Y_size * 1.01, World_XZ_size * 1.01);
    lWorld = new G4LogicalVolume(sWorld,        // shape
                                 vac,      // material
                                 "lWorld");     // name

    pWorld = new G4PVPlacement(0,               // no rotation
                               G4ThreeVector(), // at (0,0,0)
                               lWorld,          // logical volume
                               "pWorld",        // name
                               0,               // mother volume
                               false,           // no boolean operation
                               0);              // copy number

    G4VisAttributes *VisAtt_vac = new G4VisAttributes(G4Colour(0.0, 0.8, 0.0, 0.2));
    lWorld->SetVisAttributes(VisAtt_vac);
    //    VisAtt_vac->SetVisibility(false);


    create_thunderstorm_electric_fields();
    lWorld->SetFieldManager(globalfieldMgr, false);


    build_air_layers();

    // always return the root volume
    //
    return pWorld;
} // DetectorConstruction::Construct


// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Air layers
void DetectorConstruction::build_air_layers()
{
    sens_det_List.clear();

    G4FieldManager *null_field = nullptr;

    G4double EF_val = 0.0;

    std::vector < G4Material * > Airs;

    altitudes = calculate_altitudes_list(0.1 * km,  World_Y_size, nb_altitudes);
    Airs = Construct_Atmos_layers_Materials_simple(altitudes);

    std::vector<G4VSolid *> solid_boxes;
    std::vector<G4LogicalVolume *>   logic_boxes;
    std::vector<G4VPhysicalVolume *> phys_boxes;

    G4VisAttributes *visAtt_reg = new G4VisAttributes(G4Colour(0.8, 0.2, 0.8, 0.2));
    visAtt_reg->SetVisibility(true);

    G4VisAttributes *visAtt_reg2 = new G4VisAttributes(G4Colour(0.8, 0.3, 0.8, 0.2));
    visAtt_reg2->SetVisibility(true);

    G4VisAttributes *visAtt_eField_region = new G4VisAttributes(G4Colour(0.49, 0.98, 0., 0.2));
    visAtt_eField_region->SetVisibility(true);

    G4VisAttributes *visAtt_record_region = new G4VisAttributes(G4Colour(0.78, 0.29, 0., 0.2));
    visAtt_eField_region->SetVisibility(true);

    // building altitude layers

    G4double alt_Min_Efield = (Settings::EFIELD_REGION_ALT_CENTER - Settings::EFIELD_REGION_LEN / 2.) * km;
    G4double alt_Max_Efield = (Settings::EFIELD_REGION_ALT_CENTER + Settings::EFIELD_REGION_LEN / 2.) * km;

    G4int id_SD = 0;

    for (unsigned int idx_alt = 0; idx_alt < altitudes.size() - 1; ++idx_alt)
        {
            const G4double delta_alt = (altitudes[idx_alt + 1] - altitudes[idx_alt]);
            const G4double mean_alt = (altitudes[idx_alt + 1] + altitudes[idx_alt]) * G4double(0.50000000);

            const G4double ymin = altitudes[idx_alt];
            const G4double ymax = altitudes[idx_alt + 1];

            const G4double position_Y = mean_alt;

            G4bool has_efield = false;
            G4bool has_record = false;

            //            G4cout << Airs[index_air]->GetDensity() / (g / cm3) << G4endl ;

            solid_boxes.push_back(new G4Box("s_layer_" + std::to_string(idx_alt), World_XZ_size, delta_alt * G4double(0.50000000), World_XZ_size));

            logic_boxes.push_back(new G4LogicalVolume(solid_boxes.back(), Airs[idx_alt], "l_layer_" + std::to_string(idx_alt)));

            const G4double densityval = Airs[idx_alt]->GetDensity() / (g / cm3);

            // small change of color to better visualize the layers
            if (idx_alt % 2 == 0) logic_boxes.back()->SetVisAttributes(visAtt_reg);
            else logic_boxes.back()->SetVisAttributes(visAtt_reg2);


            // setting sensitive dector in record layers

            for (uint idx_r_alt = 0; idx_r_alt < Settings::RECORD_ALTITUDES.size(); ++idx_r_alt) // loop over all requested record altitudes
                {
                    G4double alt_record = (Settings::RECORD_ALTITUDES[idx_r_alt]) * km;

                    if (ymin == alt_record)
                        {
                            if (Settings::USE_MAX_STEP_FOR_RECORD) logic_boxes.back()->SetUserLimits(stepLimit_record); // maybe be overridden if this volume has also an electric field

                            logic_boxes.back()->SetVisAttributes(visAtt_record_region);
                            has_record = true;

                            sens_det_List.push_back(new SensitiveDet("sens_det_" + std::to_string(id_SD), id_SD, alt_record / km));
                            id_SD++;
                            logic_boxes.back()->SetSensitiveDetector(sens_det_List.back());
                            G4SDManager::GetSDMpointer()->AddNewDetector(sens_det_List.back());
                            Settings::SIZE_RECORD_LAYER = std::min(Settings::SIZE_RECORD_LAYER, delta_alt * G4double(0.50000000));
                        }
                }

            if (Settings::current_efield_status == Settings::BELOW_PLANE)
                {
                    // setting E-field in some regions
                    if ((ymin >= alt_Min_Efield && ymax <= alt_Max_Efield))
                        {
                            logic_boxes.back()->SetRegion(EFIELD_Region);
                            EFIELD_Region->AddRootLogicalVolume(logic_boxes.back());

                            if (Settings::USE_MAX_STEP_FOR_EFIELD) logic_boxes.back()->SetUserLimits(stepLimit_efield);

                            logic_boxes.back()->SetFieldManager(localfieldMgr, true);

                            logic_boxes.back()->SetVisAttributes(visAtt_eField_region);

                            has_efield = true;
                            EF_val = (Settings::POTENTIAL_VALUE * megavolt) / (Settings::EFIELD_REGION_LEN * km);
                        }
                    else
                        {
                            logic_boxes.back()->SetFieldManager(null_field, true);
                            has_efield = false;
                            EF_val = 0.0;
                        }
                }

            G4ThreeVector position_box {0, position_Y, 0};

            phys_boxes.push_back(new G4PVPlacement(rotation_null, position_box, logic_boxes.back(),
                                                   "atmosLayer_" + std::to_string(idx_alt) + "_" + std::to_string(mean_alt / km),
                                                   lWorld, false, 0, 0));

            if (Settings::ATMOS_LAYERS_OUTPUT_TO_FILE)
                {
                    asciiFile <<  idx_alt << " " << ymin / km << " " << ymax / km << " " << delta_alt *double(0.50) / km << " " << densityval << " " << has_efield << " " << EF_val << " " << has_record << G4endl;
                }
        }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::create_thunderstorm_electric_fields()
{
    G4double MinStep = 1. * micrometer; // minimal step, micrometer

    G4double minEps = pow(10., -5); //   Minimum & value for smallest steps
    G4double maxEps = pow(10., -5);

    // set up global field (that is zero / null) and sub volumes should overrride it
    //    G4ElectroMagneticField *myEfield = new Custom_E_Field;
    //    G4ElectroMagneticField *myEfield = new Custom_E_Field;

    G4ElectroMagneticField *myEfield = new G4UniformElectricField(G4ThreeVector(0., 0., 0.));
    G4EqMagElectricField *Equation = new G4EqMagElectricField(myEfield);
    G4int nvar = 8; // to integrate time and energy
    G4MagIntegratorStepper *EStepper = new G4DormandPrince745(Equation, nvar);
    // Get the global field manager
    globalfieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    // Set this field to the global field manager
    globalfieldMgr->SetDetectorField(myEfield);
    G4MagInt_Driver *EIntgrDriver = new G4MagInt_Driver(MinStep, EStepper, EStepper->GetNumberOfVariables());
    G4ChordFinder *EChordFinder = new G4ChordFinder(EIntgrDriver);
    globalfieldMgr->SetChordFinder(EChordFinder);
    globalfieldMgr->SetMinimumEpsilonStep(minEps);
    globalfieldMgr->SetMaximumEpsilonStep(maxEps);
    //    globalfieldMgr->SetDeltaOneStep(0.5e-3 * mm);


    // set up local field (non null)
    G4double EFIELD_VAL = (Settings::POTENTIAL_VALUE * megavolt) / (Settings::EFIELD_REGION_LEN * km) ;

    // NOT DONE HERE ANYMORE : check if above RREA threshold to see if it will use Stacking Action or not
    // should not be done here, as this will be excecuted after stacking action is initilized.
    //    G4double EFIELD_VAL_KV_PER_M = std::abs(EFIELD_VAL / (kilovolt / meter));
    //    G4double atmos_scale = get_scale(Settings::EFIELD_REGION_ALT_CENTER); // should be > 1

    //    Settings::DELTA_T = 27.3 / (EFIELD_VAL_KV_PER_M * atmos_scale - 277.0) * atmos_scale; // in micro second for the stacking action that set up the fake time oriented
    //    double DELTA_T_tmp = Settings::DELTA_T ; // variable for debug mode check

    //    if (atmos_scale < 1.0)
    //        {
    //            G4cout << "ERROR in create_thunderstorm_electric_fields : atmos_scale is <1. Aborting" << G4endl;
    //            std::abort();
    //        }

    //    // if below RREA threshold, no need to have the stacking action for fake time oriented
    //    if (Settings::DELTA_T < 0.0)
    //        {
    //            Settings::USE_STACKING_ACTION = false;
    //        }
    //    else
    //        {
    //            Settings::USE_STACKING_ACTION = true;
    //        }

    /////////////////////////////////////

    G4ElectroMagneticField *myEfield_local = new G4UniformElectricField_timeCut(G4ThreeVector(-EFIELD_VAL * std::sin(Settings::TILT), -EFIELD_VAL * std::cos(Settings::TILT), 0.));
    G4EqMagElectricField *Equation_local = new G4EqMagElectricField(myEfield_local);
    G4MagIntegratorStepper *EStepper_local = new G4DormandPrince745(Equation_local, nvar);
    localfieldMgr = new G4FieldManager(myEfield_local);
    G4MagInt_Driver *EIntgrDriver_local = new G4MagInt_Driver(MinStep, EStepper_local, EStepper_local->GetNumberOfVariables());
    G4ChordFinder *EChordFinder_local = new G4ChordFinder(EIntgrDriver_local);
    localfieldMgr->SetChordFinder(EChordFinder_local);
    localfieldMgr->SetMinimumEpsilonStep(minEps);
    localfieldMgr->SetMaximumEpsilonStep(maxEps);
    //    localfieldMgr->SetDeltaOneStep(0.5e-3 * mm);


    // NOT USED ANYMORE, USE OF STEP MAX IS THE PHYSICS LIST (and by region) IS PREFERED
    //    G4double MAX_STEP = 1 * meter;
    G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetLargestAcceptableStep(1.0 * meter);
}

std::vector<G4double> DetectorConstruction::calculate_altitudes_list(G4double alt_min, G4double alt_max_construction, G4int nb_altitudes_)
// fills the vector altitudes
{
    // defining the altitude vector
    for (G4int idx_alt = 0; idx_alt < nb_altitudes_; idx_alt++)   // geocentric altitudes
        {
            //            G4double small_rand_variation = 1.0 + (G4UniformRand() - 0.5) / 500.0;

            altitudes.push_back(
                std::exp(std::log(alt_min) + (std::log(alt_max_construction) - std::log(alt_min)) * G4double(idx_alt) / G4double(nb_altitudes_ - 1))
            );

            //            G4cout << altitudes.back() << G4endl;

            if (std::isnan(altitudes.back()))
                {
                    G4cout << "Error, one altitude is NaN" << G4endl;
                    std::abort();
                }
        }


    //    // defining the altitude vector
    //    for (G4int idx_alt = 0; idx_alt < nb_altitudes_; idx_alt++)   // geocentric altitudes
    //        {
    //            altitudes.push_back(
    //                (alt_min) + ((alt_max_construction) - (alt_min)) * double(idx_alt) / double(nb_altitudes_ - 1)
    //            );

    //            //            G4cout << altitudes.back() << G4endl;

    //            if (std::isnan(altitudes.back()))
    //                {
    //                    G4cout << "Error, one altitude is NaN" << G4endl;
    //                    std::abort();
    //                }
    //        }

    // adding altitudes of E-field region

    G4double alt_Min_Efield = (Settings::EFIELD_REGION_ALT_CENTER - Settings::EFIELD_REGION_LEN / G4double(2.00000000)) * km;
    G4double alt_Max_Efield = (Settings::EFIELD_REGION_ALT_CENTER + Settings::EFIELD_REGION_LEN / G4double(2.00000000)) * km;
    G4double alt_Mid_Efield = Settings::EFIELD_REGION_ALT_CENTER * km;

    if (!contains(altitudes, alt_Mid_Efield)) altitudes.push_back(alt_Mid_Efield);

    if (!contains(altitudes, alt_Min_Efield)) altitudes.push_back(alt_Min_Efield);

    if (!contains(altitudes, alt_Max_Efield)) altitudes.push_back(alt_Max_Efield);

    // adding altitudes of record region (narrow slice)

    for (uint kk = 0; kk < Settings::RECORD_ALTITUDES.size(); ++kk)
        {
            G4double to_add = Settings::RECORD_ALTITUDES[kk] * km;

            if (!contains(altitudes, to_add)) altitudes.push_back(to_add);

            to_add = Settings::RECORD_ALTITUDES[kk] * km + maxStep_record * G4double(10.0000000);

            if (!contains(altitudes, to_add)) altitudes.push_back(to_add);
        }

    if (!contains(altitudes, 13.0 * km)) altitudes.push_back(13.0 * km);

    // check if there are duplicates on the list (that should not happen), if so -> crash

    if (hasDuplicates(altitudes))
        {
            G4cout << "ERROR : There are duplicates values in the altitude list. Aborting." << G4endl;
            std::abort();
        }

    // sorting in increasing value
    std::sort(altitudes.begin(), altitudes.end());

    G4cout << "Successfully created altitude list" << G4endl;

    return altitudes;
}

// check if sorting is good
//for (uint var = 0; var < altitudes.size() - 1; ++var)
//    {
//        if (altitudes[var] >= altitudes[var + 1])
//            {
//                G4cout << "ERROR : altitude lsit is not strictly increasing. Aborting." << G4endl;
//                std::abort();
//            }
//    }

bool DetectorConstruction::hasDuplicates(const std::vector<G4double> &arr)
{
    for (uint i = 0; i < arr.size(); ++i)
        {
            for (uint j = i + 1; j < arr.size(); ++j)
                {
                    if (arr[i] == arr[j])
                        return true;
                }
        }

    return false;
}

// Caltulating materials of the atmopsheric layers, based on the MSIS C++ model integrated to this code
// ref : https://ccmc.gsfc.nasa.gov/modelweb/atmos/nrlmsise00.html
// simplified version using only the total mass density and same ratio of N2 O2 and AR
// Caltulating materials of the atmopsheric layers, based on the MSIS C++ model integrated to this code
// ref : https://ccmc.gsfc.nasa.gov/modelweb/atmos/nrlmsise00.html
// simplified version using only the total mass density and same ratio of N2 O2 and AR
std::vector<G4Material *> DetectorConstruction::Construct_Atmos_layers_Materials_simple(const std::vector<G4double>altitudes_)
{
    std::vector < G4Material * > Airs;

    // Vaccum
    G4NistManager *man = G4NistManager::Instance();
    G4Material *vaccum = man->FindOrBuildMaterial("G4_Galactic");

    //    elHe = new G4Element(name = "Helium", symbol = "He", z = 2., He_molarMass);
    //    elH  = new G4Element(name = "Hydrogen", symbol = "H", z = 1., H_molarMass);

    for (uint idx_alt = 0; idx_alt < altitudes_.size() - 1; idx_alt++)
        {
            const double innerAlt = altitudes_[idx_alt];
            const double outerAlt = altitudes_[idx_alt + 1];
            const double mid_altitude_in_km = (innerAlt + outerAlt) / 2. / km; // geocentric altitude

            //            G4cout << altitude_in_km << G4endl;

            if (mid_altitude_in_km > World_Y_size / km)
                {
                    Airs.push_back(vaccum);
                }
            else
                {
                    INTEGER input_iyd    = 172; // IYD - YEAR AND DAY AS YYDDD
                    REAL input_sec    = 29000.0;
                    REAL input_alt    = (REAL) mid_altitude_in_km;
                    REAL input_g_lat  = (REAL) Settings::latitude;
                    REAL input_g_long = (REAL) Settings::longitude;
                    REAL input_lst    = 16.0;
                    REAL input_f107A  = 150.0;
                    REAL input_f107   = 150.0;
                    REAL input_ap     = 4.0;
                    INTEGER input_mass     = 48;
                    REAL output_D[9];
                    REAL output_T[2];

                    gtd7_(input_iyd, input_sec, input_alt, input_g_lat, input_g_long, input_lst, input_f107A, input_f107, input_ap, input_mass, output_D, output_T); // MSIS, fortran function call

                    if (std::isnan(output_D[5]) || std::isinf(isnan(output_D[5])))
                        {
                            G4cout << "ERROR : density from gtd7_ is NaN. Aborting" << G4endl;
                            std::abort();
                        }

                    G4double density_air = output_D[5] * g / cm3; // getting density and converting it to the GEANT4 system of

                    //                    G4cout << input_alt << " " << output_D[5] << G4endl;

                    G4String name;
                    G4int ncomponents = 2;
                    G4double fractionmass;

                    Airs.push_back(new G4Material(name = "Air_" + std::to_string(idx_alt), density_air, ncomponents = 2));

                    Airs.back()->AddElement(elN, fractionmass = 0.7615);
                    Airs.back()->AddElement(elO, fractionmass = 0.2385);
                }
        }

    G4cout << "Successfully created air materials list" << G4endl;
    return Airs;
}

//======================================================================
// Returns interpolated value at x from parallel arrays ( xData, yData )
//   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
//   boolean argument extrapolate determines behaviour beyond ends of array (if needed)
G4double DetectorConstruction::interpolate(std::vector<G4double> &xData, std::vector<G4double> &yData, G4double x, bool extrapolate)
{
    int size = xData.size();

    int i = 0;                                                                  // find left end of interval for interpolation

    if (x >= xData[size - 2])                                                   // special case: beyond right end
        {
            i = size - 2;
        }
    else
        {
            while (x > xData[i + 1]) i++;
        }

    double xL = xData[i], yL = yData[i], xR = xData[i + 1], yR = yData[i + 1];  // points on either side (unless beyond ends)

    if (!extrapolate)                                                           // if beyond ends of array and not extrapolating
        {
            if (x < xL) yR = yL;

            if (x > xR) yL = yR;
        }

    double dydx = (yR - yL) / (xR - xL);                                        // gradient

    return yL + dydx * (x - xL);                                                // linear interpolation
}

////////////////////////////////////////////////////////////////////////////////////////////////////


G4double DetectorConstruction::contains(std::vector<G4double> v, G4double x)
{
    if (std::find(v.begin(), v.end(), x) != v.end())
        {
            return true; /* v contains x */
        }
    else
        {
            return false; /* v does not contain x */
        }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

G4double DetectorConstruction::get_scale(G4double alt)
{
    // returns the atmospheric relative scale compared to sea level (>1)
    std::vector<G4double> alt_list =
    {
        0,
        0.500000000000000,
        1.000000000000000,
        1.500000000000000,
        2.000000000000000,
        2.500000000000000,
        3.000000000000000,
        3.500000000000000,
        4.000000000000000,
        4.500000000000000,
        5.000000000000000,
        5.500000000000000,
        6.000000000000000,
        6.500000000000000,
        7.000000000000000,
        7.500000000000000,
        8.000000000000000,
        8.500000000000000,
        9.000000000000000,
        9.500000000000000,
        10.000000000000000,
        10.500000000000000,
        11.000000000000000,
        11.500000000000000,
        12.000000000000000,
        12.500000000000000,
        13.000000000000000,
        13.500000000000000,
        14.000000000000000,
        14.500000000000000,
        15.000000000000000,
        15.500000000000000,
        16.000000000000000,
        16.500000000000000,
        17.000000000000000,
        17.500000000000000,
        18.000000000000000,
        18.500000000000000,
        19.000000000000000,
        19.500000000000000,
        20.000000000000000,
    };

    std::vector<G4double> scale_list =
    {
        1.000000000000000,
        1.059301380991064,
        1.121238177128117,
        1.184377838328792,
        1.249042145593870,
        1.317304778260430,
        1.388415672913118,
        1.463688404983724,
        1.543377914546100,
        1.628371628371629,
        1.719862833025587,
        1.818435364663227,
        1.925291598996014,
        2.041647095663066,
        2.168634624979211,
        2.307556184746062,
        2.460377358490566,
        2.628502318081032,
        2.813376483279396,
        3.018518518518519,
        3.245395719263315,
        3.496916063287745,
        3.774240231548481,
        4.080100125156446,
        4.414353419092755,
        4.781811514484781,
        5.180770758839889,
        5.613430908308223,
        6.082089552238807,
        6.589186457806973,
        7.133479212253830,
        7.720544701006513,
        8.348271446862997,
        9.024221453287199,
        9.753178758414361,
        10.541632983023444,
        11.398601398601400,
        12.325141776937620,
        13.334696799263730,
        14.431164231961047,
        15.627996164908916,
    };

    if ((alt > 20.) || (alt < 0.))
        {
            G4cout << "ERROR in get_scale : altitude is not between 0 and 20. Aborting.";
            std::abort();
        }

    return interpolate(alt_list, scale_list, alt, false);

}

// Caltulating materials of the atmopsheric layers, based on the MSIS C++ model integrated to this code
// ref : https://ccmc.gsfc.nasa.gov/modelweb/atmos/nrlmsise00.html

//std::vector<G4Material *> DetectorConstruction::Construct_Atmos_layers_Materials(const std::vector<G4double>altitudes_)
//{

//    std::vector < G4Material * > Airs;

//    // Vaccum
//    G4NistManager *man = G4NistManager::Instance();
//    G4Material *vaccum = man->FindOrBuildMaterial("G4_Galactic");

//    struct nrlmsise_flags flags;

//    for (int jj = 1; jj < 24; jj++)
//        {
//            flags.switches[jj] = 1;
//        }

//    G4double N_molarMass  = 14.0067 * g / mole;
//    G4double O_molarMass  = 15.9994 * g / mole;
//    G4double Ar_molarMass = 39.948 * g / mole;
//    //    G4double He_molarMass = 4.002602 * g / mole;
//    //    G4double H_molarMass  = 1.00794 * g / mole;
//    G4Element  *elA = 0;
//    G4Element  *elN = 0;
//    G4Element  *elO = 0;
//    //    G4Element  *elHe = 0;
//    //    G4Element  *elH = 0;
//    G4String name, symbol;
//    G4double z;
//    G4int    ncomponents, natoms;
//    elN  = new G4Element(name = "Nitrogen", symbol = "N", z = 7., N_molarMass);
//    elO  = new G4Element(name = "Oxygen", symbol = "O", z = 8., O_molarMass);
//    elA  = new G4Element(name = "Argon", symbol = "Ar", z = 18., Ar_molarMass);
//    //    elHe = new G4Element(name = "Helium", symbol = "He", z = 2., He_molarMass);
//    //    elH  = new G4Element(name = "Hydrogen", symbol = "H", z = 1., H_molarMass);

//    for (unsigned int idx_alt = 0; idx_alt < altitudes_.size() - 1; idx_alt++)
//        {
//            const double innerAlt = altitudes_[idx_alt];
//            const double outerAlt = altitudes_[idx_alt + 1];
//            const double altitude_in_km = (innerAlt + outerAlt) / 2. / km; // geocentric altitude

//            //            G4cout << altitude_in_km << G4endl;

//            if (altitude_in_km > World_Y_size / km)
//                {
//                    Airs.push_back(vaccum);
//                }
//            else
//                {
//                    struct geocentric_coords geocentric;
//                    struct geodetic_coords   geodetic;
//                    struct ecef_coords ecef;

//                    // The MSIS routine takes in argument geodetic coordinates
//                    // we know geodetic coordinates of the source (chosen)
//                    // but we use geocentric coordinates to make the atmospheric layers

//                    geodetic.lat = Settings::latitude;
//                    geodetic.lon = Settings::longitude;
//                    geodetic.alt = 15.;

//                    // determining geocentric lat and lon corresponding to geodetic lat and lon
//                    geodetic_to_ecef(&geodetic, &ecef);
//                    ecef_to_geocentric(&ecef, &geocentric);

//                    geocentric.alt = altitude_in_km; // putting in the (geocentric) altitude of the layer

//                    // determining the geodetic altitude corresponding to the geocentric altitude
//                    geocentric_to_ecef(&geocentric, &ecef);
//                    ecef_to_geodetic_olson(&ecef, &geodetic);

//                    //IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T

//                    int input_iyd    = 172; // IYD - YEAR AND DAY AS YYDDD
//                    double input_sec    = 29000.0;
//                    double input_alt    = geodetic.alt;
//                    double input_g_lat  = Settings::latitude;
//                    double input_g_long = Settings::longitude;
//                    double input_lst    = 16.0;
//                    double input_f107A  = 150.0;
//                    double input_f107   = 150.0;
//                    double input_ap     = 4.0;
//                    double input_mass     = 48.0;
//                    double output_D[9];
//                    double output_T[2];

//                    for (int ii = 0; ii < 9; ++ii)
//                        {
//                            output_D[ii] = 0.0;
//                        }

//                    for (int ii = 0; ii < 2; ++ii)
//                        {
//                            output_T[ii] = 0.0;
//                        }

//                    // G4cout << altitude << G4endl;

//                    gtd7_(&input_iyd, &input_sec, &input_alt, &input_g_lat, &input_g_long, &input_lst, &input_f107A, &input_f107, &input_ap, &input_mass, &output_D, &output_T); // MSIS, C function call

//                    G4double density_air = output_D[5] * g / cm3; // getting density and converting it to the GEANT4 system of unit

//                    // look how many elements in Air with non-null densities and compute total and proportions

//                    G4double nb_density_total = 0.;

//                    G4int ncomponents2 = 0;

//                    //                    if (output.d[0] > 1.e-24)
//                    //                        {
//                    //                            nb_density_total = nb_density_total + output.d[0];
//                    //                            ncomponents2++;
//                    //                        }

//                    if (output.d[1] > 1.e-24)
//                        {
//                            nb_density_total = nb_density_total + output.d[1];
//                            ncomponents2++;
//                        }

//                    if (output.d[2] > 1.e-24)
//                        {
//                            nb_density_total = nb_density_total + output.d[2];
//                            ncomponents2++;
//                        }

//                    if (output.d[3] > 1.e-24)
//                        {
//                            nb_density_total = nb_density_total + output.d[3];
//                            ncomponents2++;
//                        }

//                    if (output.d[4] > 1.e-24)
//                        {
//                            nb_density_total = nb_density_total + output.d[4];
//                            ncomponents2++;
//                        }

//                    //                    if (output.d[6] > 1.e-24)
//                    //                        {
//                    //                            nb_density_total = nb_density_total + output.d[6];
//                    //                            ncomponents2++;
//                    //                        }

//                    if (output.d[7] > 1.e-24)
//                        {
//                            nb_density_total = nb_density_total + output.d[7];
//                            ncomponents2++;
//                        }

//                    G4double proporN2 = output.d[2] / nb_density_total;
//                    G4double proporO2 = output.d[3] / nb_density_total;
//                    G4double proporAr = output.d[4] / nb_density_total;
//                    G4double proporO  = output.d[1] / nb_density_total;
//                    //                    G4double proporHe = output.d[0] / nb_density_total;
//                    //                    G4double proporH  = output.d[6] / nb_density_total;
//                    G4double proporN  = output.d[7] / nb_density_total;

//                    G4double densityN2 = proporN2 * density_air;
//                    G4double densityO2 = proporO2 * density_air;
//                    G4double densityAr = proporAr * density_air;
//                    G4double densityO  = proporO * density_air;
//                    //                    G4double densityHe = proporHe * density_air;
//                    //                    G4double densityH  = proporH * density;
//                    G4double densityN  = proporN * density_air;

//                    // to remove warning messages
//                    // the vales modified by this if statement will not be used
//                    if (densityN2 <= 1.e-24) densityN2 = 1.;

//                    if (densityO2 <= 1.e-24) densityO2 = 1.;

//                    if (densityAr <= 1.e-24) densityAr = 1.;

//                    if (densityO <= 1.e-24) densityO = 1.;

//                    //                    if (densityHe <= 1.e-24) densityHe = 1.;

//                    //                    if (densityH <= 1.e-24) densityH = 1.;

//                    if (densityN <= 1.e-24) densityN = 1.;

//                    G4Material *N2 = 0;
//                    G4Material *O2 = 0;
//                    G4Material *Ar = 0;
//                    G4Material *O = 0;
//                    G4Material *N = 0;
//                    //                    G4Material *H = 0;
//                    //                    G4Material *He = 0;

//                    N2 = new G4Material(name = "N2_" + std::to_string(idx_alt), densityN2, ncomponents = 1);
//                    N2->AddElement(elN, natoms = 2);

//                    O2 = new G4Material(name = "O2_" + std::to_string(idx_alt), densityO2, ncomponents = 1);
//                    O2->AddElement(elO, natoms = 2);

//                    Ar = new G4Material(name = "Ar_" + std::to_string(idx_alt), densityAr, ncomponents = 1);
//                    Ar->AddElement(elA, natoms = 1);

//                    O = new G4Material(name = "O_" + std::to_string(idx_alt), densityO, ncomponents = 1);
//                    O->AddElement(elO, natoms = 1);

//                    //                    He = new G4Material(name = "He_" + std::to_string(idx_alt), densityHe, ncomponents = 1);
//                    //                    He->AddElement(elHe, natoms = 1);

//                    //                    H = new G4Material(name = "H_" + std::to_string(idx_alt), densityH, ncomponents = 1);
//                    //                    H->AddElement(elH, natoms = 1);

//                    N = new G4Material(name = "N_" + std::to_string(idx_alt), densityN, ncomponents = 1);
//                    N->AddElement(elN, natoms = 1);

//                    Airs.push_back(new G4Material(name = "air_" + std::to_string(idx_alt), density_air, ncomponents2));

//                    //                    if (output.d[0] > 1.e-24)
//                    //                        {
//                    //                            Airs.back()->AddMaterial(He, proporHe);
//                    //                        }

//                    if (output.d[1] > 1.e-24)
//                        {
//                            Airs.back()->AddMaterial(O, proporO);
//                        }

//                    if (output.d[2] > 1.e-24)
//                        {
//                            Airs.back()->AddMaterial(N2, proporN2);
//                        }

//                    if (output.d[3] > 1.e-24)
//                        {
//                            Airs.back()->AddMaterial(O2, proporO2);
//                        }

//                    if (output.d[4] > 1.e-24)
//                        {
//                            Airs.back()->AddMaterial(Ar, proporAr);
//                        }

//                    //                    if (output.d[6] > 1.e-24)
//                    //                        {
//                    //                            Airs.back()->AddMaterial(H, proporH);
//                    //                        }

//                    if (output.d[7] > 1.e-24)
//                        {
//                            Airs.back()->AddMaterial(N, proporN);
//                        }
//                }
//        }

//    G4cout << "Successfully created air materials list" << G4endl;
//    return Airs;
//} // TGFDetectorConstruction::Construct_Atmos_layers_Material

