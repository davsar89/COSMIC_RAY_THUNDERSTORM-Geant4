////////////////////////////////////////////////////////////////////////////////

// // ********************************************************************
// // * License and Disclaimer                                           *
// // *                                                                  *
// // * The  Geant4 software  is  copyright of the Copyright Holders  of *
// // * the Geant4 Collaboration.  It is provided  under  the terms  and *
// // * conditions of the Geant4 Software License,  included in the file *
// // * LICENSE and available at  http://cern.ch/geant4/license .  These *
// // * include a list of copyright holders.                             *
// // *                                                                  *
// // * Neither the authors of this software system, nor their employing *
// // * institutes,nor the agencies providing financial support for this *
// // * work  make  any representation or  warranty, express or implied, *
// // * regarding  this  software system or assume any liability for its *
// // * use.  Please see the license in the file  LICENSE  and URL above *
// // * for the full disclaimer and the limitation of liability.         *
// // *                                                                  *
// // * This  code  implementation is the result of  the  scientific and *
// // * technical work of the GEANT4 collaboration.                      *
// // * By using,  copying,  modifying or  distributing the software (or *
// // * any work based  on the software)  you  agree  to acknowledge its *
// // * use  in  resulting  scientific  publications,  and indicate your *
// // * acceptance of all terms of the Geant4 Software license.          *
// // ********************************************************************
////////////////////////////////////////////////////////////////////////////////

#include "G4SystemOfUnits.hh"
#include "E_Field.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Custom_E_Field::Custom_E_Field()
{
    EFIELD0_AMBIENT = Settings::POTENTIAL_VALUE * megavolt / (Settings::EFIELD_REGION_LEN * km);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Custom_E_Field::~Custom_E_Field() {}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// global E-field is made with a component from the plane (1/r decrease) and a component with a constant E-field

void Custom_E_Field::GetFieldValue(const G4double Point[4], G4double *field) const
{

    // field is really field[6]: Bx,By,Bz,Ex,Ey,Ez
    // Point[0],Point[1],Point[2] are x-, y-, z-cordinates, Point[3] is time

    field[0] = field[1] = field[2] = field[3] = field[4] = field[5] = 0.0;

    // protect against Geant4 bug that calls us with point[] NaN.
    if (Point[0] != Point[0]) return;


    if (Point[1] >= negative_layer_altitude && Point[1] <= positive_layer_altitude)
        {
            field[4] = EField_sign * EFIELD0_AMBIENT;
        }

} // EarthMagField::GetFieldValue
