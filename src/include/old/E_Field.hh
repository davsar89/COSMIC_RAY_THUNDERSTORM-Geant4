////////////////////////////////////////////////////////////////////////////////

// /* GEANT4 code for propagation of gamma-rays, electron and positrons in Earth's atmosphere */
//
// //
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
#pragma once

#include <vector>

#include "globals.hh"

#include "G4ElectroMagneticField.hh"
#include "Settings.hh"

class Custom_E_Field: public G4ElectroMagneticField
{
    public:

        Custom_E_Field();
        ~Custom_E_Field();

        /// DoesFieldChangeEnergy() returns true.
        virtual G4bool DoesFieldChangeEnergy() const
        {
            return true;
        }

        void GetFieldValue(const G4double Point[4],
                           G4double *field) const;

    private:

        G4double EFIELD0_AMBIENT = 0.0e5 * volt / meter;  // value for initialization

        G4double EField_sign = -1;


        const G4double density_ground = 1.340E-03;
        const G4double density_20km = 8.573E-05 ;
        const G4double density_10km =  4.128E-04;
        const G4double density_11km =  3.550E-04;
        const G4double density_12km =  3.035E-04;
        const G4double density_13km =  2.587E-04;

        const G4double RREA_thres_STP = 2.84e5; // V/m

        const G4double positive_layer_altitude = Settings::EFIELD_REGION_ALT_CENTER *km + Settings::EFIELD_REGION_LEN *km / 2.;
        const G4double negative_layer_altitude = Settings::EFIELD_REGION_ALT_CENTER *km - Settings::EFIELD_REGION_LEN *km / 2.;

};
