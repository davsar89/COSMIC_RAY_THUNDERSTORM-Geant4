//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************

// E-field that will be set to zero the number of secondaries electrons with a given energy threshold is too high

#include "G4UniformElectricField_timeCut.hh"
#include "G4PhysicalConstants.hh"

G4UniformElectricField_timeCut::G4UniformElectricField_timeCut(const G4ThreeVector FieldVector)
{
    fFieldComponents[0] = 0.0;
    fFieldComponents[1] = 0.0;
    fFieldComponents[2] = 0.0;
    fFieldComponents[3] = FieldVector.x();
    fFieldComponents[4] = FieldVector.y();
    fFieldComponents[5] = FieldVector.z();
}

G4UniformElectricField_timeCut::G4UniformElectricField_timeCut(G4double vField,
        G4double vTheta,
        G4double vPhi)
{
    if ((vField < 0) || (vTheta < 0) || (vTheta > pi) || (vPhi < 0) || (vPhi > twopi))
        {
            G4Exception("G4UniformElectricField_timeCut::G4UniformElectricField_timeCut()",
                        "GeomField0002", FatalException, "Invalid parameters.");
        }

    fFieldComponents[0] = 0.0;
    fFieldComponents[1] = 0.0;
    fFieldComponents[2] = 0.0;
    fFieldComponents[3] = vField * std::sin(vTheta) * std::cos(vPhi) ;
    fFieldComponents[4] = vField * std::sin(vTheta) * std::sin(vPhi) ;
    fFieldComponents[5] = vField * std::cos(vTheta) ;
}

G4Field *G4UniformElectricField_timeCut::Clone() const
{
    return new G4UniformElectricField_timeCut(G4ThreeVector(fFieldComponents[3],
            fFieldComponents[4],
            fFieldComponents[5]));
}

G4UniformElectricField_timeCut::~G4UniformElectricField_timeCut()
{
}

G4UniformElectricField_timeCut::G4UniformElectricField_timeCut(const G4UniformElectricField_timeCut &p)
    : G4ElectricField(p)
{
    for (G4int i = 0; i < 6; i++)
        {
            fFieldComponents[i] = p.fFieldComponents[i];
        }
}

G4UniformElectricField_timeCut &
G4UniformElectricField_timeCut::operator = (const G4UniformElectricField_timeCut &p)
{
    if (&p == this) return *this;

    G4ElectricField::operator=(p);

    for (G4int i = 0; i < 6; i++)
        {
            fFieldComponents[i] = p.fFieldComponents[i];
        }

    return *this;
}

// ------------------------------------------------------------------------
void G4UniformElectricField_timeCut::GetFieldValue(const G4double pos[4],
        G4double *fieldBandE) const
{

    //     G4double radius_squared = (pos[0]*pos[0]+pos[2]*pos[2])/(km*km);

    if (std::abs(pos[0]) > Settings::EFIELD_XY_HALF_SIZE * km ||
        std::abs(pos[2]) > Settings::EFIELD_XY_HALF_SIZE * km ||
        Settings::current_efield_status == Settings::OFF)
        {
            fieldBandE[0] = 0.0;
            fieldBandE[1] = 0.0;
            fieldBandE[2] = 0.0;
            fieldBandE[3] = 0.0;
            fieldBandE[4] = 0.0;
            fieldBandE[5] = 0.0;
        }
    else
        {
            fieldBandE[0] = 0.0;
            fieldBandE[1] = 0.0;
            fieldBandE[2] = 0.0;
            fieldBandE[3] = fFieldComponents[3] ;
            fieldBandE[4] = fFieldComponents[4] ;
            fieldBandE[5] = fFieldComponents[5] ;
        }
}



