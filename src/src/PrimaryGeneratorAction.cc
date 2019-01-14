// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fParticleGun(nullptr) {

    cosmic_ray_gene_p = new Cosmic_Ray_Generator_PARMA();

    G4int n_particle = 1;

    fParticleGun = new G4ParticleGun(n_particle);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {

    out_cosmic = cosmic_ray_gene_p->generate_Cosmic_ray();

    fParticleGun->SetParticleEnergy(out_cosmic.energy);
    fParticleGun->SetParticlePosition(out_cosmic.position_ini);
    fParticleGun->SetParticleMomentumDirection(out_cosmic.momentum_ini);
    fParticleGun->SetParticleTime(0.0);
    fParticleGun->SetParticleDefinition(out_cosmic.g4_particle);
    fParticleGun->GeneratePrimaryVertex(anEvent);

//    int PDG = out_cosmic.g4_particle->GetPDGEncoding();

} // PrimaryGeneratorAction::GeneratePrimaries
