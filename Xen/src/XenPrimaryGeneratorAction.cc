#include "XenPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "EnergyManager.hh"
#include "XenAnalysis.hh"

///@author Carlos Olguin
///@version 0.8

XenPrimaryGeneratorAction::XenPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{
    myEventCounter =0;
    fParticleGun= new G4ParticleGun(1);//Number of particles

    // default particle kinematic
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle= particleTable->FindParticle(particleName="neutron");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
//    fParticleGun->SetParticlePolarization(G4ThreeVector(0.,0.,-1.));//This is not part of original code

    //Static energy
    fParticleGun->SetParticleEnergy(5.0e-3*eV);



    G4cout<<"Particle used: "<<particleName<<G4endl;

}

XenPrimaryGeneratorAction::~XenPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void XenPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

    //Changed to Gaussian distribution
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    double _wl=EnergyManager::getEnergy();
    analysisManager->FillH1(4, _wl);
    double _energy=(0.841/(_wl*_wl));
    fParticleGun->SetParticleEnergy((_energy/1000)*eV);
    //fParticleGun->SetParticleEnergy(5.0e-3*eV);
    G4cout<<"[E]Energy:"<<fParticleGun->GetParticleEnergy()<<"meV"<<G4endl;
	/***
	 * This is not part of original source
	 */

  G4double envSizeXY = 0;
  G4double envSizeZ = 0;

  G4LogicalVolume* envLV= G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box* envBox = NULL;
  if ( envLV ) envBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  if ( envBox ) {
    envSizeXY = envBox->GetXHalfLength()*2.;
    envSizeZ = envBox->GetZHalfLength()*2.;
  }
  else  {
    G4cerr << "Envelope volume of box shape not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }

//--------------------------------------------------------------------------
    G4double size = 50;
    G4double x0 = size*(G4UniformRand()-0.5);
    G4double y0 = size*(G4UniformRand()-0.5);
    //G4double z0 = -0.5 * envSizeZ;
    //z0=4.7*cm;//This would be the position just at the beginning of the right window.
    G4double z0=200*cm;
    G4cout<<"[GUN]:("<<x0<<","<<y0<<","<<z0<<")"<<G4endl;
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
