
#ifndef XenPhysicsList_h
#define XenPhysicsList_h 1
// #include "G4VModularPhysicsList.hh"
//#include "HadronPhysicsLHEP_PRECO_HP.hh"
//#include "LHEP_PRECO_HP.hh"
#include "G4FastSimulationManagerProcess.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4VUserPhysicsList.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// class N3HeSimPhysicsList: public G4VModularPhysicsList

class XenPhysicsList: public G4VUserPhysicsList{
    public:
        XenPhysicsList();
        virtual ~XenPhysicsList();



    protected:
        void ConstructParticle();
        void ConstructProcess();
        void ConstructHad();


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

