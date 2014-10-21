#include "XenRunAction.hh"
#include "XenPrimaryGeneratorAction.hh"
#include "XenEventAction.hh"
#include "XenSteppingAction.hh"
#include "XenSteppingAction.hh"
#include "XenAnalysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "TFile.h"
#include "TH1F.h"

#include "CellManager.hh"
#include "EnergyManager.hh"

///@author Carlos Olguin
///@version 0.8

double XenRunAction::_energy=0;
double XenRunAction::_bEnergy=0;

XenRunAction::XenRunAction()
: G4UserRunAction()
{
;
}

XenRunAction::~XenRunAction()
{}

void XenRunAction::BeginOfRunAction(const G4Run* aRun)
{
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    G4String colNam;

    G4AnalysisManager* man =
    G4AnalysisManager::Instance();
    // Open an output file
    man->OpenFile("Histograms");
    man->SetFirstHistoId(1);
    // Create histogram(s)
    man->CreateH1("ionDist","Ionization distribution", 18, -.8*cm, 2.8*cm);//TODO: change to the appropiate size
    man->CreateH1("asym","Asymmetry", 20, -1, 1);
    man->CreateH1("asymDeg","Asymmetry in rad", 50, 0, 3.14);
    man->CreateH1("wlDist","Wave Length Distribution", 100, 0, 1);
    man->CreateH1("xMomentum","Initial Momentum along #hatx", 200, -1, 1);
    man->CreateH1("yMomentum","Initial Momentum along #haty", 200, -1, 1);
    man->CreateH1("zMomentum","Initial Momentum along #hatz", 200, -1, 1);

    man->CreateH1("rnd","Random Distribution", 200, -1, 1);
    man->CreateH1("sinn","Sin(rnd)", 200, -1, 1);

    man->CreateH2("GeometricFactor","Geometric Factor",16,0,16,9,0,9);
    man->CreateH2("dilution","Dilution Factor",16,0,16,9,0,9);

    EnergyManager::printDist();

}


void XenRunAction::EndOfRunAction(const G4Run* aRun)
{
    G4AnalysisManager* man =G4AnalysisManager::Instance();

	for(int i=0; i<9;i++)
			for(int j=0;j<16;j++)
				man->FillH2(1, i,j,CellManager::getGFactor(j,i));

    for(int i=0;i<=16;i++)
            for(int j=0;j<=9;j++){
            	G4double _d= CellManager::getDFactor(j,i);
                man->FillH2(2, i,j,_d);
            }

    G4int nofEvents = aRun->GetNumberOfEvent();
    if (nofEvents == 0) return;

    G4cout << aRun->GetRunID() << G4endl;

    if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
    }
    CellManager::print();
    std::cout<<"RUN___________________________"<<std::endl;

    man->Write();
    man->CloseFile();

    //START Writing to file

	std::ofstream ofs ("wholeDilutionFactor.dat");
	ofs <<"Energy:"<<XenRunAction::_energy<<", BEnergy:"<<XenRunAction::_bEnergy<<std::endl;
	ofs.close();
	//END Writing

}
