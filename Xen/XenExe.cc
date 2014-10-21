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
//
// $Id$
//
///@author Carlos Olguin
///@version 0.8

/// \file XenExe.cc
/// \brief Main program

#include "XenDetectorConstruction.hh"
#include "XenPrimaryGeneratorAction.hh"
#include "XenRunAction.hh"
#include "XenEventAction.hh"
#include "XenSteppingAction.hh"
#include "XenTrackingAction.hh"
#include "EnergyManager.hh"
#include "MomentumManager.hh"

#include "CellManager.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "QBBC.hh"
#include "QGSP_BERT_HP.hh"
#include "G4HadronElasticProcess.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPElastic.hh"
#include "G4ProcessManager.hh"
#include "XenPhysicsList.hh"
#include "G4StepLimiterBuilder.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "Randomize.hh"

#include "G4GDMLParser.hh"
#include "GDMLDetectorConstruction.hh"
using namespace std;

class CellManager;
int main(int argc,char** argv)
{
    // Construct the default run manager
    //
    G4RunManager * runManager = new G4RunManager;

    // Set mandatory initialization classes
    string _inv=("invisible\n");
    if(argc!=1 && _inv.compare(argv[1])!=0){

        runManager->SetUserInitialization(new XenDetectorConstruction(false));
    }
    else{
        //CellManager init
        CellManager::init();

        runManager->SetUserInitialization(new XenDetectorConstruction());}


    // Physics list
    //QGSP_BERT_HP physics list contains reference to the necessary data sets and calculations for cold neutrons
    G4VModularPhysicsList* physicsList= new QGSP_BERT_HP;
    physicsList->RegisterPhysics(new G4StepLimiterBuilder());
    physicsList->SetVerboseLevel(1);

    runManager->SetUserInitialization(physicsList);

    //Energy manager initialization - helps with the Gaussian distribution
    EnergyManager::init();
    MomentumManager::init();
    //EnergyManager::setRangeWaveLength(0,1);
    EnergyManager::setRangeWaveLength(.32,.69);


    // Primary generator action
    runManager->SetUserAction(new XenPrimaryGeneratorAction());

    // Set user action classes
    //
    // Stepping action
    runManager->SetUserAction(new XenSteppingAction());
    runManager->SetUserAction(new XenTrackingAction());

    // Event action
    runManager->SetUserAction(new XenEventAction());

    // Run action
    runManager->SetUserAction(new XenRunAction());



    runManager->Initialize();


#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  if (argc!=1&&_inv.compare(argv[1])==0) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode : define UI session
#ifdef G4UI_USE
    argv[1]="";//TODO: Correct in case in line commmands used.
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac");
#else
    UImanager->ApplyCommand("/control/execute init.mac");
#endif
    ui->SessionStart();
    delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}
