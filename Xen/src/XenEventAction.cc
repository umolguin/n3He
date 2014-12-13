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
/// \file XenEventAction.cc
/// \brief Implementation of the XenEventAction class

#include "XenEventAction.hh"

#include "XenRunAction.hh"
#include "XenSteppingAction.hh"
#include "XenAnalysis.hh"
  // use of stepping action to get and reset accumulated energy

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4RadioactiveDecay.hh"
#include "TFile.h"
#include "TH1F.h"

#include <cmath>

///@author Carlos Olguin
///@version 0.8

G4double XenEventAction::originZ=0;
G4double XenEventAction::originX=0;
G4double XenEventAction::originY=0;
G4double XenEventAction::accumDelta=0;
G4int XenEventAction::lastCellID=0;

G4double XenEventAction::lastTrackID=0;
G4ThreeVector XenEventAction::initialMomentum;
double XenEventAction::lastZMomentum=0;

std::vector<int> XenEventAction::al28IDs;

XenEventAction* XenEventAction::Instance()
{

}

XenEventAction::XenEventAction()
: G4UserEventAction()
{

}



XenEventAction::~XenEventAction()
{


}


void XenEventAction::BeginOfEventAction(const G4Event* event)
{
    G4int eventNb = event->GetEventID();

    G4cout << "\n---> This is an event with ID: " << eventNb << G4endl;

    al28IDs.clear();
    originX=0;
    originY=0;
    originX=0;
    lastTrackID=0;
    accumDelta=0;
    lastCellID=-1;
    lastZMomentum=0;

//    G4RadioactiveDecay *theRadioactiveDecay = new G4RadioactiveDecay();
//
//   G4ProcessManager* pmanager =
//   G4ParticleTable::GetParticleTable()->FindParticle("GenericIon")->GetProcessManager();

//   pmanager ->AddProcess(theRadioactiveDecay);
//   pmanager ->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
//   pmanager ->SetProcessOrdering(theRadioactiveDecay, idxAtRest);



}
void XenEventAction::EndOfEventAction(const G4Event* event)
{
    //re-initialize
    originX=0;
    originY=0;
    originX=0;
    accumDelta=0;
    lastTrackID=0;
    lastZMomentum=0;


}
