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
// * regarding  this  software system or assume any liaility for its *
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
/// \file XenSteppingAction.cc
/// \brief Implementation of the XenSteppingAction class

#include "XenSteppingAction.hh"
#include "XenEventAction.hh"
#include "XenRunAction.hh"
#include "XenDetectorConstruction.hh"
#include "MomentumManager.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "XenAnalysis.hh"
#include <cmath>

#include "CellManager.hh"
#include "Logger.hh"

///@author Carlos Olguin
///@version 0.8

bool XenSteppingAction::isNewParticle;
bool XenSteppingAction::isFirstProton;

//XenSteppingAction* XenSteppingAction::Instance()
//{
//
//}
XenSteppingAction::XenSteppingAction(
                      XenDetectorConstruction* detectorConstruction,
                      XenEventAction* eventAction)
  : G4UserSteppingAction(),
    fEventAction(eventAction)
{
}


XenSteppingAction::XenSteppingAction()
: G4UserSteppingAction()

{
}


XenSteppingAction::~XenSteppingAction()
{

}
bool XenSteppingAction::_isParentAl28(int parentID)
{
//	//Printing --
//	Logger::log("____________Electron with ParentID:");
//	Logger::log(std::to_string(parentID));
//	Logger::log("_________________GAMMAsIDs:");
//	for (auto _i = XenEventAction::gammaIDs.begin (); _i != XenEventAction::gammaIDs.end (); ++_i) {
//		Logger::log(std::to_string(*_i));
//	}
//	///-- Printing
	bool _isON=false;
	for (auto _i = XenEventAction::al28IDs.begin (); _i != XenEventAction::al28IDs.end (); ++_i) {
		if(*_i==parentID){_isON=true;break;}
	}
	return _isON;
}
void XenSteppingAction::_addAl28TrackID(int id)
{
//	//Printing --
//	Logger::log("____________GAMMA IDs_________________");
//	for (auto _i = XenEventAction::gammaIDs.begin (); _i != XenEventAction::gammaIDs.end (); ++_i)
//	{
//		Logger::log(std::to_string(*_i));
//	}
//	///-- Printing

	bool _isON=false;
	for (auto _i = XenEventAction::al28IDs.begin (); _i != XenEventAction::al28IDs.end (); ++_i) {
	    if(*_i==id){_isON=true;break;}
	}
	if(!_isON){
		XenEventAction::al28IDs.push_back(id);
		Logger::log("StepAddAl28TrackID",std::to_string(id));
	}
}
void XenSteppingAction::_printStepSpecs(const G4Step* step)
{
	G4int _cellNo = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
	Logger::log("Step","__________________Start STEP_____________________________");
	std::stringstream sstm;
	G4Track* track     = step->GetTrack();
	const G4ParticleDefinition* part = track->GetDefinition();
	sstm << "[T]Track #"
	<< track->GetTrackID() << " of " << part->GetParticleName()
	<< " E(BestUnit)= " << G4BestUnit(track->GetKineticEnergy(),"Energy")
	<< " produced by Track ID= " << track->GetParentID();
	Logger::log("Step",sstm.str());sstm.clear();

	// Kinetic energy
	G4double ken = track->GetKineticEnergy()/MeV;

	// energy deposit
	G4double edep = step->GetTotalEnergyDeposit();
	//edep/MeV

	sstm<<"ProcessName: "<<step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()<<", Delta:"<<step->GetDeltaMomentum()<<", with edep:"<<edep/MeV<<", on cellID:"<<_cellNo<<", ("<<track->GetPosition().x()<<","<<track->GetPosition().y()<<","<<track->GetPosition().z()<<")"<<std::endl;
	const std::vector<const G4Track*>* secondary = step->GetSecondaryInCurrentStep();
	size_t nbtrk = (*secondary).size();
	if (nbtrk) {
		sstm <<"\n    :----- List of secondaries ----------------" << G4endl;
		sstm.precision(4);
		for (size_t lp=0; lp<(*secondary).size(); lp++) {
		 sstm << "   "
				<< std::setw(13)
				<< (*secondary)[lp]->GetDefinition()->GetParticleName()
				<< ":  energy ="
				<< std::setw(6)
				<< G4BestUnit((*secondary)[lp]->GetKineticEnergy(),"Energy")
				<< "  time ="
				<< std::setw(6)
				<< G4BestUnit((*secondary)[lp]->GetGlobalTime(),"Time");
		 sstm << G4endl;

		}

		sstm << "    :------------------------------------------\n" << G4endl;
	}
		 Logger::log("Step",sstm.str());

	////////////////
}
void XenSteppingAction::_cleanUserAction(const G4Step* step)
{
	/*Load local variables*/
	G4Track* track     = step->GetTrack();
	G4int _cellNo = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
	G4double _delta=std::abs((step->GetDeltaEnergy()/MeV));
	G4double _edep = (step->GetTotalEnergyDeposit())/MeV;
	G4double _totalEnergy=std::abs((step->GetTotalEnergyDeposit()/MeV));
	/*End*/

	/*Step 1. Print specs*/
	_printStepSpecs(step);

	/*Step 2. Check type of particle*///TODO: Should be normalized to Tritons?(i.e. 3m_t=m_p)...
	if(!track -> GetDefinition()->GetParticleName().compare("Al28")){
			_addAl28TrackID(track->GetTrackID());
			//Logger::log("StepAddingName",track -> GetDefinition()->GetParticleName());
	}
	if(track -> GetDefinition()==G4Proton::Definition())
			CellManager::addEnergy(_delta,_cellNo,track->GetVertexMomentumDirection () ,false);
	if(track -> GetDefinition()==G4Triton::Definition())
			CellManager::addEnergy(_delta,_cellNo,track->GetVertexMomentumDirection () ,false);

	if(track -> GetDefinition()==G4Electron::Definition())
			if(_isParentAl28(track->GetParentID()))
			{
				Logger::log("Adding Beta");
				CellManager::addEnergy(_edep,_cellNo,track->GetVertexMomentumDirection () ,true);
			}

	/*Step3. Draw tracks*/
	DrawTracks(true, track,fpSteppingManager);
}
void XenSteppingAction::_fillHistos(const G4Step* step)
{
	if(!step->GetTrack() -> GetDefinition()->GetParticleName().compare("Al28")){
		//for now I only need to fill H1 Id:5 Beta energies
		const std::vector<const G4Track*>* secondary = step->GetSecondaryInCurrentStep();
		size_t nbtrk = (*secondary).size();
		G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
		analysisManager->FillH1(7,step->GetTrack()->GetKineticEnergy()/MeV);
		if (nbtrk) {
			for (size_t lp=0; lp<(*secondary).size(); lp++) {
				if((*secondary)[lp]->GetDefinition()==G4Electron::Definition())
				{
					//child is electron
					G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
					analysisManager->FillH1(5,((*secondary)[lp]->GetKineticEnergy())/MeV);
				}
			}
		}
	}
	if(step->GetTrack()->GetDefinition()==G4Electron::Definition()){
		G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
		analysisManager->FillH1(6,step->GetTrack()->GetKineticEnergy()/MeV);
	}
}
void XenSteppingAction::UserSteppingAction(const G4Step* step)
{
	_fillHistos(step);
	_cleanUserAction(step);
	return;
	//TODO: Clean the code down on.
	Logger::log("Step","__________________Start STEP_____________________________");
	std::stringstream sstm;
    G4Track* track     = step->GetTrack();
    const G4ParticleDefinition* part = track->GetDefinition();
    sstm << "[T]Track #"
    << track->GetTrackID() << " of " << part->GetParticleName()
    << " E(MeV)= " << track->GetKineticEnergy()/MeV
    << " produced by Track ID= " << track->GetParentID();
    Logger::log("Step",sstm.str());sstm.clear();

    // Kinetic energy
    G4double ken = track->GetKineticEnergy()/MeV;

    // energy deposit
    G4double edep = step->GetTotalEnergyDeposit();


    sstm<<"ProcessName: "<<step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()<<", Delta:"<<step->GetDeltaMomentum()<<", with edep:"<<edep/MeV<<", ("<<track->GetPosition().x()<<","<<track->GetPosition().y()<<","<<track->GetPosition().z()<<")"<<std::endl;
    const std::vector<const G4Track*>* secondary = step->GetSecondaryInCurrentStep();
	size_t nbtrk = (*secondary).size();
	if (nbtrk) {
		sstm <<"\n    :----- List of secondaries ----------------" << G4endl;
		sstm.precision(4);
		for (size_t lp=0; lp<(*secondary).size(); lp++) {
		 sstm << "   "
				<< std::setw(13)
				<< (*secondary)[lp]->GetDefinition()->GetParticleName()
				<< ":  energy ="
				<< std::setw(6)
				<< G4BestUnit((*secondary)[lp]->GetKineticEnergy(),"Energy")
				<< "  time ="
				<< std::setw(6)
				<< G4BestUnit((*secondary)[lp]->GetGlobalTime(),"Time");
		 sstm << G4endl;

		}

		sstm << "    :------------------------------------------\n" << G4endl;
	}
         Logger::log("Step",sstm.str());

    ////////////////

      G4StepPoint* preStepPoint = step->GetPreStepPoint();
      G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
      G4int _copyNo = theTouchable->GetCopyNumber();
      G4cout<<"[T][D]CopyNo:"<<_copyNo<<G4endl;
      G4cout<<"[T][A]Accum:"<<XenEventAction::accumDelta<<", lastCellID: "<<XenEventAction::lastCellID<<G4endl;

    /////////////////////
    // step length
    G4double stepLength = 0.;

    stepLength=track->GetPosition().z();

    G4double _path=0;
    G4double realEDep=step->GetTotalEnergyDeposit() - step->GetNonIonizingEnergyDeposit();

    G4double _weight=.01;
    G4int _nIons=0;

    G4double momentumDirZ=track->GetMomentumDirection().z();
    if(step -> GetTrack() -> GetDefinition() == G4Electron::Definition())
	{
		if(_isParentAl28(track->GetParentID()))Logger::log("Stepping","Product of Gamma");
	}
    if(step -> GetTrack() -> GetDefinition() == G4Gamma::Definition())
    {
    	_addAl28TrackID(track->GetTrackID());
    }
    if(step -> GetTrack() -> GetDefinition() == G4Neutron::Definition()){
        if(ken==0)//The end of the neutron track
        {
            XenEventAction::originZ=track->GetPosition().z();
            XenEventAction::originX=track->GetPosition().x();
            XenEventAction::originY=track->GetPosition().y();
            XenEventAction::accumDelta=0;
            XenEventAction::initialMomentum=track->GetMomentumDirection();
            XenEventAction::lastCellID=_copyNo;

            G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
            _path=std::abs(stepLength);
            XenSteppingAction::isNewParticle=true;




        }
        else{G4cout<<"neutron Ken: "<<ken<<", edep: "<<edep<<", stpLength: "<<G4endl; }
    }
    else
    {


        if(step -> GetTrack() -> GetDefinition() == G4Proton::Definition())
        {

            G4cout<<"[PT]: Proton"<<G4endl;
            if(XenSteppingAction::isFirstProton)
            {
                G4cout<<"[M]Overriding momentum lastProton"<<G4endl;
                //OverrideMomentumDirection(step);
                XenSteppingAction::isFirstProton=false;
                //XenSteppingAction::isNewParticle=false;
                track=step->GetTrack();

            }



            //This is a proton
            if(XenEventAction::originZ!=0){


                G4double currentX=track->GetPosition().x();
                G4double currentY=track->GetPosition().y();
                G4double currentZ=track->GetPosition().z();

                G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

                _path=std::abs(stepLength);
                G4double _dist=std::sqrt(std::pow(XenEventAction::originX-currentX,2)+std::pow(XenEventAction::originY-currentY,2)+std::pow(XenEventAction::originZ-currentZ,2));


                G4double _delta=std::abs((step->GetDeltaEnergy()/MeV));
                G4double _totalEnergy=std::abs((step->GetTotalEnergyDeposit()/MeV));

                _nIons=_delta/_weight;
               // G4int _discreet=rounDrawTracks(true, track,fpSteppingManager);d(_dist);
                XenRunAction::_energy+=_delta;
                if(XenEventAction::lastCellID!=_copyNo)
                {
                    G4cout<<"[A]EnergyDep:"<<XenEventAction::accumDelta<<",delta:"<<_delta<<G4endl;
                    //man->FillH2(1, i,j,CellManager::getGFactor(j,i)*10);
                    CellManager::addEnergy(XenEventAction::accumDelta,XenEventAction::lastCellID,track->GetVertexMomentumDirection () ,false);
                    XenEventAction::accumDelta=0;
                }
                else
                {
                    G4cout<<"[A]LastCell=copyNO"<<_copyNo<<G4endl;
                    XenEventAction::accumDelta+=_delta;

                }
                XenEventAction::lastCellID=_copyNo;


//                analysisManager->FillH1(1, _discreet,_nIons);
//                G4double _theta=track->GetMomentumDirection().theta();
//                G4cout<<"theta:"<<_theta<<", cos():"<<cos(_theta)<<G4endl;
//                analysisManager->FillH1(2,roundf(cos(_theta)*100)/100);
//                analysisManager->FillH1(3,_theta);

            }

        }


        else if(step -> GetTrack() -> GetDefinition() == G4Triton::Definition())
        {
            XenSteppingAction::isFirstProton=true;
            if(XenSteppingAction::isNewParticle)
            {
                //OverrideMomentumDirection(step);
                XenSteppingAction::isNewParticle=false;
//                track=step->GetTrack();
//                double _currentZMomentum=track->GetMomentumDirection().z();
//                G4cout<<"[M]triton: zchanGEd:"<<_currentZMomentum<<G4endl;
//                G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//                analysisManager->FillH1(5,_currentZMomentum);

            }

            G4cout<<"[PT]Is triton!"<<G4endl;
            G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
            G4double currentX=track->GetPosition().x();
            G4double currentY=track->GetPosition().y();
            G4double currentZ=track->GetPosition().z();
            _path=std::abs(stepLength);

            G4double _dist=std::sqrt(std::pow(XenEventAction::originX-currentX,2)+std::pow(XenEventAction::originY-currentY,2)+std::pow(XenEventAction::originZ-currentZ,2));
            G4double _delta=std::abs((step->GetDeltaEnergy()/MeV));
            G4double _totalEnergy=std::abs((step->GetTotalEnergyDeposit()/MeV));
            XenRunAction::_energy+=_delta;
            G4int _discreet=std::abs(round(_dist))*-1;
            _nIons=_delta/_weight;
            analysisManager->FillH1(1, _discreet,_nIons);

            if(XenEventAction::lastCellID!=_copyNo)
            {
                G4cout<<"[T][A]EnergyDep:"<<XenEventAction::accumDelta<<",delta:"<<_delta<<G4endl;
                //man->FillH2(1, i,j,CellManager::getGFactor(j,i)*10);
                CellManager::addEnergy(XenEventAction::accumDelta,XenEventAction::lastCellID,track->GetVertexMomentumDirection () ,false);
                XenEventAction::accumDelta=0;
            }
            else
            {
                G4cout<<"[T][A]LastCell=copyNO"<<_copyNo<<G4endl;
                XenEventAction::accumDelta+=_delta;

            }
            XenEventAction::lastCellID=_copyNo;
        }
        else if(step -> GetTrack() -> GetDefinition() == G4Electron::Definition())
		{
        	if(_isParentAl28(track->GetParentID()))
        	{
        		Logger::log("Stepping","Product of Gamma");
				G4double _delta=std::abs((step->GetDeltaEnergy()/MeV));
				XenRunAction::_bEnergy+=_delta;
				if(XenEventAction::lastCellID!=_copyNo)
				{
					G4cout<<"[A]EnergyDep:Beta"<<XenEventAction::accumDelta<<",delta:"<<_delta<<G4endl;
					//man->FillH2(1, i,j,CellManager::getGFactor(j,i)*10);
					CellManager::addEnergy(XenEventAction::accumDelta,XenEventAction::lastCellID,track->GetVertexMomentumDirection (),true );
					XenEventAction::accumDelta=0;
				}
				else
				{
					G4cout<<"[A]LastCell=copyNO_BETA"<<_copyNo<<G4endl;
					XenEventAction::accumDelta+=_delta;

				}
        	}
			XenEventAction::lastCellID=_copyNo;
        	//std::cout<<"I'm an electron, delta:"<<_delta<<std::endl;
		}
    }

    DrawTracks(true, track,fpSteppingManager);
    return;
}

void XenSteppingAction::DrawTracks(G4bool drawFlag, G4Track *theTrack,G4SteppingManager* pSM)
{



      G4ParticleDefinition *particleType = theTrack->GetDefinition();

      G4Colour red      ( 255/255.,   0/255.,   0/255.);
      G4Colour blue     (   0/255.,   0/255., 255/255.);
      G4Colour green    (   0/255., 255/255.,   0/255.);
      G4Colour yellow   ( 255/255., 255/255.,   0/255.);

      G4Colour white    ( 255/255., 255/255., 255/255.);

      G4Colour orange   ( 255/255., 127/255.,   0/255.);
      G4Colour magenta  ( 237/255., 173/255., 255/255.);
      G4Colour magenta1 ( 104/255.,  49/255.,  94/255.);

      G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

      if (pVVisManager) {
          // Declare begininng of visualization

        G4Polyline polyline;
        G4Colour colour;
        G4cout<<"ParType:"<<particleType->GetParticleType()<<G4endl;
        if( particleType == G4Gamma::GammaDefinition()){
          colour = green;
          G4cout<<"GREEN_GAMMA!!"<<G4endl;
        }
        if( particleType == G4Triton::TritonDefinition()){
          colour = yellow;
          G4cout<<"TRITON_YELLOW!!"<<G4endl;
        }
        if( particleType == G4Neutron ::NeutronDefinition())            {
            G4cout<<"NEUTRON_WHITE!!"<<G4endl;
          colour = white; }
        if( particleType == G4Proton  ::ProtonDefinition())
          { G4cout<<"PROTON_BLUE!!"<<G4endl; colour = blue; }
        if( particleType == G4Electron  ::ElectronDefinition())
          { colour = red; G4cout<<"ELECTRON_RED!!"<<G4endl;}

        G4VisAttributes attribs(colour);
        polyline.SetVisAttributes(attribs);
        polyline.push_back(pSM->GetStep()->GetPreStepPoint()->GetPosition());
        polyline.push_back(pSM->GetStep()->GetPostStepPoint()->GetPosition());
        pVVisManager -> Draw(polyline);


      }

}

void XenSteppingAction::Reset()
{
  fEnergy = 0.;
}
void XenSteppingAction::OverrideMomentumDirection(G4Step* &step)
{

//    double _preX=fpSteppingManager->GetStep()->GetPreStepPoint()->GetPosition().x();
//    double _preY=fpSteppingManager->GetStep()->GetPreStepPoint()->GetPosition().y();
//    double _preZ=fpSteppingManager->GetStep()->GetPreStepPoint()->GetPosition().z();
//
//    double _postX=fpSteppingManager->GetStep()->GetPostStepPoint()->GetPosition().x();
//    double _postY=fpSteppingManager->GetStep()->GetPostStepPoint()->GetPosition().y();
//    double _postZ=fpSteppingManager->GetStep()->GetPostStepPoint()->GetPosition().z();
//
//    double _originalLength=sqrt(pow(_postX-_preX,2)+pow(_postY-_preY,2)+pow(_postZ-_preZ,2));
//
//    //Generates the new direction
//
//    double* _momentum;
//
//     if(step -> GetTrack() -> GetDefinition() == G4Triton::Definition())
//     {
//         _momentum=MomentumManager::getRandomMomentumDirection();
//         _momentum[0]*=-1;_momentum[1]*=-1;_momentum[2]*=-1;
//          std::cout<<"[M]Overriding momentum new"<<std::endl;
//     }
//     else{
//         _momentum=MomentumManager::getLastMomentumDirection();
//        std::cout<<"[M]Overriding momentum last"<<std::endl;
//     }
//
//    for(int i=0; i<3; i++)
//    {
//       // _momentum[i]=0;
//        std::cout<<"[M]Overriding momentum tooo, "<<i<<":"<<_momentum[i]<<std::endl;
//
//    }
//    //_momentum[2]=-1;
//
    G4Track* track= step->GetTrack();
//    track->SetMomentumDirection(G4ThreeVector(_momentum[0],_momentum[1],_momentum[2]));
//
//    track->SetVertexMomentumDirection(G4ThreeVector(_momentum[0],_momentum[1],_momentum[2]));
//
//    track->SetPosition(G4ThreeVector(_preX+_originalLength*_momentum[0],_preY+_originalLength*_momentum[1],_preZ+_originalLength*_momentum[2]));
//
//    step->SetTrack(track);
//
//    G4StepPoint* _point=step->GetPreStepPoint();
//    _point->SetMomentumDirection(G4ThreeVector(_momentum[0],_momentum[1],_momentum[2]));
//    step->SetPreStepPoint(_point);
//
//
//    G4StepPoint* _point2=step->GetPostStepPoint();
//    _point2->SetMomentumDirection(G4ThreeVector(_momentum[0],_momentum[1],_momentum[2]));
//    _point2->SetPosition(G4ThreeVector(_preX+_originalLength*_momentum[0],_preY+_originalLength*_momentum[1],_preZ+_originalLength*_momentum[2]));
//    step->SetPostStepPoint(_point2);
//
//
//    double _posX2=fpSteppingManager->GetStep()->GetPostStepPoint()->GetPosition().x();
//    double _posY2=fpSteppingManager->GetStep()->GetPostStepPoint()->GetPosition().y();
//    double _posZ2=fpSteppingManager->GetStep()->GetPostStepPoint()->GetPosition().z();
////    std::cout<<"distance/prePoint x:"<<_preX<<",y:"<<_preY<<",z:"<<_preZ<<"postPoint/ x:"<<_postX<<",y:"<<_postY<<",z:"<<_postZ<<_originalLength<<std::endl;
////    std::cout<<"changing to, x:"<<_preX+_originalLength*_momentum[0]<<", y:"<<_preY+_originalLength*_momentum[1]<<",z:"<<_preZ+_originalLength*_momentum[2]<<std::endl;
////    std::cout<<"pre/ x:"<<_postX<<",y:"<<_postY<<",z:"<<_postZ<<",x2:"<<_posX2<<",y2:"<<_posY2<<",z2:"<<_posZ2<<", originalLength:"<<_originalLength<<std::endl;
//
////    //Histograms
//
    double _currentZMomentum=track->GetVertexMomentumDirection().z();

    G4cout<<"[M]: zchanGEd:"<<_currentZMomentum<<G4endl;
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(5,track->GetMomentumDirection().x());
    analysisManager->FillH1(6,track->GetMomentumDirection().y());
    analysisManager->FillH1(7,_currentZMomentum);

    _currentZMomentum=MomentumManager::getRandomNumber();
    analysisManager->FillH1(8,_currentZMomentum);
}



