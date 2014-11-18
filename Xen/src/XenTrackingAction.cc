#include "XenTrackingAction.hh"
#include "G4ThreeVector.hh"
#include "G4RandomDirection.hh"
#include "G4Proton.hh"
#include "G4Triton.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

XenTrackingAction::XenTrackingAction()
: G4UserTrackingAction(){}
XenTrackingAction::~XenTrackingAction()
{
}
void XenTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
	 if(aTrack != NULL){
	     if(aTrack->GetParentID() == 0){
	       G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
	       if(secondaries)
	     {
	       G4int nSeco = secondaries->size();
	       if(nSeco>0)
	         {

	           G4ThreeVector dir =  G4RandomDirection();
	           G4ThreeVector Odir(-dir.x(),-dir.y(),-dir.z());

	           for(G4int i=0;i<nSeco;i++)
	         {
	           if((*secondaries)[i]->GetDefinition() ==
	G4Proton::ProtonDefinition())
	             (*secondaries)[i]->SetMomentumDirection(dir.unit());

	           if((*secondaries)[i]->GetDefinition() ==
	G4Triton::TritonDefinition())
	             (*secondaries)[i]->SetMomentumDirection(Odir.unit());

	         }
	         }
	     }
	     }
	   }
}
void XenTrackingAction::PreUserTrackingAction(const G4Track* track)
{
   //track->SetMomentumDirection(G4ThreeVector(1,0,0.0));
}


