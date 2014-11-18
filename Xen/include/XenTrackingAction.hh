
#ifndef XenTrackingAction_h
#define XenTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "G4TrackVector.hh"
#include "globals.hh"

#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4UserTrackingAction.hh"
#include "G4VUserTrackInformation.hh"
#include "G4TrackingManager.hh"
#include "G4TrackVector.hh"


class XenTrackingAction : public G4UserTrackingAction
{
  public:

    XenTrackingAction();
    virtual ~XenTrackingAction();


    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);



};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
