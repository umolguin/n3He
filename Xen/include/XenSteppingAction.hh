#ifndef XenSteppingAction_h
#define XenSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4UserTrackingAction.hh"
#include "G4VUserTrackInformation.hh"
#include "G4TrackingManager.hh"
#include "G4TrackVector.hh"

#include "G4Colour.hh"

#include "G4Polyline.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"


#include "G4Gamma.hh"
#include "G4Triton.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Electron.hh"

#include "XenEventAction.hh"

//class XenEventAction;
class G4LogicalVolume;
class XenDetectorConstruction;

/// Stepping action class
///
/// It holds data member fEnergy for accumulating the energy deposit
/// in a selected volume step by step.
/// The selected volume is set from  the detector construction via the
/// SetVolume() function. The accumulated energy deposit is reset for each
/// new event via the Reset() function from the event action.

class XenSteppingAction : public G4UserSteppingAction
{
  public:
  XenSteppingAction(XenDetectorConstruction* detectorConstruction,
                    XenEventAction* eventAction);
    XenSteppingAction();
    virtual ~XenSteppingAction();

    // static access method
    static XenSteppingAction* Instance();

    // method from the base class
    virtual void UserSteppingAction(G4Step*);

    // reset accumulated energy
    void Reset();

    // set methods
    void SetVolume(G4LogicalVolume* volume) { fVolume = volume; }

    // get methods
    G4LogicalVolume* GetVolume() const { return fVolume; }
    G4double GetEnergy() const { return fEnergy; }
    static bool isNewParticle;
    static bool isFirstProton;

  private:
    static XenSteppingAction* fgInstance;

    G4LogicalVolume* fVolume;
    G4double  fEnergy;
    void DrawTracks(G4bool, G4Track*, G4SteppingManager*);
    void OverrideMomentumDirection(G4Step* &step);

    G4TrackVector* fSecondary;
    XenEventAction* fEventAction;
    XenDetectorConstruction* fDetConstruction;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
