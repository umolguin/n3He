#ifndef B1RunAction_h
#define B1RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#include "XenTrackerHit.hh"
#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"

class G4Run;


/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.
extern int fCollID;
class XenRunAction : public G4UserRunAction
{

  public:

    XenRunAction();
    virtual ~XenRunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    static double _energy;
    static double _bEnergy;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
