#ifndef XenDetectorConstruction_h
#define XenDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;

/// Detector construction class to define materials and geometry.

class XenDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    XenDetectorConstruction();
    XenDetectorConstruction(bool visible);
    virtual ~XenDetectorConstruction();



  public:
    virtual G4VPhysicalVolume* Construct();
    bool _visible;
  private:
    G4VPhysicalVolume* _GDML();
    G4VPhysicalVolume* _Air();
    G4VPhysicalVolume* _CellsOnly();
    G4VPhysicalVolume* _CellsOnlyGarfield();
    G4VPhysicalVolume* _CellsOnlyAl();
    G4VPhysicalVolume* _SmallBox();
    G4VPhysicalVolume* _FilledWorld();


};
#endif
