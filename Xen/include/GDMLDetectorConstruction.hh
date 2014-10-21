#ifndef _GDMLDETECTORCONSTRUCTION_H_
#define _GDMLDETECTORCONSTRUCTION_H_

#include "G4VUserDetectorConstruction.hh"

class GDMLDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    GDMLDetectorConstruction(G4VPhysicalVolume *setWorld = 0)
    {
      World = setWorld;
    }

    G4VPhysicalVolume *Construct()
    {
      return World;
    }

  private:

    G4VPhysicalVolume *World;
};

#endif
