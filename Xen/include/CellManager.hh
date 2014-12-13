#ifndef CellManager_h
#define CellManager_h 1

#include "globals.hh"
#include "G4Track.hh"

class CellManager
{
    public:

        static void init();
        static void setOrigin(G4double x, G4double y, G4double z);
        static void addEnergy(G4double energy, G4ThreeVector pos, G4ThreeVector momentum);
        static void addEnergy(G4double energy, G4int arrayPos, G4ThreeVector momentum, bool isBeta);
        static G4double getGFactor(G4int row, G4int column);
        static G4double getDFactor(G4int row, G4int column);
        static G4double getEnergyFromBeta(G4int row, G4int column);
        static G4double getEnergyFromP_T(G4int row, G4int column);
        static double* sumCosEnergy;
        static double* sumEnergy;
        static double sumNonCellEnergy;
        static double sumEnergyFrontWindow;
        static double sumEnergyBackWindow;
        static double sumCosEnergyFrontWindow;
		static double sumCosEnergyBackWindow;
		//TODO:EWWWWWWWW this is awful, fix it!
        static double* sumCosBEnergy;
        static double* sumBEnergy;
        static double sumNonCellBEnergy;
        static double sumBEnergyFrontWindow;
        static double sumBEnergyBackWindow;
        static double sumCosBEnergyFrontWindow;
		static double sumCosBEnergyBackWindow;
        static void print();

    private:
        static G4double nFrames;
        static G4double nWires;
        static G4double spcBtnFrms;
        static G4double spcBtnWires;
        static G4double spcFrontWin_1stFrm;
        static G4double spcBckWin_lastFrm;
        static G4double initX;
        static G4double initY;
        static G4double initZ;


};

#endif
