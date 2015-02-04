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
        static void addEnergy(G4double energy, G4int arrayPos, G4ThreeVector momentum, bool isBeta, G4String source);
        static G4double getGFactor(G4int row, G4int column);
        static G4double getDFactor(G4int row, G4int column, G4String source);
        static G4double getEnergyFromBeta_al(G4int row, G4int column);
        static G4double getEnergyFromBeta_cu(G4int row, G4int column);
        static G4double getEnergyFromP_T(G4int row, G4int column);
        static void logEvent(G4double index, G4double energy, G4int &x, G4int &y, G4int &z);
        static void printLoggedEvents();
        static double* sumCosEnergy;
        static double* sumEnergy;
        static double sumNonCellEnergy;
        static double sumEnergyFrontWindow;
        static double sumEnergyBackWindow;
        static double sumCosEnergyFrontWindow;
		static double sumCosEnergyBackWindow;
		//TODO:EWWWWWWWW this is awful, fix it!
        static double* sumCosBEnergy_al;
        static double* sumBEnergy_al;
        static double sumNonCellBEnergy_al;
        static double sumBEnergyFrontWindow_al;
        static double sumBEnergyBackWindow_al;
        static double sumCosBEnergyFrontWindow_al;
		static double sumCosBEnergyBackWindow_al;

        static double* sumCosBEnergy_cu;
        static double* sumBEnergy_cu;
        static double sumNonCellBEnergy_cu;
        static double sumBEnergyFrontWindow_cu;
        static double sumBEnergyBackWindow_cu;
        static double sumCosBEnergyFrontWindow_cu;
		static double sumCosBEnergyBackWindow_cu;

        static void print();
        static void g_genIndices(G4int &x, G4int &y, G4int &z,G4double);

    private:
        static G4double nFrames;
        static G4double nWires;
        static G4double spcBtnFrms;
        static G4double spcBtnWires;
        static G4double wireLength;
        static G4double spcFrontWin_1stFrm;
        static G4double spcBckWin_lastFrm;
        static G4double initX;
        static G4double initY;
        static G4double initZ;
        static G4double GarfieldCellSize;
        static G4int offsetCellId;
        static double g_transGFactorX;
        static double g_transGFactorY;
        static double g_transGFactorZ;


};

#endif
