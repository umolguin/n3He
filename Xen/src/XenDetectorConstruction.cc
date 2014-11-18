#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UserLimits.hh"

#include "XenDetectorConstruction.hh"
 #include "GDMLDetectorConstruction.hh"
#include "XenSteppingAction.hh"

#include "G4GDMLParser.hh"
#include "G4UnitsTable.hh"

#include "CellManager.hh"

using namespace std;

///@author Carlos Olguin
///@version 0.8.1
///
// Notes:
// - When the world is less than 300 cm neutrons interact incorrectly and differently between
// different rows.
// - Array of "real" cells was tried against cloning the logical volumes. Bug previously described is
// still present.
//
//



XenDetectorConstruction::XenDetectorConstruction(): G4VUserDetectorConstruction()
{
    _visible=true;
}
XenDetectorConstruction::XenDetectorConstruction(bool visible): G4VUserDetectorConstruction()
{
    _visible=visible;
}
XenDetectorConstruction::~XenDetectorConstruction(){}

G4VPhysicalVolume* XenDetectorConstruction::Construct()
{
	/// IDs
	const int _frontWindowID=3;
	const int _backWindowID=9;
	///
    //Materials
    G4double density1 = 2.700*g/cm3;
    G4double a = 26.98*g/mole;
    G4Material* Al = new G4Material("Aluminum", 13., a, density1);
//
//
//    G4Isotope* He3Iso = new G4Isotope("3HeIso", 2,3,3.016023*g/mole);
//    G4Element* elHe3 = new G4Element("Element_Helium_3","3He",1);
//    elHe3->AddIsotope(He3Iso,100.0*perCent);
//    G4Material* matHelium3Gas = new G4Material("He3Gas",0.1340*mg/cm3,1,kStateGas,273.15*kelvin,1.*atmosphere);
//    matHelium3Gas->AddElement(elHe3,100.0*perCent);

    G4int protons=2, neutrons=1, nucleons=protons+neutrons;
    G4double atomicMass = 3.016*g/mole;
    G4Isotope* he3 = new G4Isotope("He3", protons, nucleons, atomicMass);
    G4int isotopes, elements;
    G4Element* He3 = new G4Element("Helium3", "He3", isotopes=1);
    He3->AddIsotope(he3, 100*perCent);

    G4double pressure =.5*atmosphere;
    G4double temperature = 273.15*kelvin;
    G4double molar_constant = Avogadro*k_Boltzmann;  //from clhep
    G4double density = (atomicMass*pressure)/(temperature*molar_constant);// Ideal gas law
    G4Material* matHelium3Gas = new G4Material("Helium3", density, elements=1, kStateGas, temperature, pressure);
    matHelium3Gas->AddElement(He3, 100*perCent);


    //_____________GDML
//    G4GDMLParser parser;
//    parser.Read("Chamber3DMAX.gdml");
////    parser.Read("xwing.gdml");
//
////    Add properties
//    G4VPhysicalVolume* W=parser.GetWorldVolume();
//    G4LogicalVolume* lv = W->GetLogicalVolume();
//    //_____________GDML

//    G4double maxStep = 2*mm;
//    Set the step length to 2 mm
//    lv->SetUserLimits(new G4UserLimits(maxStep));

    //Create a cylinder filled with the gas. Target.
//    G4Tubs* _cylSol= new G4Tubs("Targett",0,32.4*cm,43*cm,0,twopi);
//    G4LogicalVolume* _cylLog= new G4LogicalVolume(_cylSol,matHelium3Gas,"Targett");

    G4Tubs* _cylTar1= new G4Tubs("Target1",0,12.7*cm,(1.25/2)*cm,0,twopi); //16.51 for dz since is half of the length of the total cells along z plus 1.27 cm space, between window and wireframe
    G4LogicalVolume* _cylTar1Log= new G4LogicalVolume(_cylTar1,matHelium3Gas,"Target1");

//    G4Tubs* _cylHousing= new G4Tubs("House",0,13.335*cm,17.018*cm,0,twopi);
//    G4LogicalVolume* _cylHLog= new G4LogicalVolume(_cylHousing,matHelium3Gas,"House");
//
//    G4SubtractionSolid* _cylHos=new G4SubtractionSolid("Hos",_cylHousing,_cylSol);
//    G4LogicalVolume* _cylHosLog= new G4LogicalVolume(_cylHos,matHelium3Gas,"Hos");

    //G4Box* worldBox = new G4Box("World", world_hx, world_hy, world_hz);
    // allocate memory
//    G4Box** _cells = static_cast<G4Box**>( ::operator new ( sizeof(G4Box) * 144 ) );
//    G4LogicalVolume** _lCells=static_cast<G4LogicalVolume**>( ::operator new ( sizeof(G4LogicalVolume) * 144 ) );

//    G4Box *_cellsDyna[144];//This implementation appears to work better with LastHopeU
//    G4LogicalVolume *_lCellsDyna[144];


    //G4Box* _cells= new G4Box("Cell",8.001*cm,.95*cm,.95*cm);
    double _cX=16.02*cm; double _cY=1.905*cm; double _cZ=1.905*cm;

    G4Box* _cells= new G4Box("Cell",_cX/2,_cY/2,_cZ/2);
    G4LogicalVolume* _lCells= new G4LogicalVolume(_cells,matHelium3Gas,"Cell");
//
//    G4Box* _cells2= new G4Box("Cell2",_cX/2,(_cY*9)/2,1.26*cm);
//    G4LogicalVolume* _lCells2= new G4LogicalVolume(_cells2,matHelium3Gas,"Cell2");

//    G4SubtractionSolid* _cylHos=new G4SubtractionSolid("Hos",_cylSol,_cells2);
//	G4LogicalVolume* _cylHosLog= new G4LogicalVolume(_cylHos,matHelium3Gas,"Hos");




    // invoke constuctors


    //Set material to the whole world
   // lv->SetMaterial(matHelium3Gas);

    //If command was sent to make the volume invisible then
//    if(!_visible)
//    {
//
//        _cylLog->SetVisAttributes (G4VisAttributes::Invisible);
////        new G4PVPlacement(0, G4ThreeVector(0,0,-35.2*cm),_cylLog,"Targett",lv,false,0,true);
//        //new G4PVPlacement(0, G4ThreeVector(0,0,-16.9164*cm),_cylLog,"Targett",lv,false,0,true);
//        new G4PVPlacement(0, G4ThreeVector(0,0,0),_cylHosLog,"Hos",lv,false,0,true);
//
//        for(int i=1; i<387;i++)
//        {
//            stringstream ss;//create a stringstream
//            ss << i;//add number to the stream
//
//            string _name = string("chamber3DMAX_") + ss.str()+string("_vol");
//            G4LogicalVolume* lv = parser.GetVolume( _name);
//            lv->SetVisAttributes (G4VisAttributes::Invisible);
//        }
//
//        G4LogicalVolume* lv1 = parser.GetVolume( "chamber3DMAX_vol");
//        lv1->SetVisAttributes (G4VisAttributes::Invisible);
//
//
//    }
//    else{
       // new G4PVPlacement(0, G4ThreeVector(0,0,-16.9164*cm),_cylTar1Log,"Targett",lv,false,0,true);
//        new G4PVPlacement(0, G4ThreeVector(0,0,0),_cylHosLog,"Hos",lv,false,0,true);
//    }
    //new G4PVPlacement(0, G4ThreeVector(0,0,0),_cylHosLog,"Hos",lv,false,0,true);
//    new G4PVPlacement(0,G4ThreeVector(0,0,(15.24+1.26)*cm),_cylTar1Log,"Tar1",lv,false,_frontWindowID,true);
//    new G4PVPlacement(0,G4ThreeVector(0,0,(-15.24-1.26)*cm),_cylTar1Log,"Tar1",lv,false,_backWindowID,true);

//    //////////////////////////??????????????????????????
    G4Material* vacuum =
         new G4Material("Vacuum",      //Name as String
    					1,		       //Atomic Number,  in this case we use 1 for hydrogen
                        1.008*g/mole,  //Mass per Mole "Atomic Weight"  1.008*g/mole for Hydrogen
    					1.e-25*g/cm3,  //Density of Vacuum  *Can't be Zero, Must be small instead
    					kStateGas,     //kStateGas for Gas
                        2.73*kelvin,   //Temperature for ga
    					1.e-25*g/cm3); //Pressure for Vacuum
    G4VSolid* worldS
      = new G4Box("World",           // its name
                   400*cm,400*cm, 400*cm); // its size

    G4LogicalVolume* lv
      = new G4LogicalVolume(
                   worldS,           // its solid
                   vacuum,  // its material
                   "World");         // its name

    G4VPhysicalVolume* W
      = new G4PVPlacement(
                   0,                // no rotation
                   G4ThreeVector(),  // at (0,0,0)
                   lv,          // its logical volume
                   "World",          // its name
                   0,                // its mother  volume
                   false,            // no boolean operation
                   0,                // copy number
                   true);  // checking overlaps

////    //////////////???????????????????????????????????????????
//	new G4PVPlacement(0,G4ThreeVector(0,0,(15.24+(1.25/2))*cm),_cylTar1Log,"Tar1",lv,false,_frontWindowID,true);
//	new G4PVPlacement(0,G4ThreeVector(0,0,(-15.24-(1.25/2))*cm),_cylTar1Log,"Tar1",lv,false,_backWindowID,true);
    for(int i=0; i<9;i++)
        for(int j=0;j<16;j++)
        {
//        	int _index=i*16+j;
//        	std::stringstream sstm;
//        	sstm << "Box" << _index;
            new G4PVPlacement(0, G4ThreeVector(0,((8.57-1.905/2)*cm)-(i*_cY),((15.24-1.905/2)*cm)-(j*_cZ)),_lCells,"CELLS",lv,false,10+(j+(i*16)),true);
//            sstm.clear();
        }

//    new G4PVPlacement(0, G4ThreeVector(0,0,(0)*cm),_lCells2,"Cells2",lv,false,0,true);
//    new G4PVPlacement(0, G4ThreeVector(0,0,(15.24+1.27)*cm),_lCells2,"Cells3",lv,false,1,true);

    CellManager::setOrigin(0,0,0);//Taken from the G4PVPlacement
    //Attach volumes

    W->SetLogicalVolume(lv);
    return W;
}
