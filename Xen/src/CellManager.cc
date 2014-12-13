#include "CellManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


G4double CellManager::nFrames=0;
G4double CellManager::nWires=0;
G4double CellManager::spcBtnWires=0;
G4double CellManager::spcBtnFrms=0;
G4double CellManager::spcFrontWin_1stFrm=0;
G4double CellManager::spcBckWin_lastFrm=0;
G4double CellManager::initX=0;
G4double CellManager::initY=0;
G4double CellManager::initZ=0;

G4double* CellManager::sumCosEnergy;
G4double* CellManager::sumEnergy;
G4double CellManager::sumNonCellEnergy;

G4double CellManager::sumCosEnergyFrontWindow;
G4double CellManager::sumEnergyFrontWindow;
G4double CellManager::sumCosEnergyBackWindow;
G4double CellManager::sumEnergyBackWindow;

G4double* CellManager::sumCosBEnergy;
G4double* CellManager::sumBEnergy;
G4double CellManager::sumNonCellBEnergy;

G4double CellManager::sumCosBEnergyFrontWindow;
G4double CellManager::sumBEnergyFrontWindow;
G4double CellManager::sumCosBEnergyBackWindow;
G4double CellManager::sumBEnergyBackWindow;

using namespace std;

void CellManager::print()
{
//    double _sum=0;
//    for(int i=0;i<16;i++)
//    {   G4cout<<"[P]:";
//        for(int j=0;j<9;j++){
//            G4int _index=i*nWires+j;
//            _sum=_sum+CellManager::getGFactor(j,i);
//            G4cout<<CellManager::getGFactor(j,i)<<"     ";
//        }
//        G4cout<<G4endl;
//    }
//    G4cout<<"[P]:144:"<<(sumCosEnergy[144]/sumEnergy[144])<<",145:"<<(sumCosEnergy[145]/sumEnergy[145])<<",total:"<<_sum<<G4endl;
//    G4cout<<"[P]:---------------"<<G4endl;
	//START Writing to file
	cout<<"PRINTING!!!!!!!!!"<<endl;
	std::ofstream ofs ("geometricFactor.dat");
	ofs <<"FrontWINDOW:"<<sumEnergyFrontWindow<<endl;
	for(int i=0;i<nWires;i++)
	{
		for(int j=0;j<nFrames;j++)
		{
			G4int _index=i*nFrames+j;
			ofs <<"	"<<sumEnergy[_index];
		}
		ofs<<endl;
	}
    ofs.close();
    //END Writing

    //START Writing to file
    std::ofstream ofs1 ("dilutionFactor.dat");
	ofs1 <<sumEnergyFrontWindow<<"	"<<sumBEnergyFrontWindow<<endl;
	for(int i=0;i<143;i++)
	{

			if((i%15)==0)ofs1<<endl;
			ofs1<<sumEnergy[i]<<"	"<<sumBEnergy[i]<<"	-	";

	}
	ofs1<<endl;

	ofs1.close();
	//END Writing
}

///ZERO INDEXED!!
G4double CellManager::getGFactor(G4int row, G4int column)//9,16
{
    G4int _index=row*nFrames+column;
    G4double _gFactor;
    if(sumEnergy[_index]>0)
        _gFactor=sumCosEnergy[_index]/sumEnergy[_index];
    else
        _gFactor=0;
    return _gFactor;
}

G4double CellManager::getDFactor(G4int row, G4int column)
{
    G4int _index=column*nFrames+row;
    G4double _dFactor;
    if(sumEnergy[_index]>0)
        _dFactor=sumBEnergy[_index]/sumEnergy[_index];
    else
        _dFactor=0;
    return _dFactor;
}
G4double CellManager::getEnergyFromP_T(G4int row, G4int column )
{
	G4int _index=row*nFrames+column;
	return sumEnergy[_index];
}
G4double CellManager::getEnergyFromBeta(G4int row, G4int column )
{
	G4int _index=row*nFrames+column;
	return sumBEnergy[_index];
}

void CellManager::init()
{
    //init for now hardcoded
    nFrames=16;
    nWires=9;
    spcBtnWires=1.9*cm;
    spcBtnFrms=1.9*cm;
    spcFrontWin_1stFrm=1.2954*cm;
    spcBckWin_lastFrm=2.1336*cm;
    initX=0;
    initY=0;
    initZ=0;

    //calculate the number of cells, flat representation of the volume
    int nCells=(nWires-1)*nFrames;
    sumCosEnergy= new G4double[static_cast<int>(nFrames*nWires+3)]; //The array contains nFrames*nWires cells, plus space between
    sumEnergy= new G4double[static_cast<int>(nFrames*nWires+3)];      //front and rear window

    sumCosBEnergy= new G4double[static_cast<int>(nFrames*nWires+3)]; //The array contains nFrames*nWires cells, plus space between
	sumBEnergy= new G4double[static_cast<int>(nFrames*nWires+3)];      //front and rear window

    for(int i=0;i<145;i++)
    {
            sumEnergy[i]=0;
            sumCosEnergy[i]=0;

            sumBEnergy[i]=0;
			sumCosBEnergy[i]=0;
    }
}
void CellManager::setOrigin(G4double x, G4double y, G4double z)
{
    initX=x;
    initY=y;
    initZ=z;
}
void CellManager::addEnergy(G4double energy, G4ThreeVector pos, G4ThreeVector momentum)
{

    G4double _posX=pos.x()-initX;
    G4double _posY=pos.y()-initY;
    G4double _posZ=pos.z()-initZ;

    //To get the position relative to the origin of the frames we add the space between the front
    //window and the first frame plane
    _posZ+=spcFrontWin_1stFrm;

    //assume wires run across x axis, but we calculate the angle between the y axis

    //cells are zero indexed!!!
    //calculate what frame I am
    G4int _columnZ=(floor((nFrames/2)+_posZ/spcBtnFrms));


    //Calculate which wire cell I am in the y axis
    G4int _rowY=(floor((nWires/2)-_posY/spcBtnWires));



    G4int _index=_columnZ*nWires+_rowY;
    if(_columnZ>=16)_index=144;// is between first frame and front window
    if(_columnZ<0)_index=145;// is between last frame and back window
    if(_rowY>=9||_rowY<0)_index=146;
    //if(_index<0)_index=145;

//    if(_index==144||_index==145)
//    {
//        G4cout<<"[P]:("<<pos.x()<<","<<pos.y()<<","<<pos.z()<<"), columnZ:"<<_columnZ<<", rowY:"<<_rowY<<G4endl;
//    }
    //Calculate theta (angle between y and x axis, starting from y)
    G4double _theta=0;
    G4double _phi=momentum.phi();

    if(_phi>=0)
    {
        if(_phi>=(pi/2))
            _theta=_phi-pi/2;
        else
            _theta=(pi/2)-_phi;
    }
    else
    {
        if(_phi<=(-1*pi/2))
        {
            _theta=(pi/2)+pi+_phi;

        }
        else
        {
            _theta=(pi/2)-_phi;

        }
    }

    sumCosEnergy[_index]+=cos(_theta)*energy;
    sumEnergy[_index]+=energy;

    //G4cout<<"[CM]:energy="<<energy<<", theta="<<_theta<<", phi="<<_phi<<", index="<<_index<<", pos=("<<pos.x()<<","<<pos.y()<<","<<pos.z()<<")"<<G4endl;

   //G4cout<<"[P][A]:theta:"<<_theta<<",energy:"<<energy<<",index:"<<_index<<"("<<pos.x()<<","<<pos.y()<<","<<pos.z()<<")"<<G4endl;
    //G4cout<<"[P]:after|cos:"<<sumCosEnergy[_index]<<",en:"<<sumEnergy[_index]<<G4endl;


}
void CellManager::addEnergy(G4double energy, G4int arrayPos, G4ThreeVector momentum, bool isBeta)
{
	//TODO: Add a constant declaration file
	/// IDs
	const int _frontWindowID=3;
	const int _backWindowID=9;
	const int _else=0;
	///



    G4cout<<"[A]AddEnergy:"<<energy<<", arrayPos"<<arrayPos<<G4endl;
    G4double _theta=0;
    G4double _phi=momentum.phi();

    if(_phi>=0)
    {
        if(_phi>=(pi/2))
            _theta=_phi-pi/2;
        else
            _theta=(pi/2)-_phi;
    }
    else
    {
        if(_phi<=(-1*pi/2))
        {
            _theta=(pi/2)+pi+_phi;

        }
        else
        {
            _theta=(pi/2)-_phi;

        }
    }
    if(!isBeta){
    	if(arrayPos==_else)sumNonCellEnergy=+energy;
		if(arrayPos==_frontWindowID)
		{
			sumCosEnergyFrontWindow+=cos(_theta)*energy;
			sumEnergyFrontWindow+=energy;
		}
		if(arrayPos==_backWindowID)
		{
			sumCosEnergyBackWindow+=cos(_theta)*energy;
			sumEnergyBackWindow+=energy;
		}
		arrayPos-=10;
		if(arrayPos<0)
		{
			G4cout<<"[T]Sorry out of cell range"<<G4endl;
			return;
		}
		sumCosEnergy[arrayPos]+=cos(_theta)*energy;
		sumEnergy[arrayPos]+=energy;
    }
    else{
    	if(arrayPos==_else)sumNonCellBEnergy=+energy;
		if(arrayPos==_frontWindowID)
		{
			sumCosBEnergyFrontWindow+=cos(_theta)*energy;
			sumBEnergyFrontWindow+=energy;
		}
		if(arrayPos==_backWindowID)
		{
			sumCosBEnergyBackWindow+=cos(_theta)*energy;
			sumBEnergyBackWindow+=energy;
		}
		arrayPos-=10;
		if(arrayPos<0)
		{
			G4cout<<"[T]Sorry out of cell range"<<G4endl;
			return;
		}
		sumCosBEnergy[arrayPos]+=cos(_theta)*energy;
		sumBEnergy[arrayPos]+=energy;
		print();
	}

}



