#include "MomentumManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "RandomGEN.hh"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>
using namespace std;
double MomentumManager::_theta=0;
double MomentumManager::_phi=0;
double MomentumManager::_x=0;
double MomentumManager::_y=0;
double MomentumManager::_z=0;

//double* MomentumManager::_distribution;
//std::uniform_real_distribution<> MomentumManager::_distribution(0,1);
//std::random_device MomentumManager::_rd;
//std::mt19937 MomentumManager::_gen(_rd());
void MomentumManager::init()
{
    RandomGEN::init(0,1, NORMAL);
    _x=_y=_z=0;
}
double MomentumManager::getRandomNumber()
{
//    std::uniform_real_distribution<> _dist(0,2*3);
//    std::random_device _randDev;
//    std::mt19937 _gene(_randDev());
//    double _rand=sin(_distribution(_gen)*2*pi);
//    //_rand=_dist(_gene);
    double _rr=RandomGEN::getRnd();
    G4cout<<"[MM]"<<_rr<<G4endl;
    return _rr;
}
double* MomentumManager::getLastMomentumDirection()
{
    static double _res[3] = {0};
//    _x=RandomGEN::getRnd();
//    _y=RandomGEN::getRnd();
//    _z=RandomGEN::getRnd();
//_x=0;
//_y=0;
//_z=0;
    _res[0]=_x;
    _res[1]=_y;
    _res[2]=_z;

    return _res;

}
double* MomentumManager::getRandomMomentumDirection()
{
    static double _res[3] = {0};
    _x=RandomGEN::getRnd();
    _y=RandomGEN::getRnd();
    _z=RandomGEN::getRnd();



    double _length=sqrt(pow(_x,2)+pow(_y,2)+pow(_z,2));
    _x=_x/_length;
    _y=_y/_length;
    _z=_z/_length;

    cout<<"normalized, x:"<<_x<<", y:"<<_y<<", z:"<<_z<<endl;

    //START Writting to file
    std::ofstream ofs ("normVec.dat", std::ofstream::app);
    ofs <<_x<<","<<_y<<","<<_z<<std::endl;
    ofs.close();
    //END Writting

    _res[0]=_x;
    _res[1]=_y;
    _res[2]=_z;
    return _res;

}
