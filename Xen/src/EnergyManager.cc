#include "EnergyManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>

double EnergyManager::_maximum=0;
double EnergyManager::_minimum=0;
//double* EnergyManager::_distribution;
std::normal_distribution<> EnergyManager::_distribution(.32,.25);
std::random_device EnergyManager::_rd;
std::mt19937 EnergyManager::_gen(_rd());

double EnergyManager::getEnergy()
{
    //Get wave length
    double _waveLength=-1;

    while(_waveLength<_minimum||_waveLength>_maximum)
    {
        _waveLength=_distribution(_gen);
       // G4cout<<"[W]TryingWaveLength: "<<_waveLength<<", max:"<<_maximum<<", min: "<<_minimum<<G4endl;
    }

   // _waveLength=std::floor(_waveLength * 100 + 0.5)/100;
    G4cout<<"[W]waveLength: "<<_waveLength<<G4endl;
    //TODO: architectually speaking is better if the transformation to energy is done here...
    //but we have to print histograms... SO TODO: change this some day
    return _waveLength;
}

void EnergyManager::init()
{
    _maximum=0;
    _minimum=0;
}

void EnergyManager::printDist()
{
    //Do stuff

    //print
    print();
}

void EnergyManager::setRangeWaveLength(double minimum, double maximum)
{
    _minimum=minimum;
    _maximum=maximum;
}

void EnergyManager::print()
{


//    std::map<int, int> hist;
//    for(int n=0; n<10000; ++n) {
//        //++hist[std::round(d(gen))];
//        std::cout<<"- "<<_distribution (_gen)<<std::endl;
//    }
//    for(auto p : hist) {
//        std::cout <<d.min()<<"/"<< std::fixed << std::setprecision(1) << std::setw(2)
//                  << p.first << ' ' << std::string(p.second/200, '*') << '\n';
//    }
}
