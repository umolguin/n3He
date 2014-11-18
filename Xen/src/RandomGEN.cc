#include "RandomGEN.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>

double RandomGEN::_maximum=0;
double RandomGEN::_minimum=0;
DistributionType RandomGEN::_distType;
//double* RandomGEN::_distribution;
//std::normal_distribution<> RandomGEN::_normalDist(0,0);
//std::uniform_real_distribution<> RandomGEN::_uniformRDist(0,0);
std::random_device RandomGEN::_rd;
std::mt19937 RandomGEN::_gen(_rd());
void RandomGEN::init(double minimum, double maximum, DistributionType distributionType)
{
    _minimum=minimum;
    _maximum=maximum;
    _distType=distributionType;
}

double RandomGEN::getRnd()
{
    if(_distType==UNIFORM)
    {

        std::uniform_real_distribution<> _uniformRDist(_minimum,_maximum);
        double _rnd=_uniformRDist(_gen);
        std::cout<<"rand:"<<_rnd<<std::endl;
        return _rnd;
    }
    if(_distType==NORMAL)
    {

        std::normal_distribution<> _normalDist(_minimum,_maximum);//this is mean and stnd_dev
        double _rnd=_normalDist(_gen);

        return _rnd;
    }
   return 0;//
}
