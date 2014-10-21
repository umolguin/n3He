#ifndef EnergyManager_h
#define EnergyManager_h 1

#include "globals.hh"
#include <map>
#include <random>
#include <cmath>

class EnergyManager
{
    public:
        static void init();
        static void setRangeWaveLength(double minimum, double maximum);
        static double getEnergy();
        static void printDist();

    private:
        static void print();
        static double _minimum;
        static double _maximum;
        //static double* _distribution;
        static std::normal_distribution<> _distribution;
        static std::random_device _rd;
        static std::mt19937 _gen;

};
#endif
