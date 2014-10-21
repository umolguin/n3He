#ifndef MomentumManager_h
#define MomentumManager_h 1

#include "globals.hh"
#include <map>
#include <random>
#include <cmath>

class MomentumManager
{
    public:
        static double* getRandomMomentumDirection();
        static double* getLastMomentumDirection();
        static double getRandomNumber();
        static void init();
    private:

        static double _theta;
        static double _phi;
        static double _x;
        static double _y;
        static double _z;

        static std::uniform_real_distribution<> _distribution;
        static std::random_device _rd;
        static std::mt19937 _gen;

};
#endif
