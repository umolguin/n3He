#ifndef RandomGEN_h
#define RandomGEN_h 1

//#include "globals.hh"
#include <map>
#include <random>
#include <cmath>

enum DistributionType{UNIFORM, NORMAL};

class RandomGEN
{
    public:
        static void init(double minimum, double maximum, DistributionType distType);
        static double getRnd();

    private:

        static double _minimum;
        static double _maximum;
        static DistributionType _distType;
        //static double* _distribution;
//        static std::normal_distribution<> _normalDist;
//        static std::uniform_real_distribution<> _uniformRDist;
        static std::random_device _rd;
        static std::mt19937 _gen;

};
#endif
