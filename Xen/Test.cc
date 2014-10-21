#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "RandomGEN.hh"

using namespace std;
int main()
{
    for(int i=0;i<10000;i++)
    {

        std::normal_distribution<> _distribution(0,1);
        std::random_device _rd;
        std::mt19937 _gen(_rd());


        double x=_distribution(_gen);
        double y=_distribution(_gen);
        double z=_distribution(_gen);

        double _length=sqrt(pow(x,2)+pow(y,2)+pow(z,2));

        x=x/_length;
        y=y/_length;
        z=z/_length;

        //START Writting to file
        std::ofstream ofs ("test.dat", std::ofstream::app);
        ofs <<x<<","<<y<<","<<z<<std::endl;
        ofs.close();
        //END Writting
    }
}
