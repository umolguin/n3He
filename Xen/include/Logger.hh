#ifndef Logger_h
#define Logger_h 1

#include "globals.hh"
#include "G4Track.hh"

using namespace std;

enum LogType{
	INFO, ERROR
};

class Logger
{
    public:

        static void init(string _fileName, int verbose);
        static void setOrigin(G4double x, G4double y, G4double z);
        static void log(LogType type, string source,string message);
        static void log(string source, string message);
        static void log(string message);//type=INFO, source=default

    private:
        static void _writeLn(string what);
        static string _fileName;
        static int _verbose;



};

#endif
