#ifndef CommonsAPI_h
#define CommonsAPI_h 1
#include <string>

class CommonsAPI
{
    public:

        static std::string getDateTime();
        static void appendLine(std::string fileName, std::string line);

};

#endif
