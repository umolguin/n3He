#include "CommonsAPI.hh"
#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>
#include <string>
#include <fstream>

using namespace std;

string CommonsAPI::getDateTime()
{
	 time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

	return buf;
}
void CommonsAPI::appendLine(string fileName, string line)
{
	std::ofstream ofs (fileName, ios::app);
	ofs <<line<<endl;
	ofs.close();
}
