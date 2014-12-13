#include "Logger.hh"
#include "CommonsAPI.hh"

using namespace std;

string Logger::_fileName="";
int Logger::_verbose=0;

void Logger::init(string fileName, int verbose)
{
	_verbose=verbose;
	_fileName=fileName;
	_writeLn("{"+CommonsAPI::getDateTime()+"}[Start]______________________________");
}
void Logger::_writeLn(string what)
{
	if(_verbose==118)
		CommonsAPI::appendLine(_fileName, what);
}
void Logger::log(LogType type,  string source,string message)
{
	string _type="INFO";
	if(type==LogType::ERROR)_type="ERROR";
	_writeLn("{"+CommonsAPI::getDateTime()+"}["+_type+"]("+source+")"+message);
}
void Logger::log(string source, string message)
{
	log(LogType::INFO,source,message);
}
void Logger::log(string message)
{
	log(LogType::INFO,"DEF",message);
}
