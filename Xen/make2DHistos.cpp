#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "TCanvas.h"
#include "TH2.h"
#include "TFile.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

using namespace std;

const int _nHistos=1;

TH2* _histos[_nHistos];

void _printHistos()
{
    for(int i=0;i<_nHistos;i++)
        _histos[i]->Write();

}
void _initHistos()
{
    _histos[0]= new TH2F("beta","Background signal from Betas",16,0,16,9,0,9);
}

///
/// Returns the Colum in the position nCol (zero index)
string _getColumn(string row, int nCol)
{
    istringstream iss(row);
    vector<string> tokens;
    copy(istream_iterator<string>(iss),
         istream_iterator<string>(),
         back_inserter<vector<string> >(tokens));
    if(tokens.size()>nCol)
        return tokens[nCol];
    else
        return "";

}
/// Returns if the string starts with the specified prefix
bool _startsWith(string line, string prefix)
{
   return !line.compare(0, prefix.size(), prefix);
}
string _getFileName(int argc, char **argv)
{
	string _res="";
	while(getopt (argc, argv, "f:")!=-1)
	{
		std::string __s(optarg);
		_res=__s;
	}
	cout<<"FileName:"<<_res<<endl;
	return _res;
}
int main(int argc,char **argv)
{
    TFile* outfile = new TFile("BetaBackground.root", "RECREATE");
    _initHistos();

	std::ifstream in_file;

	in_file.open(_getFileName(argc,argv));

	//read file
	int i=0;
	while(!in_file.eof())
	{
		string _line="";

		getline(in_file, _line);
		if(_line.compare(""))
		for(int j=0; j<16;j++){
			string _column=_getColumn(_line,j);
			cout<<_column<<"	";
			if(_column.compare(""))
				_histos[0]->Fill(j,i,std::stod(_column));
		}
		cout<<endl;
		//cout<<"line:"<<_line<<endl;
		i++;

	}
    in_file.close();
    _printHistos();
    outfile->Close();
    return 0;
}
