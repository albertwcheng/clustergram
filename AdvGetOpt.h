#ifndef _ADV_GET_OPT_H
#define _ADV_GET_OPT_H

#include <vector>
#include <string>
using namespace std;

class OptStruct{
public:
	string opname;
	string opvalue;
	inline OptStruct(const string& _opname,const string& _opvalue):opname(_opname),opvalue(_opvalue){}
};

string argv2vectorOfString(vector<string>& vectorOfString,int argc,char* argv[]); //return programName
bool preprocessFileLoadableArgs(vector<string>& args,vector<string>& processedArgs);
bool getopt(vector<OptStruct>& out_opts,vector<string>& out_args,vector<string> &in_args, string options,vector<string>* long_options=NULL);

#endif /*_ADV_GET_OPT_H*/