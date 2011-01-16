#include "AdvGetOpt.h"
#include <string.h>
#include <map>
#include <iostream>
#include <fstream>
using namespace std;

#define FILE_BUFFER_LENGTH 10240
#define FILE_ARGS_SPLIT "\t"

string preprocessFileLoadableArgs_inner_formkeystring(char* key,string extraPrefix="")
{
	
string prefix="";
	
if(key[0]!='-')
{
	if (strlen(key)==1)
	{
		prefix="-";
	}
	else
	{
		prefix="--";
	}
}
	
char* ckeyptr=key;
while(*ckeyptr!='\0'){
	if(*ckeyptr=='_')
		*ckeyptr='-';
	else if(*ckeyptr==' ')
		*ckeyptr='-';
	
	ckeyptr++;
}
	
	return prefix+extraPrefix+key;
}

string argv2vectorOfString(vector<string>& vectorOfString,int argc,char* argv[])
{
	for(int i=1;i<argc;i++)
		vectorOfString.push_back(string(argv[i]));		
	
	return string(argv[0]);
}

bool preprocessFileLoadableArgs(vector<string>& args,vector<string>& processedArgs)

{
	//cerr<<"size="<<args.size()<<endl;
	char buffer[FILE_BUFFER_LENGTH];
	
	int i=0;
	processedArgs.clear();
	int argLength=args.size();
	while(i<argLength)
	{
		string& avalue=args[i];
		if (avalue=="--@import-args")
		{
			if(i<argLength-1)
			{
				string& filename=args[i+1];
				ifstream fil(filename.c_str());
				while(fil.good())
				{	
					strcpy(buffer,"");
					fil.getline(buffer,FILE_BUFFER_LENGTH);
					if(strlen(buffer)<1)
						break; //no more to load
					
					char* pch;
					pch=strtok(buffer,FILE_ARGS_SPLIT);
					int j=0;
					while(pch!=NULL)
					{
						if(j==0){
							string argname=preprocessFileLoadableArgs_inner_formkeystring(pch);
							processedArgs.push_back(argname);
						}else {
							processedArgs.push_back(string(pch));
						}

						j++;
						pch=strtok(NULL,FILE_ARGS_SPLIT);
					}
					
					
				}
				
				fil.close();
			}else {
				cerr<<"error: --@import-args not followed by a filename"<<endl;
				return false;
			}
			
			i++;

		}
		else if(avalue=="--@import-cfg")
		{
			//not supported yet
			cerr<<"--@import-cfg not supported yet"<<endl;
			i++;
			return false;
		}else{
			processedArgs.push_back(avalue);
		}
		
		i++;
	}
				
	return true;
	
}



bool getopt(vector<OptStruct>& out_opts,vector<string>& out_args,vector<string> &in_args, string options,vector<string>* long_options)
{
	//build arg list	
	map<string,bool> optionlist; //optionlist[optname]=<whether it consumes an argvalue>
	//process short option first "ab:c" -> "-a",F "-b",T "c",F
	
	unsigned int optionLength=options.length();
 	for(unsigned int i=0;i<optionLength;i++){
		string optname=string("-")+options.substr(i,1);
		bool requireArgValue=false;
		if(i<optionLength-1 && options[i+1]==':'){
			requireArgValue=true;
			i++;
		}
		optionlist.insert(map<string,bool>::value_type(optname,requireArgValue));
		
	}
	
	//now process long options "a-1=", "a-2" => "--a-1",T "--a-2",F
	if(long_options)
	{
		for(vector<string>::iterator i=long_options->begin();i!=long_options->end();i++)
		{
			string optname=string("--")+(*i);
			int optnameLength=optname.length();
			bool requireArgValue=false;
			if(optname[optnameLength-1]=='='){
				requireArgValue=true;
				optname.erase(optnameLength-1);
			}
			
			optionlist.insert(map<string,bool>::value_type(optname,requireArgValue));
		}
	}
	
	int i=0;
	int argLength=in_args.size();
	//cerr<<"size="<<argLength<<endl;
	//now process in_args
	while(i<argLength){
		string& avalue=in_args[i];
		if(avalue[0]=='-')
		{
			map<string,bool>::iterator findOptI=optionlist.find(avalue);
			if (findOptI==optionlist.end()) {
				cerr<<"Error: Option "<<avalue<<" not found"<<endl;
				return false;
			}
			
			string& optname=avalue;
			string optvalue="";
			
			bool requireArgValue=findOptI->second;
			
			if(requireArgValue)
			{
				if(i>=argLength-1){
					cerr<<"Error: Option "<<optname<<" not followed by a required optvalue"<<endl;
					return false;
				}
				optvalue=in_args[i+1];
				i++;
			}
			
			out_opts.push_back(OptStruct(optname,optvalue));
		}else{
			out_args.push_back(avalue);
		}
		i++;
	}
		   
	return true;
	
}
