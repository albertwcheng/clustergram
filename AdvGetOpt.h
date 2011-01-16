/***************************************************************************
 Copyright 2011 Wu Albert Cheng <albertwcheng@gmail.com>
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 *******************************************************************************/ 


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