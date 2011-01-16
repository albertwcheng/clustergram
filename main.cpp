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

#include <iostream>
using namespace std;

#include <sam.h>
#include "AdvGetOptCpp/AdvGetOpt.h"

int fetch_func(const bam1_t *b, void *data)  
{  
	cout<<bam_format1((bam_header_t*)data,b)<<endl;
	return 0;
}  

int main(int argc,char*argv[])
{
	/*if(argc<3){
		cerr<<"Usage: samTry sortedBamFile region[chr:start-end]"<<endl;
		return 1;
	}*/
	
	vector<OptStruct> opts;
	vector<string> args;
	vector<string> preprocessed_in_args;
	vector<string> processed_in_args;
	vector<string> long_options;
	long_options.push_back("long-1");
	long_options.push_back("long-2=");
	
	bool success;
	
	string programName=argv2vectorOfString(preprocessed_in_args,argc,argv);
	success=preprocessFileLoadableArgs(preprocessed_in_args,processed_in_args);
	if(!success){
		cerr<<"preprocessing failed"<<endl;
		return 1;
	}
	success=getopt(opts,args,processed_in_args, "ab:c",&long_options);
	if(!success){
		cerr<<"getopt failed"<<endl;
		return 1;
	}
	
	for(vector<OptStruct>::iterator i=opts.begin();i!=opts.end();i++)
	{
		cerr<<"option "<<i->opname<<" : "<<i->opvalue<<endl;
	}
	
	for(unsigned int i=0;i<args.size();i++)
	{
		cerr<<"args "<<i<<" : "<<args[i]<<endl;
	}
	
	return 0;
	
	
	
	
	char *filename=argv[1];
	char *region=argv[2];
	samfile_t* bamfile=samopen(filename,"rb",0);
	if(!bamfile)
	{
		cerr<<"Fail to open BAM file "<<filename<<endl;
		return 1;
	}
	
	bam_index_t*idx=bam_index_load(filename);
	if(!idx){
		cerr<<"Fail to load index\n"<<endl;
		return 1;
	}
	
	int ref_id;
	int begin;
	int end;
	
	bam_parse_region(bamfile->header,region,&ref_id,&begin,&end);
	
	bam_fetch(bamfile->x.bam,idx,ref_id,begin,end,bamfile->header,fetch_func);
	
	bam_index_destroy(idx);
	
	samclose(bamfile);
	
	return 0;
}