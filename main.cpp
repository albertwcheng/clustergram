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
#include <fstream>
#include <set>
using namespace std;

#include <sam.h>

#include "AdvGetOptCpp/AdvGetOpt.h"



int fetch_func(const bam1_t *b, void *data)  
{  
	cout<<bam_format1((bam_header_t*)data,b)<<endl;
	return 0;
}  


class fetch_count_data{
public:
	bam_header_t* bamHeader;
	int count;
	int start0Constraint;
	int end1Constraint;
	inline void resetCount(){
		count=0;
	}
	inline fetch_count_data(bam_header_t* _bamHeader,int _start0Constraint=-1,int _end1Constraint=-1,int _count=0):bamHeader(_bamHeader),count(_count),start0Constraint(_start0Constraint),end1Constraint(_end1Constraint){}
	
};

int fetch_count(const bam1_t *b, void *data)
{
	fetch_count_data* fdata=(fetch_count_data*)data;
	if(fdata->start0Constraint!=-1 && b->core.pos<fdata->start0Constraint)
		return 0;
	if(fdata->end1Constraint!=-1 && b->core.pos>=fdata->end1Constraint)
		return 0;
	fdata->count++;
	return 0;
}

#define TABS "\n  "

void outArgsHelp(string argname,string help){
	cerr<<argname<<TABS<<help<<endl;
}



#define DEFAULT_BIN_SIZE "100"
#define DEFAULT_UPSTREAM_NUM_BINS "40"
#define DEFAULT_DOWNSTREAM_NUM_BINS "40"
#define DEFAULT_UNNAMEDFEATUREPREFIX "feature"
				
typedef struct __clustergram_opts{
	string bamfile;
	string bedfile;
	int binSize;
	bool useGenomicStrand;
	bool readPerMillion;
	int upstreamNumBins;
	int downstreamNumBins;
	bool alignStartBinCenter;
	string unnamedFeaturePrefix;
	bool useBedStartOnly;
	bool bedStartIsBase1;
	bool useTranscriptStartOnly;
	int totalNumOfReads;
	bool getTotalNumOfReadsFromBam;
} clustergram_opts;


int countTotalNumOfReadsInBam(string bamfile){
	int total=0;
	samfile_t* bf=samopen(bamfile.c_str(),"rb",0);
	
	
	
	if(!bf){
		cerr<<"bam file "<<bamfile<<" cannot be open for counting"<<endl;
		return -1;
	}
	
	//now do the counting
	bam1_t *dummy=bam_init1();
	
	
	while(samread(bf,dummy)>=0){
		total++;
		if(total%1000000==1){
			cerr<<"passing through read "<<total<<endl;
		}
	}
	
	bam_destroy1(dummy);
	samclose(bf);
	return total;
}


#define STRAND_FORWARD '+'
#define STRAND_REVERSE '-'

class BedEntry{
public:
	string chrom;
	int start0;
	int end1;
	string name;
	double score;
	char strand;
	inline static string getGenomeBrowserCoordinate(const string&chrom,int start0,int end1){
		char buf[256];
		sprintf(buf,"%s:%d-%d",chrom.c_str(),start0+1,end1);
		return buf;
	}
	
	inline string getGenomeBrowserCoordinate(){
		return BedEntry::getGenomeBrowserCoordinate(chrom,start0,end1);
	}
	
	inline void reset(){
		chrom="";
		start0=-1;
		end1=-1;
		name="";
		score=0.0;
		strand=STRAND_FORWARD;
		
	}
	
	inline int getLength(){
		return end1-start0;
	}
	inline void parseBedLine(const string& line,const string& defaultName){
		reset();
		char *tmp=new char[line.length()+1];
		strcpy(tmp,line.c_str());
		
		char* pch;
		pch=strtok(tmp,"\t");
		chrom=pch;
		int fieldNo=0;
		while(pch!=NULL){
			switch (fieldNo) {
				case 1:
					//start0
					start0=atoi(pch);
					break;
				case 2:
					//end1
					end1=atoi(pch);
					break;
				case 3:
					//name
					name=pch;
					break;
				case 4:
					//score
					score=atof(pch);
					break;
				case 5:
					strand=pch[0];
					break;

				default:
					break;
			}
			pch=strtok(NULL,"\t");
			fieldNo++;
		}
		delete[] tmp;
		
		if(end1==-1){
			//not set
			end1=start0+1;
		}
		
		if(name==""){
			//not set
			name=defaultName;
		}
	}
	
	inline BedEntry(){
		reset();
	}
	
	inline BedEntry(const string& line,const string& defaultName){
		parseBedLine(line,defaultName);
	}
	
};

#define FILE_BUFFER_SIZE 10240

void output_vector(ostream& os,vector<string>& voutput,string delimiter){	
	vector<string>::iterator i=voutput.begin();
	if(i==voutput.end()){
		return;
	}
	os<<*i;
	i++;
	while(i!=voutput.end()){
		os<<delimiter<<*i;
		i++;
	}
	os<<endl;
}

int fetch_count_region(samfile_t* bamfile,bam_index_t*idx,const string&chrom, int start0, int end1,bool constrainStart){
	string region=BedEntry::getGenomeBrowserCoordinate(chrom,start0,end1);
	int ref_id;
	int begin;
	int end;
	fetch_count_data counter(bamfile->header,(constrainStart?start0:-1));
	
	bam_parse_region(bamfile->header,region.c_str(),&ref_id,&begin,&end);
	bam_fetch(bamfile->x.bam,idx,ref_id,begin,end,&counter,fetch_count);
	
	return counter.count;
}


void fetch_bin_values(vector<double>& voutput,samfile_t* bamfile,bam_index_t*idx,const string&chrom, int ultimateStart0, int ultimateNumBins, int binSize,double normalization){
	int binStart0=ultimateStart0;
	for(int i=0;i<ultimateNumBins;i++){
		int countOfThisBin=fetch_count_region(bamfile,idx,chrom,binStart0,binStart0+binSize,true);
		voutput.push_back(countOfThisBin/normalization);
		binStart0+=binSize;
											  
	}
}

inline bool isChromInBam(set<string>& chromsInBam,const string& chrom){
	return chromsInBam.find(chrom)!=chromsInBam.end();
}


int runClustergram(clustergram_opts& clustergramOpts){
	

	
	if(clustergramOpts.getTotalNumOfReadsFromBam){
	   cerr<<"couting total number of reads"<<endl;
	   clustergramOpts.totalNumOfReads=countTotalNumOfReadsInBam(clustergramOpts.bamfile);
		cerr<<"total num of reads="<<clustergramOpts.totalNumOfReads<<endl;
	}
	
	cerr<<"load bamfile "<<clustergramOpts.bamfile<<endl;
	
	samfile_t* bamfile=samopen(clustergramOpts.bamfile.c_str(),"rb",0);
	if(!bamfile)
	{
		cerr<<"Fail to open BAM file "<<clustergramOpts.bamfile<<endl;
		return 1;
	}
	
	cerr<<"load bamfile index"<<endl;
	
	bam_index_t*idx=bam_index_load(clustergramOpts.bamfile.c_str());
	if(!idx){
		cerr<<"Fail to load index\n"<<endl;
		return 1;
	}
	

	set<string> chromsInBam;
	
	for(int i=0;i<bamfile->header->n_targets;i++){
		cerr<<"bam file has reference "<<bamfile->header->target_name[i]<<endl;
		chromsInBam.insert(bamfile->header->target_name[i]);
	}
	
	int unamedNo=1;
	
	cerr<<"load bedfile "<<clustergramOpts.bedfile<<endl;
	
	ifstream bedfil(clustergramOpts.bedfile.c_str(),ifstream::in );
	
	if(!bedfil.good())
	{
		cerr<<"bed file "<<clustergramOpts.bedfile<<" not good"<<endl;
		return false;
	}
	
	char buff[FILE_BUFFER_SIZE];
	
	double normalization=1.0;
	if(clustergramOpts.readPerMillion){
		normalization=1e6/double(clustergramOpts.totalNumOfReads);
	}
	
	
	int lino=0;
	
	bool firstOutputLine=true;
	
	int bedEntryNo=0;
	
	
	while(bedfil.good()){
		lino++;
		strcpy(buff,"");
		bedfil.getline(buff,FILE_BUFFER_SIZE);
		if(strlen(buff)==0)
			continue;
		
		if(string(buff).substr(0,5)=="track"){
			continue;
		}
		
		BedEntry bedEntry(buff,"");
		if(bedEntry.name==""){
			//unnamed
			//buff is now useless, so use it
			sprintf(buff,"%s%d",clustergramOpts.unnamedFeaturePrefix.c_str(),unamedNo);
			bedEntry.name=buff;
			unamedNo++;
		}
		
		if(clustergramOpts.bedStartIsBase1){
			bedEntry.start0--;
		}
		
		if(clustergramOpts.useBedStartOnly){
			bedEntry.end1=bedEntry.start0+1;
		}
		
		if(clustergramOpts.useTranscriptStartOnly){
			if(bedEntry.strand==STRAND_REVERSE){
				bedEntry.start0=bedEntry.end1-1;
			}else{
				
			}
			
			bedEntry.end1=bedEntry.start0+1;
		}
		
		
		bedEntryNo++;
		
		if(bedEntryNo%10000==1){
			cerr<<"processing bed entry "<<bedEntryNo<<" ..."<<endl;
		}
		
		if(!isChromInBam(chromsInBam, bedEntry.chrom)){
			cerr<<"bedEntry "<<bedEntry.name<<" @"<<bedEntry.chrom<<" rejected becuase "<<bedEntry.chrom<<" not found in bam file"<<endl;
			continue;
		}
		
		//now the real deal here
		//HERE:
		
		int ultimateStart0;
		int ultimateNumBins;
		bool reverseData=false;
		
		vector<string> voutput;
		
		if(clustergramOpts.alignStartBinCenter){
			int center0=bedEntry.getLength()/2+bedEntry.start0;
			int startBinStart0=center0-clustergramOpts.binSize/2;
			if(bedEntry.strand==STRAND_REVERSE && !clustergramOpts.useGenomicStrand){
				//use reverse
				reverseData=true;
				ultimateStart0=startBinStart0-clustergramOpts.downstreamNumBins*clustergramOpts.binSize;
				
			}else {
				//assume forward
				ultimateStart0=startBinStart0-clustergramOpts.upstreamNumBins*clustergramOpts.binSize;
			}

			ultimateNumBins=clustergramOpts.upstreamNumBins+clustergramOpts.downstreamNumBins+1;
			
		}
		else {
			
			if(bedEntry.strand==STRAND_REVERSE && !clustergramOpts.useGenomicStrand){
				reverseData=true;
				ultimateStart0=bedEntry.start0-(clustergramOpts.downstreamNumBins-1)*clustergramOpts.binSize;
			}
			else{
				ultimateStart0=bedEntry.start0-clustergramOpts.upstreamNumBins*clustergramOpts.binSize;
			}
			ultimateNumBins=clustergramOpts.upstreamNumBins+clustergramOpts.downstreamNumBins;
		}
		
		if(firstOutputLine){
			voutput.push_back("GENEID");
			voutput.push_back("NAME");
			
			int binStart0=-clustergramOpts.upstreamNumBins*clustergramOpts.binSize;
			for(int i=0;i<ultimateNumBins;i++)
			{
				sprintf(buff,"%d",binStart0);
				voutput.push_back(buff);
				binStart0+=clustergramOpts.binSize;
			}
			
			output_vector(cout,voutput,"\t");
			
			voutput.clear();
			firstOutputLine=false;
		}
		
		vector<double> binvalues;
		
		fetch_bin_values(binvalues,bamfile,idx,bedEntry.chrom, ultimateStart0, ultimateNumBins, clustergramOpts.binSize,normalization);
		
		sprintf(buff,"%s_%d",bedEntry.chrom.c_str(),bedEntry.start0+1);
		voutput.push_back(buff); //chr1_12512
		voutput.push_back(bedEntry.name); //GeneName
		
		
		if(bedEntry.strand==STRAND_REVERSE && !clustergramOpts.useGenomicStrand){
			for(vector<double>::reverse_iterator i=binvalues.rbegin();i!=binvalues.rend();i++)
			{
				sprintf(buff,"%.5f",*i);
				voutput.push_back(buff);
			}
		}
		else{
			
			for(vector<double>::iterator i=binvalues.begin();i!=binvalues.end();i++)
			{
				sprintf(buff,"%.5f",*i);
				voutput.push_back(buff);
			}
		}
		
		output_vector(cout,voutput,"\t");
		
	}
	
	
	bedfil.close();
	
	

	
	
	
	
	//end: clean up
	bam_index_destroy(idx);
	samclose(bamfile);
	return 0;
	
}

void printUsage(string programName){
	
	cerr<<"Usage: "<<programName<<" [options]"<<endl;
	cerr<<"preprocessor options:"<<endl;
	outArgsHelp("--@import-args","filename load arguments from tab delimited file");
	cerr<<endl;
	cerr<<"required options:"<<endl;
	outArgsHelp("--target-read-file bamfile","specify the bam file");
	outArgsHelp("--bedfile befile","specify the bed file for the regions to plot");
	cerr<<endl;
	cerr<<"optional options:"<<endl;
	outArgsHelp("--reads-per-million","output read per millions. Default is not normalization");
	outArgsHelp("--total-num-reads numReads","total number of mapped reads. Required if --reads-per-million is on");
	outArgsHelp("--use-genomic-strand","use genomic strand to define upstream and downstream instead of feature strand specified in bed file. Default is to use feature strand");
	outArgsHelp("--bin-size binSize","specify bin size. Default is 100");
	outArgsHelp("--upstream-num-bins binNum","specify the number of bin upstream. Default is 40");
	outArgsHelp("--downstream-num-bins binNum","specify the number of bin downstream. Default is 40");
	outArgsHelp("--align-start-bin-center","specify whether to align the start bin to center. Default is not");
	outArgsHelp("--unamed-feature-prefix prefix","specify the prefix for the name of an unnamed feature from bed file. Default is feature");
	outArgsHelp("--use-bed-start","Use only the start0 field of bed file assuming end1=start0+1. Default is no");
	outArgsHelp("--bed-start-is-base1","Start field of bed file is encoded in base 1 in contrary to the standard. Default is no");
	outArgsHelp("--use-transcript-start-only","Use only the transcript start, dependent on the strand. if strand is postive, start=start0=field2ofbed; if strand is negative, start=end1-1=field3ofbed-1. Default is no");
	outArgsHelp("--get-total-number-of-reads-from-bam","get the total number of reads from bam file. Warning: takes some time. Default is no");
}

int main(int argc,char*argv[])
{
	
	vector<string> long_options;
	//required:
	long_options.push_back("target-read-file=");
	long_options.push_back("bedfile=");
	vector<string> required_opts=long_options;
	
	//optional:
	long_options.push_back("reads-per-million");
	long_options.push_back("use-genomic-strand");
	long_options.push_back("bin-size=");
	long_options.push_back("upstream-num-bins=");
	long_options.push_back("downstream-num-bins=");
	long_options.push_back("align-start-bin-center");
	long_options.push_back("unamed-feature-prefix");
	long_options.push_back("use-bed-start");
	long_options.push_back("bed-start-is-base1");
	long_options.push_back("use-transcript-start-only");
	long_options.push_back("total-num-reads=");
	long_options.push_back("get-total-number-of-reads-from-bam");
	
	map<string,string> optmap;

	EasyAdvGetOptOut argsFinal=easyAdvGetOpt(argc,argv,"",&long_options);
	if(argsFinal.success){
		argsFinal.print(cerr);
		
		//cerr<<"here"<<endl;
		parseOptsIntoMap(argsFinal.opts,optmap);
		if(!checkRequiredOpts(optmap,required_opts)){
			printUsage(argsFinal.programName);
			return 1;
		}
		
	}
	else{
		printUsage(argsFinal.programName);
		return 1;
	}
	
	clustergram_opts clustergramOpts;
	clustergramOpts.totalNumOfReads=-1;
	
	clustergramOpts.bamfile=getOptValue(optmap,"--target-read-file");
	clustergramOpts.bedfile=getOptValue(optmap,"--bedfile");
	clustergramOpts.binSize=atoi(getOptValue(optmap,"--bin-size",DEFAULT_BIN_SIZE).c_str());
	clustergramOpts.useGenomicStrand=hasOpt(optmap,"--use-genomic-strand");
	clustergramOpts.readPerMillion=hasOpt(optmap,"--reads-per-million");
	clustergramOpts.upstreamNumBins=atoi(getOptValue(optmap,"--upstream-num-bins").c_str());
	clustergramOpts.downstreamNumBins=atoi(getOptValue(optmap,"--downstream-num-bins").c_str());
	clustergramOpts.alignStartBinCenter=hasOpt(optmap,"--align-start-bin-center");
	clustergramOpts.unnamedFeaturePrefix=getOptValue(optmap,"--unamed-feature-prefix",DEFAULT_UNNAMEDFEATUREPREFIX);
	clustergramOpts.useBedStartOnly=hasOpt(optmap,"--use-bed-start");
	clustergramOpts.bedStartIsBase1=hasOpt(optmap,"--bed-start-is-base1");
	clustergramOpts.useTranscriptStartOnly=hasOpt(optmap,"--use-transcript-start-only");
	clustergramOpts.getTotalNumOfReadsFromBam=hasOpt(optmap,"--get-total-number-of-reads-from-bam");
	if(clustergramOpts.readPerMillion){
		if(!hasOpt(optmap,"--total-num-reads") && !	clustergramOpts.getTotalNumOfReadsFromBam){
			cerr<<"total number of reads not specifed when --reads-per-million is on"<<endl;
			return 1;
		}
		
	}
	
	if(hasOpt(optmap,"--total-num-reads"))
		clustergramOpts.totalNumOfReads=atoi(getOptValue(optmap,"--total-num-reads").c_str());
	
	return runClustergram(clustergramOpts);
	
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