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

#define TABS "\n  "

void outArgsHelp(string argname,string help){
	cerr<<argname<<TABS<<help<<endl;
}



#define UNIT_READS 0
#define UNIT_RPKM 1
#define UNIT_RPM 2

#define PHASE_FIXED 0
#define PHASE_FLEXIBLE 1

#define MODE_USE_START 0
#define MODE_PILEUP 1


class bamwig_opts{
public:
	string bamfile;
	int binSize;
	int unitMode;
	int totalNumOfReads;
	bool getTotalNumOfReadsFromBam;
	int phaseMode;
	bool perChrom;
	string perChromDir;
	string lastChrom; //for checking if sorted
	uint32_t binCoord;
	uint32_t lastCoord;
	int countInBin;
	samfile_t* pbamfile;
	ostream* os;
	bool centerBinCoord;
	int countMode;
	
	
	inline int startBin(uint32_t newCoord){
		countInBin=0;
		//now calculate new binCoord according to phase mode
		switch (phaseMode) {
			case PHASE_FLEXIBLE:
				binCoord=newCoord;
				break;
			case PHASE_FIXED:
				binCoord=(newCoord/binSize)*binSize;
				break;
			default:
				cerr<<"Error: unknown phase mode";
				break;
		}
		
		return binCoord;
	}	
	
	inline void resetForNewChrom(const string& _chr,uint32_t newCoord){
		lastCoord=-1;
		lastChrom=_chr;
		startBin(newCoord);
		if(perChrom){
			ofstream* fos=(ofstream*)os;
			if(fos){
				fos->close();
				delete fos;
			}
			
			os=new ofstream((perChromDir+_chr+".wig").c_str());
		}
		
		
		//print starter line for chromosome
		(*os)<<"variableStep chrom="<<_chr<<endl;
	}
	
	
	~bamwig_opts()
	{
		if(perChrom){
			ofstream* fos=(ofstream*)os;
			if(fos){
				fos->close();
				delete fos;
			}
			
		}
		
	}
	
	
	bool needANewBin(uint32_t newCoord){
		//do we need a new bin?
		return newCoord>binCoord+binSize;
		
	}
	
	void printWigLine(){
		if(binCoord>=0 && countInBin>0){
			(*os)<<((centerBinCoord)?(binCoord+binSize/2):binCoord)<<"\t"<<countInBin<<endl;
		}
	}
	
	string getChromName(uint32_t tid){
		return pbamfile->header->target_name[tid];
	}
	
	int getNumOfChroms(){
		return pbamfile->header->n_targets;
	}
	
};


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

/*
 typedef struct { 
 bam1_t *b; 
 int32_t qpos; 
 int indel, level; 
 uint32_t is_del:1, is_head:1, is_tail:1; 
 } bam_pileup1_t;  

 
 
 typedef int ( *bam_pileup_f)(
 uint32_t tid,
 uint32_t pos,
 int n,
 const bam_pileup1_t *pl,
 void *data);  
 

typedef struct { 
    bam1_core_t core; 
    int l_aux, data_len, m_data; 
    uint8_t *data; 
} bam1_t;  

 typedef struct { 
 int32_t tid; 
 int32_t pos; 
 uint32_t bin:16, qual:8, l_qname:8; 
 uint32_t flag:16, n_cigar:16; 
 int32_t l_qseq; 
 int32_t mtid; 
 int32_t mpos; 
 int32_t isize; 
 } bam1_core_t;  
 
 */


int countFunc(uint32_t tid,uint32_t pos,int n,const bam_pileup1_t* pl,void *data){
	bamwig_opts* bamWigOpts=(bamwig_opts*)data;
	//now bamWigOpts remember our current state and options
	
	
	string chrom=bamWigOpts->getChromName(tid);
	if(chrom!=bamWigOpts->lastChrom){ //see if new chrom
		bamWigOpts->resetForNewChrom(chrom,pos);
	}
	
	if(bamWigOpts->needANewBin(pos)){
		bamWigOpts->printWigLine();
		bamWigOpts->startBin(pos);
	}
	
	bamWigOpts->countInBin++;
	
	return 1;
	
}


void sampass(samfile_t* bf,bam_pileup_f func,void* data){
	
	bam1_t* bamt=bam_init1();
	
	bam_pileup1_t pl;
	pl.b=bamt;
	pl.qpos=0;
	pl.indel=0;
	pl.level=0;
	
	while(samread(bf,bamt)>=0){
		
		//now pass to the func
		func(bamt->core.tid,bamt->core.pos,1,&pl,data);
	}
	
	bam_destroy1(bamt);
	
}


int runBamWig(bamwig_opts& bamWigOpts){
	

	
	if(bamWigOpts.getTotalNumOfReadsFromBam){
	   cerr<<"couting total number of reads"<<endl;
	   bamWigOpts.totalNumOfReads=countTotalNumOfReadsInBam(bamWigOpts.bamfile);
		cerr<<"total num of reads="<<bamWigOpts.totalNumOfReads<<endl;
	}
	
	cerr<<"load bamfile "<<bamWigOpts.bamfile<<endl;
	
	bamWigOpts.pbamfile=samopen(bamWigOpts.bamfile.c_str(),"rb",0);
	if(!bamWigOpts.pbamfile)
	{
		cerr<<"Fail to open BAM file "<<bamWigOpts.bamfile<<endl;
		return 1;
	}
	
	cerr<<"load bamfile index"<<endl;
	
	bam_index_t*idx=bam_index_load(bamWigOpts.bamfile.c_str());
	if(!idx){
		cerr<<"Fail to load index\n"<<endl;
		return 1;
	}
	

	set<string> chromsInBam;
	
	for(int i=0;i<bamWigOpts.pbamfile->header->n_targets;i++){
		string chrom=bamWigOpts.getChromName(i);
		cerr<<"bam file has reference "<<chrom<<endl;
		chromsInBam.insert(chrom);
	}
	
	
	
	
	switch (bamWigOpts.countMode) {
		case MODE_USE_START:
			sampileup(bamWigOpts.pbamfile,-1,countFunc,&bamWigOpts);
			break;
		case MODE_PILEUP:
			sampass(bamWigOpts.pbamfile,countFunc,&bamWigOpts);
			break;
		default:
			cerr<<"unknown count mode. abort"<<endl;
			return 1;
			
	} 
	
	//print the last line
	bamWigOpts.printWigLine();
	
	
	
	//end: clean up
	bam_index_destroy(idx);
	samclose(bamWigOpts.pbamfile);
	return 0;
	
}

void printUsage(string programName){
	
	cerr<<"Usage: "<<programName<<" [options] <bamfile> > wigFile"<<endl;
	cerr<<"Preprocessor options:"<<endl;
	outArgsHelp("--@import-args","filename load arguments from tab delimited file");
	cerr<<endl;
	cerr<<"Unit options:"<<endl;
	outArgsHelp("[default] --rpkm-auto","output read per kilobase per million reads. Automatically find the total number of reads");
	outArgsHelp("--rpm-auto","output read per million reads. Automatically find the total number of reads");
	outArgsHelp("--rpkm <totalNumOfReads>","total number of mapped reads. Required if --reads-per-million is on");
	outArgsHelp("--rpm <totalNumOfReads>","output read per millions using the provided total number of reads");
	outArgsHelp("--center-bin-coord","output the center of bin as bin coordinate");
	cerr<<"Operational:"<<endl;
	outArgsHelp("--bin-size <binSize>","specify bin size. Default is 20");
	outArgsHelp("[default] --use-start or --pileup","--use-start: count if start of read fall in bin. --pileup: pile up reads");
	outArgsHelp("[default] --flexible-phase or --fixed-phase","--flexible-phase: start a sequence of bin on reaching the first read of that sequence of bins. --fixed-phase: start a sequence of bin on a fixed interval from chrom start");
	
}

int main(int argc,char*argv[])
{
	
	vector<string> long_options;
	//required:
	vector<string> required_opts=long_options;
	
	//optional:
	long_options.push_back("rpm=");
	long_options.push_back("rpkm=");
	long_options.push_back("rpm-auto");
	long_options.push_back("rpkm-auto");
	long_options.push_back("pile-up");
	long_options.push_back("use-start");
	long_options.push_back("fixed-phase");
	long_options.push_back("flexible-phase");
	long_options.push_back("per-chr=");
	long_options.push_back("bin-size=");
	long_options.push_back("center-bin-coord");
	
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
	
	
	if(argsFinal.args.size()<1){
		cerr<<"No bam file provided"<<endl;
		printUsage(argsFinal.programName);
		return 1;
		
		
	}
	
	
	
	bamwig_opts bamwigOpts;
	
	//defaults:
	bamwigOpts.unitMode=UNIT_RPKM;
	bamwigOpts.getTotalNumOfReadsFromBam=true;
	bamwigOpts.totalNumOfReads=-1;
	bamwigOpts.phaseMode=PHASE_FLEXIBLE;
	bamwigOpts.perChrom=false;
	bamwigOpts.perChromDir="";
	bamwigOpts.countMode=MODE_USE_START;
	//get args
	bamwigOpts.bamfile=argsFinal.args[0];
	bamwigOpts.centerBinCoord=false;
	//get opts
	bamwigOpts.binSize=atoi(getOptValue(optmap,"--bin-size","20").c_str());
	if(hasOpt(optmap,"--rpkm-auto")){
		bamwigOpts.unitMode=UNIT_RPKM;
		bamwigOpts.getTotalNumOfReadsFromBam=true;
	}else if(hasOpt(optmap,"--rpm-auto")){
		bamwigOpts.unitMode=UNIT_RPM;
		bamwigOpts.getTotalNumOfReadsFromBam=true;
	}else if(hasOpt(optmap,"--rpkm")){
		bamwigOpts.unitMode=UNIT_RPKM;
		bamwigOpts.totalNumOfReads=atoi(getOptValue(optmap,"--rpkm").c_str());
	}else if(hasOpt(optmap,"--rpm")){
		bamwigOpts.unitMode=UNIT_RPM;
		bamwigOpts.totalNumOfReads=atoi(getOptValue(optmap,"--rpm").c_str());
	}
	
	
	if(hasOpt(optmap,"--use-start")){
		bamwigOpts.countMode=MODE_USE_START;
	}else if(hasOpt(optmap,"--pile-up")){
		bamwigOpts.countMode=MODE_PILEUP;
	}
	
	
	if(hasOpt(optmap,"--per-chrom")){
		bamwigOpts.perChrom=true;
		bamwigOpts.perChromDir=getOptValue(optmap,"--per-chrom",".");
		bamwigOpts.os=NULL;
	}else {
		bamwigOpts.os=&cout;
	}

	
	if(hasOpt(optmap,"--flexible-phase")){
		bamwigOpts.phaseMode=PHASE_FLEXIBLE;
	}else if(hasOpt(optmap,"--fixed-phase")){
		bamwigOpts.phaseMode=PHASE_FIXED;
	}
	
	if(hasOpt(optmap,"--center-bin-coord")){
		bamwigOpts.centerBinCoord=true;
	}
	
	return runBamWig(bamwigOpts);
	
	return 0;
	
}