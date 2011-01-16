#include <iostream>
using namespace std;

#include <sam.h>

int fetch_func(const bam1_t *b, void *data)  
{  
	cout<<bam_format1((bam_header_t*)data,b)<<endl;
	return 0;
}  

int main(int argc,char*argv[])
{
	if(argc<3){
		cerr<<"Usage: samTry sortedBamFile region[chr:start-end]"<<endl;
		return 1;
	}
	
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