#!/usr/bin/env python

from sys import *
from BedSeqUtil import *

def printUsageAndExit(programName):
	print >> stderr,"Usage:",programName,"younglabfile chrsizelist seqDir readLength > samfile"
	exit()

if __name__=="__main__":
	programName=argv[0]
	args=argv[1:]
	
	try:
		ylf,chrsizes,seqDir,readLength=args
	except:
		printUsageAndExit(programName)
	
	bedSeqClient=BedSeqClient(seqDir,"bed")
	
	seqQual="<"
	readLength=int(readLength)
	seqQualString=seqQual*readLength
	
	#first write header
	print >> stdout,"@HD\tVN:1.0"
	
	chrsiz=dict()
	
	fil=open(chrsizes)
	for lin in fil:
		fields=lin.rstrip("\r\n").split("\t")
		chrsiz[fields[0]]=int(fields[1])
		print >> stdout,"@SQ\tSN:"+fields[0]+"\tLN:"+fields[1]
		
	fil.close()
	
	curChr=""
	readNo=1
	
	fil=open(ylf)
	for lin in fil:
		lin=lin.strip()
		if len(lin)<1:
			continue
		
		if lin[0]=="#":
			continue
		if lin[0]==">":				
			curChr="chr"+lin[1:]
			curChromSize=chrsiz[curChr]
			continue
		
		ylfpos=int(lin)
		if ylfpos<0:
			pos0=-ylfpos-readLength
		else:
			pos0=ylfpos
		
		if pos0+readLength>=curChromSize:
			continue #skip this read!
		
		
		#now get seq
		readSeq=bedSeqClient.getBedSeq(curChr+"\t"+str(pos0)+"\t"+str(pos0+readLength)).split("\t")[3].upper()
		readName="read_"+str(readNo)
		
		fieldsOut=[readName,"0",curChr,str(pos0+1),"255",str(readLength)+"M","=","0","0",readSeq,seqQualString]
		print >> stdout,"\t".join(fieldsOut)
		readNo+=1
		
	fil.close()
	bedSeqClient.close()