

if [ $# -ge 1 ];then
	if [[ $1 == "-d" ]]; then
		echo "trying to download samtools from http://cdnetworks-us-1.dl.sourceforge.net/project/samtools/samtools/0.1.12/samtools-0.1.12a.tar.bz2"	fi
		wget http://cdnetworks-us-1.dl.sourceforge.net/project/samtools/samtools/0.1.12/samtools-0.1.12a.tar.bz2
		tar -xvf samtools-0.1.12a.tar.bz2
		mv samtools-0.1.12a samtools
	fi
fi


if [ ! -e samtools ]; then
  echo "samtools not exist. downlaod samtools and put inside this folder. abort. Or automatically download a version dated 2/17/2011 by specifying -d when running this make.sh"
  exit
fi

if [ ! -e samtools/libbam.a ]; then
	echo "samtools not built. Building will proceed now"
	cd samtools
	make
	cd ..
fi

if [ ! -e samtools/libbam.a ]; then 
	echo "samtools build failed. abort"
	exit
fi

g++ -Wall -O3 -o clustergram -I./samtools/ -L./samtools/ -lbam -lm -lz clustergram_main.cpp samtools/libbam.a AdvGetOptCpp/AdvGetOpt.cpp
g++ -Wall -O3 -o bam2Wig -I./samtools/ -L./samtools/ -lbam -lm -lz bamwig_main.cpp samtools/libbam.a AdvGetOptCpp/AdvGetOpt.cpp