

if [ ! -e samtools ]; then
  echo "samtools not exist. downlaod samtools and put inside this folder. abort"
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

g++ -Wall -o clustergram -I./samtools/ -L./samtools/ -lbam -lm -lz clustergram_main.cpp samtools/libbam.a AdvGetOptCpp/AdvGetOpt.cpp
