YLF2Sam.py s_5_sequence.mapview.q50.m1.txt /lab/jaenisch_albert/genomes/mm8/mm8_nr.sizes /lab/jaenisch_albert/genomes/mm8/seq/ 36 > s_5_sequence.mapview.q50.m1.sam
bsub samtools view -bt /lab/jaenisch_albert/genomes/mm8/mm8_nr.sizes -o s_5_sequence.mapview.q50.m1.bam s_5_sequence.mapview.q50.m1.sam
samtools sort s_5_sequence.mapview.q50.m1.bam s_5_sequence.mapview.q50.m1.sorted
samtools index s_5_sequence.mapview.q50.m1.sorted.bam
./clustergram --@import-args setting.txt > s_5_sequence.mapview.q50.m1.cdt