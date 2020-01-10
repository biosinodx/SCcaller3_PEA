CC=gcc
CFLAGS=-g -W -Wall -O3 -finline-functions -D_FILE_OFFSET_BITS=64
GLIBS=-lm -lbam -lz
INCLUDE=-Isamtools
LIBPATH=-Lsamtools
GENERIC_SRC= string.h vector.h file_reader.h file_reader.c hashset.h list.h dna.h sort.h 

all: samtools/libbam.a novoBreak nbasm

samtools/libbam.a:
	(cd samtools; make libbam.a)

novoBreak: $(GENERIC_SRC) filter_kmer.c sr_aln.h sr_aln.c aln_cigar.h simp_asm.h stdaln.h stdaln.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBPATH) -o $@ file_reader.c sr_aln.c stdaln.c filter_kmer.c $(GLIBS)

nbasm: file_reader.c asm_R2.c  localasm_main.c
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -f samtools/*.o samtools/libbam.a *.o *.gcda *.gcno *.gcov gmon.out novoBreak nbasm

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out
