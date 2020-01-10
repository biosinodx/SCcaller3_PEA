/*
 * 
 * Copyright (c) 2013, Zechen Chong <chongzechen@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "filter_kmer.h"
#define MAX_NM 1

int check_bam_alignment(bam1_t *b){
        uint32_t *cigars;
        uint32_t i, n_cigar;
        cigars = bam1_cigar(b);
        n_cigar = b->core.n_cigar;
		if ((b->core.flag & 4) != 0)
		    return 0;
        //for(i=0;i<n_cigar;i++){ if((cigars[i] & 0xFU) != 0) return 0; }
        for(i=0;i<n_cigar;i++){
            if((cigars[i] & 0xFU) == 3 || (cigars[i] & 0xFU) == 4)
                return 0;
        } // spliced or softclipped
        /*nm_aux = bam_aux_get(b, "NM");
        if(nm_aux && nm_aux[0] == 'i'){
                nm = bam_aux2i(nm_aux);
                if(nm <= MAX_NM) return 1;
        }*/
        return 1;
}


BitVec* fillin_bitvec(bamFile bamin, uint32_t ksize, BitVec *bt, u64hash *refhash, uint64_t *idx) {
	uint64_t KMER;
	BitVec* ret = bt;
	uint64_t k, r, kmask;
	uint32_t i, j, len, rid, lowq_num, flag;
	uint8_t *seq;
	char *quality;
	rid = 0;
	KMER = 0;
//	kmask = (1LLU << (2 * ksize)) - 1;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32-ksize)*2);
	seq = NULL;
	bam1_t *b;
	b = (bam1_t*)calloc(1, sizeof(bam1_t));
	while (bam_read1(bamin, b) >= 0) {
		flag = b->core.flag;
		rid ++;
//		if ((rid & 0xFFFFFU) == 0) {
//			fprintf(stdout, "[%s] parsed %10u reads\n", __FUNCTION__, rid);
//			fflush(stdout);
//		}
		k = 0;
		//fprintf(stdout, "rid = %u %s\n", rid, seq->seq.string); fflush(stdout);
		if (flag >= 256)
		    continue; // skip secondary alignments and other exceptions

		if(check_bam_alignment(b))
		    continue; //skip (nearly) identical reads to reference

		len = b->core.l_qseq;
		if (len < ksize+3)
		    continue;

		quality = (char*)bam1_qual(b);
		lowq_num = 0;
		for (j = 0; j < len; j++) {
			if (quality[j] + 33 <= '#') lowq_num ++;
		}
		if ((float)lowq_num/len >= 0.5)
		    continue; //skip reads with too many low quality bases

		//encap_bitvec(bt, 1024);
		seq = bam1_seq(b);
		for (i = 0; i < ksize-1; i++) {
			k = (k << 2) | bit4_bit_table[bam1_seqi(seq, i)];
		}
		for (i = 0; i <= len-ksize; i++) {
			if (quality[i+ksize-1] + 33 <= '#' && i+ksize < len && quality[i+ksize] + 33 <= '#') break; // trim ends with low qualities
			k = ((k << 2) | bit4_bit_table[bam1_seqi(seq, i+ksize-1)])  & kmask;
			//if (i + 1 < ksize) continue;
			r = dna_rev_seq(k, ksize);
			if (r < k) {
				KMER = r;
			} else {
				KMER = k;
			}
			if (exists_u64hash(refhash, KMER)) {
				one2bitvec(ret);
			} else {
				zero2bitvec(ret);
			}
			*idx = *idx+1;
		}
	}
	bam_destroy1(b);
	fprintf(stdout, "[%s] processed %10u reads\n", __FUNCTION__, rid);
	fflush(stdout);
	return ret;
}

kmerhash* build_readshash(bamFile bamin, uint32_t ksize, kmerhash *hash, BitVec *bt, uint64_t *idx, CBF *occ_table, uint32_t mincnt) {
	kmer_t KMER, *kmer;
	uint64_t k, r, kmask;
	uint32_t i, j, len, rid, lowq_num, occ, flag;
	int exists;
	char *quality;
	rid = 0;
	KMER.kmer = 0;
	KMER.cnt = 0;
	KMER.cnt2 = 0;
	occ = 0;
	uint8_t *seq;
	bam1_t *b;
	b = (bam1_t*)calloc(1, sizeof(bam1_t));
//	kmask = (1LLU << (2 * ksize)) - 1;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32-ksize)*2);
	seq = NULL;
	while (bam_read1(bamin, b) >= 0) {
		flag = b->core.flag;
		rid ++;
//		if ((rid & 0xFFFFFU) == 0) {
//			fprintf(stdout, "[%s] parsed %10u reads\n", __FUNCTION__, rid);
//			fflush(stdout);
//		}
		k = 0;
		if (flag >= 256) continue; // skip secondary alignments and other exceptions
		if(check_bam_alignment(b)) continue; //skip (nearly) identical reads to reference
		len = b->core.l_qseq;
		if (len < ksize+3) continue;
		quality = (char*)bam1_qual(b);
		lowq_num = 0;
		for (j = 0; j < len; j++) {
			if (quality[j] + 33 <= '#') lowq_num ++;
		}
		if ((float)lowq_num/len >= 0.5) continue; // too many low quality bases
		//encap_bitvec(bt, 1024);
		seq = bam1_seq(b);
		for (i = 0; i < ksize-1; i++) {
			k = (k << 2) | bit4_bit_table[bam1_seqi(seq, i)];
		}
		for (i = 0; i <= len-ksize; i++) {
			if (quality[i+ksize-1] + 33 <= '#' && i+ksize < len && quality[i+ksize] + 33 <= '#') break; // trim ends with low qualities
			k = ((k << 2) | bit4_bit_table[bam1_seqi(seq, i+ksize-1)])  & kmask;
			//if (i + 1 < ksize) continue;
			r = dna_rev_seq(k, ksize);
			if (r < k) {
				KMER.kmer = r;
			} else {
				KMER.kmer = k;
			}
			if (!get_bitvec(bt, *idx)) {
				occ = get_cbf(occ_table, &KMER.kmer, sizeof(uint64_t));
				if (occ + 1 >= mincnt) {
					kmer = prepare_kmerhash(hash, KMER, &exists);
					if (exists) {
							if (kmer->cnt < UNIQ_KMER_MAX_CNT)
									kmer->cnt ++;
					} else {
							kmer->kmer = KMER.kmer;
							kmer->cnt = occ+1;
							kmer->cnt2 = 0;
					}
				} else {
					put_cbf(occ_table, &KMER.kmer, sizeof(uint64_t));
				}
			}
			*idx = *idx+1;
		}
	}
	bam_destroy1(b);
	fprintf(stdout, "[%s] processed %10u reads\n", __FUNCTION__, rid);
	fflush(stdout);
	fprintf(stdout, "[%s]\n", date()); fflush(stdout);
	return hash;
}

kmerhash* build_refkmerhash(FileReader *fr, bamFile bamin, uint32_t ksize, kmerhash* hash, uint32_t mincnt, uint64_t off) {
	uint64_t KMER, *kmer;
	uint64_t k, r, kmask;
	uint32_t i, len, n_bit, rdlen = 101, rdnum = 0;
	int exists;
	BitVec *bt;
	CBF *occ_table;
	uint64_t idx = 0;
	u64hash *refhash = init_u64hash(1023);
	Sequence *seq;
	KMER = 0;
//	kmask = (1LLU << (2 * ksize)) - 1;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32-ksize)*2);
	seq = NULL;
//	fprintf(stdout, "Building reference kmer...\n");
//	fflush(stdout);
	while (fread_fasta(&seq, fr)) {
		fprintf(stdout, "[%s]Building reference [%s] kmer\n", __FUNCTION__, seq->name.string);
		fflush(stdout);
		k = 0;
		len = seq->seq.size;

		for (i = 0; i < ksize-1; i++) {
			k = (k << 2) | base_bit_table[(int)seq->seq.string[i]];
		}
		for (i = 0; i <= len-ksize; i++) {
			k = ((k << 2) | base_bit_table[(int)seq->seq.string[i+ksize-1]])  & kmask;
			r = dna_rev_seq(k, ksize);
			if (r < k) {
				KMER = r;
			} else {
				KMER = k;
			}
			kmer = prepare_u64hash(refhash, KMER, &exists);
			if (!exists) {
				*kmer = KMER;
			} 
		}
	}
	fprintf(stdout, "[%s]Calculating number of paired reads ...\n",__FUNCTION__);
	fflush(stdout);
	rdnum = 200000llu;
	bt = init_bitvec(2llu*rdnum*(rdlen-ksize+1));
	fillin_bitvec(bamin, ksize, bt, refhash, &idx);
	index_bitvec(bt);
	fprintf(stdout, "[%s]There are [%llu] reference kmers and [%llu] non-reference kmers\n", __FUNCTION__,
	        (unsigned long long)bt->n_ones,
	        (unsigned long long)bt->n_bit - bt->n_ones);
	fflush(stdout);
	fprintf(stdout, "[%s]Finished reference kmer building\n",__FUNCTION__);
	fflush(stdout);
	free_u64hash(refhash);
	fprintf(stdout, "[%s][%s]\n",__FUNCTION__, date());
	fflush(stdout);
	fprintf(stdout, "[%s]Freed reference hash and begin building reads hash table...\n",__FUNCTION__);
	fflush(stdout);
	fprintf(stdout, "[%s][%s]\n",__FUNCTION__, date());
	fflush(stdout);

	bam_seek(bamin, off, SEEK_SET);
	
	//TODO: process bitvec bt using CBF
	for(n_bit=2;n_bit<4 && (1U<<n_bit)<(mincnt+1);n_bit++);
	//occ_table = init_cbf(20 * 3000 * 1024 * 1024llu, n_bit, 3); // This is human genome, gsize=3000Mb, n_bit, nseed=3 TODO
	occ_table = init_cbf(60*1024*1024*1024llu, n_bit, 3); // This is human genome, gsize=3000Mb, n_bit, nseed=3 TODO
	//occ_table = init_cbf(60*1024*1024llu, n_bit, 3); // This is human genome, gsize=3000Mb, n_bit, nseed=3 TODO
	idx = 0;
	hash = build_readshash(bamin, ksize, hash, bt, &idx, occ_table, mincnt);
	/*
	index_bitvec(bt);
	fprintf(stdout, "%llu reference kmers out of total %llu kmers in reads\n", (unsigned long long)bt->n_ones, (unsigned long long)hash->size);
	fflush(stdout);
	*/
	free_cbf(occ_table);
	free_bitvec(bt);
	return hash;
}

uint32_t count_readnum(FileReader *readfr, int is_fq) {
	uint32_t ret = 0;
	Sequence *seq = NULL;
	while (is_fq?fread_fastq_adv(&seq, readfr, FASTQ_FLAG_NO_NAME):fread_fasta_adv(&seq, readfr, FASTA_FLAG_NO_NAME)) {
		ret ++;
	}
	return ret;
}

uint64_t filter_ref_kmers(kmerhash *hash, FileReader *fr, uint32_t ksize) {
	Sequence *seq;
	kmer_t KMER;
	uint64_t k, r, kmask, ret = 0;
	uint32_t i, len;
	seq = NULL;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32-ksize)*2);
	while (fread_fasta(&seq, fr)) {
		fprintf(stdout, "Filtering %s\r", seq->name.string);
		fflush(stdout);
		k = 0;
		len = seq->seq.size;
		if (len < ksize + 3) continue;
		for (i = 0; i < ksize-1; i++)
			k = (k << 2) | base_bit_table[(int)seq->seq.string[i]];
		for (i = 0; i <= len-ksize; i++) {
			k = ((k << 2) | base_bit_table[(int)seq->seq.string[i+ksize-1]])  & kmask;
			r = dna_rev_seq(k, ksize);
			if (r < k) {
				KMER.kmer = r;
			} else {
				KMER.kmer = k;
			}
			ret += remove_kmerhash(hash, KMER);
		}
	}
	return ret;
}

void  loadkmerseq(kmerhash *hash, uint32_t ksize, uint32_t mincnt, uint32_t maxcnt2, bamFile bamin, int is_somatic, chash *somanames, chash *germnames) {
	kmer_t KMER, *ret;
	uint32_t rid, len, j, lowq_num, flag; //, len_head;
	char *quality;
	uint64_t k, kmask = (1LLU << (2 * ksize)) - 1, r;
	uint8_t *seq;
	bam1_t *b;
	b = (bam1_t*)calloc(1, sizeof(bam1_t));

	rid = 0;
	KMER.cnt = 0;
	
	while (bam_read1(bamin, b) >= 0) {
		flag = b->core.flag;
		rid ++;
//		if ((rid & 0xFFFFFU) == 0) {
//			fprintf(stdout, "[%s] parsed %10u reads\n", __FUNCTION__, rid);
//			fflush(stdout);
//		}
		if (flag >= 256)
		    continue; // skip secondary alignments and other exceptions

		if(check_bam_alignment(b))
		    continue; //skip (nearly) identical reads to reference

		len = b->core.l_qseq;
		if (len < ksize+3)
		    continue;

		quality = (char*)bam1_qual(b);
		lowq_num = 0;
		for (j = 0; j < len; j++) {
			if (quality[j] + 33 <= '#') lowq_num ++;
		}
		if ((float)lowq_num/len >= 0.5)
		    continue; // too many low quality bases
		k = 0;
		seq = bam1_seq(b);
		//fprintf(stderr, "seq%d\t%s\n", i, seq);
		for (j = 0; j < len; j++) {
			k = ((k << 2) | bit4_bit_table[bam1_seqi(seq, j)]) & kmask;
			if (j + 1 < ksize) continue;
			if (quality[j] + 33 <= '#' && quality[j+1] + 33 <= '#') break; // trim ends with low qualities
			r = dna_rev_seq(k, ksize);
			if (r < k) {
				KMER.kmer = r;
			} else {
				KMER.kmer = k;
			}
			ret = get_kmerhash(hash, KMER);
			if (ret == NULL || ret->cnt < mincnt)  {
				continue;
			} 
			if (ret->cnt2 < maxcnt2 && (float)ret->cnt2/(ret->cnt+ret->cnt2) < 0.01) { //TODO magic number
				if (!exists_chash(somanames, bam1_qname(b))) 
					put_chash(somanames, strdup(bam1_qname(b)));
			}
			else if (ret->cnt2 > maxcnt2 * 2) { //TODO magic number
				if (is_somatic) continue;
				if (!exists_chash(germnames, bam1_qname(b)))
					put_chash(germnames, strdup(bam1_qname(b)));
			}
		}
	}
	bam_destroy1(b);
	return;
}
/*
static inline int cmp_pair_func(pair_t p1, pair_t p2, void *obj) {

	if (p1.pid == p2.pid) {
		if (p1.mt == p2.mt)
			return 0;
		else 
			return p1.mt-p2.mt;
	}


	return p1.pid-p2.pid;
	obj = obj;
}

define_quick_sort(sort_pairs, pair_t, cmp_pair_func);
*/

void dedup_pairs(samfile_t *somaout, samfile_t *germout, bamFile bamin, chash *somanames, chash *germnames, uint32_t ksize) {
	uint32_t rid = 0, len;
	//char pre1[256] = "", pre2[256] = "";
	//mut_type pre_t = SOMATIC;
	uint32_t flag;

	uint8_t *seq;
	bam1_t *b;
	b = (bam1_t*)calloc(1, sizeof(bam1_t));
//	kmask = (1LLU << (2 * ksize)) - 1;
	seq = NULL;
	//bam_seek(bamin, 0, SEEK_SET);
	while (bam_read1(bamin, b) >= 0) {
		flag = b->core.flag;
		rid ++;
//		if ((rid & 0xFFFFFU) == 0) {
//			fprintf(stdout, "[%s] parsed %10u reads\n", __FUNCTION__, rid);
//			fflush(stdout);
//		}
		if ((flag & 0x100) != 0) continue; // skip secondary alignments and other exceptions
		//if(check_bam_alignment(b)) continue; //skip (nearly) identical reads to reference
		len = b->core.l_qseq;
		if (len < ksize+3) continue;
		if (exists_chash(somanames, bam1_qname(b))) { 
			samwrite(somaout, b);
		}
		if (exists_chash(germnames, bam1_qname(b))) { 
			samwrite(germout, b);
		}
	}
	bam_destroy1(b);
}

void cal_ctrl_kmers(kmerhash *hash, samfile_t *ctrl, uint32_t ksize) {
	kmer_t KMER, *ret;
	uint64_t k, r, kmask;
	uint32_t i, len, rid, flag;
//	kmask = (1LLU << (2 * ksize)) - 1;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32-ksize)*2);
	ret = 0;
	rid = 0;
	uint8_t *seq;
	seq = NULL;
	bam1_t *b;
	b = (bam1_t*)calloc(1, sizeof(bam1_t));
	while (samread(ctrl, b) > 0) {
		flag = b->core.flag;
		rid ++;
//		if ((rid & 0xFFFFFU) == 0) {
//			fprintf(stdout, "[%s] parsed %10u reads\n", __FUNCTION__, rid);
//			fflush(stdout);
//		}
		k = 0;
		if (flag >= 256) continue; // skip secondary alignments and other exceptions
		len = b->core.l_qseq;
		seq = bam1_seq(b);
		if (len < ksize+3) continue;
		for (i = 0; i < ksize-1; i++)
			k = (k << 2) | bit4_bit_table[bam1_seqi(seq, i)];
		for (i = 0; i <= len-ksize; i++) {
			k = ((k << 2) | bit4_bit_table[bam1_seqi(seq, i+ksize-1)])  & kmask;
			//if (i + 1 < ksize) continue;
			r = dna_rev_seq(k, ksize);
			if (r < k) {
				KMER.kmer = r;
			} else {
				KMER.kmer = k;
			}
			ret = get_kmerhash(hash, KMER);
			if (ret) {
				if (ret->cnt2 < UNIQ_KMER_MAX_CNT)
					ret->cnt2 ++;
			}
		}
	}
	bam_destroy1(b);
	fprintf(stdout, "[%s] processed %10u reads\n", __FUNCTION__, rid);
	fflush(stdout);
	//return ret;
}

inline int cmp_kmer(const void *e1, const void *e2) {
	kmer_t *k1, *k2;
	k1 = (kmer_t*)e1;
	k2 = (kmer_t*)e2;

	if (k1->cnt > k2->cnt) return -1;
	if (k1->cnt < k2->cnt) return 1;
	return 0;
}


int usage() {
	printf("novoBreak - a tool for discovering somatic sv breakpoints\n"
		   "Auther: Zechen Chong <zchong@mdanderson.org> \n"
		   "Version: 1.1 (r20151007)\n"
		   "Usage:\n"
		   "  novoBreak -i <tumorbam> -c <normalbam> -r <reference> -o <output.kmer> [options]\n"
		   "Options:\n"
		   "  -h             This help\n"
		   "  -i <string>    Tumor bam file\n"
		   "  -c <string>    Normal bam file\n"
		   "  -r <string>    Reference file in fasta format\n"
           "  -k <int>       Kmer size, <=31 [31]\n"
		   "  -o <string>    Output kmer\n"
		   "  -g <int>       Output germline events [0]\n"  
		   "  -m <int>       Minimum kmer count regarded as novo kmers [3]\n"
		   );

	return 1;
}

int main(int argc, char **argv) {
	kmerhash *khash;
	kmer_t KMER;
	bamFile inf;
	samfile_t *ctrlf;
	FileReader *reff;
	chash *somanames, *germnames;
	FILE *out;
	samfile_t *somaout, *germout;
	char *infile; // Tumor bam file
    char *outfile; // Output kmer
    char *ctrlfile; // Normal bam file
    char *reffile; // Reference file in fasta format
	int c;
	int is_somatic = 1; // 输入参数 1 somatic  0 germline
	uint32_t ksize = 31, i, maxcnt2 = 3;
	uint32_t mincnt = 3; // 输入参数
	uint64_t ret;
	infile = outfile = ctrlfile = reffile = NULL;
	bam_header_t *header = NULL;

	while ((c = getopt(argc, argv, "hi:c:k:o:m:r:g")) != -1) {
		switch (c) {
			case 'h':
			    return usage();
			case 'i':
			    infile = optarg;
			    break;
			case 'c':
			    ctrlfile = optarg;
			    break;
			case 'k':
			    ksize = atoi(optarg);
			    break;
			case 'o':
			    outfile = optarg;
			    break;
			case 'm':
			    mincnt = atoi(optarg);
			    break;
			case 'r':
			    reffile = optarg;
			    break;
			case 'g':
			    is_somatic = 0;
			    break;
			default:
			    return usage();
		}
	}
	if((!access(outfile,0)) && (!access("somaticreads.bam",0)) && (!access("germlinereads.bam",0)))
    {
	    fprintf(stdout,"%s, somaticreads.bam, germlinereads.bam already exist. Quit.\n", outfile);
        return 0;
    }



	if(reffile == NULL || infile == NULL || outfile == NULL)
	    return usage();
	if (ksize > 31 || ksize < 11)
	    return usage();
	if ((reff = fopen_filereader(reffile)) == NULL) {
		fprintf(stderr, " -- Cannot open reference file in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	} else{
	    fprintf(stderr, "[%s]reffile=[%s]\n",__FUNCTION__, reffile);
	    fflush(stderr);
	}

	if ((inf = bam_open(infile, "r")) == NULL) {
		fprintf(stderr, " -- Cannot open input file in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	} else{
        fprintf(stderr, "[%s]infile=[%s]\n",__FUNCTION__,infile);
        fflush(stderr);
	}

	if (ctrlfile) {
		if ((ctrlf = samopen(ctrlfile, "rb", NULL)) == NULL) {
			fprintf(stderr, " -- Cannot open input file in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
			abort();
		} else {
            fprintf(stderr, "[%s]ctrlfile=[%s]\n",__FUNCTION__,ctrlfile);
            fflush(stderr);
        }
	} else ctrlf = NULL;
	if ((out = fopen(outfile, "w")) == NULL) {
		fprintf(stderr, " -- Please provide output file -o [output] --\n");
		fflush(stderr);
		abort();
	}else{
	    fprintf(stderr, "[%s]outfile=[%s]\n",__FUNCTION__,outfile);
	    fflush(stderr);
	}
	header = bam_header_read(inf);

	if ((somaout = samopen("somaticreads.bam", "wb", header)) == NULL) {
		fprintf(stderr, " cannot open file to write\n");
		fflush(stderr);
		abort();
	}

	if ((germout = samopen("germlinereads.bam", "wb", header)) == NULL) {
		fprintf(stderr, " cannot open file to write\n");
		fflush(stderr);
		abort();
	}
	fprintf(stdout, "[%s][%s]\n", __FUNCTION__, date());
	fflush(stdout);
	khash = init_kmerhash(1023);
	//refhash = init_kmerhash(1023);
	uint64_t off = bgzf_tell(inf);
	bam_header_destroy(header);
	khash = build_refkmerhash(reff, inf, ksize, khash, mincnt, off);
	fprintf(stdout, "[%s] [%llu] kmers loaded\n\n", __FUNCTION__,(unsigned long long)count_kmerhash(khash));
	fflush(stdout);
	
	fprintf(stdout, "[%s][%s]\n", __FUNCTION__, date());
	fflush(stdout);
	fprintf(stdout, "[%s] Calculating kmers from control...\n",__FUNCTION__);
	fflush(stdout);
	cal_ctrl_kmers(khash, ctrlf, ksize );

	ret = 0;
	fprintf(stdout, "[%s][%s]\n", __FUNCTION__, date());
	fflush(stdout);
	fprintf(stdout, "[%s] Loading sequences...\n",__FUNCTION__);
	fflush(stdout);
	somanames = init_chash(1023);
	germnames = init_chash(1023);
	bam_seek(inf, off, SEEK_SET);
	loadkmerseq(khash, ksize, mincnt, maxcnt2, inf, is_somatic, somanames, germnames); //TODO
	//fprintf(stdout, "There are %llu pairs loaded\n\n", (unsigned long long)count_pairv(pairs));
	fprintf(stdout, "[%s][%s]\n", __FUNCTION__, date());
	fflush(stdout);
	fprintf(stdout, "[%s]Remove duplicate sequences...\n",__FUNCTION__);
	fflush(stdout);
	bam_seek(inf, off, SEEK_SET);
	dedup_pairs(somaout, germout, inf, somanames, germnames, ksize);
	fprintf(stdout, "[%s][%s]\n\n", __FUNCTION__, date());
	fflush(stdout);
	reset_iter_kmerhash(khash);
	while (iter_kmerhash(khash, &KMER)) {
		//if (KMER.cnt < mincnt || KMER.cnt2 > maxcnt2) continue;
		if (KMER.cnt < mincnt) continue;
		if (KMER.cnt2 < maxcnt2 && (float)KMER.cnt2/(KMER.cnt+KMER.cnt2) < 0.01) {
			ret ++;
			for (i = 0; i < ksize; i++) {
				fprintf(out, "%c", bit_base_table[(KMER.kmer >> ((ksize-1-i) << 1)) & 0x03]);
			}
			fprintf(out, "\t%llu\t%llu\t%s\n", (unsigned long long)KMER.cnt, (unsigned long long)KMER.cnt2, "SOMATIC");
		} else if (KMER.cnt2 > maxcnt2*2) {
			ret ++;
			for (i = 0; i < ksize; i++) {
				fprintf(out, "%c", bit_base_table[(KMER.kmer >> ((ksize-1-i) << 1)) & 0x03]);
			}
			fprintf(out, "\t%llu\t%llu\t%s\n", (unsigned long long)KMER.cnt, (unsigned long long)KMER.cnt2, "GERMLINE");
		}
	}
	
	fprintf(stdout, "[%s][%llu] kmers passed the minimum frequency cutoff in tumor (%u) and maximum frequency cutoff in normal (%u)\n", __FUNCTION__, (unsigned long long)ret, mincnt, maxcnt2);
	fflush(stdout);
	free_kmerhash(khash);
	bam_close(inf);
	samclose(ctrlf);
	fclose_filereader(reff);
	fclose(out);
	samclose(somaout);
	samclose(germout);
	for (i = 0; i < somanames->size; i++) {
		if(exists_entity(somanames->flags, i)) { //TODO: use a list to hold the names instead of accessing hashset source code
			free(somanames->array[i]);
		}
	}
	for (i = 0; i < germnames->size; i++) {
		if(exists_entity(germnames->flags, i)) {
			free(germnames->array[i]);
		}
	}
	free_chash(somanames);
	free_chash(germnames);
	fprintf(stderr, "[%s]Program exit normally\n",__FUNCTION__);
	return 0;
}

