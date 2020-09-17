# python 3
# require, bwa, samtools
import os
import re
import argparse
import sys
import string

#Parse commandline arguments
parser=argparse.ArgumentParser(description="sc3_rewrite.py\nBased on manta vcf output")
parser.add_argument("--cellid",type=str, required=True, help="Cell ID. Required.")
parser.add_argument("--vcf",type=str, required=True, help="Input vcf file. Required.")
parser.add_argument("--wkdir",type=str, default="./", required=False, help="working dir. default: ./")
parser.add_argument("--cbam",type=str, required=True, help="Full pwd to cell bam. Required.")
parser.add_argument("-g","--genome",type=str, required=True, help="Full pwd to the reference genome. Required.")
parser.add_argument("-s","--start",type=int, help="First SV candidate to process. Required", required=True)
parser.add_argument("-e","--end", type=int, help="Last SV candidate to process. Required", required=True)
parser.add_argument("--flank", type=int, default=100, help="No. base pairs of flanking regions of breakpoints; default: 100.")
parser.add_argument("--report",type=bool, default=False, required=False, help="Reporting; default: False")
parser.add_argument("--limit", type=int, help="Max. no. read for considering a SV. default:1000", default=1000)
args=parser.parse_args()

cell=args.cellid
delly_vcf=open(args.vcf, 'r')
ref_fasta=open(args.genome, 'r')
cell_bam=args.cbam
wkdir=args.wkdir

# midresult=False; # report middle result yes or no
midresult=args.report
flankingregion=args.flank
processing_start=args.start
processing_end=args.end

readlimit=args.limit

cmd="mkdir -p "+ wkdir +"/sv_temp"+"_"+str(processing_start)+"_"+str(processing_end); os.system(cmd)
os.chdir(wkdir+"/sv_temp"+"_"+str(processing_start)+"_"+str(processing_end))

##### I. processing the reference genome #####
genome = ref_fasta.read()
ref_fasta.close()
chromo = genome.split('>')
chromo_index = list()
i = 1
while i < len(chromo):
	tmp = chromo[i][0:200]
	chromo_index.append((tmp.split('\n')[0]).split(' ')[0])
	i = i + 1

basepair = list()
i = 1
while i < len(chromo):
	tmp = chromo[i]
	tmp = tmp.split('\n')
	tmp[0] = '0'
	tmp1 = str()
	for element in tmp:
		tmp1 = tmp1 + element
	basepair.append(tmp1)
	i = i + 1

### function returns the tri-nucleotide
def triple(a_chr, b_pos, c_start=-500, d_end=500):
	if int(b_pos) + int(c_start) < 1:
		return "pos.exceedchrlength"
	a_chr = str(a_chr) # a_chr, string
	#if 'e+' in b_pos:
	#       b_pos = float(b_pos.split('e+')[0])*10**float(b_pos.split('e+')[1])
	b_pos = int(b_pos) # b_pos, int
	i=0
	pick = -1
	while i < len(chromo_index):
		if a_chr == chromo_index[i]:
			pick = i
			break
		i = i + 1
	if pick == -1:
		return "chr.notmatching"
	if b_pos + d_end + 1 >= len(basepair[pick]):
		return "pos.exceedchrlength"
	#result_triple = [basepair[i][b_pos-1], basepair[i][b_pos], basepair[i][b_pos+1]]
	#result_triple = []
	result_triple=''
	j = c_start
	while j <= d_end:
		# result_triple.append(basepair[pick][b_pos+j]) # it would return a list, i.g. ['A','T','C'], corresponding to -1, 0, +1 position of your mutation
		result_triple = result_triple + basepair[pick][b_pos+j]
		j = j + 1
	return result_triple

### example of using it, in python
#triple('chrX',123)
# note, 'chrX' should be provided in a proper way, you can check the variable 'chromo_index' of how to write it.


##### II. function(s) to decode manta vcf file #####
def format_mantaline(fd_dellyline):
	# input a line in vcf of delly
	# output a list of strings: SVid, refA_chr, refA_pos, refB_chr, refB_pos, refSV_type, refSV_CT, refSV_consensus, refA_sequence, refB_sequence, mateid
	# skip those not passing filter and/or without concensus sequence; -> return None
	fd_dellyline=fd_dellyline.split('\n')[0]
	fd_dellyline=fd_dellyline.split('\t')
	if fd_dellyline[6]!="PASS":
		return(None)
	fd_refA_chr=fd_dellyline[0]
	fd_refA_pos=fd_dellyline[1]
	fd_tmp=fd_dellyline[7].split(';')
	fd_svid=fd_dellyline[2]
	fd_refB_chr=''
	fd_refB_pos=''
	fd_refSV_consensus=''
	fd_refSV_ct=''
	fd_refA_sequence=''
	fd_refB_sequence=''
	fd_mateid=''
	# fd_refSV_type=fd_dellyline[4]
	# fd_refSV_type=(fd_refSV_type.split('<')[1]).split('>')[0]
	for fd_i in fd_tmp:
		if "CHR2=" in fd_i:
			fd_refB_chr=fd_i.split("=")[1]
		if "END2=" == fd_i[0:5]:
			fd_refB_pos=fd_i.split("=")[1]
		if "REFSV=" in fd_i:
			fd_refSV_consensus=fd_i.split("=")[1]
		if "CT=" in fd_i:
			fd_refSV_ct=fd_i.split("=")[1]
		if "SIMPLE_TYPE=" in fd_i:
			fd_refSV_type=fd_i.split("=")[1]
		if "REFA=" in fd_i:
			fd_refA_sequence=fd_i.split("=")[1]
		if "REFB=" in fd_i:
			fd_refB_sequence=fd_i.split("=")[1]
		if "MATEID=" in fd_i:
			fd_mateid=fd_i.split("=")[1]
	if fd_refSV_consensus=="":
		return(None)
	return([fd_svid, fd_refA_chr, fd_refA_pos, fd_refB_chr, fd_refB_pos, fd_refSV_type, fd_refSV_ct, fd_refSV_consensus, fd_refA_sequence, fd_refB_sequence, fd_mateid])


###
def sam2fq(samfile, fqfile):
	s2f_filein = open(samfile)
	s2f_fileout = open(fqfile, 'w')
	while True:
		tmp=s2f_filein.readline()
		if len(tmp)==0:
			break
		tmp=tmp.split("\t")
		if len(tmp)<11:
			continue
		s2f_fileout.write("@"+tmp[0]+'\n'+tmp[9]+'\n+\n'+tmp[10]+'\n')
	s2f_filein.close()
	s2f_fileout.close()


##### search support SV type #####
def readsam(rs_sam):
	sp_sam0file=open(rs_sam)
	sp_sam0lines=sp_sam0file.read()
	sp_sam0file.close()
	sp_sam0lines=sp_sam0lines.split('\n')
	sp_sam0lines_clean=[]
	for item in sp_sam0lines:
		if len(item)!=0:
			if item[0]!="@" and item[1]!="4": # filter out header, and unaligned reads
				sp_sam0lines_clean.append(item.split('\t'))
	return(sp_sam0lines_clean)

def supportref(read, sp_sam0in, sp_bp):
	sp_perfectalignment=True
	sp_passbp=False
	for item in sp_sam0in:
		if item[0]==read:
			# re module
			# https://www.cnblogs.com/huxi/archive/2010/07/04/1771073.html
			# perfect alignment: re.match(r'^\d+M$', '150M')!=None
			# NOT perfect alignment: re.match(r'^\d+M$', '150M')==None
			if re.match(r'^\d+M$', item[5])==None: # none perfect alignment
				sp_perfectalignment=False
			elif float(sp_bp)>float(item[3])+10 and float(sp_bp)<float(item[3])+float(item[5][0:(len(item[5])-1)])-10: # read with perfectalignment passes bp +- 10bp
				sp_passbp=True
	# along as sp_passbp==True: meaning reads support refX not refSV
	# if sp_perfectalignment==True, and sp_passbp==False: meaning we do not know if the read support refX or refSV
	# if sp_perfectalignment==False, and sp_passbp==False: meaning we should go to enhanced alignment, to check
	if sp_passbp==True:
		return(["refX", 1.0])
	elif sp_perfectalignment==True and sp_passbp==False:
		return(["exclude", 0])
	else:
		return(["checkfurther", 0])

# supportref("A00197:231:H57GVDSXY:3:1610:8657:18364")
# supportref("A00197:231:H57GVDSXY:3:1243:24162:24893")
# supportref("A00197:231:H57GVDSXY:3:2619:30074:5885")
# supportref("A00197:231:H57GVDSXY:3:2371:26603:15342")

def supportenhance(read, sp_samXin, sp_samSVin, sp_minoverlap=0.7, sp_mindiff=0.1):
	#read="A00197:231:H57GVDSXY:3:2619:30074:5885"
	sp_refX_overlap=[0.0]
	sp_refX_mlen=[0.0]
	for item in sp_samXin:
		if item[0]==read and item[1]!='4':
			m=re.search(r'^\d+M$', item[5])
			if m!=None: # perfect alignment
				sp_refX_mlen.append(float(item[5][m.start():(m.end()-1)]))
				sp_refX_overlap.append(1.0)
			else: # none perfect alignment
				m=re.search(r'\d+M', item[5])
				m_len=float(item[5][m.start():(m.end()-1)]); # M lenght
				# M left length
				m_l=re.search(r'^\d+\S', item[5])
				if m_l==None:
					m_r_len=0.0
				elif item[5][m_l.end()-1]!="M":
					m_l_len=float(item[5][m_l.start():(m_l.end()-1)]);
				else:
					m_l_len=0.0
				# M right length
				m_r=re.search(r'\d+\S$', item[5])
				if m_r==None:
					m_r_len=0.0
				elif item[5][m_r.end()-1]!="M":
					m_r_len=float(item[5][m_r.start():(m_r.end()-1)]);
				else:
					m_r_len=0.0
				# three intervals: [start-m_l_len, start]?, [start, start+m_len]M, [start+m_len, start+m_len+m_r_len]?; start: item[3]
				# refX or refSV intervals: [1, len(refSV)]
				# find intersect of the above two as dinominator
				sp_refX_mlen.append(m_len)
				sp_refX_overlap.append(m_len / (min([len(refSV), float(item[3])+m_len+m_r_len]) - max([1, float(item[3])-m_l_len])))
	for i in range(0, len(sp_refX_overlap)):
		if sp_refX_mlen[i]==max(sp_refX_mlen):
			sp_refX_overlap=sp_refX_overlap[i]
			break
	#
	sp_refSV_overlap=[0.0]
	sp_refSV_mlen=[0.0]
	for item in sp_samSVin:
		if item[0]==read and item[1]!='4':
			m=re.search(r'^\d+M$', item[5])
			if m!=None: # perfect alignment
				sp_refSV_mlen.append(float(item[5][m.start():(m.end()-1)]))
				sp_refSV_overlap.append(1.0)
			else: # none perfect alignment
				m=re.search(r'\d+M', item[5])
				m_len=float(item[5][m.start():(m.end()-1)]); # M lenght
				# M left length
				m_l=re.search(r'^\d+\S', item[5])
				if m_l==None:
					m_r_len=0.0
				elif item[5][m_l.end()-1]!="M":
					m_l_len=float(item[5][m_l.start():(m_l.end()-1)]);
				else:
					m_l_len=0.0
				# M right length
				m_r=re.search(r'\d+\S$', item[5])
				if m_r==None:
					m_r_len=0.0
				elif item[5][m_r.end()-1]!="M":
					m_r_len=float(item[5][m_r.start():(m_r.end()-1)]);
				else:
					m_r_len=0.0
				# three intervals: [start-m_l_len, start]?, [start, start+m_len]M, [start+m_len, start+m_len+m_r_len]?; start: item[3]
				# refX or refSV intervals: [1, len(refSV)]
				# find intersect of the above two as dinominator
				sp_refSV_mlen.append(m_len)
				sp_refSV_overlap.append(m_len / (min([len(refSV), float(item[3])+m_len+m_r_len]) - max([1, float(item[3])-m_l_len])))
	for i in range(0, len(sp_refSV_overlap)):
		if sp_refSV_mlen[i]==max(sp_refSV_mlen):
			sp_refSV_overlap=sp_refSV_overlap[i]
			break
	#sp_refSV_overlap=max(sp_refSV_overlap);
	#
	if ((sp_refSV_overlap - sp_refX_overlap) >= sp_mindiff) and sp_refSV_overlap > sp_minoverlap:
		if sp_refSV_overlap > 1:
			sp_refSV_overlap=1.0
		return(["refSV", sp_refSV_overlap])
	elif ((sp_refX_overlap - sp_refSV_overlap) >= sp_mindiff) and sp_refX_overlap > sp_minoverlap:
		if sp_refX_overlap > 1:
			sp_refX_overlap=1.0
		return(["refX", sp_refX_overlap])
	else:
		return(["unknown", 0])

def support(sp_fastq, sp_sam0, sp_samX, sp_samSV, sp_ref, sp_minoverlap=0.7, sp_mindiff=0.1):
	#sp_fastq="DEL00001068_refA.fq"
	#sp_sam0="DEL00001068_refA.0.sam"
	#sp_samX="DEL00001068_refA_refA.sam"
	#sp_samSV="DEL00001068_refA_refSV.sam"
	#sp_ref="A"
	#sp_minoverlap=0.7
	#sp_mindiff=0.1
	#sp_fastq="DEL00001068_refB.fq"
	#sp_sam0="DEL00001068_refB.0.sam"
	#sp_samX="DEL00001068_refB_refB.sam"
	#sp_samSV="DEL00001068_refB_refSV.sam"
	#sp_ref="B"
	if sp_ref=="A":
		sp_bp=x[2]
	if sp_ref=="B":
		sp_bp=x[4]
	# collect read names
	sp_reads=[]
	sp_fastqfile=open(sp_fastq)
	sp_fastq_index=0
	while True:
		sp_fastqline=sp_fastqfile.readline()
		if len(sp_fastqline)==0:
			break
			#return()###
		if sp_fastqline[0]=="@" and sp_fastq_index % 4 ==0:
			sp_reads.append(sp_fastqline[1:].split("\n")[0])
		sp_fastq_index=sp_fastq_index+1
	sp_fastqfile.close()
	#
	sp_reads_uniq=[]
	for item in sp_reads:
		if not item in sp_reads_uniq:
			sp_reads_uniq.append(item)
	# load the original sam file
	sp_sam0in = readsam(sp_sam0)
	# load the refA/B sam file
	sp_samXin = readsam(sp_samX)
	# load the refSV sam file
	sp_samSVin = readsam(sp_samSV)
	#
	sp_readname=[]
	sp_readclass=[]
	sp_readalignq=[]
	for read in sp_reads_uniq:
		# check status from original bam alignment
		read_status=supportref(read, sp_sam0in=sp_sam0in, sp_bp=sp_bp)
		if read_status[0]=="refX":
			sp_readname.append(read)
			sp_readclass.append(read_status[0])
			sp_readalignq.append(read_status[1])
		elif read_status=="exclude":
			sp_readname.append(read)
			sp_readclass.append(read_status[0])
			sp_readalignq.append(read_status[1])
		else: # for those "checkfurther"
			# check status from refA/B(X) vs refSV enhanced alignment
			sp_readname.append(read)
			read_status=supportenhance(read, sp_samXin=sp_samXin, sp_samSVin=sp_samSVin, sp_minoverlap=0.7, sp_mindiff=0.1)
			sp_readclass.append(read_status[0])
			sp_readalignq.append(read_status[1])
	sp_out=[]
	for i in range(0, len(sp_readname)):
		sp_out.append([sp_readname[i], sp_readclass[i], sp_readalignq[i]])
	return(sp_out)

# support(sp_fastq="DEL00001068_refA.fq", sp_sam0="DEL00001068_refA.0.sam", sp_samX="DEL00001068_refA_refA.sam", sp_samSV="DEL00001068_refA_refSV.sam", sp_ref="A")
# support(sp_fastq="DEL00001068_refB.fq", sp_sam0="DEL00001068_refB.0.sam", sp_samX="DEL00001068_refB_refB.sam", sp_samSV="DEL00001068_refB_refSV.sam", sp_ref="B")


##### MAIN #####
### get a pool of bed files and grep reads using python scripts together for 200 SVs a batch ###
tmpout=open(cell+'_regions.bed','w')
svlist_id=[]
svlist_chrA=[]
svlist_posA=[]
svlist_chrB=[]
svlist_posB=[]
xlist=[]

delly_line_index=-1
while True:
	delly_line=delly_vcf.readline()
	if len(delly_line)==0:
		break
	if delly_line[0]=="#":
		continue
	delly_line_index=delly_line_index+1
	if delly_line_index < processing_start:
		continue
	if delly_line_index > processing_end:
		break
	x=format_mantaline(delly_line)
	if x==None:
		continue
	# to save time; with samtools greping all the reads, can be done first by all SVs in one bed file; then grep again with only one SV in one bed file.
	tmp=x[1] + '\t' + str(int(float(x[2]))-flankingregion) + '\t' + str(int(float(x[2]))+flankingregion) + '\t' + x[0] + '\n'
	tmpout.write(tmp)
	tmp=x[3] + '\t' + str(int(float(x[4]))-flankingregion) + '\t' + str(int(float(x[4]))+flankingregion) + '\t' + x[0] + '\n'
	tmpout.write(tmp)
	xlist.append(x)
	svlist_id.append(x[0])
	svlist_chrA.append(x[1]); svlist_posA.append(int(float(x[2])))
	svlist_chrB.append(x[3]); svlist_posB.append(int(float(x[4])))
	#if "DEL00001068" in delly_line:
	#       break

tmpout.close()
# my extraction getting a little more reads than only using samtools view; but I checked, it is correct
# https://blog.csdn.net/windone0109/article/details/8895875

cmd = "samtools view -@ 4 -F 256 -q 1 -L " + cell + '_regions.bed ' + cell_bam
osopen=os.popen(cmd)
while True:
	osline=osopen.readline()
	if len(osline)==0:
		break
	tmp=osline.split('\t')
	tmp_chr=tmp[2]; # chr
	tmp_start=int(tmp[3]) # left start
	# tmp[5]; # cigar
	tmpm=re.search(r'\d+M', tmp[5])
	if tmpm==None:
		continue	
	tmp_end=int(tmp[3]) + int(tmp[5][tmpm.start():(tmpm.end()-1)]); # right end
	for i in range(0, len(svlist_chrA)):
		if svlist_chrA[i]==tmp_chr and tmp_start-flankingregion<=svlist_posA[i] and tmp_end+flankingregion>=svlist_posA[i]:
			tmpfile=open(svlist_id[i]+'_refA.0.sam','a')
			tmpfile.write(osline)
			tmpfile.close()
		if svlist_chrB[i]==tmp_chr and tmp_start-flankingregion<=svlist_posB[i] and tmp_end+flankingregion>=svlist_posB[i]:
			tmpfile=open(svlist_id[i]+'_refB.0.sam','a')
			tmpfile.write(osline)
			tmpfile.close()

osopen.close()

###
for sv_index in range(0, len(xlist)):
	x=xlist[sv_index]
	refA = x[8]; # triple(x[1],x[2], c_start=-int(round(len(x[7])/2)), d_end=int(round(len(x[7])/2)))
	refB = x[9]; # triple(x[3],x[4], c_start=-int(round(len(x[7])/2)), d_end=int(round(len(x[7])/2)))
	refSV = x[7]; # just take the consensus reported by the software tools
	tmpout=open(x[0]+'_refA.fa','w'); tmpout.write(">"+x[0]+"|refA\n"+refA+"\n"); tmpout.close()
	tmpout=open(x[0]+'_refB.fa','w'); tmpout.write(">"+x[0]+"|refB\n"+refB+"\n"); tmpout.close()
	tmpout=open(x[0]+'_refSV.fa','w'); tmpout.write(">"+x[0]+"|refSV\n"+refSV+"\n"); tmpout.close()
	# use bwa to do alignment
	# bwa mem ref.fa read1.fq read2.fq
	cmd="bwa index "+x[0]+"_refA.fa"; os.system(cmd)
	cmd="bwa index "+x[0]+"_refB.fa"; os.system(cmd)
	cmd="bwa index "+x[0]+"_refSV.fa"; os.system(cmd)
	# https://seqome.com/convert-bam-file-fastq/
	# convert sam to fastq file
	# cat ${a}.0.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > ${a}.fq
	if os.path.exists(x[0] + "_refA.0.sam"):
		num_lines = sum(1 for line in open(x[0] + "_refA.0.sam"))
		if num_lines > readlimit:
			continue
		sam2fq(x[0] + "_refA.0.sam", x[0] + "_refA.fq")
		# bwa the refA region to refA.fa and refSV.fa
		cmd="bwa mem -M "+x[0]+"_refA.fa "+x[0]+"_refA.fq > "+x[0]+"_refA_refA.sam"; os.system(cmd)
		cmd="bwa mem -M "+x[0]+"_refSV.fa "+x[0]+"_refA.fq > "+x[0]+"_refA_refSV.sam"; os.system(cmd)
		# for every read in refA region, check if supporting refA or refSV
		# exclude: for unuseful reads; often do not cover breakpoints
		# unknown: cannot ditermine from refA or refSV
		# refX: refA or refB
		# refSV: refSV
		suppA = support(sp_fastq=x[0]+"_refA.fq", sp_sam0=x[0]+"_refA.0.sam", sp_samX=x[0]+"_refA_refA.sam", sp_samSV=x[0]+"_refA_refSV.sam", sp_ref="A")
		tmpopen=open("pool.support.txt","a")
		for ele in suppA:
			tmpopen.write(x[0] + "\trefA\t" + ele[0]+"\t"+ele[1]+"\t" + str(ele[2])+ "\n")
		tmpopen.close()
	if os.path.exists(x[0] + "_refB.0.sam"):
		num_lines = sum(1 for line in open(x[0] + "_refB.0.sam"))
		if num_lines > readlimit:
			continue
		sam2fq(x[0] + "_refB.0.sam", x[0] + "_refB.fq")
		# bwa the refB region to refB.fa and refSV.fa
		cmd="bwa mem -M "+x[0]+"_refB.fa "+x[0]+"_refB.fq > "+x[0]+"_refB_refB.sam"; os.system(cmd)
		cmd="bwa mem -M "+x[0]+"_refSV.fa "+x[0]+"_refB.fq > "+x[0]+"_refB_refSV.sam"; os.system(cmd)
		suppB = support(sp_fastq=x[0]+"_refB.fq", sp_sam0=x[0]+"_refB.0.sam", sp_samX=x[0]+"_refB_refB.sam", sp_samSV=x[0]+"_refB_refSV.sam", sp_ref="B")
		tmpopen=open("pool.support.txt","a")
		for ele in suppB:
			tmpopen.write(x[0] + "\trefB\t" + ele[0]+"\t"+ele[1]+"\t" + str(ele[2])+ "\n")
		tmpopen.close()
	cmd="rm "+ x[0]+"*.fa.*";os.system(cmd)
	if midresult==False:
		cmd="rm "+ x[0]+"_*"; os.system(cmd)


