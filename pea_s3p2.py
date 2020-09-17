# python 3
# require, bwa, samtools
import os
import re
import argparse
import sys
import string
from scipy.stats import binom

#Parse commandline arguments
parser=argparse.ArgumentParser(description="sc3_rewrite.py\nBased on delly vcf output")
parser.add_argument("--cellid",type=str, required=True, help="Cell ID. Required.")
parser.add_argument("--vcf",type=str, required=True, help="Input vcf file. Required.")
parser.add_argument("--wkdir",type=str, default="./", required=False, help="working dir. default: ./")
parser.add_argument("--cbam",type=str, required=True, help="Full pwd to cell bam. Required.")
parser.add_argument("-g","--genome",type=str, required=True, help="Full pwd to the reference genome. Required.")
parser.add_argument("-s","--start",type=int, help="First SV candidate to process. Required", required=True)
parser.add_argument("-e","--end", type=int, help="Last SV candidate to process. Required", required=True)
parser.add_argument("--flank", type=int, default=100, help="No. base pairs of flanking regions of breakpoints; default: 100.")
parser.add_argument("--report",type=bool, default=False, required=False, help="Reporting; default: False")

# additional to part1
parser.add_argument("--hsnp", type=str, help="Bulk hsnp vcf file. Required", required=True)
parser.add_argument("--minmq", type=int, help="Min. mapQ. to consider hSNP in cell bam. default:30", default=30)
parser.add_argument("--lam", type=int, help="Lambda for sccaller. default:10000", default=10000)
parser.add_argument("--theta", type=float, help="Default theta for sccaller. default:0.75", default=0.75)
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

# bulkhsnp_vcf=open('/data/x001/xdong/2019-sccaller3/vcf/G12878_B.hg38.hsnp.biallelic.dbsnp.vcf', 'r')
bulkhsnp_vcf=open(args.hsnp, 'r')

# ref='/data/x001/xdong/reference/homo_sapiens/GRCh38.broad/references-hg38-v0-Homo_sapiens_assembly38.fasta'
ref=args.genome
# minmq_hsnp=30; # min mapq for hsnp; default=30
minmq_hsnp=args.minmq
# lam=10000 # default lambda=10000
lam=args.lam
# theta=0.75; # default theta=0.75
theta=args.theta
# readlimit=1000; # maximum limit on no. reads, default=1000
readlimit=args.limit

os.chdir(wkdir+"/sv_temp"+"_"+str(processing_start)+"_"+str(processing_end))

##### II. function(s) to decode delly vcf file #####
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

def format_dellyline(fd_dellyline):
	# input a line in vcf of delly
	# output a list of strings: SVid, refA_chr, refA_pos, refB_chr, refB_pos, refSV_type, refSV_CT, refSV_consensus
	fd_dellyline=fd_dellyline.split('\n')[0]
	fd_dellyline=fd_dellyline.split('\t')
	fd_refA_chr=fd_dellyline[0]
	fd_refA_pos=fd_dellyline[1]
	fd_tmp=fd_dellyline[7].split(';')
	fd_svid=fd_dellyline[2]
	fd_refB_chr=''
	fd_refB_pos=''
	fd_refSV_consensus=''
	fd_refSV_ct=''
	# fd_refSV_type=fd_dellyline[4]
	# fd_refSV_type=(fd_refSV_type.split('<')[1]).split('>')[0]
	for fd_i in fd_tmp:
		if "CHR2=" in fd_i:
			fd_refB_chr=fd_i.split("=")[1]
		if "END=" == fd_i[0:4]:
			fd_refB_pos=fd_i.split("=")[1]
		if "CONSENSUS=" in fd_i:
			fd_refSV_consensus=fd_i.split("=")[1]
		if "CT=" in fd_i:
			fd_refSV_ct=fd_i.split("=")[1]
		if "SVTYPE=" in fd_i:
			fd_refSV_type=fd_i.split("=")[1]
	return([fd_svid, fd_refA_chr, fd_refA_pos, fd_refB_chr, fd_refB_pos, fd_refSV_type, fd_refSV_ct, fd_refSV_consensus])


def readmpileup(mpline):
	mpline_ref=0
	mpline_alt=0
	if len(mpline)==0:
		return([])
	mpline=mpline.split('\n')[0]
	mpline=mpline.split('\t')
	for i in mpline[4]:
		if i=="." or i==",":
			mpline_ref=mpline_ref+1
		else:
			mpline_alt=mpline_alt+1
	return([mpline_ref, mpline_alt])

def calculate_theta(inmatrix):
	if len(inmatrix)==0:
		return(None)
	kk=[]; maf=[]
	for iii in inmatrix:
		maf.append(max(iii[0],iii[1])/(iii[0]+ iii[1]))
		tt = iii[2]/lam
		kk.append(3/4*(1-tt**2))
	kkk=0;kks=0
	for iii in range(0, len(kk)):
		kkk=kkk + kk[iii] * maf[iii]
		kks=kks + kk[iii]
	if kks==0:
		return(None)
	else:
		return(kkk/kks)

##### MAIN #####
### read delly in ###
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
	xlist.append(x)
	svlist_id.append(x[0])
	svlist_chrA.append(x[1]); svlist_posA.append(int(float(x[2])))
	svlist_chrB.append(x[3]); svlist_posB.append(int(float(x[4])))
	#if "DEL00001068" in delly_line:
	#       break

### read bulk hSNP in ###
hsnp_chr=[]
hsnp_pos=[]
while True:
	bulkhsnp_line=bulkhsnp_vcf.readline()
	if len(bulkhsnp_line)==0:
		break
	if bulkhsnp_line[0]=='#':
		continue
	tmp=bulkhsnp_line.split("\t")
	hsnp_chr.append(tmp[0])
	hsnp_pos.append(tmp[1])

### for each sv refA or refB identify related hSNP ###
svlist_refAhsnp=[]
svlist_refBhsnp=[]
for i in range(0, len(svlist_id)):
	tmpA=[]
	tmpB=[]
	for j in range(0, len(hsnp_chr)):
		if svlist_chrA[i] == hsnp_chr[j]:
			if float(svlist_posA[i])-lam <= float(hsnp_pos[j]) and float(svlist_posA[i])+lam >= float(hsnp_pos[j]):
				tmpA.append([hsnp_chr[j], hsnp_pos[j]])
		if svlist_chrB[i] == hsnp_chr[j]:
			if float(svlist_posB[i])-lam <= float(hsnp_pos[j]) and float(svlist_posB[i])+lam >= float(hsnp_pos[j]):
				tmpB.append([hsnp_chr[j], hsnp_pos[j]])
	svlist_refAhsnp.append(tmpA)
	svlist_refBhsnp.append(tmpB)


### find number of SV supporting reads in the enhanced alignment ###
### only load data here ###
file_enhancealign=open("pool.support.txt", 'r')
dat_enhancealign = []
while True:
	tmp_line=file_enhancealign.readline()
	if len(tmp_line)==0:
		break
	tmp_line=tmp_line.split("\n")[0]
	tmp_line=tmp_line.split("\t")
	if len(tmp_line)!=5:
		continue
	dat_enhancealign.append(tmp_line)

file_enhancealign.close()

fileout_final=open("pool.vcf","a")

###
for svi in range(0, len(svlist_id)):
	# svi=4
	tmpA=[]
	if len(svlist_refAhsnp)!=0:
		for ii in range(0, len(svlist_refAhsnp[svi])):
			tmp=svlist_refAhsnp[svi][ii][0] + ":" + svlist_refAhsnp[svi][ii][1] + "-" + svlist_refAhsnp[svi][ii][1]
			cmd="samtools mpileup -C50 -q " + str(minmq_hsnp) + " -f " + ref + " -r " + tmp + " " + cell_bam
			tmp=os.popen(cmd); tmp1=tmp.read(); tmp.close()
			tmp2=readmpileup(tmp1)
			if len(tmp2)==0:
				continue
			tmp2.append(abs(int(svlist_refAhsnp[svi][ii][1]) - int(svlist_posA[svi]))); # tmp2: No. ref read; No. alt read; abs distance to SV
			tmpA.append(tmp2)

	tmpB=[]
	if len(svlist_refBhsnp)!=0:
		for ii in range(0, len(svlist_refBhsnp[svi])):
			tmp=svlist_refBhsnp[svi][ii][0] + ":" + svlist_refBhsnp[svi][ii][1] + "-" + svlist_refBhsnp[svi][ii][1]
			cmd="samtools mpileup -C50 -q " + str(minmq_hsnp) + " -f " + ref + " -r " + tmp + " " + cell_bam
			tmp=os.popen(cmd); tmp1=tmp.read(); tmp.close()
			tmp2=readmpileup(tmp1)
			if len(tmp2)==0:
				continue
			tmp2.append(abs(int(svlist_refBhsnp[svi][ii][1]) - int(svlist_posB[svi]))); # tmp2: No. ref read; No. alt read; abs distance to SV
			tmpB.append(tmp2)
	thetaA=calculate_theta(tmpA)
	thetaB=calculate_theta(tmpB)
	if thetaA==None:
		thetaA=theta
	if thetaB==None:
		thetaB=theta
	# [thetaA, thetaB]
	###
	probA0=thetaA*1/8
	probA1=thetaA
	if probA1>0.95: # avoid absolete value
		probA1=0.95
	probB0=thetaB*1/8
	probB1=thetaB
	if probB1>0.95: # avoid absolete value
		probB1=0.95
	#[probA0, probA1, probB0, probB1]
	### find number of SV supporting reads in the enhanced alignment ###
	refAA_rdcount=0
	refASV_rdcount=0
	refBB_rdcount=0
	refBSV_rdcount=0
	for ii in range(0, len(dat_enhancealign)):
		if dat_enhancealign[ii][3]=='exclude' or dat_enhancealign[ii][3]=='unknown':
			continue
		if dat_enhancealign[ii][0]==svlist_id[svi]:
			if dat_enhancealign[ii][1]=='refA' and dat_enhancealign[ii][3]=='refX':
				refAA_rdcount = refAA_rdcount + float(dat_enhancealign[ii][4])
			elif dat_enhancealign[ii][1]=='refA' and dat_enhancealign[ii][3]=='refSV':
				refASV_rdcount = refASV_rdcount + float(dat_enhancealign[ii][4])
			elif dat_enhancealign[ii][1]=='refB' and dat_enhancealign[ii][3]=='refX':
				refBB_rdcount = refBB_rdcount + float(dat_enhancealign[ii][4])
			elif dat_enhancealign[ii][1]=='refB' and dat_enhancealign[ii][3]=='refSV':
				refBSV_rdcount = refBSV_rdcount + float(dat_enhancealign[ii][4])
	result="0/0"
	if (refAA_rdcount + refASV_rdcount + refBB_rdcount + refBSV_rdcount > readlimit) or refAA_rdcount + refASV_rdcount == 0 or refBB_rdcount + refBSV_rdcount == 0:
		result="./."
		outline=[svlist_chrA[svi], str(svlist_posA[svi]), svlist_id[svi], ".", xlist[svi][5], ".", "readlimit", "CHR2="+svlist_chrB[svi]+";"+"END="+str(svlist_posB[svi]), "GT", result]
	else:
		LA0=binom.pmf(k=round(refASV_rdcount), n=round(refAA_rdcount)+round(refASV_rdcount), p=probA0)
		LA1=binom.pmf(k=round(refASV_rdcount), n=round(refAA_rdcount)+round(refASV_rdcount), p=probA1)
		LB0=binom.pmf(k=round(refBSV_rdcount), n=round(refBB_rdcount)+round(refBSV_rdcount), p=probB0)
		LB1=binom.pmf(k=round(refBSV_rdcount), n=round(refBB_rdcount)+round(refBSV_rdcount), p=probB1)
		#[LA0, LA1, LB0, LB1]
		### criteria: LA1*LB1 > LA0*LB0*3
		if LA1*LB1 > LA0*LB0*3:
			result='0/1'
		outline=[svlist_chrA[svi], str(svlist_posA[svi]), svlist_id[svi], ".", xlist[svi][5], ".", "PASS", "CHR2="+svlist_chrB[svi]+";"+"END="+str(svlist_posB[svi]), "GT:thetaA:thetaB:refAAcount:refASVcount:refBBcount:refBSVcount:LA0:LA1:LB0:LB1"]
		outline.append(result+":"+str(round(thetaA, 3))+":"+str(round(thetaB, 3))+":"+str(round(refAA_rdcount))+":"+str(round(refASV_rdcount))+":"+str(round(refBB_rdcount))+":"+str(round(refBSV_rdcount))+":"+str(LA0)+":"+str(LA1)+":"+str(LB0)+":"+str(LB1))
	ooo=outline[0]
	for iii in outline[1:]:
		ooo = ooo + "\t" + str(iii)
	ooo = ooo+'\n'
	fileout_final.write(ooo)

fileout_final.close()


