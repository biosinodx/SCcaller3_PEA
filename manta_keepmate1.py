# python 3
import argparse
parser=argparse.ArgumentParser(description="PEA formating data\nBased on manta vcf output")
parser.add_argument("-i","--input",type=str, required=True, help="Input.")
parser.add_argument("-o","--output",type=str, required=True, help="Output.")
args=parser.parse_args()

filein=open(args.input, "r")
fileout=open(args.output, "w")

pastid=[]
while True:
	line=filein.readline()
	if len(line)==0:
		break
	if line[0]=="#":
		continue
	if not ("PASS" in line):
		continue
        # if "IMPRECISE" in line: ### not in use
        #        continue
	if not("MATEID=" in line):
		fileout.write(line)
		continue
	id=line.split("\t")[2]
	tmp=line.split("MATEID=")[1]
	mateid=tmp.split(";")[0]
	if mateid in pastid:
		continue
	pastid.append(id)
	fileout.write(line)

filein.close()
fileout.close()

