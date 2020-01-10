# Copyright (C) 2020  Dong X, et. al.
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero General Public License for more details.

# This package also include pre-compiled binary and source code of novoBreak; and pre-compiled binary of ssake.
# Refer to the licenses for the two software tools separately.

import sys
from subprocess import Popen, PIPE
import os
import time
main_path = os.path.split(os.path.realpath(__file__))[0]
sc3_py = os.path.join(main_path, "scripts/sc3.0.0.py")
sc2_py = os.path.join(main_path, "scripts/sccaller2.py")
bwa = os.path.join(main_path, "bin/bwa")


def main():
    if len(sys.argv) == 1:
        print """
Author: Xiao Dong, biosinodx@gmail.com, xiao.dong@einsteinmed.org; Yujue Wang, spsc83@gmail.com
Description: single cell SNP and indel (si), structure variation (sv) caller
Version: 3.0.0
Usage: 
    python sccaller.py <command> <parameters>
Dependence:
    sv: python modules sys, logging, os, re, subprocess, time, copy, struct, itertools, operator, multiprocessing, socket, Queue, collections, numpy, math, datetime
        samtools v.1.9+ (Other versions not tested)
        bwa v.0.7.10-r806-dirty
        SSAKE v3.8.3
    si: python modules os, argparse, sys, subprocess, re, collections, itertools, logging, time, functiontools, random, string, math, numpy, multiprocessing, pysam(0.15.1)
        samtools v.1.9+ (Other versions not tested)
Commands:
    -- Structure variation calling (PEA)
        sv  <file_ref> <cpu_num> <job_size> 
            <single_cell_bam> <single_cell_fq1> <single_cell_fq2> 
            <bulk_bam> <hSNP> <r1> 
            <r2> <output> <ifq_dir>
    -- SNVs and short INDELs calling
        si  -d, --wkdir: work dir. Default: ./
            -l, --lamb: lambda for bias estimation. Default=10000
            --bias: default theta (bias) for SNVs whose theta cannot be estimated. Default=0.75
            --minvar: min. num. variant supporting reads. Default: 4
            --mapq: min. mapQ. Default: 40
            --min_depth: min. reads. Default: 10
            --RD: min. read depth of known heterogous SNP called from bulk when choosing -t dbsnp. Default: 20. Recommand: 10,15,20, depending on average read depth
            --null: min. allelic fraction considered. Default=0.03
            -e, --engine: pileup engine. samtools or pysam. Default: pysam
            -w, --work_num: num. splits per chromosome for multi-process computing. Default: 100
            -n, --cpu_num: num. processes. Default: 1
            --head: first chromosome as sorted as in fasta file to analyze (1-based). Default: the first chr. in the fasta
            --tail: last chromosome as sorted as in fasta file to analyze (1-based). Default: the last chr. in the fasta
            --format: output file format. bed or vcf. Default: vcf
            --bulk: bam file of bulk DNA sequencing
            --bulk_min_var: min. num. variant supporting reads for bulk. Default: 1
            --bulk_min_mapq: min. mapQ for bulk. Default: 20
            --bulk_min_depth: min. reads for bulk. Default: 20
            -h, --help: Help
    """
        exit(0)
    org_dir = os.getcwd()
    # print type(main_path)
    # print "main_path=[{}]".format(main_path)
    if sys.argv[1] in ["SI", "si"]:
        # print sys.argv
        cmd_str = "python {0} {1}".format(sc2_py, " ".join(sys.argv[2:]))
        # print cmd_str
        ph = Popen([cmd_str], shell=True)
        ph.wait()
    elif sys.argv[1] in ["SV", "sv"]:
        # print sc3_py
        if len(sys.argv) == 2:
            cmd_str = "python {}".format(sc3_py)
            ph = Popen([cmd_str], shell=True)
            ph.wait()
            exit(0)
        file_ref = sys.argv[2]
        cpu_num = sys.argv[3]
        job_size = sys.argv[4]
        single_cell_bam = sys.argv[5]
        fq1 = sys.argv[6]
        fq2 = sys.argv[7]
        bulk_bam = sys.argv[8]
        hSNP = sys.argv[9]
        r1 = sys.argv[10]
        r2 = sys.argv[11]
        output = sys.argv[12]
        ifq_dir = sys.argv[13]

        output_dir = "{0}_{1}".format(os.path.splitext(os.path.split(bulk_bam)[1])[0],
                                      os.path.splitext(os.path.split(single_cell_bam)[1])[0])
        ph = Popen(["cd {}".format(main_path)], shell=True)
        ph.wait()

        # novoBreak
        cmd_str = "{0} {1} {2} {3} {4} {5} {6}".format(os.path.join(main_path, "scripts/run_novoBreak.sh"),
                                   main_path, file_ref, single_cell_bam, bulk_bam, cpu_num, output_dir)
        ph = Popen([cmd_str], shell=True)
        ph.wait()

        # sc3
        cmd_str = "python {0} sc3procedure_multiple {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15}" \
                  "".format(sc3_py,
                            bwa,
                            os.path.join(main_path, output_dir, "ssake/sv2alignments"),
                            file_ref,
                            os.path.join(main_path, output_dir, "novoBreak.pass.flt.vcf"),
                            cpu_num,
                            job_size,
                            single_cell_bam,
                            fq1,
                            fq2,
                            hSNP,
                            r1,
                            r2,
                            os.path.join(main_path, output_dir, output),
                            os.path.join(main_path, output_dir, "sc3"),
                            ifq_dir)
        ph = Popen([cmd_str], shell=True)
        ph.wait()
    elif sys.argv[1] == "install":
        print "installing novoBreak ..."
        # ph = Popen(["mkdir -p {}".format(os.path.join(main_path, "tmp"))], shell=True)
        # ph.wait()
        # os.chdir(os.path.join(main_path, "tmp"))
        # ph = Popen(["git clone https://github.com/czc/novobreak_src.git"], shell=True)
        # ph.wait()
        ph = Popen(["cd {} && make".format(os.path.join(main_path, "scripts/novobreak_src"))], shell=True)
        ph.wait()
        ph = Popen(["mv {0} {1}".format(os.path.join(main_path, "scripts/novobreak_src/novoBreak"),
                                        os.path.join(main_path, "bin"))], shell=True)
        ph.wait()
        # ph = Popen(["rm -rf novobreak_src"], shell=True)
        # ph.wait()
    else:
        print "wrong argument"


if __name__ == "__main__":
    main()
