# **coding: utf-8**
import sys
import logging
import os
import re
from subprocess import Popen, PIPE
import time
import copy
import struct
from itertools import compress
import operator
from multiprocessing.managers import BaseManager
import socket
import Queue
import multiprocessing
from itertools import combinations

from collections import defaultdict
import numpy as np
import math
from datetime import datetime

SAMSTEP = 100000  # do not modify this number
INDEXSTEP = 1000  # do not modify this number
TAG_SAMPLE = 0
TAG_SAME_NUM = 1
TAG_DIFF = 2
TAG_DIFFHEX = 3

EDGE_LEN = 500
# files
ENHANCE_REF_FILE = "enhanced_reference.fasta"
OVERLAPPED_SEQ_FILE = "overlaped.fq"
ENHANCE_ALIGNMENT_FILE = "enhanced_alignment.sam"
FILTERED_MAX_MAPQ_ALIGNMENT_FILE = "filtered_max_mapQ.sam"
SORTED_FILTERED_MAX_MAPQ_ALIGNMENT_FILE = "sorted_filtered_max_mapQ.sam"
FILTERED_MAX_MAPQ_ALIGNMENT_PILEUP_FILE = "filtered_max_mapQ.pileup"
UNUSUAL_SV_FILE = "unusual"
INDEL = "\+[0-9]+[ACGTNacgtn]+|\-[0-9]+[ACGTNacgtn]+"

NEED_LOG = False

# list the support assembler
assembler_list = ["ssake"]
assembler_set = set(assembler_list)

# list the support mapper
mapper_list = ["bwa"]
mapper_set = set(mapper_list)


class tools:

    def __init__(self, ):
        pass

    @staticmethod
    def clip_right(cigar_str, s_list):
        if re.findall(r"(\d+)M{}[SH]".format(max(s_list)), cigar_str):
            return True
        else:
            return False

    @staticmethod
    def is_hard_clip(cigar_str, is_right_clip):
        if is_right_clip:
            if re.findall(r"(\d+)M(\d+)H", cigar_str):
                return True
            elif re.findall(r"(\d+)M(\d+)S", cigar_str):
                return False
            else:
                raise RuntimeError("The cigar should be right clip. But it does not look like right clip.")

        else:
            if re.findall(r"(\d+)H(\d+)M", cigar_str):
                return True
            elif re.findall(r"(\d+)S(\d+)M", cigar_str):
                return False
            else:
                raise RuntimeError("The cigar should be left clip. But it does not look like left clip.")

    @staticmethod
    # check cigar1 and cigar2, make sure they are correct
    def check_cigar(cigar, org_cigar):
        org_cigar_info = re.findall(r"(\d+)([SMIDH])", org_cigar)
        l = 0
        t = ""
        new_cigar = ""
        match = 0
        for i in org_cigar_info:
            l = int(i[0])
            t = i[1]
            if t == "S" or t == "H":
                if match == 0:
                    new_cigar = "{0}{1}{2}".format(new_cigar, l, t)
                else:
                    new_cigar = "{0}{1}M{2}{3}".format(new_cigar, match, l, t)
                    match = 0
            match += l if t == "M" else 0
            if t == "I" or t == "D":
                if l > 20:
                    new_cigar = "{0}{1}M{2}{3}".format(new_cigar, match, l, t)
                    match = 0
                else:
                    match += (l if t == "I" else 0)
        if match > 0:
            new_cigar = "{0}{1}M".format(new_cigar, match)
        if cigar in new_cigar:
            return cigar
        return -1

    @staticmethod
    def is_reverse(flag):
        assert type(flag) == str
        return int(flag) & 0x10 == 0x10

    @staticmethod
    def is_forward(flag):
        assert type(flag) == str
        return not tools.is_reverse(flag)

    @staticmethod
    def get_error_probability(qual_in):
        if type(qual_in) == str:
            if qual_in.isdigit():
                qual = int(qual_in)
            else:
                qual = float(qual_in)
        else:
            qual = qual_in
        return math.pow(10, float(qual) / 10 * (-1))

    @staticmethod
    def format_region_str(region_str):
        chr, pos_str = region_str.split(":")
        pos1, pos2 = [int(i) for i in pos_str.split("-")]
        return "{0}:{1}-{2}".format(chr, min(pos1, pos2), max(pos1, pos2))

    @staticmethod
    def load_ghdbsnp(file_gsnp):
        ghdbsnp_list = []
        cmd_str = """awk '$3!="."&&(match($10,"0/1")||match($10,"1/0"))&&length($4)==length($5)&&!match($5,",")' {0} | cut -f 1,2,4,5""".format(
            file_gsnp)
        ph = Popen([cmd_str], shell=True, stdout=PIPE)
        while True:
            if not ph.stdout:
                break
            line = ph.stdout.readline()
            if not line:
                break
            else:
                ghdbsnp_list.append(line.strip().split("\t"))
        if NEED_LOG:
            logging.debug("total number of ghdbsnp is {}".format(len(ghdbsnp_list)))
        return ghdbsnp_list

    @staticmethod
    def get_vcf_key(vcf_line_list):
        contig_name = "{0}|{1}|{2}|read{3}|{4}".format(vcf_line_list[10], vcf_line_list[11],
                                                       vcf_line_list[12],
                                                       vcf_line_list[13], vcf_line_list[14])
        chr1 = vcf_line_list[0]
        pos1 = vcf_line_list[1]
        chr2 = re.findall(r"CHR2=(.+?);END=", vcf_line_list[7])[0]
        pos2 = re.findall(r";END=(.+?)\;SVLEN=", vcf_line_list[7])[0]
        vcf_key = "_".join([contig_name, chr1, pos1, chr2, pos2])
        vcf_key2 = "_".join([contig_name, chr2, pos2, chr1, pos1])
        return [vcf_key, vcf_key2]

    @staticmethod
    def get_region_list_from_enhance_ref(enhance_ref_list):
        """
        get the list of region from enhance ref. The region in the list follow the format as below:
        "chr:min_pos-max_pos"
        :param enhance_ref_list:
        :return:
        """
        region_list = []
        for seq_name in enhance_ref_list[::2]:
            region_list.extend(re.findall(r"\|([A-Za-z0-9_.]+:\d+-\d+)", seq_name))
        if len(region_list) != 6:
            logging.debug("enhance ref = {}".format(enhance_ref_list))
        assert len(region_list) == 6
        return list(set(map(lambda x: tools.format_region_str(x), region_list)))

    @staticmethod
    def get_related_ghdbsnp_dict(ghdbsnp_list, enhance_ref_list):
        # type: (list[list[str, str, str, str]], list[str,str,str,str,str,str]) -> dict[str,list[str, str, str, int, str, int, int, int]]
        """
        get the dict in which indicate the related ghdbsnp's information
        :type ghdbsnp_list: list[list[str, str, str, str]]
        :type enhance_ref_list: list[str,str,str,str,str,str]
        :param enhance_ref_list: [name_p, seq_p, namea, seq_a, nameb, seq_b]
                                name:>contig name|[consensuse,refA,refB]|chr1:pos1-pos11|chr2:pos2-pos22
        :param ghdbsnp_list: [[chr, pos, ref, alt], ]
        :return:mut_dict chr pos ref alt : [H, h, chr1, pos1, chr2, pos2, bpos1, bpos2]
        """

        def variance_in_one_region(variance, one_region):
            """

            :rtype: bool
            :type one_region: list[str, str, str]
            :type variance: list[str, str, str, str]
            :param variance: [chr, pos, ref, alt]
            :param one_region: [chr, pos1, pos2]
            :return:
            """
            if variance[0] != one_region[0]:
                return False
            if int(one_region[2]) >= int(one_region[1]):  # 1....2
                if int(variance[1]) < int(one_region[1]) or int(variance[1]) > int(one_region[2]):
                    return False
                return True
            else:  # 2....1
                if int(variance[1]) < int(one_region[2]) or int(variance[1]) > int(one_region[1]):
                    return False
                return True

        def variance_in_regions(variance, regions):
            """

            :type regions: list[list[str, str, str]]
            :type variance: list[str, str, str, str]
            :param variance: [chr, pos, ref, alt]
            :param regions: [[chr, pos1, pos2],]
            :return:
            """
            for region in regions:
                if variance_in_one_region(variance, region):
                    return True
            return False

        def variance_only_in_one_region(variance, regions):
            """
            decide if the variance is only in one of the two regions
            :param variance: [chr, pos, ref, alt]
            :param regions: [[chr1, pos11, pos12],[chr2, pos21, pos22]]
            :return: bool
            """
            return len(filter(lambda x: variance_in_one_region(variance, x), regions)) == 1

        def map_pos(region_list, snp, len_seq):
            # type: (list[list[str, str, str]], list[str, str, str, str], int) -> int
            """
            map the position to enhance reference
            :param region_list:[[chr1,pos1,pos11],[chr2,pos2,pos22]]
            :param snp:[chr,pos,ref,alt]
            :return:
            """
            if variance_in_one_region(snp, region_list[0]):
                ret = abs(int(snp[1]) - int(region_list[0][1])) + 1
            elif variance_in_one_region(snp, region_list[1]):
                # ret = abs(int(snp[1]) - int(region_list[1][1])) + 1 + abs(
                #     int(region_list[0][2]) - int(region_list[0][1])) + 1
                ret = len_seq - abs(int(region_list[1][2]) - int(snp[1]))
            else:
                raise RuntimeError("snp {0} and regions {1} have no relationship".format(snp, region_list))
            return ret

        def need_complement(snp, region_p):
            if variance_in_one_region(snp, region_p[0]):
                return int(region_p[0][1]) > int(region_p[0][2])
            elif variance_in_one_region(snp, region_p[1]):
                return int(region_p[1][1]) > int(region_p[1][2])
            else:
                RuntimeError("snp[{0}] should be in only one region of {1}".format(snp, region_p))

        mut_dict = {}  # type: dict[str, list[str, str, str, int, str, int, int]]

        # [[chr1,pos1,pos11],[chr2,pos2,pos22]]
        region_p = [re.split("[:-]", i) for i in
                    re.findall(r"\|([A-Za-z0-9_.]+:\d+-\d+)", enhance_ref_list[0])]  # type: list[list[str, str, str]]
        if len(region_p) != 2:
            logging.debug("enhance_ref_list = {}".format(enhance_ref_list))
            logging.debug("region_p = {}".format(region_p))
        assert len(region_p) == 2
        related_ghdbsnp_list = filter(lambda x: variance_only_in_one_region(x, region_p),
                                      ghdbsnp_list)  # type: list[list[str, str, str, str]]
        if not related_ghdbsnp_list:
            return mut_dict
        region_a = [re.split("[:-]", i) for i in re.findall(r"\|([A-Za-z0-9_.]+:\d+-\d+)", enhance_ref_list[2])]
        region_b = [re.split("[:-]", i) for i in re.findall(r"\|([A-Za-z0-9_.]+:\d+-\d+)", enhance_ref_list[4])]
        bpos1 = abs(int(region_p[0][2]) - int(region_p[0][1])) + 1
        bpos2 = len(enhance_ref_list[1]) - (abs(int(region_p[1][2]) - int(region_p[1][1])) + 1) + 1

        chr1 = enhance_ref_list[0]
        for related_ghdbsnp in related_ghdbsnp_list:
            key = "\t".join(related_ghdbsnp)  # chr pos ref alt
            need_com = need_complement(related_ghdbsnp, region_p)
            H = complement_seq(related_ghdbsnp[3]) if need_com else related_ghdbsnp[3]  # alt
            h = complement_seq(related_ghdbsnp[2]) if need_com else related_ghdbsnp[2]  # ref
            if variance_in_regions(related_ghdbsnp, region_a):
                chr2 = enhance_ref_list[2]
                region2 = region_a
            elif variance_in_regions(related_ghdbsnp, region_b):
                chr2 = enhance_ref_list[4]
                region2 = region_b
            else:
                continue

            pos1 = map_pos(region_p, related_ghdbsnp, len(enhance_ref_list[1]))  # type: int
            pos2 = map_pos(region2, related_ghdbsnp, len(enhance_ref_list[1]))  # type: int
            # if "E144904|contig1|size477|read3|cov1.33|consensus|15:79988581-79988801|15:79982575-79982827" in chr1:
            #     logging.debug("snp={0} pos1={1} region_p={2} ".format(related_ghdbsnp, pos1, region_p))
            mut_dict[key] = [H, h, chr1.strip(">"), pos1, chr2.strip(">"), pos2, bpos1, bpos2]
            # if key =="2\t114349599\tC\tT":
            #    logging.debug("pos1 = {0} pos2={1}".format(pos1, pos2))
        return mut_dict


def intersection(file_fq1, file_fq2, selected_file1, selected_file2):
    '''
    get the intersection of two fq files
    :param file_fq1:
    :param file_fq2:
    :param selected_file1:
    :param selected_file2:
    :return:
    '''
    with open(file_fq1, "r") as fp1, open(file_fq2, "r") as fp2:
        list1 = fp1.readlines()
        list2 = fp2.readlines()
    list1 = filter(lambda x: len(x) > 1, list1)
    list2 = filter(lambda x: len(x) > 1, list2)

    name1list = [i.split(" ")[0] for i in list1[0::4]]
    name2index1 = {}
    for i in xrange(len(name1list)):
        name2index1[name1list[i]] = i
    first_line1 = list1[0::4]
    seq1list = list1[1::4]
    info1list = list1[2::4]
    quality1list = list1[3::4]

    name2list = [i.split(" ")[0] for i in list2[0::4]]
    name2index2 = {}
    for i in xrange(len(name2list)):
        name2index2[name2list[i]] = i
    first_line2 = list2[0::4]
    seq2list = list2[1::4]
    info2list = list2[2::4]
    quality2list = list2[3::4]

    name_set = set(name1list) & set(name2list)
    name_list = filter(lambda x: x in name_set, name1list)
    with open(selected_file1, "w") as sfp1, open(selected_file2, "w") as sfp2:
        len_name_list = len(name_list)
        counter = 0
        for name in name_list:
            # index1 = name1list.index(name)
            # index2 = name2list.index(name)
            index1 = name2index1[name]
            index2 = name2index2[name]
            sfp1.write("".join([first_line1[index1], seq1list[index1],
                                info1list[index1], quality1list[index1]]))
            sfp2.write("".join([first_line2[index2], seq2list[index2],
                                info2list[index2], quality2list[index2]]))
            if counter % 10000 == 0:
                print "{0} / {1}".format(counter, len_name_list)
            counter += 1


def read_order(flag_in):
    flag_in = int(flag_in)
    if flag_in & 0x40 == 0x40:  # first in pair
        return "+"
    if flag_in & 0x80 == 0x80:  # second in pair
        return "-"
    return ""


def get_reads(sam_file, sv_file, ssake_file, contig2reads_file, out_support_sam):
    """
    根据raw和group文件，把sv中位置处支持以及不支持的read保存在sam文件里面
    """
    # contig name ---> contig data
    print "building contig name to contig data dict"
    with open(ssake_file, "r") as fp:
        contig_data = [i.strip("\n") for i in fp.readlines()]
    contig_dict = dict(zip([i.strip(">") for i in contig_data[0::2]], contig_data[1::2]))  # name : seq
    del contig_data

    # contig name ---> read names
    print "building contig2reads dict"
    contig2reads_dict = {}
    with open(contig2reads_file, 'r') as fp:
        contig2reads_data = fp.read()
    contig2reads_data = contig2reads_data.split(">")
    for i in contig2reads_data:
        tmp = i.split("\n")
        if len(tmp) > 2:
            read_name_list = filter(lambda y: len(y) > 0, map(lambda x: x.split(",")[0], tmp[1:]))
            # print "contig name: {0} read names:{1}".format(tmp[0], read_name_list)
            contig2reads_dict[tmp[0]] = []
            for j in read_name_list:
                contig2reads_dict[tmp[0]].append(j)  # contig name : [read name +/-]
    del contig2reads_data

    # read name ---> read data
    print "building read name to read data dict"
    with open(sam_file, "r") as fp:
        sam_data = fp.readlines()
    sam_head = "".join(filter(lambda x: x.startswith("@"), sam_data))
    sam_data = filter(lambda x: not x.startswith("@"), sam_data)
    read_dict = dict(zip(["{0}{1}".format(i.split("\t")[0], read_order(i.split("\t")[1])) for i in sam_data], sam_data))
    del sam_data

    with open(sv_file, "r") as fp:
        sv_data = fp.readlines()
    sv_data = filter(lambda x: not x.startswith("#"), sv_data)

    sv_line = sv_data[0]
    for sv_line in sv_data:
        sv_list = sv_line.split("\t")
        contig_name = "{0}|{1}|{2}|read{3}|{4}".format(sv_list[10], sv_list[11], sv_list[12], sv_list[13], sv_list[14])
        # group_num = sv_list[10].strip("E")
        # print "contig_name = {}".format(contig_name)
        # print "group_num = {}".format(group_num)
        # read_group_list = filter(lambda x: x[0] == group_num, group_list)  # type: list[list]
        with open(out_support_sam, "w") as fp:
            if contig_name in contig2reads_dict:
                fp.write(sam_head)
                read_name_list = contig2reads_dict[contig_name]
                for read_name in read_name_list:
                    fp.write(read_dict[read_name])
                # handle the reads here
            else:
                logging.warnings("contig {0} don't have supporting reads information".format(contig_name))


def analyze(file_in, threshold, output):
    """
    look for the sv which shorter then threshold
    :param file_in: sv file
    :param threshold:
    :param output:
    :return:
    """
    with open(file_in, "r") as fp_in, open(output, "w") as fp_out:
        while True:
            line = fp_in.readline()
            if not line:
                break
            if line.startswith("#"):
                continue
            line_list = line.strip("\n").split("\t")
            if int(line_list[7].split(";")[-1].split("=")[1]) < int(threshold):
                fp_out.write(line)


def generator_seq_from_ref(fasta_file_name):
    chr_name = "1"
    start = 0
    length = 0
    with open(fasta_file_name, "r") as fp:
        fasta_data = fp.read()
    fasta_list = fasta_data.split(">")
    del fasta_data
    fasta_list = [i.split("\n") for i in fasta_list]
    fasta_dict = dict(
        zip([i[0].split(" ")[0] for i in fasta_list], ["".join(i[1:]) for i in fasta_list]))  # chromesome: sequence
    del fasta_list
    while True:
        assert type(chr_name) == str
        assert type(start) == int
        assert type(length) == int
        [chr_name, start, length] = (yield fasta_dict[chr_name][start - 1:start - 1 + length])


def cigar_reverse(cigar):
    return "".join(re.findall(r"\d+[SMH=]", cigar)[::-1])


def cigar_total_len(cigar):
    return sum([int(i) for i in re.findall(r"(\d+)[SHM=]", cigar)])


def cigar_match_num(cigar):
    return sum([int(i) for i in re.findall(r"(\d+)[M=]", cigar)])


def cigar_clip_num(cigar):
    sh_list = re.findall(r"(\d+)[SH]", cigar)
    if not sh_list:
        return 0
    return max([int(i) for i in sh_list])


def cigar_opposite_clip_info(cigar):
    sh_list = re.findall(r"(\d+)([SH])", cigar)
    if len(sh_list) < 2:
        return [0, ""]
    else:
        assert len(sh_list) == 2
        if int(sh_list[0][0]) < int(sh_list[1][0]):
            return [int(sh_list[0][0]), sh_list[0][1]]
        else:
            return [int(sh_list[1][0]), sh_list[1][1]]


def cigar_forward_clip_info(cigar):
    """
    get the clip side clip number and clip chr(S or H)
    :param cigar:
    :return: [iclip_number, str_clip_chr]
    """
    sh_list = re.findall(r"(\d+)([SH])", cigar)
    if len(sh_list) == 0:
        return [0, ""]
    assert len(sh_list) <= 2
    if len(sh_list) == 1:
        return [int(sh_list[0][0]), sh_list[0][1]]
    elif len(sh_list) == 2:
        if int(sh_list[0][0]) > int(sh_list[1][0]):
            return [int(sh_list[0][0]), sh_list[0][1]]
        else:
            return [int(sh_list[1][0]), sh_list[1][1]]


def seq_opposite_cut_num(cigar):
    n, t = cigar_opposite_clip_info(cigar)
    if not t:
        return 0
    elif t == "H":
        return 0
    elif t == "S":
        return n
    else:
        RuntimeError("cigar type error")


def cigar_is_unusual(cigar):
    assert type(cigar) == str
    if len(re.findall("[IDNPX]", cigar)) > 0:
        return 1
    len_sh = len(re.findall(r"[SH]", cigar))
    if len_sh > 2 or len_sh == 0:
        return 2
    len_m = len(re.findall(r"[M=]", cigar))
    if len_m > 1 or len_m == 0:
        return 3
    return 0


def cigar_is_right_clip(cigar):
    if cigar_is_unusual(cigar) > 0:
        return None
    s_list = map(lambda x: int(x), re.findall(r"(\d+)[SH]", cigar))
    m_list = map(lambda x: int(x), re.findall(r"(\d+)[M=]", cigar))

    s = max(s_list)  # type: int
    m = m_list[0]
    return True if re.findall(r"{0}M{1}[SH]".format(m, s), cigar) else False


def cigar_is_left_clip(cigar):
    if cigar_is_unusual(cigar) > 0:
        return None
    return not cigar_is_right_clip(cigar)


def alignment_reverse(pos, flag, cigar, seq):
    # type: (int, str, str, str) -> list[int, int, str, str]
    r_cigar = cigar_reverse(cigar)
    r_seq = reverse_complement(seq)
    iflag = int(flag)
    r_flag = (iflag - (iflag & 16)) + ((~(iflag & 16)) & 16)
    r_pos = pos + cigar_match_num(cigar) - 1
    return [r_pos, r_flag, r_cigar, r_seq]


def seq_clipside_cut_num(cigar):
    """
    get the seq cut number of the clip side
    hard clip should not cut the sequence. Only soft clip should cut the seq.
    :param cigar:
    :return:
    """
    n, t = cigar_forward_clip_info(cigar)
    if not t:
        return 0
    elif t == "H":
        return 0
    elif t == "S":
        return n
    else:
        RuntimeError("cigar type error")


def parse_alignments2enhance_ref(alignment_info, ref_seq_generator, result_list_out):
    """
    parse "three" alignments into enhanced ref
    :param alignment_info:[contig name, chr1, pos1, flag1, seq1, cigar1, orgcigar1, vcf_pos1,
                                0        1     2      3     4      5        6          7
                           chr2, pos2, flag2, seq2, cigar2, orgcigar2, vcf_pos2
                            8      9    10     11     12       13         14
                           chr_p, pos_p, flag_p, seq_p, cigar_p, orgcigar_p]}
                            15     16      17     18     19        20
    :return:
    -2 cigar can not match, data might be wrong
    -3 there are IDNPX in cigar_p, can not handle
    -4 there are more than 2 SH in cigar_p, can not handle
    -5 there are more than 1 M= in cigar_p, can not handle
    -6 there are IDNPX in cigar_a, can not handle
    -7 there are more than 2 SH in cigar_a, can not handle
    -8 there are more than 1 M= in cigar_a, can not handle
    -9 cigar_a and cigar_p have different total number
    -10 the two cigar can not match with the other
    -11 the alignment might be wrong
    """
    assert len(alignment_info) == 21
    del result_list_out[:]
    is_reversed = False
    contig_name = alignment_info[0]
    [chr_p, pos_p, flag_p, seq_p, cigar_p, orgcigar_p, vcf_pos_p] = alignment_info[1:8] if alignment_info[
                                                                                           1:7] == alignment_info[
                                                                                                   15:] else alignment_info[
                                                                                                             8:15]
    [chr_a, pos_a, flag_a, seq_a, cigar_a, orgcigar_a, vcf_pos_a] = alignment_info[8:15] if alignment_info[
                                                                                            1:7] == alignment_info[
                                                                                                    15:] else alignment_info[
                                                                                                              1:8]
    # logging.debug("alignment info = {}".format(alignment_info))
    assert type(pos_a) == type(pos_p) == str
    pos_a = int(pos_a)
    pos_p = int(pos_p)
    # data might be wrong
    if tools.check_cigar(cigar_p, orgcigar_p) == -1 or tools.check_cigar(cigar_a, orgcigar_a) == -1:
        return -2
    # the cigar can not handle
    ret = cigar_is_unusual(orgcigar_p)
    if ret > 0:
        # logging.debug("orgcigar_p={}".format(orgcigar_p))
        return -2 - ret
    ret = cigar_is_unusual(orgcigar_a)
    if ret > 0:
        # logging.debug("orgcigar_a={}".format(orgcigar_a))
        return -5 - ret
    # unusual case
    if cigar_total_len(orgcigar_p) != cigar_total_len(orgcigar_a):
        return -9

    # impossible case
    if tools.is_forward(flag_p) ^ cigar_is_right_clip(orgcigar_p) == tools.is_forward(flag_a) ^ cigar_is_right_clip(
            orgcigar_a):
        return -10
    if tools.is_forward(flag_p) ^ tools.is_forward(flag_a):
        pos_a, flag_a, orgcigar_a, seq_a = alignment_reverse(pos_a, flag_a, orgcigar_a, seq_a)
        is_reversed = True
    if cigar_is_right_clip(orgcigar_p):
        if not is_reversed:
            ref_name_p = ">{0}|consensus|{1}:{2}-{3}|{4}:{5}-{6}" \
                         "".format(contig_name,
                                   chr_p, pos_p, pos_p + cigar_match_num(orgcigar_p) - 1,
                                   chr_a, pos_a, pos_a + cigar_match_num(orgcigar_a) - 1)
        else:
            if pos_a - cigar_match_num(orgcigar_a) + 1 < 0:
                return -11
            ref_name_p = ">{0}|consensus|{1}:{2}-{3}|{4}:{5}-{6}" \
                         "".format(contig_name,
                                   chr_p, pos_p, pos_p + cigar_match_num(orgcigar_p) - 1,
                                   chr_a, pos_a, pos_a - cigar_match_num(orgcigar_a) + 1)
        ref_seq_p = seq_p[seq_opposite_cut_num(orgcigar_p):-seq_opposite_cut_num(orgcigar_a)] if seq_opposite_cut_num(
            orgcigar_a) > 0 else seq_p[seq_opposite_cut_num(orgcigar_p):]
        ref_name_a = ">{0}|refA|{1}:{2}-{3}|{4}:{5}-{6}" \
                     "".format(contig_name,
                               chr_p, pos_p, pos_p + cigar_match_num(orgcigar_p) - 1,
                               chr_p, pos_p + cigar_match_num(orgcigar_p), pos_p + len(ref_seq_p) - 1)
        ref_seq_a = "{0}{1}".format(
            seq_p[seq_opposite_cut_num(orgcigar_p):seq_opposite_cut_num(orgcigar_p) + cigar_match_num(orgcigar_p)],
            ref_seq_generator.send([chr_p,
                                    pos_p + cigar_match_num(orgcigar_p),
                                    len(ref_seq_p) - cigar_match_num(orgcigar_p)]))
        if not is_reversed:
            if pos_a - len(ref_seq_p) + cigar_match_num(orgcigar_a) < 0 or pos_a - 1 < 0:
                return -11
            ref_name_b = ">{0}|refB|{1}:{2}-{3}|{4}:{5}-{6}" \
                         "".format(contig_name,
                                   chr_a, pos_a - len(ref_seq_p) + cigar_match_num(orgcigar_a), pos_a - 1,
                                   chr_a, pos_a, pos_a + cigar_match_num(orgcigar_a) - 1)
            ref_seq_b = "{0}{1}" \
                        "".format(ref_seq_generator.send([chr_a,
                                                          pos_a - len(ref_seq_p) + cigar_match_num(orgcigar_a),
                                                          len(ref_seq_p) - cigar_match_num(orgcigar_a)]),
                                  seq_a[
                                  seq_clipside_cut_num(orgcigar_a):seq_clipside_cut_num(orgcigar_a) + cigar_match_num(
                                      orgcigar_a)])
        else:
            if pos_a - cigar_match_num(orgcigar_a) + 1 < 0:
                return -11
            ref_name_b = ">{0}|refB|{1}:{2}-{3}|{4}:{5}-{6}" \
                         "".format(contig_name,
                                   chr_a, pos_a + len(ref_seq_p) - cigar_match_num(orgcigar_a), pos_a + 1,
                                   chr_a, pos_a, pos_a - cigar_match_num(orgcigar_a) + 1)
            ref_seq_b = "{0}{1}" \
                        "".format(reverse_complement(ref_seq_generator.send([chr_a,
                                                                             pos_a + 1,
                                                                             len(ref_seq_p) - cigar_match_num(
                                                                                 orgcigar_a)])),
                                  seq_a[
                                  seq_clipside_cut_num(orgcigar_a): seq_clipside_cut_num(orgcigar_a) + cigar_match_num(
                                      orgcigar_a)])

    else:  # p is left clip, a is right clip
        if not is_reversed:
            ref_name_p = ">{0}|consensus|{1}:{2}-{3}|{4}:{5}-{6}" \
                         "".format(contig_name,
                                   chr_a, pos_a, pos_a + cigar_match_num(orgcigar_a) - 1,
                                   chr_p, pos_p, pos_p + cigar_match_num(orgcigar_p) - 1)
        else:
            if pos_a - cigar_match_num(orgcigar_a) + 1 < 0:
                return -11
            ref_name_p = ">{0}|consensus|{1}:{2}-{3}|{4}:{5}-{6}" \
                         "".format(contig_name,
                                   chr_a, pos_a, pos_a - cigar_match_num(orgcigar_a) + 1,
                                   chr_p, pos_p, pos_p + cigar_match_num(orgcigar_p) - 1)
        ref_seq_p = seq_p[seq_opposite_cut_num(orgcigar_a):-seq_opposite_cut_num(orgcigar_p)] if seq_opposite_cut_num(
            orgcigar_p) > 0 else seq_p[seq_opposite_cut_num(orgcigar_a):]

        if pos_p - len(ref_seq_p) + cigar_match_num(orgcigar_p) < 0 or pos_p - 1 < 0:
            return -11
        ref_name_a = ">{0}|refA|{1}:{2}-{3}|{4}:{5}-{6}" \
                     "".format(contig_name,
                               chr_p, pos_p - len(ref_seq_p) + cigar_match_num(orgcigar_p), pos_p - 1,
                               chr_p, pos_p, pos_p + cigar_match_num(orgcigar_p) - 1)
        ref_seq_a = "{0}{1}".format(ref_seq_generator.send([chr_p,
                                                            pos_p - len(ref_seq_p) + cigar_match_num(orgcigar_p),
                                                            len(ref_seq_p) - cigar_match_num(orgcigar_p)]),
                                    seq_p[
                                    seq_clipside_cut_num(orgcigar_p):seq_clipside_cut_num(orgcigar_p) + cigar_match_num(
                                        orgcigar_p)])
        if not is_reversed:
            ref_name_b = ">{0}|refB|{1}:{2}-{3}|{4}:{5}-{6}" \
                         "".format(contig_name,
                                   chr_a, pos_a, pos_a + cigar_match_num(orgcigar_a) - 1,
                                   chr_a, pos_a + cigar_match_num(orgcigar_a), len(ref_seq_p) + pos_a - 1)
            ref_seq_b = "{0}{1}".format(
                seq_a[seq_opposite_cut_num(orgcigar_a):seq_opposite_cut_num(orgcigar_a) + cigar_match_num(orgcigar_a)],
                ref_seq_generator.send(
                    [chr_a, pos_a + cigar_match_num(orgcigar_a), len(ref_seq_p) - cigar_match_num(orgcigar_a)]))
        else:
            if pos_a - cigar_match_num(orgcigar_a) + 1 < 0 or pos_a - len(ref_seq_p) + 1 < 0:
                return -11
            ref_name_b = ">{0}|refB|{1}:{2}-{3}|{4}:{5}-{6}" \
                         "".format(contig_name,
                                   chr_a, pos_a, pos_a - cigar_match_num(orgcigar_a) + 1,
                                   chr_a, pos_a - cigar_match_num(orgcigar_a), pos_a - len(ref_seq_p) + 1)
            ref_seq_b = "{0}{1}".format(
                seq_a[seq_opposite_cut_num(orgcigar_a):seq_opposite_cut_num(orgcigar_a) + cigar_match_num(orgcigar_a)],
                reverse_complement(ref_seq_generator.send([chr_a,
                                                           pos_a - len(ref_seq_p) + 1,
                                                           len(ref_seq_p) - cigar_match_num(orgcigar_a)]))
            )

    result_list_out.extend([ref_name_p, ref_seq_p,
                            ref_name_a, ref_seq_a,
                            ref_name_b, ref_seq_b])
    return 0


def parse_alignments2enhance_ref_org(alignment_info, ref_seq_generator, result_list_out):
    """
    parse "three" alignments into enhanced ref
    :param alignment_info:[contig name, chr1, pos1, flag1, seq1, cigar1, orgcigar1, vcf_pos1,
                                0        1     2      3     4      5        6          7
                           chr2, pos2, flag2, seq2, cigar2, orgcigar2, vcf_pos2
                            8      9    10     11     12       13         14
                           chr_p, pos_p, flag_p, seq_p, cigar_p, orgcigar_p]}
                            15     16      17     18     19        20
    :return:
    """
    assert len(alignment_info) == 21
    contig_name = alignment_info[0]
    chr1 = alignment_info[1]
    pos1 = alignment_info[2]
    flag1 = alignment_info[3]
    seq1 = alignment_info[4]
    cigar1 = alignment_info[5]
    orgcigar1 = alignment_info[6]
    vcf_pos1 = alignment_info[7]
    chr2 = alignment_info[8]
    pos2 = alignment_info[9]
    flag2 = alignment_info[10]
    seq2 = alignment_info[11]
    cigar2 = alignment_info[12]
    orgcigar2 = alignment_info[13]
    vcf_pos2 = alignment_info[14]
    chr_p = alignment_info[15]
    pos_p = alignment_info[16]
    flag_p = alignment_info[17]
    seq_p = alignment_info[18]
    cigar_p = alignment_info[19]
    orgcigar_p = alignment_info[20]
    del result_list_out[:]
    # logging.debug("2 alignments info {}".format(alignment_info))
    # make sure the cigars are correct
    if tools.check_cigar(cigar1, orgcigar1) == -1 or tools.check_cigar(cigar2, orgcigar2) == -1:
        return -1
    groups = [[chr_p, pos_p, flag_p, seq_p, cigar_p, orgcigar_p]]
    if [chr1, pos1, flag1, seq1, cigar1, orgcigar1] not in groups:
        groups.append([chr1, pos1, flag1, seq1, cigar1, orgcigar1])
    if [chr2, pos2, flag2, seq2, cigar2, orgcigar2] not in groups:
        groups.append([chr2, pos2, flag2, seq2, cigar2, orgcigar2])
    # should be only 2 groups
    if len(groups) != 2:
        return -2

    # they should be different alignments
    if groups[0][:-2] == groups[1][:-2] and groups[0][-1] == groups[1][-1]:
        return -3

    s1_list = map(lambda x: int(x), re.findall(r"(\d+)[SH]", groups[0][4]))
    s2_list = map(lambda x: int(x), re.findall(r"(\d+)[SH]", groups[1][4]))
    m1_list = map(lambda x: int(x), re.findall(r"(\d+)M", groups[0][4]))
    m2_list = map(lambda x: int(x), re.findall(r"(\d+)M", groups[1][4]))
    # keep single SH and M
    if len(s1_list) != 1 or len(s2_list) != 1 or len(m1_list) != 1 or len(m2_list) != 1:
        return -4
    s1 = max(s1_list)  # type: int
    m1 = max(m1_list)  # type: int

    s2 = max(s2_list)  # type: int
    m2 = max(m2_list)  # type: int
    is_left_clip1 = True if re.findall(r"{0}[SH]{1}M".format(s1, m1), groups[0][4]) else False
    is_left_clip2 = True if re.findall(r"{0}[SH]{1}M".format(s2, m2), groups[1][4]) else False
    cut_num1 = 0  # the small clip in the other side of group[0]
    cut_num2 = 0  # the small clip in the other side of group[1]
    is_same_strand = (tools.is_reverse(groups[0][2]) == tools.is_reverse(groups[1][2]))
    if (sum(s1_list) + sum(m1_list)) != (sum(s2_list) + sum(m2_list)):
        cut_list = re.findall(r"(\d+)[SH]$" if is_left_clip1 else r"^(\d+)[SH]", groups[0][5])
        if cut_list:
            cut_num1 = int(cut_list[0])  # type: int
        cut_list = re.findall(r"(\d+)[SH]$" if is_left_clip2 else r"^(\d+)[SH]", groups[1][5])
        if cut_list:
            cut_num2 = int(cut_list[0])  # type: int
        if cut_num1 + sum(s1_list) + sum(m1_list) != sum(s2_list) + sum(m2_list) + cut_num2:
            return -7


    # if tools.clip_right(groups[0][4], s1_list) and not tools.is_reverse(groups[0][2]):  # p: right clip and +
    #     if tools.clip_right(groups[1][4], s2_list) and not tools.is_reverse(groups[1][2]):  # right clip and +
    #         return -11
    #     elif not tools.clip_right(groups[1][4], s2_list) and tools.is_reverse(groups[1][2]):  # left clip and -
    #         return -12
    #     elif tools.clip_right(groups[1][4], s2_list) and tools.is_reverse(groups[1][2]):  # right clip and -
    #         cut_num_left = cut_num1
    #         cut_num_right = cut_num2
    #         if cut_num_right + cut_num_left > 0 and NEED_LOG:
    #             logging.debug("p: right clip and +")
    #             logging.debug("right clip and -")
    #         tmp2 = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) + m2 - 1, int(groups[1][1]) + cut_num_right)
    #         tmpB = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) + m2 + s2 - 1 - cut_num_left,
    #                                     int(groups[1][1]) + m2)
    #         ref_seq_B = reverse_seq(
    #             "{0}{1}".format(groups[1][3] if tools.is_hard_clip(groups[1][4], True) else groups[1][3][:-1 * s2],
    #                             ref_seq_generator.send([groups[1][0], int(groups[1][1]) + m2, s2])))
    #     else:  # left clip and +
    #         cut_num_left = cut_num1
    #         cut_num_right = cut_num2
    #         if cut_num_right + cut_num_left > 0 and NEED_LOG:
    #             logging.debug("p: right clip and +")
    #             logging.debug("left clip and +")
    #         tmp2 = "{0}:{1}-{2}".format(groups[1][0], groups[1][1], int(groups[1][1]) + m2 - 1 - cut_num_right)
    #         tmpB = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) - s2 + cut_num_left, int(groups[1][1]) - 1)
    #
    #         ref_seq_B = "{0}{1}".format(ref_seq_generator.send([groups[1][0], int(groups[1][1]) - s2, s2]),
    #                                     groups[1][3] if tools.is_hard_clip(groups[1][4], False) else groups[1][3][s2:])
    #
    #     tmp1 = "{0}:{1}-{2}".format(groups[0][0], int(groups[0][1]) + cut_num_left, int(groups[0][1]) + m1 - 1)
    #     ref_name_p = ">{0}|consensus|{1}|{2}".format(contig_name, tmp1, tmp2)
    #     ref_name1 = ">{0}|refA|{1}|{2}:{3}-{4}".format(contig_name, tmp1, groups[0][0], int(groups[0][1]) + m1,
    #                                                    int(groups[0][1]) + m1 + s1 - 1 - cut_num_right)
    #     ref_name2 = ">{0}|refB|{1}|{2}".format(contig_name, tmpB, tmp2)
    #     ref_seq_p = groups[0][3]
    #     ref_seq_A = "{0}{1}".format(groups[0][3] if tools.is_hard_clip(groups[0][4], True) else groups[0][3][:-1 * s1],
    #                                 ref_seq_generator.send([groups[0][0], int(groups[0][1]) + m1, s1]))

    # elif tools.clip_right(groups[0][4], s1_list) and tools.is_reverse(groups[0][2]):  # p: right clip and -
    #     if tools.clip_right(groups[1][4], s2_list) and tools.is_reverse(groups[1][2]):  # right clip and -
    #         return -13
    #     elif not tools.clip_right(groups[1][4], s2_list) and not tools.is_reverse(groups[1][2]):  # left clip and +
    #         return -14
    #     elif not tools.clip_right(groups[1][4], s2_list) and tools.is_reverse(groups[1][2]):  # left clip and -
    #         cut_num_left = cut_num2
    #         cut_num_right = cut_num1
    #         if cut_num_right + cut_num_left > 0 and NEED_LOG:
    #             logging.debug("p: right clip and -")
    #             logging.debug("left clip and -")
    #         tmp2 = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) + m2 - 1 - cut_num_left, groups[1][1])
    #         tmpB = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) - 1, int(groups[1][1]) - s2 + cut_num_right)
    #         ref_seq_B = reverse_seq("{0}{1}".format(ref_seq_generator.send([groups[1][0], int(groups[1][1]) - s2, s2]),
    #                                                 groups[1][3] if tools.is_hard_clip(groups[1][4], False) else
    #                                                 groups[1][3][s2:]))
    #     else:  # right clip and +
    #         cut_num_left = cut_num2
    #         cut_num_right = cut_num1
    #         if cut_num_right + cut_num_left > 0 and NEED_LOG:
    #             logging.debug("p: right clip and -")
    #             logging.debug("right clip and +")
    #         tmp2 = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) + cut_num_left, int(groups[1][1]) + m2 - 1)
    #         tmpB = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) + m2,
    #                                     int(groups[1][1]) + m2 + s2 - 1 - cut_num_right)
    #         ref_seq_B = "{0}{1}".format(
    #             groups[1][3] if tools.is_hard_clip(groups[1][4], True) else groups[1][3][:-1 * s2],
    #             ref_seq_generator.send([groups[1][0], int(groups[1][1]) + m2, s2]))
    #     tmp1 = "{0}:{1}-{2}".format(groups[0][0], int(groups[0][1]) + m1 - 1, int(groups[0][1]) + cut_num_right)
    #     ref_name_p = ">{0}|consensus|{1}|{2}".format(contig_name, tmp2, tmp1)
    #     ref_name1 = ">{0}|refA|{1}:{2}-{3}|{4}".format(contig_name, groups[0][0],
    #                                                    int(groups[0][1]) + m1 + s1 - 1 - cut_num_left,
    #                                                    int(groups[0][1]) + m1,
    #                                                    tmp1)
    #     ref_name2 = ">{0}|refB|{1}|{2}".format(contig_name, tmp2, tmpB)
    #     ref_seq_p = reverse_seq(groups[0][3])
    #     ref_seq_A = reverse_seq(
    #         "{0}{1}".format(groups[0][3] if tools.is_hard_clip(groups[0][4], True) else groups[0][3][:-1 * s1],
    #                         ref_seq_generator.send([groups[0][0], int(groups[0][1]) + m1, s1])))
    elif tools.clip_right(groups[0][4], s1_list) and tools.is_reverse(groups[0][2]):  # p: right clip and -
        if tools.clip_right(groups[1][4], s2_list) and tools.is_reverse(groups[1][2]):  # right clip and -
            return -13
        elif not tools.clip_right(groups[1][4], s2_list) and not tools.is_reverse(groups[1][2]):  # left clip and +
            return -14
        elif not tools.clip_right(groups[1][4], s2_list) and tools.is_reverse(groups[1][2]):  # left clip and -
            cut_num_left = cut_num1
            cut_num_right = cut_num2
            if cut_num_right + cut_num_left > 0 and NEED_LOG:
                logging.debug("p: right clip and -")
                logging.debug("left clip and -")
            tmp2 = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) + m2 - 1 - cut_num_left, groups[1][1])
            tmpB = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) - 1, int(groups[1][1]) - s2 + cut_num_right)
            ref_seq_B = reverse_seq("{0}{1}".format(ref_seq_generator.send([groups[1][0], int(groups[1][1]) - s2, s2]),
                                                    groups[1][3] if tools.is_hard_clip(groups[1][4], False) else
                                                    groups[1][3][s2:]))
            ref_name2 = ">{0}|refB|{1}:{2}-{3}|{4}:{5}-{6}".format(contig_name, groups[1][0], int(groups[1][1]) - s2)
        else:  # right clip and +
            cut_num_left = cut_num2
            cut_num_right = cut_num1
            if cut_num_right + cut_num_left > 0 and NEED_LOG:
                logging.debug("p: right clip and -")
                logging.debug("right clip and +")
            tmp2 = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) + cut_num_left, int(groups[1][1]) + m2 - 1)
            tmpB = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) + m2,
                                        int(groups[1][1]) + m2 + s2 - 1 - cut_num_right)
            ref_seq_B = "{0}{1}".format(
                groups[1][3] if tools.is_hard_clip(groups[1][4], True) else groups[1][3][:-1 * s2],
                ref_seq_generator.send([groups[1][0], int(groups[1][1]) + m2, s2]))
        tmp1 = "{0}:{1}-{2}".format(groups[0][0], int(groups[0][1]) + m1 - 1, int(groups[0][1]) + cut_num_right)
        ref_name_p = ">{0}|consensus|{1}|{2}".format(contig_name, tmp2, tmp1)
        ref_name1 = ">{0}|refA|{1}:{2}-{3}|{4}".format(contig_name, groups[0][0],
                                                       int(groups[0][1]) + m1 + s1 - 1 - cut_num_left,
                                                       int(groups[0][1]) + m1,
                                                       tmp1)

        ref_seq_p = reverse_seq(groups[0][3])
        ref_seq_A = reverse_seq(
            "{0}{1}".format(groups[0][3] if tools.is_hard_clip(groups[0][4], True) else groups[0][3][:-1 * s1],
                            ref_seq_generator.send([groups[0][0], int(groups[0][1]) + m1, s1])))
    elif not tools.clip_right(groups[0][4], s1_list) and not tools.is_reverse(groups[0][2]):  # p: left clip and +
        if tools.clip_right(groups[1][4], s2_list) and tools.is_reverse(groups[1][2]):  # right clip and -
            return -15
        elif not tools.clip_right(groups[1][4], s2_list) and not tools.is_reverse(groups[1][2]):  # left clip and +
            return -16
        elif tools.clip_right(groups[1][4], s2_list) and not tools.is_reverse(groups[1][2]):  # right clip and +
            cut_num_left = cut_num2
            cut_num_right = cut_num1
            if cut_num_right + cut_num_left > 0 and NEED_LOG:
                logging.debug("p: left clip and +")
                logging.debug("right clip and +")
            tmp2 = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) + cut_num_left, int(groups[1][1]) + m2 - 1)
            tmpB = "{0}:{1}-{2}" \
                   "".format(groups[1][0], int(groups[1][1]) + m2, int(groups[1][1]) + m2 + s2 - 1 - cut_num_right)
            ref_seq_B = "{0}{1}".format(
                groups[1][3] if tools.is_hard_clip(groups[1][4], True) else groups[1][3][:-1 * s2],
                ref_seq_generator.send([groups[1][0], int(groups[1][1]) + m2, s2]))
        else:  # left clip and -
            cut_num_left = cut_num2
            cut_num_right = cut_num1
            if cut_num_right + cut_num_left > 0 and NEED_LOG:
                logging.debug("p: left clip and +")
                logging.debug("left clip and -")
            tmp2 = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) + m2 - 1 - cut_num_left, groups[1][1])
            tmpB = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) - 1, int(groups[1][1]) - s2 + cut_num_right)
            ref_seq_B = reverse_seq("{0}{1}".format(ref_seq_generator.send([groups[1][0], int(groups[1][1]) - s2, s2]),
                                                    groups[1][3] if tools.is_hard_clip(groups[1][4], False) else
                                                    groups[1][3][s2:]))
        tmp1 = "{0}:{1}-{2}".format(groups[0][0], groups[0][1], int(groups[0][1]) + m1 - 1 - cut_num_right)
        ref_name_p = ">{0}|consensus|{1}|{2}".format(contig_name, tmp2, tmp1)
        ref_name1 = ">{0}|refA|{1}:{2}-{3}|{4}".format(contig_name, groups[0][0],
                                                       int(groups[0][1]) - s1 + cut_num_left, int(groups[0][1]) - 1,
                                                       tmp1)
        ref_name2 = ">{0}|refB|{1}|{2}".format(contig_name, tmp2, tmpB)
        ref_seq_p = groups[0][3]
        ref_seq_A = "{0}{1}".format(ref_seq_generator.send([groups[0][0], int(groups[0][1]) - s1, s1]),
                                    groups[0][3] if tools.is_hard_clip(groups[0][4], False) else groups[0][3][s1:])

    else:  # p: left clip and -
        if tools.clip_right(groups[1][4], s2_list) and not tools.is_reverse(groups[1][2]):  # right clip and +
            return -17
        elif not tools.clip_right(groups[1][4], s2_list) and tools.is_reverse(groups[1][2]):  # left clip and -
            return -18
        elif tools.clip_right(groups[1][4], s2_list) and tools.is_reverse(groups[1][2]):  # right clip and -
            cut_num_left = cut_num1
            cut_num_right = cut_num2
            if cut_num_right + cut_num_left > 0 and NEED_LOG:
                logging.debug("p: left clip and -")
                logging.debug("right clip and -")
            tmp2 = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) + m2 - 1, int(groups[1][1]) + cut_num_right)
            tmpB = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) + m2 + s2 - 1 - cut_num_left,
                                        int(groups[1][1]) + m2)
            ref_seq_B = reverse_seq(
                "{0}{1}".format(groups[1][3] if tools.is_hard_clip(groups[1][4], True) else groups[1][3][:-1 * s2],
                                ref_seq_generator.send([groups[1][0], int(groups[1][1]) + m2, s2])))

        else:  # left clip and +
            cut_num_left = cut_num1
            cut_num_right = cut_num2
            if cut_num_right + cut_num_left > 0 and NEED_LOG:
                logging.debug("p: left clip and -")
                logging.debug("left clip and +")
            tmp2 = "{0}:{1}-{2}".format(groups[1][0], groups[1][1], int(groups[1][1]) + m2 - 1 - cut_num_right)
            tmpB = "{0}:{1}-{2}".format(groups[1][0], int(groups[1][1]) - s2 + cut_num_left, int(groups[1][1]) - 1)
            ref_seq_B = "{0}{1}".format(ref_seq_generator.send([groups[1][0], int(groups[1][1]) - s2, s2]),
                                        groups[1][3] if tools.is_hard_clip(groups[1][4], False) else groups[1][3][s2:])

        tmp1 = "{0}:{1}-{2}".format(groups[0][0], int(groups[0][1]) + m1 - 1 - cut_num_left, groups[0][1])
        ref_name_p = ">{0}|consensus|{1}|{2}".format(contig_name, tmp1, tmp2)
        ref_name1 = ">{0}|refA|{1}|{2}:{3}-{4}".format(contig_name, tmp1, groups[0][0], int(groups[0][1]) - 1,
                                                       int(groups[0][1]) - s1 + cut_num_right)
        ref_name2 = ">{0}|refB|{1}|{2}".format(contig_name, tmpB, tmp2)
        ref_seq_p = reverse_seq(groups[0][3])
        ref_seq_A = reverse_seq("{0}{1}"
                                "".format(ref_seq_generator.send([groups[0][0], int(groups[0][1]) - s1, s1]),
                                          groups[0][3] if tools.is_hard_clip(groups[0][4], False) else
                                          groups[0][3][s1:]))
    # make sure the break point are in the same position
    if m1 != s2 or m2 != s1:
        return -8
    if cut_num_left + cut_num_right > 0:
        logging.debug("cut_num_left={0} cut_num_right={1}".format(cut_num_left, cut_num_right))
    ref_seq_p = ref_seq_p[cut_num_left:-1 * cut_num_right if cut_num_right > 0 else len(ref_seq_p)]
    ref_seq_A = ref_seq_A[cut_num_left:-1 * cut_num_right if cut_num_right > 0 else len(ref_seq_A)]
    ref_seq_B = ref_seq_B[cut_num_left:-1 * cut_num_right if cut_num_right > 0 else len(ref_seq_B)]

    result_list_out.extend([ref_name_p, ref_seq_p,
                            ref_name1, ref_seq_A,
                            ref_name2, ref_seq_B])
    if cut_num_left + cut_num_right > 0:
        logging.debug("result_list_out = {}".format(result_list_out))
        logging.debug("-----------")
    return 0


def get_overlaped_seq(bam_file, region_list, fq1_source,
                      fq2_source, overlaped_read_file):
    """
    get the overlaped alignments from bam_file, then get the seqs from fq files according to the overlaped alignments
    :param bam_file: input
    :param region_list: enhanced reference's regions, input. "chr:min_pos-max_pos"
    :param fq1_source: input
    :param fq2_source: iput
    :param overlaped_read_file: output
    :return:
    """
    # print "get_overlaped_seq region_list = {}".format(region_list)
    sam_name_list1 = []
    sam_name_list2 = []
    sam_str_list = []  # type: list[str]
    for region in region_list:
        data_source = data_generator_samtools(bam_file, region)
        while True:
            sam_list = data_source.next()
            if not sam_list:
                break
            sam_str_list.append("\t".join(sam_list))

            # sam_list = sam_list.split("\t")
            if int(sam_list[1]) & 64:  # first in pair
                sam_name_list1.append(sam_list[0])
            elif int(sam_list[1]) & 128:  # second in pair
                sam_name_list2.append(sam_list[0])
            else:
                logging.debug(
                    "neither first read nor second read. name={0} flag={1}".format(sam_list[0], sam_list[1]))
    # reduce duplication
    sam_name_list1 = list(set(sam_name_list1))
    sam_name_list2 = list(set(sam_name_list2))
    # sam_str_list = list(set(sam_str_list))
    if (len(sam_name_list1) + len(sam_name_list2)) == 0:
        return -1
    with open(overlaped_read_file, "w") as fp:
        tmp = ""
        for name in sam_name_list1:
            tmp = "{0}{1}".format(tmp, fq1_source.send(name))
            # fp.write(fq1_source.send(name))
        fp.write(tmp)
        tmp = ""
        for name in sam_name_list2:
            tmp = "{0}{1}".format(tmp, fq2_source.send(name))
            # fp.write(fq2_source.send(name))
        fp.write(tmp)
    return 0

    # # handle the sam head
    # cmd_str = "rm -f {}".format(overlaped_alignment_file)
    # process_rm = Popen([cmd_str], shell=True)
    # process_rm.wait()
    # cmd_str = "samtools view -H {0} > {1}".format(bam_file, overlaped_alignment_file)
    # process_samtools = Popen([cmd_str], shell=True)
    # process_samtools.wait()
    # # handle the content
    # with open(overlaped_alignment_file, "a") as fp:
    #     fp.write("\n".join(sam_str_list))


def format_read_base(read_base):
    # remove '$' and '^'
    ret = re.sub("\^\S{1}", "", re.sub("(?<!\^)\$", "", read_base))
    # remove indel
    indel_list = re.findall(INDEL, ret)
    for i in xrange(len(indel_list)):
        indel_len = int(re.findall("[\+-](\d+)[ACGTNacgtn*]", indel_list[i])[0])
        indel_list[i] = indel_list[i][:1 + len(str(indel_len)) + indel_len]
        ret = re.sub("\\" + indel_list[i], "", ret)
    return ret.upper()


# def get_p_value(w, n, p):
#     def binomial_distribution(k, n, p):
#         return math.factorial(n) / math.factorial(k) / math.factorial(n - k) * math.pow(p, k) * math.pow((1 - p),
#                                                                                                          (n - k))
#
#     ret = 1.0
#     for i in xrange(w):
#         ret -= binomial_distribution(i, n, p)
#     return ret

def binomial_distribution(k, n, p):
    ret = 1.0
    for i in xrange(1, n + 1, 1):
        if i <= k:
            ret *= (n - i + 1) * p / i
        else:
            ret *= (n - i + 1) * (1 - p) / (n - i + 1)
    return ret


def get_p_value(w, n, p):
    ret = 0.0
    for i in xrange(w):
        ret += binomial_distribution(i, n, p)
    return ret


def parse_pileup(pileup_file, chr_sv, pos_sv, chr2, pos2, ref, alt, filtered_max_mapQ_alignments_list, r1, r2,
                 ret_list):
    """
    :param r2: min number of reads do not support snp
    :param filtered_max_mapQ_alignments_list: list[max mapQ alignment list]
    :param r1: min number of reads support snp
    :type filtered_max_mapQ_alignments_list: list[list]
    """

    def average_map_false_positive(alignments_list, ch):
        error_list = map(lambda k: tools.get_error_probability(k[4]), filter(lambda x: x[2] == ch, alignments_list))
        if len(error_list) > 0:
            return float(sum(error_list)) / len(error_list)
        else:
            RuntimeError("there is no {} in alignments_list".format(ch))

    def average_base_false_positive(selected_pileup_dict, target_base):
        """

        :type target_base: str
        """

        def get_related_mapq_baseq(selected_pileup_dict, target_chr, target_base):
            """
            get the related mapq list and related baseq list
            :rtype: list[list, list]
            """
            readbase_str = format_read_base(selected_pileup_dict[target_chr][4])
            if NEED_LOG:
                logging.debug("formatted read base = {}".format(readbase_str))
                logging.debug("target_chr={0} target_base={1}".format(target_chr, target_base))
                logging.debug(
                    "len baseq list = {0} len mapq list = {1}".format(len(selected_pileup_dict[target_chr][5]),
                                                                      len(selected_pileup_dict[target_chr][6])))
            baseq_list = [ord(i) - 33 for i in selected_pileup_dict[target_chr][5]]
            mapq_list = [ord(i) - 33 for i in selected_pileup_dict[target_chr][6]]

            assert len(readbase_str) == len(baseq_list) == len(mapq_list)
            filter_list = [i == target_base for i in readbase_str]
            baseq_list = list(compress(baseq_list, filter_list))
            mapq_list = list(compress(mapq_list, filter_list))
            return [baseq_list, mapq_list]

        target_base = target_base.upper()
        baseq_sv_list, mapq_sv_list = get_related_mapq_baseq(selected_pileup_dict, "sv", target_base)
        baseq_list, mapq_list = get_related_mapq_baseq(selected_pileup_dict, "a", target_base)
        baseq_list.extend(baseq_sv_list)
        mapq_list.extend(mapq_sv_list)
        assert len(baseq_list) == len(mapq_list)
        if NEED_LOG and len(baseq_list) == 0:
            logging.debug("selected_pileup_dict={}".format(selected_pileup_dict))
            logging.debug("baseq_list={0} mapq_list={1}".format(baseq_list, mapq_list))
        if len(baseq_list) > 0:
            return sum(
                [1 - (1 - tools.get_error_probability(baseq_list[i])) * (1 - tools.get_error_probability(mapq_list[i]))
                 for
                 i in xrange(len(baseq_list))]) / len(baseq_list)
        else:
            return None

    def predicted_proportion(x, y, n, pa, pb):
        """
        equation (1)
        :param x:
        :param y:
        :param n:
        :param pa:
        :param pb:
        :return:
        """
        return (pb * x / n + pa * y / n - pa * pb) / (1 - pa - pb)

    selected_pileup_dict = {}
    with open(pileup_file, "r") as fp:
        while True:
            pileup_line = fp.readline().strip()
            if not pileup_line:
                break
            pileup_list = pileup_line.split("\t")
            if pileup_list[0] == chr_sv and pileup_list[1] == str(pos_sv):
                if "sv" not in selected_pileup_dict:
                    selected_pileup_dict["sv"] = pileup_list
                else:
                    RuntimeError("unusual pileup file, there should be only one {0}_{1} in pileup.".format(chr_sv,
                                                                                                           pos_sv))
                continue
            if pileup_list[0] == chr2 and pileup_list[1] == str(pos2):
                if "a" not in selected_pileup_dict:
                    selected_pileup_dict["a"] = pileup_list
                else:
                    RuntimeError("unusual pileup file, there should be only one {0}_{1} in pileup".format(chr2,
                                                                                                          pos2))
                continue
    if len(selected_pileup_dict) != 2:
        if NEED_LOG:
            logging.debug("chr_sv={0} pos_sv={1} chr2={2} pos2={3}".format(chr_sv, pos_sv, chr2, pos2))
            logging.debug("selected_pileup_dict = {}".format(selected_pileup_dict))
        return -1
    ref = ref.upper()  # type: str
    alt = alt.upper()  # type: str
    n_H = max([selected_pileup_dict[i][4].upper().count(alt) for i in selected_pileup_dict])
    n_h = max([selected_pileup_dict[i][4].upper().count(ref) for i in selected_pileup_dict])
    if NEED_LOG:
        logging.debug("n_H={0} n_h={1} r1={2} r2={3}".format(n_H, n_h, r1, r2))
        logging.debug("pileup_dict={}".format(selected_pileup_dict))
    # if n_H == n_h:
    #     return -2
    if n_H < r1:
        return -3
    if n_h < r2:
        return -4
    n = len(filtered_max_mapQ_alignments_list)
    if NEED_LOG:
        logging.debug("n={}".format(n))
    if n_H == selected_pileup_dict["sv"][4].upper().count(alt):
        if NEED_LOG:
            logging.debug("case a: H is in ref_sv")
            logging.debug("need to testify h and ref_sv are mutually exclusive")
            logging.debug("A:read does't support snp, B:read is mapped to ref_sv")

        x1 = selected_pileup_dict["a"][4].upper().count(ref)
        y1 = n_H
        w1 = selected_pileup_dict["sv"][4].upper().count(ref)
        false_positive_B1 = average_map_false_positive(filtered_max_mapQ_alignments_list, chr_sv)
        false_positive_A1 = average_base_false_positive(selected_pileup_dict, ref)
        if false_positive_A1 is None:
            return -5
        proportion1 = predicted_proportion(x1, y1, n, false_positive_A1, false_positive_B1)
        if NEED_LOG:
            logging.debug("x1={0} y1={1} w1={2} fpB1={3} fpA1={4} proportion1={5}".format(x1, y1, w1,
                                                                                          false_positive_B1,
                                                                                          false_positive_A1,
                                                                                          proportion1))
        p_value1 = get_p_value(w1, n, proportion1)
        if NEED_LOG:
            logging.debug("p_value1={}".format(p_value1))
            logging.debug("need to testify H and ref_a are mutually exclusive")
            logging.debug("A: read support snp, B: read is mapped to ref_a")
        x2 = n_H
        y2 = selected_pileup_dict["a"][4].upper().count(ref)
        w2 = selected_pileup_dict["a"][4].upper().count(alt)
        false_positive_B2 = average_map_false_positive(filtered_max_mapQ_alignments_list, chr2)
        false_positive_A2 = average_base_false_positive(selected_pileup_dict, alt)
        if false_positive_A2 is None:
            return -5
        proportion2 = predicted_proportion(x2, y2, n, false_positive_A2, false_positive_B2)
        if NEED_LOG:
            logging.debug("x2={0} y2={1} w2={2} fpB2={3} fpA2={4} proportion2={5}".format(x2, y2, w2,
                                                                                          false_positive_B2,
                                                                                          false_positive_A2,
                                                                                          proportion2))
        p_value2 = get_p_value(w2, n, proportion2)
        if NEED_LOG:
            logging.debug("p_value2={}".format(p_value2))
    elif n_H == selected_pileup_dict["a"][4].upper().count(alt):
        if NEED_LOG:
            logging.debug("case b: H is in ref_a")
            logging.debug("need to testify H and ref_sv are mutually exclusive")
            logging.debug("A:read support snp, B: read is mapped to ref_sv")

        x1 = n_H
        y1 = selected_pileup_dict["sv"][4].upper().count(ref)
        w1 = selected_pileup_dict["sv"][4].upper().count(alt)
        false_positive_B1 = average_map_false_positive(filtered_max_mapQ_alignments_list, chr_sv)
        false_positive_A1 = average_base_false_positive(selected_pileup_dict, alt)
        if false_positive_A1 is None:
            return -5
        proportion1 = predicted_proportion(x1, y1, n, false_positive_A1, false_positive_B1)
        if NEED_LOG:
            logging.debug("x1={0} y1={1} w1={2} fpB1={3} fpA1={4} proportion1={5}".format(x1, y1, w1,
                                                                                          false_positive_B1,
                                                                                          false_positive_A1,
                                                                                          proportion1))
        p_value1 = get_p_value(w1, n, proportion1)
        if NEED_LOG:
            logging.debug("p_value1={}".format(p_value1))
            logging.debug("need to testify h and ref_a are mutually exclusive")
            logging.debug("A: read doesn't support snp, B: read is mapped to ref_a")
        x2 = selected_pileup_dict["sv"][4].upper().count(ref)
        y2 = n_H
        w2 = selected_pileup_dict["a"][4].upper().count(ref)
        false_positive_B2 = average_map_false_positive(filtered_max_mapQ_alignments_list, chr2)
        false_positive_A2 = average_base_false_positive(selected_pileup_dict, ref)
        if false_positive_A2 is None:
            return -5
        proportion2 = predicted_proportion(x2, y2, n, false_positive_A2, false_positive_B2)
        if NEED_LOG:
            logging.debug("x2={0} y2={1} w2={2} fpB2={3} fpA2={4} proportion2={5}".format(x2, y2, w2,
                                                                                          false_positive_B2,
                                                                                          false_positive_A2,
                                                                                          proportion2))
        p_value2 = get_p_value(w2, n, proportion2)
        if NEED_LOG:
            logging.debug("p_value2={}".format(p_value2))
    else:
        logging.error("n_H={0} selected_pileup_dict={1}".format(n_H, selected_pileup_dict))
        RuntimeError("unsual n_H")
    if NEED_LOG:
        logging.debug("selected_pileup_dict={}".format(selected_pileup_dict))
    if p_value1 < 0 or p_value1 > 1:
        p_value1 = 1
    if p_value2 < 0 or p_value2 > 1:
        p_value2 = 1
    ret = map(lambda k: str(k),
              [n, x1, y1, w1, false_positive_A1, false_positive_B1, p_value1, x2, y2, w2, false_positive_A2,
               false_positive_B2, p_value2, max(p_value1, p_value2)])
    if NEED_LOG:
        logging.debug("parse_pileup ret = {}".format(ret))
    del ret_list[:]
    ret_list.extend(ret)
    return 0


def position_in_alignment(chr, pos, alignment):
    # type: (str, int, list) -> bool
    assert type(chr) == str
    assert type(pos) == int
    assert type(alignment) == list
    # logging.debug("chr=[{0}] pos=[{1}] alignment={2}".format(chr, pos, alignment))
    a_chr = alignment[2]
    if chr != a_chr:
        # logging.debug("not in alignment different chr")
        return False
    a_pos1 = int(alignment[3])
    cigar = alignment[5]  # type: str
    a_pos2 = a_pos1 + sum([int(i[0]) for i in re.findall("(\d+)(M|D|N)", cigar)]) - 1
    if a_pos1 <= pos <= a_pos2:
        return True
    # logging.debug("not in alignment not in region. pos={0} pos1={1} pos2={2}".format(pos, a_pos1, a_pos2))
    return False


def extend_edge(ref_seq_generator, enhance_ref_list, edge_num):
    assert type(edge_num) == int
    region_p = re.findall(r"consensus\|(.+?):(\d+)-(\d+)\|(.+?):(\d+)-(\d+)", enhance_ref_list[0])
    if len(region_p) != 1:
        RuntimeError("illegal enhance ref")
    region_p = list(region_p[0])
    # print "region_p={}".format(region_p)
    assert len(region_p) == 6
    region_a = list(re.findall(r"refA\|(.+?):(\d+)-(\d+)\|(.+?):(\d+)-(\d+)", enhance_ref_list[2])[0])
    assert len(region_a) == 6
    region_b = list(re.findall(r"refB\|(.+?):(\d+)-(\d+)\|(.+?):(\d+)-(\d+)", enhance_ref_list[4])[0])
    assert len(region_b) == 6
    rpl = int(region_p[1]) - int(region_p[2])
    rpr = int(region_p[4]) - int(region_p[5])
    ral = int(region_a[1]) - int(region_a[2])
    rar = int(region_a[4]) - int(region_a[5])
    rbl = int(region_b[1]) - int(region_b[2])
    rbr = int(region_b[4]) - int(region_b[5])
    # overlap_len = abs(rpl) + 2 + abs(rpr) - len(enhance_ref_list[1])
    left_extand_num = edge_num - abs(rpl) - 1
    right_extand_num = edge_num - abs(rpr) - 1
    # print "overlap_len={}".format(overlap_len)
    # handle ref_p
    if left_extand_num > 0:
        if rpl > 0:
            enhance_ref_list[1] = "{0}{1}".format(reverse_complement(ref_seq_generator.send([region_p[0],
                                                                                             int(region_p[1]) + 1,
                                                                                             left_extand_num])),
                                                  enhance_ref_list[1])
            region_p[1] = int(region_p[1]) + left_extand_num
        else:
            enhance_ref_list[1] = "{0}{1}".format(ref_seq_generator.send([region_p[0],
                                                                          int(region_p[1]) - left_extand_num,
                                                                          left_extand_num]),
                                                  enhance_ref_list[1])
            region_p[1] = int(region_p[1]) - left_extand_num
        if region_p[1] < 0:
            return -1
        enhance_ref_list[0] = re.sub(r":(\d+)-", ":{}-".format(region_p[1]), enhance_ref_list[0], 1)

    if right_extand_num > 0:
        if rpr > 0:
            enhance_ref_list[1] = "{0}{1}" \
                                  "".format(enhance_ref_list[1],
                                            reverse_complement(ref_seq_generator.send([region_p[3],
                                                                                       int(region_p[
                                                                                               5]) - right_extand_num,
                                                                                       right_extand_num])))
            region_p[5] = int(region_p[5]) - right_extand_num
        else:
            enhance_ref_list[1] = "{0}{1}".format(enhance_ref_list[1],
                                                  ref_seq_generator.send([region_p[3],
                                                                          int(region_p[5]) + 1,
                                                                          right_extand_num]))
            region_p[5] = int(region_p[5]) + right_extand_num
        if region_p[5] < 0:
            return -1
        enhance_ref_list[0] = re.sub(r"-(\d+)$", "-{}".format(region_p[5]), enhance_ref_list[0], 1)

    # handle ref_A
    if left_extand_num > 0:
        if ral > 0:
            enhance_ref_list[3] = "{0}{1}".format(reverse_complement(ref_seq_generator.send([region_a[0],
                                                                                             int(region_a[1]) + 1,
                                                                                             left_extand_num])),
                                                  enhance_ref_list[3])
            region_a[1] = int(region_a[1]) + left_extand_num
        else:
            enhance_ref_list[3] = "{0}{1}".format(ref_seq_generator.send([region_a[0],
                                                                          int(region_a[1]) - left_extand_num,
                                                                          left_extand_num]),
                                                  enhance_ref_list[3])
            region_a[1] = int(region_a[1]) - left_extand_num
        if region_a[1] < 0:
            return -1
        enhance_ref_list[2] = re.sub(r":(\d+)-", ":{}-".format(region_a[1]), enhance_ref_list[2], 1)

    if right_extand_num > 0:
        if rar > 0:
            enhance_ref_list[3] = "{0}{1}" \
                                  "".format(enhance_ref_list[3],
                                            reverse_complement(ref_seq_generator.send([region_a[3],
                                                                                       int(region_a[
                                                                                               5]) - right_extand_num,
                                                                                       right_extand_num])))
            region_a[5] = int(region_a[5]) - right_extand_num
        else:
            enhance_ref_list[3] = "{0}{1}".format(enhance_ref_list[3],
                                                  ref_seq_generator.send([region_a[3],
                                                                          int(region_a[5]) + 1,
                                                                          right_extand_num]))
            region_a[5] = int(region_a[5]) + right_extand_num
        if region_a[5] < 0:
            return -1
        enhance_ref_list[2] = re.sub(r"-(\d+)$", "-{}".format(region_a[5]), enhance_ref_list[2], 1)

    # handle ref_B
    if left_extand_num > 0:
        if rbl > 0:
            enhance_ref_list[5] = "{0}{1}".format(reverse_complement(ref_seq_generator.send([region_b[0],
                                                                                             int(region_b[1]) + 1,
                                                                                             left_extand_num])),
                                                  enhance_ref_list[5])
            region_b[1] = int(region_b[1]) + left_extand_num
        else:
            enhance_ref_list[5] = "{0}{1}".format(ref_seq_generator.send([region_b[0],
                                                                          int(region_b[1]) - left_extand_num,
                                                                          left_extand_num]),
                                                  enhance_ref_list[5])
            region_b[1] = int(region_b[1]) - left_extand_num
        if region_b[1] < 0:
            return -1
        enhance_ref_list[4] = re.sub(r":(\d+)-", ":{}-".format(region_b[1]), enhance_ref_list[4], 1)
    if right_extand_num > 0:
        if rbr > 0:
            enhance_ref_list[5] = "{0}{1}" \
                                  "".format(enhance_ref_list[5],
                                            reverse_complement(ref_seq_generator.send([region_b[3],
                                                                                       int(region_b[
                                                                                               5]) - right_extand_num,
                                                                                       right_extand_num])))
            region_b[5] = int(region_b[5]) - right_extand_num
        else:
            enhance_ref_list[5] = "{0}{1}".format(enhance_ref_list[5],
                                                  ref_seq_generator.send([region_b[3],
                                                                          int(region_b[5]) + 1,
                                                                          right_extand_num]))
            region_b[5] = int(region_b[5]) + right_extand_num
        if region_b[5] < 0:
            return -1
        enhance_ref_list[4] = re.sub(r"-(\d+)$", "-{}".format(region_b[5]), enhance_ref_list[4], 1)
    return 0


def sc3procedure(bwa, file_sv2a, file_ref, file_sv_in,
                 bam_file_sample, sample_fq1, sample_fq2, file_gsnp,
                 r1, r2, output, workdir,
                 ifq_dir):
    """

    :param file_ghdbsnp: germline heterozygous SNPs which is in snp database in vcf format
    :param file_sv2a:
    :param file_ref:
    :param file_sv_in:
    :param enhance_ref_out:
    :param bam_file_sample:
    :param bam_file_control:
    :return:

    """
    prefix = os.path.splitext(os.path.basename(file_sv_in))[0]
    ENHANCE_REF_FILE_ = "{0}_{1}".format(prefix, ENHANCE_REF_FILE)
    OVERLAPPED_SEQ_FILE_ = "{0}_{1}".format(prefix, OVERLAPPED_SEQ_FILE)
    ENHANCE_ALIGNMENT_FILE_ = "{0}_{1}".format(prefix, ENHANCE_ALIGNMENT_FILE)
    FILTERED_MAX_MAPQ_ALIGNMENT_FILE_ = "{0}_{1}".format(prefix, FILTERED_MAX_MAPQ_ALIGNMENT_FILE)
    SORTED_FILTERED_MAX_MAPQ_ALIGNMENT_FILE_ = "{0}_{1}".format(prefix, SORTED_FILTERED_MAX_MAPQ_ALIGNMENT_FILE)
    FILTERED_MAX_MAPQ_ALIGNMENT_PILEUP_FILE_ = "{0}_{1}".format(prefix, FILTERED_MAX_MAPQ_ALIGNMENT_PILEUP_FILE)
    UNUSUAL_SV_FILE_ = "{0}_{1}".format(prefix, UNUSUAL_SV_FILE)
    output = "{0}_{1}".format(output, prefix)

    org_dir = os.getcwd()
    os.chdir(workdir)
    log_format = "%(asctime)s %(levelname)-8s %(process)d [%(funcName)s]%(message)s(%(filename)s line%(lineno)d)"
    # logging.basicConfig(filename=os.path.splitext(os.path.basename(sys.argv[5]))[0] + "_sccaller3.log",
    #                     level=logging.DEBUG,
    #                     format=log_format,
    #                     filemode="w")
    logging.basicConfig(filename=prefix + "_sccaller3.log",
                        level=logging.DEBUG, format=log_format, filemode="w")
    logging.debug("building fq index...")
    r1 = int(r1)
    r2 = int(r2)

    def init_fq_source(sample_fq1, sample_fq2, ifq_dir):  # , control_fq1, control_fq2):
        def fq_index_name(fq_name, suffix, ifq_dir):
            """
            get the ifq and iifq file name
            :type ifq_dir: str
            :param fq_name:
            :param suffix:
            :return: /aa/bb/cc.dd -> ifq_dir/cc.suffix
            """
            file_name = os.path.splitext(os.path.basename(fq_name))[0] + suffix
            if not ifq_dir.endswith("/"):
                ifq_dir = ifq_dir + "/"
            return ifq_dir + file_name
            # return os.path.splitext(fq_name)[0]

        handle_fq_index(sample_fq1,
                        fq_index_name(sample_fq1, ".ifq", ifq_dir),
                        fq_index_name(sample_fq1, ".iifq", ifq_dir))
        handle_fq_index(sample_fq2,
                        fq_index_name(sample_fq2, ".ifq", ifq_dir),
                        fq_index_name(sample_fq2, ".iifq", ifq_dir))
        # handle_fq_index(control_fq1, fq_index_name(control_fq1, ".ifq"), fq_index_name(control_fq1, ".iifq"))
        # handle_fq_index(control_fq2, fq_index_name(control_fq2, ".ifq"), fq_index_name(control_fq2, ".iifq"))
        sample_fq1_source = generator_fq(sample_fq1, fq_index_name(sample_fq1, ".ifq", ifq_dir))
        sample_fq2_source = generator_fq(sample_fq2, fq_index_name(sample_fq2, ".ifq", ifq_dir))

        # control_fq1_source = generator_fq(control_fq1, fq_index_name(control_fq1, ".ifq"))
        # control_fq2_source = generator_fq(control_fq2, fq_index_name(control_fq2, ".ifq"))
        sample_fq1_source.next()
        sample_fq2_source.next()
        # control_fq1_source.next()
        # control_fq2_source.next()
        # return [sample_fq1_source, sample_fq2_source, control_fq1_source, control_fq2_source]
        return [sample_fq1_source, sample_fq2_source]

    def build_enhance_ref(sv_line_list, sv2alignments_dict, fp_unusual,
                          ref_seq_generator, enhance_ref_list, enhance_ref_file_name):
        """

        :param sv_line_list:
        :param sv2alignments_dict:
        :param fp_unusual:
        :param ref_seq_generator:
        :param enhance_ref_list:
        :param enhance_ref_file_name:
        :return:
            -1 can not find original alignments
            -2 cigar can not match, data might be wrong
            -3 there are IDNPX in cigar_p, can not handle
            -4 there are more than 2 SH in cigar_p, can not handle
            -5 there are more than 1 M= in cigar_p, can not handle
            -6 there are IDNPX in cigar_a, can not handle
            -7 there are more than 2 SH in cigar_a, can not handle
            -8 there are more than 1 M= in cigar_a, can not handle
            -9 cigar_a and cigar_p have different total number
            -10 the two cigar can not match with the other
            -11 the alignment might be wrong
        """

        [vcf_key, vcf_key2] = tools.get_vcf_key(sv_line_list)
        key_num = 1
        if vcf_key not in sv2alignments_dict:
            key_num = 2
            if vcf_key2 not in sv2alignments_dict:
                fp_unusual.write("{}\n".format("\t".join(sv_line_list)))
                logging.debug(
                    "{0} not found in sv2alignments_dict\n vcf data = {1}".format(vcf_key, sv_line_list))
                return -1

        # key is in the dict. Get the original alignment info. Parse the alignment into enhanced reference.
        sv2alignments_info = sv2alignments_dict[vcf_key if key_num == 1 else vcf_key2]
        # {contig name_chr1_pos1_chr2_pos2:[contig name, chr1, pos1, flag1, seq1, cigar1, orgcigar1, vcf_pos1,
        #                                        0        1     2      3     4      5        6          7
        #                                   chr2, pos2, flag2, seq2, cigar2, orgcigar2, vcf_pos2
        #                                    8      9    10     11     12       13         14
        #                                   chr_p, pos_p, flag_p, seq_p, cigar_p, orgcigar_p]}
        #                                     15     16      17     18     19        20
        ret = parse_alignments2enhance_ref(sv2alignments_info, ref_seq_generator, enhance_ref_list)
        if ret < 0:
            # if NEED_LOG:
            #     logging.debug("ret = {0} vcf={1}".format(ret, sv_line_list))
            fp_unusual.write("{}\n".format("\t".join(sv_line_list)))
            return ret
        if NEED_LOG:
            logging.debug("enhance ref before extend={}".format(enhance_ref_list))
        copy_enhance_ref_list = copy.copy(enhance_ref_list)
        ret = extend_edge(ref_seq_generator, copy_enhance_ref_list, EDGE_LEN)
        if ret == 0:
            for i in xrange(len(copy_enhance_ref_list)):
                enhance_ref_list[i] = copy_enhance_ref_list[i]
        # got an enhanced reference
        with open(enhance_ref_file_name, "w") as fp_out:
            fp_out.write("{}\n".format("\n".join(enhance_ref_list)))
        if NEED_LOG:
            logging.debug(
                "\ncontig name=[{0}]\nchr1=[{1}]\npos1=[{2}]\nflag1=[{3}]\nseq1=[{4}]\ncigar1=[{5}]\norgcigar1=[{6}]\n"
                "vcf_pos1=[{7}]\nchr2=[{8}]\npos2=[{9}]\nflag2=[{10}]\nseq2=[{11}]\ncigar2=[{12}]\norgcigar2=[{13}]\n"
                "vcf_pos2=[{14}]\nchr_p=[{15}]\npos_p=[{16}]\nflag_p=[{17}]\nseq_p=[{18}]\ncigar_p=[{19}]\n"
                "orgcigar_p=[{20}]\n".format(sv2alignments_info[0], sv2alignments_info[1],
                                             sv2alignments_info[2],
                                             sv2alignments_info[3], sv2alignments_info[4],
                                             sv2alignments_info[5],
                                             sv2alignments_info[6], sv2alignments_info[7],
                                             sv2alignments_info[8],
                                             sv2alignments_info[9], sv2alignments_info[10],
                                             sv2alignments_info[11],
                                             sv2alignments_info[12], sv2alignments_info[13],
                                             sv2alignments_info[14],
                                             sv2alignments_info[15], sv2alignments_info[16],
                                             sv2alignments_info[17],
                                             sv2alignments_info[18], sv2alignments_info[19],
                                             sv2alignments_info[20]))
            logging.debug("enhance ref list = {}".format(enhance_ref_list))
        return 0

    def build_read2max_mapQ_alignment_dict(sam_in):
        # type: (str) -> dict[str, list[list, float]]
        """
        build a dict in which each read indicate only one alignment with the maximum modified mapQ
        :rtype: dict[str:[list,float]]
        :type read_alignments_dict: dict[str:list[list]]
        :param read_alignments_dict: read name: [alignment1_list, alignment2_list]
        :return: dict  read name: max mapQ alignment list  str sam head
        """

        def build_read_name2alignments_dict(sam_in):
            """
            build dict to indicate all the alignments of one read in a sam file
            :type sam_in: str
            :param sam_in: sam file name
            :return: [dict[read name: [alignment1_list, alignment2_list]] ,  sam_head_str]
            alignment_list is sam data line splitted by tab
            """
            read_alignments_dict = {}
            sam_head_str = ""
            with open(sam_in, "r") as fp:
                while True:
                    line = fp.readline()
                    if not line:
                        break
                    if line.startswith("@"):
                        sam_head_str = "{0}{1}".format(sam_head_str, line)
                        continue
                    if not line.strip():
                        continue
                    data_list = line.strip().split("\t")

                    # not mapped
                    if int(data_list[1]) & 4 or data_list[2] == "*":
                        continue
                    if data_list[0] not in read_alignments_dict:
                        read_alignments_dict[data_list[0]] = [data_list]
                    else:
                        read_alignments_dict[data_list[0]].append(data_list)
            return [read_alignments_dict, sam_head_str]

        def get_max_mapQ_alignment(alignments):
            # type: (list[list]) -> list
            """
            get the alignment with the max mapQ from the list of alignments
            :param alignments: [alignment1_list, alignment2_list]
            :return: the alignment list with the maximum mapQ
            """
            mapq_list = map(lambda alignment: int(alignment[4]), alignments)
            index_list = filter(lambda x: mapq_list[x] == max(mapq_list), xrange(len(mapq_list)))
            if len(index_list) > 1:
                return None
            return alignments[index_list[0]]

        def get_modified_error_possibility(max_mapQ_alignment, alignments):
            # type: (list, list[list]) -> float
            """
            modify the possibility with weight algorithm if the sum of all the alignment possibilities is bigger than 1
            :rtype: float
            :param max_mapQ_alignment:
            :param alignments:
            :return:read name : [max mapQ alignment, ]
            """
            p1 = sum(
                map(lambda alignment: 1 - tools.get_error_probability(int(alignment[4])), alignments))  # type: float
            if p1 <= 1:
                return tools.get_error_probability(int(max_mapQ_alignment[4]))  # type: float
            else:
                return 1 - (1 - tools.get_error_probability(int(max_mapQ_alignment[4]))) / p1  # type: float

        read2alignments_dict, sam_head_str = build_read_name2alignments_dict(sam_in)
        read_name2max_mapQ_alignment_dict = {}  # type:dict[str, list[list,float]]
        for read_name in read2alignments_dict:
            alignments_list = read2alignments_dict[read_name]
            max_mapQ_alignment_list = get_max_mapQ_alignment(alignments_list)
            if not max_mapQ_alignment_list:
                continue  # can not get the maximum mapQ alignment, go to next read
            modified_error_p = get_modified_error_possibility(max_mapQ_alignment_list, alignments_list)
            max_mapQ_alignment_list[4] = str(int(round(-10 * math.log10(modified_error_p))))
            # read_max_mapQ_alignment_dict[key] = [max_mapQ_alignment_list, modified_error_p]
            read_name2max_mapQ_alignment_dict[read_name] = max_mapQ_alignment_list
        return [read_name2max_mapQ_alignment_dict, sam_head_str]

    def filter_max_mapq_alignment(alignment, chr1, pos1, chr2, pos2, bp1, bp2):
        """
        Make sure the alignment go through snp and break point
        One snp could be in one of the two positions in enhanced reference.
        :param alignment:
        :param chr1:
        :param pos1:
        :param chr2:
        :param pos2:
        :param bp1:
        :param bp2:
        :return:
        """
        if position_in_alignment(chr1, pos1, alignment):
            if position_in_alignment(chr1, bp1, alignment) and position_in_alignment(chr1, bp2, alignment):
                return True
            else:
                return False
        if position_in_alignment(chr2, pos2, alignment):
            if position_in_alignment(chr2, bp1, alignment) and position_in_alignment(chr2, bp2, alignment):
                return True
            else:
                return False
        return False

    logging.debug("initiating fq source...")
    # print "initiating fq source..."
    [sample_fq1_source, sample_fq2_source] = init_fq_source(sample_fq1, sample_fq2, ifq_dir)
    logging.debug("loading ghdbsnp...")
    # print "loading ghdbsnp..."
    ghdbsnp_list = tools.load_ghdbsnp(file_gsnp)  # type: list[list[str, str, str, str]]  #[[chr, pos, ref, alt], ]

    logging.debug("loading reference...")
    ref_seq_generator = generator_seq_from_ref(file_ref)
    ref_seq_generator.next()
    logging.debug("loading sv2alignment file...")
    # print "loading sv2alignment file..."
    with open(file_sv2a, "r") as fp:
        sv2alignments_list = filter(lambda x: len(x) == 21, [i.strip("\n").split("\t") for i in fp.readlines()])

    # {contig name_chr1_pos1_chr2_pos2:[contig name, chr1, pos1, flag1, seq1, cigar1, orgcigar1, vcf_pos1,
    #                                        0        1     2      3     4      5        6          7
    #                                   chr2, pos2, flag2, seq2, cigar2, orgcigar2, vcf_pos2
    #                                    8      9    10     11     12       13         14
    #                                   chr_p, pos_p, flag_p, seq_p, cigar_p, orgcigar_p]}
    #                                     15     16      17     18     19        20
    sv2alignments_dict = dict(
        zip(["{0}_{1}_{2}_{3}_{4}".format(i[0], i[1], i[7], i[8], i[14]) for i in sv2alignments_list],
            sv2alignments_list))
    del sv2alignments_list

    logging.debug("handling sv...")
    # print "handling sv..."
    with open(file_sv_in, "r") as fp:
        sv_lines = fp.readlines()
    # vcf_head = "".join(filter(lambda x: x.startswith("#"), sv_lines))
    sv_lines = filter(lambda x: not x.startswith("#"), sv_lines)  # type: list[str]
    sv_list = [i.strip().split("\t") for i in sv_lines]  # type: list[list]
    logging.debug("sv_list len = {}".format(len(sv_list)))
    enhance_ref_list = []  # [name_p, seq_p, namea, seq_a, nameb, seq_b]
    icounter = 0
    icounter_1 = 0  # can not build enhance ref
    icounter_1_1 = 0  # can't find original alignments to build enhance ref
    icounter_1_2 = 0  # the cigar doesn't match, the data might be wrong
    icounter_1_3 = 0  # there are IDNPX in cigar_p, can not handle
    icounter_1_4 = 0  # there are more than 2 SH in cigar_p, can not handle
    icounter_1_5 = 0  # there are more than 1 M= in cigar_p, can not handle
    icounter_1_6 = 0  # there are IDNPX in cigar_a, can not handle
    icounter_1_7 = 0  # there are more than 2 SH in cigar_a, can not handle
    icounter_1_8 = 0  # there are more than 1 M= in cigar_a, can not handle
    icounter_1_9 = 0  # cigar_a and cigar_p have different total number
    icounter_1_10 = 0  # the two cigar can not match with the other
    icounter_1_11 = 0  # the alignment might be wrong
    icounter_2 = 0  # no related ghdbsnp
    icounter_3 = 0  # no overlapped alignments with enhanced ref
    icounter_4 = 0
    icounter_5 = 0
    icounter_6 = 0
    icounter_7 = 0
    icounter_8 = 0
    icounter_9 = 0
    with open(UNUSUAL_SV_FILE_, "w") as fp_unusual, open(output, "w") as fp_out:
        for index in xrange(len(sv_list)):
            icounter += 1
            # if icounter % 100 == 0 and icounter > 0:
            logging.debug("handling {0} ({2:.2%})sv {1}".format(icounter, sv_list[index][4], float(icounter)/len(sv_list)))
            # print "handling {0} sv {1}".format(icounter, sv_list[index][4])
            # print "step1 build enhanced ref"
            # step 1: build enhanced ref
            ret = build_enhance_ref(sv_list[index], sv2alignments_dict, fp_unusual,
                                    ref_seq_generator, enhance_ref_list, ENHANCE_REF_FILE_)
            # print "build_enhance_ref=[{}]".format(ret)
            if ret != 0:
                icounter_1 += 1
                if ret == -1:
                    icounter_1_1 += 1  # can't find original alignments to build enhance ref
                    sv_list[index][6] = "case1.1"
                elif ret == -2:
                    icounter_1_2 += 1  # the cigar doesn't match, the data might be wrong
                    sv_list[index][6] = "case1.2"
                elif ret == -3:
                    icounter_1_3 += 1  # there are IDNPX in cigar_p, can not handle
                    sv_list[index][6] = "case1.3"
                elif ret == -4:
                    icounter_1_4 += 1  # there are more than 2 SH in cigar_p, can not handle
                    sv_list[index][6] = "case1.4"
                elif ret == -5:
                    icounter_1_5 += 1  # there are more than 1 M= in cigar_p, can not handle
                    sv_list[index][6] = "case1.5"
                elif ret == -6:
                    icounter_1_6 += 1  # there are IDNPX in cigar_a, can not handle
                    sv_list[index][6] = "case1.6"
                elif ret == -7:
                    icounter_1_7 += 1  # there are more than 2 SH in cigar_a, can not handle
                    sv_list[index][6] = "case1.7"
                elif ret == -8:
                    icounter_1_8 += 1  # there are more than 1 M= in cigar_a, can not handle
                    sv_list[index][6] = "case1.8"
                elif ret == -9:
                    icounter_1_9 += 1  # cigar_a and cigar_p have different total number
                    sv_list[index][6] = "case1.9"
                elif ret == -10:
                    icounter_1_10 += 1  # the two cigar can not match with the other
                    sv_list[index][6] = "case1.10"
                elif ret == -11:
                    icounter_1_11 += 1  # the alignment might be wrong
                    sv_list[index][6] = "case1.11"
                else:
                    RuntimeError("unexpected ret({}) of build_enhance_ref".format(ret))
                # print sv_list[index][6]
                fp_out.write("{}\n".format("\t".join(sv_list[index])))
                if NEED_LOG:
                    logging.debug(sv_list[index][6])
                continue  # go to next sv
            # print "step2 begin get_related_ghdbsnp_dict"
            # step2: check related ghdbsnp number in mut_dict(chr pos ref alt : [H, h, chr1, pos1, chr2, pos2, bpos1, bpos2])
            mut_dict = tools.get_related_ghdbsnp_dict(ghdbsnp_list,
                                                      enhance_ref_list)  # type: dict[str,list[str, str, str, int, str, int, int, int]]
            if len(mut_dict) == 0:
                icounter_2 += 1
                sv_list[index][6] = "case2"
                fp_out.write("{}\n".format("\t".join(sv_list[index])))
                if NEED_LOG:
                    logging.debug("no related ghdbsnp case2")
                # print "case2"
                continue  # go to next sv
            # print "{} related GHDBSNPs".format(len(mut_dict))
            if NEED_LOG:
                logging.debug("related ghdbsnp num = [{}]".format(len(mut_dict)))
            # print "step3 get overlaped read file ---> overlaped.fq"
            # step3: get overlaped read file ---> overlaped.fq
            if get_overlaped_seq(bam_file_sample, tools.get_region_list_from_enhance_ref(enhance_ref_list),
                                 sample_fq1_source, sample_fq2_source, OVERLAPPED_SEQ_FILE_) != 0:
                icounter_3 += 1
                sv_list[index][6] = "case3"
                fp_out.write("{}\n".format("\t".join(sv_list[index])))
                if NEED_LOG:
                    logging.debug("no overlapped sequence, go to next sv {}".format(sv_list[index][6]))
                # print "case3"
                continue  # no overlapped sequence, go to next sv
            # print "step4 enhanced alignment"
            # step4: enhanced alignment
            shell_cmd = "{0} index {1}".format(bwa, ENHANCE_REF_FILE_)
            bwa_process = Popen([shell_cmd], shell=True, stdout=PIPE, stderr=PIPE)
            bwa_process.wait()
            shell_cmd = "{0} mem -M {1} {2} > {3}".format(bwa, ENHANCE_REF_FILE_,
                                                          OVERLAPPED_SEQ_FILE_, ENHANCE_ALIGNMENT_FILE_)
            bwa_process = Popen([shell_cmd], shell=True, stdout=PIPE, stderr=PIPE)
            bwa_process.wait()
            # print "begin step5 build read2max_mapQ_alignment_dict"
            # step5: build read2max_mapQ_alignment_dict
            # read name -> max mapQ alignment list
            read2max_mapQ_alignment_dict, sam_head_str = build_read2max_mapQ_alignment_dict(ENHANCE_ALIGNMENT_FILE_)
            sv_has_result = False
            tmp_counter = 0
            for mut in mut_dict:
                tmp_counter += 1
                H, h, chr1, pos1, chr2, pos2, bpos1, bpos2 = mut_dict[mut]
                assert type(H) == type(h) == type(chr1) == type(chr2) == str
                assert type(pos1) == type(pos2) == type(bpos1) == type(bpos2) == int
                if NEED_LOG:
                    logging.debug("handling mut {0} : {1}".format(mut, mut_dict[mut]))
                    logging.debug("len max mapq alignment = {}".format(len(read2max_mapQ_alignment_dict)))
                # print "handling snp{}".format(tmp_counter)
                # print "begin step6 filter the max mapQ alignment lists"
                # step 6: filter the max mapQ alignment lists
                # the enhanced alignment should go through both the break point and the snp
                filtered_max_mapQ_alignments_list = [value for key, value in read2max_mapQ_alignment_dict.items()
                                                     if filter_max_mapq_alignment(value, chr1, pos1,
                                                                                  chr2, pos2, bpos1, bpos2
                                                                                  )]  # type: list[list,float]
                if len(filtered_max_mapQ_alignments_list) == 0:
                    sv_list[index][6] = "case4"
                    fp_out.write("{}\n".format("\t".join(sv_list[index])))
                    if NEED_LOG:
                        logging.debug("no alignments go across both snp and break point. {} go to next snp".format(
                            sv_list[index][6]))
                    icounter_4 += 1
                    # print "case4"
                    continue  # go to next snp
                if NEED_LOG:
                    logging.debug(
                        "filtered max mapq alignments number = {}".format(len(filtered_max_mapQ_alignments_list)))
                # print "after step6 got {} max mapQ alignments".format(len(filtered_max_mapQ_alignments_list))
                with open(FILTERED_MAX_MAPQ_ALIGNMENT_FILE_, "w") as fp_filtered_max_mapq:
                    fp_filtered_max_mapq.write(sam_head_str)
                    fp_filtered_max_mapq.write("\n".join(["\t".join(i) for i in filtered_max_mapQ_alignments_list]))

                # step 7: pileup the filtered alignments
                # print "step7 pileup the filtered alignments"
                shell_cmd = "samtools sort -O sam {0}>{1}".format(FILTERED_MAX_MAPQ_ALIGNMENT_FILE_,
                                                                  SORTED_FILTERED_MAX_MAPQ_ALIGNMENT_FILE_)
                samtools_process = Popen([shell_cmd], shell=True)
                samtools_process.wait()
                shell_cmd = "samtools mpileup {0} -s > {1}".format(SORTED_FILTERED_MAX_MAPQ_ALIGNMENT_FILE_,
                                                                   FILTERED_MAX_MAPQ_ALIGNMENT_PILEUP_FILE_)
                samtools_process = Popen([shell_cmd], shell=True, stderr=PIPE)
                samtools_process.wait()
                # print "step8 parse the pileup to p_value"
                # step8: parse the pileup to p_value
                ret_list = []
                ret = parse_pileup(FILTERED_MAX_MAPQ_ALIGNMENT_PILEUP_FILE_, chr1, pos1, chr2, pos2, h, H,
                                   filtered_max_mapQ_alignments_list, r1, r2, ret_list)
                # print "parse_pileup = [{}]".format(ret)
                if ret == 0:
                    sv_has_result = True
                    fp_out.write("{0}\t{1}\t{2}\n".format(sv_lines[index].strip(), mut, "\t".join(ret_list)))
                    if NEED_LOG:
                        logging.debug("got 1 result")
                elif ret == -1:
                    icounter_5 += 1  # enhance alignments don't cover the two position in pileup file
                    sv_list[index][6] = "case5"
                    fp_out.write("{}\n".format("\t".join(sv_list[index])))
                    if NEED_LOG:
                        logging.debug("enhance alignments don't cover the two position in pileup file. {}".format(
                            sv_list[index][6]))
                elif ret == -3:
                    icounter_6 += 1  # max number of reads support H is less then the threshold r1
                    sv_list[index][6] = "case6"
                    fp_out.write("{}\n".format("\t".join(sv_list[index])))
                    if NEED_LOG:
                        logging.debug(
                            "max number of reads support H is less then the threshold r1. {}".format(sv_list[index][6]))
                elif ret == -4:
                    icounter_7 += 1  # max number of reads support h is less then the threshold r2
                    sv_list[index][6] = "case7"
                    fp_out.write("{}\n".format("\t".join(sv_list[index])))
                    if NEED_LOG:
                        logging.debug(
                            "max number of reads support h is less then the threshold r2. {}".format(sv_list[index][6]))
                elif ret == -5:
                    icounter_8 += 1  # can not calculate false positive
                    sv_list[index][6] = "case8"
                    fp_out.write("{}\n".format("\t".join(sv_list[index])))
                    if NEED_LOG:
                        logging.debug("can not calculate false positive. {}".format(sv_list[index][6]))
                else:
                    RuntimeError("unexpected ret({}) of parse_pileup".format(ret))
            if sv_has_result:
                icounter_9 += 1
            # break
    logging.debug("handled {0} sv in total".format(icounter))
    logging.debug("case 1 can not build enhance ref ----{0} SVs {1:.2%}"
                  "".format(icounter_1, icounter_1 / float(icounter)))
    logging.debug("case 1.1 can't find original alignments to build enhance ref ----{0} SVs {1:.2%}"
                  "".format(icounter_1_1, icounter_1_1 / float(icounter)))
    logging.debug("case 1.2 the cigar doesn't match, the data might be wrong ----{0} SVs {1:.2%}"
                  "".format(icounter_1_2, icounter_1_2 / float(icounter)))
    logging.debug("case 1.3 there are IDNPX in cigar_p, can not handle ----{0} SVs {1:.2%}"
                  "".format(icounter_1_3, icounter_1_3 / float(icounter)))
    logging.debug("case 1.4 there are more than 2 SH in cigar_p, can not handle ----{0} SVs {1:.2%}"
                  "".format(icounter_1_4, icounter_1_4 / float(icounter)))
    logging.debug("case 1.5 there are more than 1 M= in cigar_p, can not handle ----{0} SVs {1:.2%}"
                  "".format(icounter_1_5, icounter_1_5 / float(icounter)))
    logging.debug("case 1.6 there are IDNPX in cigar_a, can not handle ----{0} SVs {1:.2%}"
                  "".format(icounter_1_6, icounter_1_6 / float(icounter)))
    logging.debug("case 1.7 there are more than 2 SH in cigar_a, can not handle ----{0} SVs {1:.2%}"
                  "".format(icounter_1_7, icounter_1_7 / float(icounter)))
    logging.debug("case 1.8 there are more than 1 M= in cigar_a, can not handle ----{0} SVs {1:.2%}"
                  "".format(icounter_1_8, icounter_1_8 / float(icounter)))
    logging.debug("case 1.9 cigar_a and cigar_p have different total number ----{0} SVs {1:.2%}"
                  "".format(icounter_1_9, icounter_1_9 / float(icounter)))
    logging.debug("case 1.10 the two cigar can not match with the other ----{0} SVs {1:.2%}"
                  "".format(icounter_1_10, icounter_1_10 / float(icounter)))
    logging.debug("case 1.11 the alignment might be wrong ----{0} SVs {1:.2%}"
                  "".format(icounter_1_11, icounter_1_11 / float(icounter)))
    logging.debug("case 2 there is no related ghdbsnp ----{0} SVs {1:.2%}"
                  "".format(icounter_2, icounter_2 / float(icounter)))
    logging.debug("case 3 there is no overlapped alignments with enhanced ref ----{0} SVs {1:.2%}"
                  "".format(icounter_3, icounter_3 / float(icounter)))
    logging.debug("case 4 no alignments go across both snp and break point ----{0} sv/snp".format(icounter_4))
    logging.debug(
        "case 5 enhance alignments don't cover the two position in pileup file ----{0} sv/snp".format(icounter_5))
    logging.debug(
        "case 6 max number of reads support H is less then the threshold r1 ----{0} sv/snp".format(icounter_6))
    logging.debug(
        "case 7 max number of reads support h is less then the threshold r2 ----{0} sv/snp".format(icounter_7))
    logging.debug("case 8 can not calculate false positive ----{0} sv/snp".format(icounter_8))
    logging.debug("there are {0} SVs {1:.2%} in reault"
                  "".format(icounter_9, icounter_9 / float(icounter)))
    logging.debug("all done")
    os.chdir(org_dir)


def sc3procedure_multiple(bwa, file_sv2a, file_ref, file_sv_in, cpu_num, job_size,
                          bam_file_sample, sample_fq1, sample_fq2, file_gsnp,
                          r1, r2, output, workdir, ifq_dir):
    """

    :param file_ghdbsnp: germline heterozygous SNPs which is in snp database in vcf format
    :param file_sv2a:
    :param file_ref:
    :param file_sv_in:
    :param enhance_ref_out:
    :param bam_file_sample:
    :param bam_file_control:
    :return:

    """

    # if not workdir.endswith("/"):
    #     workdir += "/"
    suffix = os.path.splitext(os.path.split(output)[-1])[0]
    # print "suffix={}".format(suffix)
    ph = Popen(["mkdir -p {}".format(workdir)], shell=True, stdout=PIPE)
    ph.wait()

    cmd_str = "split -l {0} {1} {2}".format(job_size, file_sv_in, os.path.join(workdir, "sc3splitvcf_"))
    ph = Popen([cmd_str], shell=True, stdout=PIPE)
    ph.wait()
    process_pool = multiprocessing.Pool(processes=int(cpu_num))
    for splited_vcf in os.listdir(workdir):
        # print splited_vcf
        if not os.path.isdir(os.path.join(workdir, splited_vcf)) and splited_vcf.startswith("sc3splitvcf_"):
            process_pool.apply_async(sc3procedure, (bwa, file_sv2a, file_ref,
                                                    splited_vcf, bam_file_sample, sample_fq1,
                                                    sample_fq2, file_gsnp, r1,
                                                    r2, suffix, workdir,
                                                    ifq_dir))
    process_pool.close()
    process_pool.join()
    ph = Popen(["cat {0} > {1}".format(os.path.join(workdir, suffix + "_*"), output)], shell=True, stdout=PIPE)
    ph.wait()

    ph = Popen(["rm {}".format(os.path.join(workdir, "*.sam"))], shell=True, stdout=PIPE)
    ph.wait()
    ph = Popen(["rm {}".format(os.path.join(workdir, "*.fasta"))], shell=True, stdout=PIPE)
    ph.wait()
    ph = Popen(["rm {}".format(os.path.join(workdir, "*.fasta.amb"))], shell=True, stdout=PIPE)
    ph.wait()
    ph = Popen(["rm {}".format(os.path.join(workdir, "*.fasta.ann"))], shell=True, stdout=PIPE)
    ph.wait()
    ph = Popen(["rm {}".format(os.path.join(workdir, "*.fasta.bwt"))], shell=True, stdout=PIPE)
    ph.wait()
    ph = Popen(["rm {}".format(os.path.join(workdir, "*.fasta.pac"))], shell=True, stdout=PIPE)
    ph.wait()
    ph = Popen(["rm {}".format(os.path.join(workdir, "*.fasta.sa"))], shell=True, stdout=PIPE)
    ph.wait()
    ph = Popen(["rm {}".format(os.path.join(workdir, "*.pileup"))], shell=True, stdout=PIPE)
    ph.wait()
    ph = Popen(["rm {}".format(os.path.join(workdir, "*.log"))], shell=True, stdout=PIPE)
    ph.wait()
    ph = Popen(["rm {}".format(os.path.join(workdir, "*_unusual"))], shell=True, stdout=PIPE)
    ph.wait()
    ph = Popen(["rm {}".format(os.path.join(workdir, "*.fq"))], shell=True, stdout=PIPE)
    ph.wait()


def data_generator_samtools(file_sample, region=None):
    def read_line_from_process(process_handle, ch, spliter, columns, data_out):
        # type: (Popen, str, str, int, list) -> int
        """
        Read data from the process's stdout, filter out the comments, split by spliter, and check the format
        :param process_handle: handle of the process
        :param ch: comment character
        :param spliter:
        :param columns: input. Greater than 0: column number of the data. Others, Unlimited
        :param data_out: output data (should clear buffer before using)
        :return: 0, Success. Others, processes have been terminated, and stdout can't read the data.
        """

        def read_line_from_file(from_file, ch, spliter, columns, data_out):
            # type: (file, str, str, int, list) -> int
            """
            Read data from file, filter out comments, and split by spliter, check format
            :param from_file: file pointer
            :param ch: comment character
            :param spliter:
            :param columns: input. Greater than 0: column number of the data. Others, Unlimited
            :param data_out: output data (should clear buffer before using)
            :return: 0, Success. Others, end of file.
            eg:
            data = []
            with open("test") as file:
            ret = read_data_from_file(file,"#","\t",8,data)
            """
            while 1:
                buf = from_file.readline()
                if len(buf) == 0:
                    return -1
                if buf.startswith(ch):
                    continue
                buf = buf.strip("\n")
                buf2 = buf.split(spliter)
                if columns > 0:
                    if len(buf2) != columns:
                        continue
                    else:
                        break
                else:
                    break
            data_out.extend(buf2)
            return 0

        while 1:
            buf = []
            if read_line_from_file(process_handle.stdout, ch, spliter, columns, buf) != 0:
                if not process_handle.poll() is None:  # samtools has been terminated
                    return -1
                else:
                    continue
            else:
                del data_out[:]
                data_out.extend(buf)
                return 0

    if region != None:
        cmd_str = "samtools view {0} {1}".format(file_sample, region)
    else:
        cmd_str = "samtools view {0}".format(file_sample)
    process_samtools = Popen([cmd_str], shell=True, stdout=PIPE)
    data_buf = []
    while 1:
        if read_line_from_process(process_samtools, "@", "\t", 0, data_buf) != 0:
            break
        else:
            yield data_buf
    while 1:
        yield None


def generator_fq(fq_file, ifq_file):
    def line2list(line):
        data_list = line.strip("\n").split("\t")
        return [data_list[0], int(data_list[1])]

    ifq_file_name_list = ifq_file.split(".")
    ifq_file_name_list[-1] = "iifq"
    iifq_file = ".".join(ifq_file_name_list)
    del ifq_file_name_list

    if os.path.exists(iifq_file):
        with open(iifq_file, "r") as fp_iifq:
            iifq_list = [line2list(i) for i in fp_iifq.readlines()]  # type: list[list[str,int]]
    else:
        iifq_list = build_second_index(ifq_file)

    iiifq_list = zip([i[0] for i in iifq_list[::100]],
                     [i * 100 for i in list(xrange(int(math.ceil(len(iifq_list) / 100.0))))])
    if iifq_list[-1][0] != iiifq_list[-1][0]:
        iiifq_list.append((iifq_list[-1][0], len(iifq_list) - 1))
    # logging.debug("len(iiifq_list)={0} iiifq_list[-1]={1}".format(len(iiifq_list),iiifq_list[-1]))

    i4fq_list = zip([i[0] for i in iiifq_list[::100]],
                    [i * 100 for i in list(xrange(int(math.ceil(len(iiifq_list) / 100.0))))])
    if iiifq_list[-1][0] != i4fq_list[-1][0]:
        i4fq_list.append((iiifq_list[-1][0], len(iiifq_list) - 1))
    # logging.debug("len(i4fq_list)={0} i4fq_list[-1]={1}".format(len(i4fq_list),i4fq_list[-1]))
    name = yield ""
    i4fq_list_len = len(i4fq_list)
    with open(ifq_file, "r") as fp_ifq, open(fq_file, "r") as fp_fq:
        while True:
            # logging.debug("len i4fq_list = {0} searching [{1}]".format(len(i4fq_list), name))
            # logging.debug("i4fq_list:{}".format(i4fq_list))
            iiifq_offset = -1
            for index in xrange(i4fq_list_len):
                # if name_bigger(name, i4fq_list[index][0]):
                if operator.gt(name, i4fq_list[index][0]):
                    continue
                if name == i4fq_list[index][0]:
                    iiifq_offset = i4fq_list[index][1]
                    # iiifq_name = i4fq_list[index][0]
                    # logging.debug("name[{0}]==iiifq_name[{1}]".format(name,iiifq_name))
                    break
                iiifq_offset = i4fq_list[index - 1][1]
                # iiifq_name = i4fq_list[index - 1][0]
                # logging.debug("iiifq_name[{0}]<name[{1}]<[{2}]".format(iiifq_name, name, i4fq_list[index][0]))
                break
            if iiifq_offset == -1:
                logging.error("Didn't find {0} in i4fq_list".format(name))
                print "Didn't find {0} in i4fq_list".format(name)
                name = yield ""
                continue
                # raise RuntimeError("Didn't find {0} in fq file [{1}]".format(name, fq_file))
            # logging.debug("iiifq_offset=[{0}] iiifq_name=[{1}]".format(iiifq_offset, iiifq_name))
            iifq_offset = -1
            # for index in xrange(len(iiifq_list)):
            for index in xrange(101):
                # if name_bigger(name, iiifq_list[iiifq_offset + index][0]):
                if operator.gt(name, iiifq_list[iiifq_offset + index][0]):
                    # logging.debug("index = [{0}] name [{1}] is bigger than [{2}]".format(index, name, iiifq_list[index][0]))
                    continue
                # else:
                #     logging.debug(
                #         "index = [{0}] name [{1}] is smaller than [{2}]".format(iiifq_offset + index, name, iiifq_list[iiifq_offset + index][0]))
                if name == iiifq_list[iiifq_offset + index][0]:
                    iifq_offset = iiifq_list[iiifq_offset + index][1]
                    # iifq_name = iiifq_list[iiifq_offset + index][0]
                    # logging.debug("got index=[{0}] iifq_offset=[{1}] iifq_name[{2}]==name[{3}] ".format(index, iifq_offset, iifq_name,name))
                    break
                iifq_offset = iiifq_list[iiifq_offset + index - 1][1]
                # iifq_name = iiifq_list[iiifq_offset + index - 1][0]
                # logging.debug("2 got index=[{0}] iifq_offset=[{1}] iifq_name[{2}] < name[{3}] < [{4}]".format(index, iifq_offset,
                #                                                                                          iifq_name, name,
                #                                                                                          iiifq_list[iiifq_offset + index][0]))
                break
            if iifq_offset == -1:
                logging.error("Didn't find {0} in iiifq_list".format(name))
                print "Didn't find {0} in iiifq_list".format(name)
                name = yield ""
                continue
                # raise RuntimeError("Didn't find {0} in fq file [{1}]".format(name, fq_file))
            # logging.debug("iifq_offset = [{0}] iifq_name = [{1}]".format(iifq_offset, iifq_name))
            # iifq_list ---> ifq_offset
            ifq_offset = -1
            for index in xrange(101):
                # logging.debug("candidate iifq_name={}".format(iifq_list[iifq_offset + index][0]))
                # if name_bigger(name, iifq_list[iifq_offset + index][0]):
                if operator.gt(name, iifq_list[iifq_offset + index][0]):
                    # logging.debug("continue")
                    continue
                if name == iifq_list[iifq_offset + index][0]:
                    ifq_offset = iifq_list[iifq_offset + index][1]
                    sample_name = iifq_list[iifq_offset + index][0]
                    # logging.debug("sample_name[{0}] == name[{1}]".format(sample_name, name))
                    break
                ifq_offset = iifq_list[iifq_offset + index - 1][1]
                sample_name = iifq_list[iifq_offset + index - 1][0]
                # logging.debug("sample_name=[{0}] < name[{1}] < [{2}]".format(sample_name, name, iifq_list[iifq_offset + index][0]))
                break
            if ifq_offset == -1:
                logging.error("Didn't find {0} in iifq_list".format(name))
                print "Didn't find {0} in iifq_list".format(name)
                name = yield ""
                continue
                # raise RuntimeError("Didn't find {0} in fq file [{1}]".format(name, fq_file))

            # logging.debug("ifq_offset = {0} sample_name = {1} name = [{2}]".format(ifq_offset, sample_name, name))
            # ifq_file ---> fq_offset
            fq_offset = -1
            fp_ifq.seek(ifq_offset, 0)
            sample_line = get_sample_line_from_ifq(fp_ifq)
            sample_list = sample_line.strip("\n").split("\t")
            if sample_name == name:
                fq_offset = int(sample_list[1])
            else:
                for i in xrange(INDEXSTEP + 1):
                    length_str = fp_ifq.read(1)
                    if not length_str:
                        break
                    length = struct.unpack("<B", length_str)[0]
                    tlv_str = fp_ifq.read(length)
                    fp_ifq.seek(1, 1)
                    ifq_line = rebuild_line_with_sample(sample_line, tlv_str)
                    ifq_list = ifq_line.strip("\n").split("\t")
                    if name == ifq_list[0]:
                        fq_offset = int(ifq_list[1])
                        break
            if fq_offset == -1:
                logging.debug("Didn't find {0} in fq file [{1}]. Return \"\".".format(name, ifq_file))
                print "Didn't find {0} in fq file [{1}]. Return \"\".".format(name, ifq_file)
                name = yield ""
                continue
                # raise RuntimeError("Didn't find {0} in fq file [{1}]".format(name, fq_file))

            fp_fq.seek(fq_offset)

            ret = fp_fq.readline()
            curr_seq_name = ret[1:].strip().split(" ")[0]
            if curr_seq_name.endswith("/1") or curr_seq_name.endswith("/2"):
                curr_seq_name = curr_seq_name[:-2]
            ret_list = [[ret, fp_fq.readline(), fp_fq.readline(), fp_fq.readline()]]
            # ret = "{0}{1}{2}{3}".format(ret, fp_fq.readline(), fp_fq.readline(), fp_fq.readline())
            while True:
                first_line = fp_fq.readline()
                if not first_line:
                    break
                first_line = first_line.strip()
                if not first_line:
                    continue
                seq_name = first_line[1:].split(" ")[0]
                if seq_name.endswith("/1") or seq_name.endswith("/2"):
                    seq_name = seq_name[:-2]
                if seq_name != curr_seq_name:
                    break
                ret_list.append(["{}\n".format(first_line), fp_fq.readline(), fp_fq.readline(), fp_fq.readline()])
                # ret = "{0}{1}\n{2}{3}{4}".format(ret, first_line, fp_fq.readline(), fp_fq.readline(), fp_fq.readline())
            ret_length_list = map(lambda x: len(x[1]), ret_list)
            if ret_length_list.count(max(ret_length_list)) != 1:
                logging.debug("Didn't find {0} in fq file [{1}]. Return \"\".".format(name, ifq_file))
                print "Didn't find {0} in fq file [{1}]. Return \"\".".format(name, ifq_file)
                name = yield ""
                continue
            name = yield "".join(ret_list[ret_length_list.index(max(ret_length_list))])


def rebuild_line_with_sample(sample, tlv_str):
    tlv_list = []
    unpack_tlv(tlv_str, tlv_list)
    ret_line = ""
    handled_num = 0
    for oper in tlv_list:
        if oper[0] == TAG_DIFF:
            ret_line += oper[1]
            handled_num += len(oper[1])
        elif oper[0] == TAG_DIFFHEX:
            tmp_str = str(int(oper[1], 16))
            ret_line += tmp_str
            handled_num += len(tmp_str)
        elif oper[0] == TAG_SAME_NUM:
            same_num = struct.unpack("<B", oper[1])[0]
            ret_line += sample[handled_num:handled_num + same_num]
            handled_num += same_num
        elif oper[0] == TAG_SAMPLE:
            ret_line += oper[1]
            handled_num += len(oper[1])
    return ret_line


def get_sample_line_from_ifq(fp_ifq):
    """

    :type fp_ifq: file
    """
    curr = fp_ifq.tell()
    length_str = fp_ifq.read(1)
    if not length_str:
        return None
    length = struct.unpack("<B", length_str)[0]
    tlv_str = fp_ifq.read(length)
    fp_ifq.seek(1, 1)
    tlv_list = []
    unpack_tlv(tlv_str, tlv_list)
    assert len(tlv_list) == 1
    assert tlv_list[0][0] == TAG_SAMPLE
    return tlv_list[0][1]


def build_second_index(index_file):
    iCounter = 0
    assert type(index_file) == str
    should_break = False
    ret = []
    with open(index_file, "r") as fp_index:
        index_list = []
        while True:
            offset = fp_index.tell()
            sample_line = get_sample_line_from_ifq(fp_index)
            if not sample_line:
                if index_list:
                    ret.append([index_list[0], pre_offset])
                break
            sample_list = sample_line.split("\t")
            ret.append([sample_list[0], offset])
            for i in xrange(INDEXSTEP - 1):
                offset = fp_index.tell()
                length_str = fp_index.read(1)
                if not length_str:
                    should_break = True
                    if i > 0:
                        ret.append([index_list[0], pre_offset])
                    break
                pre_offset = offset
                length = struct.unpack("<B", length_str)[0]
                tlv_str = fp_index.read(length)
                fp_index.seek(1, 1)
                index_list = rebuild_line_with_sample(sample_line, tlv_str).split("\t")
            if should_break:
                break
            iCounter += INDEXSTEP
            if iCounter % 100000 == 0:
                logging.debug("handled {}".format(iCounter))
    return ret


def handle_fq_index(fq_file, index_file, index2_file):
    index_file_ok = False  # type: bool
    index2_file_ok = False  # type: bool
    if os.path.exists(index_file) and os.path.getmtime(index_file) > os.path.getmtime(fq_file):
        index_file_ok = True
        logging.debug("{} is fine.".format(index_file))

    if os.path.exists(index2_file) and os.path.getmtime(index2_file) > os.path.getmtime(fq_file):
        index2_file_ok = True
        logging.debug("{} is fine.".format(index2_file))

    if index_file_ok and index2_file_ok:
        return

    if not index_file_ok:
        logging.debug("now building {}".format(index_file))
        print "now building {}".format(index_file)
        # build first level index
        tmp_file = index_file + ".tmp"
        logging.debug("tmp_file = [{}]".format(tmp_file))
        print "tmp_file = [{}]".format(tmp_file)
        with open(tmp_file, "w") as fp_index, open(fq_file, "r") as fp_fq:
            curr_name = ""
            iread_times = 0
            while True:
                offset = fp_fq.tell()
                name_line = fp_fq.readline()
                if not name_line:
                    break
                name_line = name_line.strip()
                if not name_line:
                    continue
                name = name_line[1:].split(" ")[0]
                if name.endswith("/1") or name.endswith("/2"):
                    name = name[:-2]
                fp_fq.readline()
                fp_fq.readline()
                fp_fq.readline()
                iread_times += 1
                if curr_name != name and curr_name != "":
                    if iread_times == 2:
                        fp_index.write("{0}\t{1}\n".format(last_data[0], last_data[1]))
                    else:
                        logging.critical(
                            "Your fastq file might be wrong. Sequence [{}] are repeated. It will not be indexed."
                            "".format(last_data[0]))
                        print "Your fastq file might be wrong. Sequence [{}] are repeated. It will not be indexed." \
                              "".format(last_data[0])
                    curr_name = name
                    iread_times = 1
                    last_data = [name, offset]
                elif curr_name != name and curr_name == "":
                    curr_name = name
                    last_data = [name, offset]
            if iread_times == 1:
                fp_index.write("{0}\t{1}\n".format(last_data[0], last_data[1]))
            else:
                logging.critical(
                    "Your fastq file might be wrong. Sequence [{}] are repeated. It will not be indexed."
                    "".format(last_data[0]))
                print "Your fastq file might be wrong. Sequence [{}] are repeated. It will not be indexed." \
                      "".format(last_data[0])

        # process_sort = Popen(["cp ./{0} ./ifq_org".format(tmp_file)], shell=True)
        # process_sort.wait()
        index_path = os.path.split(index_file)[0] if os.path.split(index_file)[0] else "."

        process_sort = Popen(["sort -T {0} {1} > {2}".format(index_path, tmp_file, index_file)], shell=True)
        process_sort.wait()
        process_rm = Popen(["rm {}".format(index_file + ".tmp")], shell=True)
        process_rm.wait()
        process_mv = Popen(["mv {0} {1}".format(index_file, index_file + ".tmp")], shell=True)
        process_mv.wait()

        compress_ifq(index_file + ".tmp", index_file)
        process_rm = Popen(["rm {}".format(index_file + ".tmp")], shell=True)
        process_rm.wait()
    logging.debug("build index1 over")
    print "build index1 over"
    # build second level index
    if not index2_file_ok:
        logging.debug("now building {}".format(index2_file))
        print "now building {}".format(index2_file)
        second_index_list = build_second_index(index_file)
        with open(index2_file, "w") as fp_index2:
            fp_index2.write("\n".join(["{0}\t{1}".format(i[0], i[1]) for i in second_index_list]))
    logging.debug("index1 and index2 are ok now.")
    print "index1 and index2 are ok now."


def reverse_seq(string):
    return string[::-1]


complement_dict = {"A": "T", "a": "t",
                   "T": "A", "t": "a",
                   "C": "G", "c": "g",
                   "G": "C", "g": "c",
                   "Y": "R", "y": "r",
                   "R": "Y", "r": "y",
                   "H": "D", "h": "d",
                   "D": "H", "d": "h",
                   "N": "N", "n": "n",
                   "K": "M", "k": "m",
                   "M": "K", "m": "k",
                   "S": "S", "s": "s",
                   "W": "W", "w": "w",
                   "B": "V", "b": "v",
                   "V": "B", "v": "b",
                   "-": "-"}


def complement_seq(string_in):
    my_list = list(string_in)
    for index in xrange(len(my_list)):
        if my_list[index] in complement_dict:
            my_list[index] = complement_dict[my_list[index]]
        else:
            logging.debug("there is unexpected charactor in seq=[{}]".format(string_in))
            raise RuntimeError("wrong character in seq")
    return "".join(my_list)


def reverse_complement(seq):
    return reverse_seq(complement_seq(seq))


def pack_tlv(tag_num, value_str):
    """
    pack tlv data
    :param tag_num:
    :param value_str:
    :return:
    """
    if tag_num >= 4:
        raise RuntimeError("tag number should be less than 8")
    tag_num = tag_num << 6
    tag_len_list = []
    if len(value_str) <= 31:
        tag_num = tag_num + len(value_str)
        tag_len_list.append(tag_num)
    else:
        tag_num = tag_num + 63
        tag_len_list.append(tag_num)
        # for i in xrange(int(math.ceil((len(value_str) - 31) / 127.0)) - 1):
        for i in xrange(((len(value_str) - 31) - (len(value_str) - 31) % 127) / 127):
            tag_len_list.append(255)
        tag_len_list.append((len(value_str) - 31) % 127)
    fmt_list = map(lambda x: struct.pack("<B", x), tag_len_list)

    return "{0}{1}".format("".join(fmt_list), value_str)

    # return str(tag_str+struct.pack("<B", len(value_str))+value_str)


def unpack_tlv(tlv_str, tlv_ret):
    """

    :type tlv_ret: list
    """
    del tlv_ret[:]
    index = 0
    while True:
        tag = struct.unpack("<B", tlv_str[index])[0] >> 6
        value_len = 0
        len_len = 0
        # parse len
        while True:
            curr_len = struct.unpack("<B", tlv_str[index + len_len])[0]
            if len_len == 0:
                # len_curr = struct.unpack("<B", tlv_str[index])[0] & 31
                # flag = struct.unpack("<B", tlv_str[0])[0] & 32
                len_curr = curr_len & 31
                flag = curr_len & 32
            else:
                # len_curr = struct.unpack("<B", tlv_str[index])[0] & 127
                # flag = struct.unpack("<B", tlv_str[0])[0] & 128
                len_curr = curr_len & 127
                flag = curr_len & 128
            len_len += 1
            value_len += len_curr
            if not flag:
                break

        value_str = tlv_str[index + len_len:index + len_len + value_len]
        tlv_ret.append([tag, value_str])
        index += (len_len + value_len)
        if index >= len(tlv_str):
            break


def build_new_line(sample, sample_len, curr_line):
    """
    compress the curr_line which is not the sample
    ps: length of curr_line should be less than 256
    :type curr_line: str
    """

    def get_same_num(sample, sample_len, curr_line, last_index):
        same_num = 0
        for index in xrange(len(curr_line) - last_index - 1):
            if last_index + index + 2 > sample_len:
                break
            if sample[last_index + 1 + index] == curr_line[last_index + 1 + index]:
                same_num += 1
            else:
                break
        return [same_num, same_num + last_index]

    def get_different_num(sample, sample_len, curr_line, last_index):
        different_num = 0
        for index in xrange(len(curr_line) - last_index - 1):
            if last_index + index + 2 > sample_len:
                different_num += 1
                continue
            if sample[last_index + 1 + index] != curr_line[last_index + 1 + index]:
                different_num += 1
            else:
                break
        return [different_num, last_index + different_num]

    last_index = -1
    curr_line_len = len(curr_line)
    # get oper list
    oper_list = []
    while True:
        curr_pos = last_index + 1
        [same_num, last_index] = get_same_num(sample, sample_len, curr_line, last_index)
        oper_list.append([0, same_num, last_index])
        [different_num, last_index] = get_different_num(sample, sample_len, curr_line, last_index)
        oper_list.append([1, different_num, last_index])
        if last_index == curr_line_len - 1:
            break
    # print oper_list

    # modify oper list
    filter = [True for i in xrange(len(oper_list))]
    for index in xrange(len(oper_list) - 1):
        if oper_list[index][0] == 0 and oper_list[index][1] <= 2:
            oper_list[index + 1][1] += oper_list[index][1]
            filter[index] = False
            if index > 0:
                oper_list[index + 1][1] += oper_list[index - 1][1]
                filter[index - 1] = False
    oper_list = list(compress(oper_list, filter))
    # print oper_list

    # operation
    new_line = ""
    for operation in oper_list:
        value = curr_line[(operation[2] + 1 - operation[1]):(operation[2] + 1)]
        if operation[0] == 0:  # 1 same as sample
            new_line += pack_tlv(TAG_SAME_NUM, struct.pack("<B", len(value)))
        else:
            if value.isdigit() and value[0] != "0":
                new_line += pack_tlv(TAG_DIFFHEX, hex(int(value))[2:])  # 3 different number(in hex format)
            else:
                new_line += pack_tlv(TAG_DIFF, value)  # 2 different string
    new_line = struct.pack("<B", len(new_line)) + new_line + struct.pack("<B", len(new_line))
    return new_line


def compress_ifq(org_ifq, c_ifq):
    """
    compress the ifq file to tlv version.
    0 sample line
    1 number of bytes same as the sample(in hex format, should copy from sample directly)
    2 different string
    3 different number(in hex format)
    :param org_ifq:
    :param c_ifq:
    :return:
    """
    curr_line_index = 0
    with open(org_ifq, "r") as fp_org, open(c_ifq, "w") as fp_c:
        while True:
            curr_line = fp_org.readline()
            curr_pos = fp_org.tell()
            next_line = fp_org.readline()
            fp_org.seek(curr_pos)

            if not curr_line:
                break
            if len(curr_line) <= 2:
                continue
            curr_line = curr_line.strip("\n")
            if curr_line_index % INDEXSTEP == 0 or not next_line:
                sample_line = curr_line
                sample_len = len(sample_line)
                new_line = pack_tlv(TAG_SAMPLE, sample_line)
                new_line = struct.pack("<B", len(new_line)) + new_line + struct.pack("<B", len(new_line))
            else:
                new_line = build_new_line(sample_line, sample_len, curr_line)
            fp_c.write(new_line)
            curr_line_index += 1


def analyze2(sv_file, snp_file, delta):
    """
    write all the sv which has snp nearby in the log file.
    :param sv_file:
    :param snp_file:
    :param delta:
    :return:
    """
    delta = int(delta)

    def filter_snp(snp, chr1, pos1, chr2, pos2, delta):
        chr_snp = snp[0]
        pos_snp = snp[1]
        if chr_snp != chr1 and chr_snp != chr2:
            return False
        if chr_snp == chr1 and pos1 - delta <= pos_snp <= pos1 + delta:
            return True
        if chr_snp == chr2 and pos2 - delta <= pos_snp <= pos2 + delta:
            return True
        return False

    log_format = "%(asctime)s %(levelname)-8s %(process)d [%(funcName)s]%(message)s(%(filename)s line%(lineno)d)"
    logging.basicConfig(filename="analyze2.log", level=logging.DEBUG, format=log_format, filemode="w")
    print "loading snps..."
    with open(snp_file, "r") as fp_snp:
        snp_data = [j.split("\t") for j in
                    filter(lambda x: (not x.startswith("#")) and len(x) > 1, [i.strip() for i in fp_snp.readlines()])]
    snp_list = [[i[0], int(i[1])] for i in snp_data]

    print "loading sv..."
    with open(sv_file, "r") as fp_sv:
        sv_data = [i.strip().split("\t") for i in
                   filter(lambda x: (not x.startswith("#")) and len(x.strip()) > 1, fp_sv.readlines())]
    pos2_list = [re.findall(r"CHR2=(.+?);END=(\d+)", i[7])[0] for i in sv_data]
    # ret = []
    icounter = 0
    print "begin handling sv"
    for i in xrange(len(sv_data)):
        icounter += 1
        if icounter % 10 == 0:
            print "handling {0} / {1} sv".format(icounter, len(sv_data))
        # ret.append([sv_data[i][0], int(sv_data[i][1]), pos2_list[i][0], int(pos2_list[i][1])])
        chr1 = sv_data[i][0]
        pos1 = int(sv_data[i][1])
        chr2 = pos2_list[i][0]
        pos2 = int(pos2_list[i][1])
        curr_snp_list = filter(lambda x: filter_snp(x, chr1, pos1, chr2, pos2, delta), snp_list)
        if len(curr_snp_list) > 0:
            logging.debug("sv={}".format(sv_data[i]))
            logging.debug("chr1={0} pos1={1} chr2={2} pos2={3}".format(chr1, pos1, chr2, pos2))
            for j in xrange(len(curr_snp_list)):
                logging.debug("snp{0}={1}".format(j + 1, curr_snp_list[j]))

    print "all done"
    logging.debug("all done")


def select_sc3_log(sc3_log, chr1, pos1, chr2, pos2):
    cmd_str = "cat {0} | grep {1} | grep consensus | grep build_enhance_ref | grep {2} | grep {3}: | grep {4}:" \
              "".format(sc3_log, pos1, pos2, chr1, chr2)
    selected_line = Popen([cmd_str], shell=True, stdout=PIPE).communicate()[0]
    selected_line = re.findall(r"consensus\|(.+?)',", selected_line)

    selected_line = filter(lambda x: str(pos1) in x and str(pos2) in x, selected_line)
    selected_line = selected_line[0]
    # print selected_line
    return selected_line


def analyze3(analyze2_log, sc3_log):
    def parse_selected_sc3_log(selected_line):
        return [re.split("[:-]", i) for i in selected_line.split("|")]

    def parse_overlap_rate(region1, region2):
        tmp1 = int(region1[1])
        tmp2 = int(region1[2])
        pos11 = min(tmp1, tmp2)
        pos12 = max(tmp1, tmp2)
        tmp3 = int(region2[1])
        tmp4 = int(region2[2])
        pos21 = min(tmp3, tmp4)
        pos22 = max(tmp3, tmp4)
        pos_max = max([tmp1, tmp2, tmp3, tmp4])
        pos_min = min([tmp1, tmp2, tmp3, tmp4])
        return (min(pos12, pos22) - max(pos11, pos21)) / float(pos_max - pos_min)

    def snp_in_region(snp_chr, snp_pos, region):
        """

        :param snp_chr: str
        :param snp_pos: int
        :param region: ['1', '1045607', '1045719']
        :return:
        """
        if snp_chr != region[0]:
            return False
        if snp_pos < min(int(region[1]), int(region[2])) or snp_pos > max(int(region[1]), int(region[2])):
            return False
        return True

    log_format = "%(asctime)s %(levelname)-8s %(process)d [%(funcName)s]%(message)s(%(filename)s line%(lineno)d)"
    logging.basicConfig(filename="analyze3.log", level=logging.DEBUG, format=log_format, filemode="w")
    icounter = 0
    with open(analyze2_log, "r") as fp_analyze2_log:
        while True:
            analyze2_log_line = fp_analyze2_log.readline()
            if not analyze2_log_line:
                break
            if "[analyze2]sv=[" in analyze2_log_line:
                continue
            if "[analyze2]chr1=" in analyze2_log_line:
                chr1, pos1, chr2, pos2 = list(
                    re.findall(r"chr1=(.+?) pos1=(\d+) chr2=(.+?) pos2=(\d+)\(", analyze2_log_line)[0])
                snp_list = []  # [(chr, pos)]
                while True:
                    analyze2_log_line = fp_analyze2_log.readline()
                    if "[analyze2]snp" in analyze2_log_line:
                        analyze2_block = re.findall(r"snp\d+=\['([A-Za-z0-9_.]+)', (\d+)]\(", analyze2_log_line)
                        if len(analyze2_block) == 0:
                            print analyze2_log_line
                        snp_list.append(analyze2_block[0])
                    else:
                        break

                selected_sc3_line = select_sc3_log(sc3_log, chr1, pos1, chr2, pos2)
                # print selected_sc3_line
                region_list = parse_selected_sc3_log(selected_sc3_line)
                # print "snp_list = {}".format(snp_list)
                assert len(region_list) == 2
                overlape_rate = parse_overlap_rate(region_list[0], region_list[1])
                # print "region_list={}".format(region_list)
                icounter += 1
                # if icounter % 10 == 0:
                print "handled {} sv ".format(icounter)
                for snp in snp_list:
                    for region in region_list:
                        if snp_in_region(snp[0], int(snp[1]), region):
                            logging.debug(
                                "this line has snp.selected_sc3 = [{0}], snp=[{1}] overlap rate={2}"
                                "".format(selected_sc3_line, snp, overlape_rate))
    print "all done"


def analyze4(sccaller_log, output):
    log_line1 = ''
    log_line2 = ''
    with open(sccaller_log, "r") as fp, open(output, "w") as fp_out:
        while True:
            log_line1 = log_line2
            log_line2 = fp.readline()
            if not log_line2:
                break
            if "case2" in log_line2:
                fp_out.write(log_line1)

    print "all done"



task_queue = Queue.Queue()
result_queue = Queue.Queue(maxsize=150)
msgc_queue = Queue.Queue()
msgm_queue = Queue.Queue()

def main():
    if len(sys.argv) == 1:
        print """Usage: 
    {0} intersection <fq1> <fq2> <fq1out> <fq2out>
    {0} get_reads <sam_file> <sv_file> <ssake_file> <contig2reads_file> <out_support_sam>
    {0} analyze <file_in> <threshold> <output>
    {0} sc3procedure <bwa> <file_sv2a> <file_ref> <file_sv_in> <bam_file_sample> <sample_fq1> <sample_fq2> <file_gsnp> <r1> <r2> <output> <workdir> <ifq_dir>
    {0} reverse_complement <seq>
    {0} handle_fq_index <fq_file> <index_file> <index2_file>
    {0} analyze2 <sv_file> <snp_file> <delta>
    {0} analyze3 <analyze2_log> <sc3_log>
    {0} analyze4 <sccaller_log> <output>
    {0} sc3procedure_multiple <bwa> <file_sv2a> <file_ref> <file_sv_in> <cpu_num> <job_size>
                          <bam_file_sample> <sample_fq1> <sample_fq2> <file_gsnp>
                          <r1> <r2> <output> <workdir> <ifq_dir>
    """.format("python " + sys.argv[0])
        exit(0)
    # logging.debug("got cmd {}".format(sys.argv[1]))
    if sys.argv[1] == "intersection":
        intersection(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == "get_reads":
        get_reads(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "analyze":
        analyze(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "sc3procedure":

        sc3procedure(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
                     sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9],
                     sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13],
                     sys.argv[14])
    elif sys.argv[1] == "reverse_complement":
        reverse_complement(sys.argv[2])
    elif sys.argv[1] == "handle_fq_index":
        handle_fq_index(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "analyze2":
        analyze2(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "analyze3":
        analyze3(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "analyze4":
        analyze4(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "sc3procedure_multiple":
        sc3procedure_multiple(sys.argv[2], sys.argv[3], sys.argv[4],
                              sys.argv[5], sys.argv[6], sys.argv[7],
                              sys.argv[8], sys.argv[9], sys.argv[10],
                              sys.argv[11], sys.argv[12], sys.argv[13],
                              sys.argv[14], sys.argv[15], sys.argv[16])
    elif sys.argv[1] == "test":
        test()
    else:
        print "wrong argument"


def test_index_fq():
    icounter = 0
    # logging.debug("handling index.")
    # handle_fq_index("./data/hs-j258-h24-r1.fq", "./data/hs-j258-h24-r1.ifq", "./data/hs-j258-h24-r1.iifq")
    #
    # start = time.time()
    # gfq = generator_fq("./data/hs-j258-h24-r1.fq", "./data/hs-j258-h24-r1.ifq")
    # gfq.next()
    # logging.debug("build generator take {}s".format(str(time.time() - start)))
    # start = time.time()

    # with open("./data/hs-j258-h24-r1.fq", "r") as fp, open("result2", "w") as fp_out:
    #     while True:
    #         if icounter % 10000 == 0:
    #             logging.debug("handled {} fq lines".format(icounter))
    #         line = fp.readline()
    #         if not line:
    #             break
    #         icounter += 4
    #         name = line.split(" ")[0][1:]
    #         fp_out.write(gfq.send(name))
    #         fp.readline()
    #         fp.readline()
    #         fp.readline()
    # logging.debug("duration = {0}s handled {1} search requests".format(str(time.time() - start), icounter/4))

    handle_fq_index("./data/hs-j258-h24-r1.fq", "./data/hs-j258-h24-r1.ifq", "./data/hs-j258-h24-r1.iifq")
    logging.debug("building generator...")
    start = time.time()
    gfq = generator_fq("./data/hs-j258-h24-r1.fq", "./data/hs-j258-h24-r1.ifq")
    gfq.next()
    logging.debug("build generator take {}s".format(str(time.time() - start)))
    start = time.time()
    # gfq.send("HWI-ST0844:0326:H5WYCBCXX:1:1101:10182:84759")
    with open("./data/testdata", "r") as fp, open("test_result", "w") as fp_out:
        while True:
            if icounter % 10000 == 0:
                logging.debug("handled {} fq lines".format(icounter))
            line = fp.readline()
            if not line:
                break
            icounter += 4
            name = line.split(" ")[0][1:]
            # if name == "HWI-ST0844:0326:H5WYCBCXX:1:1101:10032:99000":
            #     name = name
            fp_out.write(gfq.send(name))
            fp.readline()
            fp.readline()
            fp.readline()
    logging.debug("duration = {0}s handled {1} search requests".format(str(time.time() - start), icounter / 4))


def test():
    # A2
    org_seq = "TGTGTGAGGCATTTGTGAGTTCCATTAGCACTGATAGGACTTTTAAATGCAGAAGCAGCAGTGAAGGGAATGGCATTTCCCAGCTGTAGGGCACGGGCCTTCTTTGAGATAAAAACCCACCCTGGCATCGATGACATGGGAACTATTCTGTGCTTACACCGGTGTGCCAGGGTTCAAACCGGTCGCGGTTTTATCGATCGACCGTCTGGTTCGATTCTACTAGTACGAAACCCTTGATACACTTCATCAAAGTGGATAGCACCGTCGAAATCGAGGGCTAGCTCCGTCAGGCCGGAATTTCAACCTCAAATGTTCAGCCAGGCTACCCATGTACTCGACCTGTCTTTTCACCCGAGAGCTCTTCAATCAACTGATGATCGCACCTGCATAACTGCTACCAGACCTGCTAAGGGGGAGCCTGGCCCAGCCATCTCTTCTTTGTGGTCACAAGCATGAACGGCCCTGGGACATCCTCTTTCCTCACTGTCCC"
    read_len = 150
    step = 10
    icounter = 0
    with open("oo", "w") as fp:
        while 1:
            if icounter * step + read_len > len(org_seq):
                break
            read_seq = org_seq[icounter * step:icounter * step + read_len]
            fp.write(">test_{0}\n{1}\n".format(icounter * step, read_seq))
            icounter += 1


if __name__ == "__main__":
    main()
