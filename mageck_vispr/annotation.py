__author__ = "Chen-Hao Chen"
__copyright__ = "Copyright 2015, Chen-Hao Chen, Liu lab"
__email__ = "hyalin1127@gmail.com"
__license__ = "MIT"

import os
import sys
import time
import csv
from collections import defaultdict, Counter
from urllib.request import urlopen
import bz2
import gzip
import io
import logging
import pandas as pd
import operator


# TODO rewrite this with Pandas


class Annotator():
    def __init__(self, library ):
        #self.customized_table = annotation_table
        self.sequence_table = library
        self.sequence_dict = defaultdict(list)
        self.sequence_set = set()
        self.seq_match_record = defaultdict(list)
        self.non_gene_match_record = defaultdict(list)
        self.value_frame = pd.DataFrame()
        self.value_frame_column=None
        self.estimated_sgrna_len=None # estimation of sgrna length

    def add_value_frame(self,args):
        if args.bedvalue is not None:
            self.value_frame=pd.read_table(args.bedvalue,index_col=0)
            self.value_frame_column=args.bedvalue_column
            if self.value_frame_column is None:
                raise SyntaxError("need to specify --bedvalue-column option.")
                # exit(1)
            if self.value_frame_column not in self.value_frame:
                raise SyntaxError(""+args.bedvalue_column+" is not in the columns of "+args.bedvalue+".")
                #exit(1)

    def annotate(self,args):
        self.sequence_table_import()
        self.custom_bed_get(args)
        self.write_output()

    def sequence_table_import(self):
        with open(self.sequence_table) as csvfile:
            # [cuiyb]++ to support gRNA library in csv and tab-separated txt format
            if self.sequence_table.upper().endswith('CSV'):
                delimiter = ','
            else:
                delimiter = '\t'
            reader = csv.reader(csvfile, delimiter=delimiter) # sgRNAid seq gene_id
            possible_sg_len={}
            for elements in reader:
                self.sequence_dict[elements[1].upper()] = [elements[0],
                        elements[2].upper()] # 
                self.sequence_set.add(elements[1].upper())
                this_sg_len=len(elements[1])
                if this_sg_len not in possible_sg_len:
                    possible_sg_len[this_sg_len]=0
                possible_sg_len[this_sg_len]=possible_sg_len[this_sg_len]+1
            # estimate the most likely sgrna length
            if len(possible_sg_len)>1: # a mixture of different sgrna length?
                logging.warning('The library file contails a mixture of sgRNAs with different lengths.')
                leninfo=','.join([str(a)+':'+str(b) for (a,b) in possible_sg_len.items()])
                logging.warning('sgRNA length and count in the library is: '+leninfo)
            #self.estimated_sgrna_len=max(possible_sg_len.items(),key=operator.itemgetter(1))[0]
            self.estimated_sgrna_len=[k for k in possible_sg_len.keys()]
            logging.info('Estimated sgRNA length:'+str(self.estimated_sgrna_len))



    def custom_bed_get(self,args):
        #library=args.library
        library=self.sequence_table

        annotation_table=args.annotation_table
        assembly=args.assembly
        sgrna_len=args.sgrna_len
        candidate_file_list=[]
        if sgrna_len=='AUTO':
            sgrna_len=None
        if annotation_table is None: 
            if assembly is None:
                logging.error('Need to specify the --assembly option if annotation table is not provided.')
            sgrna_len_candidate=[]
            if sgrna_len is None: # the sgrna_len is not specified, get from the library
                #sgrna_len_candidate=self.estimated_sgrna_len
                for sg_i in self.estimated_sgrna_len: 
                    if sg_i!=19 and sg_i!=20:
                        logging.warning('Unsupported sgRNA length: '+str(sg_i)+'. These sgRNAs will not be annotated.')
                        #exit(0)
                    else:
                        sgrna_len_candidate+=[sg_i]
            else:
                sgrna_len_candidate=[sgrna_len] # the sgrna_len is specified

            for sg_c in sgrna_len_candidate:
                if args.annotation_table_folder is not None:
                    annotation_table_file = (
                        "sgrna_annotation_{assembly}_exome_{len}bp.txt.bz2"
                    ).format(assembly=assembly,len=sg_c)
                    annotation_table=os.path.join(args.annotation_table_folder,annotation_table_file)
                    logging.info("Using local annotation library: "+annotation_table)
                else:
                    annotation_table = (
                        "https://bitbucket.org/liulab/mageck-vispr/"
                        "downloads/sgrna_annotation_{assembly}_exome_{len}bp.txt.bz2"
                    ).format(assembly=assembly,len=sg_c)
                    logging.info("Downloading files from bitbucket:"+annotation_table)
                candidate_file_list+=[annotation_table]
        else:
            # the annotation table is specified
            logging.info("Using existing annotation library:"+annotation_table)
            candidate_file_list=[annotation_table]

        # self.customized_table=annotation_table


        for candidate_file in candidate_file_list:
            if candidate_file.startswith("http"):
                file = urlopen(candidate_file)
            else:
                file = open(candidate_file, "rb")
            if candidate_file.endswith(".bz2"):
                file = io.BufferedReader(bz2.open(file))
            elif candidate_file.endswith(".gz"):
                file = io.BufferedReader(gzip.open(file))

            for i, line in enumerate(file):
                chr, chrstart, chrend, gene, score, strand, seq = line.decode(
                ).strip().split("\t")
                gene = gene.upper()
                seq = seq.upper()
                try:
                    int(chrstart)
                    int(chrend)
                    float(score)

                except ValueError:
                    raise SyntaxError(
                        "Error parsing line {} in annotation table.".format(i))

                if (seq in self.sequence_set):
                    # sequence_dict[seq]=(sgRNA_id,gene_id)
                    library_sg_id=self.sequence_dict[seq][0]
                    library_gene_id=self.sequence_dict[seq][1]
                    if self.value_frame_column is not None:
                        if library_sg_id in self.value_frame.index:
                            score=self.value_frame[self.value_frame_column][library_sg_id]
                        else:
                            score=0
                        score=str(score)
                    insert = [chr, chrstart, chrend, library_sg_id,
                              score, strand]
                    self.seq_match_record[seq].append(
                        insert + [gene, library_gene_id, seq])

            # end for i, line in enumerate(file):
        # end for candidate_file in candidate_file_list:

    def write_output(self):
        for values in self.seq_match_record.values():
            record = 0
            for i in values:
                if i[6] == i[7]:
                    print("{0}".format("\t".join(i[:6])))
                    record = 1

            if record == 0:
                print("{0}".format("\t".join(i[:6])))
                logging.warning("{0}".format("\t".join(
                    ["Warning: gene not matched"] + [i[k] for k in [3, 8, 7]] +
                    ["|"] + [i[k] for k in [0, 1, 2, 6]])))

        for j in [i for i in self.sequence_dict.keys()
            if i not in self.seq_match_record.keys()]:
                temp = self.sequence_dict[j]
                logging.warning("{0}".format("\t".join(["Warning: sequence not found", temp[0],j, temp[1]])))



