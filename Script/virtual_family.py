from __future__ import print_function
import ast    
import pysam
import collections 
import sys
import re
import os
import pandas as pd
import pickle
import numpy as np
import argparse
from virtual_family_library import xopen
from virtual_family_library import Ds_filter

def get_family(samfile,chromosome,position,minquality=30,mq=30,withtag=''):
    templen=''
    start=''
    tag=''
    virtual_family=collections.defaultdict(list)
    virtual_strand =collections.defaultdict(list)
    qname= collections.defaultdict(list)
    
    for pileupcolumn in samfile.pileup(chromosome, position-1, position,max_depth=100000,truncate=True):
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_proper_pair:
                if pileupread.alignment.query_qualities[pileupread.query_position]< minquality:
                    continue
                if pileupread.alignment.mapping_quality <mq:
                    continue
                qname[pileupread.alignment.qname].append(pileupread.alignment.query_sequence[pileupread.query_position])

    for pileupcolumn in samfile.pileup(chromosome, position-1, position,max_depth=100000,truncate=True):
        for pileupread in pileupcolumn.pileups:
            if withtag:
                templen=abs(pileupread.alignment.template_length)
                start=min(pileupread.alignment.reference_start,pileupread.alignment.next_reference_start)
                tag=pileupread.alignment.query_name.split('|')[0]
            else:
                start=min(pileupread.alignment.reference_start,pileupread.alignment.next_reference_start)
                templen=abs(pileupread.alignment.template_length)
            if pileupread.alignment.is_reverse is pileupread.alignment.is_read1:
                base_strand='-'
            else:
                base_strand='+'
            VUMI=str(start) +','+str(templen)
            if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_proper_pair:
                if pileupread.alignment.query_qualities[pileupread.query_position]< minquality:
                    continue
                if pileupread.alignment.mapping_quality <mq:
                    continue
                if len(set(qname[pileupread.alignment.qname]))==1:
                    read_base = pileupread.alignment.query_sequence[pileupread.query_position]
                    virtual_family[VUMI].append(read_base)
                    virtual_strand[VUMI].append(base_strand)
                    del qname[pileupread.alignment.qname]
                else:
                    continue
    return [virtual_family,virtual_strand,position-1]

    
def get_var_info(virtual_info,alt,cutoff=1.0):
    family_number=0
    alt_family_size=[]
    alt_number=0 
    var_raw =0
    depth_raw=0
    f_value=0.0
    f_1_number=0
    f_1_size=[]
    virtual_strand=0
    virtual_strand_list=[]
    virtual_strand_f_1=[]
    template_length_list=[]
    variant_singleton=0
    plus_strand_number=0
    minus_strand_number=0
    virtual_duplex=0
    plus_strand_f1=0
    minus_strand_f1=0
    virtual_duplex_f1=0
    mean_f1_size=0
    std_f1_size=0
    mean_template=0
    std_template=0
    query_position = virtual_info[2]
    for family_tag,baselist in virtual_info[0].items():
        start_position,template_length = [int(x) for x in family_tag.split(",")]
        #template_length = int(family_tag.split(",")[0])
        family_number = family_number+1
        alt_number = baselist.count(alt)
        var_raw = var_raw + alt_number
        family_size=len(baselist)  
        depth_raw =depth_raw+len(baselist)
        f_value = alt_number/float(family_size)

        if set(virtual_info[1][family_tag]).issubset(["+"]):
            virtual_strand=1
            virtual_strand_list.append(1)
        elif set(virtual_info[1][family_tag]).issubset(["-"]):
            virtual_strand=2
            virtual_strand_list.append(2)
        else:
            virtual_strand=3
            virtual_strand_list.append(3)
    
        if f_value >= cutoff:

            Ds_value = Ds_filter(start_position,start_position+template_length,query_position)
            if Ds_value<3 or Ds_value>=149:
                continue

            template_length_list.append(template_length)
            f_1_number =f_1_number+1
            f_1_size.append(family_size)
            virtual_strand_f_1.append(virtual_strand)

    plus_strand_number = virtual_strand_list.count(1)
    minus_strand_number = virtual_strand_list.count(2)
    virtual_duplex = virtual_strand_list.count(3)
    plus_strand_f1=virtual_strand_f_1.count(1)
    minus_strand_f1=virtual_strand_f_1.count(2)
    virtual_duplex_f1=virtual_strand_f_1.count(3)
    mean_f1_size = np.mean(f_1_size)
    std_f1_size = np.std(f_1_size)
    mean_template=np.mean(template_length_list)
    std_template = np.std(template_length_list)

    return [depth_raw,var_raw,family_number,plus_strand_number,minus_strand_number,virtual_duplex,f_1_number,plus_strand_f1,minus_strand_f1,virtual_duplex_f1,mean_f1_size,std_f1_size,mean_template,std_template] 


def main():
    ap = argparse.ArgumentParser()
     
    ap.add_argument("-bq","--basequality",help="mininum base quality;default is 30",nargs="?", const=30, type=int,default=30)
    ap.add_argument("-mq","--mappingquality",help="mininum base quality;default is 30",nargs="?", const=30, type=int,default=30)
    ap.add_argument("-fratio","--familyratio",help="mininum family ratio;default is 1.0",nargs="?",const=1.0,type=float,default=1.0)
    ap.add_argument("-tag","--UMI",help="UMI data;default is false",type=ast.literal_eval,nargs="?",const=False,default=False)
    ap.add_argument("bamfilepath",help="bamfile_path",type=str)
    if "-I" in sys.argv[1:] or "--Infile" in sys.argv[1:]:
        ap.add_argument("-I","--Infile",help="VCF format file",type=str)
        args=ap.parse_args()
    elif "-ih" in sys.argv[1:] or "--ihgvs" in sys.argv[1]:
        ap.add_argument("-ih","--ihgvs",help="file following HGVS standards:7:55259515T>C",type=str)
        args = ap.parse_args()
    else:
        ap.add_argument("-is","--isite", help = "single hgvs site",type=str)
        args = ap.parse_args()

    samfile = pysam.AlignmentFile(args.bamfilepath, "rb" )
    print ("\t".join(["#chr","start","end","ref","alt","depth_raw","var_raw","family_number","+","-","+/-","f1 number","f1+","f1-","f+/-","mean_size","std_size","mean_length","std_length","singleton_ratio","fdr_s_ratio","P_s_ratio"]))

    if 'Infile' in vars(args).keys():
       with xopen(args.Infile,'r') as f:
           for line in f:
               if line.startswith("#"):
                   continue
               data = line.strip().split("\t")
               #print (data)
               family_info = get_family(samfile,data[0],int(data[1]),args.basequality,args.mappingquality,args.UMI)
               var_info = get_var_info(family_info,data[4],args.familyratio)
               var_info.extend(data[6:9])
               all_info = data[0:5]
               all_info.extend(var_info)
               print (*all_info,sep="\t")  
             
    elif "ihgvs" in vars(args).keys():
       pattern = re.compile("[0-9]+|[ATCG]$")
       with xopen(args.ihgvs,'r') as f:
           for line in f:
               in_data =pattern.findall(line.strip())
               family_info = get_family(samfile,data[0],int(data[1]),args.basequality,args.mappingquality,args.UMI)
               var_info = get_var_info(family_info,data[4],args.familyratio)
               all_info = data[0:5].extend(var_info)
               print (*all_info,sep="\t")  
    else:
       pattern = re.compile("[0-9]+|[ATCG]$")
       in_data =pattern.findall(args.isite)
       family_info = get_family(samfile,data[0],int(data[1]),args.basequality,args.mappingquality,args.UMI)
       var_info = get_var_info(family_info,data[4],args.familyratio)
       all_info = data[0:5].extend(var_info)
       print (*all_info,sep="\t")  
    samfile.close()
    
if __name__=='__main__':
    main()
