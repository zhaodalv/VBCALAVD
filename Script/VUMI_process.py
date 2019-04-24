#!/usr/bin/env python


from __future__ import print_function
import pysam
import collections
import sys
import ast
import argparse
from virtual_family_library import Ds_filter
from virtual_family_library import xopen



def get_family(pileupcolumn,input_ref,alt_family_cutoff=2.0,f_cutoff=1.0,with_tag=''):
    #print (alt_family_cutoff)
    templen=''
    start=''
    base_strand=''
    quali_family=0
    alt_base=[]
    all_family=collections.defaultdict(int)
    alt_family = collections.defaultdict(list)
    result=[]
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_proper_pair:
            start=min(pileupread.alignment.reference_start,pileupread.alignment.next_reference_start)
            templen=abs(pileupread.alignment.template_length)
            if pileupread.alignment.is_reverse is pileupread.alignment.is_read1:
                base_strand='-'
            else:
                base_strand='+'            
            VUMI=str(start)+","+str(templen)+","+base_strand
            all_family[VUMI]=all_family[VUMI]+1    
            read_base = pileupread.alignment.query_sequence[pileupread.query_position]
            if read_base != input_ref:
                Ds_value = Ds_filter(start,start+templen,pileupcolumn.pos)
                if abs(Ds_value) <3 or Ds_value >=149:
                      continue
                alt_family[VUMI].append(read_base)
    candidate = get_candidate(all_family,alt_family,f_cutoff)
    if candidate:
        for key,value in candidate[1].items():
            if value >= alt_family_cutoff:
                result.append([key,value,candidate[0][key]])
        return result


def get_candidate(all_f,alt_f,f_cutoff):
    #print (f_cutoff)
    f_value=0.0
    qualified_alt_family=collections.defaultdict(int)
    alt_family = collections.defaultdict(int)
    v_singleton=collections.defaultdict(int)
    singleton_ratio={}
    alt_family_number=0
    for key,value in alt_f.items():
        if len(set(value))==1:
            f_value = float(len(value))/all_f[key]
            if f_value >= f_cutoff:
                if len(value)==1:
                    v_singleton[value[0]]= v_singleton[value[0]]+1.0
                else:
                    if len(value)>=2: #min_size is 2                                          
                        qualified_alt_family[value[0]] = qualified_alt_family[value[0]]+1.0
                    #alt_family[value[0]] = alt_family[value[0]]+1.0
    if len(qualified_alt_family)==0:
        return False
    else:
       #for key,value in alt_family.items():
        for key,value in qualified_alt_family.items():
            if key in v_singleton:
                singleton_ratio[key] = v_singleton[key]/value
            else:
                singleton_ratio[key] =0.0
    return [singleton_ratio,qualified_alt_family]


def get_chromose_params(Infile):
    data=""
    input_position =[] 
    input_ref ={} 
    with xopen(Infile,'r') as f:
        for line in f:
            if line.startswith("#"):
                   continue
            data = line.strip().split("\t")
            base_0_Pos = int(data[1])-1
            input_position.append(base_0_Pos)
            input_ref[base_0_Pos]=data[3]
    input_position = sorted(input_position)
    return data[0],input_position,input_ref

        

def main():
    family_info= collections.defaultdict(int)
    input_ref ={} 
    input_position = []
    candidate_info=[]
    ap = argparse.ArgumentParser()
    ap.add_argument("-I","--Infile",help="VCF format file",type=str)
    ap.add_argument("-alt_support","--alt_family_cutoff",help="Minimum alt family support;default is 2.0",nargs="?",type=float,default=2.0)
    ap.add_argument("-fratio","--f_cutoff",help="f value cutoff;Default is 1.0",nargs="?",const=1.0,type=float,default=1.0)
    ap.add_argument("-UMI","--UMI_tag",help="UMI_tag_provided;Default is false",type=ast.literal_eval,nargs="?",const=False,default=False)
    ap.add_argument("bamfilepath",help="bamfile_path",type=str)
    args=ap.parse_args()

    samfile = pysam.AlignmentFile(args.bamfilepath, "rb" )
    chro,input_position,input_ref = get_chromose_params(args.Infile)    

    for pileupcolumn in samfile.pileup(chro,input_position[0],input_position[-1]+1,max_depth=100000,truncate=True,nofilter=True):
       if pileupcolumn.pos not in input_position:
           continue
       else:
           ref_base = input_ref[pileupcolumn.pos]
           candidate_info = get_family(pileupcolumn,ref_base,args.alt_family_cutoff,args.f_cutoff)
           if candidate_info:
               for single_info in candidate_info:
                   out_data = [chro,pileupcolumn.pos+1,".",ref_base,single_info[0],single_info[1],single_info[2]]
                   print (*out_data,sep="\t")
    samfile.close()
    

if __name__=="__main__":
    main()
