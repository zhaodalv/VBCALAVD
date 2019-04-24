from __future__ import print_function
import pandas as pd
import scipy.stats as st
import argparse
import numpy as np
from virtual_family_library import get_fdr, get_Pval,FDR_adjust 

#def singleton_ratio(ratio_array):
#    loc,scale = st.norm.fit(ratio_array,cutoff)
#    return (loc,scale)

def calculate_cutoff(ratio_array,cutoff_p):
    row=ratio_array
    param = st.johnsonsu.fit(row)
    arg = param[:-2]
    loc = param[-2]
    scale = param[-1]
    cutoff_value = st.johnsonsu.ppf(1-cutoff_p, loc=loc, scale=scale,*arg)
    r=st.probplot(ratio_array,dist=st.johnsonsu,sparams=param)[-1][-1]
    print ("Fit r value is: ",r,"&& cutoff value is: ",cutoff_value)
    return cutoff_value

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-candidate","--Candidate_file",help = "Result from VUMI preprocessing",type=str)
    ap.add_argument("-mutation","--Mutation_file",help = "Mutation level file from VD",type=str)
    ap.add_argument("-cutoff_P","--cutoff_p_value",help="cutoff probability",type=float)
    ap.add_argument("out_path",help = "out_candidate_file",type=str)
    ap.add_argument("filtered_path",help = "filtered_site info",type=str)
    args = ap.parse_args()
    mutation = pd.read_csv(args.Mutation_file,header=None,sep="\t",low_memory=False)
    candidate =  pd.read_csv(args.Candidate_file,header=None,sep="\t",low_memory=False)
    mutation.columns = ["#chrom","position","end","ref_base","alt","depth","variant"]
    mutation.loc[:,"VAF_raw"] = mutation["variant"]/mutation["depth"]
    candidate.columns = ["#chrom","position","end","ref_base","alt","family_number","singleton_ratio"]
    mutation=mutation.astype({"#chrom":str,"position":int,"ref_base":str,"alt":str})
    candidate = candidate.astype({"#chrom":str,"position":int,"ref_base":str,"alt":str})
    candidate_merge = pd.merge(candidate,mutation[["#chrom","position","ref_base","alt","VAF_raw"]],on=["#chrom","position","ref_base","alt"]).copy()

    param = st.johnsonsu.fit(candidate_merge.loc[candidate_merge["VAF_raw"]>=0.05,"singleton_ratio"])
    r=st.probplot(candidate_merge.loc[candidate_merge["VAF_raw"]>=0.05,"singleton_ratio"],dist=st.johnsonsu,sparams=param)[-1][-1]
    print ("Fit r value is: ",r)
    length = len(candidate_merge)
    candidate_merge.loc[:,"P_val"] = candidate_merge[["singleton_ratio"]].apply(lambda row:get_Pval(row[0],param),axis=1)
    candidate_merge.sort_values(["P_val"],inplace=True,ascending=False)
    candidate_merge.loc[:,"queue_num"] = np.arange(len(candidate_merge),0,-1)
    #candidate_merge.loc[:,"fdr"] = candidate_merge[["P_val","queue_num"]].apply(lambda row:get_fdr(row[0],row[1],length),axis=1)
    candidate_merge.loc[:,"fdr_ad"] = FDR_adjust(candidate_merge["P_val"])
    candidate_merge = candidate_merge[["#chrom","position","position","ref_base","alt","family_number","singleton_ratio","fdr_ad","P_val","VAF_raw"]]
    #FDR_adjust(candidate_merge["P_val"])
    


    #candidate.loc[:,"VAF_raw"] = mutation.loc[mutation.index.isin(candidate.index),"VAF_raw"].copy()
    #loc,scale = singleton_ratio(candidate["singleton_ratio"])
   # cutoff_value = calculate_cutoff(candidate_merge.loc[candidate_merge["VAF_raw"]>=0.05,"singleton_ratio"],args.cutoff_p_value)

   
     #candidate.loc[:,"P_value"] = candidate[["singleton_ratio"]].apply(lambda row:st.norm.sf(row[0],loc,scale),axis=1)


    out_candidate_profile = candidate_merge.loc[candidate_merge["fdr_ad"]>args.cutoff_p_value].copy()
    out_candidate_profile.to_csv(args.out_path,sep="\t",index=False,header=False)
    filtered_candidate_profile = candidate_merge.loc[candidate_merge["fdr_ad"]<=args.cutoff_p_value].copy()
    filtered_candidate_profile.to_csv(args.filtered_path,sep="\t",index=False,header=False)
    
    

if __name__ == '__main__':
    main()
