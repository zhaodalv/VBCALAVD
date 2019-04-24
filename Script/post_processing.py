from __future__ import print_function
import pandas as pd
import argparse
import scipy.stats as st


#this file is used to generate real SNV info from candidate
#1.genmerline anno
#2.minmun template,template length assement and polishing assessment if possible
#3.virtual_duplex and strand bias assement


def fisher_exact_test(sample_control):
    odds_ratio,p_value = st.fisher_exact(sample_control,alternative='two-sided')
    if p_value < 0.05:
        if odds_ratio > 1.0:
            return "Somatic"
        else:
            return "Germline"
    else:
        return "Germline"

#"template_indicator","length_indicator","mean_size"
def calling_confidence_assement(row):
    if row[0]=="LC":
        return "REJECT" # reject by low templates
    elif row[0] =="MC":
        if row[1]=="HC" and row[2]>2:
            return "PASS_H" #pass by size and template length
        elif row[1]=="MC" or row[2]>2:
            return "PASS_L"
        else:
            return "REJECT" #reject by extreme low mean size
    else:
        #if row[2]>2:
        return "PASS_H"
        #else:
        #    return "PASS_L"
    

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-control","--vardict_control",help="Germline site annotation",type=str)
    ap.add_argument("-polish","--polishing",help="Polishing site annotation",type=str)
    ap.add_argument("-COSMIC","--COSMIC_file",help="COSIMIC site annotation",type=str)
    ap.add_argument("candidate_file",help="candidate from virtual family.py",type=str)
    ap.add_argument("Out_file",help="Result out",type=str)
    args = ap.parse_args()


    result_file = pd.read_csv(args.candidate_file,sep="\t",header=0)
    result_file.dropna(inplace=True)
    result_file[["#chr"]] = result_file["#chr"].astype("str")
    result = result_file.copy()
    result.loc[:,"AF_raw"] = result["var_raw"]/result["depth_raw"]
    result.loc[:,"AF_template"] = result["f1 number"]/result["family_number"]
    result.loc[:,"length_indicator"] = result[["mean_length","std_length"]].apply(lambda row: "HC" if (row[0]+row[1])<=150 else ("MC" if row[0]<=150 else "LC"),axis=1)   
    result.loc[:,"template_indicator"] = result[["f1 number"]].apply(lambda row: "HC" if row[0] >=5 else ("LC" if row[0]<=3 else "MC"),axis=1) 
    result=result.astype({"#chr":str,"start":int,"ref":str,"alt":str})


    if args.vardict_control:
        control_file = pd.read_csv(args.vardict_control,sep="\t",header=None,low_memory=False)
        HWT_info = control_file[[2,3,4,5,6,29,28,22]].copy()
        HWT_info.columns = ['#chr','start','end1','ref','alt','depth','alt_depth','AF']
        HWT_anno = HWT_info.loc[HWT_info["alt_depth"]>4].copy()

        HWT_anno=HWT_anno.astype({"#chr":str,"start":int,"ref":str,"alt":str})
        result = pd.merge(result,HWT_anno,on=["#chr","start","ref","alt"],how="left",indicator="Germline")
        result.loc[:,"genotype"] = result[['f1 number','family_number','alt_depth','depth','Germline']].apply(lambda row:fisher_exact_test(row[0:4].values.reshape(2,2)) if row[4]=="both" else "Somatic",axis=1)
       

    if args.COSMIC_file:
        COSMIC = pd.read_csv(args.COSMIC_file,sep="\t",header=None,low_memory=False)
        COSMIC.columns = ['#chr','start','end1','ref','alt']
        COSMIC=COSMIC.astype({"#chr":str,"start":int,"ref":str,"alt":str})
        COSMIC=COSMIC.drop_duplicates()
        result = pd.merge(result,COSMIC,on=["#chr","start","ref","alt"],how="left",indicator="COSMIC_sites")
        result.loc[:,"COSMIC"] = result["COSMIC_sites"]=="both"


    if args.polishing:
        polishing_file = pd.read_csv(args.polishing,sep="\t",header=None,low_memory=False)
        polishing_file.columns = ["#chr","start","end2","ref","alt","cutoff"]
        polishing_file = polishing_file.astype({"#chr":str,"start":int,"ref":str,"alt":str})
        result = pd.merge(result,polishing_file,on=["#chr","start","ref","alt"],how="left")
        result.loc[:,"stere_noise"] =result[["AF_raw","cutoff"]].apply(lambda row:1 if row[0]<= row[1] else 0,axis=1)
        result = result.loc[result["stere_noise"]==0,:]
        result.loc[:,"somatic_info"]= result[["template_indicator","length_indicator","mean_size"]].apply(lambda row: calling_confidence_assement(row),axis=1)
        result[["#chr","start","end","ref","alt","depth_raw","var_raw","family_number","f1 number","AF_raw","AF_template","mean_size","std_size","mean_length","std_length","singleton_ratio","genotype","COSMIC","somatic_info"]].to_csv(args.Out_file,sep="\t",index = False)
    else:
        result[["#chr","start","end","ref","alt","depth_raw","var_raw","family_number","f1 number","AF_raw","AF_template","mean_size","std_size","mean_template length","std_length","singleton_ratio","genotype","COSMIC","somatic_info"]].to_csv(args.Out_file,sep="\t",index = False)


if __name__=='__main__':
    main()

