Name:Virtual-Barcode-Calling-Algorithm-for-Low-Allelic-Variants-Detection(VBCALAVD)
This algorithm is panel-wide calling algorithm for low AF variants detection from high-depth cfDNA data.
This pilple is consist of three part: Preprocessing,Virtual-family+get-variant-info and Postprocessing

Input file: this file is from Vardict pileup file (code is "VarDict -G ${HUMAN_REF} -p -z 1 -c 1 -S 2 -E 3 -g 4 ${PANEL_BED} -b ${BAM} > ${OUT_FILE}")

Preprocessing is composed of 4 steps: 1) Remove variants less than 4 high quality read supports and Transform variant level data to genimic site level data. 2) Vitual-Unique-Molecular-Identifier(VUMI) remove read-level random error and false family for every chromosome. 3) Variant singleton ratio filter to remove sites with imbalaced variant singleton ratio. 4)Remove candidates in the backlist region such as LC region,STR region.

Virtual-family+get-variant-info: Candidates from Preprocessing step were as input of this step and get detailed virtual family infomation.

Postprocessing: 1) Polishing the steretypical noise. 2) Identify whether a variant is germline or not. 3) Annotation a variant is a COSMIC site. 4) Required 5 variant ctDNA template. For polishing step, we provide the polishing data based on our panel and you can build your own polishing database using polishing_step.py file.

 File output format(tab delimited format):
1. #chr    
2. start   
3. end     
4. ref     
5. alt     
6. depth_raw       # depth
7. var_raw         # variant read numbers
8. family_number   # total virtual family numbers
9. f1 number       # f=1.0 variant virtual family numbers
10. AF_raw         # AF = var_raw/depth_raw
11. AF_template    # AF = f1 number/family_number
12. mean_size      # mean size of f=1.0 virtual family 
13. std_size       # standard deviation of size
14. mean_length    # mean length of f=1.0 virtual family(mean length of ctDNA)
15. std_length     # standard deviation of length 
16. singleton_ratio # variant singleton ratio = singleton/f1 number(size>2)
17. genotype        # germline or Somatic
18. COSMIC         # annotation of a site with COSMIC database
19. somatic_info   # give confidence level of a candidate being a real SNV (REJECT,PASS_L,PASS_H).

Code constitution:
Main code: VBCALAVD.sh
All support Python module files and sh file is in the "Script" directory
Support files: COSMIC annotation file, polishing file and back genomic list region are in the Script/bed_file. if you build polishing database, you should use the polishing-code/polishing_step.py and put the database into the Script/bed_file/polishing
 

Requirement:
Python2.7.11 above; python modules: ast,pysam,collections,pandas,pickle,numpy,scipy,argparse


Parameters and Usage:
Parameters:
-h, --help
    Display this usage message and exit.

-I <val>, --INPUT <val>, --INPUT=<val>
  Input file from pipleup command and Input file is Tab delimited; Recommanded VarDict; 

-O <val>, --OUTP <val>, --OUTP=<val>
  Path for result

-S <val>, --SCRIPT <val>, --SCRIPT=<val>
  Directory of all python modules and sh scripts for calling (Script/)

-V <val>,  --VFS <val>, --VFS=<val>
  Minimun f=1.0 virtual family support; if not present in the option; default is 2.0

-F <val>, --FV <val>, --FV=<val>
  For every non-allelic virtual family,the f value cutoff;if not present in the option; default is 1.0

-U <val>, --UMI <val>, --UMI=<val>
  Work for ctDNA capture NGS data with UMI tag,dedault is False

-P <val>, --PVAL <val>, --PVAL=<val>
  Cutoff FDR value for Variant singleton ratio in the whole panel. Removing outliners with extremely high varint singleton ratio. Default is 0.01

-B <val>, --BASEQ <val>, --BASEQ=<val>
  Minimum basequality for pileupreads. if not present in the option; default is 30

-M <val>, --MAPPINGQ <val>, --MAPPINGQ=<val>
  Minimum mapping quality for pileupreads. if not present in the option; default is 30

-E <val>, --SIZE <val>, --SIZE <val>
  Minimum mean f=1.0 virtual family size.Default is 2.0

-C <val>, --PBL <val>, --PBL=<val>
  Pileup file path for control such as PBL data. The data format is the same with --INPUT file.

-D <val>, --DATAP <val>, --DATAP=<val>
  Path for polishing data; 

-R <val>, --RCOSMIC <val>, --RCOSMIC=<val>
  Path for reference COSMIC sites.

--'
  Treat the remaining arguments as BAM file name.

Example:
./VBCALAMD.sh -I path_Vardict_pileup.txt -C path_PBL_control_pileup.txt -D ./Script/bed_file/polishing/polishing_oncosmart2.txt -R ./Script/bed_file/COSMIC/hg19_cosmic80.txt -S ./Script -O path_for_result -P 0.01 -- path_for_bam 





