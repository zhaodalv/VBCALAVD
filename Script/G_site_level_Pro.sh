#!/bin/env sh

PREPROCESSING()

{
echo "Processing start"
OUT_RESULT=$1
VD_FILE=$2
BAM=$3
ALT_SUPPORT=$4
F_VALUE=$5
UMI=$6
CUTOFF_P_VALUE=$7
#BEDTOOLS_VESION=$8   #eg. v2.24.0 and need  2.21.0 later : alias bedtools=bedtools.v2.24.0

#echo ${OUT_RESULT},${VD_FILE},${BAM},${ALT_SUPPORT},${F_VALUE},${UMI},${CUTOFF_P_VALUE}
if [ ! -d ${OUT_RESULT}/Split_params ]; then

mkdir -p ${OUT_RESULT}/Split_params
mkdir -p ${OUT_RESULT}/Candidate

fi

grep 'SNV' ${VD_FILE} | awk '$29>4 {print $3,$4,$5,$6,$7,$30,$29}' OFS="\t" >${OUT_RESULT}/Mutation_input

python ${SCRIPT_PATH}/G_site_parameter.py ${OUT_RESULT}/Mutation_input >${OUT_RESULT}/Genomic_sites_input 

awk -v var="${OUT_RESULT}/Split_params/" '{print > var $1}' ${OUT_RESULT}/Genomic_sites_input


for file in ${OUT_RESULT}/Split_params/*; do

python ${SCRIPT_PATH}/VUMI_process.py -I ${file} ${BAM} -alt_support ${ALT_SUPPORT} -fratio ${F_VALUE} -UMI ${UMI} > ${OUT_RESULT}/Candidate/$(basename ${file}).candidate&

sleep 2s
done

wait

cat ${OUT_RESULT}/Candidate/*.candidate > ${OUT_RESULT}/Candidate/candidate.txt

#add a control flow
python ${SCRIPT_PATH}/singleton_ratio_fdr.py -mutation ${OUT_RESULT}/Mutation_input -candidate ${OUT_RESULT}/Candidate/candidate.txt -cutoff_P ${CUTOFF_P_VALUE} ${OUT_RESULT}/Candidate/candidate_for_virtual_family ${OUT_RESULT}/Candidate/Vsingleton_filter_candidates


BACKLIST_tar=(${SCRIPT_PATH}/bed_file/Backlist_region/*)
if [ ! -d ${SCRIPT_PATH}/bed_file/Backlist_region/region_files ] && [ -f $BACKLIST_tar ]; then 
mkdir -p ${SCRIPT_PATH}/bed_file/Backlist_region/region_files
tar -xzvf ${BACKLIST_tar} -C ${SCRIPT_PATH}/bed_file/Backlist_region/region_files
fi

BACKLIST=(${SCRIPT_PATH}/bed_file/Backlist_region/region_files/*)

bedtools intersect -a ${OUT_RESULT}/Candidate/candidate_for_virtual_family -b ${BACKLIST[@]} -v >${OUT_RESULT}/Candidate/candidate_after_process

bedtools intersect -a ${OUT_RESULT}/Candidate/candidate_for_virtual_family -b ${BACKLIST[@]} >${OUT_RESULT}/Candidate/candidate_filtered_by_BACKLIST

 
echo "Fished!"

}

