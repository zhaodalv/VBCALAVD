#!/usr/bin/env sh

VIRTUR_FAMILY()
{
BASE_Q=$1
MAPPING_Q=$2
CONTROL_PBL=$3
POLISHING=$4
COSMIC=$5


echo "Virtual family start"
if [ ! -d ${RESULT_PATH}/Result ];then
mkdir -p ${RESULT_PATH}/Result
fi

python ${SCRIPT_PATH}/virtual_family.py -I ${RESULT_PATH}/Candidate/candidate_after_process -bq ${BASE_Q} -mq ${MAPPING_Q} ${BAM} >${RESULT_PATH}/Result/HC_candidate

python ${SCRIPT_PATH}/post_processing.py -control ${CONTROL_PBL} -polish ${POLISHING} -COSMIC ${COSMIC} ${RESULT_PATH}/Result/HC_candidate ${RESULT_PATH}/Result/Calling_info.txt

echo "Virtual family finshed!"
}
