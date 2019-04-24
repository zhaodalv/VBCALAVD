#/bin/env sh


BAM_PATH=""
VD_INPUT=""
RESULT_PATH=""
SCRIPT_PATH=""

ALT_SUPPORT=2.0
F_VALUE=1.0
UMI=False
CUTOFF_P_VALUE=0.01
BASE_Q=30
MAPPING_Q=30

CONTROL_PBL=False
DATABASE=False
SIZE=2


usage() {
    cat <<EOF
Usage: $0 [options] [--] [file...]

Arguments:

  -h, --help
    Display this usage message and exit.

  -I <val>, --INPUT <val>, --INPUT=<val>
    Input file from pipleup command and Input file is Tab delimited; Recommanded VarDict; 
    Format:chr start end ref alt.

  -O <val>, --OUTP <val>, --OUTP=<val>
    Path for result

  -S <val>, --SCRIPT <val>, --SCRIPT=<val>
    Derectory of all python scripts for calling
  
  -V <val>,  --VFS <val>, --VFS=<val>
    Minimun f=1.0 virtual family support; if not present in the option; default is 2.0
 
  -F <val>, --FV <val>, --FV=<val>
    For every non-allelic virtual family,the f value cutoff;if not present in the option; default is 1.0

  -U <val>, --UMI <val>, --UMI=<val>

    Work for ctDNA capture NGS data with UMI tag

  -P <val>, --PVAL <val>, --PVAL=<val>

    Cutoff P value for Variant singleton ratio in the whole panel. Removing outliners with extremely high varint singleton ratio.

  -B <val>, --BASEQ <val>, --BASEQ=<val>
 
    Minimum basequality for pileupreads. if not present in the option; default is 30

  -M <val>, --MAPPINGQ <val>, --MAPPINGQ=<val>

    Minimum mapping quality for pileupreads. if not present in the option; default is 30
  
  -E <val>, --SIZE <val>, --SIZE <val>
    
    Minimum mean f=1.0 virtual family size. 

  -C <val>, --PBL <val>, --PBL=<val>

    Pileup file path for control such as PBL data. The data format is the same with --INPUT file.

  -D <val>, --DATAP <val>, --DATAP=<val>
    
    Path for polishing data; 

  -R <val>, --RCOSMIC <val>, --RCOSMIC=<val>

    Path for reference COSMIC sites.

  --'

    Treat the remaining arguments as BAM file name.
    file name might begin with '-'.

  file...
    Optional list of file names.  If the first file name in the list
    begins with '-', it will be treated as an option unless it comes
    after the '--' option.
EOF
}


log() { printf '%s\n' "$*"; }
error() { log "ERROR: $*" >&2; }
fatal() { error "${arg} duplicated"; exit 1; }
usage_fatal() { error "$*"; usage >&2; exit 1; }


while [ "$#" -gt 0 ]; do
    arg=$1
    case $1 in
        # convert "--opt=the value" to --opt "the value".
        # the quotes around the equals sign is to work around a
        # bug in emacs' syntax parsing
        --*'='*) shift; set -- "${arg%%=*}" "${arg#*=}" "$@"; continue;;
        -I|--INPUT)  shift; if [ x${VD_INPUT} != x ]; then fatal; else VD_INPUT=$1; fi ;;
        -O|--OUPT)   shift; if [ x${RESULT_PATH} != x ]; then fatal; else RESULT_PATH=$1; fi ;;
        -S|--SCRIPT) shift; if [ x${SCRIPT_PATH} != x ]; then fatal; else SCRIPT_PATH=$1; fi ;;
        -V|--VFS)    shift; ALT_SUPPORT=$1;;
        -F|--FV)     shift; F_VALUE=$1;;
        -U|--UMI)    shift; UMI=$1;;
        -P|--PVAL)   shift; CUTOFF_P_VALUE=$1;;
        -B|--BASEQ)  shift; BASE_Q=$1;;
        -M|--MAPPINGQ) shift; MAPPING_Q=$1;;
        -E|--SIZE)   shift; SIZE=$1;;
        -C|--PBL)    shift; CONTROL_PBL=$1 ;;
        -D|--DATAP) shift; DATAP=$1 ;;
        -R|--RCOSMIC)  shift; RCOSMIC=$1 ;;
        -h|--help) usage; exit 0;;
        --) shift; BAM_PATH=$1;break;;
        -*) usage_fatal "unknown option: '$1'";;
        *) break;; # reached the list of file names
    esac
    shift || usage_fatal "option '${arg}' requires a value"
done

#echo ${VD_INPUT},${RESULT_PATH},${SCRIPT_PATH},${ALT_SUPPORT},${F_VALUE},${UMI},${CUTOFF_P_VALUE},${BASE_Q},${MAPPING_Q},${SIZE},${CONTROL_PBL},${DATAP},${RCOSMIC},${BAM_PATH}


python ${SCRIPT_PATH}/find_library.py ${SCRIPT_PATH}


source ${SCRIPT_PATH}/G_site_level_Pro.sh

#first step is VUMI propresssing:VUMI + Ds filter + variant singleton ratio + genome backlist
PREPROCESSING ${RESULT_PATH} ${VD_INPUT} ${BAM_PATH} ${ALT_SUPPORT} ${F_VALUE} ${UMI} ${CUTOFF_P_VALUE}


source ${SCRIPT_PATH}/calling_and_post_processing.sh

#get virtual family info: get virtual family and get variant info for every high confient variants from preprocessing
VIRTUR_FAMILY ${BASE_Q} ${MAPPING_Q} ${CONTROL_PBL} ${DATAP} ${RCOSMIC}




