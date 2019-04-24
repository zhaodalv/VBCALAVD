import sys
import scipy.stats as st


def xopen(filename, mode='r'):
    assert isinstance(filename, basestring)

    if filename == '-':
        if 'r' in mode:
            return sys.stdin
        elif 'w' in mode or 'a' in mode:
            return sys.stdout
        else:
            raise ValueError('Wrong Mode Type: {0}'.format(mode))
    elif filename == '@' and ('w' in mode or 'a' in mode):
        return sys.stderr

    if filename.endswith('.gz'):
        import gzip
        return gzip.open(filename, mode)
    elif filename.endswith('.bz') or filename.endswith('.bz2'):
        import bz2
        return bz2.BZ2File(filename, mode=mode)

    return open(filename, mode)

def Ds_filter(start,end,query_position):
    Ds = query_position-start
    De = end-query_position
    return min(Ds,De)


def get_fdr(Pval,queue_num,length):
    fdr=Pval*(length/queue_num)
    if fdr <=1.0:
        return fdr
    else:
        return 1.0


def get_Pval(ratio,params,dist=st.johnsonsu):
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]
    Pval =dist.sf(ratio,loc=loc,scale=scale,*arg)
    return Pval


def calculate_P_value(x,param):
    arg = param[:-2]
    loc = param[-2]
    scale = param[-1]
    dist = st.johnsonsu(*arg, loc=loc, scale=scale) 
    p = dist.sf(x)
    return p

def FDR_adjust(sorted_P_values):
    V = len(sorted_P_values)
    out=[]
    prev=1.0
    queue_num=V
    for value in sorted_P_values:
        fdr = value*(V/float(queue_num))
        #print fdr
        if fdr<prev:
            prev=fdr
        out.append(prev)
        queue_num=queue_num-1
    return out

    
