#!/usr/bin/env python
import numpy as np
import collections
import scipy.stats as st
import warnings
import operator
import pandas as pd
import argparse

def best_fit_distribution(data, bins=200, b_value=0.3):
    sse_value = collections.defaultdict(list)
    fit_param = collections.defaultdict(list)
    y, x = np.histogram(data, bins=bins, density=True)
    x = (x + np.roll(x, -1))[:-1] / 2.0
    b_best=b_value
    distname=""
    out_param =()
    DISTRIBUTIONS = [st.norm, st.johnsonsu,st.beta,st.dweibull,st.weibull_min,st.alpha,st.gamma,st.lognorm,st.exponnorm,st.nct]
    for distribution in DISTRIBUTIONS:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            params = distribution.fit(data)
            arg = params[:-2]
            loc = params[-2]
            scale = params[-1]
            pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
            sse = np.sum(np.power(y - pdf, 2.0))
            sse_value[distribution].append(sse)
            fit_param[distribution].append(params)
    sorted_sse = sorted(sse_value.items(),key = operator.itemgetter(1))[0:3]
    for distribution,see in sorted_sse:
        params = distribution.fit(data)
        a,b = st.probplot(data,sparams =params , dist=distribution)
        if b[-1] > b_best :
            b_best=b[-1]
            distname = distribution
    return (distname,b_best)

def calculate_probability(row_data,distname):
    row=row_data
    length=len(row)
    param = distname.fit(row)
    arg = param[:-2]
    loc = param[-2]
    scale = param[-1]
    sig_greater_than_999 = distname.ppf(1-(0.05/length), loc=loc, scale=scale,*arg) 
    return sig_greater_than_999,param

if __name__=='__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument("AF_file",help="background AF sample",type=str)
    args = ap.parse_args()
    data = pd.read_csv(args.AF_file,sep="\t",index_col=0)
    for index,value in data.iterrows():
        data_in = value[~pd.isnull(value)]
        length = len(data_in)
        name,b_best= best_fit_distribution(np.array(data_in),bins=200,b_value=0)
        if name:
            cutoff,par =calculate_probability(np.array(data_in),name)
            print index,"\t",cutoff,"\t",max(data_in)
        else:
            print "NA"
