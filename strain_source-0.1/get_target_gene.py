 # -*- coding: UTF-8 -*- 
import pandas as pd
import re
import sys
import operator 
import os

def locafun(item):
    i=item['MSTree_Branch'].split('_')[1]
    branch_list = ['N01','N02','N03','N04','N05','N06','N07','N08','N09','N10','4.ANTa'] 
    if i in branch_list:
        return '2 ' #2=branch
    else:
        string=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s']
        a=i[-1]
        #result= re.match(str,i)
        if a in string:
            # print (i,'S') 
             return '0' ##0=strain
        else:
            # print (i,'P')
             return '1' #1=phygrou

def func_ALT_apply(item):
    i=item['ALT'].replace('<*>',item["REF"])
    AD=item['INFO'].split(';')[1].replace('AD=','')
    ADarr=AD.split(',')
    if int(ADarr[0])== 0:
         i=i[:-2]    
    return i
    
def typefun(item):
    rate=item['RATIO'].split(',')
    maxrate_index,max_rate=max(enumerate(rate),key=operator.itemgetter(1))
    #print (item['REF'],'ref')
    #print (maxrate_index,max_rate)
    alt=item['ALT'].split(',')
    #print (alt)
    if maxrate_index == 0:
        ALT_type = item['REF']
    else:
        index = maxrate_index-1
        alt_type=operator.itemgetter(index)
        ALT_type=alt_type(alt)
    #print (ALT_type,'alt')
    if ALT_type == item['ANC_NT']:
         istype='0'
    elif ALT_type == item['DER_NT']:
         istype='1'
    else:
         istype='?'
#    print (ALT_type,item['REF'],istype)
    return istype	
    
   # if (re.match(',',item['ALT'])):
   #     print(item)
   # else:
   #     if item['ALT']==item['Anc_Nt']:
   #         istype='1'
   #     elif item['ALT']==item['Der_Nt']:
   #         istype='0'
   #     else:
   #         istype=item['ALT']
   # return istype

def qualfun(item):
    arr=''
    if len(item['qual'])<2:
        n=int(ord(item['qual'])-33)
        if n>20:
            arr=str(n)+','
        else:
            arr=str(0)+','	
    else:
        for i in item['qual']:
            n =str(ord(i)-33)
            arr=arr+n+','
    tmp=arr[:-1]
    return arr
def all_dbfun(item):
    aLL_DP=item['INFO'].split(';')[0][3:]
    #print (aLL_DP,'a')
    return aLL_DP    

def ratefun(item):
    AD=item['INFO'].split(';')[1][3:]
    ADarr=AD.split(',')
    arr2=[ int(x) for x in ADarr]
    sum_AD=sum(arr2)
    rate=''
    for i in ADarr:
        if int(sum_AD)==0:
            r=0
        else:
            r = round(int(i)/int(sum_AD),2)
       	rate=rate+str(r)+','
    tmp= rate[:-1]
    return tmp

def ref_fun(item):
    AD=item['INFO'].split(';')[1][3:]
    ref_DP=AD.split(',')[0]
    return ref_DP

def filt_dpfun(item):
    AD=item['INFO'].split(';')[1][3:]
    ADarr=AD.split(',')
    arr2=[ int(x) for x in ADarr]
    n=sum(arr2)
    return n
    

def alt_fun(item):
    AD=item['INFO'].split(';')[1][3:]
    alt_DP=AD.split(',')[1:]
    return alt_DP

def maxrate(item):
    arr=item['RATIO'].split(',')
    maxnum=max(arr)
    a=float(maxnum)
    return a

def main(a,c,d):
    bcf_file=c+'.tmp.bcf'
    bcf_path=os.path.join(a,bcf_file)
    size_b = os.path.getsize(bcf_path)
    if size_b == 0:
        print (bcf_file+" input file is empty")
    else:

    	YStype=pd.read_csv('./node_file.txt',header=None,names=['POS','DER_NT','ANC_NT','NODE'],sep=' ',index_col=0,encoding='unicode_escape')
    	snp=pd.read_table(bcf_path,header=None,names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','MQ'],sep='\t',index_col=1)
    	fina=pd.merge(snp,YStype,how='inner',on='POS')
         
    	if len(fina) == 0:
            print(c +"isnot posistion")
    	else:
            stat_name=c+'.txt'
            star_file=os.path.join(a,stat_name)
            fina.to_csv(star_file,sep='\t')    
        
            fina['ALL_DP']=fina.apply(all_dbfun,axis=1)
            fina['FILT_DP']=fina.apply(filt_dpfun,axis=1)
            fina['REF_DP']=fina.apply(ref_fun,axis=1)
            fina['ALT_DP']=fina.apply(alt_fun,axis=1)
            fina['RATIO']=fina.apply(ratefun,axis=1)
            fina['TYPE']=fina.apply(typefun,axis=1)
            fina['MAX_HETERNUM']=fina.apply(maxrate,axis=1)
           
            fina.drop(columns=['INFO','FORMAT','QUAL','FILTER','ID','MQ'],axis=1,inplace=True)
            filt_name=c+'.filt.txt'
            filt_file=os.path.join(a,filt_name)
            fina=fina.drop(fina[fina['MAX_HETERNUM']<0.90].index)
            d=int(d)
            fina=fina.drop(fina[fina['FILT_DP']<d].index)
            fina.to_csv(filt_file,sep='\t')

if __name__=='__main__':

    main(sys.argv[2],sys.argv[1],sys.argv[3])
