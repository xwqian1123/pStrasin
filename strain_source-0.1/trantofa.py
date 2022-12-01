 # -*- coding: UTF-8 -*- 
import pandas as pd
import sys



def maxrate(item):
	return item.str.cat()

def main(a,b):	
	data = pd.read_table(a,sep='\t')
	transposed = data.T
	transposed = transposed.drop(transposed.index[0])
	transposed['maxnum']=transposed.apply(maxrate,axis=1)
	transposed.to_csv('a.txt',sep='\t')
	transposed = transposed.rename(index=lambda s:'>'+s)
	newfa = transposed['maxnum']
	newfa.to_csv(b,sep='\n',header=None)			




if __name__=='__main__':
	main(sys.argv[1],sys.argv[2])
