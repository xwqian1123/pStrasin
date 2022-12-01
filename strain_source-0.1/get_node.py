import os
import argparse
import sys
#import pandas as pd


def main(treefile,matrix,outgrup):
    os.system('/mnt/beegfs/user/wuyr/01.Software/bin/R --slave --file=/mnt/beegfs/user/qianxw/03.qianxw/Project/06.sim/01.script/ggtree_node_table.r --args {} {} {} '.format(treefile,matrix,outgrup))
    isfile = os.path.exists('result.csv')
    if isfile:
        with open ('get_all.sh','w') as d:
            with open ('get_new.sh','w') as n:
                with open ('result.csv','r') as f:
                    line=f.readlines()
                    for i in line:
                        arr=i.strip().split('\t')
                        d.write('awk -F "\\t" '+" '{if("+arr[0]+"){print "+arr[1]+"}}' "+matrix+"|sed -n '2,$p' > " + 'all_node/'+arr[2]+'_'+arr[4]+'.node.txt'+'\n')
                        n.write('cat '+ 'all_node/'+arr[2]+'_'+arr[4]+'.node.txt '+ 'all_node/'+arr[2]+'_'+arr[3]+'.node.txt '+'| sort | uniq -u > '+'new_node/'+arr[3]+'_'+arr[4]+'.node.txt'+'\n')
    get_awk_file = os.path.exists('get_all.sh')
    get_cat_file = os.path.exists('get_new.sh')
    if get_cat_file & get_awk_file:
        os.system('sh {}'.format('/mnt/beegfs/user/qianxw/03.qianxw/Project/06.sim/01.script/get_sh.sh'))
    os.system('rm -rf {}'.format('result.csv'))
    os.system('rm -rf {}'.format('get_all.sh'))
    os.system('rm -rf {}'.format('get_new.sh'))
    os.system('rm -rf {}'.format('parent.sh'))
    os.system('rm -rf {}'.format('new_node'))
    os.system('rm -rf {}'.format('empty_node'))
    os.system('rm -rf {}'.format('all_node'))


if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3])
