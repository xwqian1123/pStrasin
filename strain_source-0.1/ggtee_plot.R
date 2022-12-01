library('ggplot2')
library('ggtree')
library('treeio')
library("ape")
library(tibble)
library(geiger)
library(tidyr)

args <- commandArgs(T)
pdf_file<-paste(args[2],'.pdf',sep='')
n_csv<-paste(args[2],'.csv',sep='')
filt_file<-paste(args[2],'.filt.txt',sep='')
pdf(pdf_file,width=20,height=30)
tree = read.tree('SNP.matrix.fa.treefile')
node_name = read.table('./node_file.txt',header=F,sep=' ')

branch<-read.csv(args[3],sep='\t',head=F)
branch<-unite(branch,V4,c(V1,V2), sep="_", remove = T)
g <- split(branch$V4,list( branch$V3))
p <- ggtree(tree)
clades <- sapply(g, function(n) MRCA(p, n))
cols=c("#580000", "#EC762F", "#C5560A", "#9E3600","#791200","#FFCB7F")
p <- groupClade(p, clades) +geom_tiplab(aes(label=label),color='#000000',size = 4)+ scale_color_manual(values = c(cols)) 
treedata=fortify(tree)
treenode<-unique(node_name[,4])

#bb<-data.frame()
get_X<-function(a,b){
    print (a)
    p<<- p+geom_point2(aes(subset=(node==a)),size=6,color='#ee6666',shape=b)
}
get_loc<-function(a){
    child_name<-child(tree, a)
    for (i in child_name){
	new_empty=length(subset(sampledataVa ,as.numeric(sampledataVa$after)==i)$after)
	new_1=length(subset(sampledataVa ,as.numeric(sampledataVa$after)==i & sampledataVa$istype ==1 )$after)
	new_0=length(subset(sampledataVa ,as.numeric(sampledataVa$after)==i & sampledataVa$istype ==0 )$after)
	#cha=new_1-new_0
#----------V-------
	if (new_0!=0){
            print ('dddd')
            get_X(i,4)
	}else{
	    get_loc(i)
        }
#---------X------
        #if(new_0!=0){
          #  get_X(i,23)
         #   get_loc(i)
	#}else{
          #  next
        #}
    }
}

outgrup=args[1]
root_name<-match(outgrup,treedata$label)
root_node<-parent(tree,root_name)
#get_loc(root_node)

node <- unique(node_name[,4])

get_pie <-function (){
    dd <- data.frame()
    for (i in node){
            num=length(subset(node_name,node_name[,4]==i)$V4)
            new=length(subset(sampledata ,sampledata$node==i)$CHROM)
            new_0=length(subset(sampledata ,sampledata$node==i & sampledata$istype == 0 )$node)
            new_1=length(subset(sampledata ,sampledata$node==i  & sampledata$istype == 1)$node)
            no=num-new
            if (new == 0){
		print (new)
	    }else{	
	    	cc<-c(i,num,new_0,new_1,no)
            	percent <- round(as.numeric(cc[3:5])/as.numeric(num)*100, 2)
         #   label <- paste((piename), "(", percent, "% )")
            # print (cc) 
            	nd<-strsplit(i,'_')[[1]][2]
            	ee<-c(as.numeric(percent),as.numeric(nd))
            	dd <-rbind(dd,ee)
	    }
                # 定制标签
        }
    names(dd)=c("ANC", "DER","no",'node')
    write.csv(dd,file=n_csv)
    return (dd)
}
#a<-get_pie()
#pies <- nodepie(a, cols = 1:3)

#p<-ggtree(tree)+geom_tiplab()
#pies <- lapply(pies, function(g) g+scale_fill_manual(values = cols))
#p2 <- p + geom_inset(pies, width = .05, height = .05,x = "branch")
#print(p2)
a<-file.exists(filt_file)
if (a=='FALSE'){
        print ('file no exist')
}else{
	sampledata <<- read.table(filt_file,header=T,sep='\t')
	sampledataVa <<- sampledata %>%separate(node, c("before", "after"), "_")
	get_loc(root_node)
	a<-get_pie()
	pies <- nodepie(a, cols = 1:3,color=c('#eb746a','#28ad40','#d8d6d6'))
	p2 <- p + geom_inset(pies, width = .05, height = .05,x = "branch")+ theme(legend.position = "right") 
	print(p2)
}

