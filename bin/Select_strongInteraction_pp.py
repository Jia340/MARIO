#!/usr/bin/env python
# Pengfei Yu yupf05@gmail.com

import sys, argparse
from time import time
from operator import attrgetter
from scipy.stats import hypergeom
from random import shuffle
import math
import pp
import copy
from data_structure import *
from xplib import DBI

from UnionFind import *
# include functions UF, find, merge, connected, count

'''
class annotated_bed():
    def __init__(self,x=None,**kwargs):
        if type(x)==type("str"):
            x=x.split("\t")
        self.chr=str(x[0]).strip()
        self.start=int(x[1])
        self.end=int(x[2])
        try:
            self.seq=str(x[3]).strip()
            self.type=str(x[4]).strip()
            self.name=str(x[5]).strip()
            self.subtype=str(x[6]).strip()
        except:
            pass
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])
    def __lt__(self, other): # for the purpose of ordering clusters
        if self.chr==other.chr:
            return ((self.start<other.start)&(self.end<other.end))
        else:
            return (self.chr<other.chr)
    def overlap(self,other): # for the purpose of finding overlap between regions
        return ((self.chr==other.chr)&(self.end>other.start)&(self.start<other.end))
    def Cluster(self,c):
        self.cluster=c
        #self.c_info=c_info
    def Update(self,S,E): # used for update or expand cluster locations
        self.start=min(self.start,S)
        self.end=max(self.end,E)
    def __str__(self):
        return "\t".join(str(f) for f in [self.chr,self.start,self.end,self.type,self.name,self.subtype,self.cluster])
    # self.cluster become number of regions within cluster for cluster_pool object
'''


##########################################################
#cluster the regions together if they are connected (overlapped) using Union-find approach
def cluster_regions(part,min_clusterS):
    N=len(part)
    cluster_loc=copy.deepcopy(part)
    chro=part[0].chr
    uf_object=UnionFind.UF(N) # union find object
    for i in range(N):
        for j in range(i+1,N):
            if part[i].overlap(part[j]):
                uf_object.merge(i,j)
                #update cluster location (expand)
                cluster_loc[uf_object.find(j)].Update(min(part[i].start,part[j].start),max(part[i].end,part[j].end))
            if part[i]<part[j]:
                break
        if i%10000==0: print >> sys.stderr, "  Merging segment for clusters %s,(%d/%d)\r"%(chro,i,N),
    
    c_pool=[] # cluster pool
    for i in range(N):
        c=uf_object.find(i)
        #c_info="%d:%d"%(cluster_loc[c].start,cluster_loc[c].end)
        part[i].Cluster(chro+".%d"%(c))
        c_pool.append(c)

    cluster_pool={}
    for c in set(c_pool):
        count=c_pool.count(c)
        if count>=min_clusterS:
            cluster_pool[chro+".%d"%(c)]=cluster_loc[c]
            cluster_pool[chro+".%d"%(c)].cluster=count
    return (cluster_pool,part,chro)


def Random_strongInteraction(part1,part2,cluster_pool1,cluster_pool2):
    global min_interaction, p_value
    ''' This is for counputing FDR using random permutation '''
    c_interaction={}
    for i in range(len(part1)):
        region1=str(part1[i])
        region2=str(part2[i])
        inter="%d--%d"%(part1[i].cluster,part2[i].cluster)
        if c_interaction.has_key(inter):
            c_interaction[inter]+=1
        else:
            c_interaction[inter]=0

    k=0 # record for strong interactions
    n=0
    for interaction in c_interaction:
        n=n+1
        count=c_interaction[interaction]
        if count<min_interaction: continue
        i=int(interaction.split("--")[0])
        j=int(interaction.split("--")[1])
        try:  # we select clusters with size no less than 5, so some interactions cannot be found in clusters
            count1=cluster_pool1[i].cluster
            count2=cluster_pool2[j].cluster
        except:
            continue
        real_p=1-hypergeom.cdf(count,len(part1),count1,count2)
        if real_p<=p_value:
            k=k+1
    return [n,k]

#parameters

def ParseArg():
    p=argparse.ArgumentParser(description="find strong interactions from paired genomic location data",epilog="need Scipy for hypergeometric distribution")
    p.add_argument("-i","--input",type=str,required=True,help="input file which is the output file of Stitch-seq-Aligner.py")
    p.add_argument("-M","--min_clusterS",type=int,default=5,help="minimum number of segments allowed in each cluster, default:5")
    p.add_argument("-m","--min_interaction",type=int,default=3,help="minimum number of interactions to support a strong interaction, default:3")
    p.add_argument('-p',"--p_value",type=float,default=0.05,help="the p-value based on hypergeometric distribution to call strong interactions, default: 0.05")
    p.add_argument('-o','--output',type=str,help="specify output file")
    p.add_argument("-P","--parallel",type=int,default=5,help="number of workers for parallel computing, default: 5")
    p.add_argument("-F","--FDR",action='store_true', help="Compute FDR if specified")
    p.add_argument('-a','--annotation',type=str,default="/home/yu68/bharat-interaction/new_lincRNA_data/all_RNAs-rRNA_repeat.txt",help='If specified, include the RNA type annotation for each aligned pair, need to give bed annotation RNA file')
    p.add_argument("-A","--annotationGenebed",dest="db_detail",type=str,default="/home/yu68/bharat-interaction/new_lincRNA_data/Ensembl_mm9.genebed",help="annotation bed12 file for lincRNA and mRNA with intron and exon")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        sys.exit(0)
    return p.parse_args()


def SingleFragment(p1,p2):
    '''
    determine if mapped pair is single RNA fragment or not.
      1. no linker
      2. same chromosome
      3. same RNA
    '''
    if (p1.end-p1.start>=78) and (p2.end-p2.start>=88) and (p1.chr==p2.chr) and p1.name==p2.name:
        return True
        



def Main():
    t1=time()
    
    global min_interaction, p_value
    args=ParseArg()
    inp = open(args.input, 'r')
    min_clusterS=args.min_clusterS
    min_interaction=args.min_interaction
    p_value=args.p_value
    output=open(args.output,'w')
    ncpus=args.parallel


    #store genomic location of part1 and part2
    part1=[]
    part2=[]


    k=0
    
    print >> sys.stderr,"# Inputing data..."

    chr_list=[]
    for line in inp.read().split('\n'):
        if line=='': continue
        line=line.strip().split('\t')
        p1=annotated_bed(line[0:10],id=k)
        p2=annotated_bed(line[11:],id=k)
        if isinstance(p1.start, list):
            p1.start=int(p1.start[0])
            p1.end=int(p1.end[-1])
        if isinstance(p2.start, list):
            p2.start=int(p2.start[0])
            p2.end=int(p2.end[-1])
        if SingleFragment(p1,p2): continue
        k+=1
        part1.append(p1)
        part2.append(p2)
        if p1.chr not in chr_list: chr_list.append(p1.chr)
        if p2.chr not in chr_list: chr_list.append(p2.chr)
        if k%20000==0: 
            print >> sys.stderr,"  Reading %d pairs of segments\r"%(k),
    print >> sys.stderr,"Get total %d pairs."%(k)
    
    if len(part1)!=len(part2):
        print >> sys.stderr, "## ERROR: number of regions in two part not match!!"
        sys.exit(0)

    # sort in genomic order, easy for clustering
    part1=sorted(part1, key=attrgetter('start'))
    part1=sorted(part1, key=attrgetter('chr'))
    part2=sorted(part2, key=attrgetter('start'))
    part2=sorted(part2, key=attrgetter('chr'))

    # for parallel computing 
    print >>sys.stderr,"# Generating clusters for two parts..."
    # tuple of all parallel python servers to connect with
    ppservers = ()
    job_server = pp.Server(ncpus, ppservers=ppservers)
    jobs1=[]
    jobs2=[]
    for chro in chr_list:
        part1_temp=filter(lambda p: p.chr==chro, part1)
        if len(part1_temp)>0:
            jobs1.append(job_server.submit(cluster_regions,(part1_temp,min_clusterS),(annotated_bed,),("UnionFind","copy",)))
        part2_temp=filter(lambda p: p.chr==chro, part2)
        if len(part2_temp)>0:
            jobs2.append(job_server.submit(cluster_regions,(part2_temp,min_clusterS),(annotated_bed,),("UnionFind","copy",)))
        

    cluster_pool1={}
    part1=[]
    for job in jobs1: 
        try:
            part1=part1+job()[1]
            cluster_pool1.update(job()[0])
        except:
            print >> sys.stderr, "Wrong in %s, part1"%(job()[2])
            continue
    cluster_pool2={}
    part2=[]
    for job in jobs2:
        try:
            part2=part2+job()[1]
            cluster_pool2.update(job()[0])
        except:
            continue


    print >>sys.stderr,"   cluster number for part1 is %d          "%(len(cluster_pool1))
    print >>sys.stderr,"   cluster number for part2 is %d          "%(len(cluster_pool2))

    # sort back to pair two parts
    part1=sorted(part1, key=attrgetter('id'))
    part2=sorted(part2, key=attrgetter('id'))

    print >> sys.stderr,"size of part1&2:",len(part1),len(part2)

    c_interaction={}
    for i in range(len(part1)):
        region1=str(part1[i])
        region2=str(part2[i])
        try:
            inter=part1[i].cluster+"--"+part2[i].cluster
        except:
            print >> sys.stderr,i,part1[i].cluster,part2[i].cluster
            sys.exit()
        if c_interaction.has_key(inter):
            c_interaction[inter]+=1
        else:
            c_interaction[inter]=1

    # annotation file
    print >> sys.stderr,"# Indexing annotation files"
    dbi_all=DBI.init(args.annotation,"bed")
    dbi_detail=DBI.init(args.db_detail,"bed")
    dbi_repeat=DBI.init("/home/yu68/bharat-interaction/new_lincRNA_data/mouse.repeat.txt","bed")


    print >> sys.stderr,"# finding strong interactions from clusters..."
    k=0 # record for strong interactions
    n=0

    # annotation file

    for interaction in c_interaction:
        n=n+1
        count=c_interaction[interaction]
        if count<min_interaction: continue
        i=interaction.split("--")[0]
        j=interaction.split("--")[1]
        try:  # we select clusters with size no less than 5, so some interactions cannot be found in clusters
            count1=cluster_pool1[i].cluster
            count2=cluster_pool2[j].cluster
        except:
            continue
        real_p=1-hypergeom.cdf(count,len(part1),count1,count2)
        if real_p<=p_value:
            k=k+1
            cluster_pool1[i].Annotate(dbi_all,dbi_detail,dbi_repeat)
            cluster_pool2[j].Annotate(dbi_all,dbi_detail,dbi_repeat)
            try:
                log_p = math.log(real_p)
            except:
                log_p = -float("Inf")
            print >> output,str(cluster_pool1[i])+'\t'+str(cluster_pool2[j])+'\t%d\t%.4f'%(count,log_p)
        if n%1000==0: print >> sys.stderr, "  Progress ( %d / %d )\r"%(n,len(c_interaction)),

    print >> sys.stderr,"# Find %d strong interactions. Cost time: %.2f s"%(k,time()-t1)

    if args.FDR:
        print >> sys.stderr, "# Permutated results:"
        for i in range(10):
            shuffle(part2)
            [n_r_I,n_r_SI]=Random_strongInteraction(part1,part2,cluster_pool1,cluster_pool2)
            print >> sys.stderr, "  ",i, n_r_I, n_r_SI, n_r_SI*1.0/n_r_I
            

if __name__=="__main__":
    Main()
