import sys,argparse,os
from xplib import DBI
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np
plt.ioff()

def ParseArg():
    p=argparse.ArgumentParser(description="plot interactions using a heatmap. information of linked pairs are stored in file '*_fragment_paired_align.txt'",epilog="Require: matplotlib, numpy")
    p.add_argument("interaction",type=str,help="Interaction file from output of 'Select_strongInteraction_pp.py'")
    p.add_argument("linkedPair",type=str,help="file for information of linked pairs, which is output of 'Stitch-seq_Aligner.py'")
    p.add_argument('-n','--name',type=str,help="give a gene name and plot the interaction heatmap new the gene region, exclusive with '-r' option")
    p.add_argument('-r',type=str,help="Choose region to plot, give region with format 'chr:start-end', exclusive with '-n' option")
    p.add_argument('-s','--start',type=int,nargs='+',default=(7,11),help='start column number of the second region in interaction file and linkedPair file, default=(7,11)')
    p.add_argument('-g','--genebed',type=str,default='/home/yu68/bharat-interaction/new_lincRNA_data/Ensembl_mm9.genebed',help='the genebed file from Ensembl, default: Ensembl_mm9.genebed')
    p.add_argument("-p","--pair_dist",type=int,default=1000,help="two interacted parts within this distance are considered as self-ligated and they are marked or eliminated (see option -s for slim mode), default: 1000bp")
    p.add_argument("-S","--Slim",action='store_true',help='set slim mode to eliminate self ligated interactions')
    p.add_argument("-t","--step",type=int, default=10, help="step or resolution or unit size of the heatmap, default=10bp")
    p.add_argument("-I","--SI",action='store_true',help='Specify to add strong interaction in the figure,default: False')
    p.add_argument('-o','--output',type=str,help="output plot file, can be format of emf, eps, pdf, png, ps, raw, rgba, svg, svgz")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        sys.exit(0)
    return p.parse_args()

class Bed:
    def __init__(self,x,**kwargs):
        self.chr=x[0]
        self.start=int(x[1])
        self.stop=int(x[2])
        if len(x)>=6:
            self.type=x[3]
            self.name=x[4]
            self.subtype=x[5]
        self.center=int((self.start+self.stop)/2)
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])
    def overlap(self,bed2,n):
        return (self.chr==bed2.chr) and (self.start+n<bed2.stop) and (bed2.start+n<self.stop)
    def str_region(self):
        return self.chr+":%d-%d"%(self.start,self.stop)

def read_region(Str):
    ''' input str has format 'chr:start-end'  '''
    chro = Str.split(':')[0]
    start = Str.split(':')[1].split("-")[0]
    end = Str.split(':')[1].split("-")[1]
    return Bed([chro,start,end,Str,".","."])

def read_interaction(File,s):
    '''
    s: start column number for second part of interaction
    '''
    a=open(File,'r')
    for l in a.read().split('\n'):
        if l.strip()=="": continue
        lsep=l.split('\t')
        if lsep[3] in ['+','-']:
            bed1=Bed(lsep[0:3],strand=lsep[3])
            bed2=Bed(lsep[s:(s+3)],strand=lsep[s+3])
        else:
            bed1=Bed(lsep[0:3])
            bed2=Bed(lsep[s:(s+3)])
        yield (bed1,bed2,l)


def addGeneToFig(gene,ax,start,end,name=0,bottom=0):

    '''
    add gene to figures
        start is the start of query region
    '''
    if name==1:
        if gene.start>start:
            ax.text((1.0*gene.start+0.0*start),bottom+0.015,gene.id.split("&")[0],fontsize=6,horizontalalignment='right')
        else:
            ax.text((1.0*gene.stop+0.0*end),bottom+0.015,gene.id.split("&")[0],fontsize=6,horizontalalignment='left')
    cds=gene.cds()

    utr5=gene.utr5()
    utr3=gene.utr3()
    if cds.stop!=cds.start:
        cds_exons=cds.Exons()
        for cds_exon in cds_exons:
            ax.bar(cds_exon.start,0.02,cds_exon.stop-cds_exon.start,facecolor="blue",edgecolor="blue",alpha=1,bottom=bottom,lw=0.1)
    if not utr3 is None:
        for utr3_exon in utr3.Exons():
            ax.bar(utr3_exon.start,0.01,utr3_exon.stop-utr3_exon.start,facecolor="blue",edgecolor="blue",alpha=1,bottom=bottom+0.005,lw=0.1)
    if not utr5 is None:
        for utr5_exon in utr5.Exons():
            ax.bar(utr5_exon.start,0.01,utr5_exon.stop-utr5_exon.start,facecolor="blue",edgecolor="blue",alpha=1,bottom=bottom+0.005,lw=0.1)
    interval=(end-start)/100
    yloc=bottom+0.01
    for intron in gene.Introns():
        ax.plot([intron.start,intron.stop],[yloc,yloc],lw=0.5,color='k') 
        for i in range((intron.stop-intron.start)/interval):
            if intron.strand=="+":
                loc=intron.start+(i+1)*interval
                x=[loc-0.3*interval,loc,loc-0.3*interval]
            else:
                loc=intron.stop-(i+1)*interval
                x=[loc+0.3*interval,loc,loc+0.3*interval] 
            y=[yloc-0.01,yloc,yloc+0.01]
            ax.plot(x,y,color='k',lw=0.5)

def Genetrack(bed,gene_dbi,ax,track_bottom):
    '''
    draw gene track for each part
       bed: define the region to be drawn
       gene_dbi: querible database for genes
       ax: the subplot
    '''
    query_gene=gene_dbi.query(bed)
    ### determine height of gene track

    bottoms=[0 for i in range(100)]
    max_index=0
    for i in gene_dbi.query(bed):
        index=0
        while(1):
            if i.start > bottoms[index]:
                addGeneToFig(i,ax,bed.start,bed.stop,1,0.03*index+track_bottom+0.05)
                bottoms[index]=i.stop
                if max_index<index: max_index=index
                break
            index+=1
    gene_track_number=max_index+1
    gene_track_height=0.03*gene_track_number+0.04

    # add frame
    rect=matplotlib.patches.Rectangle((bed.start,track_bottom+0.03),bed.stop-bed.start, gene_track_height, edgecolor='black', lw=0.5,fill=False)
    ax.add_patch(rect)

    return gene_track_height+track_bottom+0.03  # return the y-axis value of gene track top

def SingleFragment(p1,p2,pair_dist):
    '''
    determine if mapped pair is single RNA fragment or not.
      1. no linker
      2. same chromosome
      3. same RNA
    '''
    if (p1.stop-p1.start==90) and (p2.stop-p2.start==100) and (p1.chr==p2.chr) and (p1.strand!=p2.strand):
        if abs(p1.center-p2.center)<pair_dist:
            return True

def PatchGen(i,j,h,step,bottom):
    ''' i,j are two loc, h is the unit hight of triangle, step is the resolution on x axis 
    bottom is the y-axis value for the bottom of heatmap  '''
    low_hight = ((j-i)/step-1)*h + bottom
    x = (i+j)*0.5
    if i!=j:
        p = Polygon([[x,low_hight],[x+0.5*step,low_hight+h],[x,low_hight+2*h],[x-0.5*step,low_hight+h]],"True",edgecolor='none')
    else:
        p = Polygon([[x+0.5*step,low_hight+h],[x,low_hight+2*h],[x-0.5*step,low_hight+h]],"True",edgecolor='none')
    return p


def Main():
    args=ParseArg()
    pair_dist=args.pair_dist
    step=args.step
    
    print "\nChecking if linkedPair file is tabixed..."
    if not os.path.isfile(args.linkedPair):
        print "LinkedPair file is not exist, please check!!"
        sys.exit(0)
    if not os.path.isfile(args.linkedPair+".tbi"):
        print "  tabix-ing..."
        os.system("sort -k1,1 -k2,2n "+args.linkedPair+" > temp_linkedPair.txt")
        os.system("bgzip temp_linkedPair.txt")
        os.system("tabix -p bed temp_linkedPair.txt.gz")
        linkedPair='temp_linkedPair.txt.gz'
    else:
        linkedPair=args.linkedPair
    print "  linkedPair file is tabixed."

    print "\nTabixing the interaction file..."
    os.system("sort -k1,1 -k2,2n "+args.interaction+" > temp_interaction.txt")
    os.system("bgzip temp_interaction.txt")
    os.system("tabix -p bed temp_interaction.txt.gz")
    print "  interaction file is tabixed."

    # start column number for second regions
    # s1 for interaction file and s2 for linkedPair file
    (s1,s2)=args.start    

    print "\nGet region information."
    if args.r:
        Region = read_region(args.r)
    elif args.name:
        os.system('grep "%s" %s > temp2.txt'%(args.name,args.genebed))
        g = open("temp2.txt").read().split('\t')
        if len(g)<2:
            print >> sys.stderr, "Error: the gene name is not found in database"
            sys.exit(0)
        s = int(g[1])
        e = int(g[2])
        Region = Bed([g[0],s-(e-s)/10,e+(e-s)/10,"region",".","."])
    else:
        print >> sys.stderr, "Error: Need to specify the region by '-r' or specify the gene name by '-n'"
        sys.exit(0)

    
    print "\n Start plot heatmaps on region: "+Region.str_region()
    fig = plt.figure(figsize=(8,6))
    ax = plt.subplot(111,frameon=False,yticks=[])
    start = Region.start
    end = Region.stop
    ax.set_xlim(start,end)

    #set x ticks withour offset
    locs=ax.get_xticks()
    ax.set_xticklabels(map(lambda x: "%i"%x, locs),fontsize=6)


    print "\nStart draw gene track"
    gene_dbi=DBI.init(args.genebed,"bed")
    print "  genebed indexed!"
    print "  Plot gene track"
    gene_top=Genetrack(Region,gene_dbi,ax,0.08)
    
    h = 1.5*step/(end-start) # unit height for triangles or polycons in heatmap

    

    print "\nQuery linkedPairs within specified region"
    os.system("tabix "+linkedPair+" %s:%i-%i > temp2.txt"%(Region.chr,Region.start,Region.stop))
    Count = {}
    for b in read_interaction("temp2.txt",s2):
        col='k'
        if args.Slim and SingleFragment(b[0],b[1],pair_dist): continue
        if Region.overlap(b[0],0) and Region.overlap(b[1],0): 
            if b[0].strand=='-':
                i = b[0].start
            else:
                i = b[0].stop
            if b[1].strand=='-':
                j = b[1].start
            else:
                j = b[1].stop       
            i = (i/step+1) * step  # approximate to the nearest central point
            j = (j/step+1) * step
            if i > j:
                temp=j
                j=i
                i=temp
            if (i,j) not in Count:
                Count[(i,j)] = 1
            else:
                Count[(i,j)] +=1
    
    print Count

    patches = []
    colors = []
    for i in range(start,end+1):
        if i%step!=0: continue
        for j in range(i,end+1):
            if j%step!=0 or (i,j) not in Count: continue
            patches.append(PatchGen(i,j,h,step,gene_top+0.01))
            colors.append(np.log(Count[(i,j)]+1))

    p = PatchCollection(patches, cmap=matplotlib.cm.Reds, alpha=0.7, edgecolor='k',linewidths=0.1)
    p.set_array(np.array(colors))
    ax.add_collection(p)

    ax.set_ylim(0,((end-start)/step+2)*h+gene_top+0.01)
    plt.colorbar(p)

    if not args.SI:
        plt.savefig(args.output)
        plt.show()
        os.system("rm temp_interaction.txt.gz*")
        if not os.path.isfile(args.linkedPair+".tbi"):
            os.system("rm temp_linkedPair.txt.gz*")
        os.system("rm temp2.txt")
        sys.exit(0)

    print "\nQuery interactions"
    os.system("tabix temp_interaction.txt.gz %s:%i-%i > temp2.txt"%(Region.chr,Region.start,Region.stop))
    print "\nList of interactions plotted: "
    k=1
    cmap=cm.get_cmap('Paired', 10)
    cmap=cmap(range(10))
    bottom = gene_top+0.01
    for b in read_interaction("temp2.txt",s1):
        if b[0].overlap(b[1],0): continue
        if Region.overlap(b[1],0):
            k+=1
            if b[1].stop > b[0].stop:
                start1 = b[0].start
                end1 = b[0].stop
                start2 = b[1].start
                end2 = b[1].stop
            else:
                start1 = b[1].start
                end1 = b[1].stop
                start2 = b[0].start
                end2 = b[0].stop
            P1=Polygon([[start1,bottom],[end1,bottom],[(end1+end2)*0.5,(end2-end1)*h/step+bottom],[(start1+end2)*0.5,(end2-start1)*h/step+bottom]],"True",facecolor='none',edgecolor=cmap[k%10],alpha=0.4,lw=0.5)
            P2=Polygon([[start2,bottom],[end2,bottom],[(start1+end2)*0.5,(end2-start1)*h/step+bottom],[(start1+start2)*0.5,(start2-start1)*h/step+bottom]],"True",facecolor='none',edgecolor=cmap[k%10],alpha=0.4,lw=0.5)
            ax.add_patch(P1)
            ax.add_patch(P2)
            print "  "+b[0].str_region()+" <-> "+b[1].str_region()


    plt.savefig(args.output)
    plt.show()

    # remove temp file
    os.system("rm temp_interaction.txt.gz*")
    if not os.path.isfile(args.linkedPair+".tbi"):
        os.system("rm temp_linkedPair.txt.gz*")
    os.system("rm temp2.txt")

if __name__=="__main__":
    Main() 
