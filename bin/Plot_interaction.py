import sys,argparse,os
from xplib import DBI
from bx.bbi.bigwig_file import BigWigFile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm
import numpy as np
plt.ioff()

def ParseArg():
    p=argparse.ArgumentParser(description="plot linked pairs around a given interaction. information of linked pairs are stored in file '*_fragment_paired_align.txt'",epilog="Require: matplotlib, numpy")
    p.add_argument("interaction",type=str,help="Interaction file from output of 'Select_strongInteraction_pp.py'")
    p.add_argument("linkedPair",type=str,help="file for information of linked pairs, which is output of 'Stitch-seq_Aligner.py'")
    p.add_argument('-n',type=int,help="Choose region to plot, it can be a number (around n-th interaction in the interaction file). This is mutually exclusive with '-r' option")
    p.add_argument('-r',type=str,nargs='+',help="Choose region to plot, give two interaction regions with format 'chr:start-end', this is mutually exclusive with '-n' option")
    p.add_argument('-s','--start',type=int,nargs='+',default=(7,11),help='start column number of the second region in interaction file and linkedPair file, default=(7,11)')
    p.add_argument('-d','--distance',type=int,default=10,help='the plus-minus distance (unit: kbp) flanking the interaction regions to be plotted, default=10')
    p.add_argument('-g','--genebed',type=str,default='/home/yu68/bharat-interaction/new_lincRNA_data/Ensembl_mm9.genebed',help='the genebed file from Ensembl, default: Ensembl_mm9.genebed')
    p.add_argument('-w','--phyloP_wig',type=str,default='/data2/sysbio/UCSD-sequencing/mouse.phyloP30way.bw',help='the bigWig file for phyloP scores,defualt: mouse.phyloP30way.bw')
    p.add_argument("-p","--pair_dist",type=int,default=1000,help="two interacted parts within this distance are considered as self-ligated and they are marked or eliminated (see option -s for slim mode), default: 1000bp")
    p.add_argument("-S","--Slim",action='store_true',help='set slim mode to eliminate self ligated interactions')
    p.add_argument('-o','--output',type=str,help="output plot file, can be format of emf, eps, pdf, png, ps, raw, rgba, svg, svgz")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        sys.exit(0)
    return p.parse_args()

class Bed:
    def __init__(self,x,**kwargs):
        self.chr=x[0]
        try:
            self.start=int(x[1])
            self.stop=int(x[2])
        except:
            self.start=int(x[1].split(",")[0])
            self.stop=int(x[2].split(",")[1])
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
        print lsep
        if lsep[3] in ['+','-']:
            bed1=Bed(lsep[0:3],strand=lsep[3])
            bed2=Bed(lsep[s:(s+3)],strand=lsep[s+3])
        else:
            bed1=Bed(lsep[0:3])
            bed2=Bed(lsep[s:(s+3)])
        yield (bed1,bed2,l)


def transform(x1,start1,end1,start2,end2):
    '''
    transfrom from scale 2 to scale 1
    '''
    return int(1.0*(x1-start1)/(end1-start1)*(end2-start2)+start2)


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
    if (p1.stop-p1.start>=78) and (p2.stop-p2.start>=88) and (p1.chr==p2.chr):
        if abs(p1.center-p2.center)<pair_dist:
            return True


def Wigtrack(bed, bw, ax, track_bottom,col):
    array=bw.summarize(bed.chr,bed.start,bed.stop,(bed.stop-bed.start)/10).sum_data
    array=np.array(array)
    Min=min(array)
    Max=max(array)
    array_n=(array-Min)*0.2/(Max-Min)+track_bottom+0.01
    ax.plot(range(bed.start+5,bed.stop-4,10),array_n,color=col,lw=0.5)
    rect=matplotlib.patches.Rectangle((bed.start,track_bottom),bed.stop-bed.start, 0.22, edgecolor='black', lw=0.5,fill=False)
    ax.add_patch(rect)
    ax.text(bed.start, track_bottom+0.16, "phyloP",horizontalalignment='left',fontsize=8,color=col)
    return track_bottom+0.22
    
    


def Main():
    args=ParseArg()
    distance=args.distance*1000
    pair_dist=args.pair_dist
    
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

    print "\nExtracting interaction information..."
    if args.n:
        Interactions=open(args.interaction,'r')
        l=Interactions.read().split('\n')[args.n-1].split('\t')
        part1=Bed(l[0:6])
        part2=Bed(l[s1:(s1+6)])
    elif len(args.r)==2:
        part1=read_region(args.r[0])
        part2=read_region(args.r[1])
    else:
        print >> sys.stderr, "need to specify two regions using '-r'"
    
    if "chr" not in part1.chr or "chr" not in part2.chr:
        print >> sys.stderr, "This program only works for genomic regions."
        exit(0)    
    
    start1=part1.start-distance
    end1=part1.stop+distance
    start2=part2.start-distance
    end2=part2.stop+distance
    # if the searched regions for part1 and part2 are overlapped, using the same regions for both part
    if part1.overlap(part2,-2*distance):
        start1=min(start1,start2)
        start2=min(start1,start2)
        end1=max(end1,end2)
        end2=max(end1,end2)

    
    # initialize figure
    print "\n Start plot interaction: "+part1.str_region()+" <-> "+part2.str_region()
    col1="#4F81BD"
    col2="#C0504D"
    fig = plt.figure(figsize=(8,4))
    ax1 = plt.subplot(111,frameon=False,yticks=[])
    plt.tick_params(axis="y",which="both",left="off",right="off",labelleft="off")  # remove y ticks
    plt.subplots_adjust(top=0.75)
    ax2 = ax1.twiny()
    ax1.set_xlim(start1,end1)
    ax2.set_xlim(start2,end2)
    ax1.set_ylim(0,1)
    ax2.set_ylim(0,1)

    #set x ticks withour offset
    locs=ax1.get_xticks()
    ax1.set_xticklabels(map(lambda x: "%i"%x, locs),fontsize=8)
    locs=ax2.get_xticks()
    ax2.set_xticklabels(map(lambda x: "%i"%x, locs),fontsize=8)
    
    
    # input bigWig file for phyloP score
    bw_phyloP = BigWigFile(open(args.phyloP_wig))

    print "\nStart draw gene track"
    gene_dbi=DBI.init(args.genebed,"bed")
    print "  genebed indexed!"
    print "  Plot gene track for Part1"
    gene1_top=Genetrack(Bed([part1.chr,start1,end1]),gene_dbi,ax1,0.08)
    wig1_top=Wigtrack(Bed([part1.chr,start1,end1]), bw_phyloP, ax1, gene1_top,col1)
    
    y_1=wig1_top+0.1
    y_2=y_1+0.2
    
    print "  Plot gene track for Part2"
    gene2_top=Genetrack(Bed([part2.chr,start2,end2]),gene_dbi,ax2,y_2+0.08)    
    wig2_top=Wigtrack(Bed([part2.chr,start2,end2]), bw_phyloP, ax2, gene2_top,col2)
    print "\nQuery interactions within +-%dkbp of interaction"%(distance/1000)
    os.system("tabix temp_interaction.txt.gz %s:%i-%i > temp2.txt"%(part1.chr,start1,end1))
    print "\nList of interactions plotted: "
    k=1
    cmap=cm.get_cmap('Paired', 10)
    cmap=cmap(range(10))
    for b in read_interaction("temp2.txt",s1):
        #if args.Slim and b[0].overlap(b[1],-pair_dist): continue
        if Bed([part2.chr,start2,end2]).overlap(b[1],0):
            k+=1
            x1_2_start=transform(b[0].start,start1,end1,start2,end2)
            x1_2_end=transform(b[0].stop,start1,end1,start2,end2)
            ax2.add_patch(matplotlib.patches.Polygon([[x1_2_start,y_1+0.04],[x1_2_end,y_1+0.04],[b[1].stop,y_2],[b[1].start,y_2]],color=cmap[k%10],alpha=0.4,lw=0.5))
            ax1.add_patch(matplotlib.patches.Rectangle((b[0].start,y_1),b[0].stop-b[0].start,0.04,color=col1,lw=0.5))
            ax2.add_patch(matplotlib.patches.Rectangle((b[1].start,y_2),b[1].stop-b[1].start,0.04,color=col2,lw=0.5))
            print "  "+b[0].str_region()+" <-> "+b[1].str_region()


    ax1.plot([start1,end1],[y_1+0.02,y_1+0.02],color=col1,linewidth=1,alpha=0.7)
    ax2.plot([start2,end2],[y_2+0.02,y_2+0.02],color=col2,linewidth=1,alpha=0.7)


    print "\nQuery linkedPairs within +-%dkbp of interaction"%(distance/1000)
    os.system("tabix "+linkedPair+" %s:%i-%i > temp2.txt"%(part1.chr,start1,end1))
    print "\nList of linked pairs plotted: "
    for b in read_interaction("temp2.txt",s2):
        col='k'
        if args.Slim and SingleFragment(b[0],b[1],pair_dist): continue
        if SingleFragment(b[0],b[1],pair_dist): col='#03C03C'
        if part1.overlap(b[0],-distance) and part2.overlap(b[1],-distance):    
            x1_2_start=transform(b[0].start,start1,end1,start2,end2)
            x1_2_end=transform(b[0].stop,start1,end1,start2,end2)
            if b[0].strand=='-':
                connect1=x1_2_start
            else:
                connect1=x1_2_end
            if b[1].strand=="-":
                connect2=b[1].start
            else:
                connect2=b[1].stop
            ax2.plot([connect1,connect2],[y_1+0.02,y_2+0.02],color=col,alpha=0.3,lw=0.5)
            ax1.plot([b[0].start,b[0].stop],[y_1+0.02,y_1+0.02],color=col,alpha=0.3,lw=0.8)
            ax2.plot([b[1].start,b[1].stop],[y_2+0.02,y_2+0.02],color=col,alpha=0.3,lw=0.8)
           # print "  "+b[0].str_region()+" <-> "+b[1].str_region()
    plt.text(0.5, 1.15, part1.str_region()+" <-> "+part2.str_region(),
         horizontalalignment='center',
         fontsize=10,
         transform = ax1.transAxes)
    plt.text(0.5, 1.10, "Distance: +-%dkbp of interaction"%(distance/1000),
         horizontalalignment='center',
         fontsize=8,
         transform = ax1.transAxes)
    ax1.text(part1.center,y_1-0.03,"|".join([part1.type,part1.name,part1.subtype]),
             verticalalignment='center', horizontalalignment='center',fontsize=8,color=col1)
    ax2.text(part2.center,y_2+0.07,"|".join([part2.type,part2.name,part2.subtype]),
             verticalalignment='center', horizontalalignment='center',fontsize=8,color=col2)
    ax1.set_ylim(0,wig2_top+0.1)
    ax2.set_ylim(0,wig2_top+0.1)
    ax1.text(start1, 0.05, part1.chr,horizontalalignment='left',fontsize=8)
    ax2.text(start2, wig2_top+0.04, part2.chr,horizontalalignment='left',fontsize=8)
    plt.savefig(args.output)
    plt.show()
     
    # remove temp file
    os.system("rm temp_interaction.txt.gz*")
    if not os.path.isfile(args.linkedPair+".tbi"):
        os.system("rm temp_linkedPair.txt.gz*")
    os.system("rm temp2.txt")

if __name__=="__main__":
    Main()
