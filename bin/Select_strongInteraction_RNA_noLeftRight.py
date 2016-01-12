import sys,argparse
from time import time
import copy,math
from data_structure import *
from scipy.stats import hypergeom
from xplib import DBI
from xplib.Annotation import Bed
from Annotation import *

def ParseArg():
    p=argparse.ArgumentParser(description="find strong interactions between RNAs from annotation of RNAs, intron regions are removed, no self interactions",epilog="need Scipy for hypergeometric distribution")
    p.add_argument("-i","--input",type=str,required=True,help="input file which is the output file of Stitch-seq-Aligner.py")
    p.add_argument("-M","--min_clusterS",type=int,default=5,help="minimum number of segments allowed in each cluster, default:5")
    p.add_argument("-m","--min_interaction",type=int,default=3,help="minimum number of interactions to support a strong interaction, default:3")
    p.add_argument('-p',"--p_value",type=float,default=0.05,help="the p-value based on hypergeometric distribution to call strong interactions, default: 0.05")
    p.add_argument('-o','--output',type=str,help="specify output file")
    p.add_argument('-O','--output_intra', type=str, help="specify output file for intra interactions")
    p.add_argument('-a','--annotation',type=str, help="Annotation file for all RNA components")
    p.add_argument('-r','--annotation_repeat',type=str, help="Annotation file for all repeat components")
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
    if (p1.end-p1.start==89) and (p2.end-p2.start==98) and (p1.chr==p2.chr) and p1.name==p2.name:
        return True

def GetAnnotationName(pAnno, hasAnno, dbi, hasRepeat, dbirepeat):
    if pAnno.source == "genome":
        if hasAnno:
            bed = Bed([pAnno.chr, pAnno.start, pAnno.end, '.', 0.0, '+'])
            for hit in dbi.query(bed):
                name_component = hit.id.split(".", 2)[2]
                id_component = hit.id.split(".", 2)[1]
                if pAnno.name == name_component or pAnno.name.split("_", 2)[0] == id_component:
                    pAnno.id = id_component
                    return ":".join(str(f) for f in [pAnno.type, pAnno.id, pAnno.proper])
                    
        if hasRepeat:
            for hit in dbirepeat.query(bed):
                tempname=hit.id.split("&")
                if pAnno.name == tempname[0]:
                    pAnno.id = "".join(str(f) for f in [pAnno.chr, pAnno.name, hit.start])
                    return ":".join(str(f) for f in [pAnno.type, pAnno.id, pAnno.proper])
        if pAnno.type == 'tRNA':
            return ":".join(str(f) for f in [pAnno.type, pAnno.name, pAnno.proper])
        if hasAnno or hasRepeat:
            raise Exception('pAnno not found! ' + pAnno.chr + ':' + pAnno.name)
        pAnno.id = "".join(str(f) for f in [pAnno.type, pAnno.chr, pAnno.name, pAnno.start])
    elif ".fa" in pAnno.source:
        pAnno.id = pAnno.chr.split("_")[0]
        try:
            pAnno.chr = pAnno.chr.split("_")[1]
        except:
            pass
    else:
        pAnno.id = pAnno.name
    return pAnno.type + ":" + pAnno.id + ":" + pAnno.proper 

def Main():
    t1=time()
    args=ParseArg()
    inp = open(args.input, 'r')
    min_clusterS=args.min_clusterS
    min_interaction=args.min_interaction
    p_value=args.p_value
    output=open(args.output,'w')
    outputIntra = open(args.output_intra, 'w')

    hasAnnotation = False
    if args.annotation:
        dbi = DBI.init(args.annotation, "bed")
        hasAnnotation = True
    else:
        dbi = False

    if args.annotation_repeat:
        dbirepeat = DBI.init(args.annotation_repeat, "bed")
        hasAnnotationRepeat = True
    else:
        dbirepeat = False        

    #store count of RNA for part1 and part2
    part={}


    k=0
    sgcount = 0 #single fragment count
    
    print >> sys.stderr,"# Inputing data..."
    interaction = {}  # store number of interactions for different RNA
    selfinteraction = {}



    #Types = ["snoRNA","protein_coding","snRNA","lincRNA","tRNA","misc_RNA","pseudogene","miRNA","antisense","sense_intronic","non_coding","processed_transcript","sense_overlapping","rRNA_repeat","rRNA"]
    for line in inp.read().split('\n'):
        if line=='': continue
        line=line.strip().split('\t')
        p1=annotated_bed_proper(line[0:10],id=k,cluster=1)
        p2=annotated_bed_proper(line[11:],id=k,cluster=1)
        if isinstance(p1.start, list):
            p1.start=int(p1.start[0])
            p1.end=int(p1.end[-1])
        if isinstance(p2.start, list):
            p2.start=int(p2.start[0])
            p2.end=int(p2.end[-1])
                
        if SingleFragment(p1,p2):
            sgcount += 1
            continue
        k+=1
        #if p1.subtype=="intron" or p2.subtype=="intron": continue
        #if p1.type in Types:
        try:
            p1_name = GetAnnotationName(p1, hasAnnotation, dbi, hasAnnotationRepeat, dbirepeat) 
            if p1_name not in part:
                part[p1_name]=1
            else:
                part[p1_name]+=1  
            #if p2.type in Types:
            p2_name = GetAnnotationName(p2, hasAnnotation, dbi, hasAnnotationRepeat, dbirepeat) 
            if not p1_name == p2_name: # count once for self-interaction
                if p2_name not in part:
                    part[p2_name]=1
                else:
                    part[p2_name]+=1
            #if p1.type in Types and p2.type in Types:
            if p1_name == p2_name:
                if p1_name not in selfinteraction:
                    selfinteraction[p1_name]=copy.deepcopy(p1)
                else:
                    selfinteraction[p1_name].Update(p1.start, p1.end)
                    selfinteraction[p1_name].Update(p2.start, p2.end)
                    selfinteraction[p1_name].cluster += 1
            else:
                if p1_name>p2_name:
                    temp = p1
                    p1 = p2
                    p2 = temp
                    tempName = p1_name
                    p1_name = p2_name
                    p2_name = tempName
                inter_name = p1_name + "--" + p2_name
                if inter_name not in interaction:
                    interaction[inter_name]=[copy.deepcopy(p1),copy.deepcopy(p2)]
                else:
                    interaction[inter_name][0].Update(p1.start,p1.end)
                    interaction[inter_name][1].Update(p2.start,p2.end)
                    interaction[inter_name][0].cluster+=1
        except Exception as e:
            print >> sys.stderr, e
        if k%20000==0: 
            print >> sys.stderr,"  Reading %d pairs of segments\r"%(k),
    print >> sys.stdout,"Get total %d pairs."%(k)
    print >> sys.stdout,"Single fragment count: %d."%(sgcount)

    print >>sys.stdout,"   number of different RNAs is %d          "%(len(part))
    
    total = k # total pairs used
    n=0
    k=0  # record number of strong interactions
    for i in interaction:
        n+=1
        count = interaction[i][0].cluster
        if count < min_interaction: continue
        p1_name = i.split("--")[0]
        p2_name = i.split("--")[1]
        P1 = interaction[i][0]
        P2 = interaction[i][1]
        P1.cluster = part[p1_name]
        P2.cluster = part[p2_name]
        if part[p1_name]<min_clusterS or part[p2_name]<min_clusterS: continue
        real_p=1-hypergeom.cdf(count,total,part[p1_name],part[p2_name])
        if real_p<=p_value:
            k=k+1
            try:
                log_p = math.log(real_p)
            except:
                log_p = -float("Inf")
            print >> output, str(P1)+'\t'+str(P2)+'\t%d\t%.4f'%(count,log_p)
        if n%500==0: print >> sys.stderr, "  Progress ( %d / %d )\r"%(n,len(interaction)),
    k1=0
    for i in selfinteraction:
        n+=1
        count = selfinteraction[i].cluster
        if count < min_interaction: continue
        p1_name = i
        P1 = selfinteraction[i]
        P1.cluster = part[p1_name]
        if part[p1_name]<min_clusterS: continue
        k1=k1+1
        print >> outputIntra, str(P1)+'\t%d'%(count)
        if n%500==0: print >> sys.stderr, "  Progress ( %d / %d )\r"%(n,len(interaction)),
    print >> sys.stdout,"# Find %d strong and %d self interactions. Cost time: %.2f s"%(k,k1,time()-t1)

if __name__=="__main__":
    Main()
