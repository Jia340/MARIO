import sys, os, argparse, gc, json, copy
import pysam
import itertools,string
from Bio import SeqIO
from xplib import TableIO
from xplib import DBI
from xplib.Annotation import Bed
from AnnoMax import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

#Annotation for miRNAs is from the .blatresult itself. The sequence is from the fastq file
#Annotation for transcripts is the reference, that is, rna.fa
#Annotation for genome is from xplib

def ParseArg():
    p=argparse.ArgumentParser( description = 'Annotate aligned RNAs. ')
    p.add_argument('input_path',type=str,help="path for the input folder, which is the output folder of Stitch-seq_Aligner.py")
    p.add_argument('file_perfix',type=str,help="perfix for the file. The program will search files in the input folder based on this argument")
    p.add_argument('--annotype',type=str,metavar='at',default="genome",help="Reference types: (miRNA, genome, transcript, other). The order should be the same as that of the bam/blatresult files",nargs="*")
    p.add_argument('--tx_anno',type=str,metavar='txa',help="Annotation file for transcriptome. It is the reference for alignment, that is, rna.fa. Should be specified if annotype contains transcript")
    p.add_argument('-a','--annotation',type=str,help='If specified, include the RNA type annotation for each aligned pair, need to give bed annotation RNA file. Should be specified if annotype contains genome')
    p.add_argument("-A","--annotationGenebed",dest="db_detail",type=str,help="annotation bed12 file for lincRNA and mRNA with intron and exon. Should be specified if annotype contains genome")
    p.add_argument("-R","--annotationRepeat",dest="db_repeat",type=str,help="annotation bed6 file from repeatMasker. Should be specificed if annotype contains genome")
    p.add_argument("-M","--mapq",type=int,default=2,help="For reads mapped to genome, reads with MAPQ lower than M will be labeled as 'Many' and will be discarded in next step")
    p.add_argument("-P","--pair",type=int, help="Specify if the RNA is RNA1/2. Should be 1 or 2")
    if len(sys.argv)==1:
        #print (p.print_help(),file=sys.stderr)
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()



class ensemblSeq:

    def getSubType(self, start, end):
        # notice that this is relative to the RNA sequence
        # therefore, just compare directly with self.utrlength5 or self.utrlength3 should be enough

        # return format: [type, name, subtype]
        # first find the midpoint of the region
        midpoint = (start + end) / 2
        if midpoint < self.utrlength5:
            return "utr5"
        elif midpoint >= self.length - self.utrlength3:
            return "utr3"
        else:
            if self.biotype == "lincRNA" or self.biotype == "protein_coding":
                return "exon"
            elif "intron" in self.biotype:
                return "intron"
            else:
                return "."

    def getAnnotation(self, start, end):
        return [self.biotype, self.genename, self.getSubType(start, end)]

    #def getGeneID(self):
    #    return self.geneID
    def __init__(self, name, seq):
        # notice that the ">" in name is chopped
        self.longname = name.strip()
        tokens = name.split("|")
        self.transID = tokens[0]
        self.strand = (int(tokens[3]) > 0)
        self.tss = int(tokens[1]) if self.strand else int(tokens[2])
        self.tes = int(tokens[2]) if self.strand else int(tokens[1])
        self.genename = tokens[8]
        if not self.genename:
            self.genename = self.transID
        #self.geneID = tokens[9]
        #if not self.geneID:
        #    self.geneID = self.transID
        self.biotype = tokens[9]

        utrstarts5 = (tokens[4] if self.strand else tokens[5]).split(";")
        utrends5 = (tokens[5] if self.strand else tokens[4]).split(";")
        utrstarts3 = (tokens[6] if self.strand else tokens[7]).split(";")
        utrends3 = (tokens[7] if self.strand else tokens[6]).split(";")

        self.utrlength5 = 0
        self.utrlength3 = 0
        if tokens[4]:
            for i in xrange(len(utrstarts5)):
                self.utrlength5 += abs(int(utrends5[i]) - int(utrstarts5[i]))

        if tokens[6]:
            for i in xrange(len(utrstarts3)):
                self.utrlength3 += abs(int(utrends3[i]) - int(utrstarts3[i]))

        self.seq = seq.strip()
        self.length = len(self.seq)
        
def get_readdict(readfilename):
    seqdict = dict()

    fseq = open(readfilename, 'r')
    seqname = ''
    oldseq = ''

    for line in fseq:
        if line.startswith('>') or line.startswith('@'):
            if oldseq and seqname:
                seqdict[seqname] = oldseq
            seqname = line.strip('>@').split(" ")[0].strip()
            oldseq = ''
        elif line.startswith('+'):
            # fastq
            if oldseq and seqname:
                seqdict[seqname] = oldseq
            seqname = ''
            oldseq = ''
        elif seqname:
            oldseq += line.strip()

    if oldseq and seqname:
        seqdict[seqname] = oldseq

    fseq.close()
    return seqdict

def blat_annotation(outputfilename, typename, readfilename, requireUnique = False, posstrand = True, strandenforced = False, mismatchthr = 2, results_dict = dict()):
    # results_dict is the dictionary for all previous results, can be none
    # key is the read name, value is the annotation string
    # this function put all qualified match in outputfilename
    # then put all unqualified match together with unmappedfile
    # return value is the new annotation dictionary
    newdict = dict()
    seqdict = get_readdict(readfilename)
   
    multiset = set()

    fres = open(outputfilename, 'r')
    # first populate the dictionary
    # read the header
    for x in xrange(5):
        fres.readline()

    for line in fres:
        tokens = line.split()
        strand = (tokens[8] == '+')
        if (not strandenforced) or strand == posstrand:
            # correct strand or strand not enforced
            if strand != posstrand:
                # incorrect strand, but not enforced
                # add the 21th column saying that strand is wrong
                strandcol = 'NonProperStrand'
            else:
                strandcol = 'ProperStrand'

            match = int(tokens[0])
            mismatch = int(tokens[1])
            readname = tokens[9]
            readlen = int(tokens[10])
            readstart = int(tokens[11])
            readend = int(tokens[12])

            mismatch += (readlen - readend) + readstart
            if mismatch <= mismatchthr:
                # this one is a match
                if readname in seqdict:
                    # first time happen
                    if True:
                        curr_anno_arr = [tokens[13], tokens[15], tokens[16], tokens[8], seqdict[readname], typename, typename, tokens[13], '.']
                        if not strandenforced:
                            curr_anno_arr.append(strandcol)
                        curr_anno = '\t'.join(curr_anno_arr)
       #             else:
       #                 curr_anno_arr = [tokens[13], tokens[15], tokens[16], tokens[8], seqdict[readname], typename]
       #                 if not strandenforced:
       #                     curr_anno_arr.append(strandcol)
       #                 curr_anno = '\t'.join(curr_anno_arr)
                    newdict[readname] = curr_anno
                    del seqdict[readname]
                elif (readname in newdict) or (readname in multiset):
                    # this is a multimatch
                    # and not first time
                    multiset.add(readname)
                    if requireUnique:
                        # then this multi-match will be discarded
                        del newdict[readname]
                    else:
                        if type(newdict[readname]) is str:
                            newdict[readname] = [newdict[readname]]
                        seq = newdict[readname][0].split()[4]
                        if True:
                            curr_anno_arr = [tokens[13], tokens[15], tokens[16], tokens[8], seq, typename, typename, tokens[13], '.', strandcol]
                            newdict[readname].append('\t'.join(curr_anno_arr))
        #                else:
        #                    curr_anno_arr = [tokens[13], tokens[15], tokens[16], tokens[8], seq, typename, strandcol]
        #                    newdict[readname].append('\t'.join(curr_anno_arr))

    fres.close()

    # merge annotation
    newanno = dict(results_dict.items() + newdict.items())

    return newanno

def otherlib_annotation(outputbam, annotationfile = 'misc', requireUnique = False, posstrand = True, strandenforced = False, results_dict = dict()):

    # outputbam is a pysam file
    # annotationfile can be:
    # "misc", "guess" from RNA name
    # or a known annotation file (for example, rna.fa from NCBI)
    

    newdict = dict()
#    funmap = open(unmapfilename, 'w')

    if annotationfile != 'misc':
        refdic = dict()
        if not 'fa' in annotationfile:
            annotationfile += ".fa"

        # this will be a dictionary of ensembl entries
        fref = open(annotationfile, 'r')
        refname = ''
        refseq = ''
        refkey = ''

        for line in fref:
            if line.startswith(">"):
                if refkey:
                    # refkey: ensembl id. there is an old ref there
                    if refkey in refdic:
                        print >> sys.stderr, refkey
                    refdic[refkey] = ensemblSeq(refname, refseq)
                    refseq = ''
                    refkey = ''
                    refname = ''
                refkey = line.strip(">").split("|")[0]
                refname = line.strip(">").strip()
            else:
                refseq += line.strip()

        if refkey:
            if refkey in refdic:
                print >> sys.stderr, refkey
            refdic[refkey] = ensemblSeq(refname, refseq)
        fref.close()

    for record in outputbam:
        name = record.qname.split(" ")[0]
        unique = False
        ##transcripts should be merged. Otherwise it's hard to know whether it's multimapped. Many/One won't be used for transcripts in link_frag.py


        if ((not requireUnique) or unique) and ((not strandenforced) or record.is_reverse != posstrand) and not record.is_unmapped:
            # this is a suitable entry
            strand = "+"
            if record.is_reverse:
                strand = "-"
            if record.is_reverse == posstrand:
                # strand is wrong
                strandcol = 'NonProperStrand'
            else:
                strandcol = 'ProperStrand'
            if True:
                if annotationfile == 'misc':
                    rnaname = outputbam.getrname(record.tid).split(" ")[0]
                    rnaid = rnaname
                    if rnaname.startswith("NONMMUT"):
                        # This is a lncRNA from NONCODE
                        rnatype = "ncRNA"
                        rnasubtype = "."
                    else:
                        # tRNA
                        rnatype = "tRNA"
                        rnasubtype = "."
                else:
                    rnaensembl = refdic[outputbam.getrname(record.tid).split("|")[0]]
                    rnaid = rnaensembl.transID
                    #geneid = rnaensembl.getGeneID()
                    [rnatype, rnaname, rnasubtype] = rnaensembl.getAnnotation(record.aend - record.alen + 1, record.aend)

                curr_anno_arr = (str(f) for f in [rnaid, record.aend - record.alen + 1, record.aend, strand, record.seq, annotationfile.split("/")[-1], rnatype, rnaname, rnasubtype, strandcol])
                if not name in newdict:
                    newdict[name] = '\t'.join(curr_anno_arr)
                else:
                    if type(newdict[name]) is str:
                        newdict[name] = [newdict[name]]
                    newdict[name].append('\t'.join(curr_anno_arr))

            else:
                curr_anno_arr = (str(f) for f in [rnaid, record.aend - record.alen + 1, record.aend, strand, record.seq, annotationfile.split("/")[-1], strandcol])
                if not name in newdict:
                    newdict[name] = '\t'.join(curr_anno_arr)
                else:
                    if type(newdict[name]) is str:
                        newdict[name] = [newdict[name]]
                    newdict[name].append('\t'.join(curr_anno_arr))


    newanno = dict(results_dict.items() + newdict.items())
    return newanno

def Included(record,RequireUnique,mapq_thred):
    # record is a pysam read
    # non-unique alignment in Bowtie2 has 'XS' tag: https://www.biostars.org/p/18965/
    if RequireUnique:
        if record.mapq >= mapq_thred:
            unique=True
        else:
            unique=False
    else:
        unique=True # not consider unique
    return (not record.is_unmapped)&unique

def Exon_junction(record):
    # Parse cigar
    # return start and end of each part
    b_list=[]
    S=record.pos
    E=S
    a_s=0
    a_e=0
    max_l=0
    for item in record.cigar:
        if item[0]==1:
            pass
        elif item[0]==3:
            if E-S+1>max_l:
                max_l=E-S+1
                a_s=S
                a_e=E
            b_list.append((S,E))
            E+=int(item[1])
            S=E
        else:
            E+=int(item[1])
    b_list.append((S,E))           
    if E-S+1>max_l:
        max_l=E-S+1
        a_s=S
        a_e=E
    return b_list,a_s,a_e
              


def genome_annotation(outputbam, annotationfile, detail, annotationRepeat, mapq_thred ,strandenforced = False, posstrand = True, requireUnique = False, results_dict = dict()):
    # annotationfile is annotation file
    # detail is db_detail file

    if annotationfile:
        dbi1=DBI.init(annotationfile,"bed")
        dbi2=DBI.init(detail,"bed")
        dbi3=DBI.init(annotationRepeat,"bed")

    newdict = dict()
#    funmap = open(unmapfilename, 'w')

    for record in outputbam:
        # print >> sys.stderr, record.qname
        if "N" not in record.cigarstring:
            anno_start=record.pos
            anno_end=record.aend
            bed_start=record.pos
            bed_end=record.aend
        else:
            bed_list,anno_start,anno_end = Exon_junction(record)
            bed_start=",".join([str(f[0]) for f in bed_list])
            bed_end=",".join([str(f[1]) for f in bed_list])
#        print anno_start,anno_end,bed_start,bed_end
        IsMapped = False

        if Included(record, requireUnique,mapq_thred):
            strandactual = ("+" if posstrand else "-")
            strand = "+"
            if record.is_reverse:
                strandactual = ("-" if posstrand else "+")
                strand = "-"
            if annotationfile:
                bed=Bed([outputbam.getrname(record.tid), anno_start, anno_end,'.',0.0,strandactual])
                [typ, name, subtype, strandcol] = annotation(bed,dbi1,dbi2,dbi3)
                if (not strandenforced) or strandcol == 'ProperStrand':
                    curr_anno_arr = (str(f) for f in [outputbam.getrname(record.tid), bed_start, bed_end, strand, record.seq, 'genome', typ, name, subtype, strandcol])
                    if not record.qname in newdict:
                        newdict[record.qname] = '\t'.join(curr_anno_arr)
                        if not Included(record, True, mapq_thred):
                            # not unique
                            newdict[record.qname] = [newdict[record.qname]]
                    else:
                        if type(newdict[record.qname]) is str:
                            newdict[record.qname] = [newdict[record.qname]]
                        newdict[record.qname].append('\t'.join(curr_anno_arr))
                    IsMapped = True
            else:
                strandcol = '.'
                curr_anno_arr = (str(f) for f in [outputbam.getrname(record.tid), record.aend - record.alen + 1, record.aend, strand, record.seq, 'genome', strandcol])
                if not record.qname in newdict:
                    newdict[record.qname] = '\t'.join(curr_anno_arr)
                    if not Included(record, True, mapq_thred):
                        # not unique
                        newdict[record.qname] = [newdict[record.qname]]
                else:
                    if type(newdict[record.qname]) is str:
                        newdict[record.qname] = [newdict[record.qname]]
                    newdict[record.qname].append('\t'.join(curr_anno_arr))
                IsMapped = True

    newanno = dict(results_dict.items() + newdict.items())
    return newanno

def Main():
    args=ParseArg()
    fin=args.input_path.rstrip('/')
    fperfix=args.file_perfix
    annotation = args.annotation
    db_detail = args.db_detail
    db_repeat = args.db_repeat
    annofile = args.tx_anno
    
    del args.annotation
    del args.db_detail
    del args.db_repeat

    strandenforced = False    
    
    unique=False
    
    if args.pair==1:
        strand=True
    elif args.pair==2:
        strand=False
    else:
        print >>sys.stderr,"Incorrect -P"

    annodictentry = dict()
    reftypelist=args.annotype  
    print >>sys.stderr,"Annotation type:%s"%('\t'.join(reftypelist))
    inputbam_name=fin+"/"+"sort_"+fperfix+".bam"
    for i in xrange(len(reftypelist)):
      gc.collect()
      reftype=reftypelist[i]
      if reftype.lower()=="mirna":
          inputbam=fin+'/'+fperfix+".blatresult"
          annodictentry = blat_annotation(inputbam, reftype, fin+'/'+fperfix+'_blat.fasta', unique, strand, strandenforced, 2, annodictentry)
      else:
          inputbam_name=inputbam_name[0:-4]+"_unmap"+".bam"
          inputbam=pysam.Samfile(inputbam_name,"rb")
          if reftype.lower()=="genome":
              annodictentry = genome_annotation(inputbam, annotation, db_detail, db_repeat, args.mapq, strandenforced, strand, unique, annodictentry)
          elif reftype.lower()=="transcript":
              annodictentry = otherlib_annotation(inputbam, annofile, unique, strand, strandenforced, annodictentry)
          inputbam.close()
    gc.collect()

    fReadOut=open(fperfix+".pairOutput","w")
    for entry, value in annodictentry.iteritems():
        if type(value) is str:
            fReadOut.write('\t'.join(str(f) for f in [value, entry, 'One']) + os.linesep)
        else:
            for item in value:
                fReadOut.write('\t'.join(str(f) for f in [item, entry, 'Many']) + os.linesep)
    fReadOut.close()

if __name__ == '__main__':
    Main()
           
    
    
    
