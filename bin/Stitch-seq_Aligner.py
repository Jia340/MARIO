#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Name:        Stitch-seq_Aligner
# Purpose:
#
# Author:      Pengfei, Xiaoyi
#
# Created:     21/08/2012
# Modified:    03/09/2013
# Copyright:   (c) Pengfei 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# Updated 20151125:
#   decoupled both strands and would take -P indicating the strand

import sys, os, argparse, gc, json, copy, errno
import pysam
import itertools,string
from Bio import SeqIO
from xplib import TableIO
from xplib import DBI
from xplib.Annotation import Bed
from Annotation import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def check_negative(value):
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue

def ParseArg():
    p=argparse.ArgumentParser( description = 'Align miRNA-mRNA pairs for Stitch-seq. print the alignable miRNA-mRNA pairs with coordinates', epilog = 'Library dependency: Bio, pysam, itertools')
    p.add_argument('inputFile',type=str,metavar='part_reads',help='a part of paired fasta/fastq file')
    #p.add_argument('input2',type=str,metavar='part2_reads',help='paired part2 fasta/fastq file')

    groupMapper = p.add_mutually_exclusive_group()
    groupMapper.add_argument('-bt', '--bowtie_path',type=str,metavar='bowtie_path',help="path for the bowtie executive")
    groupMapper.add_argument('-bt2', '--bowtie2_path', type=str, metavar="bowtie2_path", help="path for bowtie2 executive")
    groupMapper.add_argument('-th', '--tophat_path', type=str, metavar="tophat_path", help="path for tophat executive")
    p.add_argument('-bl', '--blat_path', type=str, metavar='blat_path', help="path for the blat executive (for miRNA only)", default='blat')
    p.add_argument('-f','--fastq_to_fasta_path', dest='f2fpath', default='fastq_to_fasta', type=str, metavar='fastq_to_fasta_path', help="path for the fastq_to_fasta executive (for miRNA only)")
    p.add_argument('-s','--samtools_path',dest='spath', type=str,metavar='samtool_path',help="path for samtools executive",default='samtools')
    p.add_argument('--ref',type=str,metavar='rr',default="mm9",help="Reference seq files, references will be used in the order specified",nargs="*")

    # notice that rna.fa file should be fetched from Ensembl with the following attributes specified in the following order:
    # Ensembl Transcript ID
    # cDNA sequences (this is the actual sequence and is not in header)
    # Transcript start (bp)
    # Transcript end (bp)
    # Strand
    # 5' UTR start
    # 5' UTR end
    # 3' UTR start
    # 3' UTR end
    # Associated gene name
    # Transcript biotype

    p.add_argument('--reftype',type=str,metavar='rt',default="genome",help="Reference types: (miRNA, genome, transcript, other)",nargs="*")
    p.add_argument("-P", '--part_num', dest='partNum', type=check_negative, default=1, help="Specify which part (1 or 2) of the paired read is being mapped (this is used for strand matching purposes).")
    p.add_argument("-nostr", "--ignore_strand", dest="nostr", action = "store_true", help="Reads mapped onto the wrong strand will be considered as not mapped by default. Set this flag to ignore strand information.")
    p.add_argument('-p','--threads', type=check_negative, dest='nthreads', default=1, help="Number of threads used in bowtie mapping.")
    p.add_argument('-r', '--resume', action="store_true", dest="recover", help="Set to let Stitch-seq recover from previous files. Parameters other than number of threads need to be exactly the same for recovery. This may be useful if Stitch-seq crashes for CPU/memory/storage reasons.")
    p.add_argument('-l', '--mirnalen', type=check_negative, dest="mirnalen", default=35, help="Set the maximum length allow for a miRNA alignment (default = 35), a higher value may recover more miRNA alignments if there are references at such length but will be slower.")
    p.add_argument('-o', '--output_path', type=str, default='output', help='The path where all output files will be written to. Input files with the same name may need different output path values to prevent overwriting. default value: \'output\'')
    p.add_argument('-q', '--mapq_threshold', dest='mapqThreshold', type=check_negative, default=2, help='The MAPQ number threshold to detect gaps. Reads with MAPQ lower than this number but has a large difference between XS may be caused by RNA junctions and therefore should be considered as unmapped (to be remapped by later transcriptome mapping procedures, if any). default value: 2')
    p.add_argument('-G', '--gtf', type = str, help='GTF file for gene annotation. (Only used when Tophat is used as mapper.)')
    p.add_argument('-i', '--transcriptome-index', dest = 'trans_index', type = str, help='Index file for transcriptome. (Only used when Tophat is used as mapper. Will override -G if the file is up-to-date.)')

    if len(sys.argv)==1:
        #print (p.print_help(),file=sys.stderr)
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

rev_table=string.maketrans('ACGTacgtN', 'TGCAtgcaN')
def revcomp(seq, rev_table):
    return seq.translate(rev_table)[::-1]

def blat_align(b_path, read, ref, fastqToFasta, recovering, unmapfile, mirnalen, outputPath):
    # this is used for miRNA mapping only
    # b_path: blat path;
    # fastqToFasta: fastq_to_fasta path;
    # 1. use blat to map reads to miRNA references
    # 2. filter out the wrong ones:
    #   a) with wrong strand or
    #   b) not completely mapped to miRNA (have overlangs longer than 2 on both sides in total)
    # 3. annotate the correct entries with miRNA
    # 4. because blat does not generate unmapped reads, generate unmap file for the next ref

    # readfile = read.split("/")[-1].split(".")[0] + ".sam"
    tmp = read.rsplit("/", 1)
    filenametmp = tmp[-1].split(".")
    pathBase = "./"
    if outputPath.strip():
        pathBase = outputPath + '/'

    if not filenametmp[-1] in ["fa","fasta"]:   # fastq needs to be converted
        fastafile = pathBase + ".".join(filenametmp[:-1]) + ".fasta"
        os.system(fastqToFasta + " -i " + read + " -o " + fastafile);
    else:
        fastafile = read

    fastatoblat = pathBase + ".".join(filenametmp[:-1]) + "_blat.fasta"

    output = pathBase + ".".join(filenametmp[:-1]) + ".blatresult"

    if not recovering or (not os.path.isfile(output)) or os.stat(output).st_size <= 0:
        # a new file is definitely required
        # first write file to blat
        fseq = open(fastafile, 'r')
        fseqblat = open(fastatoblat, 'w')
        fsequnmap = open(unmapfile, 'w')
        seqname = ''
        oldseq = ''

        for line in fseq:
            if line.startswith('>') or line.startswith('@'):
                if oldseq and seqname:
                    if len(oldseq) <= mirnalen:
                        fseqblat.write(">" + seqname + os.linesep)
                        fseqblat.write(oldseq + os.linesep)
                    else:
                        fsequnmap.write(">" + seqname + os.linesep)
                        fsequnmap.write(oldseq + os.linesep)
                seqname = line.strip('>@').split(" ")[0].strip()
                oldseq = ''
            elif line.startswith('+'):
                # fastq
                if oldseq and seqname:
                    if len(oldseq) <= mirnalen:
                        fseqblat.write(">" + seqname + os.linesep)
                        fseqblat.write(oldseq + os.linesep)
                    else:
                        fsequnmap.write(">" + seqname + os.linesep)
                        fsequnmap.write(oldseq + os.linesep)
                seqname = ''
                oldseq = ''
            elif seqname:
                oldseq += line.strip()

        if oldseq and seqname:
            if len(oldseq) <= mirnalen:
                fseqblat.write(">" + seqname + os.linesep)
                fseqblat.write(oldseq + os.linesep)
            else:
                fsequnmap.write(">" + seqname + os.linesep)
                fsequnmap.write(oldseq + os.linesep)

        fseq.close()
        fseqblat.close()
        fsequnmap.close()
        print >> sys.stderr, 'Start mapping ...',
        #os.system("awk 'BEGIN {RS = \">\" ; ORS = \"\\n\"; FS = \"\\n\"; OFS=\"\\n\"} length($2) >= " + str(mirnalen) + " && NR > 1 {split($1, name, \" \"); print \">\"name[1], $2}' " + fastafile + " > " + unmapfile)
        # then write file to unmap
        os.system(b_path+ " -stepSize=2 -minScore=15 -tileSize=6 "+ref+" "+fastatoblat+" "+output+" >> /dev/null")
        print >> sys.stderr, 'done.'
        # otherwise just return the file name and use the old file
        # notice that this function will not handle whether the parameters are the same
    else:
        print >> sys.stderr, 'Old file exists, recovery in process.'

    return (output, fastatoblat)

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

def blat_annotation(outputfilename, typename, readfilename, unmapfilename, anno = True, requireUnique = False, posstrand = True, strandenforced = False, mismatchthr = 2, results_dict = dict()):
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
                    if anno:
                        curr_anno_arr = [tokens[13], tokens[15], tokens[16], tokens[8], seqdict[readname], typename, typename, tokens[13], '.']
                        if not strandenforced:
                            curr_anno_arr.append(strandcol)
                        curr_anno = '\t'.join(curr_anno_arr)
                    else:
                        curr_anno_arr = [tokens[13], tokens[15], tokens[16], tokens[8], seqdict[readname], typename]
                        if not strandenforced:
                            curr_anno_arr.append(strandcol)
                        curr_anno = '\t'.join(curr_anno_arr)
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
                        if anno:
                            curr_anno_arr = [tokens[13], tokens[15], tokens[16], tokens[8], seq, typename, typename, tokens[13], '.', strandcol]
                            newdict[readname].append('\t'.join(curr_anno_arr))
                        else:
                            curr_anno_arr = [tokens[13], tokens[15], tokens[16], tokens[8], seq, typename, strandcol]
                            newdict[readname].append('\t'.join(curr_anno_arr))

    fres.close()

    # merge annotation
    newanno = dict(results_dict.items() + newdict.items())

    # write (append to) unmapped reads
    funmap = open(unmapfilename, 'a')
    for key in seqdict:
        funmap.write(">" + key + os.linesep)
        funmap.write(seqdict[key] + os.linesep)

    funmap.close()

    return newanno


def bowtie_align(b_path,read,ref,s_path,bowtie2,numOfThreads,nOffrate,reftype,recovering,unmap,outputPath,mapqThreshold):
    # b_path: bowtie path;
    # s_path: samtools path;
    # bowtie2: logic, true/false
    # offrate is not used due to bowtie2 bug
    readPath = '/'.join(read.rsplit("/", 1)[:-1])
    readFileStem = read.rsplit("/", 1)[-1].rsplit(".", 1)[0]
    readPathSlash = './'

    if outputPath.strip():
        readPathSlash = outputPath + '/'

    sam = readPathSlash + readFileStem + ".sam"
    bam = readPathSlash + readFileStem + ".bam"

    hasFile = False

    if recovering and os.path.isfile(readPathSlash + "sort_" + readFileStem + ".bam"):
        # old file exists
        try:
            align = pysam.Samfile(readPathSlash + "sort_" + readFileStem + ".bam", "rb")
            hasFile = True
        except:
            hasFile = False

    if (not recovering) or (not hasFile):
        
        print >> sys.stderr, 'Start mapping.'

        if read.split(".")[-1] in ["fa","fasta"]:   # allow fasta and fastq for read
            foption=" -f"
        else:
            foption=""

        #if reftype == 'genome':
            # is genome, use temporary result file and no unmap output
        if bowtie2:
            samFileTarget = readPathSlash + "preResults.sam"
            bamFileTarget = readPathSlash + "preResults.bam"
        else:
            samFileTarget = sam
            bamFileTarget = bam

        if ref.split(".")[-1] in ["fa","fasta"]:
            base=ref.split("/")[-1].split(".")[0]
            os.system("rm " + readPathSlash + readFileStem + ".log")
            os.system(b_path+"-build "+ref+" "+base+" >> " + readPathSlash + readFileStem + ".log 2>&1")
        else:
            base = ref
        #if not bowtie2:
        #    os.system(b_path+ foption+" -a --best --strata -n 1 -l 15 -e 200 -p " + str(numOfThreads) + " --un " + unmap + " -S "+base+" "+read+" "+sam+" >> " + readPathSlash + readFileStem + ".log 2>&1")
        #else:
        #    os.system(b_path+ " -x "+base+foption+" -U "+read+ " -p " + str(numOfThreads) + " -i S,1,0.50 0R 3 -L 15 -D 20 -t " + ("-a " if reftype != "genome" else "") + " --un " + unmap + " -S "+sam+" >> " + readPathSlash + readFileStem + ".log 2>&1")
        #else:
        if not bowtie2:
            os.system(b_path+ foption+" -a --best --strata -n 1 -l 15 -e 200 -p " + str(numOfThreads) + " --un " + unmap + " -S "+base+" "+read+" "+samFileTarget+" >> " + readPathSlash + readFileStem + ".log 2>&1")
        else:
            os.system(b_path+ " -x " + base + foption + " -U " + read + " -p " + str(numOfThreads) + " -L 15 -D 20 -t " + (("-a") if reftype != "genome" else '') + " -S " + samFileTarget + " >> " + readPathSlash + readFileStem + ".log 2>&1")

        os.system(s_path + ' view -Sb -o ' + bamFileTarget + ' ' + samFileTarget)
        os.system("rm " + samFileTarget)
        
        # add genome filtering here, after that the output files will become sam and bam
        if bowtie2:
        #if reftype == 'genome':
            print >> sys.stderr, 'Filtering mapping results...'
            align = pysam.Samfile(bamFileTarget, 'rb')
            alignOutput = pysam.Samfile(bam, 'wb', template=align)

            unmapDict = dict()
            mappedDict = dict()
            
            for record in align:
                isMapped = not record.is_unmapped
                hasAS = False

                if isMapped:
                    try:
                        AS=record.opt("AS")
                        hasAS = True
                        # first filter based on previous best score
                        if '[AS=' in record.qname:
                            # has old best AS in last mapping
                            oldAS = int(record.qname.partition('[AS=')[-1].partition(']')[0])
                            # remove old AS in qname
                            record.qname = record.qname.partition('[AS=')[0] + record.qname.partition('[AS=')[-1].partition(']')[-1]
                            
                            if oldAS >= AS:
                                # old AS not worse than this new AS
                                AS = oldAS
                                isMapped = False
                                raise
                            
                        # then do genome specific stuff
                        if reftype == 'genome':
                            if record.mapq < mapqThreshold: # if mapq > mapqthreshold, may be unmapped (fall onto downstream mapping)
                                try:
                                    XS=record.opt("XS")
                                    if AS-XS >= 3:
                                        # wrongly mapped, consider as unmapped (to be rescued in lower mapping)
                                        isMapped = False
                                except:
                                    isMapped = False
                    except:
                        pass

                if isMapped:
                    alignOutput.write(record)
                    if record.qname in unmapDict:
                        del unmapDict[record.qname]
                    mappedDict[record.qname] = True
                else:
                    if record.qname not in mappedDict:
                        if record.qname not in unmapDict:
                            unmapDict[record.qname] = dict()
                        
                        if 'AS' not in unmapDict[record.qname] or (hasAS and unmapDict[record.qname]['AS'] < AS):
                            seq = record.seq
                            if record.is_reverse:
                                seq = revcomp(record.seq, rev_table)
                            if hasAS:
                                unmapDict[record.qname]['AS'] = AS
                            unmapDict[record.qname]['rec'] = SeqRecord(Seq(seq, IUPAC.unambiguous_dna), id = record.qname + (('[AS=' + str(AS) + ']') if hasAS else '') , description='')

            align.close()
            alignOutput.close()

            funmap = open(unmap, 'w')
            for unmapStuff in unmapDict.itervalues():
                SeqIO.write(unmapStuff['rec'], funmap, "fasta")
            funmap.close()

        os.system(s_path+ " sort "+ bam + " " + readPathSlash + "sort_" + readFileStem)
        os.system("rm " + bam)
        align=pysam.Samfile(readPathSlash + "sort_" + readFileStem + ".bam", "rb")
        print >> sys.stderr, 'Mapping completed.'

    else:
        print >> sys.stderr, 'Old file exists, recovery in process.'
    
    return align

def tophat_align(t_path,read,ref,s_path,gtf,trans_index,numOfThreads,reftype,recovering,unmap,outputPath,mapqThreshold):
    # t_path: tophat path;
    # s_path: samtools path;
    readPath = '/'.join(read.rsplit("/", 1)[:-1])
    readFileStem = read.rsplit("/", 1)[-1].rsplit(".", 1)[0]
    readPathSlash = './'

    if outputPath.strip():
        readPathSlash = outputPath + '/'

    sam = readPathSlash + readFileStem + ".sam"
    bam = readPathSlash + readFileStem + ".bam"

    hasFile = False

    if recovering and os.path.isfile(readPathSlash + "sort_" + readFileStem + ".bam"):
        # old file exists
        try:
            align = pysam.Samfile(readPathSlash + "sort_" + readFileStem + ".bam", "rb")
            hasFile = True
        except:
            hasFile = False

    if (not recovering) or (not hasFile):
        
        print >> sys.stderr, 'Start mapping.'

        #if reftype == 'genome':
            # is genome, use temporary result file and no unmap output
        if ref.split(".")[-1] in ["fa","fasta"]:
            print >> sys.stderr, 'Reference should be bowtie2-indexed, please run bowtie2-build first.'
            return
        else:
            base = ref
        #if not bowtie2:
        #    os.system(b_path+ foption+" -a --best --strata -n 1 -l 15 -e 200 -p " + str(numOfThreads) + " --un " + unmap + " -S "+base+" "+read+" "+sam+" >> " + readPathSlash + readFileStem + ".log 2>&1")
        #else:
        #    os.system(b_path+ " -x "+base+foption+" -U "+read+ " -p " + str(numOfThreads) + " -i S,1,0.50 0R 3 -L 15 -D 20 -t " + ("-a " if reftype != "genome" else "") + " --un " + unmap + " -S "+sam+" >> " + readPathSlash + readFileStem + ".log 2>&1")
        #else:
        os.system(t_path + " --b2-L 15 --b2-D 20" + " -p " + str(numOfThreads) + " -o " + readPathSlash + "tophat_out" + ((" -G " + gtf) if gtf else '') + ((' --transcriptome-index ' + trans_index) if trans_index else '') + ' ' + base + ' ' + read + " >> " + readPathSlash + readFileStem + ".log 2>&1")

        os.system("mv "+ readPathSlash + "tophat_out/accepted_hits.bam " + readPathSlash + "sort_" + readFileStem + '.bam')
        os.system("mv "+ readPathSlash + "tophat_out/unmapped.bam " + unmap + '.bam')
        align=pysam.Samfile(readPathSlash + "sort_" + readFileStem + ".bam", "rb")
        print >> sys.stderr, 'Mapping completed.'

    else:
        print >> sys.stderr, 'Old file exists, recovery in process.'
    
    return align

def createPath(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def Main():
    args=ParseArg()

    inRecovery = False

    nThreads = args.nthreads
    needRecover = args.recover
    #nOffrate = args.offrate
    nOffrate = 5

    del args.nthreads
    del args.recover
    #del args.offrate
    paramDict = vars(args)

    newpar = json.dumps(paramDict, sort_keys = True)

    paramDictClean = copy.deepcopy(paramDict)
    del paramDictClean['ref']
    del paramDictClean['reftype']

    newparClean = json.dumps(paramDictClean, sort_keys = True)

    oldRef = None
    oldRefType = None

    if needRecover and os.path.isfile('.StitchSeqParameters'):
        fpar = open('.StitchSeqParameters', 'r')
        oldpar = fpar.read()

        oldDict = json.loads(oldpar)
        oldDictClean = copy.deepcopy(oldDict)
        del oldDictClean['ref']
        del oldDictClean['reftype']

        oldparClean = json.dumps(oldDictClean, sort_keys = True)

        if oldparClean == newparClean:
            print >> sys.stderr, "Same old parameters found, resuming ..."
            oldRef = oldDict['ref']
            oldRefType = oldDict['reftype']
                        
        fpar.close()
        inRecovery = True

    fpar = open('.StitchSeqParameters', 'w')
    fpar.write(newpar)
    fpar.close()

    if args.bowtie2_path:
        args.bowtie_path = args.bowtie2_path
        args.bowtie2 = True
    
    refList = args.ref
    refTypeList = args.reftype
    #annodictlist = [dict(), dict()]
    strands = (args.partNum == 1)
    strandenforced = not args.nostr

    inRecoveryFrag = inRecovery

    strand = (args.partNum == 1)

    annodictentry = dict()
    if not args.output_path.rstrip('/') == '.':
        output = args.output_path.rstrip('/')
        createPath(output)
    
    gc.collect()
    # mapping for every single ends
    # also put the annotations to the file
    # then merge all the annotations from every single step
    # then try to match all the annotations

    readfile = args.inputFile
    print >> sys.stderr, 'Mapping ' + readfile + ' ...' + (' (Recovery possible)' if inRecovery else '')

    for iref in xrange(len(refList)):
        reffile = refList[iref]
        reftype = refTypeList[iref]

        if inRecoveryFrag:
            inRecoveryFrag = False
            try:
                if reffile == oldRef[iref] and reftype == oldRefType[iref]:
                    inRecoveryFrag = True
            except:
                pass

        # unmapped read file
        print >> sys.stderr, 'Mapping ' + readfile + ' with ref: ' + reftype + (' (Recovering from old reads)' if inRecoveryFrag else '')
        readPath = '/'.join(readfile.rsplit("/", 1)[:-1])
        readPathSlash = './'
        readFileStem = readfile.rsplit("/", 1)[-1].rsplit(".", 1)[0]
        if output.strip():
            readPathSlash = output + '/'
        unmap_read = readPathSlash + readFileStem + "_unmap." + readfile.rsplit("/", 1)[-1].rsplit(".", 1)[-1]

        if reftype.lower() == "mirna":
            outputfile, readused = blat_align(args.blat_path, readfile, reffile, args.f2fpath, inRecoveryFrag, unmap_read, args.mirnalen, output)
            annodictentry = blat_annotation(outputfile, reftype, readused, unmap_read, annotation, False, strand, strandenforced, 2, annodictentry)
        elif args.tophat_path:
            outputbam = tophat_align(args.tophat_path, readfile, reffile, args.spath, args.gtf, args.trans_index, nThreads, reftype.lower(), inRecoveryFrag, unmap_read, output, args.mapqThreshold)
        else:
            outputbam = bowtie_align(args.bowtie_path, readfile, reffile, args.spath, args.bowtie2, nThreads, nOffrate, reftype.lower(), inRecoveryFrag, unmap_read, output, args.mapqThreshold)

        readfile = unmap_read
        gc.collect()

                    
if __name__ == '__main__':
    Main()




