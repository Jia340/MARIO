#!/usr/bin/env python
# Pengfei Yu yupf05@gmail.com

import sys,os,argparse
from multiprocessing import *
from functools import partial
from Bio.Blast import NCBIXML
import itertools
from Bio import SeqIO
from time import time


'''
dir(hsp):
['__class__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand']
'''

class Shared_queue:
    def __init__(self):
        m=Manager()
        self.qout=m.Queue()
        self.qpair1=m.Queue()
        self.qpair2=m.Queue()
        self.qf=m.Queue()  ##front only
        self.qb=m.Queue()  ##back only
        self.qn=m.Queue()  ##no linker

    def queue_wait(self):
        self.qout.put("wait")
        self.qpair1.put("wait")
        self.qpair2.put("wait")
        self.qf.put("wait")
        self.qb.put("wait")
        self.qn.put("wait")


    def queue_end(self):
        self.qout.put("kill")
        self.qpair1.put("kill")
        self.qpair2.put("kill")
        self.qf.put("kill")
        self.qb.put("kill")
        self.qn.put("kill")

    def queue_join(self):
        self.qout.join()
        self.qpair1.join()
        self.qpair2.join()
        self.qf.join()
        self.qb.join()
        self.qn.join()


class Shared_count:
    def __init__(self):
        self.ntotal=0
        self.nnolinker=Value("i",0)
        self.npair=Value("i",0)
        self.nfront=Value("i",0)
        self.nback=Value("i",0)


def ParseArg():
    
    p=argparse.ArgumentParser( description = 'DESCRIPTION: Run BLAST, find linker sequences and split two parts connected by linkers', epilog='Library dependency: Bio, itertools')
    p.add_argument("input",type=str,help="the input fasta file containing fragment sequences of type1 and type2")
    p.add_argument("type3_1",type=str,help="read_1 for evenlong (type3) fastq file")
    p.add_argument("type3_2",type=str,help="read_2 for evenlong (type3) fastq file")
    p.add_argument("-e","--evalue",dest="evalue",type=float,default=0.00001,help="cutoff evalues, only choose alignment with evalue less than this cutoffs (default: 1e-5).")
    p.add_argument("-d","--evalue_end",dest="evalue_end",type=float,default=0.001,help="cutoff evalues at end, only choose alignment at the end with evalue less than this cutoffs (default: 1e-3).")
    p.add_argument("--linker_db",dest="linker_db",type=str,help="BLAST database of linker sequences",default="/home/yu68/Stitch-seq/blast_db/linker.fa")
    p.add_argument("--blast_path",dest="blast_path",type=str,help="path for the local blast program",default="/home/yu68/Softwares/ncbi-blast-2.2.27+/bin/blastn")
    p.add_argument("-o","--output",dest="output",type=str,help="output file containing sequences of two sepatated parts")
    p.add_argument("-t","--trim",type=int,default=10,help="trim off the first this number of nt as index, default:10")
    p.add_argument("-b","--batch",type=int,default=200000,help="batch this number of fragments for BLAST at a time. default: 200000")
    p.add_argument("-r","--release",action='store_true',help="set to allow released criterion for Paired fragment in Type 3, include those ones with no linker in two reads")
    p.add_argument("-l","--length",type=int,default=15,help="shortest length to be considered for each part of the pair, default: 15")
    p.add_argument("-k","--linker_trim_length",type=int,default=10,help="linker length to be truncated for Paired fragment in Type 3, default: 10")
    p.add_argument("-p","--parallel",type=int,default=5,help="Number of threads")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def blast_align(fasta,blast_path,linker_db):
    fasta_name=fasta.split(".")[0]
    os.system(blast_path+" -task blastn -outfmt 5 -num_threads 6 -evalue 0.1 -db "+linker_db+" -query ./temp/"+fasta+" > ./temp/"+fasta_name+"_blast_linker.xml")
    linker_records=NCBIXML.parse(open("./temp/"+fasta_name+"_blast_linker.xml"))
#    os.system("rm ./temp/"+fasta)
#    os.system("rm ./temp/"+fasta_name+"_blast_linker.xml")
    return (linker_records)

def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch
            
    
    
def batch_iterator2(iterator1,iterator2,batch_size):
    """ batch iterator for paired-end files"""        
    entry1 = True #Make sure we loop once
    while entry1 :
        batch1 = []
        batch2 = []
        while len(batch1) < batch_size :
            try :
                entry1 = iterator1.next()
                entry2 = iterator2.next()
            except StopIteration :
                entry1 = None
            if entry1 is None :
                #End of file
                break
            batch1.append(entry1)
            batch2.append(entry2)
        if batch1 :
            yield batch1, batch2

def RecoverType12(filename, q_obj, Args, types):
    ##Put the list which needs to be mapped at first
    '''        self.qout=m.Queue()
        self.qpair1=m.Queue()
        self.qpair2=m.Queue()
        self.qf=m.Queue()  ##front only
        self.qb=m.Queue()  ##back only
        self.qn=m.Queue()  ##no linker
    '''
    lout=[]
    lpair1=[]
    lpair2=[]
    lfront=[]
    lback=[]
    lnolinker=[]


    record_count=0

    trim_n=Args.trim
    min_l=Args.length

    blast_path=Args.blast_path
    linker_db=Args.linker_db

    # E-values
    evalue=float(Args.evalue)

    # End E-values
    evalueEnd = float(Args.evalue_end)
    
    # Linker trim length
    linkerTrimLen = int(Args.linker_trim_length)

    linker_records = blast_align(filename,blast_path,linker_db)
    batch_temp=SeqIO.parse(open('./temp/'+filename,"rU"),types)


    for linker_record, fragment in itertools.izip(linker_records,batch_temp):
        record_count+=1
        line=''
        start=len(fragment) # record the start of all linker alignment
        end=0 # record the end of all linker alignment
        pos={} # position of all alignment for one fragment
        counts={} # counts of linker alignment for one fragment
        Types="None"
        for alignment in linker_record.alignments:
            expect=evalue
            for hsp in alignment.hsps:
                if hsp.expect < expect:
                    start=min(hsp.query_start - 1, hsp.query_end - 1, start)
                    end=max(hsp.query_end,hsp.query_start,end)
                    if alignment.hit_def in counts.keys():
                        counts[alignment.hit_def]+=1
                    else:
                        counts[alignment.hit_def]=1
                    pos[",".join (str(f) for f in [hsp.query_start,hsp.query_end])]=alignment.hit_def
        stat=";".join (f+':%d'%(counts[f]) for f in sorted(counts.iterkeys()))
        pos_order=";".join (str(f) for f in sorted(pos.iterkeys()))
        order=";".join (str(pos[f]) for f in sorted(pos.iterkeys()))

        if start>end:
            if len(fragment)-trim_n>min_l:    
                #XXXXXXXXXX11111111222222222 or XXXXXXXXXX11111111111111111 or XXXXXXXXXX2222222222222222
                lnolinker.append(fragment[trim_n:])
                Types="NoLinker"
        elif (start>trim_n+min_l) and (end<len(fragment)-min_l):
            #XXXXXXXXXX111111LLLLLLL2222222
            lpair1.append(fragment[trim_n:start])
            lpair2.append(fragment[end:].reverse_complement(fragment.id))
            Types="Paired"
        elif (start>trim_n+min_l):
            #XXXXXXXXXX1111111LLLLLLL
            lfront.append(fragment[trim_n:start])
            Types="FrontOnly"
        elif (end<len(fragment)-min_l):
            #XXXXXXXXXXLLLLLL22222222
            lback.append(fragment[end:].reverse_complement(fragment.id))
            Types="BackOnly"

        lout.append("\t".join([fragment.id,stat,pos_order,Types,order]))
    q_obj.qout.put(lout)
    q_obj.qpair1.put(lpair1)
    q_obj.qpair2.put(lpair2)
    q_obj.qf.put(lfront)
    q_obj.qb.put(lback)
    q_obj.qn.put(lnolinker)
    return record_count


def RecoverType3(filename, q_obj, Args, types):
    '''
    self.qpair1=m.Queue()
        self.qpair2=m.Queue()
        self.qf=m.Queue()  ##front only
        self.qb=m.Queue()  ##back only
        self.qn=m.Queue()  ##no linker
    '''
#    print filename1
#    print filename2
    lout=[]
    lpair1=[]
    lpair2=[]
    lfront=[]
    lback=[]
    record_count=0
    nolinker_count=0

    trim_n=Args.trim
    min_l=Args.length

    blast_path=Args.blast_path
    linker_db=Args.linker_db

    # E-values
    evalue=float(Args.evalue)

    # End E-values
    evalueEnd = float(Args.evalue_end)
    
    # Linker trim length
    linkerTrimLen = int(Args.linker_trim_length)

    release=Args.release

    filename1=filename[0]
    filename2=filename[1]
    linker_records1 = blast_align(filename1,blast_path,linker_db)
    linker_records2 = blast_align(filename2,blast_path,linker_db)

    #print "BLAST aligned for %s.and %s" % (filename1, filename2)

    batch_temp1=SeqIO.parse(open('./temp/'+filename1,"rU"),"fasta")
    batch_temp2=SeqIO.parse(open('./temp/'+filename2,"rU"),"fasta")
 #   print batch_temp1
    for linker_record1, linker_record2 in itertools.izip(linker_records1,linker_records2):
        read1=batch_temp1.next()
        read2=batch_temp2.next()
 #       print read1,read2

        record_count+=1
        if read1.id!=read2.id:
            print "ERROR!! ID not match for paired type3 reads"
            sys.exit(0)
        start1=len(read1)
        for alignment in linker_record1.alignments:
            expect=evalue
            for hsp in alignment.hsps:
                if max(hsp.query_start, hsp.query_end) >= len(read1):
                    expect = evalueEnd
                else:
                    expect = evalue
                if hsp.expect < expect:
                    start1=min(hsp.query_start - 1, hsp.query_end - 1, start1)
        start2=len(read2)
        for alignment in linker_record2.alignments:
            expect=evalue
            for hsp in alignment.hsps:
                if max(hsp.query_start, hsp.query_end) >= len(read2):
                    expect = evalueEnd
                else:
                    expect = evalue
                if hsp.expect < expect:
                    start2=min(hsp.query_start - 1, hsp.query_end - 1, start2)
      
        if (start1==len(read1) and start2==len(read2)):
            nolinker_count+=1
            if not release: continue
            start1 = len(read1) - linkerTrimLen
            start2 = len(read2) - linkerTrimLen
        if (start1>trim_n+min_l) and (start2>min_l):
            lpair1.append(read1[trim_n:start1])
            lpair2.append(read2[:start2])
            Types="Paired"
        elif start1>trim_n+min_l:
            lfront.append(read1[trim_n:start1])
            Types="FrontOnly"
        elif start2>min_l:
            lback.append(read2[:start2])
            Types="BackOnly"

    q_obj.qpair1.put(lpair1)
    q_obj.qpair2.put(lpair2)
    q_obj.qf.put(lfront)
    q_obj.qb.put(lback)
    return (record_count,nolinker_count)

def Writer(q_obj,c_obj,Input,fout,Type):
    '''
    The process which write the output files.
            self.qout=m.Queue()
        self.qpair1=m.Queue()
        self.qpair2=m.Queue()
        self.qf=m.Queue()  ##front only
        self.qb=m.Queue()  ##back only
        self.qn=m.Queue()  ##no linker
                self.ntotal=Value("i",0)
        self.nnolinker=Value("i",0)
        self.npair=Value("i",0)
        self.nfront=Value("i",0)
        self.nback=Value("i",0)

    '''
    Input_name=Input.split("/")[-1]
    fout=open(fout,"w")
    fpair1=open("Paired1_"+Input_name,"w")
    fpair2=open("Paired2_"+Input_name,"w")
    ffront=open("frontOnly_"+Input_name,"w")
    fback=open("backOnly_"+Input_name,"w")
    fnolinker=open("NoLinker_"+Input_name,"w")
    Switch=[1,1,1,1,1,1]
    Q=[q_obj.qout, q_obj.qf, q_obj.qb, q_obj.qn, q_obj.qpair1, q_obj.qpair2]
    F=[fout,ffront,fback,fnolinker,fpair1,fpair2]
    COUT=[0, c_obj.nfront, c_obj.nback, c_obj.nnolinker, c_obj.npair, 0]
    C=[0,0,0,0,0,0]
    while sum(Switch)>0:
        for i in xrange(6):
            if Switch[i]!=0:
                q=Q[i]
                f=F[i]
                if not q.empty():
                    item=q.get()
                    if item=="wait":
                        if COUT[i]!=0:
                            COUT[i].value=C[i]
                        C[i]=0
                        q.task_done()
                        continue
                    if item=="kill":
                        if COUT[i]!=0:
                            COUT[i].value=C[i]
                        q.task_done()
                        Switch[i]=0
                        continue
                    if q==q_obj.qout:
                        print >> fout, "\n".join(item)
                    else:
                        SeqIO.write(item,f,Type)
                        C[i]+=len(item)
                    q.task_done()
    fout.close()
    ffront.close()
    fback.close()
    fnolinker.close()
    fpair1.close()
    fpair2.close()
    return


def main():

    os.system("mkdir temp") # create temp folder
    
    args=ParseArg()
    
    name=args.output.split(".")[0].split("_")[0]

    blast_path=args.blast_path
    linker_db=args.linker_db

    # create blast database if not yet
    if not os.path.isfile(linker_db+".nhr"):
        print >> sys.stderr, "  ## blast database not detected, making blast db for seq"
        blastdir,_ = os.path.split(blast_path)
        os.system(blastdir+"/makeblastdb -in "+linker_db+" -dbtype 'nucl' -title "+os.path.splitext(linker_db)[0])

  
    # determine input sequence file type
    types="fastq"
    if args.input.split(".")[-1] in ["fa","fasta"]:
        types="fasta"

    seq_file=SeqIO.parse(args.input,types)
    
    # determine input type3 sequence files type
    types2="fastq"
    if args.type3_1.split(".")[-1] in ["fa","fasta"]:
        types2="fasta"
    
    seq_type3_1=SeqIO.parse(args.type3_1,types2)
    seq_type3_2=SeqIO.parse(args.type3_2,types2)

    num_thread=args.parallel

    Qs=Shared_queue()
    Cs=Shared_count()
   
    ###################################
    ##        multiprocessing        ## 
    ###################################
    '''
    This loop is for recovered fragment (type1&2)
    '''
    print >>sys.stderr, "---- split partners for recovered fragments (Type1&2) ----"
    t0=time()
    filename_list=[]
    for i, batch in enumerate(batch_iterator(seq_file, args.batch)):
        filename=name+"group_%i.fasta" % (i+1)
        handle=open("./temp/"+filename, "w")
        count=SeqIO.write(batch,handle,"fasta")
        handle.close()
        filename_list.append(filename)

    #start a process for writing
    pw=Process(target=Writer,args=(Qs,Cs,args.input,args.output,types,))
    pw.start()

    #start the other processes for blast
    p=Pool(num_thread)
    Cs.ntotal=sum(p.map(partial(RecoverType12, q_obj=Qs, Args=args, types=types),filename_list))

    Qs.queue_wait()

    print >>sys.stderr, "---- Done blasting ----"

    p.close()
    p.join()
    Qs.queue_join()
    filename_list=[]
    print >>sys.stderr, "---- Done writing files ----"

    t1=time()
    print >>sys.stderr,"Type1-2, with_linker: (%i/%i). P:%i B:%i F:%i. Time: %.2f min."%(Cs.ntotal-Cs.nnolinker.value,Cs.ntotal, Cs.npair.value,Cs.nback.value,Cs.nfront.value,(t1-t0)/60) 

    '''
    This loop is for evenlong(type3) reads to find "Paired"/"BackOnly"/"FrontOnly"    
    '''
    i=0
    print >>sys.stderr, "\n---- split partners for type3 paired reads ----"
    t0=time()

    for batch1, batch2 in batch_iterator2(seq_type3_1, seq_type3_2, args.batch):
        filename1="type3_group_%i_1.fasta" % (i+1)
        filename2="type3_group_%i_2.fasta" % (i+1)
        handle1=open("./temp/"+filename1, "w")
        count1=SeqIO.write(batch1,handle1,"fasta")
        handle1.close()
        handle2=open("./temp/"+filename2, "w")
        count2=SeqIO.write(batch2,handle2,"fasta")
        handle2.close()
        i=i+1
        filename_list.append([filename1,filename2])

    #start the other processes for blast
    p=Pool(num_thread)
    type3_count=p.map(partial(RecoverType3, q_obj=Qs, Args=args, types=types),filename_list)

    Cs.ntotal=sum([x[0] for x in type3_count])
    Cs.nnolinker=sum(x[1] for x in type3_count)

    Qs.queue_end()

    print >>sys.stderr, "---- Done blasting ----"

    p.close()
    p.join()
    Qs.queue_join()
    pw.join()
    print >>sys.stderr, "---- Done writing files ----"

    t1=time()
    print >>sys.stderr,"Type3, with_linker: (%i/%i). P:%i B:%i F:%i. Time: %.2f min."%(Cs.ntotal-Cs.nnolinker,Cs.ntotal, Cs.npair.value,Cs.nback.value,Cs.nfront.value,(t1-t0)/60)  

    os.system("rm -r ./temp")


if __name__=="__main__":
    main()
