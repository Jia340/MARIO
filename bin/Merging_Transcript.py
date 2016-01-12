import argparse,sys
import collections

Anno={}
mdead=[]
   
class ens_transcript():
  def __init__(self,x,source,**kwargs):
    if source=="ens":
      self.txID=x[0].strip()
      self.geneID=x[1].strip()
      self.chr='chr'+x[2].strip()
      if x[3]=="1":self.strand="+"
      else:self.strand="-"
      self.name=x[4]
      self.type=x[5]
      self.exon=[]
    else:
      pass
    for key in kwargs.keys():
      setattr(self,key,kwargs[key])

  def GenomicPosition(self,p):      
    strand=self.strand
    g=0
    for exon in self.exon:
      exon_l=exon[1]-exon[0]+1
      if p>exon_l:
        p=p-exon_l    
      else:
        if strand=="+":
          g=p+exon[0]-1
        elif strand=="-":
          g=exon[1]-p+1
        return g
    if strand=="+":
      return self.exon[-1][1]+p
    elif strand=="-":
      return self.exon[-1][0]-p
             
  def GenomicPosition_junction(self,s,e):
    starts=[s]
    ends=[]
    if self.strand=="-":
      exons=sorted(self.exon)
    else:
      exons=self.exon
    for exon in exons:
      if exon[1]>=s and exon[1]<=e:
        ends.append(exon[1])
      if exon[0]>=s and exon[0]<=e:
        starts.append(exon[0])
    ends.append(e)
    if len(starts)!=len(ends) or len(exons)==1:
      return (str(s),str(e))
    return (','.join([str(x) for x in starts]),','.join([str(y) for y in ends]))

def ParseArg():
  p=argparse.ArgumentParser(description="Merge fragments mapped to different transcripts of the same gene")
  p.add_argument("-i","--input",type=str,required=True,help="paired fragments. MUST BE SORTED BY LINKER NAME")
  p.add_argument("-md","--miRNA_dead",type=str,default="/data2/jialu/RNA_RNA/new_mapping/miRNA.dead")
  p.add_argument("-T","--annotation_trans",type=str,required=True,help="Ensembl annotation file for transcript")
  p.add_argument("-o","--Output",type=str,default="output.txt",help="Output file name")
  if len(sys.argv)==1:
      print >>sys.stderr,p.print_help()
      sys.exit(0)
  return p.parse_args()


def dead_miRNA(miRNA):
  '''
  half[7]
  '''
  if "Mir" in miRNA:
    if "mir-"+miRNA.split("Mir")[1] in mdead:
      return 1
  elif "mi" in miRNA:
    try:
      if "mir-"+miRNA.split("miR")[1].strip("-") in mdead:
        print >>sys.stderr,"mir-"+miRNA.split("miR")[1].strip("-")
        return 1
    except:
      if miRNA in mdead:
        return 1
  return 0

def Overlap(tuple1,tuple2):

  if "," in tuple1[0] or "," in tuple2[0]:
    return 0
  if abs(int(tuple1[0])-int(tuple2[0]))<=10:
    return 1
  return 0

    
def Merge(halfs):
  '''
  Merge many. halfs is a list with 10 elements.
  '''
  seq_l=len(halfs[0][4])
  
  if "miRNA" in halfs[0][5]:
    #find unique fragments and return them
    halfs_set=set()
    H=halfs[0]
    name=halfs[0][0]
    for h in halfs:
      if h[6]=="miRNA":  #check miRNA
        if dead_miRNA(h[7])==1:
  #        h[0]=h[0]+"_dead"
          h[8]="dead"
        if h[0]<name:
          H=h
          name=h[0]
    halfs_set.add('\t'.join(H))
      
    return halfs_set
#usually the multi-align reads are mapped to the same miRNA family. So we simply randomly choose one to report.
 
  AlignType=halfs[0][9]

  #unique
  halfs=[y.split("\t") for y in set(['\t'.join(x[0:9]) for x in halfs])]
  
  #Merging
  pc_type_dic={}
  pc_subtype_dic={}
  other_subtype_dic={}
  merged=set()
 
  for h in halfs:
    #check miRNA
#    if h[6]=="miRNA" and dead_miRNA(h[7])==1:
#      h[0]=h[0]+"_dead"
#      merged.add('\t'.join(h)+'\t'+AlignType)
#      continue
    #check tRNA

    if h[0] not in Anno:
      print "Not found: %s"%(h[0])
      merged.add('\t'.join(h)+'\t'+AlignType)
      continue

    if "LRG" in h[0]:
      continue

    g_id=Anno[h[0]].geneID
    g_chr=Anno[h[0]].chr
    new_pos=(Anno[h[0]].GenomicPosition(int(h[1])),Anno[h[0]].GenomicPosition(int(h[2])))
    g_s=min(new_pos)
    g_e=max(new_pos)
#    g_s=Genes[g_id].GenicPosition(Anno[h[0]].GenomicPosition(int(h[1])))
#    g_e=Genes[g_id].GenicPosition(Anno[h[0]].GenomicPosition(int(h[2])))  
    if (g_e-g_s+1)!=seq_l:
      new_pos_junction=Anno[h[0]].GenomicPosition_junction(g_s,g_e)
      h[1]=new_pos_junction[0]
      h[2]=new_pos_junction[1]
    else:    
      h[1]=str(g_s)
      h[2]=str(g_e)
    
    #annotation:key is geneID+start+end
    if Anno[h[0]].type=="protein_coding":   #if the gene is protein_coding gene
      h[0]=g_id+'_'+g_chr
      
      k='\t'.join(h[0:6])+"/"+h[7]  
      if k not in pc_type_dic:
        pc_type_dic[k]=set([h[6]])
      else:
        pc_type_dic[k].add(h[6])

      if k not in pc_subtype_dic:
        pc_subtype_dic[k]=set([h[8]])
      else:      
        pc_subtype_dic[k].add(h[8])  
    else:
      if h[6]=="miRNA" and dead_miRNA(h[7])==1:
        h[8]="dead"
      h[6]=Anno[h[0]].type
      h[0]=g_id+'_'+g_chr
      k='\t'.join(h[0:8])
      if k not in other_subtype_dic:
        other_subtype_dic[k]=set([h[8]])
      else:
        other_subtype_dic[k].add(h[8])
  
  for key in pc_type_dic.keys():
    if "protein_coding" in pc_type_dic[key]:
      try:
        pc_subtype_dic[key].remove(".")
      except:
        pass
      merged.add(key.split("/")[0]+"\tprotein_coding\t"+key.split("/")[1]+"\t"+"|".join(pc_subtype_dic[key])+"\t"+AlignType)
    else:
      merged.add(key.split("/")[0]+"\tprotein_coding-noncoding\t"+key.split("/")[1]+"\t"+"|".join(pc_subtype_dic[key])+"\t"+AlignType)
  for key in other_subtype_dic.keys():
    merged.add(key+"\t"+"|".join(other_subtype_dic[key])+"\t"+AlignType)  
  return list(merged)
    

def MergeTX(frags,output):
  '''
  frags is a list of fragments with same linker ID.
  '''
  #for genome, we only check if miRNA is a dead miRNA.
  if frags[0][5]=="genome":
    frag_set=set()
    for f in frags:
      if f[6]=="miRNA" and dead_miRNA(f[7])==1:
        f[8]="dead"
      frag_set.add('\t'.join(f))
    for m in frag_set:
      print >>output,m
    return

  linker=frags[0][10]
  Proper=[]
  NonProper=[]
  for f in frags:
    if f[9]=="ProperStrand":
      Proper.append(f[0:10])
    elif f[9]=="NonProperStrand":
      NonProper.append(f[0:10])

  if len(Proper)!=0:
    merged_f=Merge(Proper)
  elif len(Proper)==0:
    merged_f=Merge(NonProper)

  if len(merged_f)==1:
    type="One"
  elif len(merged_f)>1:
    type="Many"
  for m in merged_f:
    print >>output,m+'\t'+linker+'\t'+type
  return


def Main():
  args=ParseArg()

  ens_file=open(args.annotation_trans,"r")
  global Anno
  global Genes
  print >>sys.stderr, "Reading annotation file"
  ens_file.readline()  
  while True:
    line=ens_file.readline()
    if line=="":break
    line=line.strip().split("\t")
    if line[0] not in Anno:
      Anno[line[0]]=ens_transcript(line,"ens")
    Anno[line[0]].exon.append((int(line[6]),int(line[7])))
  for key in Anno.keys():
    if Anno[key].strand=="+":
      Anno[key].exon=sorted(Anno[key].exon)
    elif Anno[key].strand=="-":
      Anno[key].exon=sorted(Anno[key].exon,reverse=True)      
      
   
  mdead_file=open(args.miRNA_dead,"r")
  global mdead
  print >>sys.stderr, "Reading miRNA_dead"
  while True:
    line=mdead_file.readline()
    if line=="":break
    else:
      mdead.append((line.strip().split("   ")[1]).split("-",1)[1])

  frag_file=open(args.input,"r")
  output=open(args.Output,"w")
  
  c=0
  prev=""
  while True:
    line=frag_file.readline()
    if line=="":
      MergeTX(frags,output)
      break 
    line=line.strip().split("\t")    
    linker=line[10]
    if linker!=prev:
      if c!=0:
        MergeTX(frags,output)
      frags=[line]
      prev=linker
    else:
      frags.append(line)
    c+=1
    if (c%1000)==0:
      print >>sys.stderr, "Merging %d fragments\r"%(c),



if __name__=='__main__':
  Main()         
    

