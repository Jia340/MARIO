"""
For the purpose of annotating RNA types for genomic regions.
"""


#from xplib import DBI
#from cogent.db.ensembl import HostAccount, Genome

def overlap(bed1,bed2):
    """
    This function compares overlap of two Bed object from same chromosome
    
    :param bed1: A Bed object from `xplib.Annotation.Bed <http://bam2xwiki.appspot.com/bed>`_ (BAM2X)
    :param bed2: A Bed object from `xplib.Annotation.Bed <http://bam2xwiki.appspot.com/bed>`_ (BAM2X)
    :returns: boolean -- True or False

    Example:
    
    >>> from xplib.Annotation import Bed
    >>> from Annotation import overlap
    >>> bed1=Bed(["chr1",10000,12000])
    >>> bed2=Bed(["chr1",9000,13000])
    >>> print overlap(bed1,bed2)
    True

    """
    try:
        return (bed1.stop>bed2.start) and (bed1.start<bed2.stop)
    except: # in case for "NonType" of bed2
        return False

def IsProperStrand(bed1, bed2):
    """
    This function determines whether the two Bed object is at the same strand
    :param bed1: A Bed object from `xplib.Annotation.Bed <http://bam2xwiki.appspot.com/bed>`_ (BAM2X)
    :param bed2: A Bed object from `xplib.Annotation.Bed <http://bam2xwiki.appspot.com/bed>`_ (BAM2X)
    :returns: boolean -- True or False

    Example:
    
    >>> from xplib.Annotation import Bed
    >>> from Annotation import overlap
    >>> bed1=Bed(["chr1",10000,12000,'-'])
    >>> bed2=Bed(["chr1",9000,13000,'+'])
    >>> print IsProperStrand(bed1,bed2)
    False

    """
    try:
        return (bed1.strand == bed2.strand) or (bed1.strand == '.') or (bed2.strand == '.')
    except:
        return True

def IsPartOf(bed1, bed2):
    """
    This function determines whether bed1 is part of bed2
    :param bed1: A Bed object from `xplib.Annotation.Bed <http://bam2xwiki.appspot.com/bed>`_ (BAM2X)
    :param bed2: A Bed object from `xplib.Annotation.Bed <http://bam2xwiki.appspot.com/bed>`_ (BAM2X)
    :returns: boolean -- True or False

    Example:
    
    >>> from xplib.Annotation import Bed
    >>> from Annotation import overlap
    >>> bed1=Bed(["chr1",10000,12000])
    >>> bed2=Bed(["chr1",9000,13000])
    >>> print IsPartOf(bed1,bed2)
    True

    This function allows N overhang nucleotides. 

    """
    N=5
    try:
        return ((bed1.stop <= bed2.stop + N) and (bed1.start >= bed2.start)) or ((bed1.stop <= bed2.stop) and (bed1.start >= bed2.start - N)) 
    except:
        return False


def Subtype(bed1,genebed,typ):
    """
    This function determines intron or exon or utr from a BED12 file.
    
    :param bed1: A Bed object defined by `xplib.Annotation.Bed <http://bam2xwiki.appspot.com/bed>`_ (BAM2X)
    :param genebed: A Bed12 object representing a transcript defined by xplib Annotaton.Bed with information of exon/intron/utr from an BED12 file
    :returns: str -- RNA subtype. "intron"/"exon"/"utr3"/"utr5"/"."

    Example:
    
    >>> from xplib.Annotation import Bed
    >>> from xplib import DBI
    >>> from Annotation import Subtype
    >>> bed1=Bed(["chr13",40975747,40975770])
    >>> a=DBI.init("../../Data/Ensembl_mm9.genebed.gz","bed")
    >>> genebed=a.query(bed1).next()
    >>> print Subtype(bed1,genebed)
    "intron"
        
    """
    N=0.85
    subtype="intron"
    if typ != "protein_coding":
        if overlap(bed1,genebed.utr3()):
            for i in genebed.utr3().Exons():
                if IsPartOf(bed1,i):
                    subtype="utr3"
        elif overlap(bed1,genebed.utr5()):
            for i in genebed.utr5().Exons():
                if IsPartOf(bed1,i):
                    subtype="utr5"
        else:
            for i in genebed.Exons():
                if IsPartOf(bed1,i):
                    subtype="exon"
                    break
    else:
        # print bed1
        # print genebed.utr3().start, genebed.utr3().stop
        if overlap(bed1,genebed.utr3()):
            # print "If passed"
            # print len(genebed.utr3().Exons())           
            for i in genebed.utr3().Exons():
                if overlap(bed1,i):
                    subtype="utr3"
        elif overlap(bed1,genebed.utr5()):
            for i in genebed.utr5().Exons():
                if overlap(bed1,i):
                    subtype="utr5"
        else:
            max_overlap=0
            flag=0
            for i in genebed.Exons():
                if IsPartOf(bed1,i):
                    subtype="exon"
                    flag=1
                    break
                elif overlap(bed1,i):
                    nt=min(bed1.stop,i.stop)-max(bed1.start,i.start)
                    if nt > max_overlap:
                        max_overlap = nt
            if flag==0 and max_overlap/float(bed1.stop-bed1.start)>=N:
                subtype="part_exon"

    return subtype

def optimize_annotation(c_dic,bed,ref_detail):
    '''
    This function will select an optimized annotation for the bed region from the genes in c_dic.

    It will select the annotation based on a list of priorities. The hypothesis is: pre-mature RNA -> mRNA+intron/nc transcript -> small functional RNA
    So the list of priorities is: exon/utr of coding transcript > exon of lincRNA > small RNA > exon/utr of nc transcript > intron of mRNA > intron of lincRNA
    ProperStrand > NonProperStrand > repeat masker (except rRNA_repeat according to the annotation files)
    '''

    #keep only one record for each type is enough
    #print c_dic
    ftyp=""
    fname=""
    fsubtype=""
    fstrandcol=""



    if "rRNA" in c_dic:
        ftyp=c_dic["rRNA"][0][0]
        fname=c_dic["rRNA"][0][1]
        fsubtype=c_dic["rRNA"][0][2]
        fstrandcol=c_dic["rRNA"][0][3]
        return [ftyp,fname,fsubtype,fstrandcol]


    if "protein_coding" in c_dic:
        for ind in xrange(len(c_dic["protein_coding"])-1,-1,-1):
            gene=c_dic["protein_coding"][ind]
      #      print gene
            flag=0
            pe_flag=0
            tmp=""
            for hit in ref_detail.query(bed):
                tempname=hit.id.split("&")
                #print tempname
                if gene[1]==tempname[0]:
                    gene[2]=Subtype(bed,hit,tempname[1])
                    if gene[2]!="intron":
                        if tempname[1]=="protein_coding" and gene[2]!="intron" and gene[2]!="part_exon":
                            #exon or utr of coding transcript
                            flag=1
                            break
                        elif tempname[1]!="protein_coding":
                            #exon or utr of non-coding transcript. Record the subtype. If the bed doesn't overlap with any exons, it wll be annotated as protein_coding-noncoding
                            tmp=gene[2]
                            flag=2
                        elif tempname[1]=="protein_coding" and gene[2]=="part_exon":
                            #the bed cover part of an exon. If it doesn't overlap with utr or exon of other transcript, it will be annotated as exon. 
                            pe_flag=1  
            # print flag,pe_flag
            # print gene[2]
            if flag==1 and gene[2]!="intron":  ##if gene type == protein_coding
                ftyp=gene[0]
                fname=gene[1]
                fsubtype=gene[2]
                fstrandcol=gene[3]
                break
            elif pe_flag==1:
                ftyp=gene[0]
                fname=gene[1]
                fsubtype="exon"
                fstrandcol=gene[3]
            elif flag==2:
                c_dic["protein_coding-noncoding"]=[["protein_coding-noncoding",gene[1],tmp,gene[3]]]
                ##it's fine to keep the same gene in the "protein_coding" list because intron has the lowest priority. All of the subtype should be intron.
            elif gene[2]==".":
                #if the bed is in intergenic region, remove this gene from the dictionary.
                c_dic["protein_coding"].remove(gene)
                if not c_dic["protein_coding"]:
                    c_dic.pop("protein_coding",None)
        del gene

    if "short_nc" in c_dic and ftyp=="":
        ftyp=c_dic["short_nc"][0][0]
        fname=c_dic["short_nc"][0][1]
        fsubtype=c_dic["short_nc"][0][2]
        fstrandcol=c_dic["short_nc"][0][3]
        return [ftyp,fname,fsubtype,fstrandcol]

 
    if "lincRNA" in c_dic and ftyp=="":
        for gene in c_dic["lincRNA"]:
            flag=0
            for hit in ref_detail.query(bed):
                if flag==0:
                    tempname=hit.id.split("&")
                    if gene[1]==tempname[0]:
                        gene[2]=Subtype(bed,hit,tempname[1])
                        if gene[2]!="intron":
                            flag=1
            if gene[2]!="intron":
                gene[2]="exon"
                ftyp=gene[0]
                fname=gene[1]
                fsubtype=gene[2]
                fstrandcol=gene[3]
                break
        del gene
        return [ftyp,fname,fsubtype,fstrandcol]

    if ftyp=="":
        if "other" in c_dic:
            gene=c_dic["other"][0]
        elif "protein_coding-noncoding" in c_dic:
            gene=c_dic["protein_coding-noncoding"][0]
        elif "protein_coding" in c_dic:
            gene=c_dic["protein_coding"][0]
        elif "lincRNA" in c_dic:
            gene=c_dic["lincRNA"][0]
        try:
            ftyp=gene[0]
            fname=gene[1]
            fsubtype=gene[2]
            fstrandcol=gene[3]
        except:
            pass

    return [ftyp,fname,fsubtype,fstrandcol]

def annotation(bed,ref_allRNA,ref_detail,ref_repeat):
    """
    This function is based on :func:`overlap` and :func:`optimize_annotation` and :func:`Subtype` functions to annotate RNA type/name/subtype for any genomic region.
    This function will first find genes with maximum overlap (+/-2 nt) with bed, and use the function optimize_annotation to select an optimized annotation for the bed. 

    >Find hits with overlaps larger than Perc_overlap of the bed region length and build dic
    >Find hits with overlaps between Perc_max-1.0 * max_overlap and build P_dic, N_dic
    >Find an annotation

    :param bed: A Bed object defined by `xplib.Annotation.Bed <http://bam2xwiki.appspot.com/bed>`_ (in BAM2X).
    :param ref_allRNA: the `DBI.init <http://bam2xwiki.appspot.com/DBI>`_ object (from BAM2X) for bed6 file of all kinds of RNA
    :param ref_detail: the `DBI.init <http://bam2xwiki.appspot.com/DBI>`_ object for bed12 file of lincRNA and mRNA with intron, exon, UTR
    :param ref_detail: the `DBI.init <http://bam2xwiki.appspot.com/DBI>`_ object for bed6 file of mouse repeat
    :returns: list of str -- [type,name,subtype, strandcolumn]

    Example:

    >>> from xplib.Annotation import Bed
    >>> from xplib import DBI
    >>> from Annotation import annotation
    >>> bed=Bed(["chr13",40975747,40975770])
    >>> ref_allRNA=DBI.init("../../Data/all_RNAs-rRNA_repeat.txt.gz","bed")
    >>> ref_detail=DBI.init("../../Data/Ensembl_mm9.genebed.gz","bed")
    >>> ref_repeat=DBI.init("../../Data/mouse.repeat.txt.gz","bed")
    >>> print annotation(bed,ref_allRNA,ref_detail,ref_repeat)
    ["protein_coding","gcnt2","intron","ProperStrand"]

    """

    Perc_overlap=0.7
    Perc_max=0.85

    flag=0
    typ = "non"
    name = "."
    subtype = "."
    strandcol = "ProperStrand"

    ftyp = ""
    fname = ""
    fsubtype = ""
    fstrandcol = ""

    bed_len=bed.stop-bed.start+1
 
    max_overlap = 0  # find annotation with largest overlap
    overlap_dic={} #key: overlap length element: list of genes
    P_dic={}  #dictionary of properstrand. Key: type
    N_dic={}  #dictionary of nonproperstrand. Key type

    ##construct overlap_dic
    for hit in ref_allRNA.query(bed):
        overlap = min(hit.stop,bed.stop)-max(hit.start,bed.start)
        if overlap >= Perc_overlap * bed_len and overlap!=0:
            name=hit.id.split(".",1)[1]
            typ=hit.id.split(".")[0]
            if not IsProperStrand(hit, bed):
                strandcol = "NonProperStrand"
            else:
                strandcol = "ProperStrand"

            if overlap not in overlap_dic:
                overlap_dic[overlap]=[]
            overlap_dic[overlap].append([typ,name,subtype,strandcol])

            if overlap > max_overlap:
                max_overlap = overlap

    ##construct P_dic and N_dic
    #print overlap_dic
    for key in overlap_dic.keys():
        # print key
        if key >= max_overlap * Perc_max:
            for gene in overlap_dic[key]:
                # print gene
                typ = gene[0]
                name = gene[1]
                subtype = gene[2]
                strandcol = gene[3]
                if strandcol == "ProperStrand":
                    if typ=="protein_coding" or typ=="lincRNA":
                        if typ in P_dic:
                            P_dic[typ].append([typ,name,subtype,strandcol])
                        else:
                            P_dic[typ]=[[typ,name,subtype,strandcol]]
                    elif typ=="rRNA_repeat" or typ=="rRNA":
                        if "rRNA" in P_dic:
                            P_dic["rRNA"].append([typ,name,subtype,strandcol])
                        else:
                            P_dic["rRNA"]=[[typ,name,subtype,strandcol]]
                    elif typ=="snoRNA" or typ=="miRNA" or typ=="snRNA":
                        if "short_nc" in P_dic:
                            P_dic["short_nc"].append([typ,name,subtype,strandcol])
                        else:
                            P_dic["short_nc"]=[[typ,name,subtype,strandcol]]
                    else:
                        if "other" in P_dic:
                            P_dic["other"].append([typ,name,subtype,strandcol])
                        else:
                            P_dic["other"]=[[typ,name,subtype,strandcol]]
                elif strandcol == "NonProperStrand":
                    if typ=="protein_coding" or typ=="lincRNA":
                        if typ in N_dic:
                            N_dic[typ].append([typ,name,subtype,strandcol])
                        else:
                            N_dic[typ]=[[typ,name,subtype,strandcol]]
                    elif typ=="rRNA_repeat" or typ=="rRNA":
                        if "rRNA" in N_dic:
                            N_dic["rRNA"].append([typ,name,subtype,strandcol])
                        else:
                            N_dic["rRNA"]=[[typ,name,subtype,strandcol]]
                    elif typ=="snoRNA" or typ=="miRNA" or typ=="snRNA":
                        if "short_nc" in N_dic:
                            N_dic["short_nc"].append([typ,name,subtype,strandcol])
                        else:
                            N_dic["short_nc"]=[[typ,name,subtype,strandcol]]
                    else:
                        if "other" in N_dic:
                            N_dic["other"].append([typ,name,subtype,strandcol])
                        else:
                            N_dic["other"]=[[typ,name,subtype,strandcol]]

    ##select optimized annotation
    if P_dic:
        [ftyp,fname,fsubtype,fstrandcol] = optimize_annotation(P_dic,bed,ref_detail)

    if ftyp=="" and N_dic:
          [ftyp,fname,fsubtype,fstrandcol] = optimize_annotation(N_dic,bed,ref_detail)

    if ftyp=="":
        max_overlap=0
        #typ=="non" try repeat masker
        #we are not using any stringent threshold here. According to the annotation, different types of repeat element, such as LINE and SINE, are (usually) exclusive. 
        #For example, if one element is annotated as LINE, it won't be SINE at the same time.
        for hit in ref_repeat.query(bed):
            overlap = min(hit.stop,bed.stop)-max(hit.start,bed.start)
            if overlap > max_overlap and overlap >= Perc_overlap * bed_len:
                max_overlap=overlap
                tempname=hit.id.split("&")
                name = tempname[0]
                typ = tempname[1]
                subtype = tempname[2]
                if not IsProperStrand(hit, bed):
                    strandcol = "NonProperStrand"
                else:
                    strandcol = "ProperStrand"
        if max_overlap>0:
            ftyp=typ
            fname=name
            fsubtype=subtype
            fstrandcol=strandcol
        else:
            ftyp="non"
            fname="."
            fsubtype="."
            fstrandcol="ProperStrand"

    return [ftyp,fname,fsubtype,fstrandcol]

        ##find a 

        # '''
        # try:
        #     tran=genome.getTranscriptByStableId(StableId=tempname).Gene
        #     typ=tran.BioType
        #     name=tran.Symbol
        # except: pass
        # '''
    # if (typ=="protein_coding" and subtype=="."):        # this means the hit is actually in a intergenic region (rare)
    #     typ = "non"
    #     name = "."

        # '''
        # try:
        #     repeats=genome.getFeatures(CoordName=bed.chr[3:], Start=bed.start, End=bed.stop, feature_types='repeat')
        #     for r in repeats:
        #        if r.RepeatClass!='dust':
        #            typ=r.RepeatType
        #            name=r.Symbol
        #            break 
        # except: pass
        # '''
    # if typ=="lincRNA" and subtype!="intron":
    #     subtype="exon"
    # if typ in ["snoRNA","snRNA","miRNA","miscRNA"]:
    #     subtype='.'
    #if typ=="pseudogene":
      #  subtype="."
    #if typ=="SINE" and subtype=="Alu":
    #    subtype="B1"

                    
           
