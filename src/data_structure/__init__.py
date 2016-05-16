from xplib.Annotation import Bed
from AnnoMax import *

class annotated_bed():
    """
    To store, compare, cluster for the genomic regions with RNA annotation information. Utilized in the program :ref:`Select_stronginteraction_pp.py<Step6>`
    """
    
    def __init__(self,x=None,**kwargs):
        """
        Initiation example:
        
        >>> str="chr13\t40975747\t40975770\t+\tATTAAG...TGA\tprotein_coding\tgcnt2\tintron"
        >>> a=annotated_bed(str)
        or
        >>> a=annotated_bed(chr="chr13",start=40975747,end=40975770,strand='+',type="protein_coding",)
        """
        self.annotated = False
        if x is not None:
            if type(x)==type("str"):
                x=x.split("\t")
            self.chr=str(x[0])
            if "," not in x[1]:
                self.start=int(x[1])
                self.end=int(x[2])
            else:
                self.start=[int(f) for f in x[1].split(",")]
                self.end=[int(f) for f in x[2].split(",")]
            self.strand=str(x[3]).strip()
            try:
                self.seq=str(x[4]).strip()
                self.source=str(x[5]).strip()
                self.type=str(x[6]).strip()
                self.name=str(x[7]).strip()
                self.subtype=str(x[8]).strip()
                self.proper=str(x[9]).strip()
                self.id=''
#                self.source=str(x[5]).strip()
            except:
                pass
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])
    def overlap(self,other): # for the purpose of finding overlap between regions
        """
        Find overlap between regions

        :param other: another :class:`annotated_bed` object
        :returns: boolean

        """
          
        return ((self.chr==other.chr)&(self.end>other.start)&(self.start<other.end))

    def __lt__(self, other): # for the purpose of ordering clusters
        """
        Compare two objects self and other when they are not **overlapped**
        
        :param other: another :class:`annotated_bed` object
        :returns: boolean -- "None" if overlapped.

        Example:
        
        >>> a=annotated_bed(chr="chr13",start=40975747,end=40975770)
        >>> b=annotated_bed(chr="chr13",start=10003212,end=10005400)
        >>> print a>b
        False

        """
        if self.overlap(other): return "None"
        if self.chr==other.chr:
            return ((self.start<other.start)&(self.end<other.end))
        else:
            return (self.chr<other.chr)

    def Cluster(self,c):
        """
        Store cluster information of self object
        
        :param c: cluster index
        
        Example:
        
        >>> a=annotated_bed(chr="chr13",start=40975747,end=40975770)
        >>> a.Cluster(3)
        >>> print a.cluster
        3

        .. Note::
        
           a.cluster will be the count information when a become a cluster object in :ref:`Select_stronginteraction_pp.py<Step6>`

        """
        self.cluster=c
        #self.c_info=c_info
    def Update(self,S,E): # used for update or expand cluster locations
        """
        Update the upper and lower bound of the cluster after adding segments using Union-Find.
        
        :param S: start loc of the newly added genomic segment
        :param E: end loc of the newly added genomic segment
        
        Example:
        
        >>> a=annotated_bed(chr="chr13",start=40975747,end=40975770)
        >>> a.Update(40975700,40975800)
        >>> print a.start, a.end
        40975700 40975800
        """

        self.start=min(self.start,S)
        self.end=max(self.end,E)
    def Annotate(self,dbi1,dbi2,dbi3): # do RNA annotation if not annotated
        """
        Update annotation.

        :param ref_allRNA: the `DBI.init <http://bam2xwiki.appspot.com/DBI>`_ object (from BAM2X) for bed6 file of all kinds of RNA
        :param ref_detail: the `DBI.init <http://bam2xwiki.appspot.com/DBI>`_ object for bed12 file of lincRNA and mRNA with intron, exon, UTR   
        :param ref_detail: the `DBI.init <http://bam2xwiki.appspot.com/DBI>`_ object for bed6 file of mouse repeat

        Example:
        
        >>> str="chr13\t40975747\t40975770\t+"
        >>> a=annotated_bed(str)
        >>> a.Cluster(3)
        >>> ref_allRNA=DBI.init("../../Data/all_RNAs-rRNA_repeat.txt.gz","bed")
        >>> ref_detail=DBI.init("../../Data/Ensembl_mm9.genebed.gz","bed")
        >>> ref_repeat=DBI.init("../../Data/mouse.repeat.txt.gz","bed")
        >>> a.Annotate(ref_allRNA,ref_detail,ref_repeat)
        >>> print a
        "chr13\t40975747\t40975770\tprotein_coding\tgcnt2\tintron\t3"
        """
        if not self.annotated:
            if "chr" in self.chr:
                bed=Bed([self.chr,self.start,self.end])
                [self.type,self.name,self.subtype,self.proper]=annotation(bed,dbi1,dbi2,dbi3)
            self.annotated = True

    def __str__(self):
        """
        Use print function to output the cluster information (chr, start, end, type, name, subtype,cluster)
        
        Example:
        
        >>> str="chr13\t40975747\t40975770\t+\tATTAAG...TGA\tprotein_coding\tgcnt2\tintron"
        >>> a=annotated_bed(str)
        >>> a.Cluster(3)
        >>> a.Update(40975700,40975800)
        >>> print a
        "chr13\t40975700\t40975800\tprotein_coding\tgcnt2\tintron\t3"
        """
        return "\t".join(str(f) for f in [self.chr,self.start,self.end,self.type,self.name,self.subtype,self.cluster])
    # self.cluster become number of regions within cluster for cluster_pool object

class annotated_bed_proper(annotated_bed):

    def __str__(self):
        """
        Use print function to output the cluster information (chr, start, end, type, name, subtype,cluster)
        
        Example:
        
        >>> str="chr13\t40975747\t40975770\t+\tATTAAG...TGA\tprotein_coding\tgcnt2\tintron"
        >>> a=annotated_bed(str)
        >>> a.Cluster(3)
        >>> a.Update(40975700,40975800)
        >>> print a
        "chr13\t40975700\t40975800\tprotein_coding\tgcnt2\tintron\t3"
        """
        return "\t".join(str(f) for f in [self.chr,self.start,self.end,self.type,self.name,self.subtype,self.id,self.proper,self.cluster])
    # self.cluster become number of regions within cluster for cluster_pool object
