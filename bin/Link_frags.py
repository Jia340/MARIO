import argparse,sys

def ParseArg():
  p=argparse.ArgumentParser( description = "Link fragment1 and 2 by linker ID. Fragment will be output only when both RNA1 and 2 are uniquely mapped." )
  p.add_argument("--RNA1",type=str,required=True,help="merged RNA1")
  p.add_argument("--RNA2",type=str,required=True,help="merged RNA2")
  p.add_argument("-o","--output",type=str,default="linkedFrags.txt",help="Output file name.")
  if len(sys.argv)==1:
    print >>sys.stderr, p.print_help()
    sys.exit(0)
  return p.parse_args()  

def Main():
  args=ParseArg()
  Frags={}
  file1=open(args.RNA1,"r")
  output=open(args.output,"w")
  while True:
    line=file1.readline()
    if line=="":break
    line=line.strip().split('\t')
    if line[11]=="One":
      if line[10] not in Frags:
        Frags[line[10]]=['\t'.join(line[0:11])]
      else:
        print >>sys.stderr,"Error: RNA1"
        print >>sys.stderr,line
        exit(0)

  file2=open(args.RNA2,"r")  
  while True:
    line=file2.readline()
    if line=="":break
    line=line.strip().split('\t')
    if line[11]=="One" and line[10] in Frags:
      Frags[line[10]].append('\t'.join(line[0:10]))

  for key in Frags:
    if len(Frags[key])==2:
      print >>output,Frags[key][0]+'\t'+Frags[key][1]
    elif len(Frags[key])>2:
      print >>sys.stderr,"Error: RNA2"
      print >>sys.stderr,Frags[key]
      exit(0) 
  
if __name__=="__main__":
  Main()
