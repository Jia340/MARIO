#!/usr/bin/env Rscript

list.of.packages <- c("argparse","RCircos", "biovizBase","rtracklayer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages,repos="http://cran.us.r-project.org");
  source("http://bioconductor.org/biocLite.R")
  biocLite(new.packages)
}
suppressPackageStartupMessages(require(rtracklayer))
suppressPackageStartupMessages(require(RCircos))
suppressPackageStartupMessages(require(argparse))

# creat parser object
parser <- ArgumentParser()

# specify options
parser$add_argument("interaction",help="the interaction file,[required]")
parser$add_argument("-g","--genome",default="mm9",help="genome information, choice: mm9/mm10/hg19 et.al., [default: %(default)s]")
parser$add_argument("part1",help="aligned BAM file for part1,[required]")
parser$add_argument("part2",help="aligned BAM file for part2,[required]")
parser$add_argument("-b","--bin",type="integer",default=100000,help="window size for the bins, [default: %(default)s]")
parser$add_argument("-o","--output",default="Interactome_view.pdf",help="output pdf file name, [default: %(default)s]")



#  modified Histogram function add color option
RCircos.Histogram.Plot<-function(hist.data, data.col, track.num, side, col)
{
  RCircos.Pos <- RCircos.Get.Plot.Positions();
  RCircos.Par <- RCircos.Get.Plot.Parameters();
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
  
  
  #       Convert raw data to plot data. The raw data will be validated
  #       first during the convertion
  #       ____________________________________________________________
  #       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  hist.data <- RCircos.Get.Plot.Data(hist.data, "plot");
  
  
  #       Each heatmap cell is centered on data point location. Make
  #       sure each one will be in the range of it chromosome
  #       ____________________________________________________________
  #       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  hist.locations <- as.numeric(hist.data[, ncol(hist.data)]);
  start <- hist.locations - RCircos.Par$hist.width;
  end   <- hist.locations + RCircos.Par$hist.width;
  
  data.chroms <- as.character(hist.data[,1]);
  chromosomes <- unique(data.chroms);
  cyto.chroms <- as.character(RCircos.Cyto$Chromosome);
  
  for(a.chr in 1:length(chromosomes))
  {
    cyto.rows <- which(cyto.chroms==chromosomes[a.chr]);
    locations <- as.numeric(RCircos.Cyto$Location[cyto.rows]);
    chr.start <- min(locations) - RCircos.Cyto$Unit[cyto.rows[1]];
    chr.end   <- max(locations);
    
    data.rows <- which(data.chroms==chromosomes[a.chr]);
    start[data.rows[start[data.rows]<chr.start]] <- chr.start;
    end[data.rows[end[data.rows]>chr.end]] <- chr.end;
  }
  
  
  #       Plot position for current track. 
  #       ___________________________________________________________
  #       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  locations <- RCircos.Track.Positions(side, track.num);
  out.pos <- locations[1];
  in.pos <- locations[2];
  
  
  #       Draw histogram [modified]
  #       ___________________________________________________________
  #       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  num.subtrack <- RCircos.Par$sub.tracks;
  RCircos.Track.Outline(out.pos, in.pos, RCircos.Par$sub.tracks, NULL);
  
  for(a.point in 1:nrow(hist.data))
  {
    hist.height <- hist.data[a.point, data.col];
    
    the.start <- start[a.point];
    the.end <- end[a.point];
    
    #       Plot rectangle with specific height for each 
    #       data point
    #       _______________________________________________
    #       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    height <- in.pos + RCircos.Par$track.height*hist.height;
    
    polygon.x<- c(RCircos.Pos[the.start:the.end,1]*height, 
                  RCircos.Pos[the.end:the.start,1]*in.pos);
    polygon.y<- c(RCircos.Pos[the.start:the.end,2]*height, 
                  RCircos.Pos[the.end:the.start,2]*in.pos);
    polygon(polygon.x, polygon.y, col=col, border=NA);
  }
  
}


# modified colorRampPalette function to include alpha channel as well
# from http://r.789695.n4.nabble.com/Alpha-channel-in-colorRamp-and-colorRampPalette-td4672229.html
colorRampPalette <- function (colors, ...) { 
  ramp <- colorRamp(colors, ...) 
  function(n) { 
    x <- ramp(seq.int(0, 1, length.out = n)) 
    rgb(x[, 1], x[, 2], x[, 3], x[, 4], maxColorValue = 255) 
  } 
} 

colorRamp <- function (colors, bias = 1, space = c("rgb", "Lab"), interpolate = c("linear", "spline")) { 
  if (bias <= 0) 
    stop("'bias' must be positive") 
  colors <- t(col2rgb(colors, alpha=T)/255) 
  space <- match.arg(space) 
  interpolate <- match.arg(interpolate) 
  if (space == "Lab") { 
    colors <- convertColor(colors, from = "sRGB", to = "Lab") 
  } 
  interpolate <- switch(interpolate, linear = stats::approxfun, spline = stats::splinefun) 
  if ((nc <- nrow(colors)) == 1L) { 
    colors <- colors[c(1L, 1L), ] 
    nc <- 2L 
  } 
  x <- seq.int(0, 1, length.out = nc)^bias 
  palette <- c(interpolate(x, colors[, 1]), interpolate(x, colors[, 2]), interpolate(x, colors[, 3]), interpolate(x, colors[, 4])) 
  roundcolor <- function(rgb) pmax(pmin(rgb, 1), 0) 
  if (space == "Lab") { 
    function(x) { 
      roundcolor(convertColor(cbind(palette[[1L]](x), palette[[2L]](x), 
                                    palette[[3L]](x), palette[[4L]](x)), from = "Lab", to = "sRGB")) * 
        255 
    } 
  } 
  else { 
    function(x) { 
      roundcolor(cbind(palette[[1L]](x), palette[[2L]](x), 
                       palette[[3L]](x), palette[[4L]](x))) * 255 
    } 
  } 
} 


#  modified link plot function change color scheme (color_list is the vector of inter-type and p-values for all interactions)
RCircos.Link.Plot<-function(link.data, track.num, color_list)
{
  RCircos.Pos <- RCircos.Get.Plot.Positions();
  RCircos.Par <- RCircos.Get.Plot.Parameters();
  
  
  #       Check chromosome names, Start, and End positions
  #       _______________________________________________________________
  #       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  link.data <- RCircos.Validate.Genomic.Data(link.data, plot.type="link");
  
  
  #       Plot position for link track.
  #       __________________________________________________________
  #       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  one.track <- RCircos.Par$track.height + RCircos.Par$track.padding;
  start <- RCircos.Par$track.in.start - (track.num-1)*one.track;
  base.positions <- RCircos.Pos*start;
  
  data.points <- matrix(rep(0, nrow(link.data)*2), ncol=2);
  colors <- cbind(rep("#FFFFFF11",16),rainbow(16))
  line.colors <- apply(color_list,1, function (x) colorRampPalette(colors[x[1],])(12)[floor((18-x[2])/4)]);
  for(a.link in 1:nrow(link.data))
  {
    if (color_list[a.link,1]==16) { next; }   #  remove "Other" types 
    data.points[a.link, 1] <- RCircos.Data.Point(
      link.data[a.link, 1], link.data[a.link, 2]);
    data.points[a.link, 2] <- RCircos.Data.Point(
      link.data[a.link, 4], link.data[a.link, 5]);
    
    if(data.points[a.link, 1]==0 || data.points[a.link, 2]==0)
    {  print("Error in chromosome locations ...");  break; }
    
  }
  
  
  #       Draw link lines for each pair of locations
  #       __________________________________________________
  #       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  for(a.link in 1:nrow(data.points))
  {  
    point.one <- data.points[a.link, 1];
    point.two <- data.points[a.link, 2];
    if(point.one > point.two)
    { 
      point.one <- data.points[a.link, 2];
      point.two <- data.points[a.link, 1];
    }
    
    P0 <- as.numeric(base.positions[point.one,]);
    P2 <- as.numeric(base.positions[point.two,]);
    links <- RCircos.Link.Line(P0, P2); 
    lines(links$pos.x, links$pos.y, type="l", 
          col=line.colors[a.link] ); 
  }
}

RCircos.Link.Line<-function(line.start, line.end)
{       
  #       Set up the points for Bezure curve
  #       ___________________________________________
  #       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  P0 <- line.start;
  P2 <- line.end;
  
  
  #       Calculate total number of points for 
  #       the Bezuer curve
  #       ___________________________________________
  #       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  bc.point.num <- 1000;
  t <- seq(0, 1, 1/bc.point.num);
  
  
  #       Calculate the point values for Bezuer curve
  #       ___________________________________________
  #       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  link.x <- (1-t)^2*P0[1] + t^2*P2[1];
  link.y <- (1-t)^2*P0[2] + t^2*P2[2];
  
  
  #       Return the coordinates
  #       ___________________________________________
  #       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  return (list(pos.x=link.x, pos.y=link.y));
}



# get the input
args <- parser$parse_args()
interaction_f <- args$interaction
genome <- args$genome
part1 <- args$part1
part2 <- args$part2
bin <- args$bin
output <- args$output

# get cytoband ideogram information
cat ("Generating cytoband information from genome version... \n")
if (genome=="mm10") {
  data(UCSC.Mouse.GRCm38.CytoBandIdeogram)
  cytoband=UCSC.Mouse.GRCm38.CytoBandIdeogram
  rm(UCSC.Mouse.GRCm38.CytoBandIdeogram)
} else if (genome=="hg19") {
  data(UCSC.HG19.Human.CytoBandIdeogram)
  cytoband=UCSC.HG19.Human.CytoBandIdeogram
  rm(UCSC.HG19.Human.CytoBandIdeogram)
} else {
  suppressPackageStartupMessages(require(biovizBase))
  cytoband=as.data.frame(getIdeogram(genome,cytobands=T))
  cytoband=cytoband[,c(1:3,6:7)]
  colnames(cytoband)=c("Chromosome","ChromStart","ChromEnd","Band","Stain")
}


# get coverage for BAM files of part1 and part2
cat("\nCounting coverage (RPKM) for part1 and part2 alignment:")
commend <- paste("bam2tab.py --bin",bin,"-s",part1,part2,"> temp_coverage.txt",sep=" ")
system(commend)
coverage=read.table("temp_coverage.txt",header=F)
coverage=coverage[coverage[,1]!="chrM",]
# for windows exceed the chromosome boundary
for (chr in unique(as.character(coverage[,1]))) {
  n=max((1:dim(coverage)[1])[coverage[,1]==chr])
  m=max((1:dim(cytoband)[1])[cytoband$Chromosome==chr])
  coverage[n,3]=min(coverage[n,3],cytoband[m,3])
}

part1_cov=coverage[,1:4]
part2_cov=coverage[,c(1:3,5)]
colnames(part1_cov)=c("Chromosome", "chromStart", "chromEnd", "Data")
colnames(part2_cov)=c("Chromosome", "chromStart", "chromEnd", "Data")
system("rm temp_coverage.txt")
part1_cov$Data=log(part1_cov$Data+1)/max(log(part1_cov$Data+1))
part2_cov$Data=log(part2_cov$Data+1)/max(log(part2_cov$Data+1))

# get interaction file ready
cat("Generating link data from interaction file... \n")
interaction<-read.table(interaction_f,sep='\t',header=F)
interaction=interaction[(!grepl("rRNA",interaction[,4]))&(!grepl("rRNA",interaction[,11])),]
interaction=interaction[(interaction[,1]!="chrM")&(interaction[,4]!="chrM"),]
#  interaction_type
Types=apply(interaction[,c(4,11)],1, function(x) paste(sort(x)[1],sort(x)[2],sep='-'))
type_n = c("snRNA","lincRNA","snoRNA","miRNA","protein_coding","pseudogene")
Types[!((interaction[,4] %in% type_n)&(interaction[,11] %in% type_n))]="Other"
Types[(as.character(interaction[,1])==as.character(interaction[,8]))&(as.character(interaction[,5])==as.character(interaction[,12]))]="Self"

Types[Types=="miRNA-snRNA"]="Other"
Types[Types=="miRNA-miRNA"]="Other"
cat("Count of different interaction types: \n")
print(table(Types))

All_types = c("lincRNA-lincRNA","lincRNA-protein_coding","lincRNA-snoRNA","lincRNA-snRNA","miRNA-protein_coding","miRNA-snoRNA","protein_coding-protein_coding","protein_coding-snoRNA","protein_coding-snRNA","snoRNA-snoRNA","snoRNA-snRNA","snRNA-snRNA","pseudogene-snoRNA","protein_coding-pseudogene","Self","Other")
Type_int = match(Types,All_types)
Type_int = Type_int[(Types!="Other")&(!is.na(Type_int))]
print(table(Type_int))
interaction = interaction[Types!="Other",]
interaction_p=interaction[,16] # for the alpha channel of color, based on -log(p-value)
interaction_p[interaction_p==-Inf]=-30.0
interaction=interaction[,c(1:3,8:10)]
data(RCircos.Link.Data)
colnames(interaction)<-colnames(RCircos.Link.Data)
rm(RCircos.Link.Data)


#  Setup RCircos core components:
RCircos.Set.Core.Components(cytoband, NULL, 2, 0);

#  Open the graphic device (here a pdf file)
cat("Open graphic device and start plot ...\n");
pdf(file=output, height=8, width=8);
RCircos.Set.Plot.Area();
title("RNA-RNA interactome");

# 	Draw chromosome ideogram
cat("Draw chromosome ideogram ...\n");
RCircos.Chromosome.Ideogram.Plot();

#  Histogram plot for part1 and part2
cat("Draw Histogram of part1/2 ...\n");
RCircos.Histogram.Plot(part1_cov, 4, 1, "in","#763a7a");
RCircos.Histogram.Plot(part2_cov, 4, 2, "in","#0288ad");

#  Link plot for interactions
cat("Draw Interaction links ...\n");
RCircos.Link.Plot(interaction,3,cbind(Type_int,interaction_p))   # Type_int for color and interaction_p for alpha

#  Close the graphic device and clear memory
dev.off()
cat("Circos interaction plotting done ...\n\n");
rm(list=ls(all=T));

