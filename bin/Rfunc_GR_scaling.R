library(parallel)
library(rtracklayer)
library(GenomicRanges)
library(megadepth)
divide_coords=function(Start, End, Bins=100){
  #Takes a start and end coordinates and scale them in [Bins] bins.
  #Outputs a 2-columns table with starts and ends coordinates of each bin.
  if(length(Start)!=1 | length(End) !=1 | length(Bins) !=1){stop("Arguments length not valids")}
  if(Start <=0 | End <=Start | Bins < 1){stop("Arguments not valids (start <0 or End <=start or Bins<1")}
  Step=(End-Start+1)/Bins
  starts=c()
  ends=c()
  end_bin=Start-1 # initializing with the first coordinates -1 because we add +1 at each iteration
  for(i in 1:Bins){
    st_bin=end_bin+1 # New stats = previous end +1
    end_bin=Start+floor(i*Step)-1 # New end = real start + nb step -1
    starts=c(starts, st_bin)
    ends=c(ends, end_bin)
  }
  return(cbind(starts, ends))
}
divide_GR = function(GRobj, Bins=100, Before=0, After=0){
  #Take a GRange object (bed like) and output it with scaled coordinates : 1 range per bin
  if(Before <0 | After < 0 | Bins < 1){stop("Arguments not valids")}
  if(length(GRobj) >1 ){stop("GRobj must be of 1 length")}
  r=divide_coords(start(GRobj)-Before, end(GRobj)+After, Bins)
  return=GRanges(
    seqnames=seqnames(GRobj),
    ranges=apply(r, 1, function(x) paste0(x[1], "-", x[2])),
    name=GRobj$name, 
    score=GRobj$score,
    strand=strand(GRobj))
}
make_scaled_coords=function(GRobj, FinalBins=100, Extention=FALSE, 
                            Ext_length=c(300, 300), Ext_value=c(10,10)){
  #Takes a single  GRange object as input (1 interval) and bin it (FinalBins)
  #Outputs a GRange object with scaled intervals, 1 interval per bin.
  #     - allows to extend the original coordinates and scale them differently
  if(!Extention){
    return(divide_GR(GRobj, Bins=FinalBins))
  }else {
    MidBins=FinalBins-sum(Ext_value)
    res=list()
    res[[1]]=divide_GR(GRobj = GRanges(seqnames=seqnames(GRobj),
                                       ranges=paste0(start(GRobj[1])-Ext_length[1], "-", start(GRobj[1])-1),
                                       name=GRobj$name, 
                                       score=GRobj$score,
                                       strand=strand(GRobj)), Ext_value[1])
    res[[2]]=divide_GR(GRobj, Bins=MidBins)
    res[[3]]=divide_GR(GRobj = GRanges(seqnames=seqnames(GRobj),
                                       ranges=paste0(end(GRobj[1])+1, "-", end(GRobj[1])+Ext_length[2]),
                                       name=GRobj$name, 
                                       score=GRobj$score,
                                       strand=strand(GRobj)), Ext_value[2])
    return(unlist(as(res, "GRangesList")))
  }
}


get_bw_score = function( bw_file, bed_file, gr_bed){
  #bw_file bigwig file from which the coverage is taken
  #bed_file = a detailed bed file, usually with multiple interval per gene/features
  #gr_bed = a GenomicRange object from the initial bed file, with 1 interval per feature.
  coverage=get_coverage(bw_file, op="mean", annotation=bed_file)
  data=c()
  for(i in 1:length(gr_bed)){
    query=gr_bed[i];
    ovl=to(findOverlaps(query=query, subject=coverage))
    if(length(ovl)>=1){
      data=rbind(data, coverage[ovl]$score)
      names(data)[length(data)]=gr_bed[i]$name
    }
  }
  #rownames(data)=gr_bed$name
  return(data)
}
