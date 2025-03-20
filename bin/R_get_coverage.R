### With the new duplicate name checking, need to verify if in the coverage per bin we still need to make the verification.
###
params=c(BedExtLengthLeft, BedExtLengthRight, BedRFinalLength, BedExtValLeft, BedExtValRight)
for(p in params){if(is.null(p) | !is.numeric(p) | p<0 | !is.finite(p) | is.na(p)){p=0}}
BedExtLengthLeft=params[1]; BedExtLengthRight=params[2]; BedRFinalLength=params[3]; BedExtValLeft=params[4]; BedExtValRight=params[5]
if(BedExtLengthLeft==0 & BedExtLengthRight==0){BedExtension=FALSE}
if(BedExtValLeft==0 & BedExtValRight==0){BedCustomScaling=FALSE}
if(BedRFinalLength==0){BedRFinalLength=10}



#Importing bed 
bed_gr=import(BedFile)

# Adding names to the intervals if not present initially
if(is.null(bed_gr$name)){
    bed_gr$name=paste0("feature_", c(1:length(bed_gr)))
} else { 
    nb_dup=sum(duplicated(bed_gr$name))
    # Editing names so that there are no duplicates
    if(nb_dup>0){
        cat(paste0("##ATTENTION, ", nb_dup, " features were carring duplicated names. Those will be renamed uniquely\n"))
        for(i in unique(bed_gr$name[duplicated(bed_gr$name)])){a=which(bed_gr$name==i); for(j in 2:length(a)){bed_gr$name[a[j]]=paste(i, j, sep="-")}} 
    }
}
# Adding a 0 score if not present initially
if(is.null(bed_gr$score)){bed_gr$score=0}


#Creating GR bed object with extended coordinates
bed_ext_gr=bed_gr
if(BedExtension){
    start(bed_ext_gr)=start(bed_gr)-1-BedExtLengthLeft
    if(sum(start(bed_ext_gr)<0)>0){paste0("Attention, at least ", sum(start(bed_ext_gr)<0), ' invervals are out of chromosomal range (start)')}
    start(bed_ext_gr)[start(bed_ext_gr)<0] = 0
    end(bed_ext_gr)=end(bed_gr)+BedExtLengthRight
}
bed_ext_gr=bed_ext_gr[start(bed_ext_gr)>0]
bed_ext_fname=paste0(Output_prefix,".ext.bed")
export.bed(bed_ext_gr,con=bed_ext_fname)

#Getting average density per interval
cl=makeCluster(Threads)
clusterEvalQ(cl, {library(GenomicRanges); library(megadepth)})
clusterExport(cl, c("bw_fnames", "bw_names", "bed_ext_fname","bed_ext_gr"))
mean_per_interval <- parLapply(cl, bw_fnames , function(x) {  
get_coverage(x, op="mean", annotation=bed_ext_fname)$score# Megadepth function to get mean coverage per interval
})
stopCluster(cl)
names(mean_per_interval)=bw_names
tab_mean_per_interval=do.call(cbind, mean_per_interval)
rownames(tab_mean_per_interval)=bed_ext_gr$name
write.table(x=tab_mean_per_interval, file=paste0(Output_prefix, '.avg_density.tsv'), quote=FALSE, row.names=TRUE, 
col.names=TRUE, sep="\t")

cat(paste("#### Bed bined start", date(), "####", "", sep="\n"))

#Creating GR bined bed object with extended coordinates according to configuration (Takes a bit of time)+
system(paste0("ln -s ",RfunctionFile, " MyFunctionsToLoad.R "))
cl=makeCluster(Threads)
clusterEvalQ(cl, source("MyFunctionsToLoad.R"))
clusterExport(cl, c("bed_gr","bed_ext_gr","BedRFinalLength","BedExtLengthLeft", "BedExtension",
                    "BedExtLengthRight", "BedExtValLeft","BedExtValRight"))
if(BedCustomScaling){
    # In case of custom scaling, we use the initial bed and scales the extentions differently
    # If there is no bed extension, the custom scaling will not change anything.
    res_bin_bed <- parLapply(cl, c(1:length(bed_gr)), function(y) 
    make_scaled_coords(bed_gr[y], FinalBins = BedRFinalLength,
                Extention = BedExtension, Ext_length = c(BedExtLengthLeft, BedExtLengthRight),
                Ext_value = c(BedExtValLeft,BedExtValRight)))
}
if(!BedCustomScaling){
    # In case of NO custom scaling, we use the extended bed file and scales the entire extended interval as one.
    res_bin_bed <- parLapply(cl, c(1:length(bed_ext_gr)), function(y) 
    make_scaled_coords(bed_ext_gr[y], FinalBins = BedRFinalLength,
                Extention = FALSE))
}
stopCluster(cl)
bed_bin_gr=unlist(as(res_bin_bed, "GRangesList"))
bed_bin_fname=paste0(Output_prefix, ".binned.bed")
export.bed(bed_bin_gr,con=bed_bin_fname)

cat(paste("#### Bed bined end", date(), "####", "", sep="\n"))

# Getting coverage from bigwig files on binned coordinates

Bins=BedRFinalLength
DScov=list()
DSdata=list()
cl=makeCluster(Threads)
clusterEvalQ(cl, {library(GenomicRanges); library(megadepth)})
clusterExport(cl, c("bw_fnames", "bw_names", "bed_bin_fname", "bed_bin_gr", "bed_ext_gr", "Bins"))
res_bin_density <- parLapply(cl, bw_fnames , function(x) {  
tt_cov=get_coverage(paste0(x), op="mean", annotation=bed_bin_fname)# Megadepth function to get mean coverage per interval
data=c();names=c()
for(i in 1:length(bed_ext_gr)){
    # retrieving the interval origin, by interval name if possible
    ovl= which(bed_bin_gr$name== bed_ext_gr[i]$name)
    # If there are several intervals with the same name, they should be on the same seq_id
    if (length(ovl) > Bins) {ovl = ovl[as.character(seqnames(bed_bin_gr)[ovl]) == as.character(seqnames(bed_ext_gr[i]))]}
    if(length(ovl)>0){
    score=as.numeric(tt_cov[ovl]$score)
    # If the strand in -, then we reverse the score order
    if(as.vector(strand(bed_ext_gr[i]))=='-'){score=rev(score)}
    data=rbind(data, score)
    names=c(names, paste(seqnames(bed_ext_gr[i]), start(bed_ext_gr[i]), end(bed_ext_gr[i]), bed_ext_gr[i]$name, bed_ext_gr[i]$score, strand(bed_ext_gr[i]), sep="_"))
    }
}
rownames(data)=names
return(data)
})
stopCluster(cl)
cat(paste("#### end", date(), "####", "", sep="\n"))
names(res_bin_density)=bw_names
bedData=res_bin_density
save(bedData, file=paste0(Output_prefix,".binned_density.R"))