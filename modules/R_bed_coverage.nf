process R_BED_COVERAGE {
    label 'multiCpu'
    //tag "$BedName:$BedExtension-$BedExtLengthLeft:$BedExtLengthRight"
    
    //publishDir "${params.outdir}/${params.name}/Coverage", mode: 'copy', //params.publish_dir_mode,
    
    input:
    tuple val(BedName), path(BedFile), val(BedExtLengthLeft), val(BedExtLengthRight), val(BedRFinalLength), val(BedExtension), val(BedExtValLeft), val(BedExtValRight)
    val(BwName)
    val(BwFile)
// path genome

    script:
    """
    echo "R --no-save --no-restore --slave <<RSCRIPT
        R.Version()
        source($params.rfunction)
        BedName="${BedName}"
        BedFile="${BedFile}"
        if(${BedExtension}=='false'){BedExtension=FALSE}
        if(${BedExtension}=='true'){BedExtension=TRUE}
        BedExtLengthLeft=${BedExtLengthLeft}
        BedExtLengthRight=${BedExtLengthRight}
        BedRFinalLength=${BedRFinalLength}
        BedExtValLeft=${BedExtValLeft}
        BedExtValRight=${BedExtValRight}
        Threads=20

        #Importing bed and filtering
        bed_gr=import(BedFile)
        bed_ext_gr=bed_gr

        #Creating GR bed object with extended coordinates
        start(bed_ext_gr)=start(bed_gr)-1-BedExtLengthLeft
        if(sum(start(bed_ext_gr)<1)>0){paste0(\"Attention, at least \", sum(start(bed_ext_gr)<1), ' invervals are out of chromosomal range (start)')}
        start(bed_ext_gr)[start(bed_ext_gr)<1] = 1
        end(bed_ext_gr)=end(bed_gr)+BedExtLengthRight

        cat(paste(\"### Bed bined start\", date(), \"####\", \"\", sep=\"\\n\"))

        #Creating GR bined bed object with extended coordinates according to configuration (Takes a bit of time)
        cl=makeCluster(Threads)
        clusterEvalQ(cl, source($params.rfunction))
        clusterExport(cl, c(\"bed_gr\",\"BedRFinalLength\",\"BedExtLengthLeft\", \"BedExtension\",
                            \"BedExtLengthRight\", \"BedExtValLeft\",\"BedExtValRight\"))
        res_bin_bed <- parLapply(cl, c(1:length(bed_gr)), function(y) 
        make_scaled_coords(bed_gr[y], FinalBins = BedRFinalLength,
                    Extention = BedExtension, Ext_length = c(BedExtLengthLeft, BedExtLengthRight),
                    Ext_value = c(BedExtValLeft,BedExtValRight)))
        stopCluster(cl)
        bed_bin_gr=unlist(as(res_bin_bed, \"GRangesList\"))
        bed_bin_fname=paste0(strsplit(basename(BedFile), split=\".bed\")[[1]][1], 
                            '_Ex', BedExtension,'L',BedExtLengthLeft, 'R', BedExtLengthRight,
                            'Vl',BedExtValLeft,'Vr',BedExtValRight, \".bed\")
        export.bed(bed_bin_gr,con=bed_bin_fname)

        cat(paste(\"#### Bed bined end\", date(), \"####\", \"\", sep=\"\\n\"))

        # Getting coverage from bigwig files on binned coordinates
        bw_fpath=$params.input_dir
        bw_fnames=c('${BwName.join('\',\'')}')
        bw_names=c('${BwFile.join('\',\'')}')
    
        Bins=BedRFinalLength
        DScov=list()
        DSdata=list()
        cl=makeCluster(Threads)
        clusterEvalQ(cl, {library(GenomicRanges); library(megadepth)})
        clusterExport(cl, c(\"bw_fpath\", \"bw_fnames\", \"bw_names\", \"bed_bin_fname\", \"bed_bin_gr\", \"bed_ext_gr\", \"Bins\"))
        res4 <- parLapply(cl, bw_fnames , function(x) {  
        tt_cov=get_coverage(paste0(bw_fpath,x), op=\"mean\", annotation=bed_bin_fname)# Roughly 4-6min
        data=c();names=c()
        for(i in 1:length(bed_ext_gr)){
            ovl= which(bed_bin_gr\$name== bed_ext_gr[i]\$name)
            if (length(ovl) > Bins) {ovl = ovl[as.character(seqnames(bed_bin_gr)[ovl]) == as.character(seqnames(bed_ext_gr[i]))]}
            if(length(ovl)>0){
            score=as.numeric(tt_cov[ovl]\$score)
            if(as.vector(strand(bed_ext_gr[i]))=='-'){score=rev(score)}
            data=rbind(data, score)
            names=c(names, paste(seqnames(bed_ext_gr[i]), start(bed_ext_gr[i]), end(bed_ext_gr[i]), bed_ext_gr[i]\$name, bed_ext_gr[i]\$score, strand(bed_ext_gr[i]), sep=\"_\"))
            }
        }
        rownames(data)=names
        return(data)
        })
        stopCluster(cl)
        cat(paste(\"#### end\", date(), \"####\", \"\", sep=\"\\n\"))
        names(res4)=bw_names
        bedData=res4
        save(bedData, file=paste0(BedName,\".binned.R\")
        "
    """
}