process R_BED_COVERAGE {
    label 'multiCpu'
    //tag "$BedName:$BedExtension-$BedExtLengthLeft:$BedExtLengthRight"
    
    //publishDir "${params.outdir}/${params.name}/Coverage", mode: 'copy', //params.publish_dir_mode,
    
    input:
    path(rfunction)
    path(rExec)
    tuple val(BedName), path(BedFile), val(BedExtLengthLeft), val(BedExtLengthRight), val(BedRFinalLength), val(BedExtension), val(BedExtValLeft), val(BedExtValRight)
    val(BwName)
    path(BwFile)
// path genome
///// NEED TO SETUP A R running script with "source" this header.
    script:
    """
    echo "R --no-save --no-restore --slave <<RSCRIPT
        R.Version()
        RfunctionFile=\"$rfunction\"
        BedName=\"${BedName}\"
        BedFile=\"${BedFile}\"
        if(${BedExtension}=='false'){BedExtension=FALSE}
        if(${BedExtension}=='true'){BedExtension=TRUE}
        BedExtLengthLeft=${BedExtLengthLeft}
        BedExtLengthRight=${BedExtLengthRight}
        BedRFinalLength=${BedRFinalLength}
        BedExtValLeft=${BedExtValLeft}
        BedExtValRight=${BedExtValRight}

        bw_fpath=\"$params.input_dir\"
        bw_fnames=c('${BwName.join('\',\'')}')
        bw_names=c('${BwFile.join('\',\'')}')
        Threads=20
        source(\"$rExec\")
    " > r_GetCoverage_$BedName
    bash  r_GetCoverage_$BedName
    """
}