process R_BED_COVERAGE {
    label 'multiCpu'
    
    tag "${BedName}_Ext${BedExtension}Scal${BedCustomScaling}FL${BedRFinalLength}L${BedExtLengthLeft}R${BedExtLengthRight}Vl${BedExtValLeft}Vr${BedExtValRight}"
    
    publishDir "${params.outdir}/${params.name}/${tag_name}/", mode: 'copy' //params.publish_dir_mode,
    
    input:
    path(rfunction)
    path(rExec)
    tuple val(BedName), path(BedFile),val(BedExtension),val(BedCustomScaling), val(BedExtLengthLeft), val(BedExtLengthRight), val(BedRFinalLength),  val(BedExtValLeft), val(BedExtValRight)
    val(BwName)
    path(BwFile)

    output:
    tuple path($tag_name".ext.bed"), path($tag_name'.avg_density.tsv'), path($tag_name".binned_density.R")
    ////// Need to define the output names.
    script:
    def tag_name="${BedName}_Ext${BedExtension}Scal${BedCustomScaling}FL${BedRFinalLength}L${BedExtLengthLeft}R${BedExtLengthRight}Vl${BedExtValLeft}Vr${BedExtValRight}"
    """
    echo "
        R.Version()
        RfunctionFile=\\"$rfunction\\"
        BedName=\\"${BedName}\\"
        BedFile=\\"${BedFile}\\"
        if('${BedExtension}'=='true'){BedExtension=TRUE}else{BedExtension=FALSE}
        if('${BedCustomScaling}'=='true'){BedCustomScaling=TRUE}else{BedCustomScaling=FALSE}
        BedExtLengthLeft=${BedExtLengthLeft}
        BedExtLengthRight=${BedExtLengthRight}
        BedRFinalLength=${BedRFinalLength}
        BedExtValLeft=${BedExtValLeft}
        BedExtValRight=${BedExtValRight}
        Output_prefix=\\"$tag_name\\"

        bw_fnames=c('${BwFile.join('\',\'')}')
        bw_names=c('${BwName.join('\',\'')}')
        Threads=20
        source(RfunctionFile)
        source(\\"$rExec\\")
    "> r_GetCoverage_$BedName
    Rscript  r_GetCoverage_$BedName
    """
}