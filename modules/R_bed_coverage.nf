process R_BED_COVERAGE {
    label 'multiCpu'
    cpus 20
    memory '60 GB'
    queue 'run16'
    tag "${BedName}_Ext${BedExtension}Scal${BedCustomScaling}FL${BedRFinalLength}L${BedExtLengthLeft}R${BedExtLengthRight}Vl${BedExtValLeft}Vr${BedExtValRight}"
    
    publishDir "${params.outdir}/${params.name}/BigWigCoverage/", mode: 'copy' //params.publish_dir_mode,
    
    input:
    path(rfunction)
    path(rExec)
    tuple val(BedName), path(BedFile),val(BedExtension),val(BedCustomScaling), val(BedExtLengthLeft), val(BedExtLengthRight), val(BedRFinalLength),  val(BedExtValLeft), val(BedExtValRight)
    val(BwName)
    path(BwFile)

    output:
    tuple val(BedName), path("*.ext.bed"), path("*.avg_density.tsv"), path("*.binned_density.R")
    
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
        tmp_dir=\\"${workflow.workDir}\\"
        bw_fnames=c('${BwFile.join('\',\'')}')
        bw_names=c('${BwName.join('\',\'')}')
        Threads=${task.cpus}
        source(RfunctionFile)
        source(\\"$rExec\\")
    "> r_GetCoverage_$BedName
    Rscript  r_GetCoverage_$BedName
    """
}