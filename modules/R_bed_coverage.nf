process R_BED_COVERAGE {
    label 'multiCpu'
    
    //tag "$BedName:$BedExtension-$BedExtLengthLeft:$BedExtLengthRight"
    
    //publishDir "${params.outdir}/${params.name}/Coverage", mode: 'copy', //params.publish_dir_mode,
    
    input:
    path(rfunction)
    path(rExec)
    tuple val(BedName), path(BedFile),val(BedExtension),val(BedCustomScaling), val(BedExtLengthLeft), val(BedExtLengthRight), val(BedRFinalLength),  val(BedExtValLeft), val(BedExtValRight)
    val(BwName)
    path(BwFile)

    tag_name = $BedName+"_Ext"+$BedExtension+"Scal"+$BedCustomScaling+"FL"+$BedRFinalLength+"L"+$BedExtLengthLeft+"R"+$BedExtLengthRight+
                "Vl"+$BedExtValLeft+'Vr'+$BedExtValRight
    tag $tag_name
/////// Need to define the output names.
    script:
    """
    echo "
        R.Version()
        RfunctionFile=\\"$rfunction\\"
        BedName=\\"${BedName}\\"
        BedFile=\\"${BedFile}\\"
        if('${BedExtension}'=='false'){BedExtension=FALSE}
        if('${BedExtension}'=='true'){BedExtension=TRUE}
        if('${BedCustomScaling}'=='false'){BedCustomScaling=FALSE}
        if('${BedCustomScaling}'=='true'){BedCustomScaling=TRUE}
        BedExtLengthLeft=${BedExtLengthLeft}
        BedExtLengthRight=${BedExtLengthRight}
        BedRFinalLength=${BedRFinalLength}
        BedExtValLeft=${BedExtValLeft}
        BedExtValRight=${BedExtValRight}
        Output_prefix=$tag_name

        bw_fnames=c('${BwFile.join('\',\'')}')
        bw_names=c('${BwName.join('\',\'')}')
        Threads=20
        source(RfunctionFile)
        source(\\"$rExec\\")
    "> r_GetCoverage_$BedName
    Rscript  r_GetCoverage_$BedName
    """
}