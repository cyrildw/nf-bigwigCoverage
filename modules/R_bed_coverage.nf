process R_BED_COVERAGE {
    label 'multiCpu'
    //tag "$BedName:$BedExtension-$BedExtLengthLeft:$BedExtLengthRight"
    
    //publishDir "${params.outdir}/${params.name}/Coverage", mode: 'copy', //params.publish_dir_mode,
    
    input:
    tuple BedName, file(BedFile), BedExtLengthLeft, BedExtLengthRight, BedRFinalLength, BedExtension, BedExtValLeft, BedExtValRight 
    tuple val(BwName),
    tuple path(BwFile)
// path genome

    script:
    """
    echo "Processing BED: ${BedName} (${bed_file})"
    """
}