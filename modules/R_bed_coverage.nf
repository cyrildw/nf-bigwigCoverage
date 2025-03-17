process R_BED_COVERAGE {
    label 'multiCpu'
    //tag "$BedName:$BedExtension-$BedExtLengthLeft:$BedExtLengthRight"
    
    //publishDir "${params.outdir}/${params.name}/Coverage", mode: 'copy', //params.publish_dir_mode,
    
    input:
    tuple val(BedName), file(BedFile), val(BedExtLengthLeft), val(BedExtLengthRight), val(BedRFinalLength), val(BedExtension), val(BedExtValLeft), val(BedExtValRight)
    tuple val(BwName),
    tuple path(BwFile)
// path genome

    script:
    """
    echo "Processing BED: ${BedName} (${bed_file})"
    """
}