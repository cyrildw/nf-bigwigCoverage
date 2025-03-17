process R_BED_COVERAGE {
    tag "$BedName:$BedExtension-$BedExtLengthLeft:$BedExtLengthRight"
    label 'multiCpu'
    publishDir "${params.outdir}/${params.name}/Coverage", mode: 'copy', //params.publish_dir_mode,
    
    input:
    tuple BedName, file(BedFile), BedExtLengthLeft, BedExtLengthRight, BedRFinalLength, BedExtension, BedExtValLeft, BedExtValRight 
    tuple val(BwName),
    tuple path(BwFile)
// path genome

    output:
    """
    echo "Processing BED: ${BedName} (${bed_file})
    ${BwName[@]}"
    """
}