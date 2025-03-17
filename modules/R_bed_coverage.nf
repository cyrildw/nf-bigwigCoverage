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
        bedName='${BedName}'
        libNames=c('${BwName.join('\',\'')}')
        libFiles=c('${BwFile.join('\',\'')}')
        sortby=order(libNames);
        libNames=libNames[sortby]; libFiles=libFiles[sortby]
        resTable=read.table(libFiles[1], head=FALSE, stringsAsFactor=FALSE)[,c(1, 2, 3, 6)]
        colnames(resTable)=c(\'chr\', \'start\', \'end\', \'strand\')
        for( i in 1:length(libNames)){
            resTable=cbind(resTable, read.table(libFiles[i], head=FALSE, stringsAsFactor=FALSE)[,c(5)])
            colnames(resTable)[dim(resTable)[2]]=libNames[i]}
        write.table(resTable, file='${BedName}.avg_TagDensity.tsv', quote=FALSE, row.names=FALSE, col.names=TRUE, sep=\'\t\')
        "
    """
}