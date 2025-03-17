nextflow.enable.dsl = 2
i=0
Channel
      .fromPath(params.bigwig_design)
      .splitCsv(header:true, sep:';')
      .map { row -> [ row.BwName, i++, file("$params.input_dir/$row.BwFile", checkIfExists: true)]}
      .set { ch_BW_design_csv}
ch_BW_design_csv.view()

ch_BW_design_csv.map { it -> [ it[0], it[3]]}
    .multiMap { it ->   labels: it[0]
                        files: it[1]}
    .set{ch_flat_BW}
ch_flat_BW.labels.view()
//ch_flat_BW.labels.view()
//ch_design_reads_csv.view()


//
// MODULE LOAD
//
//include { R_BED_COVERAGE as R_BED_COVERAGE     } from './modules/R_bed_coverage'


Channel
    .fromPath(params.bed_design)
    .splitCsv(header:true, sep:";")
    .map { row -> [row.BedName,
		file("$params.input_dir/$row.BedFile", checkIfExists: true),
        //file("$params.input_dir/$row.BedGroupFile"),
		//row.BedDTlength, //Used by Deeptools for -m option in plot heatmap
        //row.BedReferencePoint, //Used by Deeptools in opposition to "scale-region"
		row.BedExtLengthLeft,  //Used by Deeptools and R to extend the bed coordinates upstream
		row.BedExtLengthRight, //Used by Deeptools and R to extend the bed coordinates downstream
		row.BedRFinalLength,    //Used by R to set the final length of the vector
		row.BedExtension,      // Should the bed coordinates be extended
		row.BedExtValLeft,     // Used by R, how much of the FinalLength should the upstream extension represent
        row.BedExtValRight]     // Used by R, how much of the FinalLength should the upstream extension represent
        }
    .set { ch_BED_design_csv }

ch_BED_bw_combined=ch_BED_design_csv.combine(ch_flat_BW)
ch_BED_bw_combined.view()
