// Pipeline variables

/** 
=====================================================
IKMB Bacteria genome assembly and annotation pipeline
=====================================================

This pipelines performs adapter trimming, assembly and automated annotation of
bacteria-sized genomes from paired-end Illumina reads.

Written by: Marc Hoeppner, m.hoeppner@ikmb.uni-kiel.de

**/

// Configurable settings
CENTRE = params.centre
OUTDIR = params.outdir 

inputFile=file(params.samples)

// Header log info 

log.info "=========================================" 
log.info "IKMB pipeline version v${VERSION}" 
log.info "Nextflow Version: 	$workflow.nextflow.version" 
log.info "Command Line: 	$workflow.commandLine" 
log.info "Running Resfinder:	${params.resfinder}"
log.info "Resfinder DB:		${RESFINDER_DB}"
log.info "=========================================" 


// Starting the workflow

Channel.fromPath(inputFile)
	.splitCsv(sep: ';', skip:1 )
	.map { sampleID, libraryID, R1, R2 -> [ libraryID, sampleID, file(R1), file(R2) ] }
	.groupTuple(by: [0,1])
       	.set { inputMerge }

process Merge {

	tag "${libraryID}"
        // publishDir("${OUTDIR}/Data/${libraryID}")

        input:
	set libraryID,sampleID,forward_reads,reverse_reads from inputMerge

        output:
        set sampleID,libraryID,file(left_merged),file(right_merged) into inputTrimgalore

        script:
        left_merged = libraryID + "_R1.fastq.gz"
        right_merged = libraryID + "_R2.fastq.gz"

        """
                zcat ${forward_reads.join(" ")} | gzip > $left_merged
		zcat ${reverse_reads.join(" ")} | gzip > $right_merged
        """
}

process runFastp {
	
	publishDir "${OUTDIR}/fastp", mode: 'copy'

	input:
	set sampleID,libraryID,file(forward),file(reverse) from inputFastp

	output:
	set sampleID,libraryID,file(foward_trimmed),file(reverse_trimmed) into trimmed_reads,inputBwa
	set file(json),file(html) into qc_reports

	script:
	forward_trimmed = forward.getSimpleName() + "_trimmed.fastq.gz"
	reverse_trimmed = reverse.getSimpleName() + "_trimmed.fastq.gz"
	json = sampleID + "-" + libraryID + "fastp.json"
	html = sampleID + "-" + libraryID + "fastp.html"
	
	"""
		fastp --in1 $forward --in2 $reverse --out1 $forward_trimmed --out2 $reverse_trimmed --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
	"""
}

process runShovill {

  publishDir "${OUTDIR}/${sampleID}/assembly", mode: 'copy'
  
  input:
  set sampleID,libraryID,file(fw),file(rev) from trimmed_reads

  output:
  set sampleID,file(assembly_fa) into inputDfast,inputAssemblyMetrics

  script:
  assembly_fa = "shovill/contigs.fa"

  """
	shovill -R1 $fw -R2 $rev --cpus ${task.cpus} --ram ${task.memory.toGiga()} --outdir shovill
  """

}

process runDfast_core {

        publishDir "${OUTDIR}/${sampleID}/annotation", mode: 'copy'

        input:
        set sampleID,file(assembly_fa) from inputDfast

        output:
        file("dfast/*") into DfastAnnotation
        file(annotation_gff) into DfastGFF
        set val(sampleID),file(annotation_fsa) into DfastFSA
	
	script:
        annotation_gff = "dfast/genome.gff"
        annotation_gbk = "dfast/genome.gbk"
        annotation_stats = "dfast/statistics.txt"
        annotation_fsa  = "dfast/genome.fna"

	"""
		dfast --genome $assembly_fa --out dfast --minimum_length 200  --locus_tag_prefix ${sampleID} --cpu ${task.cpus} --center_name ${CENTRE}
	"""
}

process runMultiQCLibrary {

    tag "ALL"
    publishDir "${OUTDIR}/MultiQC", mode: 'copy'

    input:
    file ('*')  from qc_reports.collect()

    output:
    file "*library_multiqc.html" into multiqc_library_report

    script:
    
    """
      cp $baseDir/config/multiqc.yaml multiqc_config.yaml
      multiqc -n library_multiqc -f . 2>&1
    """	

}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}

//#############################################################################################################
//#############################################################################################################
//
// FUNCTIONS
//
//#############################################################################################################
//#############################################################################################################


// ------------------------------------------------------------------------------------------------------------
//
// Read input file and save it into list of lists
//
// ------------------------------------------------------------------------------------------------------------
def logParams(p, n) {
  File file = new File(n)
  file.write "Parameter:\tValue\n"

  for(s in p) {
     file << "${s.key}:\t${s.value}\n"
  }
}

