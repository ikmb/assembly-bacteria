// Pipeline variables

/** 
=====================================================
IKMB Bacteria genome assembly and annotation pipeline
=====================================================

This pipelines performs adapter trimming, assembly and automated annotation of
bacteria-sized genomes from paired-end Illumina reads.

Written by: Marc Hoeppner, m.hoeppner@ikmb.uni-kiel.de

**/

// Help message
helpMessage = """
===============================================================================
IKMB Bacteria assembly and annotation pipeline | version ${workflow.manifest.version}
===============================================================================
Usage: nextflow run ikmb/assembly-bacteria --reads '*_R{1,2}_001.fastq.gz' --email 'you@somehwere.com'

Required parameters:
--reads                        Regexp to define input files (must be paired-end Illumina reads)
Optional parameters:
--email                        Email address to send reports to (enclosed in '')
--skip_multiqc                 Don't attached MultiQC report to the email.
Output:
--outdir                       Local directory to which all output is written (default: results)
"""

params.help = false

// Configurable settings
CENTRE = params.centre
OUTDIR = params.outdir 


// Header log info 
log.info "=========================================" 
log.info "IKMB pipeline version v${workflow.manifest.version}" 
log.info "Nextflow Version: 	$workflow.nextflow.version" 
log.info "Command Line: 	$workflow.commandLine" 
log.info "=========================================" 


// Starting the workflow

Channel.fromFilePairs(params.reads, flat: true )
	.ifEmpty { exit 1, "Did not find any reads matching your argument --reads" }
	.set { reads }

inputMerge = reads.groupTuple()

process Merge {

        input:
	set libraryID,file(forward_reads),file(reverse_reads) from inputMerge

        output:
        set libraryID,file(left_merged),file(right_merged) into inputFastp

        script:
        left_merged = libraryID + "_R1.fastq.gz"
        right_merged = libraryID + "_R2.fastq.gz"

	if (forward_reads.size() > 1 && forward_reads.size() < 1000) {

	        """
        	        zcat ${forward_reads.join(" ")} | gzip > $left_merged
			zcat ${reverse_reads.join(" ")} | gzip > $right_merged
	        """
	} else {

		"""
			cp $forward_reads $left_merged
			cp $reverse_reads $right_merged
		"""
	}
}

process runFastp {
	
	publishDir "${OUTDIR}/fastp", mode: 'copy'

	input:
	set libraryID,file(forward),file(reverse) from inputFastp

	output:
	set libraryID,file(forward_trimmed),file(reverse_trimmed) into trimmed_reads
	set file(json),file(html) into qc_reports

	script:
	forward_trimmed = forward.getSimpleName() + "_trimmed.fastq.gz"
	reverse_trimmed = reverse.getSimpleName() + "_trimmed.fastq.gz"
	json =  libraryID + "fastp.json"
	html =  libraryID + "fastp.html"
	
	"""
		fastp --in1 $forward --in2 $reverse --out1 $forward_trimmed --out2 $reverse_trimmed --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
	"""
}

process runShovill {

  publishDir "${OUTDIR}/${libraryID}/assembly", mode: 'copy'
  
  input:
  set libraryID,file(fw),file(rev) from trimmed_reads

  output:
  set libraryID,file(assembly_fa) into inputDfast,inputAssemblyMetrics

  script:
  assembly_fa = "shovill/contigs.fa"

  """
	shovill -R1 $fw -R2 $rev --cpus ${task.cpus} --ram ${task.memory.toGiga()} --outdir shovill
  """

}

process runDfast_core {

        publishDir "${OUTDIR}/${libraryID}/annotation", mode: 'copy'

        input:
        set libraryID,file(assembly_fa) from inputDfast

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
		dfast --genome $assembly_fa --out dfast --minimum_length 200  --locus_tag_prefix ${libraryID} --cpu ${task.cpus} --center_name ${CENTRE}
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

