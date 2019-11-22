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

	label 'assembly'
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
	label 'assembly'

	publishDir "${OUTDIR}/fastp", mode: 'copy'

	input:
	set libraryID,file(forward),file(reverse) from inputFastp

	output:
	set libraryID,file(forward_trimmed),file(reverse_trimmed) into trimmed_reads
	set file(json),file(html) into qc_reports

	script:
	forward_trimmed = forward.getSimpleName() + "_trimmed.fastq.gz"
	reverse_trimmed = reverse.getSimpleName() + "_trimmed.fastq.gz"
	json =  libraryID + ".fastp.json"
	html =  libraryID + ".fastp.html"
	
	"""
		fastp --in1 $forward --in2 $reverse --out1 $forward_trimmed --out2 $reverse_trimmed --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
	"""
}

process runUnicycler {
	label 'assembly'

	publishDir "${OUTDIR}/${libraryID}/", mode: 'copy'

	input:
	set libraryID, file(fq1), file(fq2) from trimmed_reads

	output:
	set libraryID, file("${libraryID}_assembly.fasta") into inputDfast,quast_ch
	set libraryID, file("${libraryIDd}_assembly.gfa") into bandage_ch
	file("${libaryID}_assembly.gfa")
	file("${libraryID}_assembly.png")
	file("${library}_unicycler.log")

	// Stolen from Alex Pelzer, nf-core/bacass
	script:
	"""

	    unicycler -1 $fq1 -2 $fq2 --threads ${task.cpus} ${params.unicycler_args} --keep 0 -o .
	    mv unicycler.log ${libraryID}_unicycler.log
	    # rename so that quast can use the name 
	    mv assembly.gfa ${libraryID}_assembly.gfa
	    mv assembly.fasta ${libraryID}_assembly.fasta
	    Bandage image ${libraryID}_assembly.gfa ${libraryID}_assembly.png
	"""
}

process runQuast {
	
	label 'assembly'
	
	publishDir "${params.outdir}/${libraryIDid}/", mode: 'copy'
  
	input:
	set libraryID, fasta from quast_ch
  
	output:
	// multiqc only detects a file called report.tsv. to avoid
	// name clash with other samples we need a directory named by sample
	file("${libraryID}_assembly_QC/") into quast_logs_ch

	script:
	"""
		quast.py -t ${task.cpus} -o ${libraryID}_assembly_QC ${fasta}
	"""
}
process runDfast_core {

	label 'dfast'

        publishDir "${OUTDIR}/${libraryID}/annotation", mode: 'copy'

        input:
        set libraryID,file(assembly_fa) from inputDfast

        output:
        file("dfast/*") into DfastAnnotation
        file(annotation_gff) into DfastGFF
        set val(libraryID),file(annotation_fsa) into DfastFSA
	
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

	label 'assembly'

	publishDir "${OUTDIR}/MultiQC", mode: 'copy'

	input:
	file ('*')  from qc_reports.collect()
	file quast_logs from quast_logs_ch.collect().ifEmpty([])

	output:
	file "*library_multiqc.html" into multiqc_library_report

	script:
    
	"""
		multiqc -n library_multiqc -f . 2>&1
	"""	

}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}

