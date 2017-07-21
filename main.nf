// Pipeline variables

CENTRE = params.centre

PROKKA = file(params.prokka)

BWA = file(params.bwa)

PICARD = file(params.picard)

SPADES = file(params.spades)

OUTDIR = params.outdir 

SAMTOOLS=file(params.samtools)

TRIMGALORE=file(params.trimgalore)

inputFile=file(params.samples)

params.saveTrimmed = true


// Logging and reporting

logParams(params, "nextflow_parameters.txt")

VERSION = "2.0" 
// Header log info 

log.info "=========================================" 
log.info "IKMB pipeline version v${VERSION}" 
log.info "Nextflow Version: $workflow.nextflow.version" 
log.info "Command Line: $workflow.commandLine" 
log.info "=========================================" 


// Starting the workflow

Channel.from(inputFile)
       .splitCsv(sep: ';', header: true)
       .into {  inputTrimgalore }

process runTrimgalore {

   tag "${sampleID}|${libraryID}"
   publishDir "${OUTDIR}/trimgalore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
	    else params.saveTrimmed ? filename : null
        }

   input:
   set sampleID,libraryID,forward,reverse from inputTrimgalore

   output:
   set sampleID,libraryID,file("*val_1.fq.gz"),file("*val_2.fq.gz") into trimmed_reads, inputReadsBwa 
   file "*trimming_report.txt" into trimgalore_results   
   file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

   script:

   trim_option = '--paired'

   """
	$TRIMGALORE --gzip $trim_option --fastqc --length 75  $forward $reverse
   """	
}

process runSpades {

  tag "${sampleID}|#{libraryID}"
  publishDir "${OUTDIR}/${sampleID}/assembly", mode: 'copy'
  
  input:
  set sampleID,libraryID,file(fw),file(rev) from trimmed_reads

  output:
  set sampleID,file(assembly_fa) into inputProkka,inputAssemblyBwa,inputAssemblyMetrics

  script:
  assembly_fa = "spades/scaffolds.fasta"

  """
	$SPADES -1 $fw -2 $rev -t 16 -m 120 -o spades
  """

}

process runProkka {
	
	tag "${sampleID}"
	publishDir "${OUTDIR}/${sampleID}/annotation", mode: 'copy'

	input:
	set sampleID,file(assembly_fa) from inputProkka

	output:
	set sampleID,file(annotation_gff),file(annotation_gbk) into ProkkaAnnotation
	file(annotation_stats) into ProkkaStats

	script:
	annotation_gff = "prokka/" + sampleID + ".gff"
	annotation_gbk = "prokka/" + sampleID + ".gbk"
	annotation_stats = "prokka/" + sampleID + ".txt"

	"""
		$PROKKA --outdir prokka --prefix ${sampleID} --addgenes --locustag ${sampleID} --mincontiglen 200 --centre ${CENTRE} --strain ${sampleID} --cpus 16 $assembly_fa
	"""

}


mergeAssemblyAndReads_by_id = inputAssemblyBwa.combine(inputReadsBwa, by: 0)

process runBwa {

	tag "${sampleID}"
        publishDir "${OUTDIR}/${sampleID}/BAM", mode: 'copy'

        input:
        set sampleID,file(assembly),file(left),file(right) from mergeAssemblyAndReads_by_id

        output:
        set sampleID,file(bam) into outputBwa

	script:
	bam = sampleID + ".bam"

	"""
		$BWA index $assembly && $BWA mem -M -t 8 ${assembly} $left $right | $SAMTOOLS sort - > $bam
	"""

}

process runMarkDuplicates {

    tag "${sampleID}"
    publishDir "${OUTDIR}/${sampleID}/MarkDuplicates", mode: 'copy'

    input:
    set sampleID,file(bam) from outputBwa
    
    output:
    set sampleID,file(outfile_bam),file(outfile_bai) into outputMarkDuplicates,inputRunCollectMultipleMetrics
    file(outfile_metrics) into runMarkDuplicatesOutput_QC
    
    script:
    outfile_bam = sampleID + ".dedup.bam"
    outfile_bai = sampleID + ".dedup.bai"

    outfile_metrics = sampleID + "_duplicate_metrics.txt"	
	        
    """
    
	java -Xmx25G -XX:ParallelGCThreads=5 -Djava.io.tmpdir=tmp/ -jar ${PICARD} MarkDuplicates \
		INPUT=${bam} \
		OUTPUT=${outfile_bam} \
		METRICS_FILE=${outfile_metrics} \
		CREATE_INDEX=true \
		TMP_DIR=tmp
	"""  
}

combinedBamAndAssembly = inputRunCollectMultipleMetrics.combine(inputAssemblyMetrics, by: 0)

process runCollectMultipleMetrics {

    tag "${sampleID}"
    publishDir "${OUTDIR}/${sampleID}/Picard_Metrics", mode: 'copy'
 	    
    input:
    set sampleID, file(bam), file(bai),file(assembly) from combinedBamAndAssembly

    output:
    file("${sampleID}*") into CollectMultipleMetricsOutput mode flatten

    script:       

    """
        java -XX:ParallelGCThreads=1 -Xmx5g -Djava.io.tmpdir=tmp/ -jar $PICARD CreateSequenceDictionary \
		R=$assembly \
		O=scaffolds.dict \
	&& \
	java -XX:ParallelGCThreads=1 -Xmx5g -Djava.io.tmpdir=tmp/ -jar $PICARD CollectMultipleMetrics \
		PROGRAM=MeanQualityByCycle \
		PROGRAM=QualityScoreDistribution \
		PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics\
		PROGRAM=CollectGcBiasMetrics \
		PROGRAM=CollectBaseDistributionByCycle \
		INPUT=${bam} \
		REFERENCE_SEQUENCE=scaffolds.fasta \
		ASSUME_SORTED=true \
		QUIET=true \
		OUTPUT=${sampleID} \
		TMP_DIR=tmp
	"""
}	


process runMultiQCSample {

    tag "ALL"
    publishDir "${OUTDIR}/MultiQC", mode: 'copy'

    input:
  //  file (fastqc:'fastqc/*') from trimgalore_fastqc_reports.collect()
 //   file ('trimgalore/*') from trimgalore_results.collect()
    file (dedup_metrics) from runMarkDuplicatesOutput_QC.collect()
    file (prokka_stats) from ProkkaStats.collect()
    file('*') from CollectMultipleMetricsOutput.flatten().toList()



    output:
    file "*assembly_multiqc.html" into multiqc_report
    file "*multiqc_data"
    // val prefix into multiqc_prefix

    script:
    //prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    """
    multiqc -n assembly_multiqc -f . 2>&1
    """




}

process runMultiQCLibrary {

    tag "ALL"
    publishDir "${OUTDIR}/MultiQC", mode: 'copy'

    input:
    file (fastqc:'fastqc/*') from trimgalore_fastqc_reports.collect()
    file ('trimgalore/*') from trimgalore_results.collect()

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

