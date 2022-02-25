
process Extract_HLA {

	input:
	tuple file (bam), file(bam_index) from Channel.fromPath("${params.path}*99.*.cram").map{ bam -> [ bam, bam + (bam.getExtension() == "bam" ? ".bai" : ".crai") ] }
	
	output:
	file "*extracted.bam" into extracted
	file "*coverage.txt" into coverage

	publishDir "extracted/", pattern: "*extracted.bam", mode: "copy"
	publishDir "coverage/", pattern: "*coverage.txt", mode: "copy"

	"""
	
	# unmapped reads w/ mapped mate
	samtools view -T ${params.reference} -f 4 -F 264 ${bam} -b -o all.1.bam

	# mapped reads w/ unmapped mate
	samtools view -T ${params.reference} -f 8 -F 260 ${bam} -b -o all.2.bam

	# both unmapped
	samtools view -T ${params.reference} -f 12 -F 256 ${bam} -b -o all.3.bam

	# mapped in HLA
	samtools view -T ${params.reference} -F 12 ${bam} chr6:25000000-35000000 -b -o hla.mapped_only.bam

	# all merged
	samtools merge ${bam.simpleName}_extracted.bam all.1.bam all.2.bam all.3.bam hla.mapped_only.bam
	# index new extracted file
	samtools index ${bam.simpleName}_extracted.bam

	# compute coverage
	samtools coverage -q 10 -Q 10 -r chr6:25000000-35000000 ${bam.simpleName}_extracted.bam | cut -f 4- > ${bam.simpleName}_coverage.txt

	"""
}
