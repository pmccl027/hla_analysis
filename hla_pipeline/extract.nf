
process Extract_HLA {

	input:
	#file cram from Channel.fromPath("/lustre03/project/rrg-vmooser/dtaliun/BQC19_Globus/Release5_DEC_2021/BQC10000-BQC10099/cram/*.cram")
	file cram from Channel.fromPath("/lustre03/project/rrg-vmooser/dtaliun/BQC19_Globus/Release5_DEC_2021/BQC10000-BQC10099/cram/*99.cram")
	
	output:
	file "*extracted.bam" into extracted

	publishDir "extracted/", pattern: "*extracted.bam", mode: "copy"

	"""
	reference=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa
	
	# unmapped reads w/ mapped mate
	samtools view -T ${reference} -f 4 -F 264 ${cram} -b -o all.1.bam

	# mapped reads w/ unmapped mate
	samtools view -T ${reference} -f 8 -F 260 ${cram} -b -o all.2.bam

	# both unmapped
	samtools view -T ${reference} -f 12 -F 256 ${cram} -b -o all.3.bam

	# mapped in HLA
	samtools view -T ${reference} -F 12 ${cram} chr6:25000000-35000000 -b -o hla.mapped_only.bam

	# all merged
	samtools merge hla_extracted.bam all.1.bam all.2.bam all.3.bam hla.mapped_only.bam

	"""
}
