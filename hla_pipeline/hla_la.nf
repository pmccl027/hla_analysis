
process HLA_typing {

	input:
	tuple file (bam), file(bam_index) from Channel.fromPath("extracted/*4889*.bam").map{ bam -> [ bam, bam + (bam.getExtension() == "bam" ? ".bai" : ".crai") ] }
	
	output:
	file "*.zip" into output
    
	publishDir "output/", pattern: "*.zip", mode: "copy"

	"""
	
	${params.HLA_path}/HLA-LA.pl --BAM ${bam} --graph ${params.graph} --sample sample1 --workingDir /home/pmccl/scratch/HLA-LA
     
	zip -r ${bam.simpleName}.zip ${params.out_path}sample1

	"""
}
