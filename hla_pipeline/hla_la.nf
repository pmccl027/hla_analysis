
process HLA_typing {

	input:
	tuple file (bam), file(bam_index) from Channel.fromPath("extracted/*.bam").map{ bam -> [ bam, bam + (bam.getExtension() == "bam" ? ".bai" : ".crai") ] }
	
	output:
	file "*.zip" into output
    
	publishDir "output/", pattern: "*.zip", mode: "copy"

	"""
	
	${params.HLA_path}/HLA-LA.pl --BAM ${bam} --graph ${params.graph} --sample sample --workingDir . > ${bam.getSimpleName()}.log
	
	cp -r sample/hla ${bam.getSimpleName()}
        cp ${bam.getSimpleName()}.log ${bam.getSimpleName()}
     
	zip -r ${bam.getSimpleName()}.zip ${bam.getSimpleName()}

	"""
}
