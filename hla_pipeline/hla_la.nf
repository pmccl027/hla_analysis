
process HLA_typing {

	input:
	tuple file (bam), file(bam_index) from Channel.fromPath("extracted/*4889*.bam").map{ bam -> [ bam, bam + (bam.getExtension() == "bam" ? ".bai" : ".crai") ] }
	
	output:
	file "*.zip" into output
	file("R1_bestquess_G.txt") into bestguess
    
	publishDir "output/", pattern: "*.zip", mode: "copy"

	"""
	
	${params.HLA_path}/HLA-LA.pl --BAM ${bam} --graph ${params.graph} --sample sample --workingDir . > ${bam.getSimpleName()}.log
	
	cp -r sample/hla ${bam.getSimpleName()}
        cp ${bam.getSimpleName()}.log ${bam.getSimpleName()}
	
	head -n1 sample/hla/R1_bestquess_G.txt > header.txt
	sed -i s/^/Sample\t/ header.txt
	
	tail -n+2 sample/hla/R1_bestquess_G.txt > data.txt
	sed -i s/^/${bam.getSimpleName()}/ data.txt
	
	cat header.txt data.txt > R1_bestquess_G.txt
     
	zip -r ${bam.getSimpleName()}.zip ${bam.getSimpleName()}

	"""
}

process merge {
	input:
	file(bestguess_files) from bestguess.collect()
}
