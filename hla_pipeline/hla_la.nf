
process HLA_typing {

	input:
	tuple file (bam), file(bam_index) from Channel.fromPath("extracted/*.bam").map{ bam -> [ bam, bam + (bam.getExtension() == "bam" ? ".bai" : ".crai") ] }
	
	output:
	file "*.zip" into output
	file("R1_bestguess_G_*.txt") into bestguess
    
	publishDir "output/", pattern: "*.zip", mode: "copy"

	"""
	
	${params.HLA_path}/HLA-LA.pl --BAM ${bam} --graph ${params.graph} --sample sample --workingDir . > ${bam.getSimpleName()}.log
	
	cp -r sample/hla ${bam.getSimpleName()}
        cp ${bam.getSimpleName()}.log ${bam.getSimpleName()}
	
	head -n1 sample/hla/R1_bestguess_G.txt > header.txt
	sed -i 's/^/Sample\t/' header.txt
	
	tail -n+2 sample/hla/R1_bestguess_G.txt > data.txt
	sed -i 's/^/${bam.getSimpleName()}\t/' data.txt
	
	cat header.txt data.txt > R1_bestguess_G_${bam.getSimpleName()}.txt
     
	zip -r ${bam.getSimpleName()}.zip ${bam.getSimpleName()}

	"""
}

process merge {
	input:
	file(bestguess_files) from bestguess.collect()

	output:
	file "merged_bestguess.txt" into merged_bestguess

	publishDir "merged_bestguess/", pattern: "merged_bestguess.txt", mode: "copy"	

	"""
	cat ${bestguess_files} | grep -vw "^Sample" > no_header.txt
	cat ${bestguess_files} | head -n1 > only_header.txt
	cat only_header.txt no_header.txt >> merged_bestguess.txt
	"""
}
