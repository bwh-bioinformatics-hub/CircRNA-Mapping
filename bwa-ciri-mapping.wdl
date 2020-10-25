## This is a script to align non-linear single-end reads with ciri
## Specifically for Nanopore data
## Inputs:
## 		sample 	- fastq read file 
##		gtf		- genome annotation file
##		fa		- genome index file
##		out 	- outfile prefix

workflow circRNA_alignment
{
	# input files
	File sample
    File gtf
    File fa
    String out

    # aligning to reference genome with ciri-full
    call alignment
    {
    	input:
          	sample=sample,
            gtf=gtf,
            fa=fa,
            out=out
    }
    
    call mapping
    {
    	input:
          	sample=sample,
            gtf=gtf,
            fa=fa,
            out=out
    }
    
    output
    {
    	Array[File] outputFiles = glob("ciri_output/*")
    }
    
}

task alignment
{
	File sample
    File gtf
    File fa
    String out
    
    command
    {
        mkdir ciri_output
		bwa index ${fa}
        
        bwa mem -T 19 -t 4 ${fa} ${sample} 1> ciri_output/${out}.sam 2> ${out}_bwa.log

	}
        
    runtime 
    {
        docker: "tveta/runciri:v1"
        memory: "60G"
        cpu: "4"
        disk: "local-disk 2000 HDD"
  	}
}

task mapping
{
	File sample
    File gtf
    File fa
    String out
    
    command
    {
        CIRI2.pl \
            -I ciri_output/${out}.sam \
            -O ciri_output/${out}.ciri \
            -F ${fa} \
            -A ${gtf} \
            -T 4
	}

    runtime 
    {
        docker: "tveta/runciri:v1"
        memory: "60G"
        cpu: "4"
        disk: "local-disk 2000 HDD"
  	}
}
