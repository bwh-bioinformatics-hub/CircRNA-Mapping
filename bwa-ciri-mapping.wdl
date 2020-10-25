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
          	sam=alignment.sam,
            gtf=gtf,
            fa=fa,
            out=out
    }
}

task alignment
{
	File sample
    File gtf
    File fa
    
    command
    {
        bwa index ${fa}
        
        bwa mem -T 19 -t 4 ${fa} ${sample} 1> aln.sam 2> aln_bwa.log

	}
    
    output
    {
        File sam = "aln.sam"
        File log = "aln_bwa.log"
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
	File sam
    File gtf
    File fa
    
    command
    {
        CIRI2.pl \
            -I ${sam} \
            -O aln.ciri \
            -F ${fa} \
            -A ${gtf} \
            -T 4
	}

    output
    {
        File out = "aln.ciri"
    }

    runtime 
    {
        docker: "tveta/runciri:v1"
        memory: "60G"
        cpu: "4"
        disk: "local-disk 2000 HDD"
  	}
}
