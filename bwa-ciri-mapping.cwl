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
    call mapping
    {
    	input:
          	sample=sample,
            gtf=gtf,
            fa=fa,
            out=out
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
        mkdir ciri_output
		bwa index ${fa}
        
        bwa mem -T 19 -t 4 ${fa} ${sample} 1> ciri_output/${out}.sam 2> ${out}_bwa.log
        
        CIRI2.pl \
            -I ciri_output/${out}.sam \
            -O ciri_output/${out}.ciri \
            -F ${fa} \
            -A ${gtf} \
            -T 4
        
        # CIRI_AS_v1.2.pl \
        # 	-S ciri_output/${out}.sam \
        #     -C ciri_output/${out}.ciri \
        #     -F ${fa} \
        #     -A ${gtf} \
        #     -O ciri_output/${out} \
        #     -D yes
        
        # java -jar /usr/sbin/CIRI-vis.jar \
        # 	-i ciri_output/${out}_jav.list \
        #     -l ciri_output/${out}_library_length.list \
        #     -r ${fa} \
        #     -min 1 \
        #     -d ciri_output \
        #     -o ${out}
        
	}
    
    output
    {
    	Array[File] outputFiles = glob("ciri_output/*")
    }
    
    runtime 
    {
        docker: "tveta/runciri:v1"
        memory: "60G"
        cpu: "4"
        disk: "local-disk 2000 HDD"
  	}
}