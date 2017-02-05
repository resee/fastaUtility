_______________________________
fastaUtility.py 

Rishi Aryal, 2017

______________________________




USASE: fastaUtility.py [options: -h,-t,-s,-a,-n,-e,-r,-l,-w] -f fastafile [-i idfile]
   
   -h: Usage help 
   
   OPTION (each of these options produce separate output file with the results)
   
   -t: tab delimited format
   
   -s: sub sequences (need to provide an idfile with title, and start and end posiotions,
                      the title should match the first word of the sequecne title)
                      
   -a: amino acids, finds the longest orf (from all 6 reading frame) and 
                     yields amino acid sequence of the orf
                     
   -n: nucleotide change- DNA to RNA and vice versa
   
   -e: extracted fasta sequence (need to proved an idfile with titles to extract
                                 the title should match the first word of the sequecne title)
                                 
   -r: reverse complement on RNA or DNA sequences
   
   -l: sequence length on each fasta title (adds sequence length info on the fasta title)
   
   -w: sequence with swapped-case for longest ORF (from all 6 reading frame) 
        eg.  ctgtatgtacatatgaactat produces ctgtATGTACATATGAactat
   
   file options:
   -f: Fastafile (need for all above options)
   -i: IDfile (need only with the options -s and -e options above)
   
   
