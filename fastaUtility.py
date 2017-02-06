#fastaUtility.py
#Rishi Aryal, 2017 (reseearyal[at]gmail[dot]com)

'''
Usage: fastaUtility.py [options: -h,-t,-s,-a,-n,-e,-r,-l,-w] -f <fastafile> -i <idfile> 
   
   -h: Usage help 
   
   OPTIONs: produced files with 
   -t: tab delimited format
   -s: sub sequences (need to provide idfile with title start and end posiotions)
   -a: amino acids
   -n: nucleotide changed- DNA to RNA and vice versa
   -e: extracted fasta sequence (need to proved idfile with titles to extract)
   -r: reverse complement or RNA or DNA sequences
   -l: sequence length on each fasta title
   -w: sequence with swapped-case for longest ORF 
        eg.  ctgtatgtacatatgaactat produces ctgtATGTACATATGAactat
   
   file options:
   -f: Fasta file ( need for all above options)
   -i: ID file (optional; need to provide only with -s and -e options above)
   
   '''


#############################################

import re, sys, getopt, textwrap 

def tab_generator(fastafile):
   '''produces title and sequence for each fasta sequence for further processing'''
   seq = []
   title=None
   for line in fastafile:
      if line.startswith('>'):
         if seq:
            yield title, ''.join(seq).upper()
         title, seq = line.strip(), []
      else:
         seq.append(line.strip())
   yield   title, ''.join(seq)   

def rev_comp_maker(sequence):
   ''' returns a reverse complement of a provided DNA or RNA sequence. Any characters 
   other than standard nucleotides, 'A', 'G', 'T', 'C' and 'U' will not be changed'''
   
   nt_comp_dict={ 'a':'t','t':'a','c':'g', 'g':'c',
                  'A':'T','T':'A','C':'G', 'G':'C'   } 
   nt_comp_dict2={ 'a':'u','u':'a','c':'g', 'g':'c',
                  'A':'U','U':'A','C':'G', 'G':'C'   } 
   if 'U' in sequence.upper():
      return ''.join(nt_comp_dict2[i] if i in nt_comp_dict2 else i for i in sequence[::-1])
   
   else:
      return ''.join(nt_comp_dict[i] if i in nt_comp_dict else i for i in sequence[::-1])

def translate(fastafile):
   '''Produces translated sequence from the longest ORF on the given fasta sequence. 
   The program will brake if there is non nucleotide character'''

   outfile=open(fastafile.split('.')[0]+'_Output_orf-aa.fa','w')
   print ('\toutput file is : ' +fastafile.split('.')[0]+'_Output_orf-aa.fa')
   fastafile=open(fastafile,'r')
   
   
   aa_dict={ 
          'ATT':'I', 'ATC':'I', 'ATA':'I', 'CTT':'L', 
          'CTC':'L', 'CTA':'L', 'CTG':'L', 'TTA':'L', 
          'TTG':'L', 'GTT':'V', 'GTC':'V', 'GTA':'V', 
          'GTG':'V', 'TTT':'F', 'TTC':'F', 'ATG':'M', 
          'TGT':'C', 'TGC':'C', 'GCT':'A', 'GCC':'A', 
          'GCA':'A', 'GCG':'A', 'GGT':'G', 'GGC':'G', 
          'GGA':'G', 'GGG':'G', 'CCT':'P', 'CCC':'P', 
          'CCA':'P', 'CCG':'P', 'ACT':'T', 'ACC':'T', 
          'ACA':'T', 'ACG':'T', 'TCT':'S', 'TCC':'S', 
          'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S', 
          'TAT':'Y', 'TAC':'Y', 'TGG':'W', 'CAA':'Q', 
          'CAG':'Q', 'AAT':'N', 'AAC':'N', 'CAT':'H', 
          'CAC':'H', 'GAA':'E', 'GAG':'E', 'GAT':'D', 
          'GAC':'D', 'AAA':'K', 'AAG':'K', 'CGT':'R', 
          'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 
          'AGG':'R', 'TAA':'*', 'TAG':'*', 'TGA':'*',
          
          "tga": "*", "tct": "S", "gaa": "E", "ccg": "P",
          "gta": "V", "ttg": "L", "gga": "G", "gct": "A", 
          "tca": "S", "ggt": "G", "gag": "E", "cga": "R", 
          "cag": "Q", "acg": "T", "gtg": "V", "ata": "I", 
          "caa": "Q", "cca": "P", "cgc": "R", "ttc": "F", 
          "acc": "T", "tta": "L", "gac": "D", "att": "I", 
          "ctc": "L", "gtc": "V", "tgt": "C", "act": "T", 
          "cta": "L", "ccc": "P", "agg": "R", "cgt": "R", 
          "cct": "P", "tat": "Y", "cgg": "R", "ctt": "L", 
          "taa": "*", "cat": "H", "atg": "M", "aat": "N", 
          "tgc": "C", "gcg": "A", "tac": "Y", "aaa": "K", 
          "ggg": "G", "tcc": "S", "gcc": "A", "tcg": "S", 
          "agt": "S", "aag": "K", "agc": "S", "aac": "N", 
          "aca": "T", "tgg": "W", "ctg": "L", "cac": "H", 
          "ggc": "G", "atc": "I", "ttt": "F", "tag": "*", 
          "aga": "R", "gca": "A", "gtt": "V", "gat": "D"
          
          }

   orf=re.compile(r'(?=(ATG(?:(?!TAA|TAG|TGA)...)+(?:TAA|TAG|TGA)))', re.I)

   for title, seq in tab_generator(fastafile):
      seq=seq.upper()
      if 'U' in seq:
         seq=seq.replace ('U', 'T')
      rev_seq=rev_comp_maker(seq)
      if seq:
         orfs=orf.findall(seq)+ orf.findall(rev_seq)
         if orfs:
            long_orf=max(orfs, key=len)  
            list_of_orf=[long_orf[i:i+3] for i in range(0,len(long_orf),3)]
            aa = (''.join([aa_dict[i] for i in list_of_orf]))
            title = title+ ' | longest orf translated '  
         else: 
            title= title + ' | no orf is present in this sequence'
            aa='no orf'
      else:
         title = title + ' | no squence present for this title'
         aa='none'
      outfile.write (title+ '\n'+ '\n'.join(textwrap.wrap(aa, 70))+'\n\n')
      
   fastafile.close()
   outfile.close()
            
def swapped(fastafile):
   '''produce fasta sequence(s) with case swapped for open reading frame. 
   It is easy way visualize 5' and 3' UTR. It also adds the ORF co-ordinate 
   to the fasta title.'''
   outfile=open(fastafile.split('.')[0]+'_Output_orf-swapped.fa','w')
   print ('\toutput file is : ' +fastafile.split('.')[0]+'_Output_orf-swapped.fa')
   fastafile=open(fastafile,'r')
   
   orf=re.compile(r'(?=(ATG(?:(?!TAA|TAG|TGA)...)+(?:TAA|TAG|TGA)))', re.I)
   orf2=re.compile(r'(?=(AUG(?:(?!UAA|UAG|UGA)...)+(?:UAA|UAG|UGA)))', re.I)
         
         
   for title, seq in tab_generator(fastafile):
      
         
      rev_seq=rev_comp_maker(seq)
      if seq:
         if 'U' in seq.upper():
            orfs=orf2.findall(seq)+ orf2.findall(rev_seq)
         else:
            orfs=orf.findall(seq)+ orf.findall(rev_seq)
         
         if orfs:
            long_orf=max(orfs, key=len)  
            swapped_seq=seq.replace(long_orf, long_orf.swapcase())
            start = seq.find(long_orf)
            end =start+len(long_orf)
            title = title+ ' | CDS: ' +str(start+1) +'-'+str(end)   
            if start == -1:
               long_orf=rev_comp_maker(long_orf)
               start = seq.find(long_orf)
               end =start+len(long_orf)
               swapped_seq=seq.replace(long_orf, long_orf.swapcase())
               title =title+ ' | CDS: ' +str(start+1) +'-'+str(end) + \
               ' | orf on reverse strand'
         else: 
            title= title + ' | no orf is present in this sequence'
            swapped_seq=seq
      else:
         title = title + ' | no squence present for this title'
         swapped_seq= seq
      outfile.write (title+ '\n'+ '\n'.join(textwrap.wrap(swapped_seq, 70))+'\n\n')
  
   fastafile.close()   
   outfile.close()

def length_on_id(fastafile):
   '''produce fasta file with each sequence length added on title'''
   outfile=open(fastafile.split('.')[0]+'_Output_with-seq-length.fa', 'w')
   print ('\toutput file is : '+ fastafile.split('.')[0]+'_Output_with-seq-length.fa')
   with open(fastafile,'r') as fastafile:
      for title, seq in tab_generator(fastafile):
         title = title + ' | length = ' + str(len(seq))
         outfile.write (title+ '\n'+'\n'.join(textwrap.wrap(seq, 70))+'\n\n')
      
   outfile.close()

def fasta_extractor(fastafile, idfile):
   '''extracts fasta sequences from a file provided the titles in a idfile
   only first section of title before white space will be matched'''
   outfile = open(fastafile.split('.')[0]+'_Output_extracred-sequences.fa','w')
   print ('\toutput file is : '+ fastafile.split('.')[0]+'_Output_extracted-sequences.fa')
   
   
   with open(idfile, 'r') as ids:
      id=[]
      for line in ids:
         if line.startswith('>'):
            id.append(line.split()[0])
         else:
            id.append('>'+line.split()[0])
   with open(fastafile, 'r') as fastafile:
      for title, seq in tab_generator(fastafile):
         i=title.split()[0]
         if i in id:
            outfile.write(title +'\n'+'\n'.join(textwrap.wrap(seq, 70)) +'\n\n')

   outfile.close()
   

def sub_seq_extractor(fastafile, idfile):
   ''' extracts sub sequences of fasta sequences specified by the coordinates in id file'''
   outfile=open(fastafile.split('.')[0]+'_Output_sub-seq-extracted.fa','w')
   print ('\toutput file is : '+ fastafile.split('.')[0]+'_Output_sub-seq-extracted.fa')
   
   
   with open(idfile,'r') as ids:
      id_info={}
      for line in ids:
         try:
            x,y,z =[line.split()[0],line.split()[-2],line.split()[-1]]
         except:
            print( 'idfile for sub-seq should contain [id start end]'\
            ' separated by any white space (tab or space) and start should be >=1')
            sys.exit(2)
         if x.startswith('>'):
            id_info.update({x:(y,z)})
         else:
            id_info.update({'>'+x:(y,z)})
   with open(fastafile, 'r') as fastafile:         
      for title, seq in tab_generator(fastafile):
         i=title.split()[0]
         if i in id_info:
            start,end =id_info[i]
            outfile.write(title + ' | subseq:' + str(int(start)-1) + '-' + str(end) +'\n' +'\n'.join(textwrap.wrap(seq [int(start)-1:int(end)], 70))  +'\n\n')

   outfile.close()

def tab_delimited(fastafile):
   ''' produces tab delimited file of fasta sequences '''
   outfile=open(fastafile.split('.')[0]+'_Output_tab-delimited.txt', 'w')
   print ('\toutput file is : '+ fastafile.split('.')[0]+'_Output_tab-delimited.txt')
   with open(fastafile,'r') as fastafile:
      for title, seq in tab_generator(fastafile):
         outfile.write (title + '\t'+ seq+'\n')
      
   outfile.close()
   
def rev_comp(fastafile):
   '''produces reverse complement of RNA or DNA sequence from a fasta sequence'''
   outfile=open(fastafile.split('.')[0]+'_Output_reverse_complement.fa','w')
   print ('\toutput file is : ' +fastafile.split('.')[0]+'_Output_reverse_complement.fa')

   with open (fastafile, 'r') as fastafile: 
      for title, seq in tab_generator(fastafile):
         outfile.write (title+ ' | reverse_complement\n'+ '\n'.join(textwrap.wrap(rev_comp_maker(seq), 70))+'\n\n')
      
   outfile.close()

def nucleotide_changer(fastafile):
   ''' converts DNA to RNA and vice versa'''
   outfile=open(fastafile.split('.')[0]+'_Output_nucleotide-changed.fa','w')
   print ('\toutput file is : '+ fastafile.split('.')[0]+'_Output_nucleotide-changed.fa')
   fastafile=open(fastafile, 'r')
   for line in fastafile:
      if not line.startswith('>'):
         if ('U'in line) or ('u'in line):
            line=line.replace('u','t')
            line=line.replace('U','T')
         else:
            line =line.replace('T', 'U')
            line =line.replace('t', 'u')
      outfile.write (line)
  
   fastafile.close()   
   outfile.close()
   

      
#############################################      
      
def main(argv):

   idfile = ''
   fastafile=''
   tab=''
   sub_seq=''
   translate_seq=''
   nt=''
   fasta_ext=''
   r_c=''
   length_id=''
   swap=''
   
   
   try:
     opts,args=getopt.getopt(argv,"htsanerlwi:f:",["help","tab","subseq","aminoacid",\
     "ntchange","fastaextract","revcomp:","length","swapped","idfile=","fastafile="])
   except getopt.GetoptError:
      print ('Usage: fastaUtility.py [options: -h,-t,-s,-a,-n,-e,-r,-l,-w]'\
      ' -f <fastafile> -i <idfile> ')
      sys.exit(2)

   for opt, arg in opts:
      if opt in ("-h", "--help"):
         print ('fastautility.py [*options] -f <fastafile> -i <idfile> ')
         sys.exit()
      elif opt in ("-f", "--fastafile"):
         fastafile = arg
         print ('\tinput fasta is :', fastafile)
      elif opt in ("-i", "--idfile"):
         idfile = arg
         print ('\tId file is     :', idfile)
      elif opt in ("-t", "--tab"):
         tab=True
      elif opt in ("-s", "--subseq"):
         sub_seq=True
      elif opt in ("-a", "--aminoacid"):
         translate_seq=True
      elif opt in ("-n", "--ntchange"):
         nt=True
      elif opt in ("-e", "--fastaextract"):
         fasta_ext=True
      elif opt in ("-r", "--revcomp"):
         r_c=True
      elif opt in ("-l", "--length"):
         length_id=True
      elif opt in ("-w", "--swapped"):
         swap=True
         
   if tab:
      tab_delimited(fastafile)
   if sub_seq:
      sub_seq_extractor(fastafile,idfile)
   if translate_seq:
      translate(fastafile)
   if nt:
      nucleotide_changer(fastafile)
   if fasta_ext:
      fasta_extractor(fastafile, idfile)
   if r_c:
      rev_comp(fastafile)
   if length_id:
      length_on_id(fastafile)
   if swap:
      swapped(fastafile)



if __name__ == "__main__":
   main(sys.argv[1:])            
