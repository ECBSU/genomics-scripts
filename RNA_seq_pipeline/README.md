RNAseq analysis pipeline. 
---

This utilises both deigo and dell workstation.

Computation will include, Busco, Eggnog, Pfam, tmHMMv2, signalp6, blastp, blastx, (edit include) targetp.

Assembly and some annotations are done on Deigo. Other annotations and collation of data is done on the Dell workstation. 

In the slurm script (RNA_seq_pipe1.slurm), modify the inputs of SPAdes and database of BUSCO

Once this is finished, 
Annotations on Dell. will require 3 files to be copied to the Dell: 
- transcripts.fasta
- emapper.annotations
- transdecoder.pep
  
Run the script with the prefix of your fasta file as the second argument 

For example my file name is transcripts.fasta
```
RNA_seq_pipe2.sh transcripts
```
