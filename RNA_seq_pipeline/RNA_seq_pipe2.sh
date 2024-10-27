#!/bin/bash
#Need 3 files: transcripts.fasta, emapper.annotations and transdecoder.pep

conda activate /home/ecbsu/mambaforge-pypy3/envs/trinotate
#Hammer to identify protein domains
echo "starting Hammer"
hmmscan --cpu 8 --domtblout $1.TrinotatePFAM.out /home/ecbsu/programs/Trinotate-Trinotate-v3.2.2/Pfam-A.hmm $1.fasta.transdecoder.pep > $1.pfam.log

#tmHMM for transmembrane proteins
echo "Starting tmHMM"
/home/ecbsu/programs/tmhmm-2.0c/bin/tmhmm --short < $1.fasta.transdecoder.pep > $1.tmhmm.out

#BLAST homologies with Trinotate
blastx -query $1.fasta -db /home/ecbsu/programs/Trinotate-Trinotate-v3.2.2/uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 > $1.uniprot.blastx
blastp -query $1.fasta.transdecoder.pep -db /home/ecbsu/programs/Trinotate-Trinotate-v3.2.2/uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 > $1.uniprot.blastp

#SignalP6 to predct signal peptides
echo "Starting Signalp"
conda activate trinotate-V4.0.2
mkdir signalp_out
signalp6 --fastafile transcripts.fasta.transdecoder.pep --format none --mode fast --model_dir /home/ecbsu/programs/signalp6/signalp6_fast/signalp-6-package/models/ -od $1.sigP6outdir

#LOAD TRANSCRIPTS AND CODING REGIONS
#Prepare mapping file
grep ">" $1.fasta | sed s"/>//" | cut -d '_' -f1,2,3,4,5,6,7 > gene_map.txt
grep ">" $1.fasta | sed s"/>//" > transcript_ids.txt
paste gene_map.txt transcript_ids.txt > gene_trans_map.txt

#Load a fresh Trinotate.sqlite from the Trinotate folder
cp /home/ecbsu/programs/Trinotate-Trinotate-v4.0.0/TrinotateBoilerplate.sqlite Trinotate.sqlite
#Load info into sqlite database
/home/ecbsu/programs/Trinotate-Trinotate-v4.0.0/Trinotate --db Trinotate.sqlite --init --gene_trans_map gene_trans_map.txt --transcript_fasta $1.fasta --transdecoder_pep $1.fasta.transdecoder.pep

echo "Prepping annotation report"
#LOAD data and create annotation report
/home/ecbsu/programs/Trinotate-Trinotate-v4.0.0/Trinotate --db Trinotate.sqlite --LOAD_pfam $1.TrinotatePFAM.out
/home/ecbsu/programs/Trinotate-Trinotate-v4.0.0/Trinotate --db Trinotate.sqlite --LOAD_swissprot_blastx $1.uniprot.blastx
/home/ecbsu/programs/Trinotate-Trinotate-v4.0.0/Trinotate --db Trinotate.sqlite --LOAD_swissprot_blastp $1.uniprot.blastp
/home/ecbsu/programs/Trinotate-Trinotate-v4.0.0/Trinotate --db Trinotate.sqlite --LOAD_tmhmmv2 $1.tmhmm.out
/home/ecbsu/programs/Trinotate-Trinotate-v4.0.0/Trinotate --db Trinotate.sqlite --LOAD_signalp $1.sigP6outdir/output.gff3
/home/ecbsu/programs/Trinotate-Trinotate-v4.0.0/Trinotate --db Trinotate.sqlite --LOAD_EggnogMapper $1.emapper.annotations
/home/ecbsu/programs/Trinotate-Trinotate-v4.0.0/Trinotate --db Trinotate.sqlite --report > ${1}_trinotate_annotation_report.tsv
echo "Annotation report ready"