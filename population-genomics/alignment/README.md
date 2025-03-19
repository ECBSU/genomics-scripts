## Sequence alignment using Bowtie2

Read mapping with Bowtie2
- local mapping (performs soft-clipping on suspected adapter sequences)
- --very-sensitive-local setting
- included loops to avoid aligning samples that had already been aligned

Post mapping data preprocessing
- sort BAM with samtools (alphabetically with collate, by genome position with sort)
- fixmate fills in missing details among paired reads
- mark duplicate reads
- index and generate flagstat report


Alignment and BAM file preprocessing based on:
James Reeve, Amin Ghane, Pierre Barry, Alfonso Balmori-de la Puente, Roger K. Butlin, Le Qin Choo, Alan Le Moan, Diego Garica Castillo, Ana-Maria Peris Tamayo, Sarah Kingston, Erica Leder, Sean Stankowski, Basile Pajot 2024. A standard pipeline for processing short-read sequencing data from Littorina snails. protocols.io
https://dx.doi.org/10.17504/protocols.io.dm6gp3m21vzp/v3
Version created by James Reeve
