lines=''
ignlist=''
#file inputs
f=open("prefix_pseudos_blast.tsv",'w')
sumprot=open("prefix_proteome.faa.blastP_output.tsv",'r')
suminter=open("prefix_intergenic.fasta.blastX_output.tsv",'r')
pseudos=open("prefix_pseudos.gff",'r')
#making list of ids
for gene in pseudos:
    if gene.startswith("#"):
        continue
    if "," in gene:
        genes=gene.split('old_locus_tag=',1)[1]
        for x in genes.split(","):
            lines+=x.rstrip()+'\n'
        continue
    lines+=gene.split('old_locus_tag=',1)[1]
lines=lines.rstrip()
lines=lines.split('\n')
for line in lines:
    if "ign" in line:
        ignlist+=line+'\n'
ignlist=ignlist.rstrip()
ignlist=ignlist.split('\n')
#extracting the blast results
for prot in sumprot.readlines():
    for line in lines:
        if line in prot:
            f.write(prot)
f.write("\n###intergenic\n")
for igngene in suminter:
    for ign in ignlist:
        if ign in igngene:
            f.write(igngene)
print('done')
