Blobtools2, https://blobtoolkit.genomehubs.org/ has been set up on deigo.
blobtools2_create.slurm is used to setup the data for analysis on the interactive viewer.

Files you have to generate beforehand: 
- hits (see formating details here: https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/adding-data-to-a-dataset/adding-hits/)
- coverage (sorted.mapped.bam)

If you have multiple hits files, edit the slurm script line 20

---

Visualisation of the data can be done on the Dell workstation 

```
conda activate btk
blobtools view --interactive file_name (edited) 

```
