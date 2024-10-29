# Python script that consolidates the KEGGstand pipeline output of multiple fastas into a single file
## Info
This script is used to consolidate the KEGGstand output of multiple input fasta files into a single tab-delimited output file, which can then be processed in excel, R, or other software. 

## Input
The script takes two mandatory inputs, which are the path to the directory containing KEGGstand output, and a name/path of the output prefix. Input is specified using the following arguments 
### Mandatory:
* __-i *path/to/dir/*__: Path to the directory containing KEGGstand output.
* __-o *prefix*__: Path/name of the output prefix. The script will add "_merged" and/or "_merged_BRITE" after this prefix.
### Optional
* __--minimum *number*__: Minimum completion values, used to avoid outputting uninformative modules. Will only output modules where at least 1 input fasta is found to be above minimum completion.
* __--average__: Changes the calculation of the minimum completion values. Instead of outputting modules where any fasta is above minimum, only output modules where the average is above the mininum completion value.
* __--KEGG_cat *"category"*__: If specified, will output KEGG BRITE or pathway categories gene counts. Will give a tab-delimited output of gene counts per fasta for the specified category, and all its subcategories.
  Can alternatively specify a comma-delimited list of categories, or input "ALL", to output every single category. __Important: Put the categories within quotation marks, otherwise python will only consider the first word__

### Examples:
Default run to generate an output containing the module completions of every fasta processed in the KEGGstand output:
```
python KEGGstand_output_comparison.py -i /random_folder/my_favorite_fastas_KEGGstand_output/ -o /output_folder/my_favorite_prefix
``` 
Same run, but only including modules that are complete in at least one fasta:
```
python KEGGstand_output_comparison.py -i /random_folder/my_favorite_fastas_KEGGstand_output/ -o /output_folder/my_favorite_prefix --minimum 1
```
A run where only modules that are on average 0.5 complete are considered:
```
python KEGGstand_output_comparison.py -i /random_folder/my_favorite_fastas_KEGGstand_output/ -o /output_folder/my_favorite_prefix --minimum 0.5 --average
```
A run where the default module output is generated, and gene counts for the bacterial motility BRITE category and all its subcategories are given:
```
python KEGGstand_output_comparison.py -i /random_folder/my_favorite_fastas_KEGGstand_output/ -o /output_folder/my_favorite_prefix -KEGG_cat "Bacterial motility proteins [BR:ko02035]"
```
## Output
Outputs a tab-delimited text file (or two if KEGG_cat is specified). The first line is a header line (preceded with # to allow for removal). The first column lists the modules or BRITE/pathway categories. 
All other columns represent module completion or BRITE gene counts per input fasta. 
