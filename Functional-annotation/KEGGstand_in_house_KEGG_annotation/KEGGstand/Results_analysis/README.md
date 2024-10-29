# Python script that consolidates the KEGGstand pipeline output of multiple fastas into a single file
## Info
These scripts are used to merge the KEGGstand output of multiple samples. These can be merged into a single .tsv file for downstream processing, or immediately visualized as a heatmap.

## Input
The scripts are highly similar, but due to their different outputs can take slightly different arguments:

### KEGGstand_tsv_maker.py
#### Mandatory
* __-i *path/to/dir/*__: Path to the directory containing KEGGstand output.
* __-o *prefix*__: Path/name of the output prefix. The script will add ".tsv" and/or "_BRITE.tsv" after this prefix.
  
#### Optional
* __--completion *number*__: Desired completion values, used to avoid outputting uninformative modules. By default will only consider modules where at least 1 sample is above this required completion.
* __--filter_method *method*__: Changes the calculation of the desired completion values. Options are "min", "max", and "avg". "min" is the default behaviour, requiring the lowest completion to be above or equal to the desired completion. "max" specifies the inverse, with the highest value needing to be equal or higher. "avg" specifies that the average module completion needs to be higher or equal to the specified completion
* __--module_search *string*__: Specifies a string that needs to be PRESENT in the module name. Modules WITHOUT this string are not considered. Can give multiple strings by delimiting with a comma. If spaces are included the string should be within quotation marks
* __--module_filter *string*__: Specifies a string that needs to be ABSENT in the module name. Modules WITH this string are not considered.  Can give multiple strings by delimiting with a comma. If spaces are included the string should be within quotation marks
* __--collapse *number*__: Specifies that the modules are instead listed by overarching categories (eg rather than showing the Citrate cycle and Pentose phosphate pathway completions as separate, their results are combined into a Central carbohydrate metabolism category)
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
