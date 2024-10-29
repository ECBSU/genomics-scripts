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
* __--collapse *number*__: Specifies that the modules are instead to be listed by their overarching categories (eg rather than showing the Citrate cycle and Pentose phosphate pathway completions as separate, their results are combined into a Central carbohydrate metabolism category). The number represents the level of detail for the categories. 1 is the broadest categories (eg energy metabolism), 2 is slightly more specific categories (eg methane metabolism). *If this is specified, --db needs to be specified as well.*
* __--db *path/to/module/database*__: Path to the module database used to generate the results. Used to find the categories to which the modules belong. (/bucket/HusnikU/Unit_Members/Arno/Scripts/KEGG_module_db on Deigo)
* __--collapse_method *method*__: Specifies the value to be given for the collapsed categories. Default is the average completion of the contained modules (avg). Alternatives are "min" (lowest value) and "max" (highest value) 
* __--category_search *string*__: Functions like the --module_search argument, but instead of only considering modules with the specified string, will only consider categories with the specified string.
* __--category_filter *string*__: Functions like the --module_filter argument, but instead of disregarding modules with the specified string, will disregard the categories with the specified string.
* __--show_module_count__: Specifies that the number of contained modules are to be listed in the name of the categories. Only does anything if --collapse is given.
#### TSV_maker only
* __--KEGG_cat *"category"*__: If specified, will output KEGG BRITE or pathway categories gene counts. Will give a tab-delimited output of gene counts per fasta for the specified category, and all its subcategories.
  Can alternatively specify a comma-delimited list of categories, or input "ALL", to output every single category. *Put the categories within quotation marks, otherwise python will only consider the first word*

### KEGGstand_graph_maker.py
__IMPORTANT:__ The KEGGstand_graph_maker.py script requires the seaborn package and its dependencies. If used on Deigo, the easiest way to acquire these is to activate a conda environment with these packages present. I have been using: 
```
conda activate /bucket/HusnikU/Conda-envs/genovi
```
#### Mandatory
* __-i *path/to/dir/*__: Path to the directory containing KEGGstand output.
* __-o *outputname*__: Path/name of the output heatmap. __The extension used here will affect the format of the heatmap__. Supported formats are pdf, png, ps, eps, and svg. Default is .png. 
#### Optional
* __--completion *number*__: Desired completion values, used to avoid outputting uninformative modules. By default will only consider modules where at least 1 sample is above this required completion.
* __--filter_method *method*__: Changes the calculation of the desired completion values. Options are "min", "max", and "avg". "min" is the default behaviour, requiring the lowest completion to be above or equal to the desired completion. "max" specifies the inverse, with the highest value needing to be equal or higher. "avg" specifies that the average module completion needs to be higher or equal to the specified completion
* __--module_search *string*__: Specifies a string that needs to be PRESENT in the module name. Modules WITHOUT this string are not considered. Can give multiple strings by delimiting with a comma. If spaces are included the string should be within quotation marks
* __--module_filter *string*__: Specifies a string that needs to be ABSENT in the module name. Modules WITH this string are not considered.  Can give multiple strings by delimiting with a comma. If spaces are included the string should be within quotation marks
* __--collapse *number*__: Specifies that the modules are instead to be listed by their overarching categories (eg rather than showing the Citrate cycle and Pentose phosphate pathway completions as separate, their results are combined into a Central carbohydrate metabolism category). The number represents the level of detail for the categories. 1 is the broadest categories (eg energy metabolism), 2 is slightly more specific categories (eg methane metabolism). *If this is specified, --db needs to be specified as well.*
* __--db *path/to/module/database*__: Path to the module database used to generate the results. Used to find the categories to which the modules belong. (/bucket/HusnikU/Unit_Members/Arno/Scripts/KEGG_module_db on Deigo)
* __--collapse_method *method*__: Specifies the value to be given for the collapsed categories. Default is the average completion of the contained modules (avg). Alternatives are "min" (lowest value) and "max" (highest value) 
* __--category_search *string*__: Functions like the --module_search argument, but instead of only considering modules with the specified string, will only consider categories with the specified string.
* __--category_filter *string*__: Functions like the --module_filter argument, but instead of disregarding modules with the specified string, will disregard the categories with the specified string.
* __--show_module_count__: Specifies that the number of contained modules are to be listed in the name of the categories. Only does anything if --collapse is given.
#### Graph_maker only
* __--dimensions *number,number*__: Desired dimensions of the created plot in inches. First value is height, the second is width.
* __--color *color_scheme*__: Desired colorscheme of the heatmap. Default is "coolwarm". Check here for options: https://www.practicalpythonfordatascience.com/ap_seaborn_palette

## Output
Tsv_maker will make a tab-delimited text file (or two if KEGG_cat is specified). The first line is a header line (preceded with # to allow for removal). The first column lists the modules or BRITE/pathway categories. 
All other columns represent module completion or BRITE gene counts per input fasta. 

Graph_maker will make a heatmap in the desired image file format. 

## Notes on collapsing modules into categories
Please be aware that collapsing the modules into their overarching categories occurs AFTER any module filtering. In other words, the values calculated for the categories will be affected by any module filtering steps. If you removed any modules below an average of 0.5 completion, these modules will be absent when calculating the averages of their overarching categories. This was done so uninteresting modules can be removed. To keep track of how filtering is affecting the modules underlying the categories, use the "--show_module_count" flag. 

## Examples:
Default run to generate an output containing the module completions of every fasta processed in the KEGGstand output:
```
python KEGGstand_tsv_maker.py -i /random_folder/my_favorite_fastas_KEGGstand_output/ -o /output_folder/my_favorite_prefix
``` 
Same run, but only including modules that are complete in at least one fasta:
```
python  KEGGstand_tsv_maker.py -i /random_folder/my_favorite_fastas_KEGGstand_output/ -o /output_folder/my_favorite_prefix --completion 1
```
A run where only modules that are on average 0.5 complete are considered:
```
python KEGGstand_tsv_maker.py -i /random_folder/my_favorite_fastas_KEGGstand_output/ -o /output_folder/my_favorite_prefix --completion 0.5 --filter_method avg
```
A run where the default module output is generated, and gene counts for the bacterial motility BRITE category and all its subcategories are given:
```
python KEGGstand_tsv_maker.py -i /random_folder/my_favorite_fastas_KEGGstand_output/ -o /output_folder/my_favorite_prefix -KEGG_cat "Bacterial motility proteins [BR:ko02035]"
```
A run for generating a heatmap of the samples, containing completion values for all modules. 
```
python KEGGstand_graph_maker.py -i /random_folder/my_favorite_fastas_KEGGstand_output/ -o /output_folder/heatmap.pdf
```
A run for generating a heatmap, considering only modulesthat are above 0.5 completion in every sample, collapsing them into their broadest categories, and removing the categories including the word Xenobiotics. 
```
python KEGGstand_graph_maker.py -i /random_folder/my_favorite_fastas_KEGGstand_output/ -o /output_folder/heatmap.pdf --collapse 1 --completion 0.5 --db Scripts/KEGG_module_db --filter_method min --show_module_count --category_filter Xenobiotics
```
