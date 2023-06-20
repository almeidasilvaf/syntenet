# syntenet 0.99.0

NEW FEATURES

* Added a `NEWS.md` file to track changes to the package.

# syntenet 0.99.5

NEW FEATURES

* Replaced distance measure for clustering of phylogenomic profiles: from
vegan::vegdist to labdsv::dsvdis 

# syntenet 0.99.7

NEW FEATURES

* Added functions fasta2AAStringSetlist and gff2GRangesList to help users
read multiple FASTA and GFF/GTF files as list of AAStringSet and
GRangesList, respectively.

* Updated vignette to instruct users on how to load FASTA and GFF/GTF files
to the R session.

# syntenet 1.1.1

NEW FEATURES

* Ward's clustering of synteny clusters is now performed prior to plotting
in `plot_profiles()`, not in `phylogenomic_profile()`. As a consequence,
`phylogenomic_profile()` now returns only a matrix of profiles, not a list
containing the matrix and an hclust object.

* Added an option to handle names in vector *cluster_species* as new
names for display in the heatmap. This way, species abbreviations can be
easily replaced with species' full names to make plots look better.

* Added parameters *dist_function* and *dist_params* to allow users to specify
function and parameters to calculate the distance matrix that will be passed
to Ward's clustering.

# syntenet 1.1.2

NEW FEATURES

* Added parameters *clust_function* and *clust_params* in `cluster_network()`
to let users pass any igraph::cluster_* function to cluster the synteny network.

* Added parameters *clust_function* and *clust_params* in `plot_profiles()` to
let users have more control on the method used to cluster the distance matrix
(columns in phylogenomic profiles).

* Updated vignette to reflect the changes mentioned above and included an FAQ
item with instructions on how to update the R PATH variable.

# syntenet 1.1.3

BUG FIXES

* Tidy evaluation with aes_() was deprecated in ggplot 3.0.0, and
testthat now returns warnings for it. Replaced aes_() with aes() and .data
from the rlang package.

# syntenet 1.1.4

BUG FIXES

* Replaced sprintf calls with snprintf calls in C++ code to address warnings
in the devel branch of Bioc


UPDATES

* Added CITATION file with reference to published paper

# syntenet 1.1.5

BUG FIXES

* Fixed species ID retrieval by adding an exported function
named `create_species_id_table()` that correctly creates unique species 
IDs (3-5 characters), even when the first 5 characters are equal.

* `intraspecies_synteny()` and `interspecies_synteny()` (originally unexported) 
now take the same output of `process_input()`, which makes them consistent
with the entire package. 


NEW FEATURES

* To make it easier for users who want to run DIAMOND from the command line,
I added the functions `export_sequences()` and `read_diamond()`, which
write processed sequences to FASTA files and read the DIAMOND output, 
respectively. An example code on how to run DIAMOND from the command line
has been added to the vignette.

* Included a section in the vignette on how to use syntenet as a synteny
detection program (i.e., to find synteny within a single genome or between
two genomes).

# syntenet 1.1.6

BUG FIXES

* Added a check for which IQ-TREE version is installed. If IQ-TREE2 is 
installed (and not IQ-TREE), arguments and call are modified accordingly,
because IQ-TREE developers changed some parameters (e.g., -bb is now -B).

NEW FEATURES

* Added a parameter `as` to `parse_collinearity()` that allows the extraction
on synteny block information from .collinearity files. The vignette was 
updated accordingly.

# syntenet 1.3.1

NEW FEATURES

* Added function `collapse_protein_ids()` to replace protein IDs in sequence
names (equivalent to FASTA headers) with gene IDs. If there are multiple
protein for the same gene, onlt the longest is kept. Vignette
was updated accordingly.

# syntenet 1.3.2

BUG FIXES

* Strand information is now preserved in the output of `process_input()` (for
users who want to plot synteny).

# syntenet 1.3.3

# syntenet 1.3.4

BUG FIXES

# syntenet 1.3.5

NEW FEATURES

* Added function `run_last()` to run alternative BLAST search. Vignette
was updated accordingly.

