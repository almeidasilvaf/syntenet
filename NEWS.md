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

# syntenet 1.0.1

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

# syntenet 1.0.2

NEW FEATURES

* Added parameters *clust_function* and *clust_params* in `cluster_network()`
to let users pass any igraph::cluster_* function to cluster the synteny network.

* Added parameters *clust_function* and *clust_params* in `plot_profiles()` to
let users have more control on the method used to cluster the distance matrix
(columns in phylogenomic profiles).

* Updated vignette to reflect the changes mentioned above and included an FAQ
item with instructions on how to update the R PATH variable.

# syntenet 1.0.3

BUG FIXES

* Tidy evaluation with aes_() was deprecated in ggplot 3.0.0, and
testthat now returns warnings for it. Replaced aes_() with aes() and .data
from the rlang package.



