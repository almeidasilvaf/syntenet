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

