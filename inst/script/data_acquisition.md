Data acquisition
================

# Data in data/

## proteomes.rda

The protein sequences of the *Ostreococcus sp.* algae species (from
primary transcripts only) were obtained from Pico-PLAZA 3.0 and stored
in a list of AAStringSet object.

``` r
links <- c(
    Olucimarinus = "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/Fasta/proteome.selected_transcript.olu.fasta.gz",
    Osp_RCC809 = "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/Fasta/proteome.selected_transcript.orcc809.fasta.gz",
    Otauri = "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/Fasta/proteome.selected_transcript.ota.fasta.gz"
)

proteomes <- lapply(links, Biostrings::readAAStringSet)

usethis::use_data(
    proteomes, compress = "xz"
)
```

## annotation.rda

``` r
links <- c(
    Olucimarinus = "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/GFF/olu/annotation.selected_transcript.all_features.olu.gff3.gz",
    Osp_RCC809 = "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/GFF/orcc809/annotation.selected_transcript.all_features.orcc809.gff3.gz",
    Otauri = "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/GFF/ota/annotation.selected_transcript.all_features.ota.gff3.gz"
)

annotation <- lapply(links, rtracklayer::import)
annotation <- GenomicRanges::GRangesList(annotation)

usethis::use_data(
    annotation, compress = "xz"
)
```

# Data in inst/extdata
