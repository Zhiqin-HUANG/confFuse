# confFuse
The confFuse script and artefact list are available for download (https://github.com/Zhiqin- HUANG/confFuse), and an annotation file from the GENCODE website is also required (http://www.gencodegenes.org/releases/19.html). The Release 19 (GRCh37.p13) comprehensive gene annotation GTF file was used in this study (ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_sc aff.annotation.gtf.gz). A folder should then be created which includes only deFuse output files (results.filtered.tsv). confFuse can parse multiple deFuse result files within a folder in parallel and finish the computation within a few minutes.

Usage: python confFuse.py path_to_deFuse_folder path_to_output_folder artefact_list.tab gencode.v19.chr_patch_hapl_scaff.annotation.gtf


The deFuse results from 7 pilocytic astrocytoma RNA-seq samples are available for download in folder test_data.



