# IGV folder

## Purpose:
Demonstrate validation of IGL gene annotations on the ferret draft genome assemblies.

## Methods:
The IGL locus genomic sequence (flanking regions were added for a total of 50kb) from the draft genome assemblies were aligned to the ferret reference genome using minimap2 (`-ax map-ont`).
Output minimap2 BAM files and the constant gene annotations on the ferret reference genome were uploaded to IGV.

## Result:
Alignments show gaps between `IGLC1*01` and `IGLC2*01` and over `IGLC5*01`. The gap completely covers the `IGLC5*01` gene, supporting its deletion in the draft assemblies. Howevever, draft assemblies cover the `IGLC1*01` coding sequence and part of the `IGLC2*01` 3'UTR. This suggests IGLC genes were potentially collapsed during assembly. Since the draft assemblies cover the `IGLC1*01` coding sequence, `IGLC1*01` is annotated in this position and `IGLC2*01` is deleted.

### Draft genome assemblies IGL locus genomic sequences aligned to the ferret reference genome
![alt text](https://github.com/ncsu-penglab/FerretIgTCR/blob/main/Annotations/IGV/DraftAssembliesIGLGenomic2Reference.png)
#### The top track is reference genome assembly, the middle tracks are the draft assemblies BAM files, and then the bottom track is the constant region annotations on the reference genome GFF3 file. Blue, green, and red blocks show different haplotypes among the draft assemblies based on IGL gene retention/deletion events.

### Zoomed in images of draft genome assemblies IGL locus genomic sequences coverage of `IGLC1*01` and `IGLC2*01`
![alt text](https://github.com/ncsu-penglab/FerretIgTCR/blob/main/Annotations/IGV/IGLC1*01_DraftAssemblies2Reference_91424.png)
![alt text](https://github.com/ncsu-penglab/FerretIgTCR/blob/main/Annotations/IGV/IGLC2*01_DraftAssemblies2Reference_91424.png)
#### The top track is reference genome assembly, the middle tracks are the draft assemblies BAM files, and then the bottom track is the constant region annotations on the reference genome GFF3 file.
