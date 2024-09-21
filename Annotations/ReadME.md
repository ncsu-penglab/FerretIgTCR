# Annotations folder
## **IG_TR_10FerretAssemblies.gff3**
- Contains Walsh et. al (2023) Ig/TCR constant region annotations on the following ferret genome assemblies: ASM985912v1, ASM985917v1, ASM985921v1, ASM985922v1, ASM985923v1, ASM985924v1, ASM985926v1, M1713, ASM1036756v1, and ASM1036757v1.
- The GFF3 file contains only exon annotations, ranging from the 5' splice site to the 3' UTR, and does not provide frame information.
- Membrane-bound constant regions contain "_MB" suffix in the ID.
- Quality of previously published draft genome assemblies and allelic variation between animals were not assessed.
- The following annotations were <ins>**excluded**</ins> from the GFF3 file for the following scenarios:
   - Low mapping quality and overlapped in alignment with another C region 
      -   IGLC6*01: ASM985921v1, ASM985926v1
   - Deleted in draft assemblies **(see IGV folder for more information)**
      -   IGLC5*01: ASM985912v1, ASM985917v1, ASM985921v1, ASM985922v1, ASM985923v1, ASM985924v1, ASM985926v1, M1713, ASM1036756v1
   - Gene was collapsed during assembly and/or has low/no CDS support in draft assemblies **(see IGV folder for more information)**
      - IGLC2*01: ASM985912v1, ASM985921v1, ASM985922v1, ASM985923v1, ASM985924v1, ASM985926v1, M1713, ASM1036756v1, ASM1036757v1
- Annotations that had high mapping quality but aligned to another contig suggest potential assembly errors.
   - TRGC1*01: ASM985912v1, ASM985923v1, ASM1036756v1
   - TRGC2*01: ASM985923v1, ASM1036756v1
   - TRGC3*01: ASM1036756v1, ASM1036757v1
