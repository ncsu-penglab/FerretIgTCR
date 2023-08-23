#/opt/R/4.1.0/bin/R

###########################################################################################################
#	Contact Dr. Peng (xpeng5@ncsu.edu)
#	
#	Objective: use most-up-to-date variable region reference to search for more variable region genes 
#		1) use most-up-to-date variable region reference to mask previously annotated positions within the genome  (create new fa refernce file)
#		2) align variable region refernce to the masked reference assebmly with relaxed alignment parameters
#		3) QC data by removing bad sequence alignments
#		4) merge filtered bed file 
#		4) extract 50 bp downstream of each new, putative sequence
#		5) score potential 28 bp RSS sequences using the probabilities calculated with the originally IDed RSSs
#		6) score potential 39 bp RSS sequences using the probabilities calculated with the originally IDed RSSs
#		7) search putative sequences for the presence of a stop codon 
#		8) join bed file, rss data, and orf data to filter sequences 
###########################################################################################################
	
	## Set width for viewing and plotting
	# options(width = Sys.getenv("COLUMNS"))
	
	## Install packages as needed
	.cran_packages <- c("tidyverse", "dplyr", "stringr", "gdata",'ggpubr','parallel','SciViews') 
	.inst <- .cran_packages %in% installed.packages()
	if(any(!.inst)) install.packages(.cran_packages[!.inst], repos='http://cran.us.r-project.org')
	
	## Load packages
	sapply(.cran_packages, require, character.only = TRUE)
	
	## Tools
	bedtools <- "/opt/bedtools/2.30.0/bin/bedtools"
	minimap2 <- "/peng_1/peng_lab/tools/seqanal/minimap2-2.22_x64-linux/minimap2"
	samtools <- "/peng_1/peng_lab/tools/seqanal/samtools-1.13/samtools"
	TRANSEQ <- "/mnt/peng_1/peng_lab/tools/seqanal/EMBOSS-6.6.0/bin/bin/transeq"
	NTHREADS <- 4
	
	## Functions 
	Fasta_edit <- function(fasta,rss_typ){ ## edit fasta file to create a sequences of length equal to rss_typ
		fasta.l <- list()
		## determine the last start position to make N sequences of length=rss_typ
		seq_length <- fasta$seq %>% nchar %>% unique
		last_start <- seq_length - rss_typ + 1
		
		for (i in 1:nrow(fasta)){
			seqid <- fasta[i,'id',drop=T] ## sequence name 
			
			for (j in 1:last_start){
				seq.fix <- fasta[i,'seq',drop=T] %>% substr(., j, (j+rss_typ-1)) ## subset of sequence with length = rss_typ
				fasta.l[[length(fasta.l)+1]] <- data.frame(id = paste0(seqid, '::', j), seq = seq.fix) ## join the sequence id (plus identifier) and sequence with length = rss_typ
			}
		}
		
		fasta.new <- bind_rows(fasta.l)
		return(fasta.new)
	}
	
	## I/O
	root <- '/peng_1/peng_lab/results/final/Ferret_germline_anno'
	outroot <- file.path(root, "search"); if (!dir.exists(outroot)) dir.create(outroot)
	outdir <- file.path(outroot, 'minimap2'); if (!dir.exists(outdir)) dir.create(outdir)
	
	refdir <- "/peng_1/peng_lab/seq_idx/Mustela_putorius_furo_PacBio"
	refID <- 'GCA_011764305.2_ASM1176430v1.1_genomic.fna'
	
	################################################################################################################################################
	## 1) use most-up-to-date variable region reference to mask previously annotated positions within the genome (create new fa refernce file)
	
	print("################################################################")
	print("## 1) use most-up-to-date variable region reference to mask previously annotated positions within the genome (create new fa refernce file)")
	
	## get variable region reference (look for most recent in output directory, if none use the original refernce)
	bedfile.name <- 'variable.sortedbyCoord.bed'
	bedfiles <- list.files(outdir, bedfile.name, recursive=T, full.names=T)
	if (length(bedfiles) > 0){
		# index <- bedfiles %>% dirname %>% basename %>% sub('search','',.) %>% as.numeric %>% sort %>% .[length(.)]
		index <- length(bedfiles)
		bed.file <- grep(paste0('search',index),bedfiles, value=T)
		org.bed.file <- bed.file ## save bed file as another object to be joined with new data later  
	} else {
		bed.file <- file.path(root, 'RSS_ferret_probs', "variable.sortedbyCoord.bed")
		org.bed.file <- bed.file ## save bed file as another object to be joined with new data later  
	}
	
	ref.fastafile <- file.path(refdir,'GCA_011764305.2_ASM1176430v1.1_genomic.fna');
	masked.output <- file.path(outdir, 'masked_fastafiles');if (!dir.exists(masked.output)) dir.create(masked.output)
	num.files <- list.files(masked.output) %>% length ## count the number of files weve already masked to keep them in order  
	masked.file <- file.path(masked.output, paste0((num.files+1), '_', 'GCA_011764305.2.fa'))
	
	## Usage: $ bedtools maskfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF> -fo <output FASTA>
	cmd <- paste(bedtools, 'maskfasta -fi', ref.fastafile, '-bed', bed.file, '-fo', masked.file)
	system(cmd, intern=F)
	
	## if we have already performed one search, then check if there are any sequences that were removed
	## mask them as well 
	if (length(bedfiles) > 0){
		index <- length(bedfiles)
		search.num <- paste0('search',index)
		remove.bedfile <- file.path(outdir,search.num,'metadata','variable2remove.bed')
		
		if (file.size(remove.bedfile) > 0L){ ## in the bed file contains any entries mask them 
			outfile <- file.path(masked.output,'tmp.bed')
			cmd <- paste(bedtools, 'maskfasta -fi', masked.file, '-bed', remove.bedfile, '-fo', outfile)
			system(cmd, intern=F)
			
			## update masked.file 
			cmd <- paste('mv',outfile, masked.file); system(cmd,intern=F)
			
		}
	}
	
	################################################################################################################################################
	## 2) align variable region refernce to the masked reference assebmly with relaxed alignment parameters
	
	print("################################################################")
	print("## 2) align variable region refernce to the masked reference assebmly with relaxed alignment parameters")
	
	fastafile.name <- 'variable.sortedbyCoord.fa'
	fastafiles <- list.files(outdir, fastafile.name, recursive=T, full.names=T)
	if (length(fastafiles) > 0){
		# index <- fastafiles %>% dirname %>% basename %>% sub('search','',.) %>% as.numeric %>% sort %>% .[length(.)]
		index <- length(fastafiles)
		query.fastafile <- grep(paste0('search',index),fastafiles, value=T)
		
	} else {
		query.fastafile <- file.path(root, 'RSS_ferret_probs', "variable.sortedbyCoord.fa")
	}
	
	## create a directory in outdir for this specific search
	num.files <- length(fastafiles)
	searchdir <- file.path(outdir, paste0('search', num.files+1)); if (!dir.exists(searchdir)) dir.create(searchdir)
	sam.file <- file.path(searchdir, 'minimap2.sam')
	cmd <- paste(minimap2, "-ax splice:hq -uf -p 0.1 -N 10 ", masked.file, query.fastafile, ">", sam.file) # PacBio 
	
	if (!file.exists(sam.file)){
		system(cmd, intern=F)
		## sort sam file 
		cmd <- paste(samtools, "sort -u -o", sub('minimap2.sam','minimap2_sort.sam',sam.file),"-O sam", sam.file); system(cmd, intern=F)
		cmd <- paste("mv", sub('minimap2.sam','minimap2_sort.sam',sam.file), sam.file); system(cmd, intern=F)
		
	}
	
	## convert to a bed file
	bedfile <- file.path(searchdir, 'minimap2.bed12')
	cmd <- paste(samtools, "view -bS", sam.file, "|", bedtools, "bamtobed -bed12 -i - >", bedfile)
	system(cmd, intern=F)
	
	if (file.size(bedfile) == 0L){ ## if no new putative sequences quit the program
		quit()
	}
	
	################################################################################################################################################
	## 3) QC data by removing bad sequence alignments
	
	print("################################################################")
	print("## 3) QC data by removing bad sequence alignments")
	
	bed.dat <- read.table(bedfile, sep='\t', header=F, stringsAsFactors=F) %>% as_tibble %>% setNames(c('chrom','start', 'end','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts'))
	
	## remove sequences with more or less than 2 exons 
	bed.dat <- bed.dat %>% filter(blockCount == 2)
	
	## remove sequences with more than 1000 bp between exon start positions 
	bed.dat <- bed.dat %>% mutate(tmp = sub('0,','',blockStarts) %>% as.numeric) %>% filter(tmp < 1000) %>% dplyr::select(-tmp)
	
	## remove sequences with first exons greater than 200 bp 
	t2 <- bed.dat %>% filter(blockCount == 2) %>% dplyr::select(name, strand, blockSizes) %>% separate_rows(blockSizes, sep=',') %>% group_by(name) %>% mutate(id = row_number()) %>% ungroup %>% filter(id %in% c(1,2))
	tmp <- expand.grid(c('+','-'), c(1,2), stringsAsFactors=F) %>% setNames(c('strand','id'))
	tmp$exon <- c('exon1','exon2','exon2','exon1')
	t2 <- full_join(t2, tmp)
	seqs2remove <- t2 %>% mutate(blockSizes = as.numeric(blockSizes)) %>% filter(exon == 'exon1' & blockSizes > 200) %>% .$name
	bed.dat <- bed.dat %>% filter(!(name %in% seqs2remove))
	
	if (nrow(bed.dat) == 0){ ## if no new putative sequences quit the program
		quit()
	}
	
	## save the filtered bed file
	filtered.bedfile <- file.path(searchdir,'minimap2_filtered.bed')
	write.table(bed.dat, file=filtered.bedfile, sep='\t', col.names=F, row.names=F, quote=F)
	
	################################################################################################################################################
	## 4) merge filtered bed file 
	
	print("################################################################")
	print("## 4) merge filtered bed file ")
	
	## use bedtools merge to collapse overlapping sequences (need to convert to bed 6 format first)
	merged.bed.file <- file.path(searchdir, 'minimap2_merged.bed')
	cmd <- paste('cat', filtered.bedfile, "|", bedtools, "bed12tobed6 -i stdin | sort -k1,1 -k2,2n - |", bedtools, "merge -s -o distinct,mean,distinct -c 4,5,6 -i - >", merged.bed.file)
	system(cmd, intern=F)
	
	## read in merged bed file to simplify naming scheme 
	bed.dat <- read.table(merged.bed.file, sep='\t', header=F, stringsAsFactors=F) %>% as_tibble %>% setNames(c("chrom", 'start', 'end', 'name', 'score', 'strand'))  %>% mutate(locus = substr(name, 1, 3))  %>% group_by(name) %>% mutate(locus = paste0(locus, cur_group_id())) %>% ungroup# %>% dplyr::select(-name) 
	
	## remove any sequences that were observed an odd number of times 
	## expect to see at least 2 (2 exons) but may see more from additional alignments
	bed.dat <- bed.dat %>% group_by(name) %>% mutate(n=n()) %>% filter(n %% 2 == 0) %>% ungroup %>% dplyr::select(-name,-n)
	
	
	## position sort the bed file (exon 1 and exon 2 of same gene will be next to each other)
	bed.dat <- bed.dat[order(bed.dat$chrom, bed.dat$start),]
	bed.dat <- bed.dat %>% mutate(rowID = row_number()) %>% mutate(seqID = ntile(x=rowID, n=nrow(.)/2)) %>% mutate(locus = paste0(locus,'_',seqID)) %>%  dplyr::select(chrom, start, end, locus, score, strand) %>% dplyr::rename(name = locus)
	
	merged.bed.file <- file.path(searchdir, 'minimap2_merged_rename.bed')
	write.table(bed.dat, file=merged.bed.file, sep='\t', col.names=F, row.names=F, quote=F)
	
	###############################################################
	## convert bed6 file to bed12 file 
	outfile <- file.path(searchdir,'minimap2_merged_12.bed')
	cmd <- paste('/peng_1/peng_lab/scripts/Ferret_germline_anno/leader_seq/bed6tobed12.sh -i', merged.bed.file, '>', outfile)
	system(cmd,intern=F)
	
	###############################################################
	## read in merged bed12 file to QC data 
	bed.dat <- read.table(outfile, sep='\t', header=F, stringsAsFactors=F) %>% as_tibble %>% setNames(c('chrom','start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'))
	## remove sequences that did not merge properly
	bed.dat <- bed.dat %>% filter(blockCount > 1)
	write.table(bed.dat, outfile, sep='\t', col.names=F, row.names=F, quote=F) ## save over the file 
	
	
	################################################################################################################################################
	## 4) extract 50 bp downstream of each new, putative sequence
	
	print("################################################################")
	print("## 4) extract 50 bp downstream of each new, putative sequence")
	rss.dat.l <- list()
	for (k in 1:nrow(bed.dat)){
		strand <- bed.dat[k, 'strand', drop=T]
		start_pos <- bed.dat[k, 'start', drop=T]; end_pos <- bed.dat[k, 'end', drop=T]
		id <- bed.dat[k,'name', drop=T]; contig <- bed.dat[k,'chrom', drop=T]
		if (strand == '-'){
			## get heptamer and nonamer coordinates 
			range.end <- start_pos; range.start <- range.end - 50 
			
			## bind data into bed format 
			.dat <- cbind(contig, range.start, range.end, id, strand) %>% as.data.frame(., stringsAsFactors=F) %>% setNames(c('chrom','start','end','name','strand'))
			rss.dat.l[[length(rss.dat.l) + 1]] <- .dat 
			
		} else {
			
			range.start <- end_pos ; range.end <- range.start + 50
			
			## bind data into bed format 
			.dat <- cbind(contig, range.start, range.end, id, strand) %>% as.data.frame(., stringsAsFactors=F) %>% setNames(c('chrom','start','end','name','strand'))
			rss.dat.l[[length(rss.dat.l) + 1]] <- .dat 
			
		}
	}
	
	## save RSS bed data
	dat.file <- file.path(searchdir, 'RSS.bed')
	bed.dat <- rss.dat.l %>% bind_rows %>% as_tibble %>% mutate(score =60) %>% dplyr::select(chrom,start,end,name,score,strand)
	write.table(bed.dat, file=dat.file, sep='\t', col.names=F, row.names=F, quote=F)
	
	## get RSS fasta 
	# p_val <- strsplit(bed.file, split='/') %>% unlist %>% .[9]
	genome.fasta <- file.path(refdir, "GCA_011764305.2_ASM1176430v1.1_genomic.fna")
	outfile <- file.path(searchdir, 'RSS.fa')
	cmd <- paste(bedtools, 'getfasta -s -nameOnly -fi', genome.fasta, '-bed', dat.file, '>', outfile)
	system(cmd, intern=F)
	
	
	################################################################################################################################################
	## 5) score potential 28 bp RSS sequences using the probabilities calculated with the originally IDed RSSs
	
	print("################################################################")
	print("## 5) score potential 28 bp RSS sequences using the probabilities calculated with the originally IDed RSSs")
	
	fastafile <- outfile; rss_typ=28
	fasta <- read.table(fastafile, sep='\t', header=F, stringsAsFactors=F)
	fasta <- do.call("cbind", split(fasta, rep(c(1, 2), length.out = nrow(fasta)))) 
	colnames(fasta) <- c('id','seq')
	fasta <- fasta %>% as.data.frame %>% as_tibble %>% mutate(seq = gsub('a','A',seq) %>% gsub('t','T',.) %>% gsub('c','C',.) %>% gsub('g','G',.))
	
	## check what the length of the input sequences are 
	seq_length <- fasta$seq %>% nchar %>% unique
	if (n_distinct(seq_length) > 1){
		stop('Stop, sequences are not of equal length. Exiting!')
	}
	
	## overwrite the rss_typ parameter if the sequence lengths are equal to one of the two types (28 or 39(
	rss_typ <- ifelse(seq_length %in% c(28,39), yes=seq_length, no=rss_typ)
	
	## edit the fasta file to sequences of length = rss_typ if they are not already that length 
	if (seq_length != rss_typ){
		fasta <- Fasta_edit(fasta, rss_typ) %>% as_tibble #%>% mutate(id = sub('_[0-9]+$','',id))
	}
	
	## record the shift from the gene of V gene sequence 
	fasta <- fasta %>% separate(id, into=c('id','shift'), sep='::') %>% mutate(shift = as.numeric(shift)-1)
	
	
	## remove sequences that are not IGK and start with CA 
	fasta <- fasta %>% filter(grepl("IGK",id) & grepl("^CA",seq)) %>% mutate(score=0)
	
	if (nrow(fasta) > 0){
		prob.file <- file.path(root,'RSS_ferret_probs','rss12_probs.rds')
		rss.probs <- readRDS(prob.file)
		for (i in 1:nrow(fasta)){
			seq. <- fasta$seq[i] %>% gsub('a','A',.) %>% gsub('t','T',.) %>% gsub('c','C',.) %>% gsub('g','G',.)
			
			## initialize a list of probabilities that will be selected 
			prob.list <- list()
			
			## split the string into individual characters 
			seq.split <- strsplit(seq.,'') %>% unlist 
			for (j in 1:length(rss.probs)){
				seq.pos <- names(rss.probs)[j] %>% strsplit(split=',') %>% unlist %>% as.numeric
				rss.prob <- rss.probs[[j]]
				
				## create the sequence of interest 
				seq.oi <- seq.split[seq.pos] %>% paste0(collapse='')
				prob.list[[length(prob.list)+1]] <- rss.prob[rss.prob$nts==seq.oi,]
			}
			
			## add score to fasta data frame 
			prob.dat <- bind_rows(prob.list) %>% mutate(ln_Freq = ln(prob))
			score = sum(prob.dat$ln_Freq)
			fasta[i,'score'] <- score
			
		}
		
		fasta <- fasta %>% group_by(id) %>% filter(score == max(score))
		rss_score_12 <- fasta 
	} else {
		rss_score_12 <- NULL
	}
	
	################################################################################################################################################
	## 6) score potential 39 bp RSS sequences using the probabilities calculated with the originally IDed RSSs
	
	print("################################################################")
	print("## 6) score potential 39 bp RSS sequences using the probabilities calculated with the originally IDed RSSs")
	
	fastafile <- outfile; rss_typ=39
	fasta <- read.table(fastafile, sep='\t', header=F, stringsAsFactors=F)
	fasta <- do.call("cbind", split(fasta, rep(c(1, 2), length.out = nrow(fasta)))) 
	colnames(fasta) <- c('id','seq')
	fasta <- fasta %>% as.data.frame %>% as_tibble %>% mutate(seq = gsub('a','A',seq) %>% gsub('t','T',.) %>% gsub('c','C',.) %>% gsub('g','G',.))
	
	## check what the length of the input sequences are 
	seq_length <- fasta$seq %>% nchar %>% unique
	if (n_distinct(seq_length) > 1){
		stop('Stop, sequences are not of equal length. Exiting!')
	}
	
	## overwrite the rss_typ parameter if the sequence lengths are equal to one of the two types (28 or 39(
	rss_typ <- ifelse(seq_length %in% c(28,39), yes=seq_length, no=rss_typ)
	
	## edit the fasta file to sequences of length = rss_typ if they are not already that length 
	if (seq_length != rss_typ){
		fasta <- Fasta_edit(fasta, rss_typ) %>% as_tibble #%>% mutate(id = sub('_[0-9]+$','',id))
	}
	
	## record the shift from the gene of V gene sequence 
	fasta <- fasta %>% separate(id, into=c('id','shift'), sep='::') %>% mutate(shift = as.numeric(shift)-1)
	
	## remove sequences that are not IGK and start with CA 
	fasta <- fasta %>% filter(!(grepl("IGK",id)) & grepl("^CA",seq)) %>% mutate(score=0)
	if (nrow(fasta) > 0){
		
		prob.file <- file.path(root,'RSS_ferret_probs','rss23_probs.rds')
		rss.probs <- readRDS(prob.file)
		for (i in 1:nrow(fasta)){
			seq. <- fasta$seq[i] %>% gsub('a','A',.) %>% gsub('t','T',.) %>% gsub('c','C',.) %>% gsub('g','G',.)
			
			## initialize a list of probabilities that will be selected 
			prob.list <- list()
			
			## split the string into individual characters 
			seq.split <- strsplit(seq.,'') %>% unlist 
			for (j in 1:length(rss.probs)){
				seq.pos <- names(rss.probs)[j] %>% strsplit(split=',') %>% unlist %>% as.numeric
				rss.prob <- rss.probs[[j]]
				
				## create the sequence of interest 
				seq.oi <- seq.split[seq.pos] %>% paste0(collapse='')
				prob.list[[length(prob.list)+1]] <- rss.prob[rss.prob$nts==seq.oi,]
			}
			
			## add score to fasta data frame 
			prob.dat <- bind_rows(prob.list) %>% mutate(ln_Freq = ln(prob))
			score = sum(prob.dat$ln_Freq)
			fasta[i,'score'] <- score
			
		}
		
		fasta <- fasta %>% group_by(id) %>% filter(score == max(score))
		rss_score_23 <- fasta 
	} else {
		rss_score_23 <- NULL
	}
	
	## classify rss sequences based on their RIC score 
	if (!is.null(rss_score_12)) rss_score_12 <- rss_score_12 %>% mutate(RSS = ifelse(score > -38.81, yes=T, no=F))
	if (!is.null(rss_score_23)) rss_score_23 <- rss_score_23 %>% mutate(RSS = ifelse(score > -58.45, yes=T, no=F))
	rss_scores <- bind_rows(rss_score_12, rss_score_23) %>% mutate(id = sub('\\(.*','',id) %>% sub('>','',.))
	
	outfile <- file.path(searchdir,'rss_data.txt')
	write.table(rss_scores, file=outfile, sep='\t', col.names=T, row.names=F, quote=F)
	
	################################################################################################################################################
	## 8) search putative sequences for the presence of a stop codon 
	
	print("################################################################")
	print("## 8) search putative sequences for the presence of a stop codon ")
	
	## save updated bed file 
	# bed.file <- file.path(searchdir, 'minimap2.filtered.bed12')
	# write.table(bed.dat, file=bed.file, sep='\t', col.names=F, row.names=F, quote=F)
	
	## get fasta file of filtered sequences 
	bed.file <- file.path(searchdir,'minimap2_merged_12.bed')
	fastaout <- sub('bed$','fa', bed.file)
	cmd <- paste(bedtools, 'getfasta -s -split -nameOnly -fi', genome.fasta, '-bed', bed.file, '>', fastaout)
	system(cmd, intern=F)
	
	## fix naming in fasta file for transeq 
	cmd <- paste("sed -i 's/(.*//'", fastaout); system(cmd, intern=F)
	# cmd <- paste("sed -i 's/|/_/g'", fastaout); system(cmd, intern=F)
	
	## run translation 
	for (frame in 1:3) {
		outfile <- file.path(searchdir, paste0('ORF_', frame, ".fasta"))
		# cmd <- paste(TRANSEQ, "-sequence", fastaout, "-outseq", outfile, "-frame", frame, "-auto -clean")
		cmd <- paste(TRANSEQ, "-sequence", fastaout, "-outseq", outfile, "-frame", frame)
		system(cmd)
	}
	
	## Read in ORFs, convert to single line FASTA and trim off sequence starting from stop codon
	orf.files <- list.files(searchdir, pattern = "^ORF_[1-3].fasta", full.names = T)
	orf.dat <- NULL; seq <- NULL
	for (orf.file in orf.files){
		orf.rl <- readLines(orf.file)
		for (i in 1:length(orf.rl)) {
			tmp <- orf.rl[i]
			if (substr(tmp, 1, 1) == ">") {
				if (i > 1) {
					orf.dat <- rbind(orf.dat, c(id, seq))
				}
				id <- tmp
				seq <- NULL
			} else {
				seq <- paste0(seq, tmp)
			}
		}
		orf.dat <- rbind(orf.dat, c(id, seq))
	}
	
	colnames(orf.dat) <- c("id", "seq")
	orf.dat <- as.data.frame(orf.dat, stringsAsFactors = F) %>% as_tibble
	orf.dat$frame <- as.factor(substr(orf.dat$id, nchar(orf.dat$id), nchar(orf.dat$id)))
	orf.dat$id <- sub("_[1-9]*$", "", orf.dat$id) # Clean up seq id
	
	orf.dat <- orf.dat %>% mutate(id = sub('^>','',id))
	
	## analyze only ORFs with a start codon (methionine)
	orf.dat <- orf.dat %>% filter(grepl('M',seq)) %>% mutate(seq = sub('.*?M','',seq))
	
	## check if a stop codon is present 
	# orf.dat$stop_codon <- grepl("X", orf.dat$seq)
	orf.dat$stop_codon <- grepl("\\*", orf.dat$seq)
	
	## trim off the X at the end of the sequence that signifies an unknown sequence (remainder of 1 or 2)
	orf.dat <- orf.dat %>% mutate(seq = sub('X$','',seq))
	
	## measure the length of the orf before and after removing the sequence after the stop codon 
	orf.dat <- orf.dat %>% mutate(seq_length = nchar(seq)) %>% mutate(orf = sub('\\*.*$','',seq)) %>% mutate(orf_length = nchar(orf)) %>% dplyr::select(id, seq, orf, seq_length, orf_length, frame, stop_codon)
	
	## pick the longest ORF per sequence 
	orf.dat <- orf.dat %>% group_by(id) %>% filter(orf_length == max(orf_length)) %>% ungroup %>% mutate(length_diff = seq_length - orf_length) %>% ungroup
	
	## use stop codon presence and position to classfiy sequences 
	orf.dat <- orf.dat %>% mutate(imgt = ifelse(length_diff == 0, yes="ORF", no=ifelse(length_diff == 1, yes='F', no='P')))
	
	## save data 
	outfile <- file.path(searchdir, 'orf_analysis.txt')
	write.table(orf.dat, file=outfile, sep='\t', col.names=T, row.names=F,quote=F)
	
	################################################################################################################################################
	## 9) join bed file, rss data, and orf data to filter sequences 
	
	print("################################################################")
	print("## 9) join bed file, rss data, and orf data to filter sequences ")
	
	## create a directory to store file output ("metadata"0
	metadir <- file.path(searchdir, 'metadata'); 
	if (!dir.exists(metadir)) dir.create(metadir)
	
	## read in bed file 
	bed.file <- file.path(searchdir, 'minimap2_merged_12.bed')
	bed.dat <- read.table(bed.file, sep='\t', header=F, stringsAsFactors=F) %>% setNames(c('chrom','start', 'end','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts')) %>% as_tibble
	bed6 <- bed.dat[,1:6]
	
	## read in rss data 
	rss.file <- file.path(searchdir,'rss_data.txt')
	rss.scores <- read.table(rss.file, sep='\t', header=T, stringsAsFactors=F)  %>% as_tibble
	
	## read in orf data 
	orf.file <-file.path(searchdir, 'orf_analysis.txt')
	orf.dat <- read.table(orf.file, sep='\t', header=T, stringsAsFactors=F) %>% as_tibble
	
	## join data and summarize by chain type (locus)
	dat <- bed6 %>% full_join(., rss.scores %>% dplyr::select(id, shift, RSS), by=c('name'='id')) %>% mutate(RSS = ifelse(is.na(RSS), yes=F, no=RSS))
	dat <- dat %>% full_join(., orf.dat %>% dplyr::select(id, imgt), by=c('name'='id')) %>% mutate(search = basename(searchdir) %>% sub('search','',.) %>% as.numeric) %>% mutate(FILTER = ifelse(RSS, yes=T, no=ifelse(imgt %in% c('F','ORF'), yes=T, no=F)))
	sum.dat <- dat %>% mutate(locus=substr(name, 1, 3)) %>% group_by(locus, FILTER) %>% summarize(n=n()) %>% mutate(search = basename(searchdir) %>% sub('search','',.) %>% as.numeric)
	
	## save output 
	write.table(dat, file=file.path(metadir, 'metadat.txt'), col.names=T, row.names=F, quote=F)
	write.table(sum.dat, file=file.path(metadir, 'summary.txt'), col.names=T, row.names=F, quote=F)
	
	## save the sequences that were removed (won't be added to reference)
	## will still mask out so we do not detect again 
	seqs2keep <- dat %>% filter(FILTER) %>% .$name
	bed.remove <- bed.dat %>% filter(!(name %in% seqs2keep))
	outfile <- file.path(metadir, 'variable2remove.bed')
	write.table(bed.remove, file=outfile, col.names=F, row.names=F, sep='\t', quote=F)
	if (basename(searchdir) != 'search1'){ ## if this is not the first search add the previous search removed sequences to mask 
		search.prev <- sub('search','',basename(searchdir)) %>% as.numeric %>% sum(.,-1) %>% paste0('search',.) %>% file.path(outdir,.)
		remove.prev <- file.path(search.prev,'metadata','variable2remove.bed')
		
		cmd <- paste("cat", remove.prev, '>>', outfile); system(cmd, intern=F)
		
	}
	
	## save the filter bed and fasta files 
	bed.dat <- bed.dat %>% filter(name %in% seqs2keep)
	
	## join bed file with previous version 
	org.bed <- read.table(org.bed.file, sep='\t', header=F, stringsAsFactors=F) %>% setNames(c('chrom','start', 'end','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts')) %>% as_tibble
	bed.dat <- bed.dat %>% bind_rows(., org.bed)
	
	## position sort the bed file and save 
	bed.dat <- bed.dat[order(bed.dat$chrom, bed.dat$start),] 
	bed.file <- file.path(metadir,'variable.sortedbyCoord.bed')
	write.table(bed.dat, bed.file, sep='\t', col.names=F, row.names=F, quote=F)
	
	fastaout <- file.path(metadir, 'variable.sortedbyCoord.fa')
	cmd <- paste(bedtools, 'getfasta -s -split -nameOnly -fi', genome.fasta, '-bed', bed.file, '>', fastaout)
	system(cmd, intern=F)
	
	
	
