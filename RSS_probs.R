#/opt/R/4.1.0/bin/R

###########################################################################################################
#	@author Evan Walsh
#	@file /peng_1/peng_lab/scripts/Ferret_germline_anno_dj/final/8_RSS_probs.R
#	
#	Objective: ensure that R code RIC calculations are equal to the RIC calcuations from RSSsite
#		1) use the sequences from the RSSsite ferret data to measure nucleotide frequences and calculate Bayesian probabilities
#		2) create an RSS search region fasta file (50 bp) downstream to look for RSS 
#		3) calculate RIC scores for IGK variable sequences based on mouse rss12 sequences 
#		4) calculate RIC scores for other (!IGK) variable sequences based on mouse rss23 sequences 
#		5) compare R code RIC scores to RSSsite RIC scores 
#		
#	Source: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2002-3-12-research0072#Tab2
###########################################################################################################
	
	## Set width for viewing and plotting
	options(width = Sys.getenv("COLUMNS"))
	
	## Install packages as needed
	.cran_packages <- c("tidyverse", "stringr", "gdata",'ggpubr', 'gtools','SciViews') 
	.inst <- .cran_packages %in% installed.packages()
	if(any(!.inst)) install.packages(.cran_packages[!.inst], repos='http://cran.us.r-project.org')
	
	## Load packages
	sapply(.cran_packages, require, character.only = TRUE)
	
	## Tools
	bedtools="/opt/bedtools/2.30.0/bin/bedtools"
	
	## Functions 
	kmer.generation <- function(k){
		bases <- c("A","T","C","G")
		t1 <- permutations(n = length(bases), v = bases, r = k, repeats.allowed = T)
		return(apply(t1,1,FUN=function(x) paste0(x,collapse='')))
	}
	
	Fasta_edit <- function(fasta,rss_typ){ ## edit fasta file to create a sequences of length equal to rss_typ
		fasta.l <- list()
		## determine the last start position to make N sequences of length=rss_typ
		seq_length <- fasta$seq %>% nchar %>% unique
		last_start <- (seq_length - rss_typ) + 1
		
		for (i in 1:nrow(fasta)){
			seqid <- fasta[i,'id',drop=T] ## sequence name 
			
			for (j in 1:last_start){
				seq.fix <- fasta[i,'seq',drop=T] %>% substr(., j, (j+rss_typ-1)) ## subset of sequence with length = rss_typ
				fasta.l[[length(fasta.l)+1]] <- data.frame(id = paste0(seqid, '_', j), seq = seq.fix) ## join the sequence id (plus identifier) and sequence with length = rss_typ
			}
		}
		
		fasta.new <- bind_rows(fasta.l)
		return(fasta.new)
	}
	
	## I/O
	root <- '/peng_1/peng_lab/results/final/Ferret_germline_anno_dj'
	indir <- file.path(root, 'RSS')
	outdir <- file.path(root, 'RSS_ferret_probs')
	if (!dir.exists(outdir)) dir.create(outdir)
	
	# refdir <- "/peng_1/peng_lab/seq_idx/Mustela_putorius_furo_PacBio"
	# refID <- 'GCA_011764305.2_ASM1176430v1.1_genomic.fna'
	
	####################################################################################################
	## 1) load the predicted RSS data and format to be analyzed in R
	
	############################################
	## joining 
	tab.file <- list.files(indir, pattern='RSSpredict', full.names=T, recursive=T) %>% sort %>% .[length(.)]
	if (!file.exists(tab.file)) stop('Cannot analyze RSS prediction data...file does not exist!')
	tab <- read.table(tab.file, sep='\t', header=F, stringsAsFactors=F, fill =T) %>% as_tibble %>% setNames(c('seqname','start','end','RSSseq','strand','score','RIC'))
	
	header.lines <- which(!complete.cases(tab)); dat.l <- list()
	for (i in 1:length(header.lines)){
		if (i < length(header.lines)){
			start.pos <- header.lines[i]; end.pos <- header.lines[i+1] -1 
		} else {
			start.pos <- header.lines[i]; end.pos <- nrow(tab)
		}
		
		typ <- tab[start.pos,'seqname',drop=T]
		print(typ)
		
		## classifiy sequences based on RIC score 
		tmp <- tab[(start.pos+1):end.pos,] %>% mutate(typ = typ) %>% mutate(RIC = ifelse(RIC == 'PASS', yes=T, no=F))
		
		## filter data to only contain highest score per seqname
		tmp <- tmp %>% group_by(seqname) %>% filter(score == max(score)) %>% ungroup
		
		dat.l[[typ]] <- tmp
	}
	
	tab <- dat.l %>% bind_rows %>% filter(RIC) %>% group_by(seqname) %>% filter(score == max(score))
	
	## filter discordant hits (12 RSS hit but 23 RSS chain type)
	tab12.j <- tab %>% filter(grepl('12',typ)) %>% filter(!grepl('IGK|IGH',seqname))
	tab23.j <- tab %>% filter(grepl('23',typ)) %>% filter(grepl('IGK|IGH',seqname))
	# tab.j <- bind_rows(tab12, tab23)
	
	############################################
	## variable region 
	tab.file <- "/peng_1/peng_lab/results/final/Ferret_germline_anno/RSS/RSS_predict.txt"
	if (!file.exists(tab.file)) stop('Cannot analyze RSS prediction data...file does not exist!')
	tab <- read.table(tab.file, sep='\t', header=F, stringsAsFactors=F, fill =T) %>% as_tibble %>% setNames(c('seqname','start','end','RSSseq','strand','score','RIC'))
	
	header.lines <- which(!complete.cases(tab)); dat.l <- list()
	for (i in 1:length(header.lines)){
		if (i < length(header.lines)){
			start.pos <- header.lines[i]; end.pos <- header.lines[i+1] -1 
		} else {
			start.pos <- header.lines[i]; end.pos <- nrow(tab)
		}
		
		typ <- tab[start.pos,'seqname',drop=T]
		print(typ)
		
		## classifiy sequences based on RIC score 
		tmp <- tab[(start.pos+1):end.pos,] %>% mutate(typ = typ) %>% mutate(RIC = ifelse(RIC == 'PASS', yes=T, no=F))
		
		## filter data to only contain highest score per seqname
		tmp <- tmp %>% group_by(seqname) %>% filter(score == max(score)) %>% ungroup
		
		dat.l[[typ]] <- tmp
	}
	
	tab <- dat.l %>% bind_rows %>% filter(RIC) %>% group_by(seqname) %>% filter(score == max(score))
	
	## filter discordant hits (12 RSS hit but 23 RSS chain type)
	tab12.v <- tab %>% filter(grepl('12',typ)) %>% filter(grepl('IGK',seqname))
	tab23.v <- tab %>% filter(grepl('23',typ)) %>% filter(!grepl('IGK',seqname))
	# tab.v <- bind_rows(tab12.v, tab23.v)
	
	tab12 <- bind_rows(tab12.j, tab12.v)
	tab23 <- bind_rows(tab23.j, tab23.v)
	
	###################################################################################################
	## 2) create probabilit tables of nulceotide combinations at specific positions of interest 
	# RIC12 = ln [P1 P2 P3,15,25 P4,5 P6,28 P7,8,19 P9,26 P10,12 P11,27 P13,14,23 P16,17,18 P20,21,22 P24]
	# RIC23 = ln [P1 P2 P3 P4,14 P5,39 P6 P7,24,25 P8,9,21 P10,16 P11,12 P13,22 P15,23 P17,18 P19,27,30,31,32,33,37 P20,26 P28,29 P34,38 P35,36].
	
	## create a matrix that separates RSS sequences (RSS position X RSS sequence)
	matrix12 <- tab12 %>% ungroup %>% dplyr::select(RSSseq) %>% mutate(RSSseq = gsub('a','A',RSSseq) %>% gsub('c','C',.) %>% gsub('t','T',.) %>% gsub('g','G',.)) %>% .$RSSseq %>% lapply(., FUN=function(x) strsplit(x, split='') %>% unlist) %>% do.call(cbind,.)
	
	rss12.list <- list(1, 2, c(3,15,25), c(4,5), c(6,28), c(7,8,19), c(9,26), c(10,12), c(11,27), c(13,14,23), c(16,17,18), c(20,21,22), 24)
	rss12.probs <- list()
	for (i in 1:length(rss12.list)){
		rss.pos <- rss12.list[[i]] ## get rss position(s) of interest
		kmer.len <- length(rss.pos) ## get the numb of positions 
		kmers <- kmer.generation(kmer.len) ## generate all possible kmers of that length 
		
		## creat a dataframe of kmers 
		dat <- kmers %>% as.data.frame(., stringsAsFactors=F) %>% setNames('nts') %>% as_tibble 
		
		## calculate the freq of each kmer in RSS data 
		quant <- matrix12[rss.pos,,drop=F] %>% apply(., 2, FUN=function(x) paste0(x, collapse='')) %>% table %>% as.data.frame %>% setNames(c('nts','Count'))
		dat <- full_join(dat, quant) %>% mutate(Count = ifelse(is.na(Count), yes=0, no=Count)) %>% mutate(Freq = Count/sum(Count))
		
		## order data frame in decsending order 
		dat <- dat[order(dat$Count, decreasing=T),]
		
		## calculate Bayesian probability
		dat <- dat %>% mutate(N = sum(Count), r=4^length(rss.pos)) %>% mutate(prob = (Count + (2/r)) / (N+2))
		rss12.probs[[paste0(rss.pos, collapse=',')]] <- dat
	}
	
	## save prob tables 
	saveRDS(rss12.probs, file=file.path(outdir,'rss12_probs.rds'))
	
	## create a matrix that separates RSS sequences (RSS position X RSS sequence)
	matrix23 <- tab23 %>% ungroup %>% dplyr::select(RSSseq) %>% mutate(RSSseq = gsub('a','A',RSSseq) %>% gsub('c','C',.) %>% gsub('t','T',.) %>% gsub('g','G',.)) %>% .$RSSseq %>% lapply(., FUN=function(x) strsplit(x, split='') %>% unlist) %>% do.call(cbind,.)
	
	rss23.list <- list(1, 2, 3, c(4,14), c(5,39), 6, c(7,24,25), c(8,9,21), c(10,16), c(11,12), c(13,22), c(15,23), c(17,18), c(19,27,30,31,32,33,37), c(20,26), c(28,29), c(34,38), c(35,36))
	rss23.probs <- list()
	for (i in 1:length(rss23.list)){
		rss.pos <- rss23.list[[i]] ## get rss position(s) of interest
		kmer.len <- length(rss.pos) ## get the numb of positions 
		kmers <- kmer.generation(kmer.len) ## generate all possible kmers of that length 
		
		## creat a dataframe of kmers 
		dat <- kmers %>% as.data.frame(., stringsAsFactors=F) %>% setNames('nts') %>% as_tibble 
		
		## calculate the freq of each kmer in RSS data 
		quant <- matrix23[rss.pos,,drop=F] %>% apply(., 2, FUN=function(x) paste0(x, collapse='')) %>% table %>% as.data.frame %>% setNames(c('nts','Count'))
		dat <- full_join(dat, quant) %>% mutate(Count = ifelse(is.na(Count), yes=0, no=Count)) %>% mutate(Freq = Count/sum(Count))
		
		## order data frame in decsending order 
		dat <- dat[order(dat$Count, decreasing=T),]
		
		## calculate Bayesian probability
		dat <- dat %>% mutate(N = sum(Count), r=4^length(rss.pos)) %>% mutate(prob = (Count + (2/r)) / (N+2))
		
		rss23.probs[[paste0(rss.pos, collapse=',')]] <- dat
	}
	
	## save prob tables 
	saveRDS(rss23.probs, file=file.path(outdir,'rss23_probs.rds'))
	
	###################################################################################################################
	## 2) create an RSS search region fasta file (50 bp) downstream to look for RSS 
	# rss.bed.file <- file.path(indir, 'RSS_joining.bed'); 
	# genome.fasta <- file.path(refdir, refID); outfile <- file.path(outdir,'RSSsearch.fa')
	# cmd <- paste(bedtools, 'getfasta -s -nameOnly -fi', genome.fasta, '-bed', rss.bed.file, '>', outfile)
	# system(cmd, intern=F)
	