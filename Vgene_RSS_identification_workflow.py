### Contact: Dr. York Ian (ite1@cdc.gov)
import Bio
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os,  os.path
import re
from subprocess import Popen, PIPE

## Output files
gf = open('V-region.fas', 'w') ### V-region nucleotide sequences for selected protein sequences
tf = open('V-region_RSS.txt', 'w') ### Tab-delimited file containing RSS heptamer and nonamer sequences for matching V-regions
tf.write("Contig ID\tStart\tEnd\tRSS-Heptamer\tRSS-Nonamer\tPutative V-region\n")

def get_vregion_with_rss_or_cdnaseq(ref_genome, gid, start_pos, end_pos, rc):
    """ Extract genome sequence for matching Vregions and include 500bps upstream of each match"""
    if re.search('\d+', str(start_pos)) and start_pos < 0:
        start_pos=0
    for gs in ref_genome:
		if gs.description.split(' ')[0] == gid:
			if rc == 'yes': ## Need to reverse complement sequence
				vh_match = SeqRecord(gs.seq[start_pos:end_pos]).reverse_complement()
				return str(vh_match.seq)
			if rc == 'no':
				vh_match = SeqRecord(gs.seq[start_pos:end_pos])
				return str(vh_match.seq)
			if rc == 'cdna': # refine cdna sequence boundaries based on splice site
				cdna = start_pos.strip()
				pc = ""
				pc = re.compile(cdna).search(str(gs.seq).upper())
				if not pc:
					rcdna = str(Bio.Seq.Seq(cdna).reverse_complement())
					pc = re.compile(rcdna).search(str(gs.seq).upper())
				if pc:
					return pc.start(), pc.end()
				else: 
					return 0, 0

def find_genome_location(prot_file_path, db_file_path):
    """ Run translated BLAST search on the input protein sequences to identify matches in the genome"""
    blast_hits = {}
    all_genome_hits = {}
    for seq_record in SeqIO.parse(prot_file_path, 'fasta'):
        
        if '|'  in seq_record.description: ### This conditional statement is customized for this input and can be changed as needed
            query_len = len(str(seq_record.seq))

            #Generate query sequence file
            handle = open('query.fas', 'w')
            SeqIO.write(seq_record, handle, "fasta")
            handle.close()
          
            tblastn_path = "../../file path for tblastn"
            
            #Run tblasn
            tblastn_cmd = "%s -query %s -db %s -qcov_hsp_perc 100 -num_threads 5 -outfmt \"6 qseqid sseqid pident length sstart send\" -max_target_seqs 5" % (tblastn_path, 'query.fas', db_file_path)
        
            blast_out = Popen([tblastn_cmd], stdout=PIPE, shell=True)
            if blast_out:
                linelist = blast_out.stdout.readlines()
                for ei, line in enumerate(linelist):
                    query_id, genome_id, perc_iden, match_len, gstart, gend=line.decode("utf-8").split('\t')
                    if perc_iden == '100.000' and match_len.strip() == str(query_len):
                        blast_hits["%d-%s=%s" % (ei, query_id, genome_id)] = "%s-%s" % (gstart, gend)
                        #print(query_id, genome_id, match_len, perc_iden, gstart, gend)

    get_genome_seqids = [all_genome_hits.update({kb.split('=')[1]:''}) for kb in blast_hits.keys()]

    #Retrieve genome matches from reference genome database
    genome_fastas = [s for s in SeqIO.parse(db_file_path,  'fasta') if s.description.split(' ')[0] in all_genome_hits.keys()]
    with open('V-region_DNA_with_RSS.fas', 'w') as fg:
        for k, v in blast_hits.items():
            qid, gid = k.split('=')
            gstart, gend = v.split('-')
            if int(gstart) > int(gend):
                start_pos = int(gend)-55
                end_pos = int(gstart)+500
                gvregion_rss = get_vregion_with_rss_or_cdnaseq(genome_fastas, gid, start_pos, end_pos, 'yes')

                ori_start_pos = int(gend)-1
                ori_end_pos = int(gstart)
                gvregion_cdna = get_vregion_with_rss_or_cdnaseq(genome_fastas, gid, ori_start_pos, ori_end_pos, 'yes')

            if int(gstart) < int(gend):
                start_pos = int(gstart)-500
                end_pos = int(gend)+55
                gvregion_rss = get_vregion_with_rss_or_cdnaseq(genome_fastas, gid, start_pos, end_pos, 'no')

                ori_start_pos = int(gstart)-1
                ori_end_pos = int(gend)
                gvregion_cdna = get_vregion_with_rss_or_cdnaseq(genome_fastas, gid, ori_start_pos, ori_end_pos, 'no')

            ## Refining cdna sequence boundaries based on splice sites
            cpos = gvregion_rss.upper().index(gvregion_cdna.upper())
            #print(cpos, len(gvregion_cdna))
            extd_cdna = gvregion_rss[cpos-10:cpos+len(gvregion_cdna)].upper()
            vcdna = extd_cdna[extd_cdna.index('AG')+2:len(extd_cdna)]
            cstart, cend = get_vregion_with_rss_or_cdnaseq(genome_fastas, gid, vcdna, 0, 'cdna')


            ### Recombination Signal Sequence (RSS) Sequence patterns for Immunoglobulin chain types. 
            ## Please select appropriate RSS extraction method,  depending on your input chain type
            ## RSS regular expressions for TRV's are yet to be included
            
            # RSS regular expression for IGKV
            rss_igk = re.compile('CAC.{17}C.{2}A')
            # RSS regular expression for IGLV and IGHV
            rss_igl = re.compile('CAC.{28}C.{2}A')
            
            hepta = " "
            nano = " "
            try:
                ### Get heptamer and nonamer sequences for IGKV ###
                #rss_region = gvregion_rss.strip()[-55:].upper()
                #find_rss = rss_igk.search(rss_region)
                #hepta = rss_region[find_rss.start():find_rss.start()+7]
                #nano = rss_region[find_rss.start()+19:find_rss.start()+28]

                ### Get heptamer and nonamer sequences for IGLV and IGHV ###
                rss_region = gvregion_rss.strip()[-55:].upper()
                find_rss = rss_igl.search(rss_region)
                hepta = rss_region[find_rss.start():find_rss.start()+7]
                nano = rss_region[find_rss.start()+30:find_rss.start()+39]
            except:
                pass


            fg.write(">%s\n%s\n" % (k, gvregion_rss.strip().upper()))
            gf.write(">%s\n%s\n" % (k, vcdna.strip().upper()))
            ## Customized tabulated output for RSS##
            tf.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (gid, str(cstart).strip(), str(cend).strip(), hepta, nano, vcdna.strip().upper()))

    gf.close()
    tf.close()


def read_files():
    print("   ================================================================ \n \
    This script takes a V-region protein sequence output from Vgenextract tool (Olivieri et al.,  2013) \n \
    and corresponding reference genome sequences in FASTA format. A Blast database of reference \n \
    genome sequence is a prerequisite for this script.\n \
    ================================================================ \n")

    for i in range(2):
        if i == 0:
            filename = input(f"Enter the name of V-region protein sequence  file: ")
        if i == 1:
            filename = input(f"Enter  the name of reference genome sequence  file: ")
            
        ### Check if the file exists
        if not os.path.isfile(filename):
            print(f"The file '{filename}' was not found.")
            break
        if i == 0:
            protein_seqs = filename.strip()
        if i == 1:
            reference_seqs = filename.strip()
    try:               
        find_genome_location(protein_seqs, reference_seqs) ### Get genome sequence for Vgene matches
    except Exception as e:
        print(f"An error occured: {e}")



if __name__ == "__main__":
    """ This script takes immunoglobulin heavy or light chain V-region protein sequence and its
        reference genome sequence as input.
        The script identifies genomic sequences that match the V-region protein input,  finding the 
        heptamer and nonamer Recombination Signal Sequence (RSS) within these matches,  and
        also extracting the 500bp sequence upstream of each V-region in genome. These sequences are
        all included in the output file. 
        The RSS sequences are included in the output file V-region_RSS.txt."""
    read_files()
        
        
