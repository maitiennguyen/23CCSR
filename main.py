from blast import *
from annotation import *

def main(argv):
	
	# BLAST
	
	seq_query = "ndc10_scerevisiae.fasta"
	ds_query = "scerevisiae.faa"
	blastp_evalue = "1e-01"  
	tblastn_evalue = "1e-01" 
	tblastx_evalue = "1e-01"
	blastx_evalue = "1e-01"
	taxIDS = ["147537"]
	download = "no"
	q_type = "prot"
	
	nucl_fasta_file = None 
	prot_fasta_file = None 
	all_specs = None
	prot_specs = None
	
	# if this is the first run and fasta files need to be downloaded	
	if download == "yes":
		# get user tax id input
		taxID_list = taxIDS

		# get protein and nucleotide fasta files
		print("Compiling nucl and aa datasets of targets into fasta files...")
		nucl_fasta_file, prot_fasta_file, all_specs, prot_specs = get_fasta_files(taxID_list)
		print("Done")

		# write all specs name to txt file
		write_list("all_specs", all_specs)

		# write all specs with protein dataset to txt file
		write_list("prot_data_specs", prot_specs)

	# if this is not the first run and/or user already have these 4 files
	else:							   	
		nucl_fasta_file = "nucl.fna"
		prot_fasta_file = "prot.faa"
		all_specs = read_list("all_specs.txt")
		prot_specs = read_list("prot_data_specs.txt")

	
	if q_type == "prot":
		
		blastp_hit_dict = {}
		tblastn_hit_dict = {}
		
		# make subject amino acid database for reciprical blasts 
		print("\n\nCompiling query aa dataset into BLAST database...")
		subj_db = get_dbs(ds_query, q_type)

		# perform blastp using provided prot seq and provided aa db to confirm qseq_id in aa db
		print("Done")
		qseq_id = get_qseq_id(seq_query, subj_db, q_type)

		
		# check if there is a protein fasta file
		if prot_fasta_file is not None:

			# get protein blast database
			print("\n\nCompiling target aa dataset into BLAST database...")
			prot_db = get_dbs(prot_fasta_file, "prot")
			print("Done")

			# perform blastp of the query sequence against protein database, get dict of seq id and its info
			print("\n\nPerforming blastp...")
			blastp_hit_dict, blastp_rslt_file = blast_aa_ds(seq_query, prot_db, q_type, blastp_evalue)
			print("Done")

			# write blast results to txt file
			write_dict("blastp1", blastp_hit_dict)

			# reverse blasto tp confirm results from 
			if len(blastp_hit_dict) > 0:
				print("\n\nPerforming reciprocal blastp...")
				blastp_rev_list = recip_blast("blastp", blastp_hit_dict, subj_db, qseq_id)
				print("Done")

				# update blastp_hit_dict after validation
				# if 0 seqs were valid, make dict empty
				if len(blastp_rev_list) == 0:
					blastp_hit_dict = {}
				# if len of valid seqs != len of dict, then dict must be updated
				elif len(blastp_rev_list) != len(blastp_hit_dict):
					update_blast_dict(blastp_hit_dict, blastp_rev_list)

				# write updated blast results to txt file
				write_dict("blastp2", blastp_hit_dict)

				# write summary txt file
				write_summary("blastp", blastp_hit_dict, prot_specs, all_specs)

				# check if there is a nucleotide fasta file
				if nucl_fasta_file is not None:
					# remove blastp hit species from nucleotide fasta file, update nucl_fasta_file
					print("\n\nRemoving species from target nucl dataset...")
					nucl_fasta_file = rm_seqs_fasta(blastp_hit_dict, nucl_fasta_file)
					print("Done")


		# check if there is a nucleotide fasta file
		if nucl_fasta_file is not None:

			# get nucleotide blast database
			print("\n\nCompiling target nucl dataset into BLAST database...")
			nucl_db = get_dbs(nucl_fasta_file, "nucl")
			print("Done")

			# perfrom tblasn of query sequence against nucleotide database, get dict of seq id and its info
			print("\n\nPerforming tblastn...")
			tblastn_hit_dict, tblastn_rslt_file = blast_nucl_ds(seq_query, q_type, nucl_db, tblastn_evalue)
			print("Done")

			# write blast results to txt file
			write_dict("tblastn", tblastn_hit_dict)

			# perform blastx on each tblastn result, keep valid results
			if len(tblastn_hit_dict) > 0:
				print("\n\nPerforming reciprocal blastx...")
				blastx_rev_list = recip_blast("blastx", tblastn_hit_dict, subj_db, qseq_id)
				print("Done")

				# update tblastn_hit_dict after validation
				# if 0 seqs were valid, make dict empty
				if len(blastx_rev_list) == 0:
					tblastn_hit_dict = {}
				# if len of valid seqs != len of dict, then dict must be updated
				elif len(blastx_rev_list) != len(tblastn_hit_dict):
					update_blast_dict(tblastn_hit_dict, blastx_rev_list)

				# write blast results to txt file
				write_dict("blastx", tblastn_hit_dict)

				# write summary txt file
				write_summary("tblastn", tblastn_hit_dict, blastp_hit_dict.values(), all_specs)
				
				
		# ANNOTATION
		
		annotation = BlastAnnot(tblastn_hit_dict)
		no_gap_dict, gap_dict = annotation.process_seqs()
		write_dict("no_gap", no_gap_dict)
		write_dict("gap", gap_dict)
		
		for seq_id in no_gap_dict.keys():
			stop_five, start, stop_three = annotation.find_codons(seq_id, nucl_db)
	
	elif q_type == "nucl":
		
		blastx_hit_dict = {}
		tblastx_hit_dict = {}
		
		# make subject nucl database for reciprical blasts
		print("\n\nCompiling query nucl dataset into BLAST database...")
		subj_db = get_dbs(ds_query, q_type)
		print("Done")

		# perform blastp using provided prot seq and provided aa db to confirm qseq_id in aa db
		qseq_id = get_qseq_id(seq_query, subj_db, q_type)
		
		# check if there is a protein fasta file
		if prot_fasta_file is not None:

			# get protein blast database
			print("\n\nCompiling target aa dataset into BLAST database...")
			prot_db = get_dbs(prot_fasta_file, "prot")
			print("Done")

			# perform blastx of the query sequence against protein database, get dict of seq id and its info
			print("\n\nPerforming blastx...")
			blastx_hit_dict, blastx_rslt_file = blast_aa_ds(seq_query, q_type, prot_db, blastx_evalue)
			print("Done")

			# write blast results to txt file
			write_dict("blastx", blastx_hit_dict)

			# reverse blasto tp confirm results from 
			if len(blastx_hit_dict) > 0:
				print("\n\nPerforming reciprocal tblastn...")
				tblastn_rev_list = recip_blast("tblastn", blastx_hit_dict, subj_db, qseq_id)
				print("Done")

				# update blastx_hit_dict after validation
				# if 0 seqs were valid, make dict empty
				if len(tblastn_rev_list) == 0:
					blastx_hit_dict = {}
				# if len of valid seqs != len of dict, then dict must be updated
				elif len(tblastn_rev_list) != len(blastx_hit_dict):
					update_blast_dict(blastx_hit_dict, tblastn_rev_list)

				# write updated blast results to txt file
				write_dict("tblastn", blastx_hit_dict)

				# write summary txt file
				write_summary("blastx", blastx_hit_dict, prot_specs, all_specs)

				# check if there is a nucleotide fasta file
				if nucl_fasta_file is not None:
					# remove blastx hit species from nucleotide fasta file, update nucl_fasta_file
					print("\n\nRemoving species from target nucl dataset...")
					nucl_fasta_file = rm_seqs_fasta(blastx_hit_dict, nucl_fasta_file)
					print("Done")
		
		# check if there is a nucleotide fasta file
		if nucl_fasta_file is not None:

			# get nucleotide blast database
			print("\n\nCompiling target nucl dataset into BLAST database...")
			nucl_db = get_dbs(nucl_fasta_file, "nucl")
			print("Done")

			# perfrom tblasx of query sequence against nucleotide database, get dict of seq id and its info
			print("\n\nPerforming tblastx...")
			tblastx_hit_dict, tblastx_rslt_file = blast_nucl_ds(seq_query, q_type, nucl_db, tblastx_evalue)
			print("Done")

			# write blast results to txt file
			write_dict("tblastx1", tblastx_hit_dict)

			# perform blastx on each tblastn result, keep valid results
			if len(tblastx_hit_dict) > 0:
				print("\n\nPerforming reciprocal tblastx...")
				tblastx_rev_list = recip_blast("tblastx", tblastx_hit_dict, subj_db, qseq_id)
				print("Done")

				# update tblastn_hit_dict after validation
				# if 0 seqs were valid, make dict empty
				if len(tblastx_rev_list) == 0:
					tblastn_hit_dict = {}
				# if len of valid seqs != len of dict, then dict must be updated
				elif len(tblastx_rev_list) != len(tblastx_hit_dict):
					update_blast_dict(tblastx_hit_dict, tblastx_rev_list)

				# write blast results to txt file
				write_dict("tblastx2", tblastx_hit_dict)

				# write summary txt file
				write_summary("tblastx", tblastx_hit_dict, blastx_hit_dict.values(), all_specs)
	
	
		

if __name__ == "__main__":
	main(sys.argv)