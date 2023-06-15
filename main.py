import sys
from blast import *
from annotation import *
from clustal import *
from process_specs import *

def main(argv):
	
	# BLAST
	
	seq_query = "skp1p_scereviseae.fasta"
	ds_query = "scerevisiae.faa"
	blastp_evalue = "1e-01"  
	tblastn_evalue = "1e-01" 
	tblastx_evalue = "1e-01"
	blastx_evalue = "1e-01"
	taxIDS = ["147537"]
	download = "no"
	q_type = "prot"
	q_spec_name = "Saccharomyces cerevisiae"
	n = 1
	
	nucl_fasta_file = None 
	prot_fasta_file = None 
	all_specs = None
	prot_specs = None
	prot_file_paths = None
	
	prot_db = None
	nucl_db = None
	prot_dict = {}
	nucl_dict = {}
	
	
	# if this is the first run and fasta files need to be downloaded	
	if download == "yes":
		# get user tax id input
		taxID_list = taxIDS

		# get protein and nucleotide fasta files, all specs txt file, specs w/ prot ds file, and dict associates specs with prot files
		print("\n\nCompiling nucl and aa datasets of targets into fasta files...\n\n")
		nucl_fasta_file, prot_fasta_file, all_specs, prot_specs, prot_file_paths = get_fasta_files(taxID_list)
		print("Done")

		# write all specs name to txt file
		write_list("all_specs", all_specs)

		# write all specs with protein dataset to txt file
		write_list("prot_data_specs", prot_specs)
		
		# write dict that ssociates specs with prot files
		write_dict("prot_files", prot_file_paths, "all")

	# if this is not the first run and/or user already have these 4 files
	else:							   	
		nucl_fasta_file = "nucl.fna"
		prot_fasta_file = "prot.faa"
		all_specs = read_list("all_specs.txt")
		prot_specs = read_list("prot_data_specs.txt")
		prot_file_paths = txt_to_dict("prot_files_all_dict.txt")

	# the number of times to run queries
	for i in range(n):
		
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
				write_dict("blastp1", blastp_hit_dict, str(i+1))

				# reverse blasto tp confirm results from 
				if len(blastp_hit_dict) > 0:
					print("\n\nPerforming reciprocal blastp...")
					blastp_rev_list = recip_blast("blastp", blastp_hit_dict, subj_db, qseq_id, prot_db)
					print("Done")

					# update blastp_hit_dict after validation
					# if 0 seqs were valid, make dict empty
					if len(blastp_rev_list) == 0:
						blastp_hit_dict = {}
					# if len of valid seqs != len of dict, then dict must be updated
					elif len(blastp_rev_list) != len(blastp_hit_dict):
						update_blast_dict(blastp_hit_dict, blastp_rev_list)

					# write updated blast results to txt file
					write_dict("blastp2", blastp_hit_dict, str(i+1))

					# write summary txt file
					write_summary("blastp", blastp_hit_dict, prot_specs, all_specs, str(i+1))

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
				write_dict("tblastn", tblastn_hit_dict, str(i+1))

				# perform blastx on each tblastn result, keep valid results
				if len(tblastn_hit_dict) > 0:
					print("\n\nPerforming reciprocal blastx...")
					blastx_rev_list = recip_blast("blastx", tblastn_hit_dict, subj_db, qseq_id, nucl_db)
					print("Done")

					# update tblastn_hit_dict after validation
					# if 0 seqs were valid, make dict empty
					if len(blastx_rev_list) == 0:
						tblastn_hit_dict = {}
					# if len of valid seqs != len of dict, then dict must be updated
					elif len(blastx_rev_list) != len(tblastn_hit_dict):
						update_blast_dict(tblastn_hit_dict, blastx_rev_list)

					# write blast results to txt file
					write_dict("blastx", tblastn_hit_dict, str(i+1))

					# write summary txt file
					write_summary("tblastn", tblastn_hit_dict, blastp_hit_dict.values(), all_specs, str(i+1))

			# add results from this query round to result dicts
			prot_dict[q_spec_name] = blastp_hit_dict
			nucl_dict[q_spec_name] = tblastn_hit_dict


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
				write_dict("blastx", blastx_hit_dict,str(i+1))

				# reverse blasto tp confirm results from 
				if len(blastx_hit_dict) > 0:
					print("\n\nPerforming reciprocal tblastn...")
					tblastn_rev_list = recip_blast("tblastn", blastx_hit_dict, subj_db, qseq_id, prot_db)
					print("Done")

					# update blastx_hit_dict after validation
					# if 0 seqs were valid, make dict empty
					if len(tblastn_rev_list) == 0:
						blastx_hit_dict = {}
					# if len of valid seqs != len of dict, then dict must be updated
					elif len(tblastn_rev_list) != len(blastx_hit_dict):
						update_blast_dict(blastx_hit_dict, tblastn_rev_list)

					# write updated blast results to txt file
					write_dict("tblastn", blastx_hit_dict, str(i+1))

					# write summary txt file
					write_summary("blastx", blastx_hit_dict, prot_specs, all_specs, str(i+1))

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
				write_dict("tblastx1", tblastx_hit_dict, str(i+1))

				# perform blastx on each tblastn result, keep valid results
				if len(tblastx_hit_dict) > 0:
					print("\n\nPerforming reciprocal tblastx...")
					tblastx_rev_list = recip_blast("tblastx", tblastx_hit_dict, subj_db, qseq_id, nucl_db)
					print("Done")

					# update tblastn_hit_dict after validation
					# if 0 seqs were valid, make dict empty
					if len(tblastx_rev_list) == 0:
						tblastn_hit_dict = {}
					# if len of valid seqs != len of dict, then dict must be updated
					elif len(tblastx_rev_list) != len(tblastx_hit_dict):
						update_blast_dict(tblastx_hit_dict, tblastx_rev_list)

					# write blast results to txt file
					write_dict("tblastx2", tblastx_hit_dict, str(i+1))

					# write summary txt file
					write_summary("tblastx", tblastx_hit_dict, blastx_hit_dict.values(), all_specs, str(i+1))
		
			# add results from this query round to result dicts
			prot_dict[q_spec_name] = blastx_hit_dict
			nucl_dict[q_spec_name] = tblastx_hit_dict
		
		# select next species w/ protein dataset to automate query
		if n - i > 1:
			# make sure there are species to select
			if len(prot_dict.keys()) == 0:
				print("\n\nNo species with protein dataset to select from.\n\n")
				break
			
			processor = CladesProcessor(q_spec_name, prot_dict[list(prot_dict.keys())[-1]], prot_db)
			
			# select a random species in a different family
			next_spec_name = processor.get_rand_spec()
			
			# make sure there is a next species
			if next_spec_name is None:
				print("\n\nNo species in other family ranks to select from.\n\n")
				break
			
			# get next species prot id
			next_id = processor.get_id(next_spec_name)
			
			# get fasta file of the next query protein
			next_id_fasta = processor.get_id_fasta(next_id)
			
			# get fasta file of the next query protein database
			
			# update current query files
			seq_query = next_id_fasta
			ds_query = prot_file_paths[next_spec_name]
			
			# update current species name
			q_spec_name = next_spec_name
			
			# reassign nucl_db
			nucl_fasta_file = "nucl.fna"
			
	
	# ANNOTATION
	print("\n\nPerforming annotation on nucleotide sequences...")
	
	man_anno_files = []
	anno_dict_list = {}
	
	for q_spec_name in nucl_dict:
		# make an annotation object
		annotation = BlastAnnot(nucl_dict[q_spec_name], q_spec_name)
		# get the inital list of seqs with gap (>= 10 aa or multiple alignments), and no gap (anything else)
		no_gap_dict, gap_dict = annotation.process_seqs()

		# find 5' and 5' stop codons, find start codon starting from 5' direction
		annotated_dict, no_start_seqs = annotation.annotate_no_gaps(no_gap_dict, nucl_db)
		anno_dict_list[q_spec_name] = annotated_dict
		
		name = q_spec_name.replace(" ", "_")
		write_dict("annotated", annotated_dict, name)

		# update no_gap_dict and gap_dict if any sequence in no_gap_dict does not have a start codon
		annotation.update_gap_dicts(no_start_seqs, no_gap_dict, gap_dict)

		# output file of nucleotide seqs that need manual annotation
		man_anno_file = annotation.get_man_annot(gap_dict)
		man_anno_files.append(man_anno_file)
		
	outfile = "all_man_anno.txt"  # Specify the name of the output file
	lines = []
	
	for manu_file in man_anno_files:
		with open(manu_file, "r") as file:
			ls = file.readlines()[1:]
			lines.append(ls)

	with open(outfile, "w") as outfile:
		outfile.write("Query_Species\tSubject_ID\tSubject_Species\tE-Value\tSubject_Start\tSubject_End\tQuery_Start\tQuery_End\tQuery_Alignment\tSubject_Alignment\tReading_Frame\n")
		for line in lines:
			outfile.writelines(line)
			
	print("Done")
	
	# CLUSTAL
	print("\n\nPerforming clustal analysis on nucleotide sequences...")
	
	# remove any dupe protein from multi queries
	final_prot_dict = {}
	final_annotated_dict = {}
	
	for query_spec in prot_dict.keys():
		for seq_id in prot_dict[query_spec].keys():
			if seq_id not in final_prot_dict.keys():
				final_prot_dict[seq_id] = prot_dict[query_spec][seq_id]
				
	for query_spec in anno_dict_list.keys():
		for seq_id in anno_dict_list[query_spec].keys():
			if seq_id not in final_annotated_dict.keys():
				final_annotated_dict[seq_id] = anno_dict_list[query_spec][seq_id]
			else:
				curr_len = len(anno_dict_list[query_spec][seq_id][4])
				poss_len = len(final_annotated_dict[seq_id][4])
				if poss_len > curr_len:
					final_annotated_dict[seq_id] = anno_dict_list[query_spec][seq_id]
				
	
	# make a clustal object
	clustal = Clustal(final_annotated_dict, final_prot_dict, prot_db)
	
	# get fasta file for all seqs fron aumotated annotation
	clustal.get_seqs_fasta()
	
	# run clustal with the result fasta file => output result file in clustal and fasta format
	clustal.run_clustal()
	print("Done")
	
	

if __name__ == "__main__":
	main(sys.argv)