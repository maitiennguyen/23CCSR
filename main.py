import sys
import os
import csv
from blast import *
from annotation import *
from clustal import *
from process_specs import *
from user_input import *
from organize_files import *

def main(argv):
	
	# process terminal arguments
	input_processor = InputProcessor(argv)
	input_processor.check_invalid_flag()
	
	# get required parameters
	seq_query = input_processor.get_qseq()
	ds_query = input_processor.get_qdb()
	q_spec_name = input_processor.get_qname()
	q_type = input_processor.get_qtype()
	taxIDS = input_processor.get_sset()
	download = input_processor.get_download()
	
	# get optional parameters
	blastp_evalue = input_processor.get_evalue()
	tblastn_evalue = input_processor.get_evalue()
	tblastx_evalue = input_processor.get_evalue()
	blastx_evalue = input_processor.get_evalue()
	
	# check if provided parameter files exist
	if os.path.exists(seq_query):
		pass
	else:
		print("\nERROR: Query sequence file not found in directory.\n")
		sys.exit()
	if os.path.exists(ds_query):
		pass
	else:
		print("\nERROR: Query nucl/prot dataset file not found in directory.\n")
		sys.exit()
		
	# save a copy of the original input fil names
	seq_query_copy = seq_query
	ds_query_copy = ds_query
	
	
	# BLAST
	
	nucl_fasta_file = None 
	prot_fasta_file = None 
	all_specs = None
	prot_specs = None
	prot_file_paths = None
	
	
	prot_db = None
	nucl_db = None
	prot_dict = {}
	nucl_dict = {}
	i = 1
	q_specs_list = [q_spec_name]
	
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
		try:
			nucl_fasta_file = "nucl.fna"
			prot_fasta_file = "prot.faa"
			all_specs = read_list("all_specs.txt")
			prot_specs = read_list("prot_data_specs.txt")
			prot_file_paths = txt_to_dict("prot_files_all_dict.txt")
		except FileNotFoundError:
			print("\nERROR: Required search set dababase and/or files not found in directory. See README.md to troubleshoot\n")
			sys.exit()

	# the number of times to run queries
	while True:
		
		print("\n\nRunning blast... Round: " + str(i) + ". Query species: " + q_spec_name)
		
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
				write_dict("blastp1", blastp_hit_dict, str(i))

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
					write_dict("blastp2", blastp_hit_dict, str(i))

					# write summary tsv file
					write_summary("blastp", blastp_hit_dict, prot_specs, all_specs, str(i))

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
				write_dict("tblastn", tblastn_hit_dict, str(i))

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
					write_dict("blastx", tblastn_hit_dict, str(i))

					# write summary tsv file
					write_summary("tblastn", tblastn_hit_dict, blastp_hit_dict.values(), all_specs, str(i))
					
			# see if there are any new results
			prot_ids = []
			for spec in prot_dict.keys():
				for seq_id in prot_dict[spec].keys():
					if seq_id not in prot_ids:
						prot_ids.append(seq_id)
					
			nucl_seqs = {}
			for spec in nucl_dict.keys():
				for scaff in nucl_dict[spec].keys():
					if scaff not in nucl_seqs.keys():
						nucl_seqs[scaff] = []
					for value in nucl_dict[spec][scaff]:
						nucl_seqs[scaff].append(value)
					
			curr_prot_rslts = set(blastp_hit_dict.keys())
			all_prot_rslts = set(prot_ids)
			
			no_new_nucl = False
			for seq in tblastn_hit_dict.keys():
				if seq in nucl_seqs.keys():
					for hit in tblastn_hit_dict[seq]:
						for hit2 in nucl_seqs[seq]:
							if hit2[1] < hit2[2]:
								start = hit2[1]
								end = hit2[2]
							else:
								start = hit2[2]
								end = hit2[1]
							if start <= hit[1] <= end or start <= hit[2] <= end:
								no_new_nucl = True
			
			# add results from this query round to result dicts
			prot_dict[q_spec_name] = blastp_hit_dict
			nucl_dict[q_spec_name] = tblastn_hit_dict
								
			# end blast if no new results found					
			if curr_prot_rslts.issubset(all_prot_rslts) and no_new_nucl:
				print("\n\nNo new results were found in this round, moving onto annotation.")
				break
				

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
				write_dict("blastx", blastx_hit_dict,str(i))

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
					write_dict("tblastn", blastx_hit_dict, str(i))

					# write summary tsv file
					write_summary("blastx", blastx_hit_dict, prot_specs, all_specs, str(i))

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
				write_dict("tblastx1", tblastx_hit_dict, str(i))

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
					write_dict("tblastx2", tblastx_hit_dict, str(i))

					# write summary tsv file
					write_summary("tblastx", tblastx_hit_dict, blastx_hit_dict.values(), all_specs, str(i))
					
			# see if there are any new results
			prot_ids = []
			for spec in prot_dict.keys():
				for seq_id in prot_dict[spec].keys():
					if seq_id not in prot_ids:
						prot_ids.append(seq_id)
					
			nucl_seqs = {}
			for spec in nucl_dict.keys():
				for scaff in nucl_dict[spec].keys():
					if scaff not in nucl_seqs.keys():
						nucl_seqs[scaff] = []
					for value in nucl_dict[spec][scaff]:
						nucl_seqs[scaff].append(value)
					
			curr_prot_rslts = set(blastx_hit_dict.keys())
			all_prot_rslts = set(prot_ids)
			
			no_new_nucl = False
			for seq in tblastx_hit_dict.keys():
				if seq in nucl_seqs.keys():
					for hit in tblastx_hit_dict[seq]:
						for hit2 in nucl_seqs[seq]:
							start = None
							end = None
							if hit2[1] < hit2[2]:
								start = hit2[1]
								end = hit2[2]
							else:
								start = hit2[2]
								end = hit2[1]
							if start <= hit[1] <= end or start <= hit[2] <= end:
								no_new_nucl = True
			
			# add results from this query round to result dicts
			prot_dict[q_spec_name] = blastx_hit_dict
			nucl_dict[q_spec_name] = tblastx_hit_dict
			
			# end blast if no new results found
			if curr_prot_rslts.issubset(all_prot_rslts) and no_new_nucl:
				print("\n\nNo new results were found with this run, moving onto annotation.")
				break
		
		# select next species w/ protein dataset to automate query
		# make sure there are species to select
		if len(prot_dict.keys()) == 0:
			print("\n\nNo species in protein results to select from, moving onto annotation.")
			break

		processor = CladesProcessor(prot_dict, prot_db, q_specs_list)

		# select a random species in a different family
		next_spec_name = processor.get_rand_spec()

		# make sure there is a next species
		if next_spec_name is None:
			print("\n\nNo species in other family ranks to select from, moving onto annotation.")
			break
			
		# add next query species to list of query spcies
		q_specs_list.append(next_spec_name)

		# get next species prot id
		next_id = processor.get_id(next_spec_name)

		# get fasta file of the next query protein
		next_id_fasta = processor.get_id_fasta(next_id)

		# update current query files
		seq_query = next_id_fasta
		ds_query = prot_file_paths[next_spec_name]
		
		# update query type if neccessary, (since next query will always be prot)
		if q_type == "nucl":
			q_type = "prot"

		# update current species name
		q_spec_name = next_spec_name

		# reassign nucl_db
		nucl_fasta_file = "nucl.fna"
		
		# increase # of run
		i += 1
	
	# write summary report of all blast results into one file
	final_report = {}
	blast1 = ''
	blast2 = ''
	if q_type == "prot":
		blast1 = "blastp_"
		blast2 = "tblastn_"
	elif q_type == "nucl":
		blast1 = "blastx_"
		blast2 = "tblastx_"
		
	for run in range(i):
		fname1 = blast1 + str(run+1) + "_summary_report.tsv"
		fname2 = blast2 + str(run+1) + "_summary_report.tsv"

		if os.path.exists(fname1):
			with open(fname1, 'r') as file1:
				reader = csv.reader(file1, delimiter="\t")
				for line_num, line in enumerate(reader):
					if line_num == 0:
						continue
                        
					species = line[0]
					
					prot_hit_ids = line[2]
					if prot_hit_ids == 'N/A':
						prot_hit_ids = []
					else:
						prot_hit_ids = prot_hit_ids.split(', ')
					prot_hit_ids = [prot_id.split()[0] for prot_id in prot_hit_ids]

					if species not in final_report.keys():
						final_report[species] = [prot_hit_ids, []]
					else:
						for prot_id in prot_hit_ids:
							if prot_id not in final_report[species][0]:
								final_report[species][0].append(prot_id)

		if os.path.exists(fname2):
			with open(fname2, 'r') as file2:
				reader = csv.reader(file2, delimiter="\t")
				for line_num, line in enumerate(reader):
					if line_num == 0:
						continue
                        
					species = line[0]
					  
					nucl_hit_ids = line[2]
					if nucl_hit_ids == 'N/A':
						nucl_hit_ids = []
					else:
						nucl_hit_ids = nucl_hit_ids.split(', ')
					nucl_hit_ids = [nucl_id.split()[0] for nucl_id in nucl_hit_ids]

					if species not in final_report.keys():
						final_report[species] = [[], nucl_hit_ids]
					else:
						for nucl_id in nucl_hit_ids:
							if nucl_id not in final_report[species][1]:
								final_report[species][1].append(nucl_id)
                                
	# rows to write to tsv file
	blast_sum_rows = []
	# headers for tsv file
	blast_sum_headers = ["Species", "Num Hit(s)", "Protein Hit ID(s)", "Nucleotide Hit ID(s)"]

	for key, value in final_report.items():
		species = key
		count = 0
		prot_hits = value[0]
		nucl_hits = value[1]

		if len(prot_hits) == 0:
			prot_hits = 'N/A'
		else:
			count += len(prot_hits)
			prot_hits = ', '.join(prot_hits)	

		if len(nucl_hits) == 0:
			nucl_hits = 'N/A'
		else:
			count += len(nucl_hits)
			nucl_hits = ', '.join(nucl_hits)

		row = {"Species": species, 
			   "Num Hit(s)": count, 
			   "Protein Hit ID(s)": prot_hits, 
			   "Nucleotide Hit ID(s)": nucl_hits}

		blast_sum_rows.append(row)
        
	# write to tsv file  
	with open("complete_blast_summary_report.tsv", 'w') as blast_sum_file:
		writer = csv.DictWriter(blast_sum_file, delimiter='\t', fieldnames=blast_sum_headers)
		writer.writeheader()
		writer.writerows(blast_sum_rows)
        

	
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
		
	outfile = "all_man_anno.txt" 
	lines = []
	
	for manu_file in man_anno_files:
		with open(manu_file, "r") as file:
			ls = file.readlines()[1:]
			lines.append(ls)

	with open(outfile, "w") as outfile:
		outfile.write("Query_Species\tSubject_ID\tSubject_Species\tE-Value\tSubject_Start\tSubject_End\tQuery_Start\tQuery_End\tQuery_Alignment\tSubject_Alignment\tReading_Frame\n")
		for line in lines:
			outfile.writelines(line)
	
	# write summary report for all manual annotation seqs
	filename = 'all_man_anno.txt' 
	man_sum = "complete_manual_anno_summary.tsv"

	manual_report = {}

	with open(filename, 'r') as file:
		for line_num, line in enumerate(file):
			if line_num == 0:
				continue
			line = line.strip()  
			columns = line.split('\t') 
			query_species = columns[0]
			subject_id = columns[1]
			subject_species = columns[2]
			subject_start = columns[4].split(', ') if columns[4] else []  
			subject_end = columns[5].split(', ') if columns[5] else []  
			posits = subject_start + subject_end
			posits = [int(posit) for posit in posits]
			min_posit = min(posits)
			max_posit = max(posits)

			if subject_species not in manual_report.keys():
				manual_report[subject_species] = [[subject_id], [min_posit], [max_posit], [query_species]]
			else:
				if subject_id not in manual_report[subject_species][0]:
					manual_report[subject_species][0].append(subject_id)
					manual_report[subject_species][1].append(min_posit)
					manual_report[subject_species][2].append(max_posit)
				else:
					index =  manual_report[subject_species][0].index(subject_id)
					new_min = min([manual_report[subject_species][1][index], min_posit])
					new_max = max([manual_report[subject_species][2][index], max_posit])
					manual_report[subject_species][1][index] = new_min
					manual_report[subject_species][2][index] = new_max
					
				if query_species not in  manual_report[subject_species][3]:
					manual_report[subject_species][3].append(query_species)
	
	# rows to write to tsv file
	man_sum_rows = []
	# headers for tsv file
	man_sum_headers = ["Species", "Scaffold ID", "Min Posit", "Max Posit", "Query Species"]

	for key, value in manual_report.items():
		species = key
		scaffold_id = ', '.join(value[0])
		min_posit = ', '.join([str(p) for p in value[1]])
		max_posit = ', '.join([str(p) for p in value[2]])
		q_spec = ', '.join(value[3])
		
		row = {"Species": species, 
			   "Scaffold ID": scaffold_id, 
			   "Min Posit": min_posit, 
			   "Max Posit": max_posit, 
			   "Query Species": q_spec}
		
		man_sum_rows.append(row)
	
	# write to tsv file
	with open(man_sum, 'w') as man_sum_file:
		writer = csv.DictWriter(man_sum_file, delimiter='\t', fieldnames=man_sum_headers)
		writer.writeheader()
		writer.writerows(man_sum_rows)

	print("Done")
	
	
	# CLUSTAL
	print("\n\nPerforming clustal analysis on annotated sequences...")
	
	# remove any dupe protein and scaffold sequences from multi queries
	final_prot_dict = {}
	final_annotated_dict = {}
	nucl_to_add = {}
	nucl_seq_q_spec = {}
	
	for query_spec in prot_dict.keys():
		for seq_id in prot_dict[query_spec].keys():
			if seq_id not in final_prot_dict.keys():
				final_prot_dict[seq_id] = prot_dict[query_spec][seq_id]
						
	for query_spec in anno_dict_list.keys():
		for seq_id in anno_dict_list[query_spec].keys():
			if seq_id not in final_annotated_dict.keys():
				final_annotated_dict[seq_id] = [anno_dict_list[query_spec][seq_id]]
			else:
				for index, value in enumerate(final_annotated_dict[seq_id]):
					curr_start = value[1]
					curr_end = value[2] 
					poss_start = anno_dict_list[query_spec][seq_id][1]
					poss_end = anno_dict_list[query_spec][seq_id][2]

					# check if seq overlaps with curr seq in final dict, if yes, get the longer one
					if curr_start <= poss_start <= curr_end or curr_start <= poss_end <= curr_end or curr_start >= poss_start >= curr_end or curr_start >= poss_end >= curr_end:
						curr_len = len(value[4])
						poss_len = len(anno_dict_list[query_spec][seq_id][4]) 
						if poss_len > curr_len:
							final_annotated_dict[seq_id][index] = anno_dict_list[query_spec][seq_id]

					# if they don't overlap, add it as another seq to clustal
					else:
						if seq_id not in nucl_to_add.keys():
							nucl_to_add[seq_id] = [anno_dict_list[query_spec][seq_id]]
						else:
							if anno_dict_list[query_spec][seq_id] not in nucl_to_add[seq_id]:
								nucl_to_add[seq_id].append(anno_dict_list[query_spec][seq_id])
								
			if seq_id not in nucl_seq_q_spec.keys():
				nucl_seq_q_spec[seq_id] = [query_spec]
			else:
				if query_spec not in nucl_seq_q_spec[seq_id]:
					nucl_seq_q_spec[seq_id].append(query_spec)

								
	for seq_id in nucl_to_add.keys():
		for alignment in nucl_to_add[seq_id]:
			final_annotated_dict[seq_id].append(alignment)
			
	# write summary report for annotated nucl seqs
	auto_sum = "complete_auto_anno_summary.tsv"
	auto_report = {}

	for key, value in final_annotated_dict.items():
		scaffold_id = key
		for hit in value:
			species = hit[0]
			start = hit[1]
			end = hit[2]
			
			if species not in auto_report.keys():
				auto_report[species] = [[scaffold_id], [start], [end], nucl_seq_q_spec[scaffold_id]]
			else:
				auto_report[species][0].append(scaffold_id)
				auto_report[species][1].append(start)
				auto_report[species][2].append(end)
				for spec in nucl_seq_q_spec[scaffold_id]:
					if spec not in auto_report[species][3]:
						auto_report[species][3].append(spec)   
	
	# rows to write to tsv file
	auto_sum_rows = []
	# headers for tsv file
	auto_sum_headers = ["Species", "Scaffold ID", "Start Posit", "End Posit", "Query Species"]
	
	for key, value in auto_report.items():
		species = key
		scaffold_id = ', '.join(value[0])
		min_posit = ', '.join([str(p) for p in value[1]])
		max_posit = ', '.join([str(p) for p in value[2]])
		q_spec = ', '.join(value[3])
		
		row = {"Species": species, 
			   "Scaffold ID": scaffold_id, 
			   "Start Posit": min_posit, 
			   "End Posit": max_posit, 
			   "Query Species": q_spec}
		
		auto_sum_rows.append(row)
	
	# write to tsv file
	with open(auto_sum, 'w') as auto_sum_file:
		writer = csv.DictWriter(auto_sum_file, delimiter='\t', fieldnames=auto_sum_headers)
		writer.writeheader()
		writer.writerows(auto_sum_rows)
		
	
	# make a clustal object
	clustal = Clustal(final_annotated_dict, final_prot_dict, prot_db)
	
	# get fasta file for all seqs fron aumotated annotation
	clustal.get_seqs_fasta()
	
	# run clustal with the result fasta file => output result file in clustal and fasta format
	clustal.run_clustal()
	print("Done")
	
	
	# ORGANIZE FILES INTO FOLDERS
	organizer = Organizer(i, taxIDS, ds_query_copy, seq_query_copy)
	organizer.organize_files()
	

if __name__ == "__main__":
	main(sys.argv)