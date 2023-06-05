import subprocess
import os
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# get query id based on the provded database
def get_qseq_id(prot_seq, db):
	
	qseq_id = None
	
	subprocess.run("blastp -query {0} -db {1} -out {2} -outfmt {3} -num_threads 8".format(prot_seq, db, "qseq_id.blasted", "6").split())
	
	with open("qseq_id.blasted", "r") as blastp_rslt:
		qseq_id = blastp_rslt.readline().split("\t")[1]
		
	return qseq_id	

# returns filea with all protein seqs & nucleotide seqs in fasta format, return None if no data exists
def get_fasta_files(taxID_list):
	
	# initialize files' name
	nucl_fasta_file = None
	prot_fasta_file = None

	# store nucleotide files' path
	nucl_file_paths = []

	# store protein files path
	prot_file_paths = []
	
	# store all species names
	all_specs = []

	# download genome and proteome and save fasta files
	# for every tax ID 
	for taxID in taxID_list:

		# download reference genome/set of reference genomes and protein, combine them into taxID.zip file
		subprocess.run("datasets download genome taxon {0} --reference --include genome,protein --filename {0}.zip".format(taxID).split())

		# unzip
		subprocess.run("unzip {0}.zip -d {0}".format(taxID).split())

		# delete zip files
		subprocess.run("rm {0}.zip".format(taxID).split())

		# get path to access assembly data report
		path = os.getcwd()+"/{0}/ncbi_dataset/data".format(taxID)

		# reformat assembly data report using dataformat, append each genome's assembly accession number and assembly info to list as tuple, since its fasta file is named based on them
		acs_asm_list_temp = subprocess.run("dataformat tsv genome --inputfile {0}/assembly_data_report.jsonl --fields accession,assminfo-name,organism-name".format(path).split(), \
										capture_output=True, \
										text=True) \
										.stdout \
										.split("\n")
		
		# remove unecessary info (headers and '') from list
		acs_asm_list_temp.pop(0)
		acs_asm_list_temp.pop(-1)

		# parse into list => [assemble accession num, assembly name, species name] 
		acs_asm_list = []
		for acs_asm in acs_asm_list_temp:
			acs_asm_tuple = acs_asm.split("\t")
			
			# exclude hybrid species
			spec_name = acs_asm_tuple[2]
			if " x " not in spec_name:
				acs_asm_list.append(acs_asm_tuple)
				
				# get all valid species names
				special_case1 = ["sp.", "aff."] 
				special_case2 = ["NRRL", "CBS", "JCM", "NYNU", "CRUB", "Ashbya", "UWO(PS)", "MTCC", "UWOPS"]
				special_case3 = ["MAG"]
				
				spec_name = spec_name.replace("[", "").replace("]", "").replace("'", "").split()
				
				if any(name.lower() == special_word.lower() for special_word in special_case3 for name in spec_name):
					spec_name = ' '.join(spec_name[2:6])
				elif len(spec_name) > 2 and any (spec_name[2].lower() == special_word.lower() for special_word in special_case2):
					spec_name = ' '.join(spec_name[:4])
				elif any (spec_name[1].lower() == special_word.lower() for special_word in special_case1):
					spec_name = ' '.join(spec_name[:3])
				else:
					spec_name = ' '.join(spec_name[:2])

				if spec_name not in all_specs:
					all_specs.append(spec_name)

		# make a list of all nucleotide file paths an a list of all protein file paths
		for acs_asm_tuple in acs_asm_list:
			nucl_file = path+"/"+acs_asm_tuple[0]+"/"+acs_asm_tuple[0]+"_"+acs_asm_tuple[1]+"_genomic.fna"
			prot_file = path+"/"+acs_asm_tuple[0]+"/"+"protein.faa"
			# check if file exists
			if os.path.exists(nucl_file):
				nucl_file_paths.append(nucl_file)
			if os.path.exists(prot_file):
				prot_file_paths.append(prot_file)
				
	if len(nucl_file_paths) > 0:
		
		# comebine all taxID nucleotide databases into one single database
		nucl_dbs = ' '.join(nucl_file_paths)                      
		subprocess.run("cat {0} >> nucl.fna".format(nucl_dbs), shell=True)

		# assign file name
		nucl_fasta_file = "nucl.fna"

	# make sure there are prot data
	if len(prot_file_paths) > 0:

		# comebine all taxID protein databases into one single database
		prot_dbs = ' '.join(prot_file_paths)                      
		subprocess.run("cat {0} >> prot.faa".format(prot_dbs), shell=True)

		# assign file name
		prot_fasta_file = "prot.faa"

	# FOR FUTURE USAGE IF DUPE ERROR OCCURS
	# remove posssible duplicates in protein fasta file
	# https://stackoverflow.com/questions/66462611/remove-duplicated-sequences-in-fasta-with-python
	# seen = set()
	# records = []
	# for record in SeqIO.parse("prot_db.faa", "fasta"):  
	# 	if record.seq not in seen:
	# 		seen.add(record.seq)
	# 		records.append(record)
	# subprocess.run("rm prot_db.faa".split())
	# SeqIO.write(records, "prot_db.faa", "fasta")
	

	# delete unecessary files and folders
	# for taxID in taxID_list:
	# 	subprocess.run("rm -r {0}".format(taxID).split())

	# return file names for nucl and prot fasta files
	return nucl_fasta_file, prot_fasta_file, all_specs


# use fasta file to make blast database
def get_dbs(fasta_file, seq_type):
	file_name = fasta_file[:-4]
	# make BLAST database
	subprocess.run("makeblastdb -in {0} -out {1} -dbtype {2} -parse_seqids".format(fasta_file, file_name, seq_type).split())
	
	return file_name

# write result of tblastn into a file for blastx
def write_seq(blast_dict, key):
	file_name = "blast_rslt_seq.txt"

	# put key_val into fasta format (seq, id, description)
	blast_seq = SeqIO.SeqRecord(Seq(blast_dict[key][4]), id=key, description=blast_dict[key][0])
	
	# write to file
	with open(file_name, "w") as fasta_file:
		SeqIO.write(blast_seq, fasta_file, "fasta")
	
	# return file name
	return file_name

# remove any seq that were not validated in blastp_rev
def update_blast_dict(blast_hit_dict, valid_seq_list):
	
	# add any seq that were not validated in blastp_rev to list
	to_rm = [seq_id for seq_id in blast_hit_dict.keys() if seq_id not in valid_seq_list]

	# remove seqs in list from dictionary
	for seq_id in to_rm:
		del blast_hit_dict[seq_id]

# blastp
def blastp(query, db, evalue):
	blast_file_name = "blastp_results.blasted"
	
	# dict of hits, key=>prot_id value=>[species_name, start, end, evalue, subsequence]
	blastp_hits = {}
	
	# run blastp 
	subprocess.run("blastp -query {0} -db {1} -out {2} -outfmt {3} -evalue {4} -num_threads 8".format(query, db, blast_file_name, "6", evalue).split())
	
	# fill out information into the dictionary
	# open blastp result file
	with open(blast_file_name, "r") as blastp_rslt:
		for rslt_line in blastp_rslt:
			col = rslt_line.strip().split("\t")
			# get protein id
			seq_id = col[1]
			# get start of alignment in protein subject
			sstart = col[8]
			# get end of alignment in protein subject
			send = col[9]
			# get e-value
			evalue = col[10]
			
			# run blastdbcmd to access database to get subsequence and species name
			db_info = subprocess.run("blastdbcmd -db {0} -entry {1} -range {2}-{3}".format(db, seq_id, sstart, send).split(), capture_output=True, text=True).stdout.split("\n")
			
			# get species name
			spec_name = None
			full_name_list = None

			if '[[' in db_info[0]:
				full_name_list = re.search(r'\[\[(.+?)\]\s(.+?)\]', db_info[0])
				full_name_list = full_name_list.group(1)+ " "+full_name_list.group(2)
				full_name_list = full_name_list.split()
			else:
				full_name_list = re.search(r'\[(.*?)\]', db_info[0]).group(1).split()
			
			special_case1 = ["sp.", "aff."] 
			special_case2 = ["NRRL", "CBS", "JCM", "NYNU", "CRUB", "Ashbya", "UWO(PS)", "MTCC", "UWOPS"]
			special_case3 = ["MAG"]
			
			if any(name.lower() == special_word.lower() for special_word in special_case3 for name in full_name_list):
				spec_name = ' '.join(full_name_list[2:6])
			elif len(full_name_list) > 2 and any (full_name_list[2].lower() == special_word.lower() for special_word in special_case2):
				spec_name = ' '.join(full_name_list[:4])
			elif any (full_name_list[1].lower() == special_word.lower() for special_word in special_case1):
				spec_name = ' '.join(full_name_list[:3])
			else:
				spec_name = ' '.join(full_name_list[:2])
			
			# get alignment subsequence
			align_seq = ''.join(db_info[1:-1])
			
			# fill out info for dictionary 
			blastp_hits[seq_id] = [spec_name, sstart, send, evalue, align_seq]
												
	# return dict results and file name 
	return blastp_hits, blast_file_name

# blastp, reverse blast
def blastp_rev(query_dict, db, id_of_interest):	
	blast_file_name = "blastp_rev_results.blasted"
	
	# list of valid hits id
	blastp_hits = []
	
	# for every seq result from blastp
	for seq_id in query_dict.keys():
		
		# write seq to txt file
		query = write_seq(query_dict, seq_id)
		
		# run blast p with xt file as query against protein of interesst species protein dataset
		subprocess.run("blastp -query {0} -db {1} -out {2} -outfmt {3} -num_threads 8".format(query, db, blast_file_name, "6").split())
		
		with open(blast_file_name, "r") as blastp_rslt:
			rslt = blastp_rslt.readline().strip('\n').split("\t")
			qseq_id = rslt[0]
			sseq_id = rslt[1]
			
			# check if the first result is id_of_interest
			if sseq_id == id_of_interest:
				blastp_hits.append(qseq_id)
				
	return blastp_hits


# remove seqs from fasta
def rm_seqs_fasta(blast_dict, fasta_file):
	out_file = "nucl_2.fna"
	
	blastp_spec = []
	# get species from blastp
	for result in blast_dict.keys():
		spec_name = blast_dict[result][0].lower()
		if spec_name not in blastp_spec:
			blastp_spec.append(spec_name)
			
	# check if species name is present in the description, write sequences w/o species name
	with open(out_file, "w") as file:
		for seq in SeqIO.parse(fasta_file, "fasta"):
			seq_des = seq.description.lower()
			if "[" in seq_des:
				seq_des = seq_des.replace('[', '').replace(']', '')
			if "'" in seq_des:
				seq_des = seq_des.replace("'", '')
			if not any(spec_name in seq_des for spec_name in blastp_spec):
				SeqIO.write(seq, file, "fasta")

	# return file name
	return out_file


# tblastn
def tblastn(query, db, evalue):
	blast_file_name = "tblastn_results.blasted"
	
	# dict of hits, key=>nucl_id value=>[species_name, start, end, evalue, subsequence]
	tblastn_hits = {}
	
	# run tblastn 
	subprocess.run("tblastn -query {0} -db {1} -out {2} -evalue {3} -num_threads 8".format(query, db, blast_file_name, evalue).split() + ["-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"])
	
	# fill out information into the dictionary
	# open tblastn result file
	with open(blast_file_name, "r") as tblastn_rslt:
		for rslt_line in tblastn_rslt:
			col = rslt_line.strip().split("\t")
			# get nucleotide id
			seq_id = col[1]
			# get start of alignment in nucl subject
			sstart = col[8]
			# get end of alignment in nucl subject
			send = col[9]
			# get e-value
			evalue = col[10]
			# get aa seq
			sseq = col[12]
			
			# initialize string to capture info
			db_info = None
			
			# if plus strand
			if int(sstart) < int(send):
				# run blastdbcmd to access database to get subsequence and species name
				db_info = subprocess.run("blastdbcmd -db {0} -entry {1} -range {2}-{3} -strand plus".format(db, seq_id, sstart, send).split(), capture_output=True, text=True).stdout.split("\n")
			
			# if mi nus strand
			else:
				# run blastdbcmd to access database to get subsequence and species name
				db_info = subprocess.run("blastdbcmd -db {0} -entry {1} -range {2}-{3} -strand minus".format(db, seq_id, send, sstart).split(), capture_output=True, text=True).stdout.split("\n")
			
			# get species name
			spec_name = None
			full_name_list = None
			full_name = db_info[0][1:]
			
			if "[" in full_name:
				full_name = full_name.replace('[', '').replace(']', '')
			if "'" in full_name:
				full_name = full_name.replace("'", '')
			
			special_case1 = ["sp.", "aff."] 
			special_case2 = ["NRRL", "CBS", "JCM", "NYNU", "CRUB", "Ashbya", "UWO(PS)", "MTCC", "UWOPS"]
			special_case3 = ["MAG"]
			
			full_name_list = full_name.split()

			if any(name.lower() == special_word.lower() for special_word in special_case3 for name in full_name_list):
				spec_name = ' '.join(full_name_list[3:7])
			elif any(full_name_list[3].lower() == special_word.lower() for special_word in special_case2):
				spec_name = ' '.join(full_name_list[1:5])
			elif any(full_name_list[2].lower() == special_word.lower() for special_word in special_case1):
				spec_name = ' '.join(full_name_list[1:4])
			else:
				spec_name = ' '.join(full_name_list[1:3])
			
			# get alignment subsequence
			align_seq = ''.join(db_info[1:-1])

			# fill out info for dictionary 
			tblastn_hits[seq_id.split("|")[1]] = [spec_name, sstart, send, evalue, align_seq, sseq]
						
	# return results dict and file name
	return tblastn_hits, blast_file_name


# blastx
def blastx_rev(query_dict, db, id_of_interest):
	blast_file_name = "blastx_rev_results.blasted"
	
	# list of valid hits id
	blastx_hits = []
	
	# dict of the num of hits to check for posssible genome dup
	num_hits_dict = {}
	
	for seq_id in query_dict.keys():
	
		query = write_seq(query_dict, seq_id)
		
		subprocess.run("blastx -query {0} -db {1} -out {2} -outfmt {3} -num_threads 8".format(query, db, blast_file_name, "6").split())
		
		with open(blast_file_name, "r") as blastx_rslt:
			rslt = blastx_rslt.readline().strip('\n').split("\t")
			qseq_id = rslt[0]
			sseq_id = rslt[1]
			if sseq_id == id_of_interest:
				blastx_hits.append(qseq_id)
				
	return blastx_hits
						

def main(argv):
	cse4_prot = "cse4_prot.fasta"
	saccer_aa = "sacc_cere.faa"
	blastp_evalue = "1e-01"  
	tblastn_evalue = "1e-01" 
	taxIDS = ["147537"]
	
	# make amino acid database for blastx subject
	subj_db = get_dbs(saccer_aa, "prot")
	
	# perform blastp using provided prot seq and provided aa db to confirm qseq_id in aa db
	qseq_id = get_qseq_id(cse4_prot, subj_db)
	
	# get user tax id input
	taxID_list = taxIDS
	
	# get protein and nucleotide fasta files
	nucl_fasta_file, prot_fasta_file, all_specs = get_fasta_files(taxID_list)

	# prot_fasta_file = "prot.faa"
	# subj_db = "sacc_cere"
	# qseq_id = "NP_012875.2"
	# nucl_fasta_file = "nucl_2.fna"
	# nucl_db = "nucl"
	# prot_db = "prot"
	
	# check if there is a protein fasta file
	if prot_fasta_file is not None:
		
		# get protein blast database
		prot_db = get_dbs(prot_fasta_file, "prot")
	
		#perform blastp of the query sequence against protein database, get dict of seq id and its info
		blastp_hit_dict, blastp_rslt_file = blastp(cse4_prot, prot_db, blastp_evalue)
		
		with open("blastp_hit_dict1.txt", "w") as file:
	 		for key, value in blastp_hit_dict.items():
	 			file.write(f"{key}: {value}\n")
		
		# reverse blasto tp confirm results from 
		if len(blastp_hit_dict) > 0:
			blastp_rev_list = blastp_rev(blastp_hit_dict, subj_db, qseq_id)
		
			# update blastp_hit_dict after validation
			# if 0 seqs were valid, make dict empty
			if len(blastp_rev_list) == 0:
				blastp_hit_dict = {}
			# if len of valid seqs != len of dict, then dict must be updated
			elif len(blastp_rev_list) != len(blastp_hit_dict):
				update_blast_dict(blastp_hit_dict, blastp_rev_list)

			with open("blastp_hit_dict2.txt", "w") as file:
				for key, value in blastp_hit_dict.items():
					file.write(f"{key}: {value}\n")
					
			with open("blastp_SUMMARY.txt", "w") as file:
				column_widths = [50, 20, 100]
				
				file.write("{:<{}} {:<{}} {:<{}}\n".format("Species", column_widths[0], "Num Hits", column_widths[1], "Hit ID(s)", column_widths[2]))
				
				added_spec = {}
				
				for key, value in blastp_hit_dict.items():
					if value[0] not in added_spec.keys():
						added_spec[value[0]] = [key]
					else:
						added_spec[value[0]].append(key)
						
				for spec in all_specs:
					if spec not in added_spec.keys():
						added_spec[spec] = []
				
				for key, value in added_spec.items():
					species = key
					count = str(len(value))
					hit_ids = ', '.join(value)
					
					file.write("{:<{}} {:<{}} {:<{}}\n".format(species, column_widths[0], count, column_widths[1], hit_ids, column_widths[2]))
					
			# check if there is a nucleotide fasta file
			if nucl_fasta_file is not None:
				# remove blastp hit species from nucleotide fasta file, update nucl_fasta_file
				nucl_fasta_file = rm_seqs_fasta(blastp_hit_dict, nucl_fasta_file)

	
	# check if there is a nucleotide fasta file
	if nucl_fasta_file is not None:
		
		# get nucleotide blast database
		nucl_db = get_dbs(nucl_fasta_file, "nucl")
	
		# perfrom tblasn of query sequence against nucleotide database, get dict of seq id and its info
		tblastn_hit_dict, tblastn_rslt_file = tblastn(cse4_prot, nucl_db, tblastn_evalue)
		
		with open("tblastn_hit_dict1.txt", "w") as file:
			for key, value in tblastn_hit_dict.items():
				file.write(f"{key}: {value}\n")
		
		# perform blastx on each tblastn result, keep valid results
		if len(tblastn_hit_dict) > 0:
			blastx_rev_list = blastx_rev(tblastn_hit_dict, subj_db, qseq_id)
		
			# update tblastn_hit_dict after validation
			# if 0 seqs were valid, make dict empty
			if len(blastx_rev_list) == 0:
				tblastn_hit_dict = {}
			# if len of valid seqs != len of dict, then dict must be updated
			elif len(blastx_rev_list) != len(tblastn_hit_dict):
				update_blast_dict(tblastn_hit_dict, blastx_rev_list)

			with open("tblastn_hit_dict2.txt", "w") as file:
				for key, value in tblastn_hit_dict.items():
					file.write(f"{key}: {value}\n")
					
			with open("tblastn_SUMMARY.txt", "w") as file:
				column_widths = [50, 20, 100]
				
				file.write("{:<{}} {:<{}} {:<{}}\n".format("Species", column_widths[0], "Num Hits", column_widths[1], "Hit ID(s)", column_widths[2]))
    
				added_spec = {}
				
				for key, value in tblastn_hit_dict.items():
					if value[0] not in added_spec.keys():
						added_spec[value[0]] = [key]
					else:
						added_spec[value[0]].append(key)	
				for spec in all_specs:
					if any(spec == info[0] for info in blastp_hit_dict.values()) and spec not in added_spec.keys():
						added_spec[spec] = "N/A"
					elif spec not in added_spec.keys():
						added_spec[spec] = []
				
				for key, value in added_spec.items():
					species = key
					count = "N/A"
					hit_ids = ''
					if value != "N/A":
						count = str(len(value))
						hit_ids = ', '.join(value)

					file.write("{:<{}} {:<{}} {:<{}}\n".format(species, column_widths[0], count, column_widths[1], hit_ids, column_widths[2]))

		


if __name__ == "__main__":
	main(sys.argv)

