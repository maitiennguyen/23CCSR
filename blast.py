import subprocess
import os
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# get query id based on the provded database
def get_qseq_id(seq, db, typ):
	qseq_id = None
	
	if typ == "prot":
		subprocess.run("blastp -query {0} -db {1} -out {2} -outfmt {3} -num_threads 8".format(seq, db, "qseq_id.blasted", "6").split())
		
	elif typ == "nucl":
		subprocess.run("blastn -query {0} -db {1} -out {2} -outfmt {3} -num_threads 8".format(seq, db, "qseq_id.blasted", "6").split())
	
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
	prot_files = []
	
	# store all species names
	all_specs = []
	
	# species w/ protein dataset
	prot_specs = []
	
	# dict to associate species with their prot file path
	prot_file_paths = {}
	

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
				special_case2 = ["NRRL", "CBS", "JCM", "NYNU", "CRUB", "Ashbya", "UWO(PS)", "MTCC", "UWOPS", "isolate"]
				special_case3 = ["MAG"]
				
				spec_name = spec_name.replace("[", "").replace("]", "").replace("'", "").split()
				
				if any(name.lower() == special_word.lower() for special_word in special_case3 for name in spec_name):
					spec_name = ' '.join(spec_name[2:6])
					
				elif any (spec_name[1].lower() == special_word.lower() for special_word in special_case1):
					if len(spec_name) > 2 and any (spec_name[2].lower() == special_word.lower() for special_word in special_case2):
						spec_name = ' '.join(spec_name[:4])
					else:
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
				prot_files.append(prot_file)
				for name in all_specs:
					if name in acs_asm_tuple[2].replace("[", "").replace("]", ""):
						prot_specs.append(name)
						prot_file_paths[name] = prot_file
				
	if len(nucl_file_paths) > 0:
		
		# comebine all taxID nucleotide databases into one single database
		nucl_dbs = ' '.join(nucl_file_paths)                      
		subprocess.run("cat {0} >> nucl.fna".format(nucl_dbs), shell=True)

		# assign file name
		nucl_fasta_file = "nucl.fna"

	# make sure there are prot data
	if len(prot_files) > 0:

		# comebine all taxID protein databases into one single database
		prot_dbs = ' '.join(prot_files)                      
		subprocess.run("cat {0} >> prot.faa".format(prot_dbs), shell=True)

		# assign file name
		prot_fasta_file = "prot.faa"
		
	# return file names for nucl and prot fasta files
	return nucl_fasta_file, prot_fasta_file, all_specs, prot_specs, prot_file_paths


# use fasta file to make blast database
def get_dbs(fasta_file, seq_type):
	file_name = fasta_file[:-4]
	# make BLAST database
	subprocess.run("makeblastdb -in {0} -out {1} -dbtype {2} -parse_seqids".format(fasta_file, file_name, seq_type).split())
	
	return file_name


# write result of tblastn and tblastx into a file for blastx
def write_seq(blast_dict, key, typ, db):
	file_name = "recip_seq.txt"
	seq = ''
	
	# get seq to perform reciprocal blast
	if typ == "prot":
		db_info = subprocess.run("blastdbcmd -db {0} -entry {1}".format(db, key).split(), capture_output=True, text=True).stdout.split("\n")
		seq = ''.join(db_info[1:-1])
		
	elif typ == "nucl":
		starts = []
		ends = []
		frames = []
		start = None
		end = None
		db_info = None
		
		for info in blast_dict[key]:
			starts.append(info[1])
			ends.append(info[2])
			frames.append(info[9])
		
		if all(int(frame) > 0 for frame in frames):
			start = min(starts)
			end = max(ends)
			db_info = subprocess.run("blastdbcmd -db {0} -entry {1} -strand plus -range {2}-{3}".format(db, key, start, end).split(), capture_output=True, text=True).stdout.split("\n")
			
		elif all(int(frame) < 0 for frame in frames):
			start = max(starts)
			end = min(ends)
			db_info = subprocess.run("blastdbcmd -db {0} -entry {1} -strand minus -range {2}-{3}".format(db, key, end, start).split(), capture_output=True, text=True).stdout.split("\n")
			
		else:
			posits = starts + ends
			posits = [int(p) for p in posits]
			start = min(posits)
			end = max(posits)
			db_info = subprocess.run("blastdbcmd -db {0} -entry {1} -strand plus -range {2}-{3}".format(db, key, start, end).split(), capture_output=True, text=True).stdout.split("\n")
			
		seq = ''.join(db_info[1:-1])
	
	# put key_val into fasta format (seq, id, description)
	blast_seq = SeqIO.SeqRecord(Seq(seq), id=key, description=blast_dict[key][0][0])
	
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


# perform blast against amino acids dataset
def blast_aa_ds(query, typ, db, evalue):
    
	# dict of hits, key=>prot_id value=>[species_name, start, end, evalue, subsequence]
	blast_hits = {}

	# if query is prot seq => perform blastp, if query is nucl => perform blastx
	if typ == "prot":
		blast_typ = "blastp"
	elif typ == "nucl":
		blast_typ = "blastx"

	blast_file_name = blast_typ + "_results.blasted"

	if blast_typ == "blastp":
		# run blastp 
		subprocess.run("blastp -query {0} -db {1} -out {2} -outfmt {3} -evalue {4} -num_threads 8".format(query, db, blast_file_name, "6", evalue).split())
	elif blast_typ == "blastx":
		# run blastx 
		subprocess.run("blastx -query {0} -db {1} -out {2} -outfmt {3} -evalue {4} -num_threads 8".format(query, db, blast_file_name, "6", evalue).split())

	# fill out information into the dictionary
	# open blast result file
	# make sure file is not empty
	if os.stat(blast_file_name).st_size != 0:
		with open(blast_file_name, "r") as blast_rslt:
			for rslt_line in blast_rslt:
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
				db_info = subprocess.run("blastdbcmd -db {0} -entry {1}".format(db, seq_id, sstart, send).split(), capture_output=True, text=True).stdout.split("\n")

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
				special_case2 = ["NRRL", "CBS", "JCM", "NYNU", "CRUB", "Ashbya", "UWO(PS)", "MTCC", "UWOPS", "isolate"]
				special_case3 = ["MAG"]

				if any(name.lower() == special_word.lower() for special_word in special_case3 for name in full_name_list):
					spec_name = ' '.join(full_name_list[2:6])
				elif any (full_name_list[1].lower() == special_word.lower() for special_word in special_case1):
					if len(full_name_list) > 2 and any (full_name_list[2].lower() == special_word.lower() for special_word in special_case2):
						spec_name = ' '.join(full_name_list[:4])
					else:
						spec_name = ' '.join(full_name_list[:3])
				else:
					spec_name = ' '.join(full_name_list[:2])

				# fill out info for dictionary 
				if seq_id not in blast_hits.keys():
					blast_hits[seq_id] = [[spec_name, sstart, send, evalue]]
				else:
					blast_hits[seq_id].append([spec_name, sstart, send, evalue])
												
	# return dict results and file name 
	return blast_hits, blast_file_name


# perform blast against nucleotides dataset
def blast_nucl_ds(query, typ, db, evalue):

	# dict of hits, key=>prot_id value=>[species_name, start, end, evalue, subsequence]
	blast_hits = {}

	# if query is prot seq => perform blastp, if query is nucl => perform blastx
	if typ == "prot":
		blast_typ = "tblastn"
	elif typ == "nucl":
		blast_typ = "tblastx"

	blast_file_name = blast_typ + "_results.blasted"

	if blast_typ == "tblastn":
		# run tblastn 
		subprocess.run("tblastn -query {0} -db {1} -out {2} -evalue {3} -num_threads 8".format(query, db, blast_file_name, evalue).split() + ["-outfmt", "6 qseqid sseqid length qstart qend sstart send evalue qseq sseq sframe"])
	elif blast_typ == "tblastx":
		# run tblastx 
		subprocess.run("tblastx -query {0} -db {1} -out {2} -evalue {3} -num_threads 8".format(query, db, blast_file_name, evalue).split() + ["-outfmt", "6 qseqid sseqid length qstart qend sstart send evalue qseq sseq sframe"])
	
    # fill out information into the dictionary
	# open blast result file
	# make sure file is not empty
	if os.stat(blast_file_name).st_size != 0:
		with open(blast_file_name, "r") as blast_rslt:
			for rslt_line in blast_rslt:
				col = rslt_line.strip().split("\t")
				# get nucleotide id
				seq_id = col[1]
				# length of alignment
				length = col[2]
				# get start of alignment in query
				qstart = col[3]
				# get end of alignment in query
				qend = col[4]
				# get start of alignment in nucl subject
				sstart = col[5]
				# get end of alignment in nucl subject
				send = col[6]
				# get e-value
				evalue = col[7]
				# get query aa seq alignment
				qseq = col[8]
				# get subject aa seq alignment
				sseq = col[9]
				# get subject frame
				sframe = col[10]
				
				# initialize string to capture info
				db_info = None

				# if plus strand
				if int(sstart) < int(send):
					# run blastdbcmd to access database to get subsequence and species name
					db_info = subprocess.run("blastdbcmd -db {0} -entry {1}".format(db, seq_id, sstart, send).split(), capture_output=True, text=True).stdout.split("\n")

				# if mi nus strand
				else:
					# run blastdbcmd to access database to get subsequence and species name
					db_info = subprocess.run("blastdbcmd -db {0} -entry {1}".format(db, seq_id, send, sstart).split(), capture_output=True, text=True).stdout.split("\n")

				# get species name
				spec_name = None
				full_name_list = None
				full_name = db_info[0][1:]

				if "[" in full_name:
					full_name = full_name.replace('[', '').replace(']', '')
				if "'" in full_name:
					full_name = full_name.replace("'", '')

				special_case1 = ["sp.", "aff."] 
				special_case2 = ["NRRL", "CBS", "JCM", "NYNU", "CRUB", "Ashbya", "UWO(PS)", "MTCC", "UWOPS", "isolate"]
				special_case3 = ["MAG"]

				full_name_list = full_name.split()

				if any(name.lower() == special_word.lower() for special_word in special_case3 for name in full_name_list):
					spec_name = ' '.join(full_name_list[3:7])
				elif any(full_name_list[2].lower() == special_word.lower() for special_word in special_case1):
					if any(full_name_list[3].lower() == special_word.lower() for special_word in special_case2):
						spec_name = ' '.join(full_name_list[1:5])
					else:
						spec_name = ' '.join(full_name_list[1:4])
				else:
					spec_name = ' '.join(full_name_list[1:3])

				# fill out info for dictionary 
				seq_id = seq_id.split("|")[1]
				if seq_id not in blast_hits.keys():
					blast_hits[seq_id] = [[spec_name, sstart, send, evalue, qstart, qend, length, qseq, sseq, sframe]]
				else:
					blast_hits[seq_id].append([spec_name, sstart, send, evalue, qstart, qend, length, qseq, sseq, sframe])
						
	# return results dict and file name
	return blast_hits, blast_file_name

# reverse blast
def recip_blast(blast_type, query_dict, db, id_of_interest, db2):	
	blast_file_name = blast_type + "_rev_results.blasted"

	# list of valid hits id
	blastp_hits = []

	# for every seq result from blastp
	for seq_id in query_dict.keys():

		# run reciprical blast with result from first blast as query against protein of interesst species protein dataset
		if blast_type == "blastp":
			# write seq to txt file
			query = write_seq(query_dict, seq_id, "prot", db2)
			subprocess.run("blastp -query {0} -db {1} -out {2} -outfmt {3} -num_threads 8".format(query, db, blast_file_name, "6").split())
		elif blast_type == "blastx":
			# write seq to txt file
			query = write_seq(query_dict, seq_id, "nucl", db2)
			subprocess.run("blastx -query {0} -db {1} -out {2} -outfmt {3} -num_threads 8".format(query, db, blast_file_name, "6").split())
		elif blast_type == "tblastn":
			# write seq to txt file
			query = write_seq(query_dict, seq_id, "prot", db2)
			subprocess.run("tblastn -query {0} -db {1} -out {2} -outfmt {3} -num_threads 8".format(query, db, blast_file_name, "6").split())
		elif blast_type == "tblastx":
			# write seq to txt file
			query = write_seq(query_dict, seq_id, "nucl", db2)
			subprocess.run("tblastx -query {0} -db {1} -out {2} -outfmt {3} -num_threads 8".format(query, db, blast_file_name, "6").split())

		# make sure file is not empty
		if os.stat(blast_file_name).st_size != 0:
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
		spec_name = blast_dict[result][0][0].lower()
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
				
	if os.stat(out_file).st_size == 0:
		out_file = None

	# return file name
	return out_file


# write blast results and and their info to txt file
def write_dict(blast_type, blast_dict, i):
	with open(blast_type + "_" + i + "_dict.txt", "w") as file:
		for key, value in blast_dict.items():
			file.write(f"{key}: {value}\n")
			
			
# write summary for blastp and/or tblastn
def write_summary(blast_type, blast_dict, special_case_list, all_specs_list, i):
		
	with open(blast_type + "_" + i + "_summary_report.txt", "w") as file:
		column_widths = [50, 30, 100]
		
		# write hearder
		file.write("{:<{}} {:<{}} {:<{}}\n".format("Species", column_widths[0], "Num Hit(s)", column_widths[1], "Hit ID(s)", column_widths[2]))
		
		# keep track of species that has been seen
		added_spec = {}
		
		# add all species adn their seq id in blast result into dict, add how mnay hits per id
		for key, value in blast_dict.items():
			spec_name = value[0][0]
			if spec_name not in added_spec.keys():
				added_spec[spec_name] = [(key, len(value))]
			else:
				added_spec[spec_name].append((key, len(value)))
		
		# check of species has either been found in protein blast round or does not have protein dataset
		for spec in all_specs_list:

			if blast_type == "blastp" or blast_type == "blastx":
				if spec not in added_spec.keys() and spec not in special_case_list:
					added_spec[spec] = 0
				elif spec not in added_spec.keys():
					added_spec[spec] = []
			elif blast_type == "tblastn" or blast_type == "tblastx":
				if len(special_case_list) > 0 and any(spec == info[0][0] for info in special_case_list) and spec not in added_spec.keys():
					added_spec[spec] = 0
				elif spec not in added_spec.keys():
					added_spec[spec] = []

		# fill in info to write to file
		for key, value in added_spec.items():
			species = key
			count = None
			
			if blast_type == "blastp" or blast_type == "blastx":
				count = "no protein dataset"
			elif blast_type == "tblastn" or blast_type == "tblastx":
				count = "blastp"
				
			hit_ids = "N/A"
			
			if value != 0:
				count = str(len(value))
				if len(value) != 0:
					hit_ids = ', '.join(f'{seq_id} ({num_hit})' for seq_id, num_hit in value)

			file.write("{:<{}} {:<{}} {:<{}}\n".format(species, column_widths[0], count, column_widths[1], hit_ids, column_widths[2]))
	
# write list to txt file
def write_list(filename, l):
	with open(filename + ".txt", "w") as file:
		for line in l:
			file.write(line + "\n")

# read txt file and return a list of lines from file
def read_list(filename):
	l = []
	with open (filename, "r") as file:
		for line in file:
			l.append(line.strip("\n"))
	return l

# read txt and return a dict
def txt_to_dict(filename):
	result_dict = {}

	with open(filename, 'r') as file:
		for line in file:
			line = line.strip()
			if line:
				parts = line.split(': ')
			
				key = parts[0]
				
				s = "'" + parts[1] + "'"
				values = eval(s)

				result_dict[key] = values

	return result_dict
	

