from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import subprocess

class Clustal():
	
	def __init__(self, anno_nucl_dict, blast_dict, prot_db):
		self.nucl_dict = anno_nucl_dict
		self.rslt_filename = ''
		self.blast_dict = blast_dict
		self.prot_db = prot_db
		
		
	# fasta files for aa and nucl seqs to go into clustal
	def get_seqs_fasta(self):
		file_name = "auto_anno_seqs.fasta"
		
		# clear file contents if already exists
		if os.path.isfile(file_name):
			with open(file_name, "w") as file:
				pass
			
		# add prot seqs to fasta file
		if self.blast_dict is not None:
			for seq_id in self.blast_dict.keys():
				db_info = subprocess.run("blastdbcmd -db {0} -entry {1}".format(self.prot_db, seq_id).split(), capture_output=True, text=True).stdout.split("\n")
				des = '\t'.join(db_info[0].split()[1:])
				aa_seq = Seq(''.join(db_info[1:-1]))

				fasta_seq = SeqIO.SeqRecord(aa_seq, id=seq_id, description=des)

				with open(file_name, "a") as fasta_file:
					SeqIO.write(fasta_seq, fasta_file, "fasta")
		
		# add nucl seqs to fasta file
		for seq_id in self.nucl_dict.keys():
			str_annotated_seqs = [str(item) for item in self.nucl_dict[seq_id][0:-1]]
			des = '\t'.join(str_annotated_seqs)
			
			nucl_seq = Seq(self.nucl_dict[seq_id][4])
			# translate nucl seq
			t_nucl_seq = nucl_seq.translate(table=1, stop_symbol="")
			
			fasta_seq = SeqIO.SeqRecord(t_nucl_seq, id=seq_id, description=des)

			with open(file_name, "a") as fasta_file:
				SeqIO.write(fasta_seq, fasta_file, "fasta")
		
		self.rslt_filename = file_name
		
	
	# add spec name next to alignment in clustal file format 
	def add_spec_names(self, clustal_file):
		
		# combine two dicts
		self.nucl_dict.update(self.blast_dict)
		anno_dicts = self.nucl_dict 
		
		# get clustal results
		with open(clustal_file, 'r') as f:
			clustal_lines = f.readlines()

		new_clustal_lines = []
		for line in clustal_lines:
			new_line = line.rstrip('\n')

			if len(new_line) > 0 and new_line[0] != ' ' and "CLUSTAL" not in new_line:
				seq_id = new_line.split()[0].strip()
				spec_name = ''
				
				if isinstance(anno_dicts[seq_id][0], list):
					spec_name = anno_dicts[seq_id][0][0]
				else:	
					spec_name = anno_dicts[seq_id][0]
					
				if spec_name:
					new_line += f'\t\t{spec_name.capitalize()}'

			new_clustal_lines.append(new_line)
		
		# write clustal results with name
		with open(clustal_file, 'w') as f:
			f.write('\n'.join(new_clustal_lines))
		
	# run clustal
	def run_clustal(self):
		
		if os.path.exists(self.rslt_filename) and os.stat(self.rslt_filename).st_size != 0:
			
			count = 0
			with open(self.rslt_filename, "r") as file:
				for record in SeqIO.parse(file, "fasta"):
					count += 1
					
			if count > 1:
		
				out_name1 = "auto_algn.clustal"
				out_name2 = "auto_algn.fasta"

				subprocess.run("clustalo -i {0} -o {1} --outfmt=clustal --resno --threads=16 --force".format(self.rslt_filename, out_name1).split())
				subprocess.run("clustalo -i {0} -o {1} --outfmt=fasta --resno --threads=16 --force".format(self.rslt_filename, out_name2).split())

				self.add_spec_names(out_name1)
				
			else:
				print("Fasta file contains 1 sequence, nothing to align.")
				
			
		else:
			print("No sequences to feed into clustal.")
		