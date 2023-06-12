from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import subprocess

class Clustal():
	
	def __init__(self, annotated_seqs, blastx_dict):
		self.seqs_dict = annotated_seqs
		self.seq_ids = annotated_seqs.keys()
		self.nucl_filename = ''
		#self.aa_filename = ''
		self.spec_dict = blastx_dict
		
		
	# fasta files for aa and nucl seqs to go into clustal
	def get_seqs_fasta(self):
		file_name1 = "auto_anno_nucl_seqs.fasta"
		#file_name2 = "auto_anno_aa_seqs.fasta"
		
		# clear file contents if already exists
		if os.path.isfile(file_name1):
			with open(file_name1, "w") as file:
				pass	
		# if os.path.isfile(file_name2):
		# 	with open(file_name2, "w") as file:
		# 		pass
		
		for seq_id in self.seq_ids:
			str_annotated_seqs = [str(item) for item in self.seqs_dict[seq_id][0:-1]]
			des = '\t'.join(str_annotated_seqs)
			
			nucl_seq = Seq(self.seqs_dict[seq_id][4])
			
			fasta_nucl_seq = SeqIO.SeqRecord(nucl_seq, id=seq_id, description=des)

			with open(file_name1, "a") as fasta_file:
				SeqIO.write(fasta_nucl_seq, fasta_file, "fasta")
				
# 			aa_seq = nucl_seq.translate(table=1, stop_symbol="")
# 			fasta_aa_seq = SeqIO.SeqRecord(aa_seq, id=seq_id, description=des)
			
# 			with open(file_name2, "a") as fasta_file:
# 				SeqIO.write(fasta_aa_seq, fasta_file, "fasta")
		
		self.nucl_filename = file_name1
		#self.aa_filename = file_name2
		
	
	# add spec name next to alignment in clustal file format 
	def add_spec_names(self, clustal_file):
		with open(clustal_file, 'r') as f:
			clustal_lines = f.readlines()

		new_clustal_lines = []
		for line in clustal_lines:
			new_line = line.rstrip('\n')

			if len(new_line) > 0 and new_line[0] != ' ' and "CLUSTAL" not in new_line:
				seq_id = new_line.split()[0].strip()
				spec_name = self.spec_dict[seq_id][0]
				if spec_name:
					new_line += f'\t{spec_name.capitalize()}'

			new_clustal_lines.append(new_line)

		with open(clustal_file, 'w') as f:
			f.write('\n'.join(new_clustal_lines))
		
	# run clustal
	def run_clustal(self):
		
		out_name1 = self.nucl_filename[:-6] + "_algn.clustal"
		out_name2 = self.nucl_filename[:-6] + "_algn.fasta"
		
		subprocess.run("clustalo -i {0} -o {1} --outfmt=clustal --resno --threads=16 --force".format(self.nucl_filename, out_name1).split())
		subprocess.run("clustalo -i {0} -o {1} --outfmt=fasta --resno --threads=16 --force".format(self.nucl_filename, out_name2).split())
		
# 		out_name3 = self.aa_filename[:-6] + "_algn.clustal"
# 		out_name4 = self.aa_filename[:-6] + "_algn.fasta"
		
# 		subprocess.run("clustalo -i {0} -o {1} --outfmt=clustal --resno --threads=16 --force".format(self.aa_filename, out_name3).split())
# 		subprocess.run("clustalo -i {0} -o {1} --outfmt=fasta --resno --threads=16 --force".format(self.aa_filename, out_name4).split())
		
		self.add_spec_names(out_name1)
		#self.add_spec_names(out_name3)
		