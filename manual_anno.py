import subprocess
import sys
import os
parent_directory_path = os.path.dirname(os.getcwd())
sys.path.append(parent_directory_path)
from blast import *

# convert all_man_anno.txt into dictionary of dictionaries, outer key = query sepcies, inner key = scaffld id, inner value = species name, s_starts, s_stops, reading frame, sequence
man_anno_file = "all_man_anno.txt"
man_anno_dict = {}

with open(man_anno_file, 'r') as file:
	for line_num, line in enumerate(file):
		if line_num == 0:
			continue
			
		line = line.strip()
		columns = line.split('\t') 
		q_spec = columns[0]
		s_id = columns[1]
		s_spec = columns[2]
		s_starts = [int(num) for num in columns[4].split(', ')]
		s_stops = [int(num) for num in columns[5].split(', ')]
		frames = [int(num) for num in columns[10].split(', ')]
		q_starts = [int(num) for num in columns[6].split(', ')]
		q_stops = [int(num) for num in columns[7].split(', ')]
		
		
		if q_spec not in man_anno_dict.keys():
			man_anno_dict[q_spec] = {s_id : [s_spec, s_starts, s_stops, frames, q_starts, q_stops]}
		else:
			man_anno_dict[q_spec][s_id] = [s_spec, s_starts, s_stops, frames,  q_starts, q_stops]

# get dict of species and thei prot dataset files
prot_file_paths = txt_to_dict("prot_files_all_dict.txt")

# q_spec and their query prot id
q_prot_ids = {"Saccharomyces cerevisiae" : "NP_010615.3", "Candida verbasci" : "CAI5756721.1"}
			
class MultiAlignAnno():
	
	def __init__(self, all_man_anno_dict, nucl_db, prot_file_paths, q_prot_ids):
		self.man_anno_dict = all_man_anno_dict
		self.nucl_db = nucl_db
		self.prot_file_paths = prot_file_paths
		self.q_prot_ids = q_prot_ids
		self.to_remove = {}
		self.to_further_anno = {}
		self.to_keep = {}
		
		
	def add_to_remove(self, q_spec, s_id):
		if q_spec not in self.to_remove.keys():
			self.to_remove[q_spec] = [s_id]
		else:
			self.to_remove[q_spec].append(s_id)
		
		
	def add_to_further_anno(self, q_spec, s_id):
		if q_spec not in self.to_further_anno.keys():
			self.to_further_anno[q_spec] = [s_id]
		else:
			self.to_further_anno[q_spec].append(s_id)
			
	def add_to_keep(self, q_spec, s_id, index):
		if q_spec not in self.to_keep.keys():
			self.to_keep[q_spec] = {s_id : index}
		else:
			self.to_keep[q_spec][s_id] = index
			
			
	def get_nucl_seq(self, seq_id, mode, start, stop):
		db_info = subprocess.run("blastdbcmd -db {0} -entry {1} -strand {2} -range {3}-{4}".format(self.nucl_db, seq_id, mode, start, stop).split(), capture_output=True, text=True).stdout.split("\n")
		return ''.join(db_info[1:-1])

	
	def alignments_overlap(self, starts, stops):
		if starts[0] <= starts[1] <= stops[0] or starts[0] <= stops[1] <= stops[0]:
			return True
		else:
			return False
	
	
	def write_query_fasta(self, q_spec, s_id, start, stop):
		seq = ''
		
		if start < stop:
			seq = self.get_nucl_seq(s_id, "plus", start, stop)
		else:
			seq = self.get_nucl_seq(s_id, "minus", stop, start)
			
		blast_seq = SeqIO.SeqRecord(Seq(seq), id=s_id, description=self.man_anno_dict[q_spec][s_id][0])

		file_name = "man_recip_q.fasta"
		with open(file_name, "w") as fasta_file:
			SeqIO.write(blast_seq, fasta_file, "fasta")

		return file_name
		
		
	def blastx(self, s_id, q_spec, q_prot_db):
		align_info = self.man_anno_dict[q_spec][s_id]
		blast_file_name = "man_recip_blastx.blasted"
		valid_hits = []
		
		for i in range(len(align_info[1])):	
			start = None
			stop = None
			if align_info[3][i] > 0:
				start = align_info[1][i]
				stop = align_info[2][i]
			elif align_info[3][i] < 0:
				start = align_info[2][i]
				stop = align_info[1][i]
			query = self.write_query_fasta(q_spec, s_id, start, stop)
			subprocess.run("blastx -query {0} -db {1} -out {2} -outfmt {3} -num_threads 8".format(query, q_prot_db, blast_file_name, "6").split())
			
			if os.stat(blast_file_name).st_size != 0:
				with open(blast_file_name, "r") as blastx_rslt:
					rslt = blastx_rslt.readline().strip('\n').split("\t")
					qseq_id = rslt[0]
					sseq_id = rslt[1]

					# check if the first result is id_of_interest
					if sseq_id == self.q_prot_ids[q_spec]:
						valid_hits.append(i)
		return valid_hits
	
	
	def process_longer(self, s_id, q_spec, q_prot_db):
		valid_hits = self.blastx(s_id, q_spec, q_prot_db)
		
		if len(valid_hits) > 1:
			max_i = None
			max_len = 0
			for i in valid_hits:
				curr_len = abs(self.man_anno_dict[q_spec][s_id][1][i] - self.man_anno_dict[q_spec][s_id][2][i])
				if curr_len > max_len:
					max_i = i
					max_len = curr_len
			
			self.add_to_keep(q_spec, s_id, max_i)
			
		elif len(valid_hits) == 0:
			self.add_to_remove(q_spec, s_id)
		
		
		
	def process_manual(self, s_id, q_spec, q_prot_db):
		valid_hits = self.blastx(s_id, q_spec, q_prot_db)
		
		if len(valid_hits) > 1:
			self.add_to_further_anno(q_spec, s_id)
		
		elif len(valid_hits) == 0:
			self.add_to_remove(q_spec, s_id)
		
		else:
			self.add_to_keep(q_spec, s_id, valid_hits[0])
			
			
	def process_all_seqs(self):
		further_anno = {}
		auto_anno = {}
		
		if len(self.to_further_anno.keys()) > 0:
			for q_spec in self.to_further_anno.keys():
				for s_id in self.to_further_anno[q_spec]:	
					if q_spec not in further_anno.keys():
						further_anno[q_spec] = {s_id : self.man_anno_dict[q_spec][s_id]}
					else:
						further_anno[q_spec][s_id] = self.man_anno_dict[q_spec][s_id]
		else:
			print("No scaffolds require further annotation.")
		
		if len(self.to_keep.keys()) > 0:
			for q_spec, hit_dict in self.man_anno_dict.items():
				for s_id, align_info in hit_dict.items():
					if q_spec in self.to_keep.keys() and s_id in self.to_keep[q_spec].keys():
						index =  self.to_keep[q_spec][s_id]
						name = self.man_anno_dict[q_spec][s_id][0]
						s_start = self.man_anno_dict[q_spec][s_id][1][index]
						s_stop = self.man_anno_dict[q_spec][s_id][2][index]
						frame = self.man_anno_dict[q_spec][s_id][3][index]
						q_start = self.man_anno_dict[q_spec][s_id][4][index]
						q_stop = self.man_anno_dict[q_spec][s_id][5][index]
						if q_spec not in auto_anno.keys():
							auto_anno[q_spec] = {s_id : [name, s_start, s_stop, frame, q_start, q_stop]}
						else:
							auto_anno[q_spec][s_id] = [name, s_start, s_stop, frame, q_start, q_stop]
		else:
			print("All scaffolds require further annotation.")
			
		if len(self.to_remove.keys()) > 1:
			print("All alignments of the following scaffolds failed reciprocal blast:")
			for q_spec in self.to_remove.keys():
				for s_id in self.to_remove[q_spec]: 
					print(s_id)
					

		return further_anno, auto_anno
	
		

annotator = MultiAlignAnno(man_anno_dict, "nucl", prot_file_paths, q_prot_ids)

for q_spec, hit_dict in annotator.man_anno_dict.items():
	
	q_prot_ds = annotator.prot_file_paths[q_spec]
	q_prot_db = get_dbs(q_prot_ds, "prot")
	
	for s_id, align_info in hit_dict.items():
		s_spec = align_info[0]
		s_starts = align_info[1]
		s_stops = align_info[2]
		frames = align_info[3]
		q_starts = align_info[4]
		q_stops = align_info[5]
		
		
		# more than one alignments on subject scaffold
		if len(s_starts) > 2:
			annotator.process_manual(s_id, q_spec, q_prot_db)
		
		# all frames are either positive or negative
		elif all(num > 0 for num in frames) or all(num < 0 for num in frames): 
			
			# overlap in subject sequence alignments and overlap in query sequence alignments
			if (annotator.alignments_overlap(s_starts, s_stops) and annotator.alignments_overlap(q_starts, q_stops)) or \
			(annotator.alignments_overlap(s_stops, s_starts) and annotator.alignments_overlap(q_stops, q_starts)):
				annotator.process_manual(s_id, q_spec, q_prot_db)
			
			# overlap in subject sequence alignments (different parts of query, overlap subject parts)
			elif annotator.alignments_overlap(s_starts, s_stops) or annotator.alignments_overlap(s_stops, s_starts):
				annotator.process_longer(s_id, q_spec, q_prot_db)		
			
			# overlap in query sequence alignments (different parts of subject, overlap query parts)
			elif annotator.alignments_overlap(q_starts, q_stops) or annotator.alignments_overlap(q_stops, q_starts):
				annotator.process_manual(s_id, q_spec, q_prot_db)
				
			# no overlap in query and subject alignments
			else:
				pass
			
		
		# mix of positive and minus frames
		else:
			annotator.process_longer(s_id, q_spec, q_prot_db)
		
further_anno, auto_anno = annotator.process_all_seqs()
for spec in further_anno.keys():
	print(len(further_anno[spec].keys()))
for spec in auto_anno.keys():
	print(len(auto_anno[spec].keys()))
