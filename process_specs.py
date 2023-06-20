import random
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import NCBITaxa

class CladesProcessor():
	
	def __init__(self, results_dict, db, prev_q_specs):	
		self.specs = []
		for q_spec in results_dict.keys():
			for seq_id in results_dict[q_spec]:
				spec = results_dict[q_spec][seq_id][0][0]
				if spec not in self.specs:
					self.specs.append(spec)
				
		self.results_dict = results_dict
		self.db = db
		self.prev_q_specs = prev_q_specs


	# get tax id for each species	
	def get_taxid(self, spec_name):
		ncbi = NCBITaxa()
		
		taxids = ncbi.get_name_translator([spec_name])
		if spec_name in taxids:
			taxid = taxids[spec_name][0]
		else:
			return None
		
		lineage = ncbi.get_lineage(taxid)
		
		family_taxid = None
		
		for tax in lineage:
			rank = ncbi.get_rank([tax])
			if rank[tax] == "family":
				family_taxid = tax
				break

		return family_taxid
	
	
	def get_query_spec_taxids(self):
		prev_q_taxids = []
		for spec in self.prev_q_specs:
			spec_taxid = self.get_taxid(spec)
			prev_q_taxids.append(spec_taxid)
		return prev_q_taxids
	
	# sort species into clades
	def process_all_specs(self):
		query_taxids = self.get_query_spec_taxids()
		query_clades = []
		other_clades = []
		diff_genus_specs = [spec for spec in self.specs for q_spec in self.prev_q_specs if spec.split()[0] != q_spec.split()[0]]
		
		for spec in diff_genus_specs:
			taxid = self.get_taxid(spec)
			if taxid in query_taxids:
				query_clades.append(spec)
			else:
				other_clades.append(spec)
		
		return query_clades, other_clades
	
	# pick random species in list
	def get_rand_spec(self):
		query_clades, other_clades = self.process_all_specs()
		
		if len(other_clades) == 0:
			return None
		else:
			choice = random.choice(other_clades)
			return choice
	
	# get seq id using species name
	def get_id(self, spec_name):
		for q_spec in self.results_dict.keys():
			for seq_id in self.results_dict[q_spec]:
				if spec_name == self.results_dict[q_spec][seq_id][0][0]:
					return seq_id
		return None
	
	# get fasta file of a given seq id
	def get_id_fasta(self, seq_id):
		db_info = subprocess.run("blastdbcmd -db {0} -entry {1}".format(self.db, seq_id).split(), capture_output=True, text=True).stdout.split("\n")
		des = db_info[0][1:]
		seq = Seq(''.join(db_info[1:-1]))
		
		fasta_seq = SeqIO.SeqRecord(seq, id=seq_id, description=des)
		
		filename = seq_id + "_query.fasta"
		with open(filename, "w") as fasta_file:
				SeqIO.write(fasta_seq, fasta_file, "fasta")
				
		return filename
		
												  
		