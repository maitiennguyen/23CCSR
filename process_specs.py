import random
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import NCBITaxa

class CladesProcessor():
	
	def __init__(self, q_spec_name, results_dict, db):
		self.query_spec = q_spec_name
		
		self.specs = []
		for seq_id in results_dict.keys():
			spec = results_dict[seq_id][0][0]
			if spec not in self.specs:
				self.specs.append(spec)
				
		self.results_dict = results_dict
		self.seq_ids = results_dict.keys()
		self.db = db


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
	
	
	def get_query_spec_taxid(self):
		return self.get_taxid(self.query_spec)
	
	# sort species into clades
	def process_all_specs(self):
		query_taxid = self.get_query_spec_taxid()
		query_clade = []
		other_clades = []
		
		for spec in self.specs:
			if spec.split()[0] != self.query_spec.split()[0]:
				taxid = self.get_taxid(spec)
				if taxid == query_taxid:
					query_clade.append(spec)
				else:
					other_clades.append(spec)
			
		return query_clade, other_clades
	
	# pick random species in list
	def get_rand_spec(self):
		query_clade, other_clades = self.process_all_specs()
		
		if len(other_clades) == 0:
			return None
		else:
			choice = random.choice(other_clades)
			return choice
	
	# get seq id using species name
	def get_id(self, spec_name):
		for seq_id in self.results_dict.keys():
			if spec_name == self.results_dict[seq_id][0][0]:
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
		
												  
		