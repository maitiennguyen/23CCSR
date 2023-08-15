import subprocess
import os
import sys

class Organizer():
	
	def __init__(self, q_num, taxIDs, ds_query, seq_query):
		self.q_num = q_num
		self.taxIDs = taxIDs
		self.ds_query = ds_query
		self.seq_query = seq_query
	
	# make folders
	def mkdir(self, fol_name):
		num = 1
		
		while os.path.exists(fol_name):
			fol_name = '{0}_{1}'.format(fol_name, num)
			num += 1
			
		subprocess.run("mkdir {0}".format(fol_name).split())
		
		return fol_name
	
		
	# move files into folders
	def mv(self, files, folder):
		files = ' '.join(files)
		
		subprocess.run("mv {0} {1}/".format(files, folder).split())
	
	
	# remove files
	def rm(self, files):
		files = ' '.join(files)
		
		subprocess.run("rm {0}".format(files).split())
		
	
	# organize files
	def organize_files(self):
		
		# add input files into one folder
		input_fol = self.mkdir('InputFiles')
		input_files = [self.seq_query, self.ds_query]
		self.mv(input_files, input_fol)
		
		
        # add required search set files into one folder (for future runs)
		search_set_fol = self.mkdir('SearchSetFiles')
		search_set_files = ['nucl.fna', 'prot.faa', 'all_specs.txt', 'prot_data_specs.txt', 'prot_files_all_dict.txt', ' '.join(self.taxIDs)]
		search_set_files2 = []
		for file in search_set_files:
			if os.path.exists(file):
				search_set_files2.append(file)
		self.mv(search_set_files2, search_set_fol)
        
		
        # remove very miscellaneous files
		dataset_files = ['.pjs', '.ptf', '.pto', '.pot', '.pdb', '.pos', '.pog', '.psq', '.phr', '.pin']
		dataset_files2 = ['.njs', '.ntf', '.nto', '.not', '.ndb', '.nos', '.nog', '.nsq', '.nhr', '.nin', '.fna']
		mis_files= []
		ds_name = self.ds_query.split('.')[0]
        
		for i in range(len(dataset_files)):
			file_name1 = ds_name + dataset_files[i]
			if os.path.exists(file_name1):
				mis_files.append(file_name1)
				
			file_name2 = 'prot' + dataset_files[i]
			if os.path.exists(file_name2):
				mis_files.append(file_name2)
		
		for i in range(len(dataset_files2)):
			file_name3 = 'nucl_2' + dataset_files2[i]
			if os.path.exists(file_name3):
				mis_files.append(file_name3)
				
			file_name4 = 'nucl' + dataset_files2[i]
			if os.path.exists(file_name4):
				mis_files.append(file_name4)
				
			file_name5 = ds_name + dataset_files2[i]
			if os.path.exists(file_name5):
				mis_files.append(file_name5)
				
		self.rm(mis_files)
		self.rm(['recip_seq.txt'])
		 
			
		# move dictionaries from all runs into one folder
		dict_fol = self.mkdir('DictReports')
		blast_dicts = ['blastp1', 'blastp2', 'tblastn', 'blastx', 'tblastx1', 'tblastx2']
		dict_files = []
		
		for i in range(len(blast_dicts)):
			for j in range(self.q_num):
				dict_file = blast_dicts[i] + '_' + str(j+1) + '_dict.txt'
				if os.path.exists(dict_file):
					dict_files.append(dict_file)
					
		self.mv(dict_files, dict_fol)
		
		
		# move summary reports from all runs into one folder
		sum_rep_fol = self.mkdir('SummaryReports')
		sum_rep_dicts = ['blastp', 'tblastn', 'blastx', 'tblastx']
		sum_rep_files = []
		
		for i in range(len(sum_rep_dicts)):
			for j in range(self.q_num):
				sum_rep_file = sum_rep_dicts[i] + '_' + str(j+1) + '_summary_report.txt'
				if os.path.exists(sum_rep_file):
					sum_rep_files.append(sum_rep_file)
					
		self.mv(sum_rep_files, sum_rep_fol)
		
		
		# move blast results from all runs into one folder
		blast_fol = self.mkdir('BlastFiles')
		subprocess.run("mv {0} {1}/".format('*.blasted', blast_fol), shell=True)
		
		
		# move main files into one folder
		main_fol = self.mkdir('MainFiles')
		main_files = ['complete_blast_summary_report.txt', 'complete_manual_anno_summary.txt', 'complete_auto_anno_summary.txt', 'auto_algn.fasta', 'auto_algn.clustal']
		main_files2 = []
		
		for file in main_files:
			if os.path.exists(file):
				main_files2.append(file)
		self.mv(main_files2, main_fol)
		
		
		# move annotation files into one folder
		anno_fol = self.mkdir('AnnotationFiles')
		subprocess.run("mv {0} {1}/".format('*.txt', anno_fol), shell=True)
		if os.path.exists('auto_anno_seqs.fasta'):
			self.mv(['auto_anno_seqs.fasta'], anno_fol)
		
		
		# move query seq files into one folder
		qseq_fol = self.mkdir('QuerySeqs')
		file_list = os.listdir('.')
		contains_fasta = any('.fasta' in file for file in file_list)
		if contains_fasta:
			subprocess.run("mv {0} {1}/".format('*.fasta', qseq_fol), shell=True)