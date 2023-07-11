import sys

class InputProcessor():
	
	def __init__(self, argv):
		self.argv = argv
		self.req_args = ["-qseq", "-qdb", "-qtype", "-qname", "-sset", "-download"]
		self.opt_args = ["-evalue", "-threads"]
		
	def valid_index(self, arg, index):
		if index == len(self.argv) - 1:
			print("\nERROR: No argument found: " + arg + "\n")
			sys.exit()
			
	def valid_arg(self, arg, para):
		if para in self.req_args or para in self.opt_args:
			print("\nERROR: Flag '" + para + "' found in place of argument: " + arg + "\n")
			sys.exit()
		
	def get_qseq(self):
		qseq = None
		
		for index, arg in enumerate(self.argv):
			if arg == "-qseq":
				self.valid_index(arg, index)
				qseq = self.argv[index+1]
				self.valid_arg(arg, qseq)
				break
				
		if qseq is None:
			print("\nERROR: Query sequence not found. Use '-qseq' to indicate query sequence file. See README.md for usage.\n")
			sys.exit()
		else:
			return qseq
		
		
	def get_qdb(self):
		qdb = None
		
		for index, arg in enumerate(self.argv):
			if arg == "-qdb":
				self.valid_index(arg, index)
				qdb = self.argv[index+1]
				self.valid_arg(arg, qdb)
				break
				
		if qdb is None:
			print("\nERROR: Query nucl/prot dataset not found. Use '-qdb' to indicate query dataset file. See README.md for usage.\n")
			sys.exit()
		else:
			return qdb
		
		
	def get_qtype(self):
		qtype = None
		
		for index, arg in enumerate(self.argv):
			if arg == "-qtype":
				self.valid_index(arg, index)
				qtype = self.argv[index+1]
				self.valid_arg(arg, qtype)
				break
				
		if qtype is None:
			print("\nERROR: Query type not found. Use '-qtype' to indicate query type. See README.md for usage.\n")
			sys.exit()
		else:
			qtype = qtype.lower()
			if qtype == "prot" or qtype == "nucl":
				return qtype
			else:
				print("\nERROR: Invalid input for qtype parameter. Use 'prot' or 'nock'. See README.md for usage.\n")
				sys.exit()
		
		
	def get_qname(self):
		qname = None
		
		for index, arg in enumerate(self.argv):
			if arg == "-qname":
				self.valid_index(arg, index)
				qname = ' '.join(self.argv[index+1].split(','))
				self.valid_arg(arg, qname)
				break
				
		if qname is None:
			print("\nERROR: Query species name not found. Use '-qname' to indicate query species name. See README.md for usage.\n")
			sys.exit()
		else:
			return qname
		
		
	def get_sset(self):
		sset = None
		
		for index, arg in enumerate(self.argv):
			if arg == "-sset":
				self.valid_index(arg, index)
				sset = self.argv[index+1].split(',')
				self.valid_arg(arg, sset)
				break
				
		if sset is None:
			print("\nERROR: Search set not found. Use '-sset' to indicate search set. See README.md for usage.\n")
			sys.exit()
		else:
			return sset
		
		
	def get_download(self):
		download = None
		
		for index, arg in enumerate(self.argv):
			if arg == "-download":
				self.valid_index(arg, index)
				download = self.argv[index+1]
				self.valid_arg(arg, download)
				break
				
		if download is None:
			print("\nERROR: No indication whether to download search set data. Use '-download' to indicate whether data should be downloaded. See README.md for usage.\n")
			sys.exit()
		else:
			download = download.lower()
			if download == "yes" or download == "no":
				return download
			else:
				print("\nERROR: Invalid input for download parameter. Use 'yes' or 'no'. See README.md for usage.\n")
				sys.exit()
		
		
	def get_evalue(self):
		evalue = "1e-01"  
		
		for index, arg in enumerate(self.argv):
			if arg == "-evalue":
				self.valid_index(arg, index)
				evalue = self.argv[index+1]
				self.valid_arg(arg, evalue)
				break
		
		try:
			float(evalue)
		except ValueError:
			print("\nERROR: Invalid input for e-value. See README.md for usage.\n")
			sys.exit()
			
		return evalue
	
	def get_threads(self):
		thread = "1" 
		
		for index, arg in enumerate(self.argv):
			if arg == "-threads":
				self.valid_index(arg, index)
				thread = self.argv[index+1]
				self.valid_arg(arg, thread)
				break
		
		try:
			int(thread)
		except ValueError:
			print("\nERROR: Invalid input for threads. See README.md for usage.\n")
			sys.exit()
			
		if int(thread) == 0:
			print("\nERROR: Invalid input for threads. See README.md for usage.\n")
			sys.exit()
			
		return thread
	
	def check_invalid_flag(self):
		for arg in self.argv:
			if "-" in arg:
				if arg not in self.req_args and arg not in self.opt_args:
					print("\nInvalid flag detected: " + arg + "\n")
					sys.exit()
