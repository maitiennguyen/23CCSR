import subprocess
import sys

def main(argv):
	filename = argv[1]
	
	outname = filename + "_algn"
	
	subprocess.run("clustalo -i {0} -o {1} --outfmt=clustal --resno --threads=16".format(filename, outname + ".clustal").split())
	subprocess.run("clustalo -i {0} -o {1} --outfmt=fasta --resno --threads=16".format(filename, outname + ".fasta").split())
	
if __name__ == "__main__":
	main(sys.argv)
	
	