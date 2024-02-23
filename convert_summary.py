import csv
import sys
import re

def main(argv):
    
    input_file = argv[1]
    output_file = 'converted_summary.csv'

    with open(input_file, 'r') as in_file, open(output_file, 'w', newline='') as out_file:
        writer = csv.writer(out_file)
        for line in in_file:
            row = re.split('\s{2,}', line.strip())
            writer.writerow(row)

    
    
if __name__ == "__main__":
	main(sys.argv)
    
