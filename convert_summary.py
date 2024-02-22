import csv
import sys

column_widths = [50, 30, 100, 100]\

headers = []
lines = []

def main(argv):
    
    filename = argv[1]

    with open(filename, 'r') as file:
        for line_num, line in enumerate(file):
            if line_num == 0:
                headers.append(line[:column_widths[0]].strip())
                headers.append(line[column_widths[1]:column_widths[1]+column_widths[0]+column_widths[1]].strip()[:10])
                headers.append(line[column_widths[0]+column_widths[1]:].strip())
            else:
                row = []
                row.append(line[:column_widths[0]].strip())
                row.append(line[column_widths[1]:column_widths[1]+column_widths[0]+column_widths[1]].strip()[:1])
                row.append(line[column_widths[0]+column_widths[1]:].strip())
                lines.append(row)

    with open('./converted_summary.csv', 'w', newline ='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows([headers])
        writer.writerows(lines)

    
    
if __name__ == "__main__":
	main(sys.argv)
    
