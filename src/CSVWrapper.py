import csv
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix

def write(name, dist_list):
	try:
		with open(name, 'a', newline='') as csvfile:
			file_write = csv.writer(csvfile)
			file_write.writerow(dist_list)
			csvfile.close()

		return 0

	except Exception as exception:
		print(exception)
		return -1

def read(name):
	try:
		f = open(name, 'r')
		csv_read = csv.reader(f, delimiter=';')
		matrix = []
		for line in csv_read:
			matrix.append(list(map(float, line[0].split(","))))
		return matrix
	except Exception as exception:
		print(exception)
		return -1

def read_sparse(name):
	try:
		f = open(name, 'r')
		csv_read = csv.reader(f, delimiter=';')
		row = []
		col = []
		data = []
		for line in csv_read:
			line_split = line[0].split(",")
			if len(line_split) == 3:
				row.append(int(line_split[0]))
				col.append(int(line_split[1]))
				data.append(float(line_split[2]))
			else:
				rows = int(line_split[0])
				cols = int(line_split[1])
		row = np.array(row)
		col = np.array(col)
		data = np.array(data)
		sparse_matrix = csr_matrix((data, (row, col)), shape=(rows, cols))
		return sparse_matrix

	except Exception as exception:
		print(exception)
		return -1

def read_sparse_row(name, rowID):
	pass

if __name__ == "__main__":
	#write('distance_matrix.csv', [0.3, 0.2, 0.4])
	print(read_sparse('../scripts/L1-Push-Out.csv')[0])