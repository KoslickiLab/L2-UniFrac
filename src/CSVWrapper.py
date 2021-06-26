################################################################################
##
##  File           : write_to_csv.py
##  Description    : Writes UniFrac data to CSV file periodically. 
##
##   Author        : *** Andrew Millward ***
##   Last Modified : *** 05/30/2021 ***
##

## Import Files
import csv

################################################################################
##
## Function     : write_spreadsheet
## Description  : Takes input data and exports it to a new spreadsheet.
##
## Inputs       : elem_list - List of elements to write from query.
## Outputs      : 0 if successful, -1 if failure.

def write(name, dist_list):

    try:
        with open(name, 'a', newline='') as csvfile:
            file_write = csv.writer(csvfile)
            file_write.writerow(dist_list)
            csvfile.close()

        return ( 0 ) # Success

    except:
        return ( -1 ) # Failure

if __name__ == "__main__":
    write('distance_matrix.csv', [0.3, 0.2, 0.4])