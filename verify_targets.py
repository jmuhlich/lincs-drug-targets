import sys
import os
import csv

verify_columns = ('drug_id', 'drug_moniker', 'protein_id', 'protein_moniker')

master_filename = os.path.join(os.path.dirname(sys.argv[0]), 'targets.tsv')
curated_filename = sys.argv[1]

master_file = open(master_filename, 'rb')
curated_file = open(curated_filename, 'rb')

master_reader = csv.DictReader(master_file, dialect='excel-tab')
curated_reader = csv.DictReader(curated_file, dialect='excel-tab')

for master_row in master_reader:
    curated_row = curated_reader.next()
    for column in verify_columns:
        master_data = master_row[column]
        curated_data = curated_row[column]
        if master_data != curated_data:
            print ("%s: Mismatch in column '%s' at line %d. "
                   "Should be '%s', got '%s' instead." %
                   (curated_filename, column, curated_reader.line_num,
                    master_data, curated_data))
            sys.exit(1)
