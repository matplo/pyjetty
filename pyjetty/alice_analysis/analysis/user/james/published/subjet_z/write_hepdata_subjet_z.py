import hepdata_lib
import pickle

input_files = ['/Users/jamesmulligan/Analysis_subjet_z/output_analysis/pp/493677-04-paper/inclusive_subjet_z/final_results/hepdata/tables.pkl',
               '/Users/jamesmulligan/Analysis_subjet_z/output_analysis/pp/493677-04-paper/leading_subjet_z/final_results/hepdata/tables.pkl',
               '/Users/jamesmulligan/Analysis_subjet_z/paper/fig_v5/hepdata/tables.pkl']
output_dir = '/Users/jamesmulligan/Analysis_subjet_z/paper/hepdata/20220318_hepdata'


# Get all tables, and store them in a list
key_list = []
table_list = []
for input_file in input_files:

    with open(input_file, 'rb') as f:
        table_dict = pickle.load(f)

        for key,table in table_dict.items():
            key_list.append(key)
            table_list.append(table)

# Sort tables
sorted_keys, sorted_tables = zip(*sorted(zip(key_list, table_list)))

# Create submission
hepdata_submission = hepdata_lib.Submission()
for table in sorted_tables:
    hepdata_submission.add_table(table)
hepdata_submission.create_files(output_dir)