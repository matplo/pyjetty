#!/usr/bin/env python3

import os

def merge_filelists():
    
    path = '/rstorage/alice/data/LHC18qr/413-414'
    n = 100

    filenames_passed = ['{}/files_passed_{}.txt'.format(path, i) for i in range(1, n+1)]
    filenames_failed = ['{}/files_failed_{}.txt'.format(path, i) for i in range(1, n+1)]

    concatenate(filenames_passed, '{}/files_passed.txt'.format(path))
    concatenate(filenames_failed, '{}/files_failed.txt'.format(path))

    [os.remove(file) for file in filenames_passed]
    [os.remove(file) for file in filenames_failed]

def concatenate(file_list, output_file):

    with open(output_file, 'w') as outfile:
        for file in file_list:
            with open(file) as infile:
                for line in infile:
                    outfile.write(line)
                
if __name__ == '__main__':
    merge_filelists()
