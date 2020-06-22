#! /bin/bash

# Use find command to generate text file with list of all ROOT files
# contained within the current directory.

# Note: Execute this from /rstorage, so that the paths will be those
#       of the cross-mounted read-only storage space

find "$PWD" -name "*.root" >files.txt

# To see number of files: wc -l files.txt
# To remove empty directories: find . -type d -empty -delete
