#!/bin/bash
# Proper header for a Bash script.

# tidy, version 2

# clean the large files that grow during MESA runs.
rm -f LOGS/*
rm -f LOGS1/*
rm -f LOGS2/*

rm -f photos/*
rm -f photos1/*
rm -f photos2/*

rm -f png/*
rm -f png1/*
rm -f png2/*


echo "cleaned photos*, LOGS*, png* directories"

exit #  The right and proper method of "exiting" from a script.
     #  A bare "exit" (no parameter) returns the exit status
     #+ of the preceding command. 
