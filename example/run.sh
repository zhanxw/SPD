#!/bin/bash

# 1. Build executable
make -C ../

# 2. Perform calculation
../main infile.txt outfile.txt result
