# (C) Copyright 2023- ECMWF.
# (C) Copyright 2023- Meteo-France.
import sys 
from sys import argv

verbose = True
#verbose = False
file_name = argv[1]
file_name = file_name.replace(".F90", "_openacc.F90")
with open(file_name, 'r') as file:
    lines = file.readlines()
new_file = []
for line in lines:
    if "acc routine" in line:
        if verbose: print(f"{line=}")
        line = line.rstrip() + " seq\n"
    new_file.append(line)
with open(file_name, 'w') as file:
    file.writelines(new_file)

