#!/bin/bash

#Generate code to test (in loki and tmp repo)
./openacc.sh
./cpg_dyn_slg.sh

#Run test
./diff_openacc.sh
./diff_cpg_dyn_slg.sh
./diff_tmp.sh
