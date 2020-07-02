A collection of NIM code, as well as a singularity collection containing the nim compiler and various nim packages

Download the container

singularity pull --name avocado_container.sif shub://J35P312/Nim_tes

Nim packages in the container

HTS-nim
xlsx
argparse

NOTE:

Remember to set the htslib path:

on uppmax:
module load bioinfo-tools htslib

else:

export LD_LIBRARY_PATH=/path/to/htslib/lib
