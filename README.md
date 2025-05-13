# pyplink
Self-contained python C interface to read plink bed files.

Usage:

First compile it with gcc `gcc -shared -o plink_reader.so -fPIC plink_reader.c`.

Define the function as below:

```
import ctypes
import numpy as np

# Define the function 
plink_reader.read_bed_file.argtypes = [
    ctypes.c_char_p,  # Path to the .bed file
    ctypes.c_int,     # Number of samples
    ctypes.c_int,     # Number of SNPs
    ctypes.c_int,     # Number of SNP array
    ctypes.c_int,     # Number of SNP sample array
    ctypes.POINTER(ctypes.c_int),  # SNP array
    ctypes.POINTER(ctypes.c_int),  # Sample array
    ctypes.POINTER(ctypes.c_float),  # Genotype matrix
]

# Load the shared library
plink_reader = ctypes.CDLL("./plink_reader.so")

# read function
def read_bed_file(bed_file, num_samples, num_snps, snp_array, sample_array):
    num_snp_array = len(snp_array)
    num_sample_array = len(sample_array)
        
    snp_array = (ctypes.c_int * num_snp_array)(*snp_array)
    sample_array = (ctypes.c_int * num_sample_array)(*sample_array)
    
    # Allocate memory for the genotype matrix
    print("Allocating memory for genotype matrix.")
    genotype_matrix = np.empty((num_snp_array, num_sample_array), dtype=np.float32)
    print("Finished allocation!")
    matrix_ptr = genotype_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    # Call the C function
    result = plink_reader.read_bed_file(
        bed_file.encode('utf-8'),
        num_samples,
        num_snps,
        num_snp_array,
        num_sample_array,
        snp_array,
        sample_array,
        matrix_ptr,
    )

    if result != 0:
        raise RuntimeError("Failed to read .bed file")

    return genotype_matrix
```

Use the following code to read the genotype matrix into python:

```
fam = pd.read_csv("XXX.fam", sep="\t", header=None)
num_fam = len(fam) # number of samples in fam
bim = pd.read_csv("XXX.bim", sep="\t", header=None)
sample_indices = [x for x in range(num_fam)] # sample indices
snp_indices = [x for x in range(num_snps)] # SNP indices
num_snps = len(snp_indices)
filename = "XXX.bed"
genotype_matrix = read_bed_file(filename, num_fam, num_snps, snp_indices, sample_indices)
```
