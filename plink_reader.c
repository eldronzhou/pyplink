#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

int read_bed_file(const char* bed_file, int num_samples, int num_snps, \
	int num_snp_array, int num_sample_array, int *snp_array, int *sample_array, \
	float* genotype_matrix) {
    
    FILE* file = fopen(bed_file, "rb");
    if (!file) {
        fprintf(stderr, "Error: Cannot open .bed file.\n");
        return -1;
    }

    // Read magic number and SNP-major mode
    uint8_t magic[3];
    fread(magic, sizeof(uint8_t), 3, file);
    if (magic[0] != 0x6C || magic[1] != 0x1B || magic[2] != 0x01) {
        fprintf(stderr, "Error: Invalid .bed file format.\n");
        fclose(file);
        return -1;
    }

    int bytes_per_snp = (num_samples + 3) / 4;
    uint8_t* buffer = (uint8_t*)malloc(bytes_per_snp);

    for (size_t snp_idx = 0; snp_idx < num_snp_array; snp_idx++) {
        
	size_t snp_index = snp_array[snp_idx];

	long long file_offset = 0;
	if (snp_idx == 0) {
	    file_offset = (long long)snp_index * bytes_per_snp;	
	}
	else {
	    file_offset = (long long)(snp_index-snp_array[snp_idx-1]-1) * bytes_per_snp;
	}

	long chunk_size = 0x7FFFFFFF;
	long long remaining_offset = file_offset;

	while (remaining_offset > chunk_size) {
	    if (fseek(file, chunk_size, SEEK_CUR) != 0) {
		perror("Error seeking in file");
		fclose(file);
		return 1;
	    }
	    remaining_offset -= chunk_size;
	}

	if (fseek(file, (long)remaining_offset, SEEK_CUR) != 0) {
            perror("Failed to seek SNP position");
            free(buffer);
            fclose(file);
            return -1;
        }

	if (fread(buffer, 1, bytes_per_snp, file) != bytes_per_snp) {
            perror("Failed to read SNP data");
            free(buffer);
            fclose(file);
            return -1;
        }

	float count = 0, n_miss = 0;
	for (size_t sample_idx = 0; sample_idx < num_sample_array; sample_idx++) {
	    int sample_index = sample_array[sample_idx];
            int byte_index = sample_index / 4;
            int bit_offset = (sample_index % 4) * 2;
            uint8_t genotype_code = (buffer[byte_index] >> bit_offset) & 0b11;

            float value;
            if (genotype_code == 0b00) {
		value = 2.0; // Homozygous first
		count += 2;
	    }
            else if (genotype_code == 0b10) {
		value = 1.0; // Heterozygous
		count++;
	    }
            else if (genotype_code == 0b11) {
		value = 0.0; // Homozygous second
	    }
            else {
		value = -1.0; // Missing
		n_miss++;
	    }
            genotype_matrix[snp_idx * num_sample_array + sample_idx] = value;
        }

	// fill missing with freq
	float freq = count / (2 * (num_sample_array - n_miss) );
	float sd = sqrt(2 * freq * (1 - freq));
	for (size_t sample_idx = 0; sample_idx < num_sample_array; sample_idx++) {
	    if (genotype_matrix[snp_idx * num_sample_array + sample_idx] == -1) {
		genotype_matrix[snp_idx * num_sample_array + sample_idx] = 0;
	    }
	    else {
		genotype_matrix[snp_idx * num_sample_array + sample_idx] -= 2*freq;
		genotype_matrix[snp_idx * num_sample_array + sample_idx] /= sd;
	    }
	}
    
    }

    free(buffer);
    fclose(file);
    return 0;
}

