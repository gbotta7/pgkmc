#ifndef CLUST_H
#define CLUST_H

#include <stdio.h>
#include <stdint.h>

#include "shtab.h"

#define NOISE_THRESHOLD 3
#define ABS_PEN_QUERY 0.1   // error when the reads show a SNP-mer that the reference sample lacks (NA in the matrix)
#define ABS_PEN_REF 0.1     // error when the reference sample has a SNP-mer that the reads do not show
                            // (a squared composition difference is at most 1). Keep equal to ABS_PEN_* in compute_cenhaps.py

typedef struct {
    uint64_t *h_snpmer;     // NULL if the SNP-mer could not be hashed
} marker_t;

// Reference panel for one chromosome: rows = SNP-mers, columns = reference samples.
// Each sample carries the genotype label of its haplotype pair (e.g. "1/2").
typedef struct {
    char *chr;
    int n_snp, n_samp, n_geno;
    char **samp_names;      // n_samp
    int *samp_geno;         // n_samp, index into geno_names; -1 = sample has no genotype (ignored)
    char **geno_names;      // n_geno distinct genotype labels
    char **snpmer;
    uint64_t **h_snpmer;
    float *p;               // n_snp * n_samp, mean expected composition of the sample (0 where missing)
    uint8_t *miss;          // n_snp * n_samp, 1 = no profile for this sample/SNP-mer
} geno_mat_t;

typedef struct {
    int nn;                     // nearest reference sample (column index, -1 if unresolved)
    int geno;                   // genotype of the nearest sample (index into geno_names, -1 if unresolved)
    int nn2;                    // nearest sample whose genotype differs from `geno` (-1 if none)
    double dist_best, dist_second; // distance to nn / nn2 (lower = closer); inf if unresolved / none
    int n_shared;               // SNP-mers used for the distance to nn
    int n_obs;                  // SNP-mers with counts > 0 in the sample
    int n_missing;              // SNP-mers flagged missing (both counts 0)
    int n_used;                 // SNP-mers scored (observed in the sample)
    uint64_t depth;
} geno_call_t;

uint64_t *precompute_snpmer_hash(pg_msht_t *h, const char *snpmer);
int snpmer_lookup(pg_msht_t *h, marker_t *mk, uint32_t *cnt1, uint32_t *cnt2);

// Reads a TSV holding one or more concatenated per-chromosome matrices
// ("snpmer<TAB>sample1<TAB>sample2...<TAB>chr", one row per SNP-mer).
// labels_fn: "sample<TAB>genotype" TSV (genotype = label such as "1/2", "?" = not genotyped)
// giving the genotype of each matrix column. Required.
// Name of a reference sample to leave out of every panel (leave-one-out: the sample being genotyped
// is itself a column of the matrices). NULL = use all reference samples.
void clust_set_exclude(const char *name);
geno_mat_t **load_geno_matrices(pg_msht_t *h, const char *fn, const char *labels_fn, int *n_out);
// Reads a text file with one matrix path per line, and a text file with one labels
// path per line (same order), and loads every chromosome block found into a single array.
geno_mat_t **load_matrices_from_list(pg_msht_t *h, const char *list_fn, const char *labels_list_fn, int *n_out);
void destroy_geno_matrices(geno_mat_t **gms, int n);

// Nearest-neighbour call: the genotype of the reference sample whose composition profile is closest
// (RMSE over all SNP-mers, absence mismatches penalised) to the composition observed in the sequenced sample.
void classify_genotype(pg_msht_t *h, const geno_mat_t *gm, geno_call_t *out);
void print_geno_call(FILE *out, const geno_mat_t *gm, const geno_call_t *c);

#endif // CLUST_H