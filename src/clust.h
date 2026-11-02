#ifndef CLUST_H
#define CLUST_H

#include <stdio.h>
#include <stdint.h>

#include "shtab.h"

#define ABS_PEN_QUERY 0.1   // error when the reads show a SNP-mer that the reference sample lacks (NA in the matrix)
#define ABS_PEN_REF 0.1     // error when the reference sample has a SNP-mer that the reads do not show

typedef struct {
    char *chr;
    int n_snps, n_haps, n_cls, n_geno;
    char **hap_names;
    char **cls_name;        // n_cls, haplotype class names ("1", "2", ...)
    int *hap_cls;           // n_haps, haplotype class from the labels file (-1 = unlabelled), index to cls_name
    char **geno_names;      // n_geno, genotype names ("1/1", "1/2", "2/2", ...)
    int *gt_map;            // converts a pair of haplotypes into a genotype id, index to geno_names
    char **snpmer;
    uint64_t **h_snpmer;
    uint8_t *miss;
    float *p;               // n_snps * n_haps, mean expected composition of the sample (0 where missing)
} geno_mat_t;

typedef struct {
    int nn_a, nn_b;             // nearest haplotype pair among the winning genotype's voters
    int gt_idx;                 // winning genotype (index into geno_names, -1 if unresolved)
    long n_pairs;               // labelled haplotype pairs compared against the query
    double dist_best;           // distance to nn_a, nn_b (inf if unresolved)
    int k_used;                 // neighbours that actually voted (min(k, labelled samples with a distance))
    int votes;                  // votes for the winning genotype
    int tie;                    // 1 if another genotype had the same number of votes (broken by distance)
    double dist_vote;           // average distance of the winning genotype's voters
} geno_call_t;

typedef struct { // map of labels tsv file
    char **hap, **geno, **chr;
    int n;
} labels_t;

typedef struct { // a block of SNP-mers for a chromosome
    geno_mat_t *gm;
    int m_snp;
} blk_t;

geno_mat_t **load_geno_matrix(pg_msht_t *h, const char *fn, const char *labels_fn, int *n_out);
void destroy_geno_matrices(geno_mat_t **gms, int n);
uint64_t *precompute_snpmer_hash(pg_msht_t *h, const char *snpmer);
int snpmer_lookup(pg_msht_t *h, const uint64_t *h_snpmer, uint32_t *cnt1, uint32_t *cnt2);

// Nearest-neighbour call: the genotype of the reference sample whose composition profile is closest
void classify_genotype(pg_msht_t *h, const geno_mat_t *gm, int k, geno_call_t *out, int use_asm);
void classify_query(const geno_mat_t *gm, const double *q, int k, geno_call_t *out);
void print_geno_call(FILE *out, const geno_mat_t *gm, const geno_call_t *c);

#endif // CLUST_H