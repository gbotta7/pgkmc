#ifndef CLUST_H
#define CLUST_H

#include <stdint.h>

#include "shtab.h"

typedef struct {
    int parent_id;
    char side[8];
    int left_child;
    int right_child;
} tnode_t;

typedef struct {
    uint64_t *h_snpmer;
    char move[16];
} marker_t;

typedef struct {
    marker_t *a;
    int n, m;
} marker_list_t;

typedef struct {
    tnode_t *nodes;          // 1-indexed, size cap+1
    marker_list_t *markers;  // 1-indexed, size cap+1
    int n_nodes;
    int cap;                 // allocated capacity (nodes may still be growing during load)
} clust_tree_t;

typedef struct {
    char *chr;
    clust_tree_t tree;
} chr_tree_t;

typedef struct {
    chr_tree_t *a;
    int n, m;
} chr_forest_t;

typedef struct {
	int left_votes, right_votes;
} vote_t;

typedef struct {
	int node_id;
	char side;
	int matches;
} clust_result_t;

typedef struct {
    clust_result_t *a;
    int n, m;
} result_list_t;

typedef struct {
    char *chr;
    clust_result_t best;
} chr_result_t;

clust_tree_t *forest_get_or_create(chr_forest_t *forest, const char *chr);
void tree_ensure_capacity(clust_tree_t *t, int node_id);
uint64_t *precompute_snpmer_hash(pg_msht_t *h, const char *snpmer);
void push_marker(pg_msht_t *h, marker_list_t *ml, const char *snpmer, const char *move);
chr_forest_t *load_clust_forest(pg_msht_t *h, const char *structure_fn, const char *snpmer_fn);
void destroy_forest(chr_forest_t *forest);
int snpmer_lookup(pg_msht_t *h, marker_t *mk, uint32_t *cnt1, uint32_t *cnt2);
void vote_branch(pg_msht_t *h, marker_t *mk, vote_t *vt);
void push_result(result_list_t *rl, int node_id, char side, int matches);
void walk_tree(pg_msht_t *h, clust_tree_t *t, int node_id, int matches, result_list_t *results);
clust_result_t pick_best(result_list_t *rl);
chr_result_t *classify_cenhaps(pg_msht_t *h, chr_forest_t *forest, int *n_out);

#endif // CLUST_H
