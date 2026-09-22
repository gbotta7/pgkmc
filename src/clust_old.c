#include <stdio.h>

#include "clust.h"
#include "khashl.h" // hash table
#include "shtab.h"
#include "utils.h"

// ---- SNP-mer decision tree for clustering ----
clust_tree_t *forest_get_or_create(chr_forest_t *forest, const char *chr)
{
    for (int i = 0; i < forest->n; ++i)
        if (strcmp(forest->a[i].chr, chr) == 0)
            return &forest->a[i].tree;

    if (forest->n == forest->m) {
        forest->m = forest->m < 8 ? 8 : forest->m + (forest->m >> 1);
        REALLOC(forest->a, forest->m);
    }
    chr_tree_t *ct = &forest->a[forest->n++];
    ct->chr = strdup(chr);
    memset(&ct->tree, 0, sizeof(clust_tree_t));
    return &ct->tree;
}

void tree_ensure_capacity(clust_tree_t *t, int node_id)
{
    if (node_id <= t->cap) return;
    int new_cap = t->cap;
    while (node_id > new_cap) new_cap = new_cap < 16 ? 16 : new_cap + (new_cap >> 1);
    REALLOC(t->nodes, new_cap + 1);
    REALLOC(t->markers, new_cap + 1);
    memset(t->nodes + t->cap + 1, 0, (new_cap - t->cap) * sizeof(tnode_t));
    memset(t->markers + t->cap + 1, 0, (new_cap - t->cap) * sizeof(marker_list_t));
    t->cap = new_cap;
}

uint64_t *precompute_snpmer_hash(pg_msht_t *h, const char *snpmer)
{
    int half = h->k >> 1;
    if ((int)strlen(snpmer) != h->k + 4) return NULL;

    char left[64], right[64];
    memcpy(left, snpmer, half);              left[half] = '\0';
    memcpy(right, snpmer + half + 5, half);  right[half] = '\0';

    const char *b = snpmer + half;
    if (b[0] != '[' || b[2] != '/' || b[4] != ']') return NULL;
    char a1 = b[1];

    uint64_t x[2] = {0, 0};
    uint64_t mask = (1ULL << (h->k * 2)) - 1;
    uint64_t shift = (uint64_t)(h->k - 1) * 2;

    for (int i = 0; i < half; ++i) {
        int c = seq_nt4_table[(uint8_t)left[i]];
        if (c >= 4) return NULL;
        x[0] = (x[0] << 2 | c) & mask;
        x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
    }

    uint64_t cb1 = seq_nt4_table[(uint8_t)a1];
    if (cb1 >= 4) return NULL;
    x[0] = (x[0] << 2 | cb1) & mask;
    x[1] = x[1] >> 2 | (uint64_t)(3 - cb1) << shift;

    for (int i = 0; i < half; ++i) {
        int c = seq_nt4_table[(uint8_t)right[i]];
        if (c >= 4) return NULL;
        x[0] = (x[0] << 2 | c) & mask;
        x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
    }

    uint64_t y = x[0] < x[1] ? x[0] : x[1];
    uint64_t y_rev = x[0] < x[1] ? x[1] : x[0];
    uint64_t flanks = (y & ((1ULL << half * 2) - 1)) | ((y >> ((half + 1) * 2)) << (half * 2));
    uint64_t rev_flanks = (y_rev & ((1ULL << half * 2) - 1)) | ((y_rev >> ((half + 1) * 2)) << (half * 2));
    if (flanks == rev_flanks) {
        fprintf(stderr, "[E::%s] SNP-mer '%s' is palindromic and cannot be resolved\n", __func__, snpmer);
        return NULL;
    }

    uint64_t hash_mask = (1ULL << ((h->k - 1) * 2)) - 1;
    uint64_t *h_flanks;
    MALLOC(h_flanks, 1);
    *h_flanks = pg_hash64(flanks, hash_mask);
    return h_flanks;
}

void push_marker(pg_msht_t *h, marker_list_t *ml, const char *snpmer, const char *move)
{
    if (ml->n == ml->m) {
        ml->m = ml->m < 8 ? 8 : ml->m + (ml->m >> 1);
        REALLOC(ml->a, ml->m);
    }
    marker_t *mk = &ml->a[ml->n];
    strncpy(mk->move, move, sizeof(mk->move) - 1);
    mk->move[sizeof(mk->move) - 1] = '\0';
    mk->h_snpmer = precompute_snpmer_hash(h, snpmer);
    ml->n++;
}

chr_forest_t *load_clust_forest(pg_msht_t *h, const char *structure_fn, const char *snpmer_fn)
{
    chr_forest_t *forest;
    CALLOC(forest, 1);

    FILE *fp = fopen(structure_fn, "r");
    if (!fp) {
        fprintf(stderr, "[E::%s] failed to open '%s'\n", __func__, structure_fn);
        free(forest);
        return NULL;
    }

    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    int64_t n_parsed;

    // structure file columns: node_id  parent_id  side  left_child  right_child  chr
    n_parsed = 0;
    while ((len = getline(&line, &cap, fp)) >= 0) {
        char chr[64], side[8];
        int node_id, parent_id, left_c, right_c;
        if (sscanf(line, "%d\t%d\t%7s\t%d\t%d\t%63s",
                   &node_id, &parent_id, side, &left_c, &right_c, chr) != 6)
            continue;

        clust_tree_t *t = forest_get_or_create(forest, chr);
        tree_ensure_capacity(t, node_id);

        t->nodes[node_id].parent_id = parent_id;
        strncpy(t->nodes[node_id].side, side, sizeof(t->nodes[node_id].side) - 1);
        t->nodes[node_id].left_child = left_c;
        t->nodes[node_id].right_child = right_c;
        if (node_id > t->n_nodes) t->n_nodes = node_id;
        n_parsed++;
    }
    fclose(fp);

    if (n_parsed == 0)
        fprintf(stderr, "[W::%s] no rows parsed from '%s' — check column order/format\n", __func__, structure_fn);

    fp = fopen(snpmer_fn, "r");
    if (!fp) {
        fprintf(stderr, "[E::%s] failed to open '%s'\n", __func__, snpmer_fn);
        free(line);
        return forest;
    }

    // snpmer file columns: SNPmer  node_id  move  chr
    n_parsed = 0;
    while ((len = getline(&line, &cap, fp)) >= 0) {
        char chr[64], snpmer[128], move[16];
        int node_id;
        if (sscanf(line, "%127s\t%d\t%15s\t%63s", snpmer, &node_id, move, chr) != 4)
            continue;

        clust_tree_t *t = forest_get_or_create(forest, chr);
        if (node_id < 1 || node_id > t->n_nodes) continue; // structure row must exist already
        push_marker(h, &t->markers[node_id], snpmer, move);
        n_parsed++;
    }
    fclose(fp);
    free(line);

    if (n_parsed == 0)
        fprintf(stderr, "[W::%s] no rows parsed from '%s' — check column order/format\n", __func__, snpmer_fn);

    return forest;
}

void destroy_forest(chr_forest_t *forest)
{
    if (!forest) return;
    for (int c = 0; c < forest->n; ++c) {
        clust_tree_t *t = &forest->a[c].tree;
        for (int i = 1; i <= t->n_nodes; ++i) {
            for (int j = 0; j < t->markers[i].n; ++j)
                free(t->markers[i].a[j].h_snpmer);
            free(t->markers[i].a);
        }
        free(t->markers);
        free(t->nodes);
        free(forest->a[c].chr);
    }
    free(forest->a);
    free(forest);
}

int snpmer_lookup(pg_msht_t *h, marker_t *mk, uint32_t *cnt1, uint32_t *cnt2)
{
    if (!mk->h_snpmer) return -1;
    int bucket = *mk->h_snpmer & ((1 << h->pre) - 1);
    uint64_t key = (*mk->h_snpmer >> h->pre) << F_VAL_INFO_BITS;

    pg_sht1_t *g = &h->h[bucket];
    khint_t k = pg_sht_get(g->h, key);
    if (k == kh_end(g->h)) { *cnt1 = 0; *cnt2 = 0; return 0; }

    uint32_t v = kh_val(g->h, k);
    *cnt1 = s_val_count1(v);
    *cnt2 = s_val_count2(v);
    return 1;
}

void vote_branch(pg_msht_t *h, marker_t *mk, vote_t *vt)
{
    uint32_t cnt1 = 0, cnt2 = 0;
    int found = snpmer_lookup(h, mk, &cnt1, &cnt2);
    if (cnt1 < NOISE_THRESHOLD && cnt2 < NOISE_THRESHOLD) found = 0;
    int total = found == 1 ? (int)(cnt1 + cnt2) : 0;
    int has_call = total > 0;

    if (strcmp(mk->move, "left_none") == 0) {
        has_call ? vt->right_votes++ : vt->left_votes++;
        return;
    }
    if (strcmp(mk->move, "right_none") == 0) {
        has_call ? vt->left_votes++ : vt->right_votes++;
        return;
    }
    if (!has_call) return; // left_zero/left_one need an actual call

    double comp = (double)cnt1 / (double)total;
	if (strcmp(mk->move, "left_one") == 0) {
        if (comp == 0.0) vt->right_votes++;
        if (comp == 1.0) vt->left_votes++;
        return;
    }
	if (strcmp(mk->move, "left_zero") == 0) {
        if (comp == 0.0) vt->left_votes++;
        if (comp == 1.0) vt->right_votes++;
        return;
    }
}

void push_result(result_list_t *rl, int node_id, char side, int matches)
{
    if (rl->n == rl->m) {
        rl->m = rl->m < 8 ? 8 : rl->m + (rl->m >> 1);
        REALLOC(rl->a, rl->m);
    }
    rl->a[rl->n++] = (clust_result_t){node_id, side, matches};
}

void walk_tree(pg_msht_t *h, clust_tree_t *t, int node_id, int matches, result_list_t *results)
{
    tnode_t *nd = &t->nodes[node_id];
    marker_list_t *ml = &t->markers[node_id];

    vote_t vt = {0, 0};
    for (int i = 0; i < ml->n; ++i)
        vote_branch(h, &ml->a[i], &vt);

    if (vt.left_votes == vt.right_votes) {
        if (nd->left_child == 0) push_result(results, node_id, 'L', matches);
        else walk_tree(h, t, nd->left_child, matches, results);
        if (nd->right_child == 0) push_result(results, node_id, 'R', matches);
        else walk_tree(h, t, nd->right_child, matches, results);
    } else if (vt.left_votes > vt.right_votes) {
        int m2 = matches + vt.left_votes;
        if (nd->left_child == 0) push_result(results, node_id, 'L', m2);
        else walk_tree(h, t, nd->left_child, m2, results);
    } else {
        int m2 = matches + vt.right_votes;
        if (nd->right_child == 0) push_result(results, node_id, 'R', m2);
        else walk_tree(h, t, nd->right_child, m2, results);
    }
}

clust_result_t pick_best(result_list_t *rl)
{
    clust_result_t best = rl->a[0];
    for (int i = 1; i < rl->n; ++i)
        if (rl->a[i].matches > best.matches)
            best = rl->a[i];
    return best;
}

chr_result_t *classify_cenhaps(pg_msht_t *h, chr_forest_t *forest, int *n_out)
{
    chr_result_t *out;
    CALLOC(out, forest->n);

    for (int c = 0; c < forest->n; ++c) {
        clust_tree_t *t = &forest->a[c].tree;
        result_list_t results = {0};
        walk_tree(h, t, 1, 0, &results);

        out[c].chr = forest->a[c].chr; // borrowed, freed with forest
        out[c].best = results.n > 0 ? pick_best(&results) : (clust_result_t){0, '?', 0};
        free(results.a);
    }

    *n_out = forest->n;
    return out;
}