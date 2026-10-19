#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "clust.h"
#include "khashl.h"
#include "shtab.h"
#include "utils.h"

// ---- SNP-mer hashing / lookup ----

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

// returns -1 if the marker could not be hashed, 0 if absent from the table, 1 if found
int snpmer_lookup(pg_msht_t *h, marker_t *mk, uint32_t *cnt1, uint32_t *cnt2)
{
    if (!mk->h_snpmer) { *cnt1 = 0; *cnt2 = 0; return -1; }
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

// ---- reference panel loading ----
// One or more concatenated blocks, one per chromosome. A block starts with a
// header line "snpmer<TAB>sample1<TAB>sample2...<TAB>chr", followed by one row per SNP-mer:
// mean expected composition count1/(count1+count2) of each reference sample, then the chromosome.
// Empty, "NA", "nan" or "NaN" = missing (stored as a flag, never as a NaN value).

static char *next_field(char **s)
{
    char *f = *s;
    if (!f) return NULL;
    char *t = strchr(f, '\t');
    if (t) { *t = '\0'; *s = t + 1; } else *s = NULL;
    return f;
}

// returns 1 if the field is missing, else 0 and sets *v
static int parse_val(const char *f, float *v)
{
    *v = 0.0f;
    if (!*f || strcmp(f, "NA") == 0 || strcmp(f, "nan") == 0 || strcmp(f, "NaN") == 0)
        return 1;
    char *end;
    double d = strtod(f, &end);
    if (end == f) return 1;
    *v = (float)d;
    return 0;
}

static void geno_free(geno_mat_t *gm)
{
    if (!gm) return;
    for (int j = 0; j < gm->n_samp; ++j) free(gm->samp_names[j]);
    for (int g = 0; g < gm->n_geno; ++g) free(gm->geno_names[g]);
    for (int i = 0; i < gm->n_snp; ++i) {
        free(gm->snpmer[i]);
        free(gm->h_snpmer[i]);
    }
    free(gm->chr);
    free(gm->samp_names);
    free(gm->samp_geno);
    free(gm->geno_names);
    free(gm->snpmer);
    free(gm->h_snpmer);
    free(gm->p);
    free(gm->miss);
    free(gm);
}

// parse a header line (modified in place); returns NULL if it has no sample columns
static geno_mat_t *geno_open(char *line)
{
    geno_mat_t *gm;
    CALLOC(gm, 1);

    char *s = line;
    next_field(&s); // skip "snpmer"
    int m = 0;
    char *f;
    while (s && (f = next_field(&s))) {
        if (strcmp(f, "chr") == 0) continue; // metadata column, not a sample
        if (gm->n_samp == m) {
            m = m < 8 ? 8 : m + (m >> 1);
            REALLOC(gm->samp_names, m);
        }
        gm->samp_names[gm->n_samp++] = strdup(f);
    }
    if (gm->n_samp == 0) { geno_free(gm); return NULL; }
    return gm;
}

// Sample -> genotype label table from a "sample<TAB>genotype" TSV (header skipped).
typedef struct { char **sample, **geno; int n; } labels_t;

static void labels_free(labels_t *lb)
{
    for (int i = 0; i < lb->n; ++i) { free(lb->sample[i]); free(lb->geno[i]); }
    free(lb->sample); free(lb->geno);
    lb->sample = lb->geno = NULL; lb->n = 0;
}

// returns 0 on success, -1 if the file cannot be read
static int read_labels(const char *fn, labels_t *lb)
{
    lb->sample = lb->geno = NULL; lb->n = 0;
    FILE *fp = fopen(fn, "r");
    if (!fp) return -1;
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    int m = 0, first = 1;
    while ((len = getline(&line, &cap, fp)) >= 0) {
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
        if (len == 0) continue;
        if (first) { first = 0; if (strncmp(line, "sample\t", 7) == 0) continue; }
        char *t = strchr(line, '\t');
        if (!t) continue;
        *t = '\0';
        if (lb->n == m) {
            m = m < 64 ? 64 : m + (m >> 1);
            REALLOC(lb->sample, m);
            REALLOC(lb->geno, m);
        }
        lb->sample[lb->n] = strdup(line);
        lb->geno[lb->n]   = strdup(t + 1);
        lb->n++;
    }
    free(line);
    fclose(fp);
    return 0;
}

static char *g_exclude = NULL;   // reference sample to ignore (leave-one-out), see clust_set_exclude()

void clust_set_exclude(const char *name)
{
    free(g_exclude);
    g_exclude = name ? strdup(name) : NULL;
}

// Attach a genotype label to every matrix column. Columns without a label (or labelled "?")
// get samp_geno = -1 and are ignored when classifying. Returns the number of labelled samples.
static int assign_labels(geno_mat_t *gm, const labels_t *lb)
{
    int n_lab = 0, m = 0;
    MALLOC(gm->samp_geno, gm->n_samp);
    for (int j = 0; j < gm->n_samp; ++j) {
        gm->samp_geno[j] = -1;
        if (g_exclude && strcmp(g_exclude, gm->samp_names[j]) == 0) {
            fprintf(stderr, "[M::%s] %s: leaving reference sample '%s' out\n", __func__, gm->chr ? gm->chr : "NA", g_exclude);
            continue;
        }
        const char *g = NULL;
        for (int k = 0; k < lb->n; ++k)
            if (strcmp(lb->sample[k], gm->samp_names[j]) == 0) { g = lb->geno[k]; break; }
        if (!g || !*g || strcmp(g, "?") == 0 || strcmp(g, "NA") == 0) continue;

        int gi = -1;
        for (int x = 0; x < gm->n_geno; ++x)
            if (strcmp(gm->geno_names[x], g) == 0) { gi = x; break; }
        if (gi < 0) {
            if (gm->n_geno == m) {
                m = m < 8 ? 8 : m + (m >> 1);
                REALLOC(gm->geno_names, m);
            }
            gm->geno_names[gm->n_geno] = strdup(g);
            gi = gm->n_geno++;
        }
        gm->samp_geno[j] = gi;
        n_lab++;
    }
    return n_lab;
}

// finish a block (labels and logs) and append it to the set; frees it if unusable
static void geno_push(geno_mat_t ***set, int *n, int *m, geno_mat_t *gm, int n_bad_hash, const labels_t *lb)
{
    if (!gm) return;
    if (gm->n_snp == 0) {
        fprintf(stderr, "[W::%s] block '%s' has no SNP-mer rows, skipping\n",
                __func__, gm->chr ? gm->chr : "NA");
        geno_free(gm);
        return;
    }
    if (!gm->chr) gm->chr = strdup("NA");

    int n_lab = assign_labels(gm, lb);
    if (n_lab == 0) {
        fprintf(stderr, "[W::%s] %s: none of the %d matrix samples has a genotype label, skipping\n",
                __func__, gm->chr, gm->n_samp);
        geno_free(gm);
        return;
    }
    if (n_lab < gm->n_samp)
        fprintf(stderr, "[W::%s] %s: %d/%d matrix samples have no genotype label and are ignored\n",
                __func__, gm->chr, gm->n_samp - n_lab, gm->n_samp);

    fprintf(stderr, "[M::%s] %s: %d SNP-mers x %d samples (%d genotypes, %d unhashable SNP-mers)\n",
            __func__, gm->chr, gm->n_snp, n_lab, gm->n_geno, n_bad_hash);

    if (*n == *m) {
        *m = *m < 8 ? 8 : *m + (*m >> 1);
        REALLOC(*set, *m);
    }
    (*set)[(*n)++] = gm;
}

geno_mat_t **load_geno_matrices(pg_msht_t *h, const char *fn, const char *labels_fn, int *n_out)
{
    *n_out = 0;
    labels_t lb;
    if (!labels_fn || read_labels(labels_fn, &lb) < 0) {
        fprintf(stderr, "[E::%s] failed to read genotype labels '%s'\n", __func__, labels_fn ? labels_fn : "(none)");
        return NULL;
    }
    FILE *fp = fopen(fn, "r");
    if (!fp) {
        fprintf(stderr, "[E::%s] failed to open '%s'\n", __func__, fn);
        labels_free(&lb);
        return NULL;
    }

    geno_mat_t **set = NULL;
    int n = 0, m = 0;
    geno_mat_t *gm = NULL;
    int m_snp = 0, n_bad_hash = 0;

    char *line = NULL;
    size_t cap = 0;
    ssize_t len;

    while ((len = getline(&line, &cap, fp)) >= 0) {
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
        if (len == 0) continue;

        // a header line starts a new chromosome block
        if (strncmp(line, "snpmer\t", 7) == 0) {
            geno_push(&set, &n, &m, gm, n_bad_hash, &lb);
            gm = geno_open(line);
            m_snp = 0;
            n_bad_hash = 0;
            if (!gm) fprintf(stderr, "[W::%s] header without sample columns, block skipped\n", __func__);
            continue;
        }
        if (!gm) continue; // rows of a skipped block

        int S = gm->n_samp;
        if (gm->n_snp == m_snp) {
            m_snp = m_snp < 64 ? 64 : m_snp + (m_snp >> 1);
            REALLOC(gm->snpmer, m_snp);
            REALLOC(gm->h_snpmer, m_snp);
            REALLOC(gm->p, (size_t)m_snp * S);
            REALLOC(gm->miss, (size_t)m_snp * S);
        }
        int i = gm->n_snp;
        char *s = line;
        char *name = next_field(&s);

        int j = 0;
        char *f;
        // read exactly S values; the field after them is the chr column
        while (s && j < S && (f = next_field(&s))) {
            float v;
            size_t x = (size_t)i * S + j++;
            gm->miss[x] = (uint8_t)parse_val(f, &v);
            gm->p[x] = v;
        }
        if (j != S) {
            fprintf(stderr, "[W::%s] row '%s' has %d/%d values, skipping\n", __func__, name, j, S);
            continue;
        }
        if (!gm->chr && s) {
            char *c = next_field(&s);
            if (c && *c) gm->chr = strdup(c);
        }
        gm->snpmer[i] = strdup(name);
        gm->h_snpmer[i] = precompute_snpmer_hash(h, name);
        if (!gm->h_snpmer[i]) ++n_bad_hash;
        gm->n_snp++;
    }
    geno_push(&set, &n, &m, gm, n_bad_hash, &lb);
    fclose(fp);
    free(line);
    labels_free(&lb);

    if (n == 0) {
        fprintf(stderr, "[E::%s] no usable matrix blocks in '%s'\n", __func__, fn);
        free(set);
        return NULL;
    }
    *n_out = n;
    return set;
}

// Read a text file with one path per line into an array of strings
static char **read_path_list(const char *fn, int *n_out)
{
    *n_out = 0;
    FILE *fp = fopen(fn, "r");
    if (!fp) return NULL;
    char **v = NULL;
    int n = 0, m = 0;
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    while ((len = getline(&line, &cap, fp)) >= 0) {
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' ||
                           line[len-1] == ' '  || line[len-1] == '\t'))
            line[--len] = '\0';
        if (len == 0 || line[0] == '#') continue;
        if (n == m) { m = m < 8 ? 8 : m + (m >> 1); REALLOC(v, m); }
        v[n++] = strdup(line);
    }
    free(line);
    fclose(fp);
    *n_out = n;
    return v;
}

// Read a text file with one matrix path per line, and a text file with one labels path per
// line (same order), and load every chromosome block found into a single array.
geno_mat_t **load_matrices_from_list(pg_msht_t *h, const char *list_fn, const char *labels_list_fn, int *n_out)
{
    *n_out = 0;
    int n_mf = 0, n_lf = 0;
    char **mat_fns = read_path_list(list_fn, &n_mf);
    if (!mat_fns) {
        fprintf(stderr, "[E::%s] failed to open '%s'\n", __func__, list_fn);
        return NULL;
    }
    char **lab_fns = read_path_list(labels_list_fn, &n_lf);
    if (!lab_fns) {
        fprintf(stderr, "[E::%s] failed to open '%s'\n", __func__, labels_list_fn);
        for (int i = 0; i < n_mf; ++i) free(mat_fns[i]);
        free(mat_fns);
        return NULL;
    }
    if (n_lf != n_mf)
        fprintf(stderr, "[W::%s] %d matrix paths but %d label paths; matrices without a label file are skipped\n",
                __func__, n_mf, n_lf);

    geno_mat_t **all = NULL;
    int n = 0, m = 0;
    for (int f = 0; f < n_mf; ++f) {
        if (f >= n_lf) {
            fprintf(stderr, "[W::%s] no labels file for '%s', skipping\n", __func__, mat_fns[f]);
            continue;
        }
        int k = 0;
        geno_mat_t **gms = load_geno_matrices(h, mat_fns[f], lab_fns[f], &k);
        if (!gms) {
            fprintf(stderr, "[W::%s] no matrices loaded from '%s', skipping\n", __func__, mat_fns[f]);
            continue;
        }
        if (n + k > m) {
            m = n + k > 2 * m ? n + k : 2 * m;
            REALLOC(all, m);
        }
        memcpy(all + n, gms, (size_t)k * sizeof(*gms));
        n += k;
        free(gms); // frees only the array; the matrices now belong to 'all'
    }
    for (int i = 0; i < n_mf; ++i) free(mat_fns[i]);
    for (int i = 0; i < n_lf; ++i) free(lab_fns[i]);
    free(mat_fns); free(lab_fns);

    if (n == 0) {
        fprintf(stderr, "[E::%s] no matrices loaded from any file listed in '%s'\n", __func__, list_fn);
        free(all);
        return NULL;
    }
    *n_out = n;
    return all;
}

void destroy_geno_matrices(geno_mat_t **gms, int n)
{
    if (!gms) return;
    for (int i = 0; i < n; ++i) geno_free(gms[i]);
    free(gms);
}

// ---- classification ----
// Nearest-neighbour call against the reference samples. For every SNP-mer i of the panel the query has
// either an observed composition q_i = c1_i / (c1_i + c2_i) (reads left after the noise filter) or nothing
// (both counts 0); reference sample j has a mean composition p_ij or nothing (NA). The error of SNP-mer i is
//   both present   : (q_i - p_ij)^2
//   both absent    : 0                         (agreeing absence is a perfect match)
//   query only     : ABS_PEN_QUERY             (reads show it, the reference lacks it)
//   reference only : ABS_PEN_REF               (the reference has it, the reads show nothing)
// and the distance is the root mean error over ALL the SNP-mers of the panel:
//   d_j = sqrt( (1/S) * sum_i error_ij )
// The call is the genotype label of the sample with the smallest d_j. The denominator is the same for every
// reference, so distances are comparable; a reference that shares no present SNP-mer with the query is skipped.

void classify_genotype(pg_msht_t *h, const geno_mat_t *gm, geno_call_t *out)
{
    int S = gm->n_snp, N = gm->n_samp;
    out->nn = out->geno = out->nn2 = -1;
    out->dist_best = out->dist_second = INFINITY; // comparison start values only
    out->n_shared = 0;
    out->n_obs = out->n_missing = out->n_used = 0;
    out->depth = 0;

    double *q;
    double *dist;
    int *shared;
    MALLOC(q, S);
    MALLOC(dist, N);
    MALLOC(shared, N);

    // 1. Observed composition per SNP-mer; NaN-free: SNP-mers with both counts 0 are marked with q = -1
    for (int i = 0; i < S; ++i) {
        marker_t mk;
        memset(&mk, 0, sizeof(mk));
        mk.h_snpmer = gm->h_snpmer[i];
        uint32_t a = 0, b = 0;
        if (snpmer_lookup(h, &mk, &a, &b) < 0) { a = b = 0; } // unhashable -> missing
        if ((int)a <= NOISE_THRESHOLD) a = 0;
        if ((int)b <= NOISE_THRESHOLD) b = 0;
        if (a == 0 && b == 0) { q[i] = -1.0; out->n_missing++; }
        else {
            q[i] = (double)a / ((double)a + (double)b);
            out->n_obs++;
            out->depth += (uint64_t)a + b;
        }
    }
    out->n_used = out->n_obs;

    if (out->n_obs == 0) { // no data: nothing to compare
        free(q); free(dist); free(shared);
        return;
    }

    // 2. Distance to every labelled reference sample; keep the nearest
    for (int j = 0; j < N; ++j) {
        dist[j] = INFINITY;
        shared[j] = 0;
        if (gm->samp_geno[j] < 0) continue;
        double sum = 0.0;
        int n = 0;                                    // SNP-mers present in both
        for (int i = 0; i < S; ++i) {
            int qp = q[i] >= 0.0;                     // observed in the query
            int rp = !gm->miss[(size_t)i * N + j];    // has a profile in this reference sample
            if (qp && rp) {
                double d = q[i] - (double)gm->p[(size_t)i * N + j];
                sum += d * d;
                ++n;
            } else if (qp) sum += ABS_PEN_QUERY;
            else if (rp) sum += ABS_PEN_REF;
            // both absent: no error
        }
        shared[j] = n;
        if (n == 0) continue;                         // nothing to compare in common
        dist[j] = sqrt(sum / S);
        if (dist[j] < out->dist_best) { out->dist_best = dist[j]; out->nn = j; out->n_shared = n; }
    }
    if (out->nn < 0) { free(q); free(dist); free(shared); return; }
    out->geno = gm->samp_geno[out->nn];

    // 3. Nearest sample with a different genotype, to measure how clear the call is
    for (int j = 0; j < N; ++j) {
        if (gm->samp_geno[j] < 0 || gm->samp_geno[j] == out->geno) continue;
        if (dist[j] < out->dist_second) { out->dist_second = dist[j]; out->nn2 = j; }
    }

    free(q);
    free(dist);
    free(shared);
}

// Output: chr, genotype, nn_sample. One line per chromosome.
// genotype = label of the nearest reference sample (e.g. "1/2", same labels as the
// genotype column of the labels file); nn_sample = that reference sample.
// No call writes NA.
void print_geno_call(FILE *out, const geno_mat_t *gm, const geno_call_t *c)
{
    if (c->nn < 0) {
        fprintf(out, "%s\tNA\tNA\n", gm->chr);
        return;
    }
    fprintf(out, "%s\t%s\t%s\n", gm->chr, gm->geno_names[c->geno], gm->samp_names[c->nn]);
}