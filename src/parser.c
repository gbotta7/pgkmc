#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parser.h"
#include "kseq.h"
#include "utils.h"
#include "clust.h"

KSEQ_INIT(gzFile, gzread)

// FASTA/FASTQ parser
struct bseq_file_s {
	gzFile fp;
	kseq_t *ks;
};

bseq_file_t *bseq_open(const char *fn)
{
	bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (f == 0) return 0;
	fp = calloc(1, sizeof(bseq_file_t));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
	return fp;
}

void bseq_close(bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

bseq1_t *bseq_read(bseq_file_t *fp, int64_t chunk_size, int keep_comment, int *n_)
{
	int m, n;
	int64_t size = 0;
	bseq1_t *seqs;
	kseq_t *ks = fp->ks;
	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0) {
		bseq1_t *s;
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = realloc(seqs, m * sizeof(bseq1_t));
		}
		s = &seqs[n];
		s->name = strdup(ks->name.s);
		s->comment = (ks->comment.s && keep_comment)? strdup(ks->comment.s) : 0;
		s->seq = strdup(ks->seq.s);
		s->qual = ks->qual.l? strdup(ks->qual.s) : 0;
		s->l_seq = ks->seq.l;
		s->aux = 0;
		size += seqs[n++].l_seq;
		if (size >= chunk_size) break;
	}
	*n_ = n;
	return seqs;
}


// ---- LOAD COMPOSITION MATRIX ----
// pulls the first column of the row (the SNP-mer) out of the line and advances the cursor past it
static char *next_field(char **s)
{
    char *f = *s;
    if (!f) return NULL;
    char *t = strchr(f, '\t');
    if (t) {
        *t = '\0';
        *s = t + 1;
    } else *s = NULL;
    return f;
}

// parse a header line (modified in place) into a list of sample names and return the number of samples
static int parse_header(char *line, char ***names_out)
{
    char **names = NULL;
    int n = 0, m = 0;
    char *s = line;
    char *f;
    next_field(&s); // skip "snpmer"
    while (s && (f = next_field(&s))) {
        if (strcmp(f, "chr") == 0) continue;
        if (n == m) {
            m = m < 8 ? 8 : m + (m >> 1);
            REALLOC(names, m);
        }
        names[n++] = strdup(f);
    }
    *names_out = names; // address of a list of strings
    return n;
}

static void names_free(char **names, int n)
{
    for (int i = 0; i < n; ++i) free(names[i]);
    free(names);
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

// find the block of each chromosome, creating it on first sight
static blk_t *blk_get(blk_t **v, int *n, int *m, const char *chr, int n_haps, char **names)
{
    for (int i = *n - 1; i >= 0; --i)  // the most recent block is usually the right one
        if (strcmp((*v)[i].gm->chr, chr) == 0) return &(*v)[i];
    if (*n == *m) {
        *m = *m < 8 ? 8 : *m + (*m >> 1);
        REALLOC(*v, *m);
    }
    blk_t *b = &(*v)[(*n)++];
    memset(b, 0, sizeof(*b));
    CALLOC(b->gm, 1);
    b->gm->chr = strdup(chr); // chromosome name
    b->gm->n_haps = n_haps; // number of haplotype columns
    MALLOC(b->gm->hap_names, n_haps);
    for (int j = 0; j < n_haps; ++j) {
        b->gm->hap_names[j] = strdup(names[j]); // haplotype names
    }
    return b;
}

static int geno_lookup(geno_mat_t *gm, const char *name)
{
    for (int x = 0; x < gm->n_geno; ++x)
        if (strcmp(gm->geno_names[x], name) == 0) return x;
    REALLOC(gm->geno_names, gm->n_geno + 1);
    gm->geno_names[gm->n_geno] = strdup(name);
    return gm->n_geno++;
}

static int to_int(const char *s, long *v)
{
    char *end;
    *v = strtol(s, &end, 10);
    return end != s && *end == '\0';
}

// "a/b" with a <= b (numeric order when both are integers)
static void pair_name(const char *x, const char *y, char *buf, size_t cap)
{
    long a, b;
    int swap = (to_int(x, &a) && to_int(y, &b)) ? a > b : strcmp(x, y) > 0;
    snprintf(buf, cap, "%s/%s", swap ? y : x, swap ? x : y);
}

static void build_pairs(geno_mat_t *gm)
{
    int n_cls = gm->n_cls;
    gm->gt_map = NULL;
    gm->geno_names = NULL;
    gm->n_geno = 0;
    if (n_cls == 0) return;

    MALLOC(gm->gt_map, (size_t)n_cls * n_cls); // all combinations (e.g. 1/2 and 2/1)
    for (int a = 0; a < n_cls; ++a)
        for (int b = a; b < n_cls; ++b) {
            char nm[128];
            pair_name(gm->cls_name[a], gm->cls_name[b], nm, sizeof(nm)); // merge haplotypes and sort them numerically
            int g = geno_lookup(gm, nm);
            gm->gt_map[(size_t)a * n_cls + b] = g;
            gm->gt_map[(size_t)b * n_cls + a] = g;   // symmetric
        }
}

static int find_label(const labels_t *lb, const int *cand, int nc, const char *name)
{
    for (int c = 0; c < nc; ++c)
        if (strcmp(lb->hap[cand[c]], name) == 0) return cand[c]; // walk indices of the candidate rows, return the right index for that haplotype
    return -1;
}

static int assign_labels(geno_mat_t *gm, const labels_t *lb)
{
    int n_lab = 0, m = 0, nc = 0;
    int *cand_chrom;
    MALLOC(cand_chrom, lb->n);
    for (int k = 0; k < lb->n; ++k) // mask labels to the chromosome of this block (or all chromosomes if no chr column is present)
        if (!lb->chr[k] || strcmp(lb->chr[k], gm->chr) == 0) cand_chrom[nc++] = k;

    MALLOC(gm->hap_cls, gm->n_haps);
    for (int j = 0; j < gm->n_haps; ++j) {
        gm->hap_cls[j] = -1;
        int hap_idx = find_label(lb, cand_chrom, nc, gm->hap_names[j]);
        const char *hap_cluster = hap_idx >= 0 ? lb->geno[hap_idx] : NULL;
        if (!hap_cluster || strcmp(hap_cluster, "?") == 0 || strcmp(hap_cluster, "NA") == 0) continue; // ignore unlabelled haplotypes or strcmp crashes if g is NULL

        hap_idx = -1;
        for (int x = 0; x < gm->n_cls; ++x)
            if (strcmp(gm->cls_name[x], hap_cluster) == 0) { // scan against clusters already recorded for this chromosome
                hap_idx = x;
                break;
            }
        if (hap_idx < 0) {
            if (gm->n_cls == m) { // grow the list of clusters
                m = m < 8 ? 8 : m + (m >> 1);
                REALLOC(gm->cls_name, m);
            }
            gm->cls_name[gm->n_cls] = strdup(hap_cluster);
            hap_idx = gm->n_cls++;
        }

        gm->hap_cls[j] = hap_idx;
        n_lab++;
    }

    free(cand_chrom);
    return n_lab; // number of labelled haplotypes
}

static void geno_free(geno_mat_t *gm)
{
    if (!gm) return;
    for (int j = 0; j < gm->n_haps; ++j) free(gm->hap_names[j]);
    free(gm->hap_names);
    for (int x = 0; x < gm->n_cls; ++x) free(gm->cls_name[x]);
    free(gm->cls_name);
    free(gm->hap_cls);
    for (int g = 0; g < gm->n_geno; ++g) free(gm->geno_names[g]);
    free(gm->geno_names);
    free(gm->gt_map);
    for (int i = 0; i < gm->n_snps; ++i) {
        free(gm->snpmer[i]);
        free(gm->h_snpmer[i]);
    }
    free(gm->snpmer);
    free(gm->h_snpmer);
    free(gm->p);
    free(gm->miss);
    free(gm->chr);
    free(gm);
}

void destroy_geno_matrices(geno_mat_t **gms, int n)
{
    if (!gms) return;
    for (int i = 0; i < n; ++i) geno_free(gms[i]);
    free(gms);
}

static void geno_push(geno_mat_t ***set, int *n, int *m, blk_t *b, const labels_t *lb)
{
    geno_mat_t *gm = b->gm;
    b->gm = NULL;
    if (gm->n_snps == 0) {
        fprintf(stderr, "[W::%s] %s has no SNP-mer rows, skipping\n", __func__, gm->chr);
        geno_free(gm);
        return;
    }

    int n_lab = assign_labels(gm, lb);
    if (n_lab == 0) {
        fprintf(stderr, "[W::%s] %s: none of the %d matrix samples has a haplotype label, skipping\n", __func__, gm->chr, gm->n_haps);
        geno_free(gm);
        return;
    }
    if (n_lab < gm->n_haps)
        fprintf(stderr, "[W::%s] %s: %d/%d matrix samples have no haplotype label and are ignored\n", __func__, gm->chr, gm->n_haps - n_lab, gm->n_haps);

    build_pairs(gm);
    fprintf(stderr, "[M::%s] %s: %d SNP-mers x %d haplotypes (%d clusters, %d possible genotypes)\n", __func__, gm->chr, gm->n_snps, n_lab, gm->n_cls, gm->n_geno);

    if (*n == *m) {
        *m = *m < 8 ? 8 : *m + (*m >> 1);
        REALLOC(*set, *m);
    }
    (*set)[(*n)++] = gm;
}

static void labels_free(labels_t *lb)
{
    for (int i = 0; i < lb->n; ++i) { free(lb->hap[i]); free(lb->geno[i]); free(lb->chr[i]); }
    free(lb->hap); free(lb->geno); free(lb->chr);
    lb->hap = lb->geno = lb->chr = NULL; lb->n = 0;
}

static int read_labels(const char *fn, labels_t *lb)
{
    lb->hap = lb->geno = lb->chr = NULL; lb->n = 0;
    FILE *fp = fopen(fn, "r");
    if (!fp) return 1;
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    int m = 0, first = 1;
    while ((len = getline(&line, &cap, fp)) >= 0) {
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
        if (len == 0) continue;
        if (first) { first = 0; if (strncmp(line, "haplotype\t", 7) == 0) continue; }
        char *s = line;
        char *samp = next_field(&s);
        char *geno = next_field(&s);
        char *chr  = next_field(&s);   // NULL if the file has only two columns
        if (!samp || !geno) continue;
        if (lb->n == m) {
            m = m < 64 ? 64 : m + (m >> 1);
            REALLOC(lb->hap, m);
            REALLOC(lb->geno, m);
            REALLOC(lb->chr, m);
        }
        lb->hap[lb->n] = strdup(samp);
        lb->geno[lb->n] = strdup(geno);
        lb->chr[lb->n] = chr && *chr ? strdup(chr) : NULL; // if a chr column is present, a label only applies to the chromosome it names; without it, a label applies to every chromosome
        lb->n++;
    }

    free(line);
    fclose(fp);

    return 0;
}

// Loads the single matrix file and the single labels file
geno_mat_t **load_geno_matrix(pg_msht_t *h, const char *fn, const char *labels_fn, int *n_out)
{      
    *n_out = 0;

    // Read the labels file
    labels_t lb;
    if (!labels_fn || read_labels(labels_fn, &lb) > 0) {
        fprintf(stderr, "[E::%s] failed to read haplotype labels '%s'\n", __func__, labels_fn);
        return NULL;
    }

    // Read the matrix composition file
    FILE *fp = fopen(fn, "r");
    if (!fp) {
        fprintf(stderr, "[E::%s] failed to open matrix'%s'\n", __func__, fn);
        labels_free(&lb);
        return NULL;
    }

    blk_t *blks = NULL;
    int nb = 0, mb = 0;
    char **names = NULL;
    int n_haps = 0;
    float *vtmp = NULL;
    uint8_t *mtmp = NULL;

    char *line = NULL;
    size_t cap = 0;
    ssize_t len;

    while ((len = getline(&line, &cap, fp)) >= 0) {
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
        if (len == 0) continue;

        if (strncmp(line, "snpmer\t", 7) == 0) { // header: sample names
            n_haps = parse_header(line, &names);
            if (n_haps == 0) {
                fprintf(stderr, "[W::%s] Header without sample columns, rows skipped\n", __func__);
                return NULL; // no header
            }
            else {
                MALLOC(vtmp, n_haps);
                MALLOC(mtmp, n_haps);
            }
            continue;
        }

        char *s = line;
        char *name = next_field(&s); // name is SNP-mer
        int j = 0;
        char *f;
        while (s && j < n_haps && (f = next_field(&s))) { // exactly n_haps values
            mtmp[j] = (uint8_t)parse_val(f, &vtmp[j]); // vtmp holds the values, mtmp holds the mask of NA (1 where missing)
            ++j;
        }
        if (j != n_haps) {
            fprintf(stderr, "[W::%s] Row '%s' has values for %d/%d samples, skipping\n", __func__, name, j, n_haps);
            continue;
        }

        char *chr = s ? next_field(&s) : NULL;   // the field after the values is the chr column
        if (!chr || !*chr) chr = "NA";

        blk_t *b = blk_get(&blks, &nb, &mb, chr, n_haps, names); // create the geno_mat for this chromosome if it doesn't exist yet and return the block to store SNP-mer compositions
        geno_mat_t *gm = b->gm;
        if (gm->n_snps == b->m_snp) {
            b->m_snp = b->m_snp < 64 ? 64 : b->m_snp + (b->m_snp >> 1); // start at 64 rows, grow by 1.5x each time
            REALLOC(gm->snpmer, b->m_snp);
            REALLOC(gm->h_snpmer, b->m_snp);
            REALLOC(gm->p, (size_t)b->m_snp * n_haps);
            REALLOC(gm->miss, (size_t)b->m_snp * n_haps);
        }
        int i = gm->n_snps;
        memcpy(gm->p + (size_t)i * n_haps, vtmp, (size_t)n_haps * sizeof(float));
        memcpy(gm->miss + (size_t)i * n_haps, mtmp, (size_t)n_haps);
        gm->snpmer[i] = strdup(name);
        gm->h_snpmer[i] = precompute_snpmer_hash(h, name);
        gm->n_snps++;
    }

    fclose(fp);
    free(line);
    names_free(names, n_haps);
    free(vtmp); free(mtmp);

    // unwrap blocks into a single array of geno_mat_t, one per chromosome, store haplotypes clusters and pairs genotypes
    geno_mat_t **set = NULL;
    int n = 0, m = 0;
    for (int i = 0; i < nb; ++i) geno_push(&set, &n, &m, &blks[i], &lb);
    free(blks);
    labels_free(&lb);

    if (n == 0) {
        fprintf(stderr, "[E::%s] no usable chromosomes in '%s'\n", __func__, fn);
        free(set);
        return NULL;
    }

    *n_out = n;
    return set;
}