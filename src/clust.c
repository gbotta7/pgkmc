#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "clust.h"
#include "khashl.h"
#include "shtab.h"
#include "utils.h"

uint64_t *precompute_snpmer_hash(pg_msht_t *h, const char *snpmer)
{
    int half = h->k >> 1;
    if ((int)strlen(snpmer) != h->k + 4) return NULL;

    char left[64], right[64];
    memcpy(left, snpmer, half); left[half] = '\0';
    memcpy(right, snpmer + half + 5, half); right[half] = '\0';

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
int snpmer_lookup(pg_msht_t *h, const uint64_t *h_snpmer, uint32_t *cnt1, uint32_t *cnt2)
{
    if (!h_snpmer) { *cnt1 = 0; *cnt2 = 0; return -1; }
    int bucket = *h_snpmer & ((1 << h->pre) - 1);
    uint64_t key = (*h_snpmer >> h->pre) << F_VAL_INFO_BITS;

    pg_sht1_t *g = &h->h[bucket];
    khint_t k = pg_sht_get(g->h, key);
    if (k == kh_end(g->h)) { *cnt1 = 0; *cnt2 = 0; return 0; }

    uint32_t v = kh_val(g->h, k);
    *cnt1 = s_val_count1(v);
    *cnt2 = s_val_count2(v);
    return 1;
}

void classify_query(const geno_mat_t *gm, const double *q, int k, geno_call_t *call)
{
    int n_snps = gm->n_snps, n_haps = gm->n_haps, n_geno = gm->n_geno, n_cls = gm->n_cls;
    call->nn_a = call->nn_b = call->gt_idx = -1;
    call->dist_best = INFINITY;
    call->dist_vote = INFINITY;
    call->k_used = call->votes = call->tie = 0;
    call->n_pairs = 0;
 
    int n_obs = 0; // SNP-mers observed in the query
    for (int i = 0; i < n_snps; ++i)
        if (q[i] >= 0.0) n_obs++;
    if (n_obs == 0 || n_geno == 0) return;
 
    if (k < 1) k = 1;
    double *td; int *ta, *tb, *tg; // save top-k sorted distance, hapA, hapB, and genotype
    MALLOC(td, k); MALLOC(ta, k); MALLOC(tb, k); MALLOC(tg, k);
    int nn = 0; // number of NN stored
 
    for (int a = 0; a < n_haps; ++a) {
        if (gm->hap_cls[a] < 0) continue; // continue if there are no labels for this haplotype
        for (int b = a + 1; b < n_haps; ++b) {
            if (gm->hap_cls[b] < 0) continue; // continue if there are no labels for this haplotype
            int ca = gm->hap_cls[a], cb = gm->hap_cls[b];
            int g = gm->gt_map[(size_t)ca * n_cls + cb]; // index of the genotype
 
            double sum = 0.0;
            int n = 0;
            for (int i = 0; i < n_snps; ++i) {
                int pa = !gm->miss[(size_t)i * n_haps + a], pb = !gm->miss[(size_t)i * n_haps + b]; // check presence in both haplotypes
                int qp = q[i] >= 0.0; // query observed
                int rp = pa || pb;    // pair has a profile
                double v = 0.0;
                if (pa && pb) v = 0.5 * ((double)gm->p[(size_t)i * n_haps + a] + (double)gm->p[(size_t)i * n_haps + b]); // average if present in both haplotypes
                else if (pa) v = gm->p[(size_t)i * n_haps + a]; // only hapA
                else if (pb) v = gm->p[(size_t)i * n_haps + b]; // only hapB
                if (qp && rp) {
                    double d = q[i] - v;
                    sum += d * d; // accumulate squared error
                    ++n;
                } else if (qp ^ rp) sum += qp ? ABS_PEN_QUERY : ABS_PEN_REF; // only one of them has data
            }
            if (n == 0) continue; // skip pairs that share no present SNP-mer with the query
 
            call->n_pairs++;
            double d = sqrt(sum / n_snps); // RMSE
 
            if (nn < k || d < td[nn-1]) { // add until kNN or once you find a closer one, equal distances keep scan order
                int pos = nn < k ? nn : k - 1;
                while (pos > 0 && td[pos-1] > d) { // find where to insert the new NN
                    td[pos] = td[pos-1]; ta[pos] = ta[pos-1]; tb[pos] = tb[pos-1]; tg[pos] = tg[pos-1];
                    --pos;
                }
                td[pos] = d; ta[pos] = a; tb[pos] = b; tg[pos] = g;
                if (nn < k) ++nn;
            }
        }
    }
 
    if (nn > 0) {
        int *votes;
        double *sumd;
        CALLOC(votes, n_geno);
        CALLOC(sumd, n_geno);
        for (int r = 0; r < nn; ++r) { // votes per genotype and sum of the distances of its voters
            votes[tg[r]]++;
            sumd[tg[r]] += td[r];
        }
        int best = -1;
        for (int g = 0; g < n_geno; ++g) {
            if (votes[g] == 0) continue;
            if (best < 0 || votes[g] > votes[best] || (votes[g] == votes[best] && sumd[g] < sumd[best])) best = g; // update best genotype
        }
        for (int g = 0; g < n_geno; ++g) { // just check if there are vote ties
            if (g != best && votes[g] == votes[best]) {
                call->tie = 1;
                break;
            }
        }
        call->gt_idx = best;
        call->k_used = nn;
        call->votes = votes[best];
        call->dist_vote = sumd[best] / votes[best]; // avg distance across the votes for that genotype
 
        for (int r = 0; r < nn; ++r) // nearest pair of the winning genotype (always a voter)
            if (tg[r] == best) {
                call->nn_a = ta[r];
                call->nn_b = tb[r];
                call->dist_best = td[r];
                break;
            }
 
        free(votes); free(sumd);
    }
    free(td); free(ta); free(tb); free(tg);
}

void classify_genotype(pg_msht_t *h, const geno_mat_t *gm, int k, geno_call_t *call, int use_asm)
{
    int n_snps = gm->n_snps;
    double *q;
    MALLOC(q, n_snps);

    int noise_threshold = 3;
    if (use_asm) noise_threshold = 0;

    // observed composition per SNP-mer; SNP-mers with both counts 0 are marked with q = -1
    for (int i = 0; i < n_snps; ++i) {
        uint32_t a = 0, b = 0;
        snpmer_lookup(h, gm->h_snpmer[i], &a, &b);
        if ((int)a <= noise_threshold) a = 0;
        if ((int)b <= noise_threshold) b = 0;
        if (a == 0 && b == 0) q[i] = -1.0;
        else q[i] = (double)a / ((double)a + (double)b);
    }
    classify_query(gm, q, k, call);
    free(q);
}

// Output: chr, genotype, nn_pair, votes, dist_vote
void print_geno_call(FILE *out, const geno_mat_t *gm, const geno_call_t *c)
{
    if (c->nn_a < 0) {
        fprintf(out, "%s\tNA\tNA\tNA\tNA\n", gm->chr);
        return;
    }
    fprintf(out, "%s\t%s\t%s+%s\t%d\t%.4f\n", gm->chr, gm->geno_names[c->gt_idx], gm->hap_names[c->nn_a], gm->hap_names[c->nn_b], c->votes, c->dist_vote);
}
