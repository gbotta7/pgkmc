#include "shtab.h"
#include "utils.h"

unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

char nt4_seq_table[5] = {
    'A', // 0
    'C', // 1
    'G', // 2
    'T', // 3
    'N'  // 4 (invalid/ambiguous)
};


void pg_opt_init(pg_opt_t *o)
{
	memset(o, 0, sizeof(pg_opt_t));
	o->k = 31;
	o->msf = 0.95;  // total allelic frequency
    o->maf = 0;     // minimum allelic frequency
    o->snp = 0;
	o->pre = 10;
	o->filt_type = 0;
	o->n_threads = 3;
	o->chunk_size = mm_parse_num("1.9g");
    o->write_info = 0;
    o->write_mko = S_COUNTER_MAX;
	o->verbose = 0;
}

int64_t mm_parse_num(const char *str)
{
    double x;
    char *p;
    x = strtod(str, &p);
    if      (*p == 'G' || *p == 'g') x *= 1e9;
    else if (*p == 'M' || *p == 'm') x *= 1e6;
    else if (*p == 'K' || *p == 'k') x *= 1e3;
    return (int64_t)(x + .499);
}

void sample_name_from_path(const char *gnm_fn, char *out, size_t out_sz)
{
    const char *bname = strrchr(gnm_fn, '/');
    bname = bname ? bname + 1 : gnm_fn;

    strncpy(out, bname, out_sz - 1);
    out[out_sz - 1] = '\0';

    // strip .gz if present
    size_t len = strlen(out);
    if (len > 3 && strcmp(out + len - 3, ".gz") == 0)
        out[len - 3] = '\0';

    // strip .fa, .fna, or .fasta
    static const char *fa_exts[] = { ".fasta", ".fna", ".fa", NULL };
    for (int i = 0; fa_exts[i]; i++) {
        len = strlen(out);
        size_t elen = strlen(fa_exts[i]);
        if (len > elen && strcmp(out + len - elen, fa_exts[i]) == 0) {
            out[len - elen] = '\0';
            break;
        }
    }
}