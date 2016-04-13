
cdef extern from "htslib/faidx.h" nogil:

    ctypedef struct faidx_t:
       pass

    int fai_build(char *fn)

    void fai_destroy(faidx_t *fai)

    faidx_t *fai_load(char *fn)

    char *fai_fetch(faidx_t *fai,
                    char *reg,
                    int *len)

    int faidx_nseq(faidx_t *fai)

    int faidx_has_seq(faidx_t *fai, const char *seq)

    char *faidx_fetch_seq(faidx_t *fai,
                         char *c_name,
                         int p_beg_i,
                         int p_end_i,
                         int *len)

    int faidx_seq_len(faidx_t *fai, const char *seq)
