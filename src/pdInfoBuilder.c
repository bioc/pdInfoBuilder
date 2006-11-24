#include "pdInfoBuilder.h"

static const R_CallMethodDef R_CallDef[] =
{
    {"PIB_25mers_to_mat", (DL_FUNC)&PIB_25mers_to_mat, 1},
    {NULL, NULL, 0},
};

void R_init_pdInfoBuilder(DllInfo *info)
{
    R_registerRoutines(info, NULL, R_CallDef, NULL, NULL);
}

/*
 * PIB_25mers_to_mat
 *
 * Return a matrix with 50 columns coding the 25mers passed as input.
 *
 * mers - a STRSXP containing DNA sequences (A, T, C, G) with 25 chars
 *        each.
 *
 * Each letter in a given 25mer is coded as two columns in the
 * resulting matrix.  The codes are:
 *
 *     A  0 0
 *     T  0 1
 *     C  1 0
 *     G  1 1
 *
 */
SEXP PIB_25mers_to_mat(SEXP mers)
{
    int n;                      /* number of 25mers */
    int i, j, k;                /* index counters */
    char *seq;
    SEXP ans;                   /* return matrix */
    SEXP dims;
    int *mat;
    const int NMER = 25;        /* number of bases in each sequence */
    const int NCOL_MER = NMER * 2; /* number of columns needed for a seq. */

    if (!isString(mers)) {
        error("mers argument must be a character vector");
    }
    n = length(mers);
    PROTECT(ans = allocVector(INTSXP, n * NCOL_MER));
    mat = INTEGER(ans);
    k = 0;
    for (i = 0; i < n; i++) {
        seq = CHAR(STRING_ELT(mers, i));
        for (j = 0; j < NMER; j++) {
            if (seq[j] == 'A') {
                mat[k++] = 0;
                mat[k++] = 0;
            } else if (seq[j] == 'T') {
                mat[k++] = 0;
                mat[k++] = 1;
            } else if (seq[j] == 'C') {
                mat[k++] = 1;
                mat[k++] = 0;
            } else if (seq[j] == 'G') {
                mat[k++] = 1;
                mat[k++] = 1;
            } else {
                error("unrecognized base: %c", seq[j]);
            }
        }
    }
    PROTECT(dims = allocVector(INTSXP, 2));
    INTEGER(dims)[0] = n;
    INTEGER(dims)[1] = NCOL_MER;
    setAttrib(ans, R_DimSymbol, dims);
    UNPROTECT(2);
    return(ans);
}
