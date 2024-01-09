//
// Created by Alexander Oleynichenko on 09.01.2024.
//

#include "mdcint.h"

#include <stdlib.h>

#include "libunf.h"
#include "mrconee.h"

void read_mdcint(char *path, mrconee_data_t *mrconee_data)
{
    unf_file_t *mdcint = unf_open(path, "r", UNF_ACCESS_SEQUENTIAL);
    if (mdcint == NULL) {
        printf(" MDCINT file not found\n");
        return;
    }

    printf(" two-electron integrals file\n");

    // read date and time, total number of Kramers pairs
    char date_time[100];
    int32_t nkr = 0;
    int64_t nkr8 = 0;
    int use_int4 = mrconee_data->dirac_int_size == 4;

    int nread = unf_read(mdcint, use_int4 ? "c18,i4" : "c18,i8", date_time, use_int4 ? &nkr : &nkr8);
    if (nread != 2 || unf_error(mdcint)) {
        perror(" error while reading mdcint");
        return;
    }

    if (!use_int4) {
        nkr = (int32_t) nkr8;
    }

    // read indices of Kramers pairs
    int32_t num_spinors = 2 * nkr;
    int32_t *kr = (int32_t *) calloc(num_spinors, sizeof(int32_t *));
    int64_t *kr8 = (int64_t *) calloc(num_spinors, sizeof(int64_t *));

    unf_backspace(mdcint);
    if (use_int4) {
        nread = unf_read(mdcint, "c18,i4,i4[i4]", date_time, &nkr, kr, &num_spinors);
    }
    else {
        nread = unf_read(mdcint, "c18,i8,i8[i4]", date_time, &nkr8, kr8, &num_spinors);
    }
    if (nread != 3 || unf_error(mdcint)) {
        perror(" error while reading mdcint");
        return;
    }

    if (!use_int4) {
        for (int i = 0; i < num_spinors; i++) {
            kr[i] = (int32_t) kr8[i];
        }
    }

    date_time[18] = '\0';
    printf(" date and time              %s\n", date_time);
    printf(" number of Kramers pairs    %d\n", nkr);
    for (int i = 0; i < nkr; i++) {
        printf("%4d%4d\n", kr[2*i], kr[2*i + 1]);
    }

    /*
     * read (luint, end = 1301, err = 1302) ikr, jkr, nonzr, &
     * (indk(inz), indl(inz), inz = 1, nonzr), &
     * (cbuf(1, inz), inz = 1, nonzr)
     */
    int32_t ikr = 0, jkr = 0, nonzr = 0;
    int64_t ikr8 = 0, jkr8 = 0, nonzr8 = 0;
    int64_t count_non_zero = 0;

    char *ind_buf = (char *) calloc(num_spinors * num_spinors, 2 * sizeof(int32_t));
    char *ind_buf8 = (char *) calloc(num_spinors * num_spinors, 2 * sizeof(int64_t));
    double *cbuf = (double *) calloc(num_spinors * num_spinors, sizeof(double));

    while (1) {
        if (use_int4) {
            nread = unf_read(mdcint, "3i4,c8[i4],r8[i4]", &ikr, &jkr, &nonzr, ind_buf, &nonzr, cbuf, &nonzr);
        }
        else {
            nread = unf_read(mdcint, "3i8,c16[i8],r8[i8]", &ikr8, &jkr8, &nonzr8, ind_buf8, &nonzr8, cbuf, &nonzr8);
        }
        if (nread != 5 || unf_error(mdcint)) {
            perror(" error while reading mdcint");
            return;
        }

        ikr = (int32_t) ikr8;
        jkr = (int32_t) jkr8;
        nonzr = (int32_t) nonzr8;

        if (ikr == 0 && jkr == 0) {
            break;
        }

        count_non_zero += nonzr;
    }

    printf(" number of non-zero ints    %lld\n\n", count_non_zero);

    // cleanup
    unf_close(mdcint);
    free(ind_buf);
    free(ind_buf8);
    free(cbuf);
    free(kr);
}
