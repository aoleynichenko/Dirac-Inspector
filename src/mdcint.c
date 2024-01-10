/*
 * Inspector of DIRAC files containing transformed molecular integrals.
 *
 * 2024 Alexander Oleynichenko
 */

#include "mdcint.h"

#include <stdbool.h>
#include <stdlib.h>
#include <sys/time.h>

#include "libunf.h"
#include "mrconee.h"

double abs_time();


void read_mdcint(char *path, mrconee_data_t *mrconee_data)
{
    if (mrconee_data == NULL) {
        printf(" MDCINT file cannot be parsed without auxuliary data from the MRCONEE file\n");
        return;
    }

    int use_int4 = mrconee_data->dirac_int_size == 4;
    int use_int8 = mrconee_data->dirac_int_size == 8;

    unf_file_t *mdcint = unf_open(path, "r", UNF_ACCESS_SEQUENTIAL);
    if (mdcint == NULL) {
        printf(" MDCINT file not found\n");
        return;
    }

    printf(" two-electron integrals file\n");

    /*
     * read date and time, total number of Kramers pairs
     */
    char date_time[100];
    int32_t nkr = 0;
    int64_t nkr_8 = 0;

    int nread;
    if (use_int4) {
        nread = unf_read(mdcint, "c18,i4", date_time, &nkr);
    }
    else {
        nread = unf_read(mdcint, "c18,i8", date_time, &nkr_8);
        nkr = (int32_t) nkr_8;
    }
    if (nread != 2 || unf_error(mdcint)) {
        perror(" error while reading MDCINT file");
        return;
    }

    /*
     * also read indices of Kramers pairs
     */
    unf_backspace(mdcint);

    int32_t num_spinors = 2 * nkr;
    int32_t *kr = (int32_t *) calloc(num_spinors, sizeof(int32_t *));
    int64_t *kr_8 = (int64_t *) calloc(num_spinors, sizeof(int64_t *));

    if (use_int4) {
        nread = unf_read(mdcint, "c18,i4,i4[i4]", date_time, &nkr, kr, &num_spinors);
    }
    else {
        nread = unf_read(mdcint, "c18,i8,i8[i4]", date_time, &nkr_8, kr_8, &num_spinors);
    }
    if (nread != 3 || unf_error(mdcint)) {
        perror(" error while reading MDCINT file");
        return;
    }

    if (use_int8) {
        for (int i = 0; i < num_spinors; i++) {
            kr[i] = (int32_t) kr_8[i];
        }
    }
    free(kr_8);

    date_time[18] = '\0';
    printf(" date and time              %s\n", date_time);
    printf(" number of Kramers pairs    %d\n", nkr);
    for (int i = 0; i < nkr; i++) {
        printf("%4d%4d\n", kr[2*i], kr[2*i + 1]);
    }

    /*
     * read chunks of non-zero two-electron integrals
     *
     * read (luint, end = 1301, err = 1302) ikr, jkr, nonzr, &
     * (indk(inz), indl(inz), inz = 1, nonzr), &
     * (cbuf(1, inz), inz = 1, nonzr)
     */
    int64_t count_non_zero = 0;

    char *ind_buf = (char *) calloc(num_spinors * num_spinors, 2 * sizeof(int32_t));
    char *ind_buf_8 = (char *) calloc(num_spinors * num_spinors, 2 * sizeof(int64_t));
    double *val_buf_real = (double *) calloc(num_spinors * num_spinors, sizeof(double));
    double *val_buf_complex = (double *) calloc(num_spinors * num_spinors, sizeof(double));

    int32_t ikr = 0;
    int32_t jkr = 0;
    int32_t nonzr = 0;
    int64_t ikr8 = 0;
    int64_t jkr8 = 0;
    int64_t nonzr8 = 0;
    double time_start = abs_time();

    while (true) {
        if (mrconee_data->group_arith == 1 || mrconee_data->is_spinfree == 1) {
            if (use_int4) {
                nread = unf_read(mdcint, "3i4,c8[i4],r8[i4]", &ikr, &jkr, &nonzr, ind_buf, &nonzr, val_buf_real,
                                 &nonzr);
            }
            else {
                nread = unf_read(mdcint, "3i8,c16[i8],r8[i8]", &ikr8, &jkr8, &nonzr8, ind_buf_8, &nonzr8, val_buf_real,
                                 &nonzr8);
                ikr = (int32_t) ikr8;
                jkr = (int32_t) jkr8;
                nonzr = (int32_t) nonzr8;
            }
        }
        else {
            if (use_int4) {
                nread = unf_read(mdcint, "3i4,c8[i4],z8[i4]", &ikr, &jkr, &nonzr, ind_buf, &nonzr, val_buf_complex,
                                 &nonzr);
            }
            else {
                nread = unf_read(mdcint, "3i8,c16[i8],z8[i8]", &ikr8, &jkr8, &nonzr8, ind_buf_8, &nonzr8,
                                 val_buf_complex, &nonzr8);
                ikr = (int32_t) ikr8;
                jkr = (int32_t) jkr8;
                nonzr = (int32_t) nonzr8;
            }
        }

        if (nread != 5 || unf_error(mdcint)) {

            perror(" error while reading MDCINT file");
            break;
        }

        if (ikr == 0 && jkr == 0) {
            break;
        }

        count_non_zero += nonzr;
    }

    double time_finish = abs_time();
    printf(" number of non-zero ints    %lld\n", count_non_zero);
    printf(" time for reading 2e ints   %.2f sec\n\n", time_finish - time_start);

    /*
     * cleanup
     */

    unf_close(mdcint);

    free(ind_buf);
    free(ind_buf_8);
    free(val_buf_real);
    free(val_buf_complex);
    free(kr);
}


/**
 * Interface to the system-dependent functions for time measurements.
 */
double abs_time()
{
    struct timeval cur_time;
    gettimeofday(&cur_time, NULL);
    return (cur_time.tv_sec * 1000000u + cur_time.tv_usec) / 1.e6;
}
