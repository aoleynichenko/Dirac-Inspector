/*
 * Inspector of DIRAC files containing transformed molecular integrals.
 *
 * 2024 Alexander Oleynichenko
 */

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "mdprop.h"
#include "mrconee.h"
#include "mdcint.h"

int main()
{
    mrconee_data_t *mrconee_data = read_mrconee("MRCONEE");
    if (mrconee_data == NULL) {
        printf(" MRCONEE file not found\n");
    }
    else {
        print_mrconee_data(stdout, mrconee_data);
    }

    read_mdprop("MDPROP", mrconee_data);
    read_mdcint("MDCINT", mrconee_data);

    return 0;
}



