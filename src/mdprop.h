/*
 * Inspector of DIRAC files containing transformed molecular integrals.
 *
 * 2024 Alexander Oleynichenko
 */

#ifndef DIRAC_INSPECTOR_MDPROP_H
#define DIRAC_INSPECTOR_MDPROP_H

#include "mrconee.h"

void read_mdprop(char *path, mrconee_data_t *mrconee_data);

#endif // DIRAC_INSPECTOR_MDPROP_H
