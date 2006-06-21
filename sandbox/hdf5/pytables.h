#ifndef __PYTABLES_COMPAT_H__
#define __PYTABLES_COMPAT_H__

#include <hdf5.h>

/* hardcoded pickle for Filter class */
#define FILTERS_P "ccopy_reg\n_reconstructor\np1\n(ctables.Leaf\nFilters\np2\nc__builtin__\nobject\np3\nNtRp4\n(dp5\nS'shuffle'\np6\nI0\nsS'complevel'\np7\nI0\nsS'fletcher32'\np8\nI0\nsS'complib'\np9\nS'zlib'\np10\nsb."


/* pytables.c */
herr_t find_attribute(hid_t loc_id, const char *attr_name);
herr_t set_attribute_string(hid_t obj_id, const char *attr_name, const char *attr_data);
herr_t make_array(hid_t loc_id, const char *dset_name, const int rank, const hsize_t *dims, hid_t type_id, const void *data);
/* cproto -D__CPROTO__ -q -p -f 3 pytables.c */

#endif
