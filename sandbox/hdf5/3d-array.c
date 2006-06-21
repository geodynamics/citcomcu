#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include "pytables.h"

int main(int argc, char *argv[])
{
    hid_t file;
    hid_t root;
    hid_t box_dset;

    int rank;
    hsize_t dims[3];

    int i,j,k;
    int n, nno;
    float *data;

    herr_t status;

    /* open file */
    file = H5Fcreate("3d-array1.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* open root group */
    root = H5Gopen(file, "/");

    /* root attributes */
    set_attribute_string(root, "TITLE", "pytables Numeric array!");
    set_attribute_string(root, "CLASS", "GROUP");
    set_attribute_string(root, "VERSION", "1.0");
    set_attribute_string(root, "FILTERS", FILTERS_P);
    set_attribute_string(root, "PYTABLES_FORMAT_VERSION", "1.5");

    /* array of shape (5,6,2) */
    rank = 3;
    dims[0] = 5;
    dims[1] = 6;
    dims[2] = 7;
    nno = dims[0]*dims[1]*dims[2];

    /* data */
    n = 0;
    data = (float *)malloc(nno*sizeof(float));
    for(j = 0; j < dims[1]; j++)
    {
        for(i = 0; i < dims[0]; i++)
        {
            for(k = 0; k < dims[2]; k++)
            {
                data[i*dims[2]*dims[1] + j*dims[2] + k] = (float)n;
                n++;
            }
        }
    }

    /* write data */
    box_dset = make_array(root, "box", rank, dims, H5T_NATIVE_FLOAT, data);

    /* free memory */
    free(data);

    /* close root group */
    status = H5Gclose(root);

    /* close file */
    status = H5Fclose(file);

    return EXIT_SUCCESS;
}
