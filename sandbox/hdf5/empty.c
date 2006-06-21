#include <stdlib.h>
#include <hdf5.h>

int main(void)
{
    hid_t file;
    herr_t status;

    /* Create a new file using default properties */
    file = H5Fcreate("empty1.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Terminate access to the file */
    status = H5Fclose(file);

    return EXIT_SUCCESS;
}
