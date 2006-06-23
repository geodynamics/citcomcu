/* hyperslab.c - Minimal example program to demonstrate the following
 *               capabilities:
 *
 *  - Creating an HDF5 file using parallel I/O (called hyperslab.h5)
 *  - Assigning each MPI process a section of 3D array (via hyperslab selection)
 *  - Adding PyTables compatibility so that data can be accessed from Python
 *
 * For simplicity, all dimensions for this program are fixed. The 3D array
 * consists of 8 x 8 x 8 values, distributed over 2 x 2 x 2 procs. Each
 * proc, then, is responsible for a local mesh of dimensions 4 x 4 x 4.
 * 
 *
 * To compile, use:
 *
 *      mpicc hyperslab.c pytables.o -o hyperslab -lhdf5
 *
 * To run, use:
 *
 *      mpirun -np 8 hyperslab
 *
 * To view contents of output file, use:
 *      
 *      h5dump hyperslab.h5
 *
 * To use from within python:
 *
 *      >> import tables
 *      >> file = tables.openFile('hyperslab.h5', 'r')
 *      >> mesh = file.root.proc_mesh
 *      >> print mesh[:,:,0]  # select xy-plane slice with z=0
 *      [[ 1.  1.  1.  1.  5.  5.  5.  5.]
 *       [ 1.  1.  1.  1.  5.  5.  5.  5.]
 *       [ 1.  1.  1.  1.  5.  5.  5.  5.]
 *       [ 1.  1.  1.  1.  5.  5.  5.  5.]
 *       [ 3.  3.  3.  3.  7.  7.  7.  7.]
 *       [ 3.  3.  3.  3.  7.  7.  7.  7.]
 *       [ 3.  3.  3.  3.  7.  7.  7.  7.]
 *       [ 3.  3.  3.  3.  7.  7.  7.  7.]]
 *      >>
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <hdf5.h>
#include "pytables.h"

struct Mesh3D
{
    int rank;
    hsize_t dims[3];
    int nno;
    float *data;
};

int main(int argc, char *argv[])
{
    /* HDF5 API definitions */
    hid_t file;
    hid_t root;
    hid_t dset;
    hid_t filespace;
    hid_t memspace;
    hid_t plist;

    hsize_t count[3];   /* hyperslab selection parameters */
    hsize_t stride[3];
    hsize_t block[3];
    hsize_t offset[3];

    herr_t status;

    /* MPI variables */
    int mpi_size, mpi_rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;

    int nprocs[3];
    int loc[3];

    /* data variables */
    int ijk;
    int i, j, k;
    struct Mesh3D mesh;
    struct Mesh3D lmesh;



    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    if(mpi_size != 8)
    {
        if(mpi_rank == 0)
            printf("This example is set up to use only 8 processes!\n");
        return EXIT_FAILURE;
    }



    /* set up file access property list with parallel I/O access */
    plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, comm, info);

    /* create a new file collectively and release property list identifier */
    file = H5Fcreate("hyperslab.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    H5Pclose(plist);

    /* create the dataspace for the dataset */
    mesh.rank = 3;
    mesh.dims[0] = 8;
    mesh.dims[1] = 8;
    mesh.dims[2] = 8;
    mesh.nno = 8*8*8;
    filespace = H5Screate_simple(mesh.rank, mesh.dims, NULL);

    lmesh.rank = 3;
    lmesh.dims[0] = 4;
    lmesh.dims[1] = 4;
    lmesh.dims[2] = 4;
    lmesh.nno = 4*4*4;
    memspace = H5Screate_simple(lmesh.rank, lmesh.dims, NULL);

    nprocs[0] = mesh.dims[0] / lmesh.dims[0];
    nprocs[1] = mesh.dims[1] / lmesh.dims[1];
    nprocs[2] = mesh.dims[2] / lmesh.dims[2];

    loc[2] =    mpi_rank % nprocs[2];
    loc[0] =  ((mpi_rank - loc[2]) / nprocs[2]) % nprocs[0];
    loc[1] = (((mpi_rank - loc[2]) / nprocs[2]) - loc[0]) / nprocs[0];

    /* create chunked dataset */
    plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist, lmesh.rank, lmesh.dims);
    dset = H5Dcreate(file, "/proc_mesh", H5T_NATIVE_FLOAT, filespace,
                     H5P_DEFAULT);
    H5Pclose(plist);
    H5Sclose(filespace);

    /* each process defines dataset in memory and writes it
     * to the hyperslab in the file.
     */
    count[0] = 1;
    count[1] = 1;
    count[2] = 1;
    stride[0] = 1;
    stride[1] = 1;
    stride[2] = 1;
    block[0] = lmesh.dims[0];
    block[1] = lmesh.dims[1];
    block[2] = lmesh.dims[2];
    offset[0] = loc[0] * block[0];
    offset[1] = loc[1] * block[1];
    offset[2] = loc[2] * block[2];

    /* select hyperslab in the file */
    filespace = H5Dget_space(dset);
    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, 
                                 offset, stride, count, block);

    /* initialize local data buffer */
    lmesh.data = (float *)malloc(lmesh.nno * sizeof(float));
    for(j = 0; j < lmesh.dims[1]; j++)
    {
        for(i = 0; i < lmesh.dims[0]; i++)
        {
            for(k = 0; k < lmesh.dims[2]; k++)
            {
                ijk = k + j*lmesh.dims[2] + i*lmesh.dims[2]*lmesh.dims[1];
                lmesh.data[ijk] = (float)mpi_rank + 1;
            }
        }
    }

    /* create property list for collective dataset write */
    plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

    /* to write dataset independently use
     * H5Pset_dxpl_mpio(plist, H5FD_MPIO_INDEPENDENT);
     */

    status = H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace,
                      plist, lmesh.data);
    free(lmesh.data);

    /* pytables compatibility */
    root = H5Gopen(file, "/");

    set_attribute_string(root, "TITLE", "PyTables Numeric array");
    set_attribute_string(root, "CLASS", "GROUP");
    set_attribute_string(root, "VERSION", "1.0");
    set_attribute_string(root, "FILTERS", FILTERS_P);
    set_attribute_string(root, "PYTABLES_FORMAT_VERSION", "1.5");

    set_attribute_string(dset, "TITLE", "Processor number in mesh");
    set_attribute_string(dset, "CLASS", "ARRAY");
    set_attribute_string(dset, "FLAVOR", "Numeric");
    set_attribute_string(dset, "VERSION", "2.3");

    /* close/release resources */
    H5Gclose(root);
    H5Dclose(dset);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist);
    H5Fclose(file);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
