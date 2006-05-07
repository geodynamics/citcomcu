/*
 * CitcomCU is a Finite Element Code that solves for thermochemical
 * convection within a three dimensional domain appropriate for convection
 * within the Earth's mantle. Cartesian and regional-spherical geometries
 * are implemented. See the file README contained with this distribution
 * for further details.
 * 
 * Copyright (C) 1994-2005 California Institute of Technology
 * Copyright (C) 2000-2005 The University of Colorado
 *
 * Authors: Louis Moresi, Shijie Zhong, and Michael Gurnis
 *
 * For questions or comments regarding this software, you may contact
 *
 *     Luis Armendariz <luis@geodynamics.org>
 *     http://geodynamics.org
 *     Computational Infrastructure for Geodynamics (CIG)
 *     California Institute of Technology
 *     2750 East Washington Blvd, Suite 210
 *     Pasadena, CA 91007
 *
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 2 of the License, or any
 * later version.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program; if not, write to the Free Software 
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

/* Most routines on parallel computing are here. Shijie Zhong
 * wrote these in 1996 initially for intel paragon and then gradually
 * ported to other platforms, T3E, and Linux clusters.
 */

#include <mpi.h>

#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"


void parallel_process_initialization(struct All_variables *E,
                                     int argc, char **argv)
{
    E->parallel.me = 0;
    E->parallel.nproc = 1;
    E->parallel.me_loc[1] = 0;
    E->parallel.me_loc[2] = 0;
    E->parallel.me_loc[3] = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &(E->parallel.me));
    MPI_Comm_size(MPI_COMM_WORLD, &(E->parallel.nproc));

    return;
}


/* ============================================ */
/* ============================================ */

void parallel_process_termination(int code)
{
    MPI_Finalize();
    exit(code);
    return;
}


/* ============================================ */
/* ============================================ */


void parallel_domain_decomp1(struct All_variables *E)
{
    int i, j, k;
    int nox, noz, noy;
    
    int me = E->parallel.me;
    
    int levmax = E->mesh.levmax;
    int levmin = E->mesh.levmin;
    
    int nprocx = E->parallel.nprocx;
    int nprocz = E->parallel.nprocz;
    int nprocy = E->parallel.nprocy;

    int nprocxz = E->parallel.nprocx * E->parallel.nprocz;
    int nprocxy = E->parallel.nprocx * E->parallel.nprocy;
    int nproczy = E->parallel.nprocz * E->parallel.nprocy;
    

    if(E->parallel.automa)
    {
    }

    if(E->control.NMULTIGRID || E->control.EMULTIGRID)
    {
        E->lmesh.elx = E->lmesh.mgunitx * (int)pow(2.0, (double)levmax);
        E->lmesh.elz = E->lmesh.mgunitz * (int)pow(2.0, (double)levmax);
        E->lmesh.ely = E->lmesh.mgunity * (int)pow(2.0, (double)levmax);
    }

    E->lmesh.elx = E->mesh.elx / E->parallel.nprocx;
    E->lmesh.elz = E->mesh.elz / E->parallel.nprocz;
    E->lmesh.ely = E->mesh.ely / E->parallel.nprocy;

    E->parallel.nprocxz = nprocxz;
    E->parallel.nprocxy = nprocxy;
    E->parallel.nproczy = nproczy;

    k = 0;
    for(j = 0; j < E->parallel.nproc; j++)
    {
        for(i = 0; i <= j; i++)
        {
            E->parallel.mst[j][i][1] = k++;
            E->parallel.mst[j][i][2] = k++;
        }
    }
    for(j = 0; j < E->parallel.nproc; j++)
    {
        for(i = 0; i <= E->parallel.nproc; i++)
            if(i > j)
            {
                E->parallel.mst[j][i][1] = E->parallel.mst[i][j][2];
                E->parallel.mst[j][i][2] = E->parallel.mst[i][j][1];
            }
    }

    /* for overlapping domain, good for e by e assemble */

    /* z direction first */
    j = me % E->parallel.nprocz;
    E->parallel.me_loc[3] = j;
    E->lmesh.nzs = j * E->lmesh.elz + 1;

    /* y direction then */
    k = (me + 1)/nprocxz - (((me + 1) % nprocxz == 0) ? 1 : 0);
    E->parallel.me_loc[2] = k;
    E->lmesh.nys = k * E->lmesh.ely + 1;

    /* x direction then */
    i = (me + 1 - k * nprocxz)/nprocz 
            - (((me + 1 - k * nprocxz) % nprocz == 0) ? 1 : 0);
    E->parallel.me_loc[1] = i;
    E->lmesh.nxs = i * E->lmesh.elx + 1;

/*fprintf(stderr,"b %d %d %d %d %d %d %d\n",E->parallel.me,E->parallel.me_loc[1],E->parallel.me_loc[2],E->parallel.me_loc[3],E->lmesh.nxs,E->lmesh.nzs,E->lmesh.nys); */

    E->lmesh.nox = E->lmesh.nnx[1] = E->lmesh.elx + 1;
    E->lmesh.noy = E->lmesh.nnx[2] = E->lmesh.ely + 1;
    E->lmesh.noz = E->lmesh.nnx[3] = E->lmesh.elz + 1;

    E->lmesh.nno = E->lmesh.noz * E->lmesh.nox * E->lmesh.noy;
    E->lmesh.nsf = E->lmesh.nox * E->lmesh.noy;
    E->lmesh.npno = E->lmesh.elz * E->lmesh.elx * E->lmesh.ely;
    E->lmesh.nel = E->lmesh.elz * E->lmesh.elx * E->lmesh.ely;
    E->lmesh.nex[1] = E->lmesh.elx;
    E->lmesh.nex[2] = E->lmesh.ely;
    E->lmesh.nex[3] = E->lmesh.elz;
    E->lmesh.nnov = E->lmesh.nno;
    E->lmesh.neq = E->lmesh.nnov * E->mesh.nsd;

    for(i = levmax; i >= levmin; i--)
    {
        if(E->control.NMULTIGRID || E->control.EMULTIGRID)
        {
            nox = E->lmesh.elx / ((int)pow(2.0, (double)(levmax - i))) + 1;
            noz = E->lmesh.elz / ((int)pow(2.0, (double)(levmax - i))) + 1;
            noy = E->lmesh.ely / ((int)pow(2.0, (double)(levmax - i))) + 1;
        }
        else
        {
            nox = E->lmesh.nox;
            noz = E->lmesh.noz;
            noy = E->lmesh.noy;
        }

        E->lmesh.ELX[i] = nox - 1;
        E->lmesh.ELZ[i] = noz - 1;
        E->lmesh.ELY[i] = max(noy - 1, 1);
        E->lmesh.NNO[i] = nox * noz * noy;
        E->lmesh.NEL[i] = (nox - 1) * (noz - 1) * max((noy - 1), 1);
        E->lmesh.NPNO[i] = E->lmesh.NEL[i];
        E->lmesh.NOX[i] = nox;
        E->lmesh.NOZ[i] = noz;
        E->lmesh.NOY[i] = noy;
        E->lmesh.NNOV[i] = E->lmesh.NNO[i];
        E->lmesh.NEQ[i] = E->mesh.nsd * E->lmesh.NNOV[i];
        E->lmesh.NXS[i] = E->parallel.me_loc[1] * E->lmesh.ELX[i] + 1;
        E->lmesh.NYS[i] = E->parallel.me_loc[2] * E->lmesh.ELY[i] + 1;
        E->lmesh.NZS[i] = E->parallel.me_loc[3] * E->lmesh.ELZ[i] + 1;
/*fprintf(E->fp,"b %d %d %d %d %d %d %d\n",E->parallel.me,E->lmesh.ELX[i],E->lmesh.ELZ[i],E->lmesh.ELY[i],E->lmesh.NNO[i],E->lmesh.NEL[i],E->lmesh.NEQ[i]); */

    }

    return;
}

/* ============================================ 
 shuffle element order and determine boundary nodes for 
 exchange info across the boundaries
 ============================================ */

void parallel_shuffle_ele_and_id(struct All_variables *E)
{
    if(E->mesh.periodic_x || E->mesh.periodic_y)
        parallel_shuffle_ele_and_id_bc2(E);
    else
        parallel_shuffle_ele_and_id_bc1(E);
    return;
}

void parallel_shuffle_ele_and_id_bc1(struct All_variables *E)
{
    //int i, ii, j, k, l, node, node1, el, elt, lnode, llnode, jj, k1, k2;
    int i, ii, j, k, node, node1, llnode;
    //int lev, elx, elz, ely, nel, nno, nox, noz, noy;
    int lev, elx, elz, ely, nel, nox, noz, noy;

    for(lev = E->mesh.levmax; lev >= E->mesh.levmin; lev--)
    {
        nel = E->lmesh.NEL[lev];
        elx = E->lmesh.ELX[lev];
        elz = E->lmesh.ELZ[lev];
        ely = E->lmesh.ELY[lev];
        nox = E->lmesh.NOX[lev];
        noy = E->lmesh.NOY[lev];
        noz = E->lmesh.NOZ[lev];

        ii = 0;

        for(i = 1; i <= 2; i++)
        {                       /* do the ZOY boundary elements first */
            ii++;
            if((i == 1 && E->parallel.me_loc[1] != 0) ||
               (i == 2 && E->parallel.me_loc[1] != E->parallel.nprocx - 1))
            {
                for(k = 1; k <= noy; k++)
                    for(j = 1; j <= noz; j++)
                    {
                        node = j + (((i==1) ? 1 : nox) - 1)*noz + (k-1)*noz*nox;
                        llnode = j + (k - 1) * noz;
                        E->parallel.NODE[lev][llnode].bound[1][i] = node;
                        E->NODE[lev][node] = E->NODE[lev][node] | OFFSIDE;
                        E->NODE[lev][node] = E->NODE[lev][node] | LIDN;
                        node1 = node + (((i == 1) ? 1 : -1)) * noz;
                        E->NODE[lev][node1] = E->NODE[lev][node1] | LIDN;
                    }
                E->parallel.NUM_NNO[lev].bound[1][i] = noy * noz;
            }
        }                       /* end for i   */

        for(j = 1; j <= 2; j++)
        {                       /* do XOY boundary elements */
            ii++;
            if((j == 1 && E->parallel.me_loc[3] != 0) ||
               (j == 2 && E->parallel.me_loc[3] != E->parallel.nprocz - 1))
            {
                for(k = 1; k <= noy; k++)
                    for(i = 1; i <= nox; i++)
                    {
                        node = ((j==1) ? 1 : noz) + (i-1)*noz + (k-1)*noz*nox;
                        llnode = i + (k - 1) * nox;
                        E->parallel.NODE[lev][llnode].bound[3][j] = node;
                        E->NODE[lev][node] = E->NODE[lev][node] | OFFSIDE;
                        E->NODE[lev][node] = E->NODE[lev][node] | LIDN;
                        node1 = node + (((j == 1) ? 1 : -1));
                        E->NODE[lev][node1] = E->NODE[lev][node1] | LIDN;
                    }
                E->parallel.NUM_NNO[lev].bound[3][j] = nox * noy;
            }
        }                       /* end for j   */

        for(k = 1; k <= 2; k++)
        {                       /* do XOZ boundary elements for 3D */
            ii++;
            if((k == 1 && E->parallel.me_loc[2] != 0) ||
               (k == 2 && E->parallel.me_loc[2] != E->parallel.nprocy - 1))
            {
                for(j = 1; j <= noz; j++)
                    for(i = 1; i <= nox; i++)
                    {
                        node = j + (i-1)*noz + (((k==1) ? 1 : noy) - 1)*noz*nox;
                        llnode = j + (i - 1) * noz;
                        E->parallel.NODE[lev][llnode].bound[2][k] = node;
                        E->NODE[lev][node] = E->NODE[lev][node] | OFFSIDE;
                        E->NODE[lev][node] = E->NODE[lev][node] | LIDN;
                        node1 = node + (((k == 1) ? 1 : -1)) * noz * nox;
                        E->NODE[lev][node1] = E->NODE[lev][node1] | LIDN;
                    }
                E->parallel.NUM_NNO[lev].bound[2][k] = nox * noz;
            }
        }                       /* end for k */


        E->parallel.num_b = ii;

    }                           /* end for level */

    return;
}

// for periodic BC 
//
void parallel_shuffle_ele_and_id_bc2(struct All_variables *E)
{
    //int i, ii, j, k, l, node, node1, el, elt, lnode, llnode, jj, k1, k2;
    int i, ii, j, k, node, node1, llnode;
    //int lev, elx, elz, ely, nel, nno, nox, noz, noy;
    int lev, elx, elz, ely, nel, nox, noz, noy;

    for(lev = E->mesh.levmax; lev >= E->mesh.levmin; lev--)
    {
        nel = E->lmesh.NEL[lev];
        elx = E->lmesh.ELX[lev];
        elz = E->lmesh.ELZ[lev];
        ely = E->lmesh.ELY[lev];
        nox = E->lmesh.NOX[lev];
        noy = E->lmesh.NOY[lev];
        noz = E->lmesh.NOZ[lev];

        ii = 0;

        for(i = 1; i <= 2; i++)
        {                       /* do the ZOY boundary elements first */
            ii++;
            for(k = 1; k <= noy; k++)
                for(j = 1; j <= noz; j++)
                {
                    node = j + (((i==1) ? 1 : nox) - 1)*noz + (k - 1)*noz*nox;
                    llnode = j + (k - 1) * noz;
                    E->parallel.NODE[lev][llnode].bound[1][i] = node;
                    E->NODE[lev][node] = E->NODE[lev][node] | OFFSIDE;
                    E->NODE[lev][node] = E->NODE[lev][node] | LIDN;
                    node1 = node + (((i == 1) ? 1 : -1)) * noz;
                    E->NODE[lev][node1] = E->NODE[lev][node1] | LIDN;
                }
            E->parallel.NUM_NNO[lev].bound[1][i] = noy * noz;
        }                       /* end for i   */

        for(j = 1; j <= 2; j++)
        {                       /* do XOY boundary elements */
            ii++;
            if((j == 1 && E->parallel.me_loc[3] != 0) ||
               (j == 2 && E->parallel.me_loc[3] != E->parallel.nprocz - 1))
            {
                for(k = 1; k <= noy; k++)
                    for(i = 1; i <= nox; i++)
                    {
                        node = ((j==1) ? 1 : noz) + (i-1)*noz + (k-1)*noz*nox;
                        llnode = i + (k - 1) * nox;
                        E->parallel.NODE[lev][llnode].bound[3][j] = node;
                        E->NODE[lev][node] = E->NODE[lev][node] | OFFSIDE;
                        E->NODE[lev][node] = E->NODE[lev][node] | LIDN;
                        node1 = node + (((j == 1) ? 1 : -1));
                        E->NODE[lev][node1] = E->NODE[lev][node1] | LIDN;
                    }
            }
            E->parallel.NUM_NNO[lev].bound[3][j] = nox * noy;
        }                       /* end for j   */

        for(k = 1; k <= 2; k++)
        {                       /* do XOZ boundary elements for 3D */
            ii++;
            for(j = 1; j <= noz; j++)
                for(i = 1; i <= nox; i++)
                {
                    node = j + (i-1)*noz + (((k==1) ? 1 : noy) - 1)*noz*nox;
                    llnode = j + (i - 1) * noz;
                    E->parallel.NODE[lev][llnode].bound[2][k] = node;
                    E->NODE[lev][node] = E->NODE[lev][node] | OFFSIDE;
                    E->NODE[lev][node] = E->NODE[lev][node] | LIDN;
                    node1 = node + (((k == 1) ? 1 : -1)) * noz * nox;
                    E->NODE[lev][node1] = E->NODE[lev][node1] | LIDN;
                }
            E->parallel.NUM_NNO[lev].bound[2][k] = nox * noz;
        }                       /* end for k */


        E->parallel.num_b = ii;

    }                           /* end for level */

    return;
}



/* ============================================ 
 determine communication routs  and boundary ID for 
 exchange info across the boundaries
 ============================================ */

void parallel_communication_routs(struct All_variables *E)
{
    if(E->mesh.periodic_x || E->mesh.periodic_y)
    {
        parallel_communication_routs2(E);
        if(E->control.composition)
            parallel_communication_routs4(E);
    }
    else
    {
        parallel_communication_routs1(E);
        if(E->control.composition)
            parallel_communication_routs3(E);
    }

    return;
}



void parallel_communication_routs1(struct All_variables *E)
{
    //int i, ii, j, k, l, node, el, elt, lnode, jj, doff;
    int i, ii, j, k, node, lnode, jj, doff;
    int lev, elx, elz, ely, nno, nox, noz, noy, p, kkk, kk;
    int me, nprocz, nprocx, nprocy;

    me = E->parallel.me;
    nprocx = E->parallel.nprocx;
    nprocz = E->parallel.nprocz;
    nprocy = E->parallel.nprocy;

    for(lev = E->mesh.levmax; lev >= E->mesh.levmin; lev--)
    {
        elx = E->lmesh.ELX[lev];
        elz = E->lmesh.ELZ[lev];
        ely = E->lmesh.ELY[lev];
        nox = E->lmesh.NOX[lev];
        noy = E->lmesh.NOY[lev];
        noz = E->lmesh.NOZ[lev];
        nno = E->lmesh.NNO[lev];

        ii = 0;
        kkk = 0;
        if(E->mesh.nsd == 2)
        {

            for(i = 1; i <= 2; i++)
            {                   /* do the OZ boundaries */

                ii++;

                E->parallel.NUM_PASS[lev].bound[1][i] = 1;
                if(E->parallel.me_loc[1] == 0 && i == 1)
                    E->parallel.NUM_PASS[lev].bound[1][i] = 0;
                else if(E->parallel.me_loc[1] == nprocx - 1 && i == 2)
                    E->parallel.NUM_PASS[lev].bound[1][i] = 0;

                for(p = 1; p <= E->parallel.NUM_PASS[lev].bound[1][i]; p++)
                {

                    kkk++;
                    /* determine the pass ID for ii-th boundary and p-th pass */

                    /* for kkk th pass, get the # of data to be passed */
                    E->parallel.NUM_NEQ[lev].pass[1][i] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[1][i] : noy) * E->mesh.nsd;
                    E->parallel.NUM_NODE[lev].pass[1][i] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[1][i] : noy);

                    /* for kkk th pass, determine the target processor */

                    E->parallel.PROCESSOR[lev].pass[1][i] = me - ((i == 1) ? 1 : -1) * nprocz;

                    /* for kkk th pass, determine the ID of equations to be passed */

                    for(node = 1; node <= E->parallel.NUM_NODE[lev].pass[1][i]; node++)
                    {
                        lnode = E->parallel.NODE[lev][node].bound[1][i];
                        E->parallel.EXCHANGE_NODE[lev][node].pass[1][i] = lnode;
                        for(doff = 1; doff <= E->mesh.nsd; doff++)
                        {
                            jj = doff + (node - 1) * E->mesh.nsd;
                            E->parallel.EXCHANGE_ID[lev][jj].pass[1][i] = E->ID[lev][lnode].doff[doff];
                        }
                    }

                }               /* end for loop p */
            }                   /* end for loop i */

            for(j = 1; j <= 2; j++)
            {                   /* do OX boundary */

                ii++;

                E->parallel.NUM_PASS[lev].bound[2][j] = 1;
                if(E->parallel.me_loc[2] == 0 && j == 1)
                    E->parallel.NUM_PASS[lev].bound[2][j] = 0;
                else if(E->parallel.me_loc[2] == nprocz - 1 && j == 2)
                    E->parallel.NUM_PASS[lev].bound[2][j] = 0;

                for(p = 1; p <= E->parallel.NUM_PASS[lev].bound[2][j]; p++)
                {
                    kkk++;

                    /* for kkk th pass, get the # of data to be passed */
                    E->parallel.NUM_NEQ[lev].pass[2][j] = E->parallel.NUM_NNO[lev].bound[2][j] * E->mesh.nsd;
                    E->parallel.NUM_NODE[lev].pass[2][j] = E->parallel.NUM_NNO[lev].bound[2][j];

                    /* for kkk th pass, determine the target processor */

                    E->parallel.PROCESSOR[lev].pass[2][j] = me - ((j == 1) ? 1 : -1);

                    /* for kkk th pass, determine the ID of equations to be passed */

                    for(node = 1; node <= E->parallel.NUM_NODE[lev].pass[2][j]; node++)
                    {
                        lnode = E->parallel.NODE[lev][node].bound[2][j];
                        E->parallel.EXCHANGE_NODE[lev][node].pass[2][j] = lnode;
                        for(doff = 1; doff <= E->mesh.nsd; doff++)
                        {
                            jj = doff + (node - 1) * E->mesh.nsd;
                            E->parallel.EXCHANGE_ID[lev][jj].pass[2][j] = E->ID[lev][lnode].doff[doff];
                        }
                    }           /* end loop for node */

                }               /* end for loop p */
            }                   /* end for loop j */

        }                       /* end for 2D */



        else if(E->mesh.nsd == 3)
        {

            for(i = 1; i <= 2; i++)
            {                   /* do YOZ boundaries & OY lines */

                ii++;
                E->parallel.NUM_PASS[lev].bound[1][i] = 1;
                if(E->parallel.me_loc[1] == 0 && i == 1)
                    E->parallel.NUM_PASS[lev].bound[1][i] = 0;
                else if(E->parallel.me_loc[1] == nprocx - 1 && i == 2)
                    E->parallel.NUM_PASS[lev].bound[1][i] = 0;

                for(p = 1; p <= E->parallel.NUM_PASS[lev].bound[1][i]; p++)
                {
                    kkk++;
                    /* determine the pass ID for ii-th boundary and p-th pass */
                    E->parallel.NUM_NEQ[lev].pass[1][i] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[1][i] : noy) * E->mesh.nsd;
                    E->parallel.NUM_NODE[lev].pass[1][i] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[1][i] : noy);

                    E->parallel.PROCESSOR[lev].pass[1][i] = me - ((i == 1) ? 1 : -1) * nprocz;

                    for(k = 1; k <= E->parallel.NUM_NODE[lev].pass[1][i]; k++)
                    {
                        lnode = k;
                        node = E->parallel.NODE[lev][lnode].bound[1][i];
                        E->parallel.EXCHANGE_NODE[lev][k].pass[1][i] = node;
                        for(doff = 1; doff <= E->mesh.nsd; doff++)
                        {
                            jj = doff + (k - 1) * E->mesh.nsd;
                            E->parallel.EXCHANGE_ID[lev][jj].pass[1][i] = E->ID[lev][node].doff[doff];
                        }
                    }           /* end for node k */

                }               /* end for loop p */
            }                   /* end for i */

            for(j = 1; j <= 2; j++)
            {                   /* do XOY boundaries & OX lines */

                ii++;
                E->parallel.NUM_PASS[lev].bound[3][j] = 1;
                if(E->parallel.me_loc[3] == 0 && j == 1)
                    E->parallel.NUM_PASS[lev].bound[3][j] = 0;
                else if(E->parallel.me_loc[3] == nprocz - 1 && j == 2)
                    E->parallel.NUM_PASS[lev].bound[3][j] = 0;

                for(p = 1; p <= E->parallel.NUM_PASS[lev].bound[3][j]; p++)
                {
                    kkk++;
                    /* determine the pass ID for ii-th boundary and p-th pass */
                    E->parallel.NUM_NEQ[lev].pass[3][j] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[3][j] : nox) * E->mesh.nsd;
                    E->parallel.NUM_NODE[lev].pass[3][j] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[3][j] : nox);

                    E->parallel.PROCESSOR[lev].pass[3][j] = me - ((j == 1) ? 1 : -1);

                    for(k = 1; k <= E->parallel.NUM_NODE[lev].pass[3][j]; k++)
                    {
                        lnode = k;
                        node = E->parallel.NODE[lev][lnode].bound[3][j];
                        E->parallel.EXCHANGE_NODE[lev][k].pass[3][j] = node;
                        for(doff = 1; doff <= E->mesh.nsd; doff++)
                        {
                            jj = doff + (k - 1) * E->mesh.nsd;
                            E->parallel.EXCHANGE_ID[lev][jj].pass[3][j] = E->ID[lev][node].doff[doff];
                        }
                    }           /* end for node k */

                }               /* end for loop p */

            }                   /* end for j */

            for(k = 1; k <= 2; k++)
            {                   /* do XOZ boundaries & OZ lines */

                ii++;
                E->parallel.NUM_PASS[lev].bound[2][k] = 1;
                if(E->parallel.me_loc[2] == 0 && k == 1)
                    E->parallel.NUM_PASS[lev].bound[2][k] = 0;
                else if(E->parallel.me_loc[2] == nprocy - 1 && k == 2)
                    E->parallel.NUM_PASS[lev].bound[2][k] = 0;

                for(p = 1; p <= E->parallel.NUM_PASS[lev].bound[2][k]; p++)
                {

                    kkk++;
                    /* determine the pass ID for ii-th boundary and p-th pass */
                    E->parallel.NUM_NEQ[lev].pass[2][k] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[2][k] : noz) * E->mesh.nsd;
                    E->parallel.NUM_NODE[lev].pass[2][k] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[2][k] : noz);

                    E->parallel.PROCESSOR[lev].pass[2][k] = me - ((k == 1) ? 1 : -1) * nprocx * nprocz;

                    for(kk = 1; kk <= E->parallel.NUM_NODE[lev].pass[2][k]; kk++)
                    {
                        lnode = kk;
                        node = E->parallel.NODE[lev][lnode].bound[2][k];
                        E->parallel.EXCHANGE_NODE[lev][kk].pass[2][k] = node;
                        for(doff = 1; doff <= E->mesh.nsd; doff++)
                        {
                            jj = doff + (kk - 1) * E->mesh.nsd;
                            E->parallel.EXCHANGE_ID[lev][jj].pass[2][k] = E->ID[lev][node].doff[doff];
                        }
                    }           /* end for node kk */

                }               /* end for loop p */

            }                   /* end for k */

        }                       /* end of dims==3 */

        E->parallel.TNUM_PASS[lev] = kkk;


/*    fprintf(E->fp," me= %d  pass  %d \n",E->parallel.me,E->parallel.TNUM_PASS[lev]);
    for (k=1;k<=E->parallel.TNUM_PASS[lev];k++)  
fprintf(E->fp,"proc %d and pass  %d to proc %d with %d eqn\n",E->parallel.me,k,E->parallel.PROCESSOR[lev].pass[k],E->parallel.NUM_NEQ[lev].pass[k]);
 */


    }                           /* end for level */

    return;
}


// periodic BC
void parallel_communication_routs2(struct All_variables *E)
{
    //int i, ii, j, k, l, node, el, elt, lnode, jj, doff;
    int i, ii, j, k, node, lnode, jj, doff;
    int lev, elx, elz, ely, nno, nox, noz, noy, p, kkk, kk;
    int me, nprocz, nprocx, nprocy;

    me = E->parallel.me;
    nprocx = E->parallel.nprocx;
    nprocz = E->parallel.nprocz;
    nprocy = E->parallel.nprocy;

    for(lev = E->mesh.levmax; lev >= E->mesh.levmin; lev--)
    {
        elx = E->lmesh.ELX[lev];
        elz = E->lmesh.ELZ[lev];
        ely = E->lmesh.ELY[lev];
        nox = E->lmesh.NOX[lev];
        noy = E->lmesh.NOY[lev];
        noz = E->lmesh.NOZ[lev];
        nno = E->lmesh.NNO[lev];

        ii = 0;
        kkk = 0;
        if(E->mesh.nsd == 2)
        {

            for(i = 1; i <= 2; i++)
            {                   /* do the OZ boundaries */

                ii++;

                E->parallel.NUM_PASS[lev].bound[1][i] = 1;

                for(p = 1; p <= E->parallel.NUM_PASS[lev].bound[1][i]; p++)
                {

                    kkk++;
                    /* determine the pass ID for ii-th boundary and p-th pass */

                    /* for kkk th pass, get the # of data to be passed */
                    E->parallel.NUM_NEQ[lev].pass[1][i] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[1][i] : noy) * E->mesh.nsd;
                    E->parallel.NUM_NODE[lev].pass[1][i] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[1][i] : noy);

                    /* for kkk th pass, determine the target processor */

                    E->parallel.PROCESSOR[lev].pass[1][i] = me - ((i == 1) ? 1 : -1) * nprocz;
                    if(i == 1 && E->parallel.me_loc[1] == 0)
                        E->parallel.PROCESSOR[lev].pass[1][i] = me + (nprocx - 1) * nprocz;
                    else if(i == 2 && E->parallel.me_loc[1] == nprocx - 1)
                        E->parallel.PROCESSOR[lev].pass[1][i] = me - (nprocx - 1) * nprocz;

                    /* for kkk th pass, determine the ID of equations to be passed */

                    for(node = 1; node <= E->parallel.NUM_NODE[lev].pass[1][i]; node++)
                    {
                        lnode = E->parallel.NODE[lev][node].bound[1][i];
                        E->parallel.EXCHANGE_NODE[lev][node].pass[1][i] = lnode;
                        for(doff = 1; doff <= E->mesh.nsd; doff++)
                        {
                            jj = doff + (node - 1) * E->mesh.nsd;
                            E->parallel.EXCHANGE_ID[lev][jj].pass[1][i] = E->ID[lev][lnode].doff[doff];
                        }
                    }

                }               /* end for loop p */
            }                   /* end for loop i */


            for(j = 1; j <= 2; j++)
            {                   /* do OX boundary */

                ii++;

                E->parallel.NUM_PASS[lev].bound[2][j] = 1;
                if(E->parallel.me_loc[2] == 0 && j == 1)
                    E->parallel.NUM_PASS[lev].bound[2][j] = 0;
                else if(E->parallel.me_loc[2] == nprocz - 1 && j == 2)
                    E->parallel.NUM_PASS[lev].bound[2][j] = 0;

                for(p = 1; p <= E->parallel.NUM_PASS[lev].bound[2][j]; p++)
                {
                    kkk++;

                    /* for kkk th pass, get the # of data to be passed */
                    E->parallel.NUM_NEQ[lev].pass[2][j] = E->parallel.NUM_NNO[lev].bound[2][j] * E->mesh.nsd;
                    E->parallel.NUM_NODE[lev].pass[2][j] = E->parallel.NUM_NNO[lev].bound[2][j];

                    /* for kkk th pass, determine the target processor */

                    E->parallel.PROCESSOR[lev].pass[2][j] = me - ((j == 1) ? 1 : -1);

                    /* for kkk th pass, determine the ID of equations to be passed */

                    for(node = 1; node <= E->parallel.NUM_NODE[lev].pass[2][j]; node++)
                    {
                        lnode = E->parallel.NODE[lev][node].bound[2][j];
                        E->parallel.EXCHANGE_NODE[lev][node].pass[2][j] = lnode;
                        for(doff = 1; doff <= E->mesh.nsd; doff++)
                        {
                            jj = doff + (node - 1) * E->mesh.nsd;
                            E->parallel.EXCHANGE_ID[lev][jj].pass[2][j] = E->ID[lev][lnode].doff[doff];
                        }
                    }           /* end loop for node */

                }               /* end for loop p */
            }                   /* end for loop j */

        }                       /* end for 2D */



        else if(E->mesh.nsd == 3)
        {

            for(i = 1; i <= 2; i++)
            {                   /* do YOZ boundaries & OY lines */

                ii++;
                E->parallel.NUM_PASS[lev].bound[1][i] = 1;

                for(p = 1; p <= E->parallel.NUM_PASS[lev].bound[1][i]; p++)
                {
                    kkk++;
                    /* determine the pass ID for ii-th boundary and p-th pass */
                    E->parallel.NUM_NEQ[lev].pass[1][i] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[1][i] : noy) * E->mesh.nsd;
                    E->parallel.NUM_NODE[lev].pass[1][i] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[1][i] : noy);

                    E->parallel.PROCESSOR[lev].pass[1][i] = me - ((i == 1) ? 1 : -1) * nprocz;
                    if(i == 1 && E->parallel.me_loc[1] == 0)
                        E->parallel.PROCESSOR[lev].pass[1][i] = me + (nprocx - 1) * nprocz;
                    else if(i == 2 && E->parallel.me_loc[1] == nprocx - 1)
                        E->parallel.PROCESSOR[lev].pass[1][i] = me - (nprocx - 1) * nprocz;

                    for(k = 1; k <= E->parallel.NUM_NODE[lev].pass[1][i]; k++)
                    {
                        lnode = k;
                        node = E->parallel.NODE[lev][lnode].bound[1][i];
                        E->parallel.EXCHANGE_NODE[lev][k].pass[1][i] = node;
                        for(doff = 1; doff <= E->mesh.nsd; doff++)
                        {
                            jj = doff + (k - 1) * E->mesh.nsd;
                            E->parallel.EXCHANGE_ID[lev][jj].pass[1][i] = E->ID[lev][node].doff[doff];
                        }
                    }           /* end for node k */

                }               /* end for loop p */
            }                   /* end for i */

            for(j = 1; j <= 2; j++)
            {                   /* do XOY boundaries & OX lines */

                ii++;
                E->parallel.NUM_PASS[lev].bound[3][j] = 1;
                if(E->parallel.me_loc[3] == 0 && j == 1)
                    E->parallel.NUM_PASS[lev].bound[3][j] = 0;
                else if(E->parallel.me_loc[3] == nprocz - 1 && j == 2)
                    E->parallel.NUM_PASS[lev].bound[3][j] = 0;

                for(p = 1; p <= E->parallel.NUM_PASS[lev].bound[3][j]; p++)
                {
                    kkk++;
                    /* determine the pass ID for ii-th boundary and p-th pass */
                    E->parallel.NUM_NEQ[lev].pass[3][j] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[3][j] : nox) * E->mesh.nsd;
                    E->parallel.NUM_NODE[lev].pass[3][j] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[3][j] : nox);

                    E->parallel.PROCESSOR[lev].pass[3][j] = me - ((j == 1) ? 1 : -1);

                    for(k = 1; k <= E->parallel.NUM_NODE[lev].pass[3][j]; k++)
                    {
                        lnode = k;
                        node = E->parallel.NODE[lev][lnode].bound[3][j];
                        E->parallel.EXCHANGE_NODE[lev][k].pass[3][j] = node;
                        for(doff = 1; doff <= E->mesh.nsd; doff++)
                        {
                            jj = doff + (k - 1) * E->mesh.nsd;
                            E->parallel.EXCHANGE_ID[lev][jj].pass[3][j] = E->ID[lev][node].doff[doff];
                        }
                    }           /* end for node k */

                }               /* end for loop p */

            }                   /* end for j */

            for(k = 1; k <= 2; k++)
            {                   /* do XOZ boundaries & OZ lines */

                ii++;
                E->parallel.NUM_PASS[lev].bound[2][k] = 1;

                for(p = 1; p <= E->parallel.NUM_PASS[lev].bound[2][k]; p++)
                {

                    kkk++;
                    /* determine the pass ID for ii-th boundary and p-th pass */
                    E->parallel.NUM_NEQ[lev].pass[2][k] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[2][k] : noz) * E->mesh.nsd;
                    E->parallel.NUM_NODE[lev].pass[2][k] = ((p == 1) ? E->parallel.NUM_NNO[lev].bound[2][k] : noz);

                    E->parallel.PROCESSOR[lev].pass[2][k] = me - ((k == 1) ? 1 : -1) * nprocx * nprocz;
                    if(k == 1 && E->parallel.me_loc[2] == 0)
                        E->parallel.PROCESSOR[lev].pass[2][k] = me + (nprocy - 1) * nprocx * nprocz;
                    else if(k == 2 && E->parallel.me_loc[2] == nprocy - 1)
                        E->parallel.PROCESSOR[lev].pass[2][k] = me - (nprocy - 1) * nprocx * nprocz;

                    for(kk = 1; kk <= E->parallel.NUM_NODE[lev].pass[2][k]; kk++)
                    {
                        lnode = kk;
                        node = E->parallel.NODE[lev][lnode].bound[2][k];
                        E->parallel.EXCHANGE_NODE[lev][kk].pass[2][k] = node;
                        for(doff = 1; doff <= E->mesh.nsd; doff++)
                        {
                            jj = doff + (kk - 1) * E->mesh.nsd;
                            E->parallel.EXCHANGE_ID[lev][jj].pass[2][k] = E->ID[lev][node].doff[doff];
                        }
                    }           /* end for node kk */

                }               /* end for loop p */

            }                   /* end for k */

        }                       /* end of dims==3 */

        E->parallel.TNUM_PASS[lev] = kkk;


/*
    fprintf(E->fp," me= %d  pass  %d \n",E->parallel.me,E->parallel.TNUM_PASS[lev]);
    for (k=1;k<=E->parallel.TNUM_PASS[lev];k++)  
fprintf(E->fp,"proc %d and pass  %d to proc %d with %d eqn\n",E->parallel.me,k,E->parallel.PROCESSOR[lev].pass[k],E->parallel.NUM_NEQ[lev].pass[k]);
 */

    }                           /* end for level */

    return;
}

void parallel_communication_routs3(struct All_variables *E)
{
    //int i, ii, j, k, l, node, el, elt, lnode, jj, doff;
    int i, j, k, el;
    //int lev, elx, elz, ely, nno, nox, noz, noy, p, kkk, kk;
    int lev, elx, elz, ely;
    int m1, m2, m3, me, nprocz, nprocx, nprocy, proc;

    me = E->parallel.me;
    nprocx = E->parallel.nprocx;
    nprocz = E->parallel.nprocz;
    nprocy = E->parallel.nprocy;


    E->parallel.no_neighbors = 0;

    E->parallel.neighbors_rev = (int *)malloc((E->parallel.nproc + 1) * sizeof(int));

    for(k = -1; k <= 1; k++)
        for(i = -1; i <= 1; i++)
            for(j = -1; j <= 1; j++)
            {

                m1 = E->parallel.me_loc[1] + i;
                m2 = E->parallel.me_loc[2] + k;
                m3 = E->parallel.me_loc[3] + j;
                if(m1 >= 0 && m1 < nprocx && m2 >= 0 && m2 < nprocy && m3 >= 0 && m3 < nprocz)
                {
                    proc = m3 + m1 * nprocz + m2 * E->parallel.nprocxz;
                    if(proc != me)
                    {
                        E->parallel.no_neighbors++;
                        E->parallel.neighbors[E->parallel.no_neighbors] = proc;
                        E->parallel.neighbors_rev[proc] = E->parallel.no_neighbors;
                    }
                }
            }


    for(i = 1; i <= E->parallel.no_neighbors; i++)
        fprintf(E->fp, "aaa %d %d %d\n", i, E->parallel.neighbors[i], E->parallel.neighbors_rev[E->parallel.neighbors[i]]);
    fflush(E->fp);


    lev = E->mesh.levmax;

    elx = E->lmesh.ELX[lev];
    elz = E->lmesh.ELZ[lev];
    ely = E->lmesh.ELY[lev];

    for(el = 1; el <= E->lmesh.nel; el++)
        E->Element[el] = 0;

    for(k = 1; k <= ely; k++)
        for(j = 1; j <= elz; j++)
        {
            el = j + (k - 1) * elx * elz;
            E->Element[el] = E->Element[el] | SIDEE;
            el = j + (elx - 1) * elz + (k - 1) * elx * elz;
            E->Element[el] = E->Element[el] | SIDEE;
        }

    for(i = 1; i <= elx; i++)
        for(j = 1; j <= elz; j++)
        {
            el = j + (i - 1) * elz;
            E->Element[el] = E->Element[el] | SIDEE;
            el = j + (i - 1) * elz + (ely - 1) * elx * elz;
            E->Element[el] = E->Element[el] | SIDEE;
        }

    for(k = 1; k <= ely; k++)
        for(i = 1; i <= elx; i++)
        {
            el = elz + (i - 1) * elz + (k - 1) * elx * elz;
            E->Element[el] = E->Element[el] | SIDEE;
            el = 1 + (i - 1) * elz + (k - 1) * elx * elz;
            E->Element[el] = E->Element[el] | SIDEE;
        }


    return;
}

void parallel_communication_routs4(struct All_variables *E)
{
    //int i, ii, j, k, l, node, el, elt, lnode, jj, doff;
    //int lev, elx, elz, ely, nno, nox, noz, noy, p, kkk, kk;
    //int me, nprocz, nprocx, nprocy;

    return;
}


void exchange_number_rec_markers(struct All_variables *E)
{
    //int target_proc, kk, e, node, i, ii, j, k, bound, type, idb, msginfo[8];
    int target_proc, k, idb;
    static int *S[27], *R[27];

    static int been_here = 0;

    MPI_Status status[100];
    MPI_Request request[100];

    if(been_here == 0)
    {
        for(k = 1; k <= E->parallel.no_neighbors; k++)
        {
            S[k] = (int *)malloc(2 * sizeof(int));
            R[k] = (int *)malloc(2 * sizeof(int));
        }
        been_here++;
    }

    idb = 0;
    for(k = 1; k <= E->parallel.no_neighbors; k++)
    {
        target_proc = E->parallel.neighbors[k];
        idb++;
        S[k][0] = E->parallel.traces_transfer_number[k];
        S[k][1] = E->parallel.me;
        MPI_Isend(S[k], 2, MPI_INT, target_proc, E->parallel.mst1[E->parallel.me][target_proc], MPI_COMM_WORLD, &request[idb - 1]);
    }                           /* for k */

    for(k = 1; k <= E->parallel.no_neighbors; k++)
    {
        target_proc = E->parallel.neighbors[k];
        idb++;
        MPI_Irecv(R[k], 2, MPI_INT, target_proc, E->parallel.mst1[E->parallel.me][target_proc], MPI_COMM_WORLD, &request[idb - 1]);
    }                           /* for k */

    MPI_Waitall(idb, request, status);

    for(k = 1; k <= E->parallel.no_neighbors; k++)
        E->parallel.traces_receive_number[k] = R[k][0];

    return;
}



/* ============================================ 
  exchange ID related information across boundaries
 ============================================ */


void exchange_markers(struct All_variables *E)
{
    //int target_proc, kk, e, node, i, ii, j, k, bound, type, idb, msginfo[8];
    int target_proc, kk, k, idb;

    //static int been_here = 0;
    //static int mid_recv, sizeofk;
    //const int levmax = E->mesh.levmax;

    MPI_Status status[100];
    MPI_Request request[100];

    idb = 0;
    for(k = 1; k <= E->parallel.no_neighbors; k++)
    {
        if(E->parallel.traces_transfer_number[k] > 0)
        {
            target_proc = E->parallel.neighbors[k];
            idb++;
            kk = E->parallel.traces_transfer_number[k] * 2 + 1;
            MPI_Isend(E->PINS[k], kk, MPI_INT, target_proc, E->parallel.mst1[E->parallel.me][target_proc], MPI_COMM_WORLD, &request[idb - 1]);
            idb++;
            kk = E->parallel.traces_transfer_number[k] * 2 * E->mesh.nsd + 1;
            MPI_Isend(E->PVV[k], kk, MPI_FLOAT, target_proc, E->parallel.mst2[E->parallel.me][target_proc], MPI_COMM_WORLD, &request[idb - 1]);
            idb++;
            MPI_Isend(E->PXX[k], kk, MPI_DOUBLE, target_proc, E->parallel.mst3[E->parallel.me][target_proc], MPI_COMM_WORLD, &request[idb - 1]);
        }
    }                           /* for k */

    for(k = 1; k <= E->parallel.no_neighbors; k++)
    {
        if(E->parallel.traces_receive_number[k] > 0)
        {
            target_proc = E->parallel.neighbors[k];
            idb++;
            kk = E->parallel.traces_receive_number[k] * 2 + 1;
            MPI_Irecv(E->RINS[k], kk, MPI_INT, target_proc, E->parallel.mst1[E->parallel.me][target_proc], MPI_COMM_WORLD, &request[idb - 1]);
            idb++;
            kk = E->parallel.traces_receive_number[k] * 2 * E->mesh.nsd + 1;
            MPI_Irecv(E->RVV[k], kk, MPI_FLOAT, target_proc, E->parallel.mst2[E->parallel.me][target_proc], MPI_COMM_WORLD, &request[idb - 1]);
            idb++;
            MPI_Irecv(E->RXX[k], kk, MPI_DOUBLE, target_proc, E->parallel.mst3[E->parallel.me][target_proc], MPI_COMM_WORLD, &request[idb - 1]);
        }
    }                           /* for k */

    MPI_Waitall(idb, request, status);


    return;
}



/* ============================================ 
  exchange ID related information across boundaries
 ============================================ */

void exchange_id_d20(struct All_variables *E, double *U, int lev)
{
    //int target_proc, kk, e, node, i, ii, j, k, bound, type, idb, msginfo[8];
    int target_proc, kk, i, j, k, idb;
    static double *S[27], *R[27];

    static int been_here = 0;
    //static int mid_recv, sizeofk;
    static int sizeofk;
    const int levmax = E->mesh.levmax;

    MPI_Status status[100];
    MPI_Request request[100];

    if(been_here == 0)
    {
        sizeofk = 0;
        for(i = 1; i <= E->mesh.nsd; i++)
            for(k = 1; k <= 2; k++)
            {
                if(E->parallel.NUM_PASS[levmax].bound[i][k])
                    sizeofk = max(sizeofk, (1 + E->parallel.NUM_NEQ[levmax].pass[i][k]) * sizeof(double));
            }

        for(k = 1; k <= 2; k++)
        {
            S[k] = (double *)malloc(sizeofk);
            R[k] = (double *)malloc(sizeofk);
        }
        been_here++;
    }

    for(i = 1; i <= E->mesh.nsd; i++)
    {
        idb = 0;
        for(k = 1; k <= 2; k++)
            if(E->parallel.NUM_PASS[levmax].bound[i][k])
            {
                for(j = 1; j <= E->parallel.NUM_NEQ[lev].pass[i][k]; j++)
                    S[k][j - 1] = U[E->parallel.EXCHANGE_ID[lev][j].pass[i][k]];

                target_proc = E->parallel.PROCESSOR[lev].pass[i][k];
                if(target_proc != E->parallel.me)
                {
                    idb++;
                    MPI_Isend(S[k], E->parallel.NUM_NEQ[lev].pass[i][k], MPI_DOUBLE, target_proc, E->parallel.mst[E->parallel.me][target_proc][k], MPI_COMM_WORLD, &request[idb - 1]);
                }
            }                   /* for k */

        for(k = 1; k <= 2; k++)
            if(E->parallel.NUM_PASS[levmax].bound[i][k])
            {

                target_proc = E->parallel.PROCESSOR[lev].pass[i][k];
                if(target_proc != E->parallel.me)
                {
                    idb++;
                    MPI_Irecv(R[k], E->parallel.NUM_NEQ[lev].pass[i][k], MPI_DOUBLE, target_proc, E->parallel.mst[E->parallel.me][target_proc][k], MPI_COMM_WORLD, &request[idb - 1]);
                }
            }                   /* for k */

        MPI_Waitall(idb, request, status);

        for(k = 1; k <= 2; k++)
            if(E->parallel.NUM_PASS[levmax].bound[i][k])
            {
                target_proc = E->parallel.PROCESSOR[lev].pass[i][k];
                if(target_proc != E->parallel.me)
                {
                    for(j = 1; j <= E->parallel.NUM_NEQ[lev].pass[i][k]; j++)
                        U[E->parallel.EXCHANGE_ID[lev][j].pass[i][k]] += R[k][j - 1];
                }
                else
                {
                    kk = 1;
                    if(k == 1)
                        kk = 2;
                    for(j = 1; j <= E->parallel.NUM_NEQ[lev].pass[i][k]; j++)
                        U[E->parallel.EXCHANGE_ID[lev][j].pass[i][k]] += S[kk][j - 1];
                }
            }                   /* for k */

    }                           /* for dim */

    return;
}


void exchange_node_f20(struct All_variables *E, float *U, int lev)
{
    //int target_proc, kk, e, node, i, ii, j, k, bound, type, idb, msginfo[8];
    int target_proc, kk, i, j, k, idb;
    static float *S[27], *R[27];

    static int been_here = 0;
    //static int mid_recv, sizeofk;
    static int sizeofk;
    const int levmax = E->mesh.levmax;

    MPI_Status status[100];
    MPI_Request request[100];

    if(been_here == 0)
    {
        sizeofk = 0;
        for(i = 1; i <= E->mesh.nsd; i++)
            for(k = 1; k <= 2; k++)
            {
                if(E->parallel.NUM_PASS[levmax].bound[i][k])
                    sizeofk = max(sizeofk, (1 + E->parallel.NUM_NODE[levmax].pass[i][k]) * sizeof(float));
            }

        for(k = 1; k <= 2; k++)
        {
            S[k] = (float *)malloc(sizeofk);
            R[k] = (float *)malloc(sizeofk);
        }
        been_here++;
    }

    for(i = 1; i <= E->mesh.nsd; i++)
    {
        idb = 0;
        for(k = 1; k <= 2; k++)
            if(E->parallel.NUM_PASS[levmax].bound[i][k])
            {
                for(j = 1; j <= E->parallel.NUM_NODE[lev].pass[i][k]; j++)
                    S[k][j - 1] = U[E->parallel.EXCHANGE_NODE[lev][j].pass[i][k]];

                target_proc = E->parallel.PROCESSOR[lev].pass[i][k];
                if(target_proc != E->parallel.me)
                {
                    idb++;
                    MPI_Isend(S[k], E->parallel.NUM_NODE[lev].pass[i][k], MPI_FLOAT, target_proc, E->parallel.mst[E->parallel.me][target_proc][k], MPI_COMM_WORLD, &request[idb - 1]);
                }
            }                   /* for k */

        for(k = 1; k <= 2; k++)
            if(E->parallel.NUM_PASS[levmax].bound[i][k])
            {

                target_proc = E->parallel.PROCESSOR[lev].pass[i][k];
                if(target_proc != E->parallel.me)
                {
                    idb++;
                    MPI_Irecv(R[k], E->parallel.NUM_NODE[lev].pass[i][k], MPI_FLOAT, target_proc, E->parallel.mst[E->parallel.me][target_proc][k], MPI_COMM_WORLD, &request[idb - 1]);
                }
            }                   /* for k */

        MPI_Waitall(idb, request, status);

        for(k = 1; k <= 2; k++)
            if(E->parallel.NUM_PASS[levmax].bound[i][k])
            {
                target_proc = E->parallel.PROCESSOR[lev].pass[i][k];
                if(target_proc != E->parallel.me)
                {
                    for(j = 1; j <= E->parallel.NUM_NODE[lev].pass[i][k]; j++)
                        U[E->parallel.EXCHANGE_NODE[lev][j].pass[i][k]] += R[k][j - 1];
                }
                else
                {
                    kk = 1;
                    if(k == 1)
                        kk = 2;
                    for(j = 1; j <= E->parallel.NUM_NODE[lev].pass[i][k]; j++)
                        U[E->parallel.EXCHANGE_NODE[lev][j].pass[i][k]] += S[kk][j - 1];
                }
            }                   /* for k */

    }                           /* for dim */

    return;
}

/* ==========================   */

double CPU_time0(void)
{
    return (double)MPI_Wtime();
}


void parallel_process_sync(void)
{
    MPI_Barrier(MPI_COMM_WORLD);
}
