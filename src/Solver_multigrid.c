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

#include "element_definitions.h"
#include "global_defs.h"
#include <math.h>


/* *INDENT-OFF* */
static double weightIJ[5][5] = {  {0.0, 0.0   , 0.0   , 0.0   , 0.0   },
                                  {0.0, 0.5625, 0.1875, 0.0625, 0.1875},
                                  {0.0, 0.1875, 0.5625, 0.1875, 0.0625},
                                  {0.0, 0.0625, 0.1875, 0.5625, 0.1875},
                                  {0.0, 0.1875, 0.0625, 0.1875, 0.5625}  };
/* *INDENT-ON* */


void set_mg_defaults(struct All_variables *E)
{
    E->solver_allocate_vars = mg_allocate_vars;
    E->build_forcing_term = assemble_forces_iterative;
    E->solve_stokes_problem = solve_constrained_flow_iterative;

    E->control.mg_cycle = 1;

    return;
}

void mg_allocate_vars(struct All_variables *E)
{
    return;
}


/* project_vector
 *	double *AU - data on upper/lower mesh
 *	double *AD - data on upper/lower mesh
 */
void project_vector(struct All_variables *E, int start_lev, double *AU, double *AD, int ic)
{
    int i, j;
    int el, node, e1;
    //int eqn1, eqn_minus1;
    //int eqn2, eqn_minus2;
    //int eqn3, eqn_minus3;
    //double amplifier, average1, average2, average3, w;
    double average1, average2, average3;
    //float time;

    const int sl_minus = start_lev - 1;
    const int neq_minus = E->lmesh.NEQ[start_lev - 1];
    const int nno_minus = E->lmesh.NNO[start_lev - 1];
    const int nels_minus = E->lmesh.NEL[start_lev - 1];
    const int dims = E->mesh.nsd;
    //const int ends = enodes[E->mesh.nsd];
    //const double weight = (double)1.0 / ends;


    /* on the lower level the average value of data in upper level
     * ELEMENTS are slid across to the nodes.  */

    for(i = 0; i <= neq_minus; i++)
        AD[i] = 0.0;

    if(3 == dims)
        for(el = 1; el <= nels_minus; el++)
            for(i = 1; i <= ENODES3D; i++)
            {
                average1 = average2 = average3 = 0.0;
                e1 = E->EL[sl_minus][el].sub[i];
                for(j = 1; j <= ENODES3D; j++)
                {
                    node = E->IEN[start_lev][e1].node[j];
                    average1 += AU[E->ID[start_lev][node].doff[1]];
                    average2 += AU[E->ID[start_lev][node].doff[2]];
                    average3 += AU[E->ID[start_lev][node].doff[3]];

                }
                node = E->IEN[sl_minus][el].node[i];

                AD[E->ID[sl_minus][node].doff[1]] += E->TWW[sl_minus][el].node[i] * average1;
                AD[E->ID[sl_minus][node].doff[2]] += E->TWW[sl_minus][el].node[i] * average2;
                AD[E->ID[sl_minus][node].doff[3]] += E->TWW[sl_minus][el].node[i] * average3;

            }
    else
        for(el = 1; el <= nels_minus; el++)
            for(i = 1; i <= ENODES2D; i++)
            {
                average1 = average2 = 0.0;
                e1 = E->EL[sl_minus][el].sub[i];
                for(j = 1; j <= ENODES2D; j++)
                {
                    node = E->IEN[start_lev][e1].node[j];
                    average1 += AU[E->ID[start_lev][node].doff[1]];
                    average2 += AU[E->ID[start_lev][node].doff[2]];
                }
                node = E->IEN[sl_minus][el].node[i];
                AD[E->ID[sl_minus][node].doff[1]] += E->TWW[sl_minus][el].node[i] * average1;
                AD[E->ID[sl_minus][node].doff[2]] += E->TWW[sl_minus][el].node[i] * average2;
            }

    exchange_id_d20(E, AD, sl_minus);

    if(3 == dims)
        for(i = 1; i <= nno_minus; i++)
        {
            AD[E->ID[sl_minus][i].doff[1]] = AD[E->ID[sl_minus][i].doff[1]] * E->MASS[sl_minus][i];
            AD[E->ID[sl_minus][i].doff[2]] = AD[E->ID[sl_minus][i].doff[2]] * E->MASS[sl_minus][i];
            AD[E->ID[sl_minus][i].doff[3]] = AD[E->ID[sl_minus][i].doff[3]] * E->MASS[sl_minus][i];
        }
    else
        for(i = 1; i <= nno_minus; i++)
        {
            AD[E->ID[sl_minus][i].doff[1]] = AD[E->ID[sl_minus][i].doff[1]] * E->MASS[sl_minus][i];
            AD[E->ID[sl_minus][i].doff[2]] = AD[E->ID[sl_minus][i].doff[2]] * E->MASS[sl_minus][i];
        }


    return;
}



/* =======================================================================================
   Interpolation from coarse grid to fine. See the appology attached to project() if you get
   stressed out by node based assumptions. If it makes you feel any better, I don't like
   it much either.
   ======================================================================================= */

/* interp_vector
 *	double *AD - data on upper/lower mesh
 *	double *AU - data on upper/lower mesh
 */
void interp_vector(struct All_variables *E, int start_lev, double *AD, double *AU)
{
    int i, j, k;
    float x1, x2;
    float n1, n2;
    int node0, node1, node2;
    int eqn0, eqn1, eqn2;

    const int level = start_lev + 1;
    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];

    const int nox = E->lmesh.NOX[level];
    const int noz = E->lmesh.NOZ[level];
    const int noy = E->lmesh.NOY[level];
    //const int high_eqn = E->lmesh.NEQ[level];

    if(start_lev == E->mesh.levmax)
        return;

    un_inject_vector(E, start_lev, AD, AU); /* transfer information from lower level */

    for(k = 1; k <= noy; k += 2)    /* Fill in gaps in x direction */
        for(j = 1; j <= noz; j += 2)
            for(i = 2; i < nox; i += 2)
            {
                node0 = j + (i - 1) * noz + (k - 1) * nox * noz;    /* this node */
                node1 = node0 - noz;
                node2 = node0 + noz;

                x1 = E->ECO[level][E->NEI[level].element[ends * (node1 - 1)]].size[1];
                x2 = E->ECO[level][E->NEI[level].element[(node2 - 1) * ends]].size[1];

                n1 = x2 / (x1 + x2);
                n2 = x1 / (x1 + x2);

                /* now for each direction */

                eqn0 = E->ID[level][node0].doff[1];
                eqn1 = E->ID[level][node1].doff[1];
                eqn2 = E->ID[level][node2].doff[1];
                AU[eqn0] = n1 * AU[eqn1] + n2 * AU[eqn2];

                eqn0 = E->ID[level][node0].doff[2];
                eqn1 = E->ID[level][node1].doff[2];
                eqn2 = E->ID[level][node2].doff[2];
                AU[eqn0] = n1 * AU[eqn1] + n2 * AU[eqn2];

                if(3 == dims)
                {
                    eqn0 = E->ID[level][node0].doff[3];
                    eqn1 = E->ID[level][node1].doff[3];
                    eqn2 = E->ID[level][node2].doff[3];
                    AU[eqn0] = n1 * AU[eqn1] + n2 * AU[eqn2];
                }

            }


    for(k = 1; k <= noy; k += 2)    /* Fill in gaps in z direction */
        for(i = 1; i <= nox; i++)
            for(j = 2; j < noz; j += 2)
            {
                node0 = j + (i - 1) * noz + (k - 1) * nox * noz;    /* this node */
                node1 = node0 - 1;
                node2 = node0 + 1;

                x1 = E->ECO[level][E->NEI[level].element[ends * (node1 - 1)]].size[3];
                x2 = E->ECO[level][E->NEI[level].element[(node2 - 1) * ends]].size[3];

                n1 = x2 / (x1 + x2);
                n2 = x1 / (x1 + x2);

                eqn0 = E->ID[level][node0].doff[1];
                eqn1 = E->ID[level][node1].doff[1];
                eqn2 = E->ID[level][node2].doff[1];
                AU[eqn0] = n1 * AU[eqn1] + n2 * AU[eqn2];

                eqn0 = E->ID[level][node0].doff[2];
                eqn1 = E->ID[level][node1].doff[2];
                eqn2 = E->ID[level][node2].doff[2];
                AU[eqn0] = n1 * AU[eqn1] + n2 * AU[eqn2];

                if(3 == dims)
                {
                    eqn0 = E->ID[level][node0].doff[3];
                    eqn1 = E->ID[level][node1].doff[3];
                    eqn2 = E->ID[level][node2].doff[3];
                    AU[eqn0] = n1 * AU[eqn1] + n2 * AU[eqn2];
                }
            }


    if(3 == dims)
        for(i = 1; i <= nox; i++)   /* Fill in gaps in y direction */
            for(j = 1; j <= noz; j++)
                for(k = 2; k < noy; k += 2)
                {
                    node0 = j + (i - 1) * noz + (k - 1) * nox * noz;    /* this node */
                    node1 = node0 - nox * noz;
                    node2 = node0 + nox * noz;

                    x1 = E->ECO[level][E->NEI[level].element[ends * (node1 - 1)]].size[2];
                    x2 = E->ECO[level][E->NEI[level].element[(node2 - 1) * ends]].size[2];

                    n1 = x2 / (x1 + x2);
                    n2 = x1 / (x1 + x2);

                    eqn0 = E->ID[level][node0].doff[1];
                    eqn1 = E->ID[level][node1].doff[1];
                    eqn2 = E->ID[level][node2].doff[1];
                    AU[eqn0] = n1 * AU[eqn1] + n2 * AU[eqn2];

                    eqn0 = E->ID[level][node0].doff[2];
                    eqn1 = E->ID[level][node1].doff[2];
                    eqn2 = E->ID[level][node2].doff[2];
                    AU[eqn0] = n1 * AU[eqn1] + n2 * AU[eqn2];

                    eqn0 = E->ID[level][node0].doff[3];
                    eqn1 = E->ID[level][node1].doff[3];
                    eqn2 = E->ID[level][node2].doff[3];
                    AU[eqn0] = n1 * AU[eqn1] + n2 * AU[eqn2];
                }
    return;

}

/* ==================================================== */

/* project_scalar_e
 * 	float *AU - data on upper/lower mesh
 * 	float *AD - data on upper/lower mesh
 */
void project_scalar_e(struct All_variables *E, int start_lev, float *AU, float *AD)
{
    //int i, j, m;
    int i;
    //int el, node, e;
    int el, e;
    //float average, w;
    float average;

    //static int been_here = 0;

    const int sl_minus = start_lev - 1;
    const int nels_minus = E->lmesh.NEL[start_lev - 1];
    //const int dims = E->mesh.nsd;
    const int ends = enodes[E->mesh.nsd];
    const double weight = (double)1.0 / ends;
    //const int vpts = vpoints[E->mesh.nsd];
    //const int n_minus = nels_minus * vpts;

    for(el = 1; el <= nels_minus; el++)
    {
        average = 0.0;
        for(i = 1; i <= ends; i++)
        {
            e = E->EL[sl_minus][el].sub[i];
            average += AU[e];
        }

        AD[el] = average * weight;
    }

    return;
}

/* ==================================================== */
/* project_scalar
 * 	float *AU - data on upper/lower mesh
 * 	float *AD - data on upper/lower mesh
 */
void project_scalar(struct All_variables *E, int start_lev, float *AU, float *AD)
{
    //int i, j, m;
    int i, j;
    int el, node, node1;
    float average, w;

    const int sl_minus = start_lev - 1;
    const int nno_minus = E->lmesh.NNO[start_lev - 1];
    const int nels_minus = E->lmesh.NEL[start_lev - 1];
    //const int dims = E->mesh.nsd;
    const int ends = enodes[E->mesh.nsd];
    const double weight = (double)1.0 / ends;

    for(i = 1; i <= nno_minus; i++)
        AD[i] = 0.0;

    for(el = 1; el <= nels_minus; el++)
        for(i = 1; i <= ends; i++)
        {
            average = 0.0;
            node1 = E->EL[sl_minus][el].sub[i];
            for(j = 1; j <= ends; j++)
            {
                node = E->IEN[start_lev][node1].node[j];
                average += AU[node];
            }

            w = weight * average;

            node = E->IEN[sl_minus][el].node[i];

            AD[node] += w * E->TWW[sl_minus][el].node[i];
        }

    exchange_node_f20(E, AD, sl_minus);

    for(i = 1; i <= nno_minus; i++)
    {
        AD[i] *= E->MASS[sl_minus][i];
    }


    return;
}




/*  ==============================================
    function to project viscosity down to all the 
    levels in the problem. (no gaps for vbcs)
    ==============================================  */

void project_viscosity(struct All_variables *E)
{
    //int lv, i, j, k, el, sl_minus;
    int lv, sl_minus;

    const int nsd = E->mesh.nsd;
    const int vpts = vpoints[nsd];

    float *viscU, *viscD;

    viscU = (float *)malloc((1 + E->lmesh.NNO[E->mesh.levmax]) * sizeof(float));
    viscD = (float *)malloc((1 + vpts * E->lmesh.NNO[E->mesh.levmax - 1]) * sizeof(float));

    for(lv = E->mesh.levmax; lv > E->mesh.levmin; lv--)
    {
        sl_minus = lv - 1;
        if(E->viscosity.smooth_cycles == 1)
        {
            visc_from_gint_to_nodes(E, E->EVI[lv], viscU, lv);
            project_scalar(E, lv, viscU, viscD);
            visc_from_nodes_to_gint(E, viscD, E->EVI[sl_minus], sl_minus);
        }
        else if(E->viscosity.smooth_cycles == 2)
        {
            visc_from_gint_to_ele(E, E->EVI[lv], viscU, lv);
            inject_scalar_e(E, lv, viscU, E->EVI[sl_minus]);
        }
        else if(E->viscosity.smooth_cycles == 3)
        {
            visc_from_gint_to_ele(E, E->EVI[lv], viscU, lv);
            project_scalar_e(E, lv, viscU, viscD);
            visc_from_ele_to_gint(E, viscD, E->EVI[sl_minus], sl_minus);
        }
        else if(E->viscosity.smooth_cycles == 0)
        {
            visc_from_gint_to_nodes(E, E->EVI[lv], viscU, lv);
            inject_scalar(E, lv, viscU, viscD);
            visc_from_nodes_to_gint(E, viscD, E->EVI[sl_minus], sl_minus);
        }
    }

    free((void *)viscU);
    free((void *)viscD);
    return;
}


/* inject_node_fvector
 *	float **AU - data on upper/lower mesh
 *	float **AD - data on upper/lower mesh
 */
void inject_node_fvector(struct All_variables *E, int start_lev, float **AU, float **AD)
{
    int i;
    //int el, ex, ey, ez, d;
    int el, ex, ey, ez;
    int node, node_minus;

    const int sl_minus = start_lev - 1;
    const int elx = E->lmesh.ELX[sl_minus];
    const int elz = E->lmesh.ELZ[sl_minus];
    const int ely = E->lmesh.ELY[sl_minus];
    const int dims = E->mesh.nsd;

    assert(start_lev != E->mesh.levmin);

    for(ey = 1; ey <= ely; ey++)
        for(ez = 1; ez <= elz; ez++)
            for(ex = 1; ex <= elx; ex++)
            {
                el = ez + (ex - 1) * elz + (ey - 1) * elz * elx;
                for(i = 1; i <= enodes[dims]; i++)
                {
                    node_minus = E->IEN[sl_minus][el].node[i];
                    node = E->IEN[start_lev][E->EL[sl_minus][el].sub[i]].node[i];

                    AD[1][node_minus] = AU[1][node];
                    AD[2][node_minus] = AU[2][node];
                    if(3 == dims)
                        AD[3][node_minus] = AU[3][node];
                }
            }
    return;
}

/* =====================================================
   Function to inject data from high to low grid (i.e.
   just dropping values not at shared grid points.
   ===================================================== */

/* inject
 * 	double *AU - data on upper/lower mesh
 * 	double *AD - data on upper/lower mesh
 */
void inject(struct All_variables *E, int start_lev, double *AU, double *AD)
{
    int i;
    int el, node_coarse, node_fine;
    int sl_minus;
    int eqn, eqn_coarse;

    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];

    if(start_lev == E->mesh.levmin)
    {
        fprintf(E->fp, "Warning, attempting to project below lowest level\n");
        return;
    }
    sl_minus = start_lev - 1;

    for(el = 1; el <= E->lmesh.NEL[sl_minus]; el++)
        for(i = 1; i <= ends; i++)
        {
            node_coarse = E->IEN[sl_minus][el].node[i];
            node_fine = E->IEN[start_lev][E->EL[sl_minus][el].sub[i]].node[i];

            eqn_coarse = E->ID[sl_minus][node_coarse].doff[1];
            eqn = E->ID[start_lev][node_fine].doff[1];
            AD[eqn_coarse] = AU[eqn];

            eqn_coarse = E->ID[sl_minus][node_coarse].doff[2];
            eqn = E->ID[start_lev][node_fine].doff[2];
            AD[eqn_coarse] = AU[eqn];

            if(3 == dims)
            {
                eqn_coarse = E->ID[sl_minus][node_coarse].doff[3];
                eqn = E->ID[start_lev][node_fine].doff[3];
                AD[eqn_coarse] = AU[eqn];
            }
        }
    return;
}


/* =====================================================
   Function to inject data from high to low grid (i.e.
   just dropping values not at shared grid points.
   ===================================================== */

/* un_inject_vector
 * 	double *AD - data on upper/lower mesh
 * 	double *AU - data on upper/lower mesh
 */
void un_inject_vector(struct All_variables *E, int start_lev, double *AD, double *AU)
{
    int i;
    int el, node, node_plus;
    int eqn1, eqn_plus1;
    int eqn2, eqn_plus2;
    int eqn3, eqn_plus3;

    const int dims = E->mesh.nsd;
    //const int ends = enodes[dims];
    const int sl_plus = start_lev + 1;
    const int neq = E->lmesh.NEQ[sl_plus];
    const int nels = E->lmesh.NEL[start_lev];

    assert(start_lev != E->mesh.levmax /* un_injection */ );

    for(i = 1; i <= neq; i++)
        AU[i] = 0.0;

    if(3 == dims)
        for(el = 1; el <= nels; el++)
            for(i = 1; i <= ENODES3D; i++)
            {
                node = E->IEN[start_lev][el].node[i];
                node_plus = E->IEN[sl_plus][E->EL[start_lev][el].sub[i]].node[i];

                eqn1 = E->ID[start_lev][node].doff[1];
                eqn2 = E->ID[start_lev][node].doff[2];
                eqn3 = E->ID[start_lev][node].doff[3];
                eqn_plus1 = E->ID[sl_plus][node_plus].doff[1];
                eqn_plus2 = E->ID[sl_plus][node_plus].doff[2];
                eqn_plus3 = E->ID[sl_plus][node_plus].doff[3];
                AU[eqn_plus1] = AD[eqn1];
                AU[eqn_plus2] = AD[eqn2];
                AU[eqn_plus3] = AD[eqn3];

            }
    else
        for(el = 1; el <= nels; el++)
            for(i = 1; i <= ENODES2D; i++)
            {
                node = E->IEN[start_lev][el].node[i];
                node_plus = E->IEN[sl_plus][E->EL[start_lev][el].sub[i]].node[i];

                eqn1 = E->ID[start_lev][node].doff[1];
                eqn2 = E->ID[start_lev][node].doff[2];
                eqn_plus1 = E->ID[sl_plus][node_plus].doff[1];
                eqn_plus2 = E->ID[sl_plus][node_plus].doff[2];
                AU[eqn_plus1] = AD[eqn1];
                AU[eqn_plus2] = AD[eqn2];
            }

    return;
}


/* inject_scalar
 * 	float *AU - data on upper/lower mesh
 * 	float *AD - data on upper/lower mesh
 */
void inject_scalar(struct All_variables *E, int start_lev, float *AU, float *AD)
{
    //int i, m, el, node_coarse, node_fine, sl_minus, eqn, eqn_coarse;
    int i, el, node_coarse, node_fine, sl_minus;

    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];

    if(start_lev == E->mesh.levmin)
    {
        fprintf(E->fp, "Warning, attempting to project below lowest level\n");
        return;
    }

    sl_minus = start_lev - 1;

    for(el = 1; el <= E->lmesh.NEL[sl_minus]; el++)
    {
        for(i = 1; i <= ends; i++)
        {
            node_coarse = E->IEN[sl_minus][el].node[i];
            node_fine = E->IEN[start_lev][E->EL[sl_minus][el].sub[i]].node[i];
            AD[node_coarse] = AU[node_fine];
        }
    }

    return;
}


/* inject_scalar_e
 * 	float *AU - data on upper/lower mesh
 * 	float *AD - data on upper/lower mesh
 */
void inject_scalar_e(struct All_variables *E, int start_lev, float *AU, float *AD)
{
    //int i, j, m;
    int i;
    //int el, node, e;
    int el, e;
    //float average, w;

    const int sl_minus = start_lev - 1;
    const int nels_minus = E->lmesh.NEL[start_lev - 1];
    //const int dims = E->mesh.nsd;
    const int ends = enodes[E->mesh.nsd];
    const int vpts = vpoints[E->mesh.nsd];

    for(el = 1; el <= nels_minus; el++)
        for(i = 1; i <= ends; i++)
        {
            e = E->EL[sl_minus][el].sub[i];
            AD[(el - 1) * vpts + i] = AU[e];
        }
    return;
}
