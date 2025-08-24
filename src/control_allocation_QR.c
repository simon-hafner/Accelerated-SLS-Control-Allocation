/*------------------------------------------------------------------------*/
/*                       control_allocation_QR.c                         */
/*------------------------------------------------------------------------*/
/*                                                                        */
/*  I N S T I T U T E   O F   F L I G H T   S Y S T E M   D Y N A M I C S */
/*                       Web: www.fsd.ed.tum.de                           */
/*                         _______                                        */
/*                            |  |    |\  /|                              */
/*                            |  |    | \/ |                              */
/*                            |  |____|    |                              */
/*                 Technische Universitaet Muenchen TUM                   */
/*                                                                        */
/*(c) 2023 by Institute of Flight System Dynamics                         */
/*                        All Rights Reserved                             */
/*------------------------------------------------------------------------*/
/*Description:  Provides control allocation functionality                 */
/*Type:         C - source file                                           */
/*Dependencies:                                                           */
/*------------------------------------------------------------------------*/
/*Author:       S. Hafner                                                 */
/*Date:         2023-11-08                                                */
/*Contact:      Institute of Flight System Dynamics                       */
/*              School of Engineering and Design                          */
/*              Technische Universitaet Muenchen                          */
/*              Boltzmannstrasse 15                                       */
/*              D-85748 Garching                                          */
/*              simon.hafner@tum.de                                       */
/*------------------------------------------------------------------------*/

#include <math.h>
#include "control_allocation_QR.h"

/* Static functions to perform the QR decomposition  */

static unsigned int H1(unsigned short lpivot, unsigned short l1, unsigned char transpose, unsigned short m, matrix_t *const U, float *const up, short C_max)
{
    unsigned int status = 0;
    if ((lpivot < l1) && (l1 < m))
    {
        if (transpose == 1U)
        {
            float cl = ABS(M_P_ENTRY(U, lpivot, lpivot));
            for (unsigned short j = l1; j < m; j++)
            {
                /* cl = max(abs(U(j)), cl) */
                if (ABS(M_P_ENTRY(U, lpivot, j)) > cl)
                {
                    cl = ABS(M_P_ENTRY(U, lpivot, j));
                }
            }
            /* cl != 0*/
            if (ABS(cl) > (EPS_F * 0.01f))
            {
                float clinv = 1.0f / cl;

                float sm = (M_P_ENTRY(U, lpivot, lpivot) * clinv) * (M_P_ENTRY(U, lpivot, lpivot) * clinv);

                for (unsigned short j = l1; j < m; j++)
                {
                    sm += (M_P_ENTRY(U, lpivot, j) * clinv) * (M_P_ENTRY(U, lpivot, j) * clinv);
                }

                cl = cl * sqrtf(sm);

                if (ABS(M_P_ENTRY(U, lpivot, lpivot)) <= EPS_F)
                {
                    cl = -cl;
                }

                *up = M_P_ENTRY(U, lpivot, lpivot) - cl;
                M_P_ENTRY(U, lpivot, lpivot) = cl;
            }
            float b = *up * M_P_ENTRY(U, lpivot, lpivot);

            if ((ABS(b) > (EPS_F * 0.01f)) && (C_max >= 0))
            {
                b = 1.0f / b;

                for (unsigned short j = l1; j < (unsigned short)C_max; j++)
                {
                    float sm = M_P_ENTRY(U, j, lpivot) * *up;

                    for (unsigned short i = l1; i < m; i++)
                    {
                        sm += M_P_ENTRY(U, j, i) * M_P_ENTRY(U, lpivot, i);
                    }

                    if (ABS(sm) > (EPS_F * 0.01f))
                    {
                        sm = sm * b;

                        M_P_ENTRY(U, j, lpivot) += sm * *up;

                        for (unsigned short i = l1; i < m; i++)
                        {
                            M_P_ENTRY(U, j, i) += sm * M_P_ENTRY(U, lpivot, i);
                        }
                    }
                }
            }
        }
        else
        {
            float cl = ABS(M_P_ENTRY(U, lpivot, lpivot));
            for (unsigned short j = l1; j < m; j++)
            {
                /* cl = max(abs(U(j)), cl) */
                if (ABS(M_P_ENTRY(U, j, lpivot)) > cl)
                {
                    cl = ABS(M_P_ENTRY(U, j, lpivot));
                }
            }
            /* cl != 0*/
            if (ABS(cl) > EPS_F)
            {
                float clinv = 1.0f / cl;

                float sm = (M_P_ENTRY(U, lpivot, lpivot) * clinv) * (M_P_ENTRY(U, lpivot, lpivot) * clinv);

                for (unsigned short j = l1; j < m; j++)
                {
                    sm += (M_P_ENTRY(U, j, lpivot) * clinv) * (M_P_ENTRY(U, j, lpivot) * clinv);
                }

                cl = cl * sqrtf(sm);

                if (ABS(M_P_ENTRY(U, lpivot, lpivot)) <= EPS_F)
                {
                    cl = -cl;
                }

                *up = M_P_ENTRY(U, lpivot, lpivot) - cl;
                M_P_ENTRY(U, lpivot, lpivot) = cl;
            }
            float b = *up * M_P_ENTRY(U, lpivot, lpivot);

            if ((ABS(b) > (EPS_F * 0.01f)) && (C_max >= 0))
            {
                b = 1.0f / b;

                for (unsigned short j = l1; j < (unsigned short)C_max; j++)
                {
                    float sm = M_P_ENTRY(U, lpivot, j) * *up;

                    for (unsigned short i = l1; i < m; i++)
                    {
                        sm += M_P_ENTRY(U, i, j) * M_P_ENTRY(U, i, lpivot);
                    }

                    if (ABS(sm) > (EPS_F * 0.01f))
                    {
                        sm = sm * b;

                        M_P_ENTRY(U, lpivot, j) += sm * *up;

                        for (unsigned short i = l1; i < m; i++)
                        {
                            M_P_ENTRY(U, i, j) += sm * M_P_ENTRY(U, i, lpivot);
                        }
                    }
                }
            }
        }
    }
    else
    {
        status |= 1U;
    }

    return (status);
}

static unsigned int H2(unsigned short lpivot, unsigned short l1, unsigned char transpose, unsigned short m, const matrix_t *const U, const float *const up, vector_t *const C)
{
    unsigned int status = 0;
    if ((lpivot < l1) && (l1 < m) && (ABS(M_P_ENTRY(U, lpivot, lpivot)) > TOL))
    {
        if (transpose == 1U)
        {
            float b = *up * M_P_ENTRY(U, lpivot, lpivot);

            if (ABS(b) > (EPS_F * 0.01f))
            {
                b = 1.0f / b;

                float sm = C->vector[lpivot] * *up;

                for (unsigned short i = l1; i < m; i++)
                {
                    sm += C->vector[i] * M_P_ENTRY(U, lpivot, i);
                }

                if (ABS(sm) > (EPS_F * 0.01f))
                {
                    sm = sm * b;

                    C->vector[lpivot] += sm * *up;

                    for (unsigned short i = l1; i < m; i++)
                    {
                        C->vector[i] += sm * M_P_ENTRY(U, lpivot, i);
                    }
                }
            }
        }
        else
        {
            float b = *up * M_P_ENTRY(U, lpivot, lpivot);

            if (ABS(b) > (EPS_F * 0.01f))
            {
                b = 1.0f / b;

                float sm = C->vector[lpivot] * *up;

                for (unsigned short i = l1; i < m; i++)
                {
                    sm += C->vector[i] * M_P_ENTRY(U, i, lpivot);
                }

                if (ABS(sm) > (EPS_F * 0.01f))
                {
                    sm = sm * b;

                    C->vector[lpivot] += sm * *up;

                    for (unsigned short i = l1; i < m; i++)
                    {
                        C->vector[i] += sm * M_P_ENTRY(U, i, lpivot);
                    }
                }
            }
        }
    }
    return (status);
}

unsigned int pivoted_QR_T(matrix_t *const B, vector_u_t *p, vector_t *h, const unsigned short m, const float tol, unsigned short *const rank, const unsigned short rank_max)
{
    float hmax = 0.0f;
    const float factor = 1e-5f;
    unsigned int status = 0U;

    unsigned short n = B->rows;
    // unsigned short m = B->columns;

    unsigned short mu = MIN(n, m);

    /* Reset Buffers */
    for (unsigned short i = 0U; i < h->rows; i++)
    {
        h->vector[i] = 0.0f;
    }

    for (unsigned short i = 0; i < n; i++)
    {
        p->vector[i] = 0U;
    }

    for (unsigned short j = 0U; j < mu; j++)
    {
        unsigned short lmax = j;
        /* update squared column lengths and find lmax */
        if (j != 0U)
        {
            for (unsigned short l = j; l < n; l++)
            {
                h->vector[l] -= M_P_ENTRY(B, l, (j - 1U)) * M_P_ENTRY(B, l, (j - 1U));
                if (h->vector[l] > h->vector[lmax])
                {
                    lmax = l;
                }
            }
            if (ABS((hmax + factor * h->vector[lmax]) - hmax) < EPS_F)
            {
                lmax = j;
                /* compute squared column lengths and find lmax */
                for (unsigned short l = j; l < n; l++)
                {
                    h->vector[l] = 0.0f;
                    for (unsigned short i = j; i < m; i++)
                    {
                        h->vector[l] += M_P_ENTRY(B, l, i) * M_P_ENTRY(B, l, i);
                    }
                    if (h->vector[l] > h->vector[lmax])
                    {
                        lmax = l;
                    }
                }
                hmax = h->vector[lmax];
            }
        }
        else
        {
            // lmax = j;
            /* compute squared column lengths and find lmax */
            for (unsigned short l = j; l < n; l++)
            {
                h->vector[l] = 0.0f;
                for (unsigned short i = j; i < m; i++)
                {
                    h->vector[l] += M_P_ENTRY(B, l, i) * M_P_ENTRY(B, l, i);
                }
                if (h->vector[l] > h->vector[lmax])
                {
                    lmax = l;
                }
            }
            hmax = h->vector[lmax];
        }

        /* lmax has been determined*/
        p->vector[j] = lmax;

        /* conduct column interchange if necessary */
        // TODO maybe get rid of actual swapping
        if (p->vector[j] != j)
        {
            for (unsigned short i = 0U; i < m; i++)
            {
                float tmp = M_P_ENTRY(B, j, i);
                M_P_ENTRY(B, j, i) = M_P_ENTRY(B, lmax, i);
                M_P_ENTRY(B, lmax, i) = tmp;
            }
            h->vector[lmax] = h->vector[j];
        }

        // print_matrix(B);
        /* Call H1 function */
        status |= H1(j, j + 1U, 1U, m, B, &h->vector[j], (short)n);
    }

    /* Determine the pseudorank, K, using the tolerance TOL */
    *rank = rank_max;
    if (ABS(M_P_ENTRY(B, 0U, 0U)) > tol)
    {
        float rel_TOL = MAX(tol, ABS(M_P_ENTRY(B, 0U, 0U)) * tol);

        for (unsigned short j = 1U; j < rank_max; j++)
        {
            if (ABS(M_P_ENTRY(B, j, j)) <= rel_TOL)
            {
                *rank = j;
                break;
            }
        }
    }
    else
    {
        *rank = 0U;
    }

    /* END of pivoted QR */
    return (status);
}

unsigned int pivoted_QR(matrix_t *const A, vector_u_t *p, vector_t *h, const unsigned short m, const unsigned short n, const float tol, unsigned short *const rank, const unsigned short rank_max)
{
    float hmax = 0.0f;
    const float factor = 1e-3f;
    unsigned int status = 0U;

    /* Reset Buffers */
    for (unsigned short i = 0U; i < h->rows; i++)
    {
        h->vector[i] = 0.0f;
    }

    for (unsigned short i = 0; i < n; i++)
    {
        p->vector[i] = 0U;
    }

    for (unsigned short j = 0U; j < n; j++)
    {
        unsigned short lmax = j;
        /* update squared column lengths and find lmax */
        if (j != 0U)
        {
            for (unsigned short l = j; l < n; l++)
            {
                h->vector[l] -= M_P_ENTRY(A, (j - 1U), l) * M_P_ENTRY(A, (j - 1U), l);
                if (h->vector[l] > h->vector[lmax])
                {
                    lmax = l;
                }
            }
            if ((ABS(hmax + (factor * h->vector[lmax])) - hmax) <= EPS_F)
            {
                lmax = j;
                /* compute squared column lengths and find lmax */
                for (unsigned short l = j; l < n; l++)
                {
                    h->vector[l] = 0.0f;
                    for (unsigned short i = j; i < m; i++)
                    {
                        h->vector[l] += M_P_ENTRY(A, i, l) * M_P_ENTRY(A, i, l);
                    }
                    if (h->vector[l] > h->vector[lmax])
                    {
                        lmax = l;
                    }
                }
                hmax = h->vector[lmax];
            }
        }
        else
        {
            lmax = j;
            /* compute squared column lengths and find lmax */
            for (unsigned short l = j; l < n; l++)
            {
                h->vector[l] = 0.0f;
                for (unsigned short i = j; i < m; i++)
                {
                    h->vector[l] += M_P_ENTRY(A, i, l) * M_P_ENTRY(A, i, l);
                }
                if (h->vector[l] > h->vector[lmax])
                {
                    lmax = l;
                }
            }
            hmax = h->vector[lmax];
        }

        /* lmax has been determined*/
        p->vector[j] = lmax;

        /* conduct column interchange if necessary */
        // TODO maybe get rid of actual swapping
        if (p->vector[j] != j)
        {
            for (unsigned short i = 0U; i < m; i++)
            {
                float tmp = M_P_ENTRY(A, i, j);
                M_P_ENTRY(A, i, j) = M_P_ENTRY(A, i, lmax);
                M_P_ENTRY(A, i, lmax) = tmp;
            }
            h->vector[lmax] = h->vector[j];
        }

        /* Call H1 function */
        status |= H1(j, j + 1U, 0U, m, A, &h->vector[j], (short)n);
    }

    /* Determine the pseudorank, K, using the tolerance TOL */
    *rank = rank_max;
    if (ABS(M_P_ENTRY(A, 0U, 0U)) > tol)
    {
        float rel_TOL = ABS(M_P_ENTRY(A, 0U, 0U)) * tol;
        if (rel_TOL < tol)
        {
            rel_TOL = tol;
        }

        for (unsigned short j = 1; j < rank_max; j++)
        {
            if (ABS(M_P_ENTRY(A, j, j)) <= rel_TOL)
            {
                *rank = j;
                break;
            }
        }
    }
    else
    {
        *rank = 0U;
    }

    /* END of pivoted QR */
    return (status);
}

unsigned int solve_QR_LS_udet(const matrix_t *const B, const vector_c_t *const nu_des, const unsigned short *const rank, vector_t *u_alloc, matrix_t *const R, const vector_u_t *p, const vector_t *h, vector_t *g)
{
    unsigned int status = 0U;

    unsigned short n = B->rows;
    unsigned short m = B->columns;

    /* Rank 0 solution */
    if (*rank == 0U)
    {
        for (unsigned short i = 0; i < m; i++)
        {
            u_alloc->vector[i] = 0.0f;
        }
    }
    else
    {
        /*Copy over nu des */
        for (unsigned short j = 0; j < n; j++)
        {
            u_alloc->vector[j] = nu_des->vector[j];
        }
        /* Required as simulink does not initialize with zeros */
        for (unsigned short j = n; j < m; j++)
        {
            u_alloc->vector[j] = 0.0f;
        }

        /* Sort according to pivoting vector */
        for (unsigned short j = 0; j < n; j++)
        {
            if (p->vector[j] != j)
            {
                float tmp = u_alloc->vector[p->vector[j]];
                u_alloc->vector[p->vector[j]] = u_alloc->vector[j];
                u_alloc->vector[j] = tmp;
            }
        }

        if (*rank != n)
        {
            /* Rank deficient Solution */
            /* Copy from B */
            for (unsigned short j = 0U; j < *rank; j++)
            {
                for (unsigned short i = j; i < n; i++)
                {
                    M_P_ENTRY(R, i, j) = M_P_ENTRY(B, i, j);
                }
            }

            /* Apply Householder QR again for R*/
            for (unsigned short i = 0U; i < *rank; i++)
            {
                status |= H1(i, i + 1U, 0U, n, R, &g->vector[i], (short)*rank);
                status |= H2(i, i + 1U, 0U, n, R, &g->vector[i], u_alloc);
            }

            /* Set entries of u_alloc to zero which cannot be allocated to achieve least distance solution*/
            for (unsigned short i = *rank; i < n; i++)
            {
                u_alloc->vector[i] = 0.0f;
            }

            /* Solve the reduced triangular system */
            for (unsigned short i = 0U; i < *rank; i++)
            {
                unsigned short ii = *rank - 1U - i;
                float sm = 0.0f;
                if (ii != (*rank - 1U))
                {
                    for (unsigned short j = ii + 1U; j < *rank; j++)
                    {
                        sm += M_P_ENTRY(R, ii, j) * u_alloc->vector[j];
                    }
                }
                u_alloc->vector[ii] = (u_alloc->vector[ii] - sm) / M_P_ENTRY(R, ii, ii);
            }
        }
        else
        {
            /* Full rank case */
            /* Solve the full triangular system */
            for (unsigned short i = 0U; i < *rank; i++)
            {
                float sm = 0.0f;
                if (i != 0U)
                {
                    for (unsigned short j = 0U; j < i; j++)
                    {
                        unsigned short jj = i - 1U - j;
                        sm += M_P_ENTRY(B, i, jj) * u_alloc->vector[jj];
                    }
                }
                u_alloc->vector[i] = (u_alloc->vector[i] - sm) / M_P_ENTRY(B, i, i);
            }
        }

        /* Multiply with Q^T*/
        for (unsigned short j = 0U; j < *rank; j++)
        {
            unsigned short jj = *rank - 1U - j;
            status |= H2(jj, jj + 1U, 1U, m, B, &h->vector[jj], u_alloc);
        }
    }
    return (status);
}

unsigned int solve_QR_LS_udet_SLS(const matrix_t *const B, const vector_t *const nu_des, const unsigned short m, const unsigned short *const rank, vector_t *u_alloc, matrix_t *const R, const vector_u_t *p, const vector_t *h, vector_t *g)
{
    unsigned int status = 0U;

    unsigned short n = B->rows;
    // unsigned short m = B->columns;

    unsigned short mu = MIN(n, m);

    /* Reset buffers */

    for (unsigned short i = 0U; i < g->rows; i++)
    {
        g->vector[i] = 0.0f;
    }

    for (unsigned short i = 0U; i < R->rows; i++)
    {
        for (unsigned short j = 0U; j < R->columns; j++)
        {
            M_P_ENTRY(R, i, j) = 0.0f;
        }
    }

    /* Rank 0 solution */
    if (*rank == 0U)
    {
        for (unsigned short i = 0; i < m; i++)
        {
            u_alloc->vector[i] = 0.0f;
        }
    }
    else
    {
        /*Copy over nu des */
        for (unsigned short j = 0; j < n; j++)
        {
            u_alloc->vector[j] = nu_des->vector[j];
        }
        /* Required as simulink does not initialize with zeros */
        for (unsigned short j = n; j < m; j++)
        {
            u_alloc->vector[j] = 0.0f;
        }

        /* Sort according to pivoting vector */
        for (unsigned short j = 0; j < mu; j++)
        {
            if (p->vector[j] != j)
            {
                float tmp = u_alloc->vector[p->vector[j]];
                u_alloc->vector[p->vector[j]] = u_alloc->vector[j];
                u_alloc->vector[j] = tmp;
            }
        }

        if (*rank != n)
        {
            /* Rank deficient Solution */
            /* Copy from B */
            for (unsigned short j = 0U; j < *rank; j++)
            {
                for (unsigned short i = j; i < n; i++)
                {
                    M_P_ENTRY(R, i, j) = M_P_ENTRY(B, i, j);
                }
            }

            /* Apply Householder QR again for R*/
            for (unsigned short i = 0U; i < *rank; i++)
            {
                status |= H1(i, i + 1U, 0U, n, R, &g->vector[i], (short)*rank);
                status |= H2(i, i + 1U, 0U, n, R, &g->vector[i], u_alloc);
            }

            /* Set entries of u_alloc to zero which cannot be allocated to achieve least distance solution*/
            for (unsigned short i = *rank; i < n; i++)
            {
                u_alloc->vector[i] = 0.0f;
            }

            /* Solve the reduced triangular system */
            for (unsigned short i = 0U; i < *rank; i++)
            {
                unsigned short ii = *rank - 1U - i;
                float sm = 0.0f;
                if (ii != (*rank - 1U))
                {
                    for (unsigned short j = ii + 1U; j < *rank; j++)
                    {
                        sm += M_P_ENTRY(R, ii, j) * u_alloc->vector[j];
                    }
                }
                u_alloc->vector[ii] = (u_alloc->vector[ii] - sm) / M_P_ENTRY(R, ii, ii);
            }
        }
        else
        {
            /* Full rank case */
            /* Solve the full triangular system */
            for (unsigned short i = 0U; i < *rank; i++)
            {
                float sm = 0.0f;
                if (i != 0U)
                {
                    for (unsigned short j = 0U; j < i; j++)
                    {
                        unsigned short jj = i - 1U - j;
                        sm += M_P_ENTRY(B, i, jj) * u_alloc->vector[jj];
                    }
                }
                u_alloc->vector[i] = (u_alloc->vector[i] - sm) / M_P_ENTRY(B, i, i);
            }
        }

        /* Multiply with Q */
        for (unsigned short j = 0U; j < *rank; j++)
        {
            unsigned short jj = *rank - 1U - j;
            status |= H2(jj, jj + 1U, 1U, m, B, &h->vector[jj], u_alloc);
        }
    }
    return (status);
}

unsigned int solve_QR_LS_odet_SLS(matrix_t *const A, const vector_t *const u_des, const unsigned short m, const unsigned short n, const unsigned short *const rank, vector_t *u_sol, const vector_u_t *p, const vector_t *h)
{
    unsigned int status = 0U;

    /* both vectors must be m long to use the same interface */
    /* can be changed if u_des is non constant */

    /* Rank 0 solution */
    if (*rank == 0U)
    {
        for (unsigned short i = 0; i < n; i++)
        {
            u_sol->vector[i] = 0.0f;
        }
    }
    else
    {
        /*Copy over nu des */
        for (unsigned short j = 0; j < m; j++)
        {
            u_sol->vector[j] = u_des->vector[j];
        }

        /* Multiply with Q^T */
        for (unsigned short j = 0U; j < n; j++)
        {
            status |= H2(j, j + 1U, 0U, m, A, &h->vector[j], u_sol);
        }

        /* Reset h vector */
        for (unsigned short j = 0U; j < n; j++)
        {
            h->vector[j] = 0.0f;
        }

        if (*rank < n)
        {
            /* Apply second QR decomposition */
            for (unsigned short i = 0U; i < *rank; i++)
            {
                unsigned short ii = *rank - 1U - i;
                status |= H1(ii, *rank, 1U, n, A, &h->vector[ii], ((short)ii) - 1);
            }

            /* Solve the full triangular system */
            for (unsigned short i = 0U; i < *rank; i++)
            {
                unsigned short ii = *rank - 1U - i;
                float sm = 0.0f;
                if (ii != (*rank - 1U))
                {
                    for (unsigned short j = ii + 1U; j < *rank; j++)
                    {
                        sm += M_P_ENTRY(A, ii, j) * u_sol->vector[j];
                    }
                }
                u_sol->vector[ii] = (u_sol->vector[ii] - sm) / M_P_ENTRY(A, ii, ii);
            }

            /* Multiply with T^T */
            for (unsigned short j = *rank; j < n; j++)
            {
                u_sol->vector[j] = 0.0f;
            }

            for (unsigned short j = 0U; j < *rank; j++)
            {
                status |= H2(j, *rank, 1U, n, A, &h->vector[j], u_sol);
            }
        }
        else
        {
            /* Full rank case */
            /* Solve the full triangular system */
            for (unsigned short i = 0U; i < *rank; i++)
            {
                unsigned short ii = *rank - 1U - i;
                float sm = 0.0f;
                if (ii != (*rank - 1U))
                {
                    for (unsigned short j = ii + 1U; j < *rank; j++)
                    {
                        sm += M_P_ENTRY(A, ii, j) * u_sol->vector[j];
                    }
                }
                u_sol->vector[ii] = (u_sol->vector[ii] - sm) / M_P_ENTRY(A, ii, ii);
            }
        }

        /* Sort according to pivoting vector */
        for (unsigned short j = 0; j < n; j++)
        {
            unsigned short jj = n - 1U - j;
            if (p->vector[jj] != jj)
            {
                float tmp = u_sol->vector[p->vector[jj]];
                u_sol->vector[p->vector[jj]] = u_sol->vector[jj];
                u_sol->vector[jj] = tmp;
            }
        }
    }
    return (status);
}

unsigned int nullspace(matrix_t *const B, vector_u_t *p, vector_t *h, const unsigned short m, const float tol, matrix_t *const Q2, unsigned short *const n_Q2)
{
    /* Declare and reset variables */
    unsigned int status = 0U;
    unsigned short rank = 0U;

    unsigned short rank_max = MIN(B->rows, B->columns);

    for (unsigned short i = 0U; i < p->rows; i++)
    {
        p->vector[i] = 0U;
    }

    for (unsigned short i = 0U; i < h->rows; i++)
    {
        h->vector[i] = 0.0f;
    }

    status |= pivoted_QR_T(B, p, h, m, tol, &rank, rank_max);

    /*Fill Q2*/
    if (m > rank)
    {
        *n_Q2 = m - rank;
        for (unsigned short i = 0U; i < Q2->rows; i++)
        {
            for (unsigned short j = 0U; j < *n_Q2; j++)
            {
                M_P_ENTRY(Q2, i, j) = 0.0f;
            }
        }

        for (unsigned short j = 0U; j < *n_Q2; j++)
        {
            M_P_ENTRY(Q2, rank + j, j) = 1.0f;
        }

        /* Multiply with Q^T*/
        for (unsigned short j = 0U; j < rank; j++)
        {
            unsigned short jj = rank - 1U - j;
            status |= H2_mat(jj, jj + 1U, 1U, m, B, &h->vector[jj], Q2, *n_Q2);
        }
    }
    else
    {
        *n_Q2 = 0U;
    }

    return (status);
}

unsigned int H2_mat(unsigned short lpivot, unsigned short l1, unsigned char transpose, unsigned short m, const matrix_t *const U, const float *const up, matrix_t *const C, unsigned short m_C)
{
    unsigned int status = 0;
    if ((lpivot < l1) && (l1 <= m) && (ABS(M_P_ENTRY(U, lpivot, lpivot)) > TOL))
    {
        if (transpose == 1U)
        {
            float b = *up * M_P_ENTRY(U, lpivot, lpivot);

            if (b < 0.0f)
            {
                b = 1.0f / b;

                for (unsigned short j = 0U; j < m_C; j++)
                {
                    float sm = M_P_ENTRY(C, lpivot, j) * *up;

                    for (unsigned short i = l1; i < m; i++)
                    {
                        sm += M_P_ENTRY(C, i, j) * M_P_ENTRY(U, lpivot, i);
                    }

                    if (ABS(sm) > EPS_F)
                    {
                        sm = sm * b;

                        M_P_ENTRY(C, lpivot, j) += sm * *up;

                        for (unsigned short i = l1; i < m; i++)
                        {
                            M_P_ENTRY(C, i, j) += sm * M_P_ENTRY(U, lpivot, i);
                        }
                    }
                }
            }
        }
        else
        {
            float b = *up * M_P_ENTRY(U, lpivot, lpivot);

            if (ABS(b) > EPS_F)
            {
                b = 1.0f / b;
                for (unsigned short j = 0U; j < m_C; j++)
                {
                    float sm = M_P_ENTRY(C, lpivot, j) * *up;

                    for (unsigned short i = l1; i < m; i++)
                    {
                        sm += M_P_ENTRY(C, i, j) * M_P_ENTRY(U, i, lpivot);
                    }

                    if (ABS(sm) > EPS_F)
                    {
                        sm = sm * b;

                        M_P_ENTRY(C, lpivot, j) += sm * *up;

                        for (unsigned short i = l1; i < m; i++)
                        {
                            M_P_ENTRY(C, i, j) += sm * M_P_ENTRY(U, i, lpivot);
                        }
                    }
                }
            }
        }
    }
    else
    {
        status |= 1U;
    }

    return (status);
}
