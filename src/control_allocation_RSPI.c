/*------------------------------------------------------------------------*/
/*                      control_allocation_RSPI.c                         */
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
/*Date:         2024-11-19                                                */
/*Contact:      Institute of Flight System Dynamics                       */
/*              School of Engineering and Design                          */
/*              Technische Universitaet Muenchen                          */
/*              Boltzmannstrasse 15                                       */
/*              D-85748 Garching                                          */
/*              simon.hafner@tum.de                                       */
/*------------------------------------------------------------------------*/

#include "control_allocation_RSPI.h"

/**
 * @brief Executes the SLS AS control allocation
 *
 */
unsigned int CA_RSPI(float *const B, float *const W_u, float *const W_nu, float *const nu_des, float *const u_min, float *const u_max,
                     float *const u_des, short *const W, float *const u_alloc, const unsigned short n, const unsigned short m, const unsigned short n_W_u,
                     const unsigned short max_iter, const float abort_tol, float *const buffer_float, unsigned short *const buffer_unsigned)
{
    unsigned int status = 0;

    /* Initialzation */
    unsigned short phase = 1U;
    unsigned short m_free = m;

    unsigned short rank_B = 0U;

    float residual = 1.0f;

    matrix_t B_i = {.rows = n, .columns = m, .matrix = B};
    vector_t W_nu_i = {.rows = n, .vector = W_nu};
    matrix_t W_u_i = {.rows = n_W_u, .columns = m, .matrix = W_u};
    vector_t nu_des_i = {.rows = n, .vector = nu_des};
    vector_t u_min_i = {.rows = m, .vector = u_min};
    vector_t u_max_i = {.rows = m, .vector = u_max};
    vector_t u_des_i = {.rows = m, .vector = u_des};
    vector_i_t W_i = {.rows = m, .vector = W};
    vector_t u_alloc_i = {.rows = m, .vector = u_alloc};

    /* Buffers */
    unsigned short counter = 0U;
    matrix_t B_free = {.rows = n, .columns = m, .matrix = &buffer_float[counter]};
    counter += n * m;

    matrix_t Q2 = {.rows = m, .columns = m, .matrix = &buffer_float[counter]};
    counter += m * m;
    matrix_t W_u_free = {.rows = m, .columns = m, .matrix = &buffer_float[counter]};
    counter += m * m;
    matrix_t C0 = {.rows = m, .columns = m, .matrix = &buffer_float[counter]};
    counter += m * m;

    matrix_t R = {.rows = n, .columns = n, .matrix = &buffer_float[counter]};
    counter += n * n;

    vector_t h_12 = {.rows = m, .vector = &buffer_float[counter]};
    counter += m;
    vector_t g_12 = {.rows = m, .vector = &buffer_float[counter]};
    counter += m;
    vector_t u_delta = {.rows = m, .vector = &buffer_float[counter]};
    counter += m;
    vector_t u_opt = {.rows = m, .vector = &buffer_float[counter]};
    counter += m;

    counter = 0U;
    vector_u_t permutation = {.rows = m, .vector = &buffer_unsigned[counter]};
    counter += m;
    vector_u_t pivot = {.rows = m, .vector = &buffer_unsigned[counter]};

    // % ||Wv(B(u0+p)-v)|| = ||A u_delta-d||
    /* Scale B with weighting matrix */
    /* nu weighting limited to diagonal with this implementation */
    for (unsigned short i = 0U; i < n; i++)
    {
        for (unsigned short j = 0U; j < m; j++)
        {
            M_ENTRY(B_i, i, j) = W_nu_i.vector[i] * M_ENTRY(B_i, i, j);
        }
    }

    /* W_nu*nu_des */
    for (unsigned short i = 0U; i < n; i++)
    {
        nu_des_i.vector[i] = W_nu_i.vector[i] * nu_des_i.vector[i];
    }

    /* W_nu*nu_des - B*u_alloc */
    for (unsigned short i = 0U; i < n; i++)
    {
        for (unsigned short j = 0U; j < m; j++)
        {
            nu_des_i.vector[i] -= M_ENTRY(B_i, i, j) * u_alloc_i.vector[j];
        }
    }

    unsigned short iter = 0U;
    /* Iterate until the optimum is found or max iterations have been reached*/
    for (iter = 0U; iter < max_iter; iter++)
    {

        /* Compute the permutation vector */
        /* [indexes of free u; indeces of saturated u]*/
        /* get number of free u (m_free)*/
        status |= get_permutations(&W_i, &permutation, &m_free);

        /* Copy and rearrange control effectiveness matrix based on free variables */
        for (unsigned short j = 0; j < m_free; j++)
        {
            /* Copy corresponding colums*/
            for (unsigned short i = 0; i < n; i++)
            {
                M_ENTRY(B_free, i, j) = M_ENTRY(B_i, i, permutation.vector[j]);
            }
        }

        if (phase == 1U)
        {
            /* ---------------------- */
            /*        Phase I         */
            /* ---------------------- */

            /* Solve the underdetermined equation system B_free u_free = nu_des */
            /* Calculate pivoted QR of B*/
            status |= pivoted_QR_T(&B_free, &pivot, &h_12, m_free, TOL, &rank_B, MIN(n, m_free));

            /* Solve Least Squares Problem */
            status |= solve_QR_LS_udet_SLS(&B_free, &nu_des_i, m_free, &rank_B, &u_delta, &R, &pivot, &h_12, &g_12);

            /* -------------------------- */
            /* Check if u_opt is feasible */
            /* -------------------------- */
            // TODO Could also do this based on the scaling alpha

            /* Feasiblility Check */
            unsigned short feasible = check_feasibility(&u_opt, &u_delta, &u_alloc_i, &u_min_i, &u_max_i, &permutation, m_free);

            if (feasible == 1U)
            {
                /* -------------------------- */
                /* Optimal point is feasible  */
                /* -------------------------- */

                /* Update the output u_alloc */
                /* u_alloc = u_opt; */
                for (unsigned short i = 0U; i < m_free; i++)
                {
                    u_alloc_i.vector[permutation.vector[i]] = u_opt.vector[permutation.vector[i]];
                }

                /* Update residual pseudo control */
                /* nu_des = nu_des - B_free * u_delta; */
                for (unsigned short i = 0U; i < n; i++)
                {
                    for (unsigned short j = 0U; j < m_free; j++)
                    {
                        nu_des_i.vector[i] = nu_des_i.vector[i] - (M_ENTRY(B_free, i, j) * u_delta.vector[j]);
                    }
                }

                /* switch to phase II */
                phase = 2U;

                /* Update u_des based on */
                /* u_opt = (u_des - u_alloc) */
                /* using u_opt as buffer */
                for (unsigned short i = 0U; i < m; i++)
                {
                    u_opt.vector[i] = u_des_i.vector[i] - u_alloc_i.vector[i];
                }

                /* u_des = W_u * u_opt */
                for (unsigned short i = 0U; i < n_W_u; i++)
                {
                    u_des_i.vector[i] = 0.0f;
                    for (unsigned short j = 0U; j < m; j++)
                    {
                        u_des_i.vector[i] = u_des_i.vector[i] + (M_ENTRY(W_u_i, i, j) * u_opt.vector[j]);
                    }
                }

                /* Reset working set */
                for (unsigned short i = 0U; i < m; i++)
                {
                    W_i.vector[i] = 0;
                }
            }
            else
            {
                /* ------------------------------ */
                /* u_delta leads to infeasible u  */
                /* ------------------------------ */

                /* Scale down u_delta by alpha such that it leads to a feasbile u*/

                /* Compute scaling factor alpha */
                /* 0 <= alpha <= 1 */
                /* for feasible u, alpha == 1*/

                float alpha = 1.0f;
                unsigned short idx_alpha = 0U;
                get_scaling(&u_delta, &u_alloc_i, &u_min_i, &u_max_i, &permutation, &W_i, m_free, &alpha, &idx_alpha);

                /* Update u_alloc and the residual nu_des*/
                /* u_opt = alpha * u_delta */
                /* u_opt is used as a buffer */
                for (unsigned short i = 0; i < m_free; i++)
                {
                    u_opt.vector[permutation.vector[i]] = alpha * u_delta.vector[i];
                }

                for (unsigned short i = m_free; i < m; i++)
                {
                    u_opt.vector[permutation.vector[i]] = 0.0f;
                }

                /* u_alloc = u_alloc + u_opt; */
                for (unsigned short i = 0; i < m; i++)
                {
                    u_alloc_i.vector[i] += u_opt.vector[i];
                }

                /* nu_des = nu_des - B * u_delta; */
                for (unsigned short i = 0U; i < n; i++)
                {
                    for (unsigned short j = 0U; j < m; j++)
                    {
                        nu_des_i.vector[i] -= (M_ENTRY(B_i, i, j) * u_opt.vector[j]);
                    }
                }

                /* Add the saturated input to the working set */
                if (u_delta.vector[idx_alpha] >= 0.0f)
                {
                    W_i.vector[permutation.vector[idx_alpha]] = 1;
                }
                else
                {
                    W_i.vector[permutation.vector[idx_alpha]] = -1;
                }

                /* Check if residual nu is so small that it can be stopped */
                /* Update residual */
                residual *= (1.0f - alpha);
                if (residual < abort_tol)
                {
                    /* Switch to phase II */
                    phase = 2;

                    /* Update u_des based on */
                    /* u_opt = (u_des - u_alloc) */
                    /* using u_opt as buffer */
                    for (unsigned short i = 0U; i < m; i++)
                    {
                        u_opt.vector[i] = u_des_i.vector[i] - u_alloc_i.vector[i];
                    }

                    /* u_des = W_u_i * u_opt */
                    for (unsigned short i = 0U; i < n_W_u; i++)
                    {
                        u_des_i.vector[i] = 0.0f;
                        for (unsigned short j = 0U; j < m; j++)
                        {
                            u_des_i.vector[i] = u_des_i.vector[i] + (M_ENTRY(W_u_i, i, j) * u_opt.vector[j]);
                        }
                    }
                }
            }
        }
        else
        {
            /* ---------------------- */
            /*        Phase II        */
            /* ---------------------- */

            if (m_free <= 1U)
            {
                // No Nullspace transition possible
                status = iter;
                return (status);
            }

            /* Copy and rearrange constraint matrix based on free variables */
            for (unsigned short j = 0; j < m_free; j++)
            {
                /* Copy corresponding colums*/
                for (unsigned short i = 0; i < m; i++)
                {
                    M_ENTRY(W_u_free, i, j) = M_ENTRY(W_u_i, i, permutation.vector[j]);
                }
            }

            /* Compute nullspace Q2 of B_free */
            unsigned short n_Q2 = 0U;
            /* Compute nullspace*/
            status |= nullspace(&B_free, &pivot, &h_12, m_free, TOL, &Q2, &n_Q2);

            /* If the nullspace is empty no nullspace transition is possible */
            /* Therefore abort phase II */
            if (n_Q2 == 0U)
            {
                status = iter;
                return (status);
            }

            /* Project the desired u_des onto the nullspace of B using the matrix W_u*/
            /* u_delta = Q2 * (W_u * Q2)^+ u_des */

            /* Replaced W_u_Q2 with C0 */
            /* W_u_Q2 = W_u_free * Q2 */
            for (unsigned short k = 0U; k < n_W_u; k++) /* left row counter               */
            {
                for (unsigned short j = 0U; j < n_Q2; j++) /*  right column counter          */
                {
                    M_ENTRY(C0, k, j) = 0.0f;
                    for (unsigned short i = 0U; i < m_free; i++) /* left col counter               */
                    {
                        M_ENTRY(C0, k, j) = M_ENTRY(C0, k, j) + /* calculate result     */
                                            (M_ENTRY(W_u_free, k, i) * M_ENTRY(Q2, i, j));
                    }
                }
            }

            /* Compute g_12 = (W_u * Q2)^+ u_des */
            status |= pivoted_QR(&C0, &pivot, &h_12, m_free, n_Q2, TOL, &rank_B, n_Q2); // Check max rank feature
            status |= solve_QR_LS_odet_SLS(&C0, &u_des_i, m_free, n_Q2, &rank_B, &g_12, &pivot, &h_12);

            /* u_delta = Q2 * g12 */
            for (unsigned short i = 0U; i < m_free; i++)
            {
                u_delta.vector[i] = 0.0f;
                for (unsigned short j = 0U; j < n_Q2; j++)
                {
                    u_delta.vector[i] += (M_ENTRY(Q2, i, j) * g_12.vector[j]);
                }
            }

            /* -------------------------- */
            /* Check if u_opt is feasible */
            /* -------------------------- */

            /* Feasiblility Check */
            unsigned short feasible = check_feasibility(&u_opt, &u_delta, &u_alloc_i, &u_min_i, &u_max_i, &permutation, m_free);

            if (feasible == 1U)
            {
                /* -------------------------- */
                /* Optimal point is feasible  */
                /* -------------------------- */

                /* Update the output u_alloc */
                /* u_alloc = u_opt; */
                for (unsigned short i = 0U; i < m_free; i++)
                {
                    u_alloc_i.vector[permutation.vector[i]] = u_opt.vector[permutation.vector[i]];
                }

                /* Update residual desired input */
                /* u_des = u_des - W_u * u_delta; */ // Could also only use free p here
                for (unsigned short i = 0U; i < n_W_u; i++)
                {
                    for (unsigned short j = 0U; j < m_free; j++)
                    {
                        u_des_i.vector[i] = u_des_i.vector[i] - (M_ENTRY(W_u_i, i, permutation.vector[j]) * u_delta.vector[j]);
                    }
                }

                /* --------------------------- */
                /* Stop RSPI NST               */
                /* --------------------------- */

                status = iter;
                return status;
            }
            else
            {
                /* ------------------------------ */
                /* u_delta leads to infeasible u  */
                /* ------------------------------ */

                /* Scale down u_delta by alpha such that it leads to a feasbile u*/

                /* Compute scaling factor alpha */
                /* 0 <= alpha <= 1 */
                /* for feasible u, alpha == 1*/

                float alpha = 1.0f;
                unsigned short idx_alpha = 0U;
                get_scaling(&u_delta, &u_alloc_i, &u_min_i, &u_max_i, &permutation, &W_i, m_free, &alpha, &idx_alpha);

                /* Update u_alloc and the residual u_des*/
                /* u_opt = alpha * u_delta */
                /* u_opt is used as a buffer */
                for (unsigned short i = 0; i < m_free; i++)
                {
                    u_opt.vector[permutation.vector[i]] = alpha * u_delta.vector[i];
                }

                for (unsigned short i = m_free; i < m; i++)
                {
                    u_opt.vector[permutation.vector[i]] = 0.0f;
                }

                /* u_alloc = u_alloc + u_opt; */
                for (unsigned short i = 0; i < m; i++)
                {
                    u_alloc_i.vector[i] += u_opt.vector[i];
                }

                /* u_des = u_des - W_u * u_delta; */
                for (unsigned short i = 0U; i < n_W_u; i++)
                {
                    for (unsigned short j = 0U; j < m; j++)
                    {
                        u_des_i.vector[i] -= (M_ENTRY(W_u_i, i, j) * u_opt.vector[j]);
                    }
                }

                /* Add the saturated input to the working set */
                if (u_delta.vector[idx_alpha] >= 0.0f)
                {
                    W_i.vector[permutation.vector[idx_alpha]] = 1;
                }
                else
                {
                    W_i.vector[permutation.vector[idx_alpha]] = -1;
                }
            }
        }
    }
    status = iter;
    return (status);
}
