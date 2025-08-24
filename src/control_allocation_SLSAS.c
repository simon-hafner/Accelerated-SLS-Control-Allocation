/*------------------------------------------------------------------------*/
/*                      control_allocation_SLSAS.c                        */
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

#include "control_allocation_SLSAS.h"

/* Static helper functions */

static unsigned short lambda_positive(const vector_t *const lambda, unsigned short m, const float tol)
{
    unsigned short optimum = 1U;
    for (unsigned short i = 0U; i < m; i++)
    {
        if (lambda->vector[i] < -tol)
        {
            optimum = 0U;
        }
    }
    return optimum;
}

static unsigned short lambda_min(const vector_t *const lambda, unsigned short m)
{
    unsigned short idx_min_lambda = 0U;
    float min_lambda = TOL;

    for (unsigned short i = 0U; i < m; i++)
    {
        if (lambda->vector[i] < min_lambda)
        {
            min_lambda = lambda->vector[i];
            idx_min_lambda = i;
        }
    }

    return idx_min_lambda;
}

/**
 * @brief Executes the SLS AS control allocation
 *
 */
unsigned int CA_SLSAS(float *const B, float *const W_u, float *const W_nu, float *const nu_des, float *const u_min, float *const u_max,
                      float *const u_des, short *const W, float *const u_alloc, unsigned short *const rank_B, const unsigned short n, const unsigned short m, const unsigned short n_W_u,
                      const unsigned short max_iter, const float abort_tol, const unsigned short warmstart, float *const buffer_float, unsigned short *const buffer_unsigned)
{
    unsigned int status = 0;

    /* Initialzation */
    unsigned short phase = 1U;
    unsigned short m_free = m;

    unsigned short B_QR_ready = 0U;

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
    matrix_t B_QR = {.rows = n, .columns = m, .matrix = &buffer_float[counter]};
    counter += n * m;

    matrix_t Q2 = {.rows = m, .columns = m, .matrix = &buffer_float[counter]};
    counter += m * m;
    matrix_t W_u_free = {.rows = m, .columns = m, .matrix = &buffer_float[counter]};
    counter += m * m;
    matrix_t C0 = {.rows = m, .columns = m, .matrix = &buffer_float[counter]};
    counter += m * m;

    matrix_t R = {.rows = n, .columns = n, .matrix = &buffer_float[counter]};
    counter += n * n;

    vector_t h_QR = {.rows = n, .vector = &buffer_float[counter]};
    counter += n;
    vector_t h_12 = {.rows = m, .vector = &buffer_float[counter]};
    counter += m;
    vector_t g_12 = {.rows = m, .vector = &buffer_float[counter]};
    counter += m;
    vector_t u_delta = {.rows = m, .vector = &buffer_float[counter]};
    counter += m;
    vector_t lambda = {.rows = m, .vector = &buffer_float[counter]};
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

    if ((warmstart == 1U) && (*rank_B == n))
    {

        /* Check if initialization point also fits with new limits */
        unsigned short feasible = 1U;
        for (unsigned short i = 0U; i < m; i++)
        {
            if ((u_alloc_i.vector[i] < (u_min_i.vector[i] - (TOL * 10.0f))) ||
                (u_alloc_i.vector[i] > (u_max_i.vector[i] + (TOL * 10.0f))))
            {
                feasible = 0U;
            }
        }

        if (feasible == 0U)
        {
            for (unsigned short i = 0; i < m; i++)
            {
                /* Initialize in the center of the u space */
                u_alloc_i.vector[i] = (u_max_i.vector[i] + u_min_i.vector[i]) * 0.5f;
            }

            for (unsigned short i = 0; i < m; i++)
            {
                /* Initialize W */
                W_i.vector[i] = 0;
            }
        }
    }
    else
    {
        for (unsigned short i = 0; i < m; i++)
        {
            /* Initialize in the center of the u space */
            u_alloc_i.vector[i] = (u_max_i.vector[i] + u_min_i.vector[i]) * 0.5f;
        }

        for (unsigned short i = 0; i < m; i++)
        {
            /* Initialize W */
            W_i.vector[i] = 0;
        }
    }

    /* Reset Rank of B */
    *rank_B = 0U;

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

    unsigned short iter = 0;
    /* Iterate until the optimum is found or max iterations have been reached*/
    for (iter = 0; iter < max_iter; iter++)
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
            status |= pivoted_QR_T(&B_free, &pivot, &h_12, m_free, TOL, rank_B, MIN(n, m_free));

            /* Solve Least Squares Problem */
            status |= solve_QR_LS_udet_SLS(&B_free, &nu_des_i, m_free, rank_B, &u_delta, &R, &pivot, &h_12, &g_12);

            /* Check if matrix for the QR was original matrix */
            if (m_free == m)
            {
                /* Store QR decomposition of inital matrix data for phase II optimality check */
                for (unsigned short i_rows = 0U; i_rows < B_free.rows; i_rows++)
                {
                    for (unsigned short i_cols = 0U; i_cols < B_free.columns; i_cols++)
                    {
                        M_ENTRY(B_QR, i_rows, i_cols) = M_ENTRY(B_free, i_rows, i_cols); /* Print each Column element in one row*/
                    }
                }

                for (unsigned short i = 0; i < n; i++)
                {
                    h_QR.vector[i] = h_12.vector[i];
                }
                B_QR_ready = 1U;
            }

            /* -------------------------- */
            /* Check if u_opt is feasible */
            /* -------------------------- */
            // TODO Could also do this based on the scaling alpha

            /* Compute scaling factor alpha */
            /* 0 <= alpha <= 1 */
            /* for feasible u, alpha == 1*/

            float alpha = 1.0f;
            unsigned short idx_alpha = 0U;
            get_scaling(&u_delta, &u_alloc_i, &u_min_i, &u_max_i, &permutation, &W_i, m_free, &alpha, &idx_alpha);

            /* Feasiblility Check */
            unsigned short feasible = 0U;

            if (alpha > (1.0f - (10.0f * TOL)))
            {
                feasible = 1U;
            }

            if (feasible == 1U)
            {
                /* -------------------------- */
                /* Optimal point is feasible  */
                /* -------------------------- */

                /* Update u_alloc and the residual nu_des*/
                /* u_delta = alpha * u_delta */
                for (unsigned short i = 0; i < m_free; i++)
                {
                    u_delta.vector[i] = alpha * u_delta.vector[i];
                }

                /* u_alloc = u_alloc + u_delta; */
                for (unsigned short i = 0; i < m_free; i++)
                {
                    u_alloc_i.vector[permutation.vector[i]] += u_delta.vector[i];
                }

                /* nu_des = nu_des - B * u_delta; */
                for (unsigned short i = 0U; i < n; i++)
                {
                    for (unsigned short j = 0U; j < m_free; j++)
                    {
                        nu_des_i.vector[i] -= (M_ENTRY(B_i, i, permutation.vector[j]) * u_delta.vector[j]);
                    }
                }

                /* ---------------------------- */
                /* Check if the optimal u_alloc */
                /* is also optimal in the sense */
                /* of the secondary goal        */
                /* ---------------------------- */

                if (*rank_B == n)
                {
                    /* if rank_B == n, we know that if the u_alloc is feasible */
                    /* that B u_alloc = nu_des and the optimum of phase 1 has been reached */
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
                }
                else
                {
                    /* If B_free was rank deficient (full nu_des might not be reached), */
                    /* optimality needs to be checked through Lagrange multipliers. */
                    /* Compute gradient of optimization criterion. */
                    /* lambda = W. B' nu_des */

                    /* g = B'*nu_des;   */
                    // status |= mat_t_vec(B, nu_des, &g_12);

                    for (unsigned short j = 0U; j < m; j++)
                    {
                        g_12.vector[j] = 0.0f;
                        for (unsigned short i = 0U; i < n; i++)
                        {
                            g_12.vector[j] = g_12.vector[j] + (M_ENTRY(B_i, i, j) * nu_des_i.vector[i]);
                        }
                    }

                    /* lambda = W.*g */
                    for (unsigned short i = 0U; i < m; i++)
                    {
                        lambda.vector[i] = (float)W_i.vector[i] * g_12.vector[i]; // TODO Could be optimized as W only contains -1, 0, 1
                    }

                    /* Check if each element of lambda is positive or zero */
                    /* by comparing to -TOL */
                    /* If this is true, the optimum has been reached */

                    unsigned short optimum = lambda_positive(&lambda, m, abort_tol);

                    if (optimum == 1U)
                    {
                        /* Optimum for phase I is found */
                        if (*rank_B < m_free)
                        {
                            /* If the system has coplanar effectors, further optimization is possible */
                            /* even though the B_free is rank_deficient. This is true if rank_B < m_free. */
                            /* We know that ||Wv(Bu-v)|| is minimal and can further optimize ||Wu(u-ud)|| */
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
                        else
                        {
                            /* No further optimization of the secondary objective */
                            /* is possible. Therefore, stop the algorithm. */
                            status = iter;
                            return (status);
                        }
                    }
                    else
                    {
                        /* --------------------------- */
                        /* Optimal point is not found  */
                        /* --------------------------- */

                        /* The saturated input resulting in the most negative lambda */
                        /* is freed (removed from the working set) */
                        unsigned short idx_min_lambda = lambda_min(&lambda, m);

                        W_i.vector[idx_min_lambda] = 0;
                    }
                }
            }
            else
            {
                /* ------------------------------ */
                /* u_delta leads to infeasible u  */
                /* ------------------------------ */

                /* Add the saturated input to the working set */
                if (u_delta.vector[idx_alpha] >= 0.0f)
                {
                    W_i.vector[permutation.vector[idx_alpha]] = 1;
                }
                else
                {
                    W_i.vector[permutation.vector[idx_alpha]] = -1;
                }

                /* Scale down u_delta by alpha such that it leads to a feasbile u*/

                /* Update u_alloc and the residual nu_des*/
                /* u_delta = alpha * u_delta */
                for (unsigned short i = 0; i < m_free; i++)
                {
                    u_delta.vector[i] = alpha * u_delta.vector[i];
                }

                /* u_alloc = u_alloc + u_delta; */
                for (unsigned short i = 0; i < m_free; i++)
                {
                    u_alloc_i.vector[permutation.vector[i]] += u_delta.vector[i];
                }

                /* nu_des = nu_des - B * u_delta; */
                for (unsigned short i = 0U; i < n; i++)
                {
                    for (unsigned short j = 0U; j < m_free; j++)
                    {
                        nu_des_i.vector[i] -= (M_ENTRY(B_i, i, permutation.vector[j]) * u_delta.vector[j]);
                    }
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

            if (m_free == 0U)
            {
                // No Nullspace transition possible
                status = iter;
                return (status);
            }

            /* Copy and rearrange constraint matrix based on free variables */
            for (unsigned short j = 0; j < m_free; j++)
            {
                /* Copy corresponding colums*/
                for (unsigned short i = 0; i < n_W_u; i++)
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

            unsigned short rank_W_Q2 = 0U;

            /* Compute g_12 = (W_u * Q2)^+ u_des */
            status |= pivoted_QR(&C0, &pivot, &h_12, n_W_u, n_Q2, TOL, &rank_W_Q2, MIN(n_Q2, n_W_u)); // Check max rank feature
            status |= solve_QR_LS_odet_SLS(&C0, &u_des_i, n_W_u, n_Q2, &rank_W_Q2, &g_12, &pivot, &h_12);

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

            /* Compute scaling factor alpha */
            /* 0 <= alpha <= 1 */
            /* for feasible u, alpha == 1*/

            float alpha = 1.0f;
            unsigned short idx_alpha = 0U;
            get_scaling(&u_delta, &u_alloc_i, &u_min_i, &u_max_i, &permutation, &W_i, m_free, &alpha, &idx_alpha);

            /* Feasiblility Check */

            if (alpha > (1.0f - (10.0f * TOL)))
            {
                /* -------------------------- */
                /* Optimal point is feasible  */
                /* -------------------------- */

                /* Update the output u_alloc */
                /* u_alloc = u_opt; */
                for (unsigned short i = 0U; i < m_free; i++)
                {
                    u_alloc_i.vector[permutation.vector[i]] += u_delta.vector[i];
                }

                /* Update residual desired input */
                /* u_des = u_des - W_u * u_delta; */ // Could also only use free p here
                for (unsigned short i = 0U; i < n_W_u; i++)
                {
                    for (unsigned short j = 0U; j < m_free; j++)
                    {
                        u_des_i.vector[i] = u_des_i.vector[i] - (M_ENTRY(W_u_free, i, j) * u_delta.vector[j]);
                    }
                }

                /* ---------------------------- */
                /* Check if the u_alloc is the  */
                /* optimum                      */
                /* ---------------------------- */
                if (B_QR_ready == 0U)
                {
                    /* QR decomposition of full matrix was not yet computed */

                    /* Copy full B to B_QR*/
                    for (unsigned short j = 0; j < m; j++)
                    {
                        /* Copy corresponding colums*/
                        for (unsigned short i = 0; i < n; i++)
                        {
                            M_ENTRY(B_QR, i, j) = M_ENTRY(B_i, i, j);
                        }
                    }

                    unsigned short rank_B_full = 0U;
                    status |= pivoted_QR_T(&B_QR, &pivot, &h_QR, m, TOL, &rank_B_full, n);
                    B_QR_ready = 1U;
                }

                /* g_12 = -W_u'*u_des; */
                for (unsigned short j = 0U; j < W_u_i.columns; j++)
                {
                    g_12.vector[j] = 0.0f;
                    for (unsigned short i = 0U; i < n_W_u; i++)
                    {
                        g_12.vector[j] = g_12.vector[j] - /*Calculate matrix vector product*/
                                         (M_ENTRY(W_u_i, i, j) * u_des_i.vector[i]);
                    }
                }

                /* Compute the lagrange multipliers related to the saturared inputs (active equality constraints)*/
                /* An effecive way to do so is by using the iterative QR decomposition */
                /* We can solve the least squares problem through the QR of (C_0 - Q_11*Q_11'*C_0)*/

                /* Constuct C0 which represents the active inequalities as equalities*/
                unsigned short m_sat = m - m_free;

                for (unsigned short i = 0U; i < C0.rows; i++)
                {
                    for (unsigned short j = 0U; j < C0.columns; j++)
                    {
                        M_ENTRY(C0, i, j) = 0.0f;
                    }
                }
                for (unsigned short i = 0U; i < m_sat; i++)
                {
                    M_ENTRY(C0, permutation.vector[i + m_free], i) = -(float)W_i.vector[permutation.vector[i + m_free]];
                }

                /* First get (C_0 - Q_11*Q_11'*C_0) */
                /* Q_11'*C_0 */
                // TODO Could be further optimized to compute only the changed one
                for (unsigned short j = 0U; j < n; j++)
                {
                    status |= H2_mat(j, j + 1U, 1U, m, &B_QR, &h_QR.vector[j], &C0, m_sat);
                }

                /* C0(n + 1 : end, :) = 0; */
                /* Necessary for the algorithm */ // Check if this could be avoided
                for (unsigned short i = n; i < m; i++)
                {
                    for (unsigned short j = 0U; j < m_sat; j++)
                    {
                        M_ENTRY(C0, i, j) = 0.0f;
                    }
                }

                /* Q_11*Q_11'*C_0 */
                for (unsigned short j = 0U; j < n; j++)
                {
                    unsigned short jj = n - 1U - j;
                    status |= H2_mat(jj, jj + 1U, 1U, m, &B_QR, &h_QR.vector[jj], &C0, m_sat);
                }

                /* C0 - Q_11*Q_11'*C_0 */
                for (unsigned short j = 0U; j < m_sat; j++)
                {
                    for (unsigned short i = 0U; i < m; i++)
                    {
                        /* Dynamically construct C0 in this operation*/
                        if (i == permutation.vector[j + m_free])
                        {
                            M_ENTRY(C0, i, j) = -(float)W_i.vector[i] - M_ENTRY(C0, i, j);
                        }
                        else
                        {
                            M_ENTRY(C0, i, j) = -M_ENTRY(C0, i, j);
                        }
                    }
                }

                /* Solve overdetermined and full-rank least-squares problem */
                /* (C_0 - Q_11*Q_11'*C_0)*lambda = g_12 */
                unsigned short rank_C0 = 0U;                                              // TODO Issue lies in W computation
                status |= pivoted_QR(&C0, &pivot, &h_12, m, m_sat, TOL, &rank_C0, m_sat); // TODO Check pivot size
                status |= solve_QR_LS_odet_SLS(&C0, &g_12, m, m_sat, &rank_C0, &lambda, &pivot, &h_12);

                /* Check if each element of lambda is positive or zero */
                /* by comparing to -TOL */
                /* If this is true, the optimum has been reached */
                unsigned short optimum = lambda_positive(&lambda, m_sat, abort_tol);

                if (optimum == 1U)
                {

                    /* --------------------------- */
                    /* Overall Optimum reached     */
                    /* --------------------------- */
                    status = iter;
                    return status;
                }
                else
                {
                    /* --------------------------- */
                    /* Optimal point is not found  */
                    /* --------------------------- */

                    /* The saturated input resulting in the most negative lambda */
                    /* is freed (removed from the working set) */

                    unsigned short idx_min_lambda = lambda_min(&lambda, m_sat);

                    W_i.vector[permutation.vector[idx_min_lambda + m_free]] = 0;
                }
            }
            else
            {
                /* ------------------------------ */
                /* u_delta leads to infeasible u  */
                /* ------------------------------ */

                /* Scale down u_delta by alpha such that it leads to a feasbile u*/

                /* Update u_alloc and the residual u_des*/
                /* u_opt = alpha * u_delta */
                /* u_opt is used as a buffer */

                if (alpha < TOL)
                {
                    alpha = 0.0f;
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

                /* u_delta = alpha * u_delta */
                for (unsigned short i = 0; i < m_free; i++)
                {
                    u_delta.vector[i] = alpha * u_delta.vector[i];
                }

                /* u_alloc = u_alloc + u_delta; */
                for (unsigned short i = 0; i < m_free; i++)
                {
                    u_alloc_i.vector[permutation.vector[i]] += u_delta.vector[i];
                }

                /* u_des = u_des - W_u * u_delta; */
                for (unsigned short i = 0U; i < n_W_u; i++)
                {
                    for (unsigned short j = 0U; j < m_free; j++)
                    {
                        u_des_i.vector[i] -= (M_ENTRY(W_u_free, i, j) * u_delta.vector[j]);
                    }
                }
            }
        }
    }
    status = iter;
    return (status);
}

unsigned int get_permutations(const vector_i_t *const W, vector_u_t *const permutation, unsigned short *m_free)
{
    unsigned int status = 0U;

    (*m_free) = 0U;
    unsigned short m_sat = 0U;

    for (unsigned short i = 0; i < W->rows; i++)
    {
        if (W->vector[i] == 0)
        {
            permutation->vector[*m_free] = i;
            (*m_free)++;
        }
    }

    for (unsigned short i = 0; i < W->rows; i++)
    {
        if (W->vector[i] != 0)
        {
            permutation->vector[(*m_free) + m_sat] = i;
            m_sat++;
        }
    }

    return status;
}

void get_scaling(const vector_t *const u_delta, const vector_t *const u_alloc, const vector_t *const u_min, const vector_t *const u_max, const vector_u_t *const permutation, const vector_i_t *const W, unsigned short m_free, float *alpha, unsigned short *idx_alpha)
{
    *alpha = 1.0f;
    *idx_alpha = 0U;

    for (unsigned short i = 0; i < m_free; i++)
    {
        if (W->vector[permutation->vector[i]] == 0) // check if input was free // TODO Check if this is necessarry
        {
            float alpha_temp = 1.0f;
            if (u_delta->vector[i] < -TOL)
            {
                alpha_temp = (u_min->vector[permutation->vector[i]] - u_alloc->vector[permutation->vector[i]]) / u_delta->vector[i];
            }
            else if (u_delta->vector[i] > TOL)
            {
                alpha_temp = (u_max->vector[permutation->vector[i]] - u_alloc->vector[permutation->vector[i]]) / u_delta->vector[i];
            }

            if (alpha_temp < *alpha) // Will only be executed if one of the cases above are active
            {
                *alpha = alpha_temp;
                //*idx_alpha = permutation->vector[i];
                *idx_alpha = i;
            }
        }
    }

    if (*alpha < TOL)
    {
        *alpha = 0.0f;
    }
}

unsigned short check_feasibility(const vector_t *const u_opt, const vector_t *const u_delta, const vector_t *const u_alloc, const vector_t *const u_min, const vector_t *const u_max, const vector_u_t *const permutation, unsigned short m_free)
{
    /* u_opt = u_alloc + u_delta */
    for (unsigned short i = 0U; i < u_opt->rows; i++)
    {
        if (i < m_free)
        {
            u_opt->vector[permutation->vector[i]] = u_alloc->vector[permutation->vector[i]] + u_delta->vector[i];
        }
        else
        {
            u_opt->vector[permutation->vector[i]] = 0.0f;
        }
    }

    unsigned short feasible = 1U;
    for (unsigned short i = 0U; i < m_free; i++)
    {
        if ((u_opt->vector[permutation->vector[i]] < (u_min->vector[permutation->vector[i]] - (TOL * 10.0f))) ||
            (u_opt->vector[permutation->vector[i]] > (u_max->vector[permutation->vector[i]] + (TOL * 10.0f))))
        {
            feasible = 0U;
        }
    }
    return (feasible);
}
