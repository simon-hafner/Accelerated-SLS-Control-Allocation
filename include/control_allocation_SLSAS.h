/*------------------------------------------------------------------------*/
/*                     control_allocation_SLSAS.h                         */
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
/*Type:         C - header file                                           */
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

#ifndef CONTROL_ALLOCATION_SLSAS_H
#define CONTROL_ALLOCATION_SLSAS_H

#define TOL 1e-5f

#include "mat_vec.h"
#include "control_allocation_QR.h"

unsigned int CA_SLSAS(float *const B, float *const W_u, float *const W_nu, float *const nu_des, float *const u_min, float *const u_max,
                      float *const u_des, short *const W, float *const u_alloc, unsigned short *const rank_B, const unsigned short n, const unsigned short m, const unsigned short n_W_u,
                      const unsigned short max_iter, const float abort_tol, const unsigned short warmstart, float *const buffer_float, unsigned short *const buffer_unsigned);

unsigned int get_permutations(const vector_i_t *const W, vector_u_t *const permutation, unsigned short *m_free);

void get_scaling(const vector_t *const u_delta, const vector_t *const u_alloc, const vector_t *const u_min, const vector_t *const u_max, const vector_u_t *const permutation, const vector_i_t *const W, unsigned short m_free, float *alpha, unsigned short *idx_alpha);

unsigned short check_feasibility(const vector_t *const u_opt, const vector_t *const u_delta, const vector_t *const u_alloc, const vector_t *const u_min, const vector_t *const u_max, const vector_u_t *const permutation, unsigned short m_free);

#endif /* CONTROL_ALLOCATION_SLSAS_H */
