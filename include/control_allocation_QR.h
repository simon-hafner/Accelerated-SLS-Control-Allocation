/*------------------------------------------------------------------------*/
/*                       control_allocation_QR.h                          */
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
/*Date:         2024-08-06                                                */
/*Contact:      Institute of Flight System Dynamics                       */
/*              School of Engineering and Design                          */
/*              Technische Universitaet Muenchen                          */
/*              Boltzmannstrasse 15                                       */
/*              D-85748 Garching                                          */
/*              simon.hafner@tum.de                                       */
/*------------------------------------------------------------------------*/

#ifndef CONTROL_ALLOCATION_QR_H
#define CONTROL_ALLOCATION_QR_H

#define TOL 1e-5f

#include "mat_vec.h"

unsigned int pivoted_QR_T(matrix_t *const B, vector_u_t *p, vector_t *h, const unsigned short m, float tol, unsigned short *const rank, const unsigned short rank_max);

unsigned int pivoted_QR(matrix_t *const A, vector_u_t *p, vector_t *h, const unsigned short m, const unsigned short n, const float tol, unsigned short *const rank, const unsigned short rank_max);

unsigned int solve_QR_LS_udet(const matrix_t *const B, const vector_c_t *const nu_des, const unsigned short *const rank, vector_t *u_alloc, matrix_t *const R, const vector_u_t *p, const vector_t *h, vector_t *g);

unsigned int solve_QR_LS_udet_SLS(const matrix_t *const B, const vector_t *const nu_des, const unsigned short m, const unsigned short *const rank, vector_t *u_alloc, matrix_t *const R, const vector_u_t *p, const vector_t *h, vector_t *g);

unsigned int solve_QR_LS_odet_SLS(matrix_t *const A, const vector_t *const u_des, const unsigned short m, const unsigned short n, const unsigned short *const rank, vector_t *u_sol, const vector_u_t *p, const vector_t *h);

unsigned int nullspace(matrix_t *const B, vector_u_t *p, vector_t *h, const unsigned short m, const float tol, matrix_t *const Q2, unsigned short *const n_Q2);

unsigned int H2_mat(unsigned short lpivot, unsigned short l1, unsigned char transpose, unsigned short m, const matrix_t *const U, const float *const up, matrix_t *const C, unsigned short m_C);

#endif /* CONTROL_ALLOCATION_QR_H */
