/*------------------------------------------------------------------------*/
/*                     control_allocation_RSPI.h                         */
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

#ifndef CONTROL_ALLOCATION_RSPI_H
#define CONTROL_ALLOCATION_RSPI_H

#define TOL 1e-5f

#define NUM_U 12U
#define NUM_NU 6U
#define BUFFER_R_F_SIZE (3U * NUM_U * NUM_U) + (NUM_NU * NUM_NU) + (NUM_U * NUM_NU) + (4U * NUM_U)
#define BUFFER_R_U_SIZE 2U * NUM_U

#include "mat_vec.h"
#include "control_allocation_SLSAS.h"

unsigned int CA_RSPI(float *const B, float *const W_u, float *const W_nu, float *const nu_des, float *const u_min, float *const u_max,
                     float *const u_des, short *const W, float *const u_alloc, const unsigned short n, const unsigned short m, const unsigned short n_W_u,
                     const unsigned short max_iter, const float abort_tol, float *const buffer_float, unsigned short *const buffer_unsigned);

#endif /* CONTROL_ALLOCATION_RSPI_H */
