/*------------------------------------------------------------------------*/
/*                               mat_vec.h                                 */
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
/*Description:  Provides basic matrix functionilities required            */
/*Type:         C - header file                                           */
/*Dependencies:                                                           */
/*------------------------------------------------------------------------*/
/*Author:       S. Hafner                                                 */
/*Date:         2023-10-24                                                */
/*Contact:      Institute of Flight System Dynamics                       */
/*              School of Engineering and Design                          */
/*              Technische Universitaet Muenchen                          */
/*              Boltzmannstrasse 15                                       */
/*              D-85748 Garching                                          */
/*              simon.hafner@tum.de                                       */
/*------------------------------------------------------------------------*/

#ifndef MATRIX_H
#define MATRIX_H

// #define MAXROWS 6
/* Specific to project and must be adapted*/

#define EPS_F 1.19e-07f

#define ABS(X) (((X) >= (float)0.0) ? (X) : -(X)) /* Absolutbetrag von X */

#define MIN(X, Y) (((X) <= (Y)) ? (X) : (Y)) /* Minimum of X and Y */
#define MAX(X, Y) (((X) >= (Y)) ? (X) : (Y)) /* Maximum of X and Y  */

#define M_P_ENTRY(mat, i, j) ((mat)->matrix[((j) * (mat)->rows) + (i)])
#define M_ENTRY(mat, i, j) ((mat).matrix[((j) * (mat).rows) + (i)])

// #include <math.h>

/*Variable definitions*/
/*(C) S. Myschik, OCT 2005*/

/**
 * @brief Matrix struct containing matrix, num rows, num columns
 *
 */
typedef struct
{
    unsigned short rows;
    unsigned short columns;
    float *matrix;
} matrix_t;

/**
 * @brief Constant Matrix struct containing matrix, num rows, num columns
 *
 */
typedef struct
{
    const unsigned short rows;
    const unsigned short columns;
    const float *matrix;
} matrix_c_t;

/*(C) S. Hafner. OCT 2023*/

/**
 * @brief Vector struct containing matrix, num rows
 *
 */
typedef struct
{
    unsigned short rows;
    float *vector;
} vector_t;

/**
 * @brief Constant Vector struct containing matrix, num rows
 *
 */
typedef struct
{
    const unsigned short rows;
    const float *vector;
} vector_c_t;

/**
 * @brief Vector struct containing matrix, num rows
 *
 */
typedef struct
{
    unsigned short rows;
    unsigned short *vector;
} vector_u_t;

/**
 * @brief Vector struct containing matrix, num rows
 *
 */
typedef struct
{
    unsigned short rows;
    short *vector;
} vector_i_t;

#endif /* MATRIX */