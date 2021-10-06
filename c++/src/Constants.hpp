/*!
 * \file Constants.hpp
 * \author S. Buchet
 * \brief definition of constants
 */

#pragma once

#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

/* set to 1 to output debugging information */
#define DEBUG_LOG 1                               

/* set to 1 to compute only few genes */
#define DEBUG_REDUCTION 0
/* number of variables to compute if DEBUG_REDUCTION == 1 */
#define NVAR_DEBUG 20

/* true if only a sub sample of the data is used */
#define SAMPLING false
/* between 0 and 1, percentage of cells to compute */
#define SAMPLE_SIZE 0.1

/* set to 1 to use multiple cores */
#define USE_OPENMP 1

/* number of processors to use */
#if USE_OPENMP == 1
  #define N_THREADS 6
  //#define N_THREADS 16
#endif

#endif
