/*------------------------------------------------------------------*/
/* XFLOW: A discontinuous Galerkin finite element software library. */
/*                                                                  */
/*                    Copyright  2007-2008                          */
/*           Krzysztof J. Fidkowski, kfid@alum.mit.edu              */
/*                                                                  */
/*                    Copyright  2008-2012                          */
/*                 The University of Michigan                       */
/*                    All rights reserved                           */
/*                                                                  */
/* This library is intended to be useful but is distributed without */
/* any warranty, not even merchantability or fitness for a          */
/* particular purpose.  It is free software: you can redistribute   */
/* it and/or modify it under the terms of the GNU Lesser General    */
/* Public License (LGPLv3).                                         */
/*                                                                  */
/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free       */
/* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.        */
/*------------------------------------------------------------------*/

#ifndef _xf_MPIStruct_h
#define _xf_MPIStruct_h 1

/*
  FILE:  xf_MPIStruct.h

  This file contains structures and types specific to the message
  passing interface.

*/


//#include <mpi.h>


/* Operation definitions for calling xf_MPI_Reduce */
enum xfe_MPI_Op { 
  xfe_MPI_MIN,
  xfe_MPI_MAX,
  xfe_MPI_SUM,
  xfe_MPI_OpLast
};

#endif // end ifndef _xf_MPIStruct_h
