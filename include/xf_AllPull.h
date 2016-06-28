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


#ifndef _xf_AllPull_h
#define _xf_AllPull_h 1

/*
  FILE:  xf_AllPull.h
 
  This file contains the headers for functions in AllPull.c
 
*/

/******************************************************************/
//   FUNCTION Definition: xf_CheckExist
int
xf_CheckExist(int num, int size, int *vector, int *index);

/*
 
 PURPOSE: 
 
 Check existance of num in vector.
 
 INPUTS:
 
 num : value to be checked
 size : size of vector
 vector : vector of integers among which num will be checked
 out: -1 -> vector does not have num; >=0 -> vector has num at index
 
 OUTPUTS: 
 None. value of index is changed
 
 RETURNS: Error code
 
 
 */


/******************************************************************/
//   FUNCTION Definition: xf_AllPull
int
xf_AllPull(xf_All *All, xf_All *All_Small, int egrp, int elem);
/*
 
  PURPOSE: 
 
  Copy the local All Structure (elements and 1st-level neighbors) to All_Small
  The central element will always have the same element group as the large mesh
  and its number is always 0
 
  INPUTS:
 
  All : All Structure
  All_Small : All_Small structure
  elem : element to pull local mesh structure from
  egrp : element group to which elem belongs
 
  OUTPUTS: 
  None: All_Small structure is modified.
 
 
  RETURNS: Error Code
 
*/



#endif // end ifndef _xf_AllPull_h
