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


#ifndef _xf_AllPullStruct_h
#define _xf_AllPullStruct_h 1

/*
  FILE:  xf_AllPullStruct.h
 
  This file contains the definition of structures used in AllPull.c
 
*/

/* Conversion table structure */
typedef struct
{
  int nItem;
  int *Item;
  
} xf_ConvTable;

/* Mesh conversion structure */
typedef struct
{
  int Dim;
  
  xf_ConvTable Nodes;
  
  xf_ConvTable IFaces;
  
  int nBFaceGroup;
  xf_ConvTable *BFaceGroups;
  
  int nElemGroup;
  xf_ConvTable *ElemGroups;
  
} xf_MeshPullInfo;


#endif // end ifndef _xf_AllPullStruct_h
