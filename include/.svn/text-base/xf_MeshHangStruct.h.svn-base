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

#ifndef _xf_MeshHangStruct_h
#define _xf_MeshHangStruct_h 1

/*
  FILE:  xf_MeshHangStruct.h

  This file contains structures for the xflow hanging face support

*/


/* Types of positions for segments */
enum xfe_SegPosType{
  xfe_SegPosNone,
  xfe_SegPosLeft,
  xfe_SegPosRight,
  xfe_SegPosLast
};

/* corresponding names */
static char *xfe_SegPosName[xfe_SegPosLast] = {
  "None",
  "Left",
  "Right"
};

/* Types of positions for quads.  Do not change this order */
enum xfe_QuadPosType{
  xfe_QuadPosNone,
  xfe_QuadPosSW,
  xfe_QuadPosSE,
  xfe_QuadPosNW,
  xfe_QuadPosNE,
  xfe_QuadPosLeft,
  xfe_QuadPosRight,
  xfe_QuadPosBottom,
  xfe_QuadPosTop,
  xfe_QuadPosLast
};

/* corresponding names */
static char *xfe_QuadPosName[xfe_QuadPosLast] = {
  "None",
  "SW",
  "SE",
  "NW",
  "NE",
  "Left",
  "Right",
  "Bottom",
  "Top",
};


/* Types of positions for triangles */
enum xfe_TriPosType{
  xfe_TriPosNone,
  xfe_TriPosLeft,
  xfe_TriPosRight,
  xfe_TriPosTop,
  xfe_TriPosCenter,
  xfe_TriPosLast
};

/* corresponding names */
static char *xfe_TriPosName[xfe_TriPosLast] = {
  "None",
  "Left",
  "Right",
  "Top",
  "Center"
};



/* Types of positions for hexes.  Do not change this order */
enum xfe_HexPosType{
  xfe_HexPosNone,

  xfe_HexPos000, // uniform, subelem adjacent to origin
  xfe_HexPos100, // uniform, subelem adjacent to point 1,0,0
  xfe_HexPos010, // ...
  xfe_HexPos110,
  xfe_HexPos001,
  xfe_HexPos101,
  xfe_HexPos011,
  xfe_HexPos111, // uniform, subelem adjacent to 1,1,1

  xfe_HexPos022, // SliceX, left
  xfe_HexPos122, // SliceX, right
  xfe_HexPos202, // SliceY, front
  xfe_HexPos212, // SliceY, back
  xfe_HexPos220, // SliceZ, bottom
  xfe_HexPos221, // SliceZ, top

  xfe_HexPos002, // SliceXY, adjacent to 0,0,0
  xfe_HexPos102, // SliceXY, adjacent to 1,0,0
  xfe_HexPos012, // SliceXY, adjacent to 0,1,0
  xfe_HexPos112, // SliceXY, adjacent to 1,1,0
  xfe_HexPos020, // SliceXZ, adjacent to 0,0,0
  xfe_HexPos120, // SliceXZ, adjacent to 1,0,0
  xfe_HexPos021, // SliceXZ, adjacent to 0,0,1
  xfe_HexPos121, // SliceXZ, adjacent to 1,0,1
  xfe_HexPos200, // SliceYZ, adjacent to 0,0,0
  xfe_HexPos210, // SliceYZ, adjacent to 0,1,0
  xfe_HexPos201, // SliceYZ, adjacent to 0,0,1
  xfe_HexPos211, // SliceYZ, adjacent to 0,1,1

  xfe_HexPosLast
};



#endif // end ifndef _xf_MeshHangStruct_h
