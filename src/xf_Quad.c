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

/*
  FILE:  xf_Quad.c

  This file contains functions for working with quadrature points.

*/


#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Basis.h"
#include "xf_QuadRule.h"
#include "xf_QuadStruct.h"
#include "xf_EqnSetHook.h"
#include "xfYu_Model.h"

/* Customized quadrature order increment for uniformly increasing
   quadrature accuracy */
int QuadOrderDelta = 0;

/* Delta added to requested quadrature orders in the presence of mesh
   motion.  This value will depend on the smoothness of the motion
   mapping.  The value below is a comporomise between accuracy and
   speed. */
int QuadOrderMeshMotionDelta = 1;

/* Custom quadrature order increment proportional to interpolation order p.
The quadrature order used is QuadOrder + p*QuadOrderProportionalDelta. */
int QuadOrderProportionalDelta = 0; 


/******************************************************************/
//   FUNCTION Definition: xf_CreateQuadData
int 
xf_CreateQuadData( xf_QuadData **pQuadData)
{
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pQuadData, 1, sizeof(xf_QuadData)));
  if (ierr != xf_OK) return ierr;

  (*pQuadData)->Type   = xfe_QuadDataGeneric;
  (*pQuadData)->Shape  = xfe_ShapeLast;
  (*pQuadData)->Order  = -1;
  (*pQuadData)->nquad  = 0;
  (*pQuadData)->qdim   = 0;
  (*pQuadData)->xquad  = NULL;
  (*pQuadData)->wquad  = NULL;  
  (*pQuadData)->nvec   = NULL;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyQuadData
int 
xf_DestroyQuadData( xf_QuadData *QuadData)
{
  if (QuadData == NULL) return xf_OK;

  xf_Release( (void *) QuadData->xquad);
  xf_Release( (void *) QuadData->wquad);
  xf_Release( (void *) QuadData->nvec);
  xf_Release( (void *) QuadData);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyGenericQuadData
int 
xf_DestroyGenericQuadData( xf_QuadData *QuadData)
{
  if (QuadData == NULL) return xf_OK;

  if (QuadData->Type == xfe_QuadDataGeneric)
    return xf_Error(xf_DestroyQuadData(QuadData));

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_GetQuadOrderElem
int 
xf_GetQuadOrderElem( const xf_Mesh *Mesh, const xf_EqnSet *EqnSet, 
		     int egrp, int Order, int *QuadOrder)
{
  int ierr;
  int GeomOrder;
  enum xfe_BasisType GeomBasis;

  if (Order < 0){
    (*QuadOrder) = -1; // max quadrature
    return xf_OK;
  }

  // call equation-set-specific function
  if (EqnSet != NULL){
    ierr = xf_Error(xf_EqnSetQuadOrderElem(EqnSet, Order, QuadOrder));
    if (ierr != xf_OK) return ierr;
  }
  else {
    ierr = xf_Error(ElemQuadOrder(Order, QuadOrder));
    if (ierr != xf_OK) return ierr;
  }
  //else (*QuadOrder) = Order;

  GeomOrder = Mesh->ElemGroup[egrp].QOrder;
  GeomBasis = Mesh->ElemGroup[egrp].QBasis;

  /* adjust for the geometry order: the geometry Jacobian will
     introduce additional order into the integrand if Q>1 or if Q1 is
     not linear. */
  
  if (GeomOrder > 1) (*QuadOrder) += Mesh->Dim*(GeomOrder-1);
  if (xf_Q1BasisNotLinear(GeomBasis)) (*QuadOrder) += Mesh->Dim;

  // Are we doing mesh motion? If so, add a hard-coded delta (defined above)
  if ((Mesh->Motion != NULL) && (Mesh->Motion->Active))
    (*QuadOrder) += QuadOrderMeshMotionDelta;

  // add any other hardcoded delta (defined at beginning of file)
  (*QuadOrder) += QuadOrderDelta;
  
  /* Add a hard-coded delta proportional to interpolation order p.
  Only used when EqnSet is provided, since otherwise the input Order 
  may not be p. */
  if (EqnSet != NULL)
    (*QuadOrder) += Order*QuadOrderProportionalDelta;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_GetQuadOrderIFace
int 
xf_GetQuadOrderIFace( const xf_Mesh *Mesh, const xf_EqnSet *EqnSet, 
		      xf_IFace IFace, int Order, int *QuadOrder)
{
  int ierr;
  int GeomOrder, GeomOrderL, GeomOrderR;
  enum xfe_BasisType GeomBasisL, GeomBasisR;
  enum xfe_Bool Q1NotLinear;

  if (Order < 0){
    (*QuadOrder) = -1; // max quadrature
    return xf_OK;
  }

  // call equation-set-specific function
  if (EqnSet != NULL){
    ierr = xf_Error(xf_EqnSetQuadOrderFace(EqnSet, Order, QuadOrder));
    if (ierr != xf_OK) return ierr;
  }
  else {
    ierr = xf_Error(FaceQuadOrder(Order, QuadOrder));
    if (ierr != xf_OK) return ierr;
  }
  //else (*QuadOrder) = Order;

  GeomOrderL = Mesh->ElemGroup[IFace.ElemGroupL].QOrder;
  GeomOrderR = Mesh->ElemGroup[IFace.ElemGroupR].QOrder;
  GeomOrder = max(GeomOrderL, GeomOrderR);

  GeomBasisL = Mesh->ElemGroup[IFace.ElemGroupL].QBasis;
  GeomBasisR = Mesh->ElemGroup[IFace.ElemGroupR].QBasis;
  Q1NotLinear = (xf_Q1BasisNotLinear(GeomBasisL) || xf_Q1BasisNotLinear(GeomBasisR));

  /* adjust for the geometry order: if GeomOrder > 1, the geometry
     Jacobian will introduce additional order into the integrand. */
  if (GeomOrder > 1) (*QuadOrder) += (Mesh->Dim-1)*(GeomOrder-1);
  if (Q1NotLinear) (*QuadOrder) += (Mesh->Dim-1);

  // Are we doing mesh motion? If so, add a hard-coded delta (defined above)
  if ((Mesh->Motion != NULL) && (Mesh->Motion->Active))
    (*QuadOrder) += QuadOrderMeshMotionDelta;

  (*QuadOrder) += QuadOrderDelta;
  
  /* Add a hard-coded delta proportional to interpolation order p.
  Only used when EqnSet is provided, since otherwise the input Order 
  may not be p. */
  if (EqnSet != NULL)
    (*QuadOrder) += Order*QuadOrderProportionalDelta;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_GetQuadOrderAcrossFace
//   simply a modified version of the above one
//   but given left and right info already
int
xf_GetQuadOrderAcrossFace( const xf_Mesh *Mesh, const xf_EqnSet *EqnSet,
                      int egrpL, int egrpR, int Order, int *QuadOrder)
{
    int ierr;
    int GeomOrder, GeomOrderL, GeomOrderR;
    enum xfe_BasisType GeomBasisL, GeomBasisR;
    enum xfe_Bool Q1NotLinear;
    
    if (Order < 0){
        (*QuadOrder) = -1; // max quadrature
        return xf_OK;
    }
    
    // call equation-set-specific function
    if (EqnSet != NULL){
        ierr = xf_Error(xf_EqnSetQuadOrderFace(EqnSet, Order, QuadOrder));
        if (ierr != xf_OK) return ierr;
    }
    else {
        ierr = xf_Error(FaceQuadOrder(Order, QuadOrder));
        if (ierr != xf_OK) return ierr;
    }
    //else (*QuadOrder) = Order;
    
    GeomOrderL = Mesh->ElemGroup[egrpL].QOrder;
    GeomOrderR = Mesh->ElemGroup[egrpR].QOrder;
    GeomOrder = max(GeomOrderL, GeomOrderR);
    
    GeomBasisL = Mesh->ElemGroup[egrpL].QBasis;
    GeomBasisR = Mesh->ElemGroup[egrpR].QBasis;
    Q1NotLinear = (xf_Q1BasisNotLinear(GeomBasisL) || xf_Q1BasisNotLinear(GeomBasisR));
    
    /* adjust for the geometry order: if GeomOrder > 1, the geometry
     Jacobian will introduce additional order into the integrand. */
    if (GeomOrder > 1) (*QuadOrder) += (Mesh->Dim-1)*(GeomOrder-1);
    if (Q1NotLinear) (*QuadOrder) += (Mesh->Dim-1);
    
    // Are we doing mesh motion? If so, add a hard-coded delta (defined above)
    if ((Mesh->Motion != NULL) && (Mesh->Motion->Active))
        (*QuadOrder) += QuadOrderMeshMotionDelta;
    
    (*QuadOrder) += QuadOrderDelta;
    
    /* Add a hard-coded delta proportional to interpolation order p.
     Only used when EqnSet is provided, since otherwise the input Order
     may not be p. */
    if (EqnSet != NULL)
        (*QuadOrder) += Order*QuadOrderProportionalDelta;
    
    return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_GetQuadOrderGeneralFace
int 
xf_GetQuadOrderGeneralFace( const xf_Mesh *Mesh, const xf_EqnSet *EqnSet, 
			    int egrp, int Order, int *QuadOrder)
{
  int ierr;
  int GeomOrder;
  enum xfe_Bool Q1NotLinear;

  if (Order < 0){
    (*QuadOrder) = -1; // max quadrature
    return xf_OK;
  }

  // call equation-set-specific function
  if (EqnSet != NULL){
    ierr = xf_Error(xf_EqnSetQuadOrderFace(EqnSet, Order, QuadOrder));
    if (ierr != xf_OK) return ierr;
  }
  else {
    ierr = xf_Error(FaceQuadOrder(Order, QuadOrder));
    if (ierr != xf_OK) return ierr;
  }
  //else (*QuadOrder) = Order;

  GeomOrder = Mesh->ElemGroup[egrp].QOrder;

  Q1NotLinear = xf_Q1BasisNotLinear(Mesh->ElemGroup[egrp].QBasis);

  /* adjust for the geometry order: if GeomOrder > 1, the geometry
     Jacobian will introduce additional order into the integrand. */
  if (GeomOrder > 1) (*QuadOrder) += (Mesh->Dim-1)*(GeomOrder-1);
  if (Q1NotLinear) (*QuadOrder) += (Mesh->Dim-1);
  
  // Are we doing mesh motion? If so, add a hard-coded delta (defined above)
  if ((Mesh->Motion != NULL) && (Mesh->Motion->Active))
    (*QuadOrder) += QuadOrderMeshMotionDelta;

  (*QuadOrder) += QuadOrderDelta;
  
  /* Add a hard-coded delta proportional to interpolation order p.
  Only used when EqnSet is provided, since otherwise the input Order 
  may not be p. */
  if (EqnSet != NULL)
    (*QuadOrder) += Order*QuadOrderProportionalDelta;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_GetQuadOrderBFace
int 
xf_GetQuadOrderBFace( const xf_Mesh *Mesh, const xf_EqnSet *EqnSet, 
		      xf_BFace BFace, int Order, int *QuadOrder)
{
  return xf_Error(xf_GetQuadOrderGeneralFace(Mesh, EqnSet, BFace.ElemGroup, 
					     Order, QuadOrder));
}




/******************************************************************/
//   FUNCTION Definition: xf_GetQuadPoints
int 
xf_GetQuadPoints( xf_QuadData *QuadData)
{
  int ierr;
  int Order;

  Order = QuadData->Order;
  
  switch (QuadData->Shape){
  case xfe_Point:
    QuadData->qdim = 1;
    QuadData->nquad = 1;
    ierr = xf_Error(xf_ReAlloc( (void **) &QuadData->xquad, 1, sizeof(real)));
    if (ierr != xf_OK) return ierr;  
    QuadData->xquad[0] = 0;
    ierr = xf_Error(xf_ReAlloc( (void **) &QuadData->wquad, 1, sizeof(real)));
    if (ierr != xf_OK) return ierr;  
    QuadData->wquad[0] = 1.0;
    break;
  case xfe_Segment:
    QuadData->qdim = 1;
    return xf_QuadLine(Order, &QuadData->nquad, &QuadData->xquad, 
		       &QuadData->wquad);
    break;
  case xfe_Triangle:
    QuadData->qdim = 2;
    return xf_QuadTriangle(Order, &QuadData->nquad, &QuadData->xquad, 
			   &QuadData->wquad);
    break;
  case xfe_Quadrilateral:
    QuadData->qdim = 2;
    return xf_QuadQuadrilateral(Order, &QuadData->nquad, &QuadData->xquad, 
				&QuadData->wquad);
    break;
  case xfe_Tetrahedron:
    QuadData->qdim = 3;
    return xf_QuadTetrahedron(Order, &QuadData->nquad, &QuadData->xquad, 
			      &QuadData->wquad);
    break;
  case xfe_Hexahedron:
    QuadData->qdim = 3;
    return xf_QuadHexahedron(Order, &QuadData->nquad, &QuadData->xquad, 
			     &QuadData->wquad);
    break;
  default:
    return xf_Error(xf_UNKNOWN_SHAPE);
  }

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_QuadElem
int 
xf_QuadElem( xf_Mesh *Mesh, int egrp, int elem, int Order, 
	     xf_QuadData **pQuadData, enum xfe_Bool *changed)
{
  int ierr;
  enum xfe_ShapeType Shape;
  xf_QuadData *QuadData;
  
  if (Mesh->ElemGroup[egrp].CutFlag){
    ierr = xf_Error(xf_DestroyGenericQuadData((*pQuadData)));
    if (ierr != xf_OK) return ierr;
    (*pQuadData) = Mesh->ElemGroup[egrp].CutElemData[elem].QuadData;
    (*changed) = xfe_True;
    return xf_OK;
  }

  /* At this point, the element is not cut */

  // if (*pQuadData) is pointing to specific quad data, set it to NULL
  if ((*pQuadData) != NULL)
    if ((*pQuadData)->Type == xfe_QuadDataSpecific) (*pQuadData) = NULL;
  
  /* Determine the shape of the element (triangle, tetrahedron, etc.) */
  ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
  if (ierr != xf_OK) return ierr;

  if ((*pQuadData) == NULL){
    (*changed) = xfe_True;
    ierr = xf_Error(xf_CreateQuadData(pQuadData));
    if (ierr != xf_OK) return ierr;
    (*pQuadData)->Type = xfe_QuadDataGeneric;
  }
  else{
    // (*pQuadData) is already pointing to some generic rule    
    if ( ((*pQuadData)->Shape == Shape) && ((*pQuadData)->Order == Order)){
      (*changed) = xfe_False;
      return xf_OK; // the existing quad rule is ok
    }
    (*changed) = xfe_True;
  }

  QuadData = (*pQuadData);
  QuadData->Shape = Shape;
  QuadData->Order = Order;

  // get generic quad points; reallocate the arrays in QuadData
  ierr = xf_Error(xf_GetQuadPoints(QuadData));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_QuadElemAtLeast
int 
xf_QuadElemAtLeast( xf_Mesh *Mesh, int egrp, int elem, int Order, int nmin, 
		    xf_QuadData **pQuadData, enum xfe_Bool *changed)
{
  int ierr, nq, i;
  int imax = 20; // number of times we'll up the order before we give up

  i = 0;
  nq = -1;
  while (nq <= nmin){
    ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, Order, pQuadData, changed));
    if (ierr != xf_OK) return ierr;
    if ((nq = (*pQuadData)->nquad) <= nmin) Order++;
    // we've maxed out on quad capability and haven't gotten enough points
    if (i++ > imax) return xf_Error(xf_NOT_SUPPORTED);
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_QuadElemJustAbove
int 
xf_QuadElemJustAbove( xf_Mesh *Mesh, int egrp, int elem, int nmin, 
		      xf_QuadData **pQuadData, enum xfe_Bool *changed)
{
  int ierr, nq, i;
  int Order;
  int imax = 20; // number of times we'll up the order before we give up

  i = 0;
  nq = -1;
  Order = 1;
  while (nq <= nmin){
    ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, Order, pQuadData, changed));
    if (ierr != xf_OK) return ierr;
    if ((nq = (*pQuadData)->nquad) <= nmin) Order++;
    // we've maxed out on quad capability and haven't gotten enough points
    if (i++ > imax) return xf_Error(xf_NOT_SUPPORTED);
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_QuadFace
int 
xf_QuadFace( xf_Mesh *Mesh, int egrp, int elem, int face, int Order, 
	      xf_QuadData **pQuadData, enum xfe_Bool *changed)
{
  int ierr;
  enum xfe_BasisType QBasis;
  enum xfe_ShapeType ElemShape, FaceShape;
  xf_Face Face;
  xf_CutFaceData *CutFaceData;
  xf_QuadData *QuadData;
  
  Face = Mesh->ElemGroup[egrp].Face[elem][face];

  if (Face.Group == xf_INTERIORFACE) // Interior face
    CutFaceData = Mesh->IFace[Face.Number].CutFaceData;
  else if (Face.Group >= 0)  // Boundary face
    CutFaceData = Mesh->BFaceGroup[Face.Group].BFace[Face.Number].CutFaceData;
  else return xf_Error(xf_NOT_SUPPORTED);

  if (CutFaceData != NULL){
    ierr = xf_Error(xf_DestroyGenericQuadData((*pQuadData)));
    if (ierr != xf_OK) return ierr;
    (*pQuadData) = CutFaceData->QuadData;
    (*changed) = xfe_True;
    return xf_OK;
  }

  /* At this point, the face is not cut */
  
  QBasis = Mesh->ElemGroup[egrp].QBasis;

  /* Determine the shape of the face (line, triangle, quad, etc.) */
  ierr = xf_Error(xf_Basis2Shape(QBasis, &ElemShape));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_FaceShape(ElemShape, face, &FaceShape));
  if (ierr != xf_OK) return ierr;


  // if (*pQuadData) is pointing to specific quad data, set it to NULL
  if ((*pQuadData) != NULL)
    if ((*pQuadData)->Type == xfe_QuadDataSpecific) (*pQuadData) = NULL;

  if ((*pQuadData) == NULL){
    (*changed) = xfe_True;
    ierr = xf_Error(xf_CreateQuadData(pQuadData));
    if (ierr != xf_OK) return ierr;
    (*pQuadData)->Type = xfe_QuadDataGeneric;
  }
  else{
    // (*pQuadData) is already pointing to some generic rule    
    if ( ((*pQuadData)->Shape == FaceShape) && ((*pQuadData)->Order == Order)){
      (*changed) = xfe_False;
      return xf_OK; // the existing quad rule is ok
    }
    (*changed) = xfe_True;
  }

  QuadData = (*pQuadData);
  QuadData->Shape = FaceShape;
  QuadData->Order = Order;

  // get generic quad points; reallocate the arrays in QuadData
  ierr = xf_Error(xf_GetQuadPoints(QuadData));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_WriteQuadDataBinary
int 
xf_WriteQuadDataBinary( xf_QuadData *QuadData, FILE *fid){
  int ierr, tot, rev, si;
  enum xfe_Bool flag;

  si = sizeof(int);

  rev = 0;  // writer revision number
  if (fwrite(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);

  ierr = xf_Error(xf_WriteStringBinary(xfe_QuadDataName[QuadData->Type], fid));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_WriteStringBinary(xfe_ShapeName[QuadData->Shape], fid));
  if (ierr != xf_OK) return ierr;
  
  ierr = 0;
  ierr += fwrite(&QuadData->Order, si, 1, fid);
  ierr += fwrite(&QuadData->nquad , si, 1, fid);
  ierr += fwrite(&QuadData->qdim  , si, 1, fid);
  if (ierr != 3) return xf_Error(xf_FILE_WRITE_ERROR);

  // xquad
  flag = (QuadData->xquad != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    tot = QuadData->qdim*QuadData->nquad;
    if (fwrite(QuadData->xquad, sizeof(real), tot, fid) != tot)
      return xf_Error(xf_FILE_WRITE_ERROR);
  }

  // wquad
  flag = (QuadData->wquad != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    tot = QuadData->nquad;
    if (fwrite(QuadData->wquad, sizeof(real), tot, fid) != tot)
      return xf_Error(xf_FILE_WRITE_ERROR);
  }

  // nvec
  flag = (QuadData->nvec != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    tot = QuadData->qdim*QuadData->nquad;
    if (fwrite(QuadData->nvec, sizeof(real), tot, fid) != tot)
      return xf_Error(xf_FILE_WRITE_ERROR);
  }
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_ReadQuadDataBinary
int 
xf_ReadQuadDataBinary( FILE *fid, xf_QuadData *QuadData){
  int ierr, tot, rev, si;
  enum xfe_Bool flag;

  si = sizeof(int);

  // read + check revision number
  rev = 0;
  if (fread(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);

  // Type of Quad data
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_QuadDataName, xfe_QuadDataLast, 
				    (int *) &QuadData->Type));
  if (ierr != xf_OK) return ierr;

  // Shape
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_ShapeName, xfe_ShapeLast, 
				    (int *) &QuadData->Shape));
  if (ierr != xf_OK) return ierr;
  
  ierr = 0;
  ierr += fread(&QuadData->Order, si, 1, fid);
  ierr += fread(&QuadData->nquad, si, 1, fid);
  ierr += fread(&QuadData->qdim , si, 1, fid);
  if (ierr != 3) return xf_Error(xf_FILE_READ_ERROR);

  // xquad
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    tot = QuadData->qdim*QuadData->nquad;
    ierr = xf_Error(xf_Alloc((void **) QuadData->xquad, tot, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    if (fread(QuadData->xquad, sizeof(real), tot, fid) != tot)
      return xf_Error(xf_FILE_READ_ERROR);
  }
  else
    QuadData->xquad = NULL;

  // wquad
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    tot = QuadData->nquad;
    ierr = xf_Error(xf_Alloc((void **) QuadData->wquad, tot, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    if (fread(QuadData->wquad, sizeof(real), tot, fid) != tot)
      return xf_Error(xf_FILE_READ_ERROR);
  }
  else
    QuadData->wquad = NULL;

  // nvec
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    tot = QuadData->qdim*QuadData->nquad;
    ierr = xf_Error(xf_Alloc((void **) QuadData->nvec, tot, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    if (fread(QuadData->nvec, sizeof(real), tot, fid) != tot)
      return xf_Error(xf_FILE_READ_ERROR);
  }
  else
    QuadData->nvec = NULL;
  return xf_OK;
}


