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
 FILE:  xf_MeshToolsElemSearch.c
 
 This file contains functions for efficiently locating an element,
 starting from global coordinates.
 
 */




/******************************************************************/
//   FUNCTION Definition: xf_BuildElemBoundBox
static int 
xf_BuildElemBoundBox(xf_All *All, xf_Vector **pElemBoundBox)
{
/*
PURPOSE:

  Creates a vector of bounding boxes for each element.

INPUTS:
 
  All   : All structure

OUTPUTS: 

  (*pBoundBox) : vector of bounding boxes
                 [xmin xmax ymin ymax (zmin zmax)]

RETURN:

  Error Code
*/

  int ierr, StateRank, dim;
  int egrp, elem, d, n, nn;
  enum xfe_BasisType QBasis;
  enum xfe_Bool PointsChanged;
  const char Title[]="ElemBoundBox";
  real *xn = NULL, *xglob = NULL;
  real *BB, xmid, xdel;
  real fac = 1.5; // hard-coded safety factor
  xf_BasisData *GeomPhiData = NULL;

  xf_Mesh *Mesh;
  Mesh = All->Mesh;
  dim = Mesh->Dim;

  /* Create a vector for storing bounding boxes */

  StateRank = 2*All->Mesh->Dim;   // StateRank is 2 * dim

  ierr = xf_Error(xf_FindVector(All, Title, xfe_LinkageGlobElem, StateRank, NULL, 0, 0, 
				NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False, 
				xfe_False, NULL, pElemBoundBox, NULL));
  if (ierr != xf_OK) return ierr;

  /* Loop over elements and create bounding boxes */
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    // basis of element group
    QBasis = Mesh->ElemGroup[egrp].QBasis;
    // Lagrange nodes in reference space
    ierr = xf_Error(xf_LagrangeNodes(QBasis, Mesh->ElemGroup[egrp].QOrder, &nn, NULL, &xn));
    if (ierr != xf_OK) return ierr;
    // re-allocate global coordinates
    ierr = xf_Error(xf_Alloc((void **) &xglob, nn*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    PointsChanged = xfe_True;

    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

      // global coordinates
      ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, 
				      PointsChanged, nn, xn, xglob));
      if (ierr != xf_OK) return ierr;
      PointsChanged = xfe_False;

      if (nn==0) return xf_Error(xf_CODE_LOGIC_ERROR);

      // take min/max to obtain bounding box
      for (d=0; d<dim; d++){
	// pointer to bounding box data
	BB = (*pElemBoundBox)->GenArray[egrp].rValue[elem]+2*d;
	BB[0] = BB[1] = xglob[d];
	for (n=1; n<nn; n++){
	  BB[0] = min(BB[0], xglob[n*dim+d]);
	  BB[1] = max(BB[1], xglob[n*dim+d]);
	} // n

	// apply safety factor to account for nonlinearity
	xmid = 0.5*(BB[0]+BB[1]);
	xdel = 0.5*(BB[1]-BB[0]);
	BB[0] = xmid-fac*xdel;
	BB[1] = xmid+fac*xdel;
	
	// print out for debugging
	//xf_printf("(%d, %d): (%.6E, %.6E)\n", egrp, elem, BB[0], BB[1]);
       
      } // d

    } // elem
  } // egrp

  /* Destroy Geometry Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  xf_Release( (void *) xn);
  xf_Release( (void *) xglob);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyElemSearchStructure
int
xf_DestroyElemSearchStructure(xf_ElemSearchStruct *ESS)
{
  int d, iBin;
  
  for (d=0; d<ESS->dim; d++){
    for (iBin=0; iBin<ESS->nBin[d]; iBin++){
      xf_Release( (void *) ESS->Elem[d][iBin]);
    }
    xf_Release( (void *) ESS->nElem[d]);
    xf_Release( (void *) ESS->xBin[d]);
    xf_Release( (void *) ESS->Elem[d]);    
  }
  xf_Release( (void *) ESS->buf);

  return xf_OK;
}

  
/******************************************************************/
//   FUNCTION Definition: xf_BuildElemSearchStructure
int
xf_BuildElemSearchStructure(xf_All *All, xf_ElemSearchStruct *ESS)
{
  int ierr, dim, d;
  int nelemtot, nbuf;
  int egrp, elem, eglob, nelem;
  int nLevel, nLevelMax = 5;
  int iBin, nBin, Delta, Break;
  int *pos = NULL;
  enum xfe_Verbosity Verbosity;
  real *xmid = NULL;
  real *BB = NULL;
  real domainBB[6];
  xf_Vector *ElemBoundBox = NULL;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  dim = Mesh->Dim;

  // determine verbosity
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
				     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
				     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;
  
  // number of elements
  ierr = xf_Error(xf_GetnElem(Mesh, NULL, &nelemtot));
  if (ierr != xf_OK) return ierr;

  // set number of levels based on dim and number of elements
  nLevel = max(1, (int) log(pow(nelemtot, 1.0/dim))/log(2));
  nLevel = min(nLevel, nLevelMax);

  if (Verbosity > xfe_VerbosityMedium)
    xf_printf("nLevel = %d\n", nLevel);

  // number of bins is 2^nLevel
  nBin = xf_PowInt(2, nLevel);

  // allocate memory and initialize ESS
  ESS->dim = dim;
  for (d=0; d<dim; d++){
    ESS->nBin[d] = nBin;
    ierr = xf_Error(xf_Alloc((void **) &ESS->xBin[d], nBin+1, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &ESS->nElem[d], nBin, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &ESS->Elem[d], nBin, sizeof(int *)));
    if (ierr != xf_OK) return ierr;
  }
  nbuf = 0; // counter for buffer size

  // create element bounding boxes  
  ierr = xf_Error(xf_BuildElemBoundBox(All, &ElemBoundBox));
  if (ierr != xf_OK) return ierr;

  // allocate temporary memory
  ierr = xf_Error(xf_Alloc((void **) &xmid, nelemtot, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) &pos, nelemtot, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // start loop over dimensions
  for (d=0; d<dim; d++){	       
    // create processor-local list of elem centroids
    domainBB[2*d+0] = 1e30; domainBB[2*d+1] = -1e30;
    for (egrp=0, eglob=0; egrp<Mesh->nElemGroup; egrp++)
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
	BB = ElemBoundBox->GenArray[egrp].rValue[elem]+2*d;
	xmid[eglob++] = 0.5*(BB[0] + BB[1]);
	domainBB[2*d+0] = min(domainBB[2*d+0], BB[0]);
	domainBB[2*d+1] = max(domainBB[2*d+1], BB[1]);
      }
    
    if (eglob != nelemtot) return xf_Error(xf_CODE_LOGIC_ERROR);

    // sort processor-local list of centroids
    ierr = xf_Error(xf_SortRealParallel(xmid, nelemtot, xfe_True, pos));
    if (ierr != xf_OK) return ierr;

    // set bin delta, breakpoint
    Delta = nelemtot/nBin;
    Break = nBin - nelemtot%nBin;

    // loop over interval bins
    eglob = 0;
    ESS->xBin[d][0] = domainBB[2*d+0];
    for (iBin=0; iBin<nBin; iBin++){
      // increment eglob = global element number
      if (iBin == Break) Delta++;
      eglob += Delta;
      //xf_printf("d = %d, iBin = %d, Delta = %d, Break = %d\n",
      //	d, iBin, Delta, Break);
      // set xBin
      if (iBin == nBin-1) ESS->xBin[d][iBin+1] = domainBB[2*d+1];
      else ESS->xBin[d][iBin+1] = xmid[eglob];
      // count number of elements whose bbox overlaps interval
      for (egrp=0, nelem=0; egrp<Mesh->nElemGroup; egrp++)
	for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
	  BB = ElemBoundBox->GenArray[egrp].rValue[elem]+2*d;
	  if ((BB[0] < ESS->xBin[d][iBin+1]) && (BB[1] > ESS->xBin[d][iBin])){
	    nelem++;
	  }
	}
      ESS->nElem[d][iBin] = nelem;
      // allocate Elem[d][bin][*]
      ierr = xf_Error(xf_Alloc((void **) ESS->Elem[d]+iBin, 2*nelem, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      // keep track of maximum for buffer
      nbuf = max(nbuf, nelem);
      // fill in SearchElem by looping over elems again
      for (egrp=0, nelem=0; egrp<Mesh->nElemGroup; egrp++)
	for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
	  BB = ElemBoundBox->GenArray[egrp].rValue[elem]+2*d;
	  if ((BB[0] < ESS->xBin[d][iBin+1]) && (BB[1] > ESS->xBin[d][iBin])){
	    ESS->Elem[d][iBin][2*nelem+0] = egrp;
	    ESS->Elem[d][iBin][2*nelem+1] = elem;
	    //xf_printf("d=%d, iBin=%d, added(%d): %d,%d\n", d, iBin, nelem,egrp,elem);
	    nelem++;
	  }
	}
      
    } // end loop over bins

    if (eglob != nelemtot){
      xf_printf("eglob = %d, nelemtot = %d\n", eglob, nelemtot);
      return xf_Error(xf_CODE_LOGIC_ERROR);
    }

  }// end loop over dimensions
	

  // allocate buffer for later use of search vectors	       
  ierr = xf_Error(xf_Alloc((void **) &ESS->buf, 2*nbuf, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // destroy element bounding boxes
  ierr = xf_Error(xf_DestroyVector(ElemBoundBox, xfe_True));
  if (ierr != xf_OK) return ierr;

  // release memory
  xf_Release( (void *) xmid);
  xf_Release( (void *) pos);

  return xf_OK;
}


  
/******************************************************************/
//   FUNCTION Definition: xf_FindElemUsingSearchStructure
//   modified by Yu to accommedate more grid points (Aug 2015)
int
xf_FindElemUsingSearchStructure(xf_All *All, const int np, real *gxpoint,
                                xf_ElemSearchStruct *ESS, int *pegrp,
                                int *pelem, real *xref)
{
  int ierr, dim, d, dmax, ip;
  int ibin, nbin, dbin, dir;
  int ielem, nelem, nelemmax, egrp, elem;
  int egrpLocal, elemLocal;
  int vBin[3]; // vector of bins
  int I[3], storefound;
  int *PE[3], *infotank;
  enum xfe_Bool done, found, agree;
  enum xfe_Bool converged, inside, OnProc;
  enum xfe_ShapeType Shape;
  real tol = 1e-8, *xpoint;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  dim = Mesh->Dim;
    
  infotank = NULL;
  ierr = xf_Error(xf_Alloc((void **) &infotank, np, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  //xf_printf("\nPoint = %.6E %.6E\n", xpoint[0], xpoint[1]);

  
for(ip=0; ip<np; ip++){
    xpoint = gxpoint + ip*dim;
    
  // locate bins containing xpoint in all dimensions
  for (d=0; d<dim; d++){
    nbin = ESS->nBin[d];  // power of 2
    ibin = nbin/2; dbin = nbin/2; 
    done = xfe_False;
    //xf_printf("nbin = %d, ibin = %d, dbin = %d\n", nbin, ibin, dbin);
    while (dbin > 0){
      if (xpoint[d] < ESS->xBin[d][ibin]) dir = -1;
      else if (xpoint[d] > ESS->xBin[d][ibin+1]) dir = 1;
      else {
	done = xfe_True;
	break; // on interval, we're done
      }
      if ((dbin == 1) && (dir == -1)) ibin--;
      dbin = dbin/2;
      ibin = ibin + dir*dbin;
    }
    if ((xpoint[d] >= ESS->xBin[d][ibin  ]) &&
	(xpoint[d] <= ESS->xBin[d][ibin+1]))
      done = xfe_True;

    OnProc = done;
      
    //this is sanity check; can be avoided
   /*
    // Reduce "done" among all processors .. error if not on any bins
    ierr = xf_Error(xf_MPI_Allreduce((int *) &done, 1, xfe_SizeInt, xfe_MPI_MAX));
    if (ierr != xf_OK) return ierr;

    if (!done){
      xf_printf("Could not find bin containing xpoint[%d] = %.6E on any processors.\n",
		d, xpoint[d]);
      return xf_NOT_FOUND;
    }
    */
      
    vBin[d] = ibin;
  } // d
    
    
  nelem = 0;
  if (OnProc){ // only do the following if point is on this processor

    // using ESS->buf, build elem list = intersection of search lists
    for (d=0; d<dim; d++) I[d] = 0; // starting point in all lists
    done = xfe_False;
    while (!done){
      // done if exceeded any element lists
      done = xfe_False;
      for (d=0; d<dim; d++) 
	if (I[d] >= ESS->nElem[d][vBin[d]]){
	  done = xfe_True;
	  break;
	}
      if (done) break;
      // pointers to elements
      for (d=0; d<dim; d++)
	PE[d] = ESS->Elem[d][vBin[d]] + 2*I[d];
      //xf_printf("I=[%d %d], (%d,%d), (%d,%d)\n", I[0], I[1],
      //      PE[0][0], PE[0][1], PE[1][0], PE[1][1]);
      // do all dim agree?
      agree = xfe_True;
      for (d=0; d<dim-1; d++){
	if ((PE[d][0] != PE[d+1][0]) ||
	    (PE[d][1] != PE[d+1][1]))
	  agree = xfe_False;
      }
      // if all dim agree, add to list, increment all dim, continue
      if (agree){
	ESS->buf[2*nelem+0] = PE[d][0];
	ESS->buf[2*nelem+1] = PE[d][1];
	nelem++;
	for (d=0; d<dim; d++) I[d]++;
	continue;
      }
      // at this point, elems are different
      // determine max elem
      dmax = 0;
      for (d=1, dmax=0; d<dim; d++)
	if ((PE[d][0] > PE[dmax][0]) ||
	    ((PE[d][0] == PE[dmax][0]) && (PE[d][1] > PE[dmax][1])))
	  dmax = d;
      // increase index of all < max elem until pass max elem or end of list (done)
      for (d=0; d<dim; d++){
	if (d == dmax) continue;
	// "catch up" index d
	while ((PE[d][0] < PE[dmax][0]) ||
	       ((PE[d][0] == PE[dmax][0]) && (PE[d][1] < PE[dmax][1]))){
	  I[d]++;
	  if (I[d] >= ESS->nElem[d][vBin[d]]){
	    done = xfe_True;
	    break;
	  }
	  PE[d] = ESS->Elem[d][vBin[d]] + 2*I[d];
	}
	if (done) break;
      } // d
      if (done) break;
    } // while ! done
  }

    
  // some processor must have flagged some elements by now
  /*
  nelemmax = nelem;
  ierr = xf_Error(xf_MPI_Allreduce(&nelemmax, 1, xfe_SizeInt, xfe_MPI_MAX));
  if (ierr != xf_OK) return ierr;

  if (nelemmax <= 0){
    xf_printf("No common elements found in the dimension-based bin search.\n");
    return xf_NOT_FOUND;
  }*/

  // apply Glob2RefElem to all elements in ESS->buf
  egrp = -1;
  elem = -1;
  found = xfe_False;
  if (nelem > 0){
    for (ielem=0; ielem<nelem; ielem++){
      egrp = ESS->buf[2*ielem+0];
      elem = ESS->buf[2*ielem+1];
    
      //xf_printf("In search: looking at (%d,%d)\n", egrp,elem);
      // determine element Shape 
      ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
      if (ierr != xf_OK) return ierr;
      // get reference coordinates
      ierr = xf_Error(xf_Glob2RefElem(Mesh, egrp, elem, xpoint, tol, xfe_False,
				      xref+ip*dim, &converged));
      if (ierr != xf_OK) return ierr;
      //xf_printf(" converged = %d\n", converged);
      if (!converged) continue; // skip elements where glob2refelem fails
      // if inside, found=True, break
      ierr = xf_Error(xf_InsideShape(Shape, xref+ip*dim, tol, &inside));
      if (ierr != xf_OK) return ierr;
      //xf_printf(" inside = %d\n", inside);
      if (inside){
	found = xfe_True;
	break;
      }
    } // ielem
  }

  d = found;
    
  // reduce "found" among all processors
  infotank[ip] = (int) found;
    
    // set element to return
    if(d){
        pegrp[ip] = egrp;
        pelem[ip] = elem;
    }
    if(!d){  //not in this decomposition
        pegrp[ip] = -1;
        pelem[ip] = -1;
    }
}//np
    
  //ierr = xf_Error(xf_MPI_Allreduce((int *) &found, 1, xfe_SizeInt, xfe_MPI_MAX));
  //if (ierr != xf_OK) return ierr;
  xf_printf("Start transmitting point-search results\n");
  ierr = xf_Error(xf_MPI_Allreduce(infotank, np, xfe_SizeInt, xfe_MPI_MAX));
  if (ierr != xf_OK) return ierr;
  xf_printf("Finish transmitting point-search results\n");

  //another loop for search result check
  for(ip=0; ip<np; ip++){
        
    
  //if (!found){
  if (!infotank[ip]){
    xf_printf("Point not inside any of the elements after Glob2RefElem:\n");
    for (d=0; d<dim; d++) xf_printf("%.6E ", xpoint[d]);
    xf_printf("\n");
    return xf_NOT_FOUND;
  }
 
  // in case there are multiple matches (e.g. on an element boundary), use unique val
  // not using MIN because some procs have egrp=elem=-1
  //ierr = xf_Error(xf_MPI_Allreduce(&egrp, 1, xfe_SizeInt, xfe_MPI_MAX));
  //if (ierr != xf_OK) return ierr;
  //ierr = xf_Error(xf_MPI_Allreduce(&elem, 1, xfe_SizeInt, xfe_MPI_MAX));
  //if (ierr != xf_OK) return ierr;
  
  }//np
    
  xf_Release( (void *) infotank);
  
  return xf_OK;
}


