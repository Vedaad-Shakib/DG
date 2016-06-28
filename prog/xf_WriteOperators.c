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
  FILE:  xf_WriteOperators.c

  Writes out custom operators (e.g. restriction or prolongation) in text format

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_MeshTools.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Math.h"
#include "xf_Quad.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_Arg.h"
#include "xf_Solver.h"




/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr;
  int i1, i2;
  int n1, n2;
  int order1, order2;
  enum xfe_BasisType basis1, basis2;
  char *ArgIn[] = {"basis1", "TriLagrange", "first basis",
		   "basis2", "TriLagrange", "second basis",
                   "order1", "0", "first order",
		   "order2", "0", "second order",
                   "outfile", "NULL", "outputfile",
		   "\0"};
  char outfile[xf_MAXSTRLEN];
  real *TT, val;
  xf_Matrix *T = NULL;
  FILE *fid;
  xf_KeyValue KeyValue;

  
  xf_printf("\n");
  xf_printf("=== Operator text writing program  ===\n");
  xf_printf("\n");
    
      
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);
  if ((ierr != xf_OK) && (ierr != xf_FORCE_QUIT)) return xf_Error(ierr);
    
  // Get basis1
  ierr = xf_Error(xf_GetKeyValueEnum(KeyValue, "basis1", xfe_BasisName, 
                                     xfe_BasisLast, (int *) &basis1));
  if (ierr != xf_OK) return ierr;

  // Get basis2
  ierr = xf_Error(xf_GetKeyValueEnum(KeyValue, "basis1", xfe_BasisName, 
                                     xfe_BasisLast, (int *) &basis2));
  if (ierr != xf_OK) return ierr;

  // Get order 1
  ierr = xf_Error(xf_GetKeyValueInt(KeyValue, "order1", &order1));
  if (ierr != xf_OK) return ierr;

  // Get order 2
  ierr = xf_Error(xf_GetKeyValueInt(KeyValue, "order2", &order2));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_GetKeyValue(KeyValue, "outfile", outfile));
  if (ierr != xf_OK) return ierr;

  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;

  
  // are we riting to a file?
  if (xf_NotNull(outfile)){
    if ((fid = fopen(outfile, "w")) == NULL) return xf_Error(xf_FILE_WRITE_ERROR);
  }
  else fid = stdout;


  /***  prolongation matrix ***/

  ierr = xf_Error(xf_FindTransferMatrix(NULL, basis1, order1, basis2, order2, &T));
  if (ierr != xf_OK) return ierr;
  
  // pull off real array of the matrix
  TT = T->GenArray->rValue[0];

  // matrix dimensions
  ierr = xf_Error(xf_Order2nNode(basis1, order1, &n1));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Order2nNode(basis2, order2, &n2));
  if (ierr != xf_OK) return ierr;

  fprintf(fid,"%%Prolongation matrix FROM (order1,basis1) = (%d,%s) TO (order2,basis2) = (%d,%s):\n",
            order1, xfe_BasisName[basis1], order2, xfe_BasisName[basis2]);
  fprintf(fid,"%%(the matrix is %d by %d)\n", n2, n1);
  for (i2=0; i2<n2; i2++){
    for (i1=0; i1<n1; i1++){
      val = TT[i2*n1+i1];
      if (fabs(val) < 100.*MEPS) val = 0.;  // round to 0
      fprintf(fid,"%20.15e ", val);
    }
    fprintf(fid,"\n");
  }


  // T was created stand-alone and must be destroyed
  ierr = xf_Error(xf_DestroyMatrix(T));
  if (ierr != xf_OK) return ierr;

  xf_printf("xf_WriteOperators finished.\n");

  return xf_OK;
}
