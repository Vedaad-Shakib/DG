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
 FILE:  xf_MeshStats.c
 
 This file contains functions for for calculating mesh statistics
 
 */

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_Mesh.h"
#include "xf_MeshTools.h"
#include "xf_Math.h"
#include "xf_MPI.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_Data.h"
#include "xf_Arg.h"

/******************************************************************/
//   FUNCTION Definition: xf_CalcElemQuality
static int
xf_CalcElemQuality(xf_All *All, xf_Vector *ElemQuality, 
                   enum xfe_ElemStatType StatMeasure, enum xfe_Bool UseLog)
{
  int ierr, egrp, elem, Dim;
  real Surface, Volume;
  xf_Vector *EG;
  
  //find or calculate the geometric properties of the elements
  ierr = xf_Error(xf_FindElemGeom(All, &EG));
  if (ierr != xf_OK) return ierr;
  
  switch (StatMeasure) {
    case xfe_ESVolume:
      for (egrp = 0; egrp < All->Mesh->nElemGroup; egrp++) {
        for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++) {
          Volume = EG->GenArray[egrp].rValue[elem][xfe_EGVolume];
          if (UseLog)
            ElemQuality->GenArray[egrp].rValue[elem][0] = log10(Volume);
          else
            ElemQuality->GenArray[egrp].rValue[elem][0] = Volume;
        }
      }
      break;
    case xfe_ESSurfArea:
      for (egrp = 0; egrp < All->Mesh->nElemGroup; egrp++) {
        for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++) {
          Surface = EG->GenArray[egrp].rValue[elem][xfe_EGSurfArea];
          if (UseLog)
            ElemQuality->GenArray[egrp].rValue[elem][0] = log10(Surface);
          else 
            ElemQuality->GenArray[egrp].rValue[elem][0] = Surface;

        }
      }
      break;
    case xfe_ESAspectRatio:
      for (egrp = 0; egrp < All->Mesh->nElemGroup; egrp++) {
        for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++) {
          Volume = EG->GenArray[egrp].rValue[elem][xfe_EGVolume];
          Surface = EG->GenArray[egrp].rValue[elem][xfe_EGSurfArea];
          Dim = All->Mesh->Dim;
          if (UseLog)
            ElemQuality->GenArray[egrp].rValue[elem][0] = log10(pow(Surface/(2.0*Dim), 
                                                                    Dim/(Dim-1))/Volume);
          else
            ElemQuality->GenArray[egrp].rValue[elem][0] = pow(Surface/(2.0*Dim), 
                                                              Dim/(Dim-1))/Volume;
        }
      }
      break;
    default:
      return xf_NOT_SUPPORTED;
  }
  
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, terr, i, len, nbin, egrp, elem, MaxOrder, MinOrder, Order, nElem, j;
  real MaxQ, MinQ, dQ, rbuf[5], frac, *EQ, *EQunsorted, AbsMaxQ, AbsMinQ, Q;
  int **hist, bin, *pos, rank_top, rank_bottom;
  int myRank, nProc, nElemGlob, top, bottom, *nSelElemInEgrp, **SelElemInEgrp;
  char *ArgIn[] = {"in", "NULL", "input mesh file (with extension)",
    "frac", "1.0", "fraction of the number of elements to include in the histogram",
    "lowest","True", "if True, the histogram will include the \"frac\" elements with lowest \"quality\"",
    "nbin", "10", "number of bins in the histogram",
    "logspace", "False","if True, log spacing is used for the bins",
    "quality", "NULL", "type of quality measure (Surface, Volume, AspectRatio)",
    "writexfa", "False", "if True, a new xfa will be written with quality measure in it",
    "\0"};
  char inFile[xf_MAXSTRLEN], outFile[xf_MAXSTRLEN], quality[xf_MAXSTRLEN];  
  char *inExt;
  enum xfe_Bool IsThereData = xfe_False, VariableOrder = xfe_False;
  enum xfe_Bool lowest, logspacing, writexfa;
  enum xfe_ElemStatType StatMeasure;
  xf_KeyValue KeyValue;
  xf_All *All;
  xf_Data *D;
  xf_Mesh *Mesh, *Mesh_Glob;
  xf_Vector *U, *ElemQuality;
  FILE *fid;
  
  /* Initialize parallel-run (no effect in serial) */
  ierr = xf_Error(xf_MPI_Init(&argc, &argv));
  if (ierr != xf_OK) return ierr;
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  xf_printf("\n");
  xf_printf("=== xf_MeshStats: Mesh Statistics  ===\n");
  xf_printf("\n");
  
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;
  
  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
  
  xf_printf("nKey = %d\n", KeyValue.nKey);
  for (i=0; i<KeyValue.nKey; i++)
    xf_printf("%d : Key = %s, Value = %s\n", i, KeyValue.Key[i], 
              KeyValue.Value[i]);
  
  /* Get inputs */
  // Get inFile
  ierr = xf_GetKeyValue(KeyValue, "in", inFile);
  if (ierr != xf_OK) return ierr;
  
  // Get quality
  ierr = xf_GetKeyValueEnum(KeyValue, "quality", xfe_ElemStatName, 
                            xfe_ESLast, (int *)&StatMeasure);
  if (ierr != xf_OK) return ierr;
  
  // Get nbin
  ierr = xf_GetKeyValueInt(KeyValue, "nbin", &nbin);
  if (ierr != xf_OK) return ierr;
  
  // Get range
  ierr = xf_GetKeyValueReal(KeyValue, "frac", &frac);
  if (ierr != xf_OK) return ierr;
  if (frac > 1.0 || frac < 0.0) {
    xf_printf("frac is larger than 1.0 or lower than 0.0\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  
  ierr = xf_GetKeyValueBool(KeyValue, "lowest", &lowest);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_GetKeyValueBool(KeyValue, "logspace", &logspacing);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_GetKeyValueBool(KeyValue, "writexfa", &writexfa);
  if (ierr != xf_OK) return ierr;
  /* Finished getting inputs */
  
  // extensions are required
  if ((len = strlen(inFile)) < 4) return xf_Error(xf_FILE_READ_ERROR);
  inExt = inFile + len - 4; // pointer to extension
  
  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;
  
  //create the entire all structure just in case the input is xfa.
  ierr = xf_Error(xf_CreateAll(&All, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  // Read in the input file
  if (strncmp(inExt, ".gri", 4) == 0){
    Mesh = All->Mesh;
    // root creates Mesh_Glob and reads the ascii file
    if (myRank == 0){
      terr = xf_Error(xf_CreateMesh(&Mesh_Glob));
      if (terr == xf_OK)
        terr = xf_Error( xf_ReadGriFile(inFile, NULL, Mesh_Glob) );
    }
    ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    if (terr != xf_OK) return xf_Error(terr);
    
    // Mesh_Glob is parallelized
    ierr = xf_Error(xf_ParallelizeMesh(Mesh_Glob, Mesh, NULL, NULL));
    if (ierr != xf_OK) return ierr;
    
    // no longer need Mesh_Glob
    ierr = xf_Error(xf_DestroyMesh(Mesh_Glob));
    if (ierr != xf_OK) return ierr;
  }
  else if (strncmp(inExt, ".xfa", 4) == 0){
    ierr = xf_Error(xf_ReadAllBinary(inFile, All));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_FindPrimalState(All->DataSet, 0, &D, NULL));
    if (ierr != xf_OK) return ierr;
    U = (xf_Vector *)D->Data;
    IsThereData = xfe_True;
    if (U->vOrder == NULL)
      VariableOrder = xfe_False;
    else
      VariableOrder = xfe_True;
  }
  else {
    xf_printf("inFile extension: %s not supported.\n",inExt);
    return xf_Error(xf_NOT_SUPPORTED);
  }
  
  //get total number of elements
  ierr = xf_Error(xf_GetnElem(All->Mesh, NULL, &nElem));
  if (ierr != xf_OK) return ierr;
  //let all processors know the global total number of elements
  nElemGlob = nElem;
  ierr = xf_Error(xf_MPI_Allreduce(&nElemGlob, 1, xfe_SizeInt, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;
  
  //create a vector for storing stats
  if (logspacing)
    sprintf(quality,"log(%s)",xfe_ElemStatName[StatMeasure]);
  else
    sprintf(quality,"%s",xfe_ElemStatName[StatMeasure]);
  ierr = xf_Error(xf_FindVector(All, quality, xfe_LinkageGlobElem, 
                                1, NULL, 0, 0, NULL, NULL, NULL, 
                                NULL, NULL, xfe_SizeReal, xfe_False, 
                                writexfa, &D, &ElemQuality, NULL));
  if (ierr != xf_OK) return ierr;
  if (writexfa)//this is this way to avoid memory violation
    D->ReadWrite = xfe_True;
  
  //array for ranking elements
  ierr = xf_Error(xf_Alloc((void **)&pos, nElem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **)&EQ, nElem, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **)&EQunsorted, nElem, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_CalcElemQuality(All, ElemQuality, StatMeasure, logspacing));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **)&nSelElemInEgrp, All->Mesh->nElemGroup, 
                           sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  //backup data for ranking
  i=0;
  MaxOrder = -1;
  MinOrder = 100;//Maybe one day we will get to this approximation order...hehe
  Mesh = All->Mesh;
  Order = 0;
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    //for now, nSelElemInEgrp will store the number of elements in egrp
    nSelElemInEgrp[egrp] = Mesh->ElemGroup[egrp].nElem;
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
      //get max and min order
      if (IsThereData){
        Order = xf_InterpOrder(U, egrp, elem);
        if (Order < MinOrder)
          MinOrder = Order;
        if (Order > MaxOrder)
          MaxOrder = Order;
      }
      EQ[i] = ElemQuality->GenArray[egrp].rValue[elem][0];
      EQunsorted[i] = ElemQuality->GenArray[egrp].rValue[elem][0];
      i++;
    }
  }
  if (VariableOrder){
    ierr = xf_Error(xf_MPI_Allreduce(&MaxOrder, 1, xfe_SizeInt, xfe_MPI_MAX));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_MPI_Allreduce(&MinOrder, 1, xfe_SizeInt, xfe_MPI_MIN));
    if (ierr != xf_OK) return ierr;
  }
  else {
    MaxOrder = MinOrder = Order;
  }

  //sort data
  ierr = xf_Error(xf_SortRealParallel(EQ, nElem, xfe_False, pos));
  if (ierr != xf_OK) return ierr;
  
  //cutoff position
  if (lowest){
    rank_top = 0;
    rank_bottom = (int)(frac*nElemGlob)-1;
    if (rank_bottom < 0)
      return xf_Error(xf_INPUT_ERROR);
  }
  else {
    rank_top = (int)((1.0-frac)*nElemGlob);
    if (rank_top == nElemGlob)
      return xf_Error(xf_INPUT_ERROR);
    rank_bottom = nElemGlob-1;
  }
  
  /*loop through elements in each processor and find range,
  and create a lookup table for the elements*/
  MaxQ = MinQ = -1.0e300;
  top = bottom = -1;
  i = 0;
  ierr = xf_Error(xf_VAlloc2((void ***)&SelElemInEgrp, 
                             Mesh->nElemGroup, nSelElemInEgrp, 
                             sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++) {
    nSelElemInEgrp[egrp] = 0;
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++) {
      if (pos[i] == rank_top) {//found top of the list
        top = i;
        MinQ = EQunsorted[i];
      }
      if (pos[i] == rank_bottom) {//found bottom of the list
        bottom = i;
        MaxQ = EQunsorted[i];
      }
      //within range
      if (pos[i] >= rank_top && pos[i] <= rank_bottom) {
        SelElemInEgrp[egrp][nSelElemInEgrp[egrp]] = elem;
        nSelElemInEgrp[egrp]++;
      }
      i++;
    }
  }

  ierr = xf_Error(xf_MPI_Allreduce(&top, 1, xfe_SizeInt, xfe_MPI_MAX));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_MPI_Allreduce(&bottom, 1, xfe_SizeInt, xfe_MPI_MAX));
  if (ierr != xf_OK) return ierr;
  
  
  rbuf[0] = MinQ;
  rbuf[1] = MaxQ;
  rbuf[2] = EQ[nElem-1];
  
  ierr = xf_Error(xf_MPI_Allreduce(rbuf, 3, xfe_SizeReal, xfe_MPI_MAX));
  if (ierr != xf_OK) return ierr;
  
  MinQ = rbuf[0];
  MaxQ = rbuf[1];
  AbsMaxQ = rbuf[2];
  
  AbsMinQ = EQ[0];
  ierr = xf_Error(xf_MPI_Allreduce(&AbsMinQ, 1, xfe_SizeReal, xfe_MPI_MIN));
  if (ierr != xf_OK) return ierr;
  
  //No need to reduce the range since it's the same on all processors
  dQ = (MaxQ-MinQ)/nbin;
  //we are ready to fill the bins.
  ierr = xf_Error(xf_Alloc2((void ***)&hist, nbin, 
                            (MaxOrder-MinOrder)+1, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  for (i = 0; i < nbin; i++)
    for (j = 0; j < MaxOrder-MinOrder+1; j++)
      hist[i][j] = 0;
  
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    for (i = 0; i < nSelElemInEgrp[egrp]; i++){
      elem = SelElemInEgrp[egrp][i];
      if (IsThereData)
        Order = xf_InterpOrder(U, egrp, elem);
      else 
        Order = MinOrder;
      //histogram
      Q = ElemQuality->GenArray[egrp].rValue[elem][0];
      bin = (int)((Q-MinQ)/dQ);
      if (bin == nbin)
        bin--;
      hist[bin][Order-MinOrder]++;
    }
  }
  //reduce
  ierr = xf_Error(xf_MPI_Allreduce(hist[0], nbin*((MaxOrder-MinOrder)+1), 
                                   xfe_SizeInt, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;

  //write out in matlab format
  //file name
  if (myRank == 0){//only root does the writing
    //file name
    inFile[strlen(inFile)-4] = '\0';
    sprintf(outFile,"%s_%s.m",inFile,xfe_ElemStatName[StatMeasure]);
    
    if ((fid = fopen(outFile, "w")) == NULL)
      return xf_Error(xf_FILE_WRITE_ERROR);
    fprintf(fid, "nbin = %d;\n",nbin);
    if (logspacing){
      fprintf(fid, "AbsMin%s = %1.10E;\n",xfe_ElemStatName[StatMeasure],pow(10.0,AbsMinQ));
      fprintf(fid, "AbsMax%s = %1.10E;\n",xfe_ElemStatName[StatMeasure],pow(10.0,AbsMaxQ));
      fprintf(fid, "Min%s = %1.10E;\n",xfe_ElemStatName[StatMeasure],pow(10.0,MinQ));
      fprintf(fid, "Max%s = %1.10E;\n",xfe_ElemStatName[StatMeasure],pow(10.0,MaxQ));
    }
    else {
      fprintf(fid, "AbsMin%s = %1.10E;\n",xfe_ElemStatName[StatMeasure],AbsMinQ);
      fprintf(fid, "AbsMax%s = %1.10E;\n",xfe_ElemStatName[StatMeasure],AbsMaxQ);
      fprintf(fid, "Min%s = %1.10E;\n",xfe_ElemStatName[StatMeasure],MinQ);
      fprintf(fid, "Max%s = %1.10E;\n",xfe_ElemStatName[StatMeasure],MaxQ);
    }

            
    fprintf(fid, "hist_%s=[",xfe_ElemStatName[StatMeasure]);
    for (i = 0; i < nbin; i++){
      for (j = 0; j < MaxOrder-MinOrder+1; j++){
        fprintf(fid, "%d ",hist[i][j]);
      }
      fprintf(fid, "%% bin %d\n",i+1);
    }
    fprintf(fid, "];\n");
        
    fprintf(fid, "cent_%s=[",xfe_ElemStatName[StatMeasure]);
    for (i = 0; i < nbin; i++){
      if (logspacing)
        fprintf(fid, "%1.12E\n",pow(10.0,MinQ+((double)i+0.5)*dQ));
      else
        fprintf(fid, "%1.12E\n",MinQ+((double)i+0.5)*dQ);
    }
    fprintf(fid, "];\n");
    
    fprintf(fid, "figure(1)\n");
    fprintf(fid, "bar(cent_%s,hist_%s,'stacked')\n",
            xfe_ElemStatName[StatMeasure],xfe_ElemStatName[StatMeasure]);
    if (logspacing)
      fprintf(fid, "xlabel('%s - log_{10} spaced')\n",xfe_ElemStatName[StatMeasure]);
    else
      fprintf(fid, "xlabel('%s')\n",xfe_ElemStatName[StatMeasure]);
    fprintf(fid, "ylabel('Occurrences')\n");
    fclose(fid);
  }
  xf_MPI_Barrier();

  //clean up
  xf_Release((void *)EQ);
  xf_Release((void *)EQunsorted);
  xf_Release((void *)pos);
  xf_Release((void *)nSelElemInEgrp);
  xf_Release2((void **)SelElemInEgrp);
  xf_Release2((void **)hist);
  if (writexfa){
    sprintf(outFile,"%s_%s.xfa",inFile,xfe_ElemStatName[StatMeasure]);
    ierr = xf_Error(xf_WriteAllBinary(All, outFile));
    if (ierr != xf_OK) return ierr;
  }
  else
    xf_DestroyVector(ElemQuality, xfe_True);
  
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;

  xf_printf("xf_MeshStats finished.\n");
  
  /* MPI finalize */
  ierr = xf_Error(xf_MPI_Finalize());
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}



