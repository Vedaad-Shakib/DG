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
 FILE:  xf_EqnSet.c
 
 This file contains functions for working with the EqnSet structure.
 
 */


#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Param.h"
#include "xf_EqnSetHook.h"
#include "xf_Output.h"


/******************************************************************/
//   FUNCTION Definition: xf_CreateResTerms
static int 
xf_CreateResTerms( xf_ResTerms **pResTerms){
  
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pResTerms, 1, sizeof(xf_ResTerms)));
  if (ierr != xf_OK) return ierr;
  
  (*pResTerms)->nResTerm = 0;
  (*pResTerms)->ResTerm  = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_AllocResTerms
static int 
xf_AllocResTerms( xf_ResTerms *ResTerms, int nResTerm){
  
  int ierr, i;
  
  ierr = xf_Error(xf_Alloc((void **) &ResTerms->ResTerm, nResTerm, sizeof(xf_ResTerm)));
  if (ierr != xf_OK) return ierr;
  
  ResTerms->nResTerm = nResTerm;
  
  for (i=0; i<nResTerm; i++){
    ResTerms->ResTerm[i].Type = xfe_ResTermUnknown;
    
    /* Initialize key-value structure */
    ierr = xf_Error(xf_InitKeyValue(&ResTerms->ResTerm[i].KeyValue));
    if (ierr != xf_OK) return ierr;
    
    ResTerms->ResTerm[i].Active = xfe_True;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyResTerm
int 
xf_DestroyResTerm( xf_ResTerm *ResTerm){
  
  int ierr;
  
  /* Destroy key-value structure */
  ierr = xf_Error(xf_DestroyKeyValue(&ResTerm->KeyValue));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyResTerms
int 
xf_DestroyResTerms( xf_ResTerms *ResTerms){
  
  int ierr, i;
  
  if (ResTerms == NULL) return xf_OK;
  
  for (i=0; i<ResTerms->nResTerm; i++){
    ierr = xf_Error(xf_DestroyResTerm(ResTerms->ResTerm+i));
    if (ierr != xf_OK) return ierr;
  }
  xf_Release((void *) ResTerms->ResTerm);
  
  xf_Release((void *) ResTerms);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CopyResTerm
int 
xf_CopyResTerm( xf_ResTerm *ResTerm1, xf_ResTerm *ResTerm2)
{
  int ierr;
  
  ResTerm2->Type = ResTerm1->Type;
  
  /* Copy key-value structure */
  ierr = xf_Error(xf_CopyKeyValue(&ResTerm1->KeyValue, &ResTerm2->KeyValue));
  if (ierr != xf_OK) return ierr;
  
  ResTerm2->Active = ResTerm1->Active;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_AddResTerm
int 
xf_AddResTerm( xf_ResTerms *ResTerms, enum xfe_ResTermType Type, char **KeyValueList)
{
  
  int ierr;
  int nResTerm;
  
  nResTerm = ResTerms->nResTerm++;
  ierr = xf_Error(xf_ReAlloc((void **) &ResTerms->ResTerm, nResTerm+1, sizeof(xf_ResTerm)));
  if (ierr != xf_OK) return ierr;
  
  ResTerms->ResTerm[nResTerm].Type = Type;
  
  /* Initialize key-value structure */
  ierr = xf_Error(xf_InitKeyValue(&ResTerms->ResTerm[nResTerm].KeyValue));
  if (ierr != xf_OK) return ierr;
  
  /* Add list to key-values */
  ierr = xf_Error(xf_AddKeyValueList(&ResTerms->ResTerm[nResTerm].KeyValue, 
                                     KeyValueList, xfe_True, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ResTerms->ResTerm[nResTerm].Active = xfe_True;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateICs
static int 
xf_CreateICs( xf_ICs **pICs){
  
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pICs, 1, sizeof(xf_ICs)));
  if (ierr != xf_OK) return ierr;
  
  (*pICs)->nIC = 0;
  (*pICs)->IC  = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_AllocICs
static int 
xf_AllocICs( xf_ICs *ICs, int nIC){
  
  int ierr, i;
  
  ierr = xf_Error(xf_Alloc((void **) &ICs->IC, nIC, sizeof(xf_IC)));
  if (ierr != xf_OK) return ierr;
  
  ICs->nIC = nIC;
  
  for (i=0; i<nIC; i++){
    ICs->IC[i].Type     = NULL;
    ICs->IC[i].Function = NULL;
    ICs->IC[i].Header   = NULL;
    ICs->IC[i].Data     = NULL;
    ICs->IC[i].AlterFunction = NULL;
    ICs->IC[i].AlterData     = NULL;
    ICs->IC[i].PriorSteadySolve = xfe_False;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyIC
static int 
xf_DestroyIC( xf_IC *IC){
  
  xf_Release(IC->Type);
  xf_Release(IC->Function);
  xf_Release(IC->Header);
  xf_Release(IC->Data);
  xf_Release(IC->AlterFunction);
  xf_Release(IC->AlterData);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyICs
static int 
xf_DestroyICs( xf_ICs *ICs){
  
  int ierr, i;
  
  if (ICs == NULL) return xf_OK;
  
  for (i=0; i<ICs->nIC; i++){
    ierr = xf_Error(xf_DestroyIC(ICs->IC+i));
    if (ierr != xf_OK) return ierr;
  }
  xf_Release((void *) ICs->IC);
  
  xf_Release((void *) ICs);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateBCs
static int 
xf_CreateBCs( xf_BCs **pBCs){
  
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pBCs, 1, sizeof(xf_BCs)));
  if (ierr != xf_OK) return ierr;
  
  (*pBCs)->nBC = 0;
  (*pBCs)->BC  = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_AllocBCs
static int 
xf_AllocBCs( xf_BCs *BCs, int nBC){
  
  int ierr, i;
  
  ierr = xf_Error(xf_Alloc((void **) &BCs->BC, nBC, sizeof(xf_BC)));
  if (ierr != xf_OK) return ierr;
  
  BCs->nBC = nBC;
  
  for (i=0; i<nBC; i++){
    BCs->BC[i].BFGTitle      = NULL;
    BCs->BC[i].Type          = NULL;
    BCs->BC[i].Function      = NULL;
    BCs->BC[i].Header        = NULL;
    BCs->BC[i].Data          = NULL;
    BCs->BC[i].OutputLinkage = NULL;
    BCs->BC[i].nBCParam      = 0;
    BCs->BC[i].BCParam       = NULL;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyBC
static int 
xf_DestroyBC( xf_BC *BC){
  
  xf_Release(BC->BFGTitle);
  xf_Release(BC->Type);
  xf_Release(BC->Function);
  xf_Release(BC->Header);
  xf_Release(BC->Data);
  xf_Release(BC->OutputLinkage);
  xf_Release(BC->BCParam);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyBCs
static int 
xf_DestroyBCs( xf_BCs *BCs){
  
  int ierr, i;
  
  if (BCs == NULL) return xf_OK;
  
  for (i=0; i<BCs->nBC; i++){
    ierr = xf_Error(xf_DestroyBC(BCs->BC+i));
    if (ierr != xf_OK) return ierr;
  }
  xf_Release((void *) BCs->BC);
  
  xf_Release((void *) BCs);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CreateOutputs
static int 
xf_CreateOutputs( xf_Outputs **pOutputs)
{
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pOutputs, 1, sizeof(xf_Outputs)));
  if (ierr != xf_OK) return ierr;
  
  (*pOutputs)->nOutput = 0;
  (*pOutputs)->Output  = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_InitOutput
static void 
xf_InitOutput( xf_Output *Output)
{
  int j;
  
  Output->Name                 = NULL;
  Output->Type                 = xfe_DomainIntegral;
  Output->DomainNorm           = xfe_DomainNormNone;
  Output->TimeNorm             = xfe_TimeNormNone;
  Output->UsesFlux             = xfe_False;
  Output->ScalarName           = NULL;
  Output->VectorName           = NULL;
  Output->Function             = NULL;
  Output->Data                 = NULL;
  Output->nFluxComponent       = 0;
  Output->FluxComponentNames   = NULL;
  Output->FluxComponentWeights = NULL;
  Output->FluxComponentMoments = NULL;
  for (j=0; j<6; j++) Output->LineCoord[j] = 0.;
  for (j=0; j<4; j++) Output->CutPlane[j] = 0.;
  Output->StartTime            = 0.;
  Output->EndTime              = 0.;
  Output->Value                = 0.;
  Output->ErrEst               = 0.;
  Output->nBFG                 = 0;
  Output->BFGTitles            = NULL;
  Output->DumpFile             = NULL;
  Output->CutPlaneIntersect    = NULL;
  Output->egrp                 = 0;
  Output->elem                 = 0;
  Output->elemLocal            = NULL;
  for (j=0; j<3; j++) Output->xref[j] = 0.;
  Output->nSumOutput           = 0;
  Output->SumOutputNames       = NULL;
  Output->SumOutputWeights     = NULL;
  Output->SumOutputErrTols     = NULL;
  Output->nSensitivity         = 0;
  Output->Sensitivity          = NULL;
}

/******************************************************************/
//   FUNCTION Definition: xf_AllocOutputs
static int 
xf_AllocOutputs( xf_Outputs *Outputs, int nOutput)
{
  int ierr, i;
  
  ierr = xf_Error(xf_Alloc((void **) &Outputs->Output, nOutput, sizeof(xf_Output)));
  if (ierr != xf_OK) return ierr;
  
  Outputs->nOutput = nOutput;
  
  for (i=0; i<nOutput; i++) xf_InitOutput(Outputs->Output+i);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroySensitivity
static void 
xf_DestroySensitivity( xf_Sensitivity *Sensitivity)
{
  xf_Release( (void *) Sensitivity->ParamName);
  xf_Release( (void *) Sensitivity->value);
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyOutput
static int 
xf_DestroyOutput( xf_Output *Output)
{
  int ierr, i;
  
  xf_Release( (void * ) Output->Name);
  xf_Release( (void * ) Output->ScalarName);
  xf_Release( (void * ) Output->VectorName);
  xf_Release( (void * ) Output->Function);
  xf_Release( (void * ) Output->Data);
  xf_Release2((void **) Output->FluxComponentNames);
  xf_Release( (void * ) Output->FluxComponentWeights);
  xf_Release( (void * ) Output->FluxComponentMoments);
  xf_Release2((void **) Output->BFGTitles);
  xf_Release2((void **) Output->DumpFile);
  xf_Release2((void **) Output->SumOutputNames);
  xf_Release( (void * ) Output->SumOutputWeights);
  xf_Release( (void * ) Output->SumOutputErrTols);
  for (i = 0; i < Output->nSensitivity; i++)
    xf_DestroySensitivity(Output->Sensitivity+i);
  xf_Release( (void * ) Output->Sensitivity);
  
  ierr = xf_Error(xf_DestroyCutPlaneIntersect(Output->CutPlaneIntersect));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void * ) Output->elemLocal);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReAllocOutputs
int 
xf_ReAllocOutputs( xf_Outputs *Outputs, int nOutput)
{
  int ierr, i;
  
  // Destroy outputs if reallocating to fewer outputs
  for (i=nOutput; i<Outputs->nOutput; i++){
    ierr = xf_Error(xf_DestroyOutput(Outputs->Output+i));
    if (ierr != xf_OK) return ierr;
  }
  
  // reallocate Outputs
  ierr = xf_Error(xf_ReAlloc((void **) &Outputs->Output, nOutput, sizeof(xf_Output)));
  if (ierr != xf_OK) return ierr;
  
  // Initialize new outputs if reallocating to more outputs
  for (i=Outputs->nOutput; i<nOutput; i++) xf_InitOutput(Outputs->Output+i);
  
  Outputs->nOutput = nOutput;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyOutputs
static int 
xf_DestroyOutputs( xf_Outputs *Outputs)
{
  int ierr, i;
  
  if (Outputs == NULL) return xf_OK;
  
  for (i=0; i<Outputs->nOutput; i++){
    ierr = xf_Error(xf_DestroyOutput(Outputs->Output+i));
    if (ierr != xf_OK) return ierr;
  }
  xf_Release((void *) Outputs->Output);
  
  xf_Release((void *) Outputs);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_InitEqnSet
static int 
xf_InitEqnSet( xf_EqnSet *EqnSet){
  
  int ierr;
  
  ierr = xf_Error(xf_InitKeyValue( &(EqnSet->KeyValue) ));
  if (ierr != xf_OK) return ierr;
  
  EqnSet->EqnSetLibrary = NULL;
  EqnSet->Dim           = 0;
  EqnSet->ResTerms      = NULL;
  EqnSet->ICs           = NULL;
  EqnSet->BCs           = NULL;
  EqnSet->Outputs       = NULL;
  EqnSet->StateRank     = 0; 
  EqnSet->StateName     = NULL;
  EqnSet->nPosInState   = 0; 
  EqnSet->PosInState    = NULL;
  EqnSet->nIParam       = 0;
  EqnSet->IParamKey     = NULL;
  EqnSet->nRParam       = 0;
  EqnSet->RParamKey     = NULL;
  EqnSet->nAuxU         = 0;
  EqnSet->AuxUNames     = NULL;
  EqnSet->PosInAuxU     = NULL;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateEqnSet
int 
xf_CreateEqnSet( xf_EqnSet **pEqnSet){
  
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pEqnSet, 1, sizeof(xf_EqnSet)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_InitEqnSet((*pEqnSet)));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_DestroyEqnSet
int 
xf_DestroyEqnSet( xf_EqnSet *EqnSet, enum xfe_Bool DestroySelf){
  
  int i, ierr;
  
  if (EqnSet == NULL) return xf_OK;
  
  xf_Release((void *)EqnSet->EqnSetLibrary);
  
  ierr = xf_Error(xf_DestroyKeyValue(&EqnSet->KeyValue));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyResTerms(EqnSet->ResTerms));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyICs(EqnSet->ICs));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyBCs(EqnSet->BCs));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyOutputs(EqnSet->Outputs));
  if (ierr != xf_OK) return ierr;
  
  xf_Release2((void **)EqnSet->StateName);
  
  xf_Release((void *)EqnSet->PosInState);
  
  xf_Release2((void **)EqnSet->IParamKey);
  xf_Release2((void **)EqnSet->RParamKey);
  
  xf_Release2((void **)EqnSet->AuxUNames);
  xf_Release((void *)EqnSet->PosInAuxU);
  
  if (DestroySelf) xf_Release((void *)EqnSet);
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_ReadEqnSetParam
static int 
xf_ReadEqnSetParam(FILE *feqn, char **InputStrings, int *piString,
                   xf_KeyValue *KeyValueParam)
{
  int ierr;
  char line0[200], *line;
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  enum xfe_Bool ReachedEnd = xfe_False;
  
  do{      
    /* read line of file */
    ierr = xf_LineFromFileOrStrings(feqn, InputStrings, piString, line0, &line);
    if (ierr == xf_END_OF_FILE) break;
    if (ierr != xf_OK) return ierr;
    
    if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line
    
    /* Check if reached end of block */
    if (strncmp(line, "ENDBLOCK", 8) == 0){
      ReachedEnd = xfe_True; break;
    }
    if (strncmp(line, "STARTBLOCK", 10) == 0) break;
    
    /* read key and value from line: every param must be of key=value form */
    ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
    if (ierr != xf_OK) return ierr;
    
    /* Store key and value in local list */
    ierr = xf_Error(xf_AddKeyValue(KeyValueParam, key, value, xfe_True));
    if (ierr == xf_OVERWROTE){
      xf_printf("Error. The key %s is assigned more than once in the EqnSet file.\n", key);
      return xf_FILE_READ_ERROR;
    }
    if (ierr != xf_OK) return ierr;
  } while (1);
  
  if (!ReachedEnd){
    xf_printf("Error. Param block not terminated with ENDBLOCK\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadEqnSetResTerm
static int 
xf_ReadEqnSetResTerms(FILE *feqn, char **InputStrings, int *piString,
                      xf_ResTerms **pResTerms)
{
  int ierr, nResTerm, iResTerm;
  char line0[200], *line;
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  enum xfe_Bool ReachedEnd = xfe_False;
  xf_ResTerms *ResTerms;
  
  ierr = xf_Error(xf_CreateResTerms(pResTerms));
  if (ierr != xf_OK) return ierr;
  
  ResTerms = (*pResTerms);
  
  nResTerm =  0;
  iResTerm = -1;
  do{      
    /* read line of file */
    ierr = xf_LineFromFileOrStrings(feqn, InputStrings, piString, line0, &line);
    if (ierr == xf_END_OF_FILE) break;
    if (ierr != xf_OK) return ierr;
    
    if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line
    
    /* Check if reached end of block */
    if (strncmp(line, "ENDBLOCK", 8) == 0){
      ReachedEnd = xfe_True; break;
    }
    if (strncmp(line, "STARTBLOCK", 10) == 0) break;
    
    if (strncmp(line, "nResTerm", 8) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      if (sscanf(value, "%d", &nResTerm) != 1) return xf_Error(xf_FILE_READ_ERROR);
      if (iResTerm != -1){
        xf_printf("Lines in ResTerm block out of order.\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      ierr = xf_Error(xf_AllocResTerms(ResTerms, nResTerm));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "TermType", 8) == 0){
      iResTerm++;  
      if (iResTerm >= nResTerm){
        xf_printf("Number of residual terms exceeds nResTerm.\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_Value2Enum(value, xfe_ResTermName, xfe_ResTermLast, 
                                    (int *) &ResTerms->ResTerm[iResTerm].Type));
      if( ierr != xf_OK ) return ierr;
      
      continue;
    }
    
    if (iResTerm < 0){
      xf_printf("Unrecognized line in residual term block.\n");
      return xf_Error(xf_FILE_READ_ERROR);
    }
    
    /* At this point, reading a key=value pair corresponding to iResTerm */
    ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
    if (ierr != xf_OK) return ierr;
    
    /* Store key and value in local list */
    ierr = xf_Error(xf_AddKeyValue(&ResTerms->ResTerm[iResTerm].KeyValue, 
                                   key, value, xfe_True));
    if (ierr == xf_OVERWROTE){
      xf_printf("Error. The key %s is assigned more than once in Residual term %d.\n", 
                key, iResTerm);
      return xf_FILE_READ_ERROR;
    }
    if (ierr != xf_OK) return ierr;
    
  } while (1);
  
  if (iResTerm != (nResTerm-1)){
    xf_printf("Number of residual terms is less than nResTerm.\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  if (!ReachedEnd){
    xf_printf("Error. ResTerm block not terminated with ENDBLOCK\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadEqnSetIC
static int 
xf_ReadEqnSetICs(FILE *feqn, char **InputStrings, int *piString, 
                 xf_ICs **pICs)
{
  
  int ierr, nIC, iIC;
  char line0[200], *line;
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  enum xfe_Bool ReachedEnd = xfe_False;
  xf_ICs *ICs;
  
  ierr = xf_Error(xf_CreateICs(pICs));
  if (ierr != xf_OK) return ierr;
  
  ICs = (*pICs);
  
  nIC =  0;
  iIC = -1;
  do{      
    /* read line of file */
    ierr = xf_LineFromFileOrStrings(feqn, InputStrings, piString, line0, &line);
    if (ierr == xf_END_OF_FILE) break;
    if (ierr != xf_OK) return ierr;
    
    if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line
    
    /* Check if reached end of block */
    if (strncmp(line, "ENDBLOCK", 8) == 0){
      ReachedEnd = xfe_True; break;
    }
    if (strncmp(line, "STARTBLOCK", 10) == 0) break;
    
    if (strncmp(line, "ICType", 6) == 0){
      
      if (nIC != 0){
        xf_printf("Please specify only one initial condition.\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      nIC = 1;
      ierr = xf_Error(xf_AllocICs(ICs, nIC));
      if (ierr != xf_OK) return ierr;
      
      iIC = 0;
      
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&ICs->IC[iIC].Type, xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (nIC != 1){
      xf_printf("Please specify ICType in eqnset initial conditions.\n");
      return xf_Error(xf_FILE_READ_ERROR);
    }
    
    if (strncmp(line, "Function", 8) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&ICs->IC[iIC].Function, xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "Header", 6) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&ICs->IC[iIC].Header, xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "Data", 4) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&ICs->IC[iIC].Data, xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "AlterFunction", 13) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&ICs->IC[iIC].AlterFunction, xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "AlterData", 9) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&ICs->IC[iIC].AlterData, xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "PriorSteadySolve", 16) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_Value2Enum(value, xfe_BoolName, xfe_BoolLast, 
                                    (int *) &ICs->IC[iIC].PriorSteadySolve));
      if( ierr != xf_OK ) return ierr;
      
      continue;
    }
    
    xf_printf("Unrecognized line in IC block.\n");
    return xf_Error(xf_FILE_READ_ERROR);
    
  } while (1);
  
  if (nIC != 1){
    xf_printf("Warning, number of initial conditions != 1.\n");
  }
  
  if (!ReachedEnd){
    xf_printf("Error. IC block not terminated with ENDBLOCK\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadEqnSetBC
static int 
xf_ReadEqnSetBCs(FILE *feqn, char **InputStrings, int *piString,
                 xf_BCs **pBCs)
{
  
  int ierr, nBC, iBC;
  char line0[200], *line;
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  enum xfe_Bool ReachedEnd = xfe_False;
  xf_BCs *BCs;
  
  ierr = xf_Error(xf_CreateBCs(pBCs));
  if (ierr != xf_OK) return ierr;
  
  BCs = (*pBCs);
  
  nBC =  0;
  iBC = -1;
  do{      
    /* read line of file */
    ierr = xf_LineFromFileOrStrings(feqn, InputStrings, piString, line0, &line);
    if (ierr == xf_END_OF_FILE) break;
    if (ierr != xf_OK) return ierr;
    
    if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line
    
    /* Check if reached end of block */
    if (strncmp(line, "ENDBLOCK", 8) == 0){
      ReachedEnd = xfe_True; break;
    }
    if (strncmp(line, "STARTBLOCK", 10) == 0) break;
    
    if (strncmp(line, "nBC", 2) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      if (sscanf(value, "%d", &nBC) != 1) return xf_Error(xf_FILE_READ_ERROR);
      if (iBC != -1){
        xf_printf("Lines in BC block out of order.\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      ierr = xf_Error(xf_AllocBCs(BCs, nBC));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "BFGTitle", 8) == 0){
      iBC++;  
      if (iBC >= nBC){
        xf_printf("Number of BFGTitles exceeds nBC.\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&BCs->BC[iBC].BFGTitle, xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (iBC < 0){
      xf_printf("Unrecognized line in BC block.\n");
      return xf_Error(xf_FILE_READ_ERROR);
    }
    
    if (strncmp(line, "BCType", 6) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&BCs->BC[iBC].Type, xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "Function", 8) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&BCs->BC[iBC].Function, xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "Header", 6) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&BCs->BC[iBC].Header, xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "Data", 4) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&BCs->BC[iBC].Data, xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "OutputLinkage", 13) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&BCs->BC[iBC].OutputLinkage, xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    xf_printf("Unrecognized line in BC block.\n");
    return xf_Error(xf_FILE_READ_ERROR);
    
  } while (1);
  
  if (iBC != (nBC-1)){
    xf_printf("Number of BC terms is less than nBC.\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  if (!ReachedEnd){
    xf_printf("Error. BC block not terminated with ENDBLOCK\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadEqnSetOutput
static int 
xf_ReadEqnSetOutputs(FILE *feqn, char **InputStrings, int *piString,
                     xf_Outputs **pOutputs)
{
  int ierr, nOutput, iOutput, i;
  int nFluxComponent, nSumOutput;
  char line0[xf_MAXLINELEN], *line, longvalue[xf_MAXLINELEN];
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  enum xfe_Bool ReachedEnd = xfe_False;
  xf_Outputs *Outputs;
  
  ierr = xf_Error(xf_CreateOutputs(pOutputs));
  if (ierr != xf_OK) return ierr;
  
  Outputs = (*pOutputs);
  
  nOutput =  0;
  iOutput = -1;
  nFluxComponent = -1;
  nSumOutput = -1;
  do{      
    /* read line of file */
    ierr = xf_LineFromFileOrStrings(feqn, InputStrings, piString, line0, &line);
    if (ierr == xf_END_OF_FILE) break;
    if (ierr != xf_OK) return ierr;
    
    if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line
    
    /* Check if reached end of block */
    if (strncmp(line, "ENDBLOCK", 8) == 0){
      ReachedEnd = xfe_True; break;
    }
    if (strncmp(line, "STARTBLOCK", 10) == 0) break;
    
    if (strncmp(line, "nOutput", 7) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      if (sscanf(value, "%d", &nOutput) != 1) return xf_Error(xf_FILE_READ_ERROR);
      if (iOutput != -1){
        xf_printf("Lines in Output block out of order.\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      ierr = xf_Error(xf_AllocOutputs(Outputs, nOutput));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (nOutput < 0){
      xf_printf("Error reading Output block.\n");
      return xf_Error(xf_FILE_READ_ERROR);
    }
    
    if (strncmp(line, "OutputName", 10) == 0){
      iOutput++;  
      if (iOutput >= nOutput){
        xf_printf("Number of OutputNames exceeds nOutput.\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&Outputs->Output[iOutput].Name, xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (iOutput < 0){
      xf_printf("Unrecognized line in Output block before first OutputName.\n");
      return xf_Error(xf_FILE_READ_ERROR);
    }
    
    
    if (strncmp(line, "OutputType", 10) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_Value2Enum(value, xfe_OutputName, xfe_OutputLast, 
                                    (int *) &Outputs->Output[iOutput].Type));
      if( ierr != xf_OK ) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "DomainNorm", 10) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_Value2Enum(value, xfe_DomainNormName, xfe_DomainNormLast, 
                                    (int *) &Outputs->Output[iOutput].DomainNorm));
      if( ierr != xf_OK ) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "TimeNorm", 8) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_Value2Enum(value, xfe_TimeNormName, xfe_TimeNormLast, 
                                    (int *) &Outputs->Output[iOutput].TimeNorm));
      if( ierr != xf_OK ) return ierr;
      
      continue;
    }
    
    
    if (strncmp(line, "UsesFlux", 7) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_Value2Enum(value, xfe_BoolName, xfe_BoolLast, 
                                    (int *) &Outputs->Output[iOutput].UsesFlux));
      if( ierr != xf_OK ) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "ScalarName", 10) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&Outputs->Output[iOutput].ScalarName, 
                                     xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "VectorName", 10) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&Outputs->Output[iOutput].VectorName, 
                                     xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "Function", 8) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&Outputs->Output[iOutput].Function, 
                                     xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    if (strncmp(line, "Data", 4) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&Outputs->Output[iOutput].Data, 
                                     xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    
    if (strncmp(line, "FluxComponentNames", 18) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, longvalue, xf_MAXLINELEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ScanXStringAlloc(longvalue, xf_MAXSTRLEN, &nFluxComponent,
                                          &Outputs->Output[iOutput].FluxComponentNames));
      if (ierr != xf_OK) return ierr;
      
      if ((Outputs->Output[iOutput].nFluxComponent > 0) && 
          (Outputs->Output[iOutput].nFluxComponent != nFluxComponent)){
        xf_printf("Mismatch in output # flux components: # weights != # names.\n");
        xf_printf(" [OutputName = %s]\n", Outputs->Output[iOutput].Name);
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      Outputs->Output[iOutput].nFluxComponent = nFluxComponent;
      continue;
    }
    
    if (strncmp(line, "FluxComponentWeights", 20) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, longvalue, xf_MAXLINELEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ScanXRealAlloc(longvalue, &nFluxComponent,
                                        &Outputs->Output[iOutput].FluxComponentWeights));
      if (ierr != xf_OK) return ierr;
      
      if ((Outputs->Output[iOutput].nFluxComponent > 0) && 
          (Outputs->Output[iOutput].nFluxComponent != nFluxComponent)){
        xf_printf("Mismatch in output # flux components: # weights != # names.\n");
        xf_printf(" [OutputName = %s]\n", Outputs->Output[iOutput].Name);
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      Outputs->Output[iOutput].nFluxComponent = nFluxComponent;
      continue;
    }
    
    if (strncmp(line, "FluxComponentMoments", 20) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, longvalue, xf_MAXLINELEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ScanXIntAlloc(longvalue, &nFluxComponent,
                                       &Outputs->Output[iOutput].FluxComponentMoments));
      if (ierr != xf_OK) return ierr;
      
      if ((Outputs->Output[iOutput].nFluxComponent > 0) && 
          (Outputs->Output[iOutput].nFluxComponent != nFluxComponent)){
        xf_printf("Mismatch in output # flux components: # moments != # names.\n");
        xf_printf(" [OutputName = %s]\n", Outputs->Output[iOutput].Name);
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      Outputs->Output[iOutput].nFluxComponent = nFluxComponent;
      continue;
    }
    
    if (strncmp(line, "LineCoord", 9) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ScanXReal(value, &i, Outputs->Output[iOutput].LineCoord));
      if (ierr != xf_OK) return ierr;
      
      if (i > 6) return xf_Error(xf_INPUT_ERROR);
      
      continue;
    }
    
    if (strncmp(line, "StartTime", 9) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      if (sscanf(value, "%lf", &Outputs->Output[iOutput].StartTime) != 1) 
        return xf_Error(xf_FILE_READ_ERROR);
      
      continue;
    }
    
    if (strncmp(line, "EndTime", 7) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      if (sscanf(value, "%lf", &Outputs->Output[iOutput].EndTime) != 1) 
        return xf_Error(xf_FILE_READ_ERROR);
      
      continue;
    }
    
    if (strncmp(line, "ErrTol", 6) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      if (sscanf(value, "%lf", &Outputs->Output[iOutput].ErrEst) != 1) 
        return xf_Error(xf_FILE_READ_ERROR);
      
      continue;
    }
    
    if (strncmp(line, "egrp", 4) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      if (sscanf(value, "%d", &Outputs->Output[iOutput].egrp) != 1) 
        return xf_Error(xf_FILE_READ_ERROR);
      
      continue;
    }
    
    if (strncmp(line, "elem", 4) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      if (sscanf(value, "%d", &Outputs->Output[iOutput].elem) != 1) 
        return xf_Error(xf_FILE_READ_ERROR);
      
      continue;
    }
    
    
    if (strncmp(line, "xref", 4) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ScanXReal(value, &i, Outputs->Output[iOutput].xref));
      if (ierr != xf_OK) return ierr;
      
      if (i > 3) return xf_Error(xf_INPUT_ERROR);
      
      continue;
    }
    
    
    if (strncmp(line, "SumOutputNames", 14) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, longvalue, xf_MAXLINELEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ScanXStringAlloc(longvalue, xf_MAXSTRLEN, &nSumOutput,
                                          &Outputs->Output[iOutput].SumOutputNames));
      if (ierr != xf_OK) return ierr;
      
      if ((Outputs->Output[iOutput].nSumOutput > 0) && 
          (Outputs->Output[iOutput].nSumOutput != nSumOutput)){
        xf_printf("Mismatch in SumOutput: # weights != # names.\n");
        xf_printf(" [OutputName = %s]\n", Outputs->Output[iOutput].Name);
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      Outputs->Output[iOutput].nSumOutput = nSumOutput;
      continue;
    }
    
    if (strncmp(line, "SumOutputWeights", 16) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, longvalue, xf_MAXLINELEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ScanXRealAlloc(longvalue, &nSumOutput,
                                        &Outputs->Output[iOutput].SumOutputWeights));
      if (ierr != xf_OK) return ierr;
      
      if ((Outputs->Output[iOutput].nSumOutput > 0) && 
          (Outputs->Output[iOutput].nSumOutput != nSumOutput)){
        xf_printf("Mismatch in SumOutput: # weights != # names.\n");
        xf_printf(" [OutputName = %s]\n", Outputs->Output[iOutput].Name);
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      Outputs->Output[iOutput].nSumOutput = nSumOutput;
      continue;
    }
    
    if (strncmp(line, "SumOutputErrTols", 16) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, longvalue, xf_MAXLINELEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ScanXRealAlloc(longvalue, &nSumOutput,
                                        &Outputs->Output[iOutput].SumOutputErrTols));
      if (ierr != xf_OK) return ierr;
      
      if ((Outputs->Output[iOutput].nSumOutput > 0) && 
          (Outputs->Output[iOutput].nSumOutput != nSumOutput)){
        xf_printf("Mismatch in SumOutput: # ErrTols != # names.\n");
        xf_printf(" [OutputName = %s]\n", Outputs->Output[iOutput].Name);
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      Outputs->Output[iOutput].nSumOutput = nSumOutput;
      continue;
    }    
    
    if (strncmp(line, "CutPlane", 8) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ScanXReal(value, &i, Outputs->Output[iOutput].CutPlane));
      if (ierr != xf_OK) return ierr;
      
      if (i > 4) return xf_Error(xf_INPUT_ERROR);
      
      continue;
    }
    
    
    if (strncmp(line, "BFGTitles", 9) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, longvalue, xf_MAXLINELEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ScanXStringAlloc(longvalue, xf_MAXSTRLEN,
                                          &Outputs->Output[iOutput].nBFG,
                                          &Outputs->Output[iOutput].BFGTitles));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    
    if (strncmp(line, "DumpFile", 8) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, longvalue, xf_MAXLINELEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&Outputs->Output[iOutput].DumpFile, 
                                     xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;
      
      continue;
    }
    
    xf_printf("Unrecognized line in Output block:\n");
    xf_printf("  --> %s\n", line);
    return xf_Error(xf_FILE_READ_ERROR);
    
  } while (1);
  
  if (iOutput != (nOutput-1)){
    xf_printf("Number of Output terms is less than nOutput.\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  if (!ReachedEnd){
    xf_printf("Error. Output block not terminated with ENDBLOCK\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadEqnSetSensitivity
static int 
xf_ReadEqnSetSensitivity(FILE *feqn, char **InputStrings, int *piString,
                         xf_EqnSet *EqnSet)
{
  int ierr, nSens,iType, iParam, iOut, sum, iSens, nSensSet;
  char line0[xf_MAXLINELEN], *line, OutputName[xf_MAXSTRLEN],value[xf_MAXSTRLEN];
  enum xfe_Bool ReachedEnd = xfe_False, SetOutput = xfe_False, SetPar = xfe_False;
  char key[xf_MAXSTRLEN], SensType[xf_MAXSTRLEN], ParamName[xf_MAXSTRLEN];
  xf_Output *Output = NULL;
  line=NULL;
  nSens = -1;
  iSens = 1;
  nSensSet = 0;
  iOut = iType = iParam = 0;
  do{
    sum = iOut+iType+iParam;
    if (iSens > 0){
      if (sum/3 == iSens){
        SetOutput = xfe_True;
        iSens++;
      }
      else{
        SetOutput = xfe_False; 
        SetPar = xfe_False;
      }
    }
    
    if (SetOutput) {
      //output should exist and should have been initialized in EqnSet
      ierr = xf_Error(xf_FindOutput(EqnSet, OutputName, &Output));
      if (ierr != xf_OK) return ierr;
      
      Output->nSensitivity++;
      ierr = xf_Error(xf_ReAlloc((void **)&Output->Sensitivity, Output->nSensitivity, 
                                 sizeof(xf_Sensitivity)));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_Value2Enum(SensType, xfe_SensitivityName, xfe_SensitivityLast, 
                                    (int*)&Output->Sensitivity[Output->nSensitivity-1].Type));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AllocString(&Output->Sensitivity[Output->nSensitivity-1].ParamName,
                                     xf_MAXSTRLEN, ParamName));
      if (ierr != xf_OK) return ierr;
      Output->Sensitivity[Output->nSensitivity-1].nVar = -1;//this means it has not been initialized
      Output->Sensitivity[Output->nSensitivity-1].value = NULL;
      ierr = xf_Error(xf_InitKeyValue(&(Output->Sensitivity[Output->nSensitivity-1].KeyValue)));
      if (ierr != xf_OK) return ierr;
      nSensSet++;
      SetOutput = xfe_False;
      SetPar = xfe_True;
    }

    /* read line of file */
    ierr = xf_LineFromFileOrStrings(feqn, InputStrings, piString, line0, &line);
    if (ierr == xf_END_OF_FILE) break;
    if (ierr != xf_OK) return ierr;
    
    /* Check if reached end of block */
    if (strncmp(line, "ENDBLOCK", 8) == 0){
      ReachedEnd = xfe_True; break;
    }
    
    if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line
    
    if (strncmp(line, "nSensitivity", 12) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      if (sscanf(value, "%d", &nSens) != 1) return xf_Error(xf_FILE_READ_ERROR);
      if (iSens != 1){
        xf_printf("Lines in Sensitivity block out of order.\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
      continue;
    }
    //Get output corresponding to this sensitivity
    if (strncmp(line, "OutputName", 10) == 0){
      iOut++;
      if (iOut > nSens){
        xf_printf("Number of OutputName exceeds nSensitivity.\n");
        xf_printf("Only one output is allowed per sensitivity.\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      ierr = xf_Error(xf_ReadKey(line, "=", key, OutputName, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      continue;
    }
    
    //Get sensitivity type
    if (strncmp(line, "SensitivityType", 15) == 0){
      iType++;  
      if (iType > nSens){
        xf_printf("Number of SensitivityType exceeds nSensitivity.\n");
        xf_printf("Only one type is allowed per sensitivity.\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      ierr = xf_Error(xf_ReadKey(line, "=", key, SensType, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      continue;
    }
    
    //Get parameter name
    if (strncmp(line, "ParamName", 9) == 0){
      iParam++;  
      if (iParam > nSens){
        xf_printf("Number of ParamName exceeds nSensitivity.\n");
        xf_printf("Only one parameter is allowed per sensitivity.\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
      
      ierr = xf_Error(xf_ReadKey(line, "=", key, ParamName, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      continue;
    }
    
    //Get other parameters (not recognized by previous if statements)
    if (SetPar){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_AddKeyValue(&Output->Sensitivity[Output->nSensitivity-1].KeyValue,
                            key, value, xfe_True);
      if (ierr == xf_OVERWROTE){
        xf_printf("Key \"%s\" repeated in eqnset file.\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
      else if (ierr != xf_OK) return xf_Error(ierr);
    }
    
  } while(1);
  
  if (nSens != nSensSet){
    xf_printf("Missing information for sensitivities.\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  if (iSens != nSens+1){
    xf_printf("Number of Sensitivity terms is less than nSensitivity.\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  if (!ReachedEnd){
    xf_printf("Error. Sensitivity block not terminated with ENDBLOCK\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadEqnSetFileSerial
static int 
xf_ReadEqnSetFileSerial(char *EqnSetFile, char **InputStrings, xf_EqnSet *EqnSet)
{
  int ierr, iString = 0;
  enum xfe_Bool ReadEqnSetLibrary = xfe_False;
  enum xfe_Bool done;
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  char line0[200], *line;
  FILE *feqn = NULL;
  
  /* Check input */
  if (((EqnSetFile == NULL) && (InputStrings == NULL)) || 
      ((EqnSetFile != NULL) && (InputStrings != NULL))) 
    return xf_Error(xf_INPUT_ERROR);
  
  if (EqnSetFile != NULL){
    /* Open file */
    if ((feqn = fopen(EqnSetFile, "r")) == NULL){
      xf_printf("Could not find Equation Set file: %s\n", EqnSetFile);
      xf_printf("Make sure the appropriate extension is included.\n");
      return xf_Error(xf_FILE_READ_ERROR);
    }
  }
  else{
    // Start at the first string
    iString = 0;
  }
  
  /* Identify EqnSetLibrary */
  done = xfe_False;
  do{
    ierr = xf_LineFromFileOrStrings(feqn, InputStrings, &iString, line0, &line);
    if (ierr == xf_END_OF_FILE) break;
    if (ierr != xf_OK) return ierr;
    
    if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line
    
    if (strncmp(line, "EqnSetLibrary", 13) == 0){
      if (ReadEqnSetLibrary){
        xf_printf("EqnSetLibrary repeated in .eqn file.\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
      ReadEqnSetLibrary = xfe_True;
      /* read key and value from line */
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_Alloc((void **)&EqnSet->EqnSetLibrary, xf_MAXSTRLEN, sizeof(char)));
      if (ierr != xf_OK) return ierr;
      
      strcpy(EqnSet->EqnSetLibrary, value);
    }
    
  } while (!done);
  
  // Rewind file (or strings) to beginning
  if (feqn != NULL)
    rewind(feqn);
  else
    iString = 0;
  
  /* Read blocks */
  done = xfe_False;
  do{
    ierr = xf_LineFromFileOrStrings(feqn, InputStrings, &iString, line0, &line);
    if (ierr == xf_END_OF_FILE) break;
    if (ierr != xf_OK) return ierr;
    
    if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line
    
    
    if (strncmp(line, "STARTBLOCK", 10) == 0){
      /* read in second token */
      ierr = xf_Error(xf_DesiredToken(line, " ", 1, value));
      if (ierr != xf_OK) return ierr;
      
      if (strncmp(value, "PARAM", 5) == 0){
        /* Read Parameter block */
        ierr = xf_Error(xf_ReadEqnSetParam(feqn, InputStrings, &iString, &EqnSet->KeyValue));
        if (ierr != xf_OK) return ierr;
      }
      else if (strncmp(value, "RESTERM", 7) == 0){
        /* Read Residual term block */
        ierr = xf_Error(xf_ReadEqnSetResTerms(feqn, InputStrings, &iString, &EqnSet->ResTerms));
        if (ierr != xf_OK) return ierr;
      }
      else if (strncmp(value, "IC", 2) == 0){
        /* Read IC block */
        ierr = xf_Error(xf_ReadEqnSetICs(feqn, InputStrings, &iString, &EqnSet->ICs));
        if (ierr != xf_OK) return ierr;
      }
      else if (strncmp(value, "BC", 2) == 0){
        /* Read BC block */
        ierr = xf_Error(xf_ReadEqnSetBCs(feqn, InputStrings, &iString, &EqnSet->BCs));
        if (ierr != xf_OK) return ierr;
      }
      else if (strncmp(value, "OUTPUT", 6) == 0){
        /* Read Output block */
        ierr = xf_Error(xf_ReadEqnSetOutputs(feqn, InputStrings, &iString, &EqnSet->Outputs));
        if (ierr != xf_OK) return ierr;
      }
      else if (strncmp(value, "SENSITIVITY", 11) == 0){
        /* Read Sensitivity block */
        ierr = xf_Error(xf_ReadEqnSetSensitivity(feqn, InputStrings, &iString, EqnSet));
        if (ierr != xf_OK) return ierr;
      }
    }
    
  } while (!done);
  
  if (feqn != NULL) fclose(feqn);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MergeEqnSet
int 
xf_MergeEqnSet(xf_EqnSet *EqnSet1, xf_EqnSet *EqnSet2)
{
  int ierr;
  
  // EqnSetLibrary
  if (EqnSet2->EqnSetLibrary != NULL){
    xf_Release((void *) EqnSet1->EqnSetLibrary);
    EqnSet1->EqnSetLibrary = EqnSet2->EqnSetLibrary;
    EqnSet2->EqnSetLibrary = NULL;
  }
  
  // Parameters: no overwrite warning
  ierr = xf_MergeKeyValue(&EqnSet1->KeyValue, EqnSet2->KeyValue, 0);
  if ((ierr != xf_OK) && (ierr != xf_OVERWROTE)) return xf_Error(ierr);
  
  // ResTerms
  if (EqnSet2->ResTerms != NULL){
    xf_DestroyResTerms(EqnSet1->ResTerms);
    EqnSet1->ResTerms = EqnSet2->ResTerms;
    EqnSet2->ResTerms = NULL;
  }
  
  // ICs
  if (EqnSet2->ICs != NULL){
    xf_DestroyICs(EqnSet1->ICs);
    EqnSet1->ICs = EqnSet2->ICs;
    EqnSet2->ICs = NULL;
  }
  
  // BCs
  if (EqnSet2->BCs != NULL){
    xf_DestroyBCs(EqnSet1->BCs);
    EqnSet1->BCs = EqnSet2->BCs;
    EqnSet2->BCs = NULL;
  }
  
  // Outputs
  if (EqnSet2->Outputs != NULL){
    xf_DestroyOutputs(EqnSet1->Outputs);
    EqnSet1->Outputs = EqnSet2->Outputs;
    EqnSet2->Outputs = NULL;
  }
  
  /* Make sure Residual terms, ICs, and BCs were read in */
  if (EqnSet1->ResTerms == NULL){
    xf_printf("Error, no ResTerms present in merged EqnSet structure.\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  
  if (EqnSet1->ICs == NULL){
    xf_printf("Error, no ICs present in merged EqnSet structure.\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  
  if (EqnSet1->BCs == NULL){
    xf_printf("Error, no BCs present in merged EqnSet structure.\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  
  /* Need > 0 ResTerms, ICs, BCs */
  if (EqnSet1->ResTerms->nResTerm == 0){
    xf_printf("Please specify > 0 residual terms.\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  
  if (EqnSet1->ICs->nIC == 0){
    xf_printf("Please specify > 0 initial conditions.\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  
  if (EqnSet1->BCs->nBC == 0){
    xf_printf("Warning 0 BCs in EqnSet; OK if all periodic. \n");
  }
  
  return xf_OK;
}


/*------------------------*/
/* Binary Reading/Writing */
/*------------------------*/  


/******************************************************************/
//   FUNCTION Definition: xf_WriteResTermBinary
int 
xf_WriteResTermBinary( xf_ResTerm *ResTerm, FILE *fid)
{
  int ierr, rev;
  
  rev = 0;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // Type
  ierr = xf_Error(xf_WriteStringBinary(xfe_ResTermName[ResTerm->Type], fid));
  if (ierr != xf_OK) return ierr;
  
  // KeyValue
  ierr = xf_Error(xf_WriteKeyValueBinary(ResTerm->KeyValue, fid));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadResTermBinary
int 
xf_ReadResTermBinary( FILE *fid, xf_ResTerm *ResTerm)
{
  int ierr, rev;
  
  rev = 0;  // writer revision number
  if (fread(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);
  
  // Type
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_ResTermName, xfe_ResTermLast, 
                                    (int *) &ResTerm->Type));
  if (ierr != xf_OK) return ierr;
  
  // KeyValue
  ierr = xf_Error(xf_ReadKeyValueBinary(fid, &ResTerm->KeyValue));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteResTermsBinary
static int 
xf_WriteResTermsBinary( xf_ResTerms *ResTerms, FILE *fid){
  int ierr, rev, i;
  
  rev = 0;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // nResTerm
  if (fwrite(&ResTerms->nResTerm, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // ResTerm
  for (i=0; i<ResTerms->nResTerm; i++){
    ierr = xf_Error(xf_WriteResTermBinary(ResTerms->ResTerm + i, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadResTermsBinary
static int 
xf_ReadResTermsBinary( FILE *fid, xf_ResTerms *ResTerms){
  int ierr, rev, i;
  
  rev = 0;  // writer revision number
  if (fread(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);
  
  // nResTerm
  if (fread(&ResTerms->nResTerm, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  // Allocate ResTerm
  ierr = xf_Error(xf_AllocResTerms(ResTerms, ResTerms->nResTerm));
  if (ierr != xf_OK) return ierr;
  
  // Read ResTerm
  for (i=0; i<ResTerms->nResTerm; i++){
    ierr = xf_Error(xf_ReadResTermBinary(fid, ResTerms->ResTerm + i));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_WriteICBinary
static int 
xf_WriteICBinary( xf_IC *IC, FILE *fid){
  int ierr, rev;
  enum xfe_Bool flag;
  
  rev = 3;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // Type
  ierr = xf_Error(xf_WriteStringBinary(IC->Type, fid));
  if (ierr != xf_OK) return ierr;
  
  // Function
  ierr = xf_Error(xf_WriteStringBinary(IC->Function, fid));
  if (ierr != xf_OK) return ierr;
  
  // Header
  ierr = xf_Error(xf_WriteStringBinary(IC->Header, fid));
  if (ierr != xf_OK) return ierr;
  
  // Data
  flag = (IC->Data != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_WriteStringBinary(IC->Data, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  // AlterFunction
  ierr = xf_Error(xf_WriteStringBinary(IC->AlterFunction, fid));
  if (ierr != xf_OK) return ierr;
  
  // AlterData
  flag = (IC->AlterData != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_WriteStringBinary(IC->AlterData, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  // PriorSteadySolve
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[IC->PriorSteadySolve], fid));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadICBinary
static int 
xf_ReadICBinary( FILE *fid, xf_IC *IC){
  int ierr, rev;
  enum xfe_Bool flag;
  
  if (fread(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev > 3) return xf_Error(xf_FILE_READ_ERROR);
  
  // Type
  ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &IC->Type));
  if (ierr != xf_OK) return ierr;
  
  // Function
  if (rev >= 1){
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &IC->Function));
    if (ierr != xf_OK) return ierr;
  }
  
  // Header
  if (rev >= 2){
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &IC->Header));
    if (ierr != xf_OK) return ierr;
  }
  
  // Data
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &IC->Data));
    if (ierr != xf_OK) return ierr;
  }
  else
    IC->Data = NULL;
  
  if (rev >= 3){
    // AlterFunction
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &IC->AlterFunction));
    if (ierr != xf_OK) return ierr;
    // AlterData
    ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
    if (ierr != xf_OK) return ierr;
    if (flag){
      ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &IC->AlterData));
      if (ierr != xf_OK) return ierr;
    }
    else
      IC->AlterData = NULL;
    // PriorSteadySolve
    ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, 
                                      (int *) &IC->PriorSteadySolve));
    if (ierr != xf_OK) return ierr;
    
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteICsBinary
static int 
xf_WriteICsBinary( xf_ICs *ICs, FILE *fid){
  int ierr, rev, i;
  
  rev = 0;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // nIC
  if (fwrite(&ICs->nIC, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // IC
  for (i=0; i<ICs->nIC; i++){
    ierr = xf_Error(xf_WriteICBinary(ICs->IC + i, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadICsBinary
static int 
xf_ReadICsBinary( FILE *fid, xf_ICs *ICs){
  int ierr, i, rev;
  
  rev = 0;  // writer revision number
  if (fread(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);
  
  // nIC
  if (fread(&ICs->nIC, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  // Allocate IC
  ierr = xf_Error(xf_AllocICs(ICs, ICs->nIC));
  if (ierr != xf_OK) return ierr;
  
  // Read IC
  for (i=0; i<ICs->nIC; i++){
    ierr = xf_Error(xf_ReadICBinary(fid, ICs->IC + i));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_WriteBCBinary
static int 
xf_WriteBCBinary( xf_BC *BC, FILE *fid){
  int ierr, rev;
  enum xfe_Bool flag;
  
  rev = 2;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // BFGTitle
  ierr = xf_Error(xf_WriteStringBinary(BC->BFGTitle, fid));
  if (ierr != xf_OK) return ierr;
  
  // Type
  ierr = xf_Error(xf_WriteStringBinary(BC->Type, fid));
  if (ierr != xf_OK) return ierr;
  
  // Function
  ierr = xf_Error(xf_WriteStringBinary(BC->Function, fid));
  if (ierr != xf_OK) return ierr;
  
  // Header
  ierr = xf_Error(xf_WriteStringBinary(BC->Header, fid));
  if (ierr != xf_OK) return ierr;
  
  // Data
  flag = (BC->Data != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_WriteStringBinary(BC->Data, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  // OutputLinkage
  ierr = xf_Error(xf_WriteStringBinary(BC->OutputLinkage, fid));
  if (ierr != xf_OK) return ierr;
  
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadBCBinary
static int 
xf_ReadBCBinary( FILE *fid, xf_BC *BC){
  int ierr, rev;
  enum xfe_Bool flag;
  
  if (fread(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev > 2) return xf_Error(xf_FILE_READ_ERROR);
  
  // BFGTitle
  ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &BC->BFGTitle));
  if (ierr != xf_OK) return ierr;
  
  // Type
  ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &BC->Type));
  if (ierr != xf_OK) return ierr;
  
  // Function
  ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &BC->Function));
  if (ierr != xf_OK) return ierr;
  
  // Header
  if (rev >= 1){
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &BC->Header));
    if (ierr != xf_OK) return ierr;
  }
  
  // Data
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &BC->Data));
    if (ierr != xf_OK) return ierr;
  }
  else
    BC->Data = NULL;
  
  // OutputLinkage
  if (rev >= 2){
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &BC->OutputLinkage));
    if (ierr != xf_OK) return ierr;
  }
  
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteBCsBinary
static int 
xf_WriteBCsBinary( xf_BCs *BCs, FILE *fid){
  int ierr, rev, i;
  
  rev = 0;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // nBC
  if (fwrite(&BCs->nBC, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // BC
  for (i=0; i<BCs->nBC; i++){
    ierr = xf_Error(xf_WriteBCBinary(BCs->BC + i, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadBCsBinary
static int 
xf_ReadBCsBinary( FILE *fid, xf_BCs *BCs){
  int ierr, i, rev;
  
  rev = 0;  // writer revision number
  if (fread(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);
  
  // nBC
  if (fread(&BCs->nBC, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  // Allocate BC
  ierr = xf_Error(xf_AllocBCs(BCs, BCs->nBC));
  if (ierr != xf_OK) return ierr;
  
  // Read BC
  for (i=0; i<BCs->nBC; i++){
    ierr = xf_Error(xf_ReadBCBinary(fid, BCs->BC + i));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_WriteOutputBinary
static int 
xf_WriteOutputBinary( xf_Output *Output, FILE *fid){
  int ierr, i, rev;
  enum xfe_Bool NullFlag;
  
  rev = 9;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // Name
  ierr = xf_Error(xf_WriteStringBinary(Output->Name, fid));
  if (ierr != xf_OK) return ierr;
  
  // Type
  ierr = xf_Error(xf_WriteStringBinary(xfe_OutputName[Output->Type], fid));
  if (ierr != xf_OK) return ierr;
  
  // DomainNorm
  ierr = xf_Error(xf_WriteStringBinary(xfe_DomainNormName[Output->DomainNorm], fid));
  if (ierr != xf_OK) return ierr;
  
  // TimeNorm
  ierr = xf_Error(xf_WriteStringBinary(xfe_TimeNormName[Output->TimeNorm], fid));
  if (ierr != xf_OK) return ierr;
  
  // UsesFlux
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[Output->UsesFlux], fid));
  if (ierr != xf_OK) return ierr;
  
  // ScalarName
  ierr = xf_Error(xf_WriteStringBinary(Output->ScalarName, fid));
  if (ierr != xf_OK) return ierr;
  
  // VectorName
  ierr = xf_Error(xf_WriteStringBinary(Output->VectorName, fid));
  if (ierr != xf_OK) return ierr;
  
  // Function
  ierr = xf_Error(xf_WriteStringBinary(Output->Function, fid));
  if (ierr != xf_OK) return ierr;
  
  // Data
  ierr = xf_Error(xf_WriteStringBinary(Output->Data, fid));
  if (ierr != xf_OK) return ierr;
  
  
  // nFluxComponent
  if (fwrite(&Output->nFluxComponent, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // FluxComponentNames + FluxComponentWeights
  for (i=0; i<Output->nFluxComponent; i++){
    ierr = xf_Error(xf_WriteStringBinary(Output->FluxComponentNames[i], fid));
    if (ierr != xf_OK) return ierr;
    
    if (fwrite(Output->FluxComponentWeights+i, sizeof(real), 1, fid) != 1) 
      return xf_Error(xf_FILE_WRITE_ERROR);
  }
  
  // LineCoord
  if (fwrite(Output->LineCoord, sizeof(real), 6, fid) != 6) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // CutPlane
  if (fwrite(Output->CutPlane, sizeof(real), 4, fid) != 4) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // nBFG
  if (fwrite(&Output->nBFG, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // egrp, elem, and xref (for PointValue)
  if (fwrite(&Output->egrp, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  if (fwrite(&Output->elem, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  if (fwrite(Output->xref, sizeof(real), 3, fid) != 3) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // StartTime and EndTime
  if (fwrite(&Output->StartTime, sizeof(real), 1, fid) != 1)  
    return xf_Error(xf_FILE_WRITE_ERROR);
  if (fwrite(&Output->EndTime, sizeof(real), 1, fid) != 1)
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // nSumOutput
  if (fwrite(&Output->nSumOutput, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // SumOutputNames == NULL?
  NullFlag = (Output->SumOutputNames == NULL);
  if (fwrite(&NullFlag, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // SumOutputNames + SumOutputWeights + SumOutputErrTols
  for (i=0; i<Output->nSumOutput; i++){
    if (!NullFlag){
      ierr = xf_Error(xf_WriteStringBinary(Output->SumOutputNames[i], fid));
      if (ierr != xf_OK) return ierr;
    }
    
    if (fwrite(Output->SumOutputWeights+i, sizeof(real), 1, fid) != 1) 
      return xf_Error(xf_FILE_WRITE_ERROR);
    
    if (fwrite(Output->SumOutputErrTols+i, sizeof(real), 1, fid) != 1) 
      return xf_Error(xf_FILE_WRITE_ERROR);
  }
  
  // BFGTitles
  for (i=0; i<Output->nBFG; i++){
    ierr = xf_Error(xf_WriteStringBinary(Output->BFGTitles[i], fid));
    if (ierr != xf_OK) return ierr;
  }
  
  // DumpFile
  ierr = xf_Error(xf_WriteStringBinary(Output->DumpFile, fid));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadOutputBinary
static int 
xf_ReadOutputBinary( FILE *fid, xf_Output *Output){
  int ierr, rev, i;
  enum xfe_Bool NullFlag;
  
  if (fread(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev > 9) return xf_Error(xf_FILE_READ_ERROR);
  
  // Name
  ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &Output->Name));
  if (ierr != xf_OK) return ierr;
  
  // Type
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_OutputName, xfe_OutputLast, 
                                    (int *) &Output->Type));
  if (ierr != xf_OK) return ierr;
  
  // DomainNorm
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_DomainNormName, xfe_DomainNormLast, 
                                    (int *) &Output->DomainNorm));
  if (ierr != xf_OK) return ierr;
  
  if (rev >= 7){
    // TimeNorm
    ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_TimeNormName, xfe_TimeNormLast, 
                                      (int *) &Output->TimeNorm));
    if (ierr != xf_OK) return ierr;
  }
  
  // UsesFlux
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &Output->UsesFlux));
  if (ierr != xf_OK) return ierr;
  
  // ScalarName
  ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &Output->ScalarName));
  if (ierr != xf_OK) return ierr;
  
  
  if (rev >= 4){
    // VectorName
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &Output->VectorName));
    if (ierr != xf_OK) return ierr;
  }
  
  if (rev >= 5){
    // Function
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &Output->Function));
    if (ierr != xf_OK) return ierr;
    // Data
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &Output->Data));
    if (ierr != xf_OK) return ierr;
  }
  
  
  // nFluxComponent
  if (fread(&Output->nFluxComponent, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  // Allocate FluxComponentNames + Weights
  ierr = xf_Error(xf_Alloc2((void ***)&Output->FluxComponentNames, 
                            Output->nFluxComponent, xf_MAXSTRLEN, sizeof(char)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **)&Output->FluxComponentWeights, 
                           Output->nFluxComponent, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  for (i=0; i<Output->nFluxComponent; i++)
    strcpy(Output->FluxComponentNames[i], "\0");
  
  // Read FluxComponentNames + Weights
  for (i=0; i<Output->nFluxComponent; i++){
    ierr = xf_Error(xf_ReadStringBinary(fid, xf_MAXSTRLEN, Output->FluxComponentNames[i], NULL));
    if (ierr != xf_OK) return ierr;
    
    if (fread(Output->FluxComponentWeights + i, sizeof(real), 1, fid) != 1) 
      return xf_Error(xf_FILE_READ_ERROR);
  }
  
  // LineCoord
  if (fread(Output->LineCoord, sizeof(real), 6, fid) != 6) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  // CutPlane
  if (rev >= 1){
    if (fread(Output->CutPlane, sizeof(real), 4, fid) != 4)
      return xf_Error(xf_FILE_READ_ERROR);
  }
  
  // nBFG
  if (fread(&Output->nBFG, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  // egrp, elem, and xref (for PointValue)
  if (rev >= 2){
    if (fread(&Output->egrp, sizeof(int), 1, fid) != 1) 
      return xf_Error(xf_FILE_READ_ERROR);
    if (fread(&Output->elem, sizeof(int), 1, fid) != 1) 
      return xf_Error(xf_FILE_READ_ERROR);
    if (fread(Output->xref, sizeof(real), 3, fid) != 3)
      return xf_Error(xf_FILE_READ_ERROR);
  }
  
  if (rev >= 8){
    // StartTime and EndTime
    if (fread(&Output->StartTime, sizeof(real), 1, fid) != 1) 
      return xf_Error(xf_FILE_READ_ERROR);
    if (fread(&Output->EndTime, sizeof(real), 1, fid) != 1) 
      return xf_Error(xf_FILE_READ_ERROR);
  }
  
  
  // SumOutputs
  if (rev >= 3){
    // nSumOutput
    if (fread(&Output->nSumOutput, sizeof(int), 1, fid) != 1) 
      return xf_Error(xf_FILE_READ_ERROR);
    
    // SumOutputNames == NULL?
    if (fread(&NullFlag, sizeof(int), 1, fid) != 1) 
      return xf_Error(xf_FILE_READ_ERROR);
    
    // Allocate SumOutputNames + Weights
    if (!NullFlag){
      ierr = xf_Error(xf_Alloc2((void ***)&Output->SumOutputNames, 
                                Output->nSumOutput, xf_MAXSTRLEN, sizeof(char)));
      if (ierr != xf_OK) return ierr;
      
      for (i=0; i<Output->nSumOutput; i++)
        strcpy(Output->SumOutputNames[i], "\0");
    }
    
    ierr = xf_Error(xf_Alloc((void **)&Output->SumOutputWeights, 
                             Output->nSumOutput, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    if (rev >= 9){
      ierr = xf_Error(xf_Alloc((void **)&Output->SumOutputErrTols, 
                               Output->nSumOutput, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
    
    // Read SumOutputNames + Weights
    for (i=0; i<Output->nSumOutput; i++){
      if (!NullFlag){
        ierr = xf_Error(xf_ReadStringBinary(fid, xf_MAXSTRLEN, Output->SumOutputNames[i], NULL));
        if (ierr != xf_OK) return ierr;
      }
      
      if (fread(Output->SumOutputWeights + i, sizeof(real), 1, fid) != 1) 
        return xf_Error(xf_FILE_READ_ERROR);
      if (rev >= 9){
        if (fread(Output->SumOutputErrTols + i, sizeof(real), 1, fid) != 1) 
          return xf_Error(xf_FILE_READ_ERROR);
      }
    }
    
  }
  
  
  // Allocate BFGTitles
  ierr = xf_Error(xf_Alloc2((void ***)&Output->BFGTitles,
                            Output->nBFG, xf_MAXSTRLEN, sizeof(char)));
  if (ierr != xf_OK) return ierr;
  
  for (i=0; i<Output->nBFG; i++)
    strcpy(Output->BFGTitles[i], "\0");
  
  // Read BFGTitles
  for (i=0; i<Output->nBFG; i++){
    ierr = xf_Error(xf_ReadStringBinary(fid, xf_MAXSTRLEN, Output->BFGTitles[i], NULL));
    if (ierr != xf_OK) return ierr;
  }
  
  if (rev >= 6){
    // DumpFile
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &Output->DumpFile));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteOutputsBinary
static int 
xf_WriteOutputsBinary( xf_Outputs *Outputs, FILE *fid){
  int ierr, rev, i;
  
  rev = 0;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // nOutput
  if (fwrite(&Outputs->nOutput, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // Output
  for (i=0; i<Outputs->nOutput; i++){
    ierr = xf_Error(xf_WriteOutputBinary(Outputs->Output + i, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadOutputsBinary
static int 
xf_ReadOutputsBinary( FILE *fid, xf_Outputs *Outputs){
  int ierr, i, rev;
  
  rev = 0;  // writer revision number
  if (fread(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);
  
  // nOutput
  if (fread(&Outputs->nOutput, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  // Allocate Output
  ierr = xf_Error(xf_AllocOutputs(Outputs, Outputs->nOutput));
  if (ierr != xf_OK) return ierr;
  
  // Read Output
  for (i=0; i<Outputs->nOutput; i++){
    ierr = xf_Error(xf_ReadOutputBinary(fid, Outputs->Output + i));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadSensitivitiesBinary
static int 
xf_ReadSensitivitiesBinary( FILE *fid, xf_EqnSet *EqnSet){
  int ierr, writerev, rev, nSens, iSens, is;
  char OutputName[xf_MAXSTRLEN], SensType[xf_MAXSTRLEN];
  char ParamName[xf_MAXSTRLEN];
  xf_Output *Output = NULL;
  
  writerev = 0;
  if (fread(&rev, sizeof(int), 1, fid) != 1)
    return xf_Error(xf_FILE_READ_ERROR);
  if (rev != writerev) return xf_Error(xf_FILE_READ_ERROR);
  
  //nSensitivity
  if(fread(&nSens, sizeof(int), 1, fid) != 1)
    return xf_Error(xf_FILE_READ_ERROR);
  
  if (nSens == 0) xf_Error(xf_CODE_LOGIC_ERROR);
  
  for (iSens = 0; iSens < nSens; iSens++) {
    //Get OutputName
    ierr = xf_Error(xf_ReadStringBinary(fid, xf_MAXSTRLEN, OutputName,
                                        NULL));
    if (ierr != xf_OK) return ierr;
    //find Output so we can link sensitivity
    ierr = xf_Error(xf_FindOutput(EqnSet, OutputName, &Output));
    if (ierr != xf_OK) return ierr;
    //sanity check
    if (Output->nSensitivity == 0 && Output->Sensitivity != NULL)
      return xf_Error(xf_CODE_LOGIC_ERROR);
    Output->nSensitivity++;
    ierr = xf_Error(xf_ReAlloc((void**)&Output->Sensitivity,
                               Output->nSensitivity,
                               sizeof(xf_Sensitivity)));
    if (ierr != xf_OK) return ierr;
    is = Output->nSensitivity-1;
    //SensitivityType
    ierr = xf_Error(xf_ReadStringBinary(fid, xf_MAXSTRLEN,
                                        SensType, NULL));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_Value2Enum(SensType, xfe_SensitivityName,
                                  xfe_SensitivityLast,
                                  (int*)&Output->Sensitivity[is].Type));
    if (ierr != xf_OK) return ierr;
    //ParamName
    ierr = xf_Error(xf_ReadStringBinary(fid, xf_MAXSTRLEN,
                                        ParamName, NULL));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_AllocString(&Output->Sensitivity[is].ParamName,
                                   xf_MAXSTRLEN, ParamName));
    if (ierr != xf_OK) return ierr;
    Output->Sensitivity[is].nVar = -1;//this means it has not been initialized
    Output->Sensitivity[is].value = NULL;
    ierr = xf_Error(xf_InitKeyValue(&Output->Sensitivity[is].KeyValue));
    if (ierr != xf_OK) return ierr;
    //KeyValue
    ierr = xf_Error(xf_ReadKeyValueBinary(fid,&Output->Sensitivity[is].KeyValue));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteSensitivitiesBinary
static int 
xf_WriteSensitivitiesBinary( xf_EqnSet *EqnSet, FILE *fid){
  int ierr, rev, s, o, nSens;
  enum xfe_SensitivityType SensType;
  xf_Outputs *Outputs = EqnSet->Outputs;
  
  rev = 0; // write revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  //count total number of sensitivities
  nSens = 0;
  for (o = 0; o < Outputs->nOutput; o++)
    for (s=0; s < Outputs->Output[o].nSensitivity; s++)
      nSens++;
  
  //nSensitivities
  if (fwrite(&nSens, sizeof(int), 1, fid) != 1)
    return xf_Error(xf_FILE_WRITE_ERROR);
  //loop over outputs and write sensitivities
  for (o = 0; o < Outputs->nOutput; o++){
    for (s=0; s < Outputs->Output[o].nSensitivity; s++){
      //OutputName
      ierr = xf_Error(xf_WriteStringBinary(Outputs->Output[o].Name, fid));
      if (ierr != xf_OK) return ierr;
      //SensitivityType
      SensType = Outputs->Output[o].Sensitivity[s].Type;
      ierr = xf_Error(xf_WriteStringBinary(xfe_SensitivityName[SensType], fid));
      if (ierr != xf_OK) return ierr;
      //ParamName
      ierr = xf_Error(xf_WriteStringBinary(Outputs->Output[o].Sensitivity[s].ParamName, fid));
      if (ierr != xf_OK) return ierr;
      //KeyValue
      ierr = xf_Error(xf_WriteKeyValueBinary(Outputs->Output[o].Sensitivity[s].KeyValue, fid));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteEqnSetBinarySerial
int 
xf_WriteEqnSetBinarySerial( xf_EqnSet *EqnSet, FILE *fid)
{
  int ierr, rev, si, i, iw, iSens;
  enum xfe_Bool flag;
  
  si = sizeof(int);
  
  rev = 0;  // writer revision number
  if (fwrite(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // EqnSetLibrary
  ierr = xf_Error(xf_WriteStringBinary(EqnSet->EqnSetLibrary, fid));
  if (ierr != xf_OK) return ierr;
  
  // Dim
  if (fwrite(&EqnSet->Dim, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // KeyValue
  ierr = xf_Error(xf_WriteKeyValueBinary(EqnSet->KeyValue, fid));
  if (ierr != xf_OK) return ierr;
  
  // ResTerms
  flag = (EqnSet->ResTerms != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_WriteResTermsBinary(EqnSet->ResTerms, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  // ICs
  flag = (EqnSet->ICs != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_WriteICsBinary(EqnSet->ICs, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  // BCs
  flag = (EqnSet->BCs != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_WriteBCsBinary(EqnSet->BCs, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  // Outputs
  flag = (EqnSet->Outputs != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_WriteOutputsBinary(EqnSet->Outputs, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  // StateRank
  if (fwrite(&EqnSet->StateRank, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // StateName
  flag = (EqnSet->StateName != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    for (i=0; i<EqnSet->StateRank; i++){
      ierr = xf_Error(xf_WriteStringBinary(EqnSet->StateName[i], fid));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  //Sensitivities
  iw=0;
  if (EqnSet->Outputs != NULL){
    //count the number of sensitivities
    for (i = 0; i < EqnSet->Outputs->nOutput; i++)
      for (iSens=0; iSens<EqnSet->Outputs->Output[i].nSensitivity; iSens++)
        iw++;
    flag = (iw > 0);
    ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
    if (ierr != xf_OK) return ierr;
    if (flag) {
      ierr = xf_Error(xf_WriteSensitivitiesBinary(EqnSet, fid));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadEqnSetBinarySerial
int 
xf_ReadEqnSetBinarySerial(  FILE *fid, xf_EqnSet *EqnSet)
{
  int ierr, rev, si, i;
  enum xfe_Bool flag;
  
  si = sizeof(int);
  
  rev = 0;  // writer revision number
  if (fread(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);
  
  // EqnSetLibrary
  ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &EqnSet->EqnSetLibrary));
  if (ierr != xf_OK) return ierr;
  
  // Dim
  if (fread(&EqnSet->Dim, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  
  // KeyValue
  ierr = xf_Error(xf_ReadKeyValueBinary(fid, &EqnSet->KeyValue));
  if (ierr != xf_OK) return ierr;
  
  // ResTerms
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_CreateResTerms(&EqnSet->ResTerms));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReadResTermsBinary(fid, EqnSet->ResTerms));
    if (ierr != xf_OK) return ierr;
  }
  else
    EqnSet->ResTerms = NULL;
  
  // ICs
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_CreateICs(&EqnSet->ICs));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReadICsBinary(fid, EqnSet->ICs));
    if (ierr != xf_OK) return ierr;
  }
  else
    EqnSet->ICs = NULL;
  
  // BCs
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_CreateBCs(&EqnSet->BCs));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReadBCsBinary(fid, EqnSet->BCs));
    if (ierr != xf_OK) return ierr;
  }
  else
    EqnSet->BCs = NULL;
  
  // Outputs
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_CreateOutputs(&EqnSet->Outputs));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReadOutputsBinary(fid, EqnSet->Outputs));
    if (ierr != xf_OK) return ierr;
  }
  else
    EqnSet->Outputs = NULL;
  
  // StateRank
  if (fread(&EqnSet->StateRank, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  
  // StateName
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_Alloc2((void ***) &EqnSet->StateName, EqnSet->StateRank, 
                              xf_MAXSTRLEN, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<EqnSet->StateRank; i++){
      ierr = xf_Error(xf_ReadStringBinary(fid, xf_MAXSTRLEN, EqnSet->StateName[i], NULL));
      if (ierr != xf_OK) return ierr;
    }
  }
  else
    EqnSet->StateName = NULL;
  
  //Sensitivities
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag) {
    ierr = xf_Error(xf_ReadSensitivitiesBinary(fid, EqnSet));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}



/*-----------------*/
/* Parallelization */
/*-----------------*/  


/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeResTerms
static int 
xf_ParallelizeResTerms( xf_ResTerms *ResTerms){
  
  int ierr, i;
  int myRank, nProc;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_MPI_Bcast((void *) &ResTerms->nResTerm, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  if (myRank > 0){
    ierr = xf_Error(xf_AllocResTerms(ResTerms, ResTerms->nResTerm));
    if (ierr != xf_OK) return ierr;
  }
  
  for (i=0; i<ResTerms->nResTerm; i++){
    ierr = xf_Error(xf_MPI_Bcast((void *) &ResTerms->ResTerm[i].Type, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ParallelizeKeyValue(&ResTerms->ResTerm[i].KeyValue));
    if (ierr != xf_OK) return ierr;
  } // i
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeICs
static int 
xf_ParallelizeICs( xf_ICs *ICs){
  
  int ierr, i;
  int myRank, nProc;
  int ibuf[7];
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_MPI_Bcast((void *) &ICs->nIC, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  if (myRank > 0){
    ierr = xf_Error(xf_AllocICs(ICs, ICs->nIC));
    if (ierr != xf_OK) return ierr;
  }
  
  for (i=0; i<ICs->nIC; i++){
    
    if (myRank == 0){
      ibuf[0] =  strlen(ICs->IC[i].Type)+1;
      ibuf[1] = ( (ICs->IC[i].Function == NULL) ? -1 : strlen(ICs->IC[i].Function)+1 );
      ibuf[2] = ( (ICs->IC[i].Data == NULL) ? -1 : strlen(ICs->IC[i].Data)+1 );
      ibuf[3] = ( (ICs->IC[i].Header == NULL) ? -1 : strlen(ICs->IC[i].Header)+1 );
      ibuf[4] = ( (ICs->IC[i].AlterFunction == NULL) ? -1 : strlen(ICs->IC[i].AlterFunction)+1 );
      ibuf[5] = ( (ICs->IC[i].AlterData == NULL) ? -1 : strlen(ICs->IC[i].AlterData)+1 );
      ibuf[6] = ICs->IC[i].PriorSteadySolve;
    }
    ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 7*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    if (myRank > 0){
      ierr = xf_Error(xf_Alloc((void **) &ICs->IC[i].Type, ibuf[0], sizeof(char)));
      if (ierr != xf_OK) return ierr;
      
      if (ibuf[1] > 0){
        ierr = xf_Error(xf_Alloc((void **) &ICs->IC[i].Function, ibuf[1], sizeof(char)));
        if (ierr != xf_OK) return ierr;
      }
      else
        ICs->IC[i].Function = NULL;
      
      if (ibuf[2] > 0){
        ierr = xf_Error(xf_Alloc((void **) &ICs->IC[i].Data, ibuf[2], sizeof(char)));
        if (ierr != xf_OK) return ierr;
      }
      else
        ICs->IC[i].Data = NULL;
      
      if (ibuf[3] > 0){
        ierr = xf_Error(xf_Alloc((void **) &ICs->IC[i].Header, ibuf[3], sizeof(char)));
        if (ierr != xf_OK) return ierr;
      }
      else
        ICs->IC[i].Header = NULL;
      
      if (ibuf[4] > 0){
        ierr = xf_Error(xf_Alloc((void **) &ICs->IC[i].AlterFunction, ibuf[4], sizeof(char)));
        if (ierr != xf_OK) return ierr;
      }
      else
        ICs->IC[i].AlterFunction = NULL;
      
      if (ibuf[5] > 0){
        ierr = xf_Error(xf_Alloc((void **) &ICs->IC[i].AlterData, ibuf[5], sizeof(char)));
        if (ierr != xf_OK) return ierr;
      }
      else
        ICs->IC[i].AlterData = NULL;
      
      // PriorSteadySolve
      ICs->IC[i].PriorSteadySolve = ibuf[6];
      
    }
    
    ierr = xf_Error(xf_MPI_Bcast((void *) ICs->IC[i].Type, ibuf[0]*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
    
    // Function
    if (ICs->IC[i].Function != NULL){
      ierr = xf_Error(xf_MPI_Bcast((void *) ICs->IC[i].Function, ibuf[1]*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // Data
    if (ICs->IC[i].Data != NULL){
      ierr = xf_Error(xf_MPI_Bcast((void *) ICs->IC[i].Data, ibuf[2]*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // Header
    if (ICs->IC[i].Header != NULL){
      ierr = xf_Error(xf_MPI_Bcast((void *) ICs->IC[i].Header, ibuf[3]*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // AlterFunction
    if (ICs->IC[i].AlterFunction != NULL){
      ierr = xf_Error(xf_MPI_Bcast((void *) ICs->IC[i].AlterFunction, ibuf[4]*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // AlterData
    if (ICs->IC[i].AlterData != NULL){
      ierr = xf_Error(xf_MPI_Bcast((void *) ICs->IC[i].AlterData, ibuf[5]*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
  } // i
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeBCs
static int 
xf_ParallelizeBCs( xf_BCs *BCs){
  
  int ierr, i;
  int myRank, nProc;
  int ibuf[6];
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_MPI_Bcast((void *) &BCs->nBC, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  if (myRank > 0){
    ierr = xf_Error(xf_AllocBCs(BCs, BCs->nBC));
    if (ierr != xf_OK) return ierr;
  }
  
  for (i=0; i<BCs->nBC; i++){
    
    if (myRank == 0){
      ibuf[0] =  strlen(BCs->BC[i].BFGTitle)+1;
      ibuf[1] =  strlen(BCs->BC[i].Type)+1;
      ibuf[2] = ( (BCs->BC[i].Function == NULL) ? -1 : strlen(BCs->BC[i].Function)+1 );
      ibuf[3] = ( (BCs->BC[i].Data == NULL) ? -1 : strlen(BCs->BC[i].Data)+1 );
      ibuf[4] = ( (BCs->BC[i].Header == NULL) ? -1 : strlen(BCs->BC[i].Header)+1 );
      ibuf[5] = ( (BCs->BC[i].OutputLinkage == NULL) ? -1 : strlen(BCs->BC[i].OutputLinkage)+1 );
    }
    ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 6*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    if (myRank > 0){
      ierr = xf_Error(xf_Alloc((void **) &BCs->BC[i].BFGTitle, ibuf[0], sizeof(char)));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_Alloc((void **) &BCs->BC[i].Type, ibuf[1], sizeof(char)));
      if (ierr != xf_OK) return ierr;
      
      if (ibuf[2] > 0){
        ierr = xf_Error(xf_Alloc((void **) &BCs->BC[i].Function, ibuf[2], sizeof(char)));
        if (ierr != xf_OK) return ierr;
      }
      else
        BCs->BC[i].Function = NULL;
      
      if (ibuf[3] > 0){
        ierr = xf_Error(xf_Alloc((void **) &BCs->BC[i].Data, ibuf[3], sizeof(char)));
        if (ierr != xf_OK) return ierr;
      }
      else
        BCs->BC[i].Data = NULL;
      
      if (ibuf[4] > 0){
        ierr = xf_Error(xf_Alloc((void **) &BCs->BC[i].Header, ibuf[4], sizeof(char)));
        if (ierr != xf_OK) return ierr;
      }
      else
        BCs->BC[i].Header = NULL;
      
      if (ibuf[5] > 0){
        ierr = xf_Error(xf_Alloc((void **) &BCs->BC[i].OutputLinkage, ibuf[5], sizeof(char)));
        if (ierr != xf_OK) return ierr;
      }
      else
        BCs->BC[i].OutputLinkage = NULL;
      
    }
    
    // Title
    ierr = xf_Error(xf_MPI_Bcast((void *) BCs->BC[i].BFGTitle, ibuf[0]*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
    
    // Type
    ierr = xf_Error(xf_MPI_Bcast((void *) BCs->BC[i].Type, ibuf[1]*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
    
    // Function
    if (BCs->BC[i].Function != NULL){
      ierr = xf_Error(xf_MPI_Bcast((void *) BCs->BC[i].Function, ibuf[2]*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // Data
    if (BCs->BC[i].Data != NULL){
      ierr = xf_Error(xf_MPI_Bcast((void *) BCs->BC[i].Data, ibuf[3]*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // Header
    if (BCs->BC[i].Header != NULL){
      ierr = xf_Error(xf_MPI_Bcast((void *) BCs->BC[i].Header, ibuf[4]*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // Header
    if (BCs->BC[i].OutputLinkage != NULL){
      ierr = xf_Error(xf_MPI_Bcast((void *) BCs->BC[i].OutputLinkage, ibuf[5]*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
  } // i
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeSensitivities
static int
xf_ParallelizeSensitivities( xf_Output *Output){
  int ierr, ibuf[3], myRank, iSens;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;
  
  //create the sensitivities structures on non-root procs
  if (myRank > 0){
    ierr = xf_Error(xf_Alloc((void **)&Output->Sensitivity,
                             Output->nSensitivity, sizeof(xf_Sensitivity)));
    if (ierr != xf_OK) return ierr;
  }
  for (iSens = 0; iSens < Output->nSensitivity; iSens++) {
    //broadcast the basic info (Type, strlen(ParamName), nVar)
    if (myRank == 0) {
      ibuf[0] = Output->Sensitivity[iSens].Type;
      ibuf[1] = strlen(Output->Sensitivity[iSens].ParamName)+1;
      ibuf[2] = Output->Sensitivity[iSens].nVar;
    }
    ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 3*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    if (myRank > 0) {
      Output->Sensitivity[iSens].Type = ibuf[0];
      ierr = xf_Error(xf_Alloc((void **)&Output->Sensitivity[iSens].ParamName,
                               ibuf[1], sizeof(char)));
      if (ierr != xf_OK) return ierr;
      Output->Sensitivity[iSens].nVar = ibuf[2];
      //allocate value
      ierr = xf_Error(xf_Alloc((void **)&Output->Sensitivity[iSens].value,
                               Output->Sensitivity[iSens].nVar,
                               sizeof(real)));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_InitKeyValue(&Output->Sensitivity[iSens].KeyValue));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_MPI_Bcast((void *)Output->Sensitivity[iSens].ParamName,
                                 ibuf[1]*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
    //parallelize KeyValue
    ierr = xf_Error(xf_ParallelizeKeyValue(&Output->Sensitivity[iSens].KeyValue));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeOutputs
static int 
xf_ParallelizeOutputs( xf_Outputs *Outputs){
  
  int ierr, i, j, totsize;
  int myRank, nProc;
  int ibuf[17];
  xf_Output *Output;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_MPI_Bcast((void *) &Outputs->nOutput, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  if (myRank > 0){
    ierr = xf_Error(xf_AllocOutputs(Outputs, Outputs->nOutput));
    if (ierr != xf_OK) return ierr;
  }
  
  for (i=0; i<Outputs->nOutput; i++){
    Output = Outputs->Output + i;
    
    if (myRank == 0){
      ibuf[0]  =  strlen(Output->Name)+1;
      ibuf[1]  =  Output->Type;
      ibuf[2]  =  Output->DomainNorm;
      ibuf[3]  =  Output->UsesFlux;
      ibuf[4]  =((Output->ScalarName == NULL) ? -1 : strlen(Output->ScalarName)+1 );
      ibuf[5]  =  Output->nFluxComponent;
      ibuf[6]  =  Output->nBFG;
      ibuf[7]  =  Output->egrp;
      ibuf[8]  =  Output->elem;
      ibuf[9]  =  Output->nSumOutput;
      ibuf[10] = (Output->FluxComponentMoments != NULL);
      ibuf[11] =((Output->VectorName == NULL) ? -1 : strlen(Output->VectorName)+1 );
      ibuf[12] =((Output->Function   == NULL) ? -1 : strlen(Output->Function  )+1 );
      ibuf[13] =((Output->Data       == NULL) ? -1 : strlen(Output->Data      )+1 );
      ibuf[14] =((Output->DumpFile   == NULL) ? -1 : strlen(Output->DumpFile  )+1 );
      ibuf[15] =  Output->TimeNorm;
      ibuf[16] =  Output->nSensitivity;
    }
    ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 17*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    if (myRank > 0){
      ierr = xf_Error(xf_Alloc((void **) &Output->Name, ibuf[0], sizeof(char)));
      if (ierr != xf_OK) return ierr;
      
      Output->Type       = ibuf[1];
      Output->DomainNorm = ibuf[2];
      Output->UsesFlux   = ibuf[3];
      
      if (ibuf[4] > 0){
        ierr = xf_Error(xf_Alloc((void **) &Output->ScalarName, ibuf[4], sizeof(char)));
        if (ierr != xf_OK) return ierr;
      }
      else
        Output->ScalarName = NULL;
      
      Output->nFluxComponent = ibuf[5];
      
      ierr = xf_Error(xf_Alloc2((void ***)&Output->FluxComponentNames, 
                                Output->nFluxComponent, xf_MAXSTRLEN, sizeof(char)));
      if (ierr != xf_OK) return ierr;
      
      
      ierr = xf_Error(xf_Alloc((void **)&Output->FluxComponentWeights, 
                               Output->nFluxComponent, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      
      Output->nBFG = ibuf[6];
      
      ierr = xf_Error(xf_Alloc2((void ***)&Output->BFGTitles,
                                Output->nBFG, xf_MAXSTRLEN, sizeof(char)));
      if (ierr != xf_OK) return ierr;
      
      for (j=0; j<Output->nBFG; j++)
        strcpy(Output->BFGTitles[j], "\0");
      
      Output->egrp = ibuf[7];
      Output->elem = ibuf[8];
      Output->nSumOutput = ibuf[9];
      
      ierr = xf_Error(xf_Alloc2((void ***)&Output->SumOutputNames,
                                Output->nSumOutput, xf_MAXSTRLEN, sizeof(char)));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_Alloc((void **)&Output->SumOutputWeights, 
                               Output->nSumOutput, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_Alloc((void **)&Output->SumOutputErrTols, 
                               Output->nSumOutput, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      
      if (ibuf[10]){
        ierr = xf_Error(xf_Alloc((void **)&Output->FluxComponentMoments, 
                                 Output->nFluxComponent, sizeof(int)));
        if (ierr != xf_OK) return ierr;
      }
      
      if (ibuf[11] > 0){
        ierr = xf_Error(xf_Alloc((void **) &Output->VectorName, ibuf[11], sizeof(char)));
        if (ierr != xf_OK) return ierr;
      }
      else
        Output->VectorName = NULL;
      
      if (ibuf[12] > 0){
        ierr = xf_Error(xf_Alloc((void **) &Output->Function, ibuf[12], sizeof(char)));
        if (ierr != xf_OK) return ierr;
      }
      else
        Output->Function = NULL;
      
      if (ibuf[13] > 0){
        ierr = xf_Error(xf_Alloc((void **) &Output->Data, ibuf[13], sizeof(char)));
        if (ierr != xf_OK) return ierr;
      }
      else
        Output->Data = NULL;
      
      if (ibuf[14] > 0){
        ierr = xf_Error(xf_Alloc((void **) &Output->DumpFile, ibuf[14], sizeof(char)));
        if (ierr != xf_OK) return ierr;
      }
      else
        Output->DumpFile = NULL;
      
      Output->TimeNorm = ibuf[15];
      
      Output->nSensitivity = ibuf[16];
    }
    
    // Name
    ierr = xf_Error(xf_MPI_Bcast((void *) Output->Name, ibuf[0]*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
    
    // ScalarName
    if (Output->ScalarName != NULL){
      ierr = xf_Error(xf_MPI_Bcast((void *) Output->ScalarName, ibuf[4]*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // VectorName
    if (Output->VectorName != NULL){
      ierr = xf_Error(xf_MPI_Bcast((void *) Output->VectorName, ibuf[11]*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // Function
    if (Output->Function != NULL){
      ierr = xf_Error(xf_MPI_Bcast((void *) Output->Function, ibuf[12]*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // Data
    if (Output->Data != NULL){
      ierr = xf_Error(xf_MPI_Bcast((void *) Output->Data, ibuf[13]*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // FluxComponentNames
    totsize = Output->nFluxComponent*xf_MAXSTRLEN;
    if (totsize > 0){
      ierr = xf_Error(xf_MPI_Bcast((void *) Output->FluxComponentNames[0], 
                                   totsize*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // FluxComponentWeights
    if (Output->nFluxComponent > 0){
      ierr = xf_Error(xf_MPI_Bcast((void *) Output->FluxComponentWeights,
                                   Output->nFluxComponent*sizeof(real), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // FluxComponentMoments
    if ((Output->nFluxComponent > 0) && (Output->FluxComponentMoments != NULL)){
      ierr = xf_Error(xf_MPI_Bcast((void *) Output->FluxComponentMoments,
                                   Output->nFluxComponent*sizeof(int), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // LineCoord
    ierr = xf_Error(xf_MPI_Bcast((void *) Output->LineCoord, 6*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;
    
    // CutPlane
    ierr = xf_Error(xf_MPI_Bcast((void *) Output->CutPlane, 4*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;
    
    
    // BFGTitles
    totsize = Output->nBFG*xf_MAXSTRLEN;
    if (totsize > 0){
      ierr = xf_Error(xf_MPI_Bcast((void *) Output->BFGTitles[0], totsize*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // DumpFile
    if (Output->DumpFile != NULL){
      ierr = xf_Error(xf_MPI_Bcast((void *) Output->DumpFile, ibuf[14]*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // xref
    ierr = xf_Error(xf_MPI_Bcast((void *) Output->xref, 3*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;
    
    // SumOutputNames
    totsize = Output->nSumOutput*xf_MAXSTRLEN;
    if (totsize > 0){
      ierr = xf_Error(xf_MPI_Bcast((void *) Output->SumOutputNames[0], 
                                   totsize*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // SumOutputWeights and SumOutputErrTols
    if (Output->nSumOutput > 0){
      ierr = xf_Error(xf_MPI_Bcast((void *) Output->SumOutputWeights,
                                   Output->nSumOutput*sizeof(real), 0));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_MPI_Bcast((void *) Output->SumOutputErrTols,
                                   Output->nSumOutput*sizeof(real), 0));
      if (ierr != xf_OK) return ierr;
    }
    
    // StartTime and EndTime
    ierr = xf_Error(xf_MPI_Bcast((void *) &Output->StartTime, 1*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_MPI_Bcast((void *) &Output->EndTime, 1*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;
    // Sensitivities
    if (Output->nSensitivity > 0){
      ierr = xf_Error(xf_ParallelizeSensitivities(Output));
      if (ierr != xf_OK) return ierr;
    }
    
    /* For safety, in case calling parallelize with internal
     structures filled in (e.g. from calls to CalculateOutput in
     serial), release all internal data. */
    ierr = xf_Error(xf_DestroyCutPlaneIntersect(Output->CutPlaneIntersect));
    if (ierr != xf_OK) return ierr;
    Output->CutPlaneIntersect = NULL;
    
    xf_Release( (void * ) Output->elemLocal);
    Output->elemLocal = NULL;
    
  } // i
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeEqnSet
static int 
xf_ParallelizeEqnSet( xf_EqnSet *EqnSet)
{
  /*
   PURPOSE:
   
   Parallelizes the EqnSet data structure.  Each proc gets a copy of the
   EqnSet data structure that initially only resides on proc 0.
   
   INPUTS:
   
   EqnSet : pointer to EqnSet.  Each proc must have this allocated.
   
   OUTPUTS: 
   
   None: EqnSet is parallelized
   
   RETURN:
   
   Error Code
   */
  
  int ierr;
  int myRank, nProc;
  int ibuf[8];
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (nProc == 1) return xf_OK; // nothing to do
  
  if (myRank == 0){
    ibuf[ 0] =  ((EqnSet->EqnSetLibrary == NULL) ? -1 : strlen(EqnSet->EqnSetLibrary)+1);
    ibuf[ 1] =  EqnSet->Dim;
    ibuf[ 2] = (EqnSet->ResTerms != NULL);
    ibuf[ 3] = (EqnSet->ICs != NULL);
    ibuf[ 4] = (EqnSet->BCs != NULL);
    ibuf[ 5] = (EqnSet->Outputs != NULL);
    ibuf[ 6] =  EqnSet->StateRank;
    ibuf[ 7] = (EqnSet->StateName != NULL);
  }
  
  ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 8*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  if (myRank > 0){
    
    // EqnSetLibrary
    if (EqnSet->EqnSetLibrary != NULL) xf_Release( (void *) EqnSet->EqnSetLibrary);
    if (ibuf[0] > 0){
      ierr = xf_Error(xf_Alloc((void **) &EqnSet->EqnSetLibrary, ibuf[0], sizeof(char)));
      if (ierr != xf_OK) return ierr;
    }
    else EqnSet->EqnSetLibrary = NULL;
    
    // Dim
    EqnSet->Dim = ibuf[1];
    
    // ResTerms
    if (EqnSet->ResTerms != NULL) xf_DestroyResTerms(EqnSet->ResTerms);
    if (ibuf[2]){ 
      ierr = xf_Error(xf_CreateResTerms(&EqnSet->ResTerms));
      if (ierr != xf_OK) return ierr;
    }
    else EqnSet->ResTerms = NULL;
    
    // ICs
    if (EqnSet->ICs != NULL) xf_DestroyICs(EqnSet->ICs);
    if (ibuf[3]){ 
      ierr = xf_Error(xf_CreateICs(&EqnSet->ICs));
      if (ierr != xf_OK) return ierr;
    }
    else EqnSet->ICs = NULL;
    
    // BCs
    if (EqnSet->BCs != NULL) xf_DestroyBCs(EqnSet->BCs);
    if (ibuf[4]){
      ierr = xf_Error(xf_CreateBCs(&EqnSet->BCs));
      if (ierr != xf_OK) return ierr;
    }
    else EqnSet->BCs = NULL;
    
    // Outputs
    if (EqnSet->Outputs != NULL) xf_DestroyOutputs(EqnSet->Outputs);
    if (ibuf[5]){ 
      ierr = xf_Error(xf_CreateOutputs(&EqnSet->Outputs));
      if (ierr != xf_OK) return ierr;
    }
    else EqnSet->Outputs = NULL;
    
    // StateRank
    EqnSet->StateRank = ibuf[6];
    
    // StateName
    if (ibuf[7]){
      if (EqnSet->StateName != NULL) xf_Release2( (void **) EqnSet->StateName);
      ierr = xf_Error(xf_Alloc2((void ***) &EqnSet->StateName, EqnSet->StateRank, 
                                xf_MAXSTRLEN, sizeof(char)));
      if (ierr != xf_OK) return ierr;
    }
    else EqnSet->StateName = NULL;
    
  } // end if myRank > 0
  
  // EqnSetLibrary
  if (ibuf[0] > 0){
    ierr = xf_Error(xf_MPI_Bcast((void *) EqnSet->EqnSetLibrary, ibuf[0]*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
  }
  
  // KeyValue
  ierr = xf_Error(xf_ParallelizeKeyValue(&EqnSet->KeyValue));
  if (ierr != xf_OK) return ierr;
  
  // ResTerms
  if (EqnSet->ResTerms != NULL){
    ierr = xf_Error(xf_ParallelizeResTerms(EqnSet->ResTerms));
    if (ierr != xf_OK) return ierr;
  }
  
  // ICs
  if (EqnSet->ICs != NULL){
    ierr = xf_Error(xf_ParallelizeICs(EqnSet->ICs));
    if (ierr != xf_OK) return ierr;
  }
  
  // BCs
  if (EqnSet->BCs != NULL){
    ierr = xf_Error(xf_ParallelizeBCs(EqnSet->BCs));
    if (ierr != xf_OK) return ierr;
  }
  
  // Outputs
  if (EqnSet->Outputs != NULL){
    ierr = xf_Error(xf_ParallelizeOutputs(EqnSet->Outputs));
    if (ierr != xf_OK) return ierr;
  }
  
  // StateName
  if (EqnSet->StateName != NULL){
    ierr = xf_Error(xf_MPI_Bcast((void *) EqnSet->StateName[0],
                                 EqnSet->StateRank*xf_MAXSTRLEN*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_WriteEqnSetBinary
int 
xf_WriteEqnSetBinary(  xf_EqnSet *EqnSet, FILE *fid)
{
  int ierr, terr, myRank;
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;
  
  // root writes out eqnset
  if (myRank == 0)
    terr = xf_Error(xf_WriteEqnSetBinarySerial(EqnSet, fid));
  
  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadEqnSetBinary
int 
xf_ReadEqnSetBinary(  FILE *fid, xf_EqnSet *EqnSet)
{
  int ierr, terr, myRank;
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;
  
  // root reads in eqnset
  if (myRank == 0)
    terr = xf_Error(xf_ReadEqnSetBinarySerial(fid, EqnSet));
  
  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);
  
  // parallelize eqnset
  ierr = xf_Error(xf_ParallelizeEqnSet(EqnSet));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_ReadEqnSetFile
int 
xf_ReadEqnSetFile(char *EqnSetFile, char **InputStrings, xf_EqnSet *EqnSet)
{
  int ierr, terr, myRank;
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;
  
  // root reads in eqnset
  if (myRank == 0)
    terr = xf_Error(xf_ReadEqnSetFileSerial(EqnSetFile, InputStrings, EqnSet));
  
  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);
  
  // parallelize eqnset
  ierr = xf_Error(xf_ParallelizeEqnSet(EqnSet));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReParallelizeEqnSet
int 
xf_ReParallelizeEqnSet(xf_EqnSet *EqnSet)
{
  int ierr, terr, myRank, nProc;
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // if not parallel, exit immediately
  if (nProc == 1) return xf_OK;
  
  // everyone but root destroys EqnSet
  if (myRank != 0){
    ierr = xf_Error(xf_DestroyEqnSet(EqnSet, xfe_False));
    if (ierr == xf_OK)
      ierr = xf_Error(xf_InitEqnSet(EqnSet));
  }
  else ierr = xf_OK;
  
  // reduce error, return if necessary
  ierr = xf_Error(xf_MPI_Allreduce(&ierr, 1, xfe_SizeInt, xfe_MPI_MAX));
  if (ierr != xf_OK) return ierr;
  if (ierr != xf_OK) return xf_Error(ierr);
  
  // parallelize eqnset
  ierr = xf_Error(xf_ParallelizeEqnSet(EqnSet));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}
