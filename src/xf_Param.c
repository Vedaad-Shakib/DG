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
  FILE:  xf_Param.c

  This file contains functions for working with parameters

*/


#include <string.h>

#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_EqnSetHook.h"


/******************************************************************/
//   FUNCTION Definition: xf_InitKeyValue
int
xf_InitKeyValue( xf_KeyValue *KeyValue){
  
  KeyValue->nKey = 0;
  KeyValue->Key   = (char **) NULL;
  KeyValue->Value = (char **) NULL;
  KeyValue->MaxStrLen = xf_MAXSTRLEN;
  KeyValue->DKey = 10;
  KeyValue->nKey0 = 0;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyKeyValue
int 
xf_DestroyKeyValue( xf_KeyValue *KeyValue){
  
  KeyValue->nKey = 0;
  xf_Release2((void **)KeyValue->Key);
  xf_Release2((void **)KeyValue->Value);
  KeyValue->nKey0 = 0;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CopyKeyValue
int
xf_CopyKeyValue( xf_KeyValue *KeyValue1, xf_KeyValue *KeyValue2 )
{
  int ierr, i, maxlen;

  KeyValue2->nKey  = KeyValue1->nKey;
  KeyValue2->DKey  = KeyValue1->DKey;
  KeyValue2->nKey0 = KeyValue1->nKey0;
  maxlen = KeyValue2->MaxStrLen = KeyValue1->MaxStrLen;

  // reallocate KeyValue2 Key and Value
  ierr = xf_Error(xf_ReAlloc2( (void ***) &KeyValue2->Key, KeyValue2->nKey0,
			       maxlen, sizeof(char)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc2( (void ***) &KeyValue2->Value, KeyValue2->nKey0,
			       maxlen, sizeof(char)));
  if (ierr != xf_OK) return ierr;

  // copy keys and values
  for (i=0; i<KeyValue2->nKey; i++){
    strncpy(KeyValue2->Key[i]  , KeyValue1->Key[i]  , maxlen);
    strncpy(KeyValue2->Value[i], KeyValue1->Value[i], maxlen);
  }


  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CreateParam
int 
xf_CreateParam( xf_Param **Param){

  int ierr;

  ierr = xf_Error(xf_Alloc((void **) Param, 1, sizeof(xf_Param)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_InitKeyValue(&((*Param)->KeyValue)));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_DestroyParam
int 
xf_DestroyParam( xf_Param *Param){
  
  int ierr;

  if (Param == NULL) return xf_OK;

  ierr = xf_Error(xf_DestroyKeyValue(&(Param->KeyValue)));
  if (ierr != xf_OK) return ierr;

  xf_Release((void *)Param);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CopyParam
int
xf_CopyParam( xf_Param *Param1, xf_Param *Param2)
{
  int ierr;

  ierr = xf_Error(xf_CopyKeyValue(&(Param1->KeyValue), &(Param2->KeyValue)));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_AddKeyValue
int 
xf_AddKeyValue( xf_KeyValue *KeyValue, const char key[], 
		const char value[], enum xfe_Bool OverWriteFlag){
  int ierr, k, len, DKey;
  char **Key;
  char **Value;

  Key = KeyValue->Key;
  Value = KeyValue->Value;

  for (k=0; k<KeyValue->nKey; k++){
    if(strcmp(Key[k], key) == 0){
      if (strlen(value) > KeyValue->MaxStrLen-1) return xf_STRING_ERROR;
      // found matching key; overwrite value if allowed to do so
      if (OverWriteFlag){
	strncpy(Value[k], value, KeyValue->MaxStrLen);
	return xf_OVERWROTE;
      }
      else return xf_OK; // key found but not overwriting
    }
  }


  /* At this point, the key was not found, so add it */

  if (KeyValue->nKey >= KeyValue->nKey0){    // need to re-allocate
    DKey = KeyValue->DKey;
    len = KeyValue->MaxStrLen;
    ierr = xf_Error(xf_ReAllocCopy2( (void ***) &KeyValue->Key, KeyValue->nKey0, len,
				     KeyValue->nKey0+DKey, len, sizeof(char)));    
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAllocCopy2( (void ***) &KeyValue->Value, KeyValue->nKey0, len,
				     KeyValue->nKey0+DKey, len, sizeof(char)));    
    if (ierr != xf_OK) return ierr;
    KeyValue->nKey0 += DKey;
  }    

  strncpy(KeyValue->Key[KeyValue->nKey], key, KeyValue->MaxStrLen);
  strncpy(KeyValue->Value[KeyValue->nKey], value, KeyValue->MaxStrLen);

  KeyValue->nKey++;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Value2Enum
int 
xf_Value2Enum(const char value[], char *EName[], int ELast, 
	      int *ival)
{
  int k;
  for (k=0; k<ELast; k++){
    if (strcmp(EName[k], value) == 0){
      *ival = k;
      return xf_OK;
    }
  }

  return xf_NOT_FOUND;
}


/******************************************************************/
//   FUNCTION Definition: xf_Line2EnumsAlloc
int 
xf_Line2EnumsAlloc(const char *line, char *EName[], int ELast, 
		   int *nval, int **pival)
{
  int ierr, k;
  char **Values = NULL;
  
  // take care of trivial case right away
  if (line == NULL){
    (*nval) = 0;
    (*pival) = NULL;
    return xf_OK;
  }

  ierr = xf_Error(xf_ScanXStringAlloc(line, xf_MAXSTRLEN, nval, &Values));
  if (ierr != xf_OK) return ierr;

  /* Allocate (*pival) */
  ierr = xf_Error(xf_Alloc((void **) pival, (*nval), sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  for (k=0; k<(*nval); k++){
    ierr = xf_Value2Enum(Values[k], EName, ELast, (*pival)+k);
    if (ierr == xf_NOT_FOUND) (*pival)[k] = -1;
    else if (ierr != xf_OK) return xf_Error(ierr);
  }

  xf_Release2( (void **) Values);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DumpKeyValue
int 
xf_DumpKeyValue(xf_KeyValue KeyValue, const char fmt[], FILE *fid)
{
  int k;
  char **Key;
  char **Value;

  Key = KeyValue.Key;
  Value = KeyValue.Value;

  for (k=0; k<KeyValue.nKey; k++)
    fprintf(fid, fmt, Key[k], Value[k]);
}

/******************************************************************/
//   FUNCTION Definition: xf_GetKeyValue
int 
xf_GetKeyValue(xf_KeyValue KeyValue, const char key[], char value[]){
  int k;
  char **Key;
  char **Value;

  Key = KeyValue.Key;
  Value = KeyValue.Value;

  for (k=0; k<KeyValue.nKey; k++){
    if(strcmp(Key[k], key) == 0){
      strncpy(value, Value[k], KeyValue.MaxStrLen);
      return xf_OK;
    }
  }
  return xf_NOT_FOUND;
}

/******************************************************************/
//   FUNCTION Definition: xf_GetKeyValueEnum
int 
xf_GetKeyValueEnum(xf_KeyValue KeyValue, const char key[], 
		   char *EnumName[], int EnumLast, int *val)
{
  int ierr;
  char value[xf_MAXSTRLEN];
  ierr = xf_GetKeyValue(KeyValue, key, value);
  if( ierr != xf_OK ) return ierr;
  ierr = xf_Error(xf_Value2Enum(value, EnumName, EnumLast, val));
  if( ierr != xf_OK ) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_GetKeyValueBool
int 
xf_GetKeyValueBool(xf_KeyValue KeyValue, const char key[], enum xfe_Bool *bval)
{
  int ierr;
  char value[xf_MAXSTRLEN];
  ierr = xf_GetKeyValue(KeyValue, key, value);
  if( ierr != xf_OK ) return ierr;
  ierr = xf_Error(xf_Value2Enum(value, xfe_BoolName, xfe_BoolLast, (int *) bval));
  if( ierr != xf_OK ) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_GetKeyValueInt
int 
xf_GetKeyValueInt(xf_KeyValue KeyValue, const char key[], int *val)
{
  int ierr;
  char value[xf_MAXSTRLEN];
  ierr = xf_GetKeyValue(KeyValue, key, value);
  if( ierr != xf_OK ) return ierr;
  ierr = sscanf(value, "%d", val);
  if (ierr != 1) return xf_STRING_ERROR;
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_GetKeyValueReal
int 
xf_GetKeyValueReal(xf_KeyValue KeyValue, const char key[], real *val)
{
  int ierr;
  char value[xf_MAXSTRLEN];
  ierr = xf_GetKeyValue(KeyValue, key, value);
  if( ierr != xf_OK ) return ierr;
  ierr = sscanf(value, "%lf", val);
  if (ierr != 1) return xf_STRING_ERROR;
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SetKeyValue
int 
xf_SetKeyValue(xf_KeyValue KeyValue, const char key[], const char value[]){
  int k;
  char **Key;
  char **Value;

  Key = KeyValue.Key;
  Value = KeyValue.Value;

  for (k=0; k<KeyValue.nKey; k++){
    if(strcmp(Key[k], key) == 0){
      strncpy(Value[k], value, KeyValue.MaxStrLen);
      return xf_OK;
    }
  }
  return xf_NOT_FOUND;
}

/******************************************************************/
//   FUNCTION Definition: xf_SetKeyValueBool
int 
xf_SetKeyValueBool(xf_KeyValue KeyValue, const char key[], enum xfe_Bool val)
{
  int ierr;
  char value[xf_MAXSTRLEN];

  sprintf(value, "%s\0", xfe_BoolName[val]);
  ierr = xf_Error(xf_SetKeyValue(KeyValue, key, value));
  if( ierr != xf_OK ) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SetKeyValueInt
int 
xf_SetKeyValueInt(xf_KeyValue KeyValue, const char key[], int val)
{
  int ierr;
  char value[xf_MAXSTRLEN];

  sprintf(value, "%d\0", val);
  ierr = xf_Error(xf_SetKeyValue(KeyValue, key, value));
  if( ierr != xf_OK ) return ierr;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_SetKeyValueReal
int 
xf_SetKeyValueReal(xf_KeyValue KeyValue, const char key[], real val)
{
  int ierr;
  char value[xf_MAXSTRLEN];

  sprintf(value, "%.16E\0", val);
  ierr = xf_Error(xf_SetKeyValue(KeyValue, key, value));
  if( ierr != xf_OK ) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SetKeyValueVerb
int 
xf_SetKeyValueVerb(xf_KeyValue KeyValue, const char key[], enum xfe_Verbosity val)
{
  int ierr;
  char value[xf_MAXSTRLEN];
  
  sprintf(value, "%s\0", xfe_VerbosityName[val]);
  ierr = xf_Error(xf_SetKeyValue(KeyValue, key, value));
  if( ierr != xf_OK ) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_AddKeyValueList
int 
xf_AddKeyValueList(xf_KeyValue *pKeyValue, char **List,
		   enum xfe_Bool OverWriteFlag,
		   enum xfe_Bool CheckFlag)
{
  int ierr, i;
  char *Key, *Value;
  char value[xf_MAXSTRLEN];
  xf_KeyValue KeyValueCheck;

  if (List == NULL) return xf_OK;

  // Initialize a KeyValue for checking
  ierr = xf_Error(xf_InitKeyValue(&KeyValueCheck));
  if (ierr != xf_OK) return ierr;

  /* Set default values */
  i = 0;
  while (strlen(List[i]) > 0){
    Key   = List[i  ];
    Value = List[i+1];
    ierr = xf_AddKeyValue(pKeyValue, Key, Value, OverWriteFlag);
    if ((ierr != xf_OK) && (ierr != xf_OVERWROTE)) return xf_Error(ierr);

    if (CheckFlag){ // Check if list contains repeats
      ierr = xf_AddKeyValue(&KeyValueCheck, Key, Value, xfe_True);
      if (ierr == xf_OVERWROTE){
	xf_printf("Error, key = %s repeated in List.\n", Key);
	return xf_Error(xf_STRING_ERROR);
      }
      if (ierr != xf_OK) return xf_Error(ierr);
    }

    i = i + 2;
    if (i > xf_ISAFETY){
      xf_printf("Error setting default parameters\n");
      xf_printf("Make sure ParamList contains a null terminating string.\n");
      return xf_Error(xf_STRING_ERROR);
    }
  }

  /* Perform a validity check if requested: warning occurs if a key in
     *pKeyValue is not present in List */
  if (CheckFlag){
    for (i = 0; i < pKeyValue->nKey; i++){
      ierr = xf_GetKeyValue(KeyValueCheck, pKeyValue->Key[i], value);
      if (ierr == xf_NOT_FOUND){
	xf_printf("Warning, key = %s does not appear in List.  Typo or old .xfa.\n",
		  pKeyValue->Key[i]);
      }
      if( ierr != xf_OK ) return ierr;    
    }
  }
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValueCheck));
  if (ierr != xf_OK) return ierr;
  

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MergeKeyValue
int
xf_MergeKeyValue(xf_KeyValue *KeyValue1, xf_KeyValue KeyValue2, 
		 int WarnOverwrite)
{
  int k, ierr, ReturnFlag;

  ReturnFlag = xf_OK;

  for (k=0; k<KeyValue2.nKey; k++){
    ierr = xf_AddKeyValue(KeyValue1, KeyValue2.Key[k], KeyValue2.Value[k], xfe_True);
    if (ierr == xf_OVERWROTE){
      ReturnFlag = xf_OVERWROTE;
      if (WarnOverwrite == 1)
	xf_printf("Warning: Key = %s is being overwritten in MergeKeyValue.\n", KeyValue2.Key[k]);
    }
    else if (ierr == xf_OK){
      if (WarnOverwrite == 2)
	xf_printf("Warning: Key = %s does not exist in default key-value list.  Possible typo?\n", 
		  KeyValue2.Key[k]);
    }
    else if (ierr != xf_OK) return xf_Error(ierr);
  }

  return ReturnFlag;
}


/******************************************************************/
//   FUNCTION Definition: xf_AllocFillRParam
int
xf_AllocFillRParam(xf_KeyValue KeyValue, int nRP, char *RPName[],
		   const real *RPDef, real **pRP, int *pnset)
{
  int ierr;
  int i, nset;
  real *RP = NULL;

  ierr = xf_Error(xf_ReAlloc( (void **) pRP, nRP, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  RP = (*pRP);

  for (i=0, nset=0; i<nRP; i++){
    ierr = xf_GetKeyValueReal(KeyValue, RPName[i], RP+i);
    if (ierr == xf_NOT_FOUND){
      if (RPDef != NULL) RP[i] = RPDef[i];
      continue;
    }
    else if (ierr != xf_OK) return xf_Error(ierr);
    nset++;
  } // i

  if (pnset != NULL) (*pnset) = nset;
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_AllocFillIParam
int
xf_AllocFillIParam(xf_KeyValue KeyValue, int nIP, char *IPName[],
		   const int *IPDef, int **pIP, int *pnset)
{
  int ierr;
  int i, nset;
  int *IP = NULL;

  ierr = xf_Error(xf_ReAlloc( (void **) pIP, nIP, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  IP = (*pIP);

  for (i=0, nset=0; i<nIP; i++){
    ierr = xf_GetKeyValueInt(KeyValue, IPName[i], IP+i);
    if (ierr == xf_NOT_FOUND){
      if (IPDef != NULL) IP[i] = IPDef[i];
      continue;
    }
    else if (ierr != xf_OK) return xf_Error(ierr);
    nset++;
  } // i

  if (pnset != NULL) (*pnset) = nset;
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteKeyValueBinary
int
xf_WriteKeyValueBinary( xf_KeyValue KeyValue, FILE *fid){
  int ierr, rev, i, si;

  si = sizeof(int);

  rev = 0;  // writer revision number
  if (fwrite(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);

  // nKey
  if (fwrite(&KeyValue.nKey, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  // MaxStrLen
  if (fwrite(&KeyValue.MaxStrLen, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  // DKey
  if (fwrite(&KeyValue.DKey, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  // nKey0
  if (fwrite(&KeyValue.nKey0, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  // Key + Value
  for (i=0; i<KeyValue.nKey; i++){
    ierr = xf_Error(xf_WriteStringBinary(KeyValue.Key[i], fid));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_WriteStringBinary(KeyValue.Value[i], fid));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadKeyValueBinary
int
xf_ReadKeyValueBinary( FILE *fid, xf_KeyValue *KeyValue){
  int ierr, rev, i, si;

  si = sizeof(int);

  // read + check revision number
  rev = 0;
  if (fread(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);

  // nKey
  if (fread(&KeyValue->nKey, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  // MaxStrLen
  if (fread(&KeyValue->MaxStrLen, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  // DKey
  if (fread(&KeyValue->DKey, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  // nKey0
  if (fread(&KeyValue->nKey0, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);


  if (KeyValue->nKey0 < KeyValue->nKey) return xf_Error(xf_OUT_OF_BOUNDS);
  
  ierr = xf_Error(xf_Alloc2((void ***) &KeyValue->Key, KeyValue->nKey0, 
			    KeyValue->MaxStrLen, sizeof(char)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc2((void ***) &KeyValue->Value, KeyValue->nKey0, 
			    KeyValue->MaxStrLen, sizeof(char)));
  if (ierr != xf_OK) return ierr;

  // Key + Value
  for (i=0; i<KeyValue->nKey; i++){
    ierr = xf_Error(xf_ReadStringBinary(fid, KeyValue->MaxStrLen, KeyValue->Key[i], NULL));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReadStringBinary(fid, KeyValue->MaxStrLen, KeyValue->Value[i], NULL));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_WriteParamBinary
int
xf_WriteParamBinary( xf_Param *Param, FILE *fid)
{
  int ierr, rev, terr;
  int myRank;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  rev = 0;  // writer revision number
  ierr = xf_Error(xf_fwrite(&rev, sizeof(int), 1, fid));
  if (ierr != xf_OK) return ierr;

  // root writes KeyValue
  if (myRank == 0)
    terr = xf_Error(xf_WriteKeyValueBinary(Param->KeyValue, fid));

  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);


  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeKeyValue
int 
xf_ParallelizeKeyValue( xf_KeyValue *KeyValue){

  int ierr, totsize;
  int myRank, nProc;
  int ibuf[4];

  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  if (nProc == 1) return xf_OK; // nothing to do

  if (myRank == 0){
    ibuf[0] = KeyValue->nKey;
    ibuf[1] = KeyValue->MaxStrLen;
    ibuf[2] = KeyValue->DKey;
    ibuf[3] = KeyValue->nKey0;
  }
  
  ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 4*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;

  if (myRank > 0){
    KeyValue->nKey      = ibuf[0];
    KeyValue->MaxStrLen = ibuf[1];
    KeyValue->DKey      = ibuf[2];
    KeyValue->nKey0     = ibuf[3];
    
    if (KeyValue->Key   != NULL) xf_Release2( (void **) KeyValue->Key);
    if (KeyValue->Value != NULL) xf_Release2( (void **) KeyValue->Value);

    ierr = xf_Error(xf_Alloc2((void ***) &KeyValue->Key, KeyValue->nKey0, 
			      KeyValue->MaxStrLen, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_Alloc2((void ***) &KeyValue->Value, KeyValue->nKey0, 
			      KeyValue->MaxStrLen, sizeof(char)));
    if (ierr != xf_OK) return ierr;    
  }
  
  totsize = KeyValue->nKey*KeyValue->MaxStrLen;
  if (totsize > 0){
    ierr = xf_Error(xf_MPI_Bcast((void *) KeyValue->Key[0], totsize*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_MPI_Bcast((void *) KeyValue->Value[0], totsize*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadParamBinary
int
xf_ReadParamBinary( FILE *fid, xf_Param *Param){
  int ierr, terr, rev;
  int myRank;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // read revision number
  ierr = xf_Error(xf_fread(fid, sizeof(int), 1, &rev));
  if (ierr != xf_OK) return ierr;

  // revision # check
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);

  // KeyValue read in by root
  if (myRank == 0)
    terr = xf_Error(xf_ReadKeyValueBinary(fid, &Param->KeyValue));

  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);

  // parallelize Param->KeyValue
  ierr = xf_Error(xf_ParallelizeKeyValue(&Param->KeyValue));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}
