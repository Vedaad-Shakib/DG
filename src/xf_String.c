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
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_Param.h"




/******************************************************************/
//   FUNCTION Definition: xf_TrimAndCheckBlank
enum xfe_Bool
xf_TrimAndCheckBlank(char **pline, int len)
{  
  int start, k;

  if (len <= 0) return xfe_True;  // blank line

  // trim leading spaces
  start = 0;
  while (((*pline)[start] == ' ') && (start < len))  start++;
  if (start >= len) 
    (*pline)[0] = '\n';
  else if (start > 0){
    for (k=start; k<len; k++) (*pline)[k-start] = (*pline)[k];
    (*pline)[len-start] = '\0';
  }

  // if blank line or comment, return True
  return ( (strncmp((*pline), "#", 1)==0) || (strncmp((*pline),"\n",1)==0) );

}

/******************************************************************/
//   FUNCTION Definition: xf_NotNull
enum xfe_Bool
xf_NotNull( const char *value)
{
  if ( (strcmp(value,"NULL")==0) ||
       (strcmp(value,"None")==0) )
    return xfe_False;
  else
    return xfe_True;
}



/******************************************************************/
//   FUNCTION Definition: xf_LineFromFileOrStrings
int 
xf_LineFromFileOrStrings( FILE *fid, char **InputStrings, int *piString,
			  char *line0, char **pline)
{  
  if (fid != NULL){
    if (fgets(line0, xf_MAXLINELEN, fid) == NULL)
      return xf_END_OF_FILE;
    (*pline) = line0;
  }
  else{
    if (InputStrings == NULL) return xf_Error(xf_INPUT_ERROR);
    (*pline) = InputStrings[(*piString)++];
    if (strlen((*pline)) <= 0) return xf_END_OF_FILE;
  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_LongLineFromFileOrStrings
int 
xf_LongLineFromFileOrStrings( FILE *fid, char **InputStrings, int *piString,
			      char *line0, char **pline)
{  
  if (fid != NULL){
    if (fgets(line0, xf_MAXLONGLINELEN, fid) == NULL)
      return xf_END_OF_FILE;
    (*pline) = line0;
  }
  else{
    if (InputStrings == NULL) return xf_Error(xf_INPUT_ERROR);
    (*pline) = InputStrings[(*piString)++];
    if (strlen((*pline)) <= 0) return xf_END_OF_FILE;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_RewindFileOrStrings
int 
xf_RewindFileOrStrings( FILE *fid, int *piString)
{  
  if (fid != NULL){
    rewind(fid);
  }
  else{
    if (piString == NULL) return xf_Error(xf_INPUT_ERROR);
    (*piString) = 0;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadKey
int
xf_ReadKey( char line[], char delim[], char key[], char value[], int maxlen)
{
  int k, startkey, endkey, startvalue, endvalue;

  startkey   = 0;
  startvalue = 0;

  /* trim line so that anything after the first '#' is ignored */
  k = 0;
  while ((line[k] != '\0') && (line[k] != '\n') && (line[k] != '#') && (k<maxlen) )
    k++;
  if (line[k] == '#') line[k] = '\0';
	 

  /* Make sure delim occurs once, and only once */
  k = 0;
  endkey = 0;
  while ((line[k] != '\0') && (line[k] != '\n')){
    if (line[k] == delim[0]){
      if (endkey != 0) {
	xf_printf("Error: > %s < Delimeter repeated.\n", line);
	return xf_Error(xf_OUT_OF_BOUNDS);
      }
      else{
	endkey = k-1;
	startvalue = k+1;
      }
    }
    k++;
  }
  endvalue = k-1;

  if (endkey == 0){
    //xf_printf("Error: > %s < No delimeter found.\n", line);
    return xf_NOT_FOUND;
  }
    

  /* make sure not exceeding maxlen for key and value */

  if (((endkey-startkey+1) > maxlen) || 
      ((endvalue-startvalue+1) > maxlen))
    return xf_Error(xf_OUT_OF_BOUNDS);


  /* trim whitespace off either side of key and value */
  
  while ((line[startkey] == ' ') && (startkey <= endkey))  startkey++;
  if (startkey > endkey) return xf_Error(xf_OUT_OF_BOUNDS); 

  while ((line[endkey] == ' ') && (startkey <= endkey))  endkey--;
  if (startkey > endkey) return xf_Error(xf_OUT_OF_BOUNDS); 

  while ((line[startvalue] == ' ') && (startvalue <= endvalue))  startvalue++;
  if (startvalue > endvalue) return xf_Error(xf_OUT_OF_BOUNDS); 

  while ((line[endvalue] == ' ') && (startvalue <= endvalue))  endvalue--;
  if (startvalue > endvalue) return xf_Error(xf_OUT_OF_BOUNDS); 
  

  /* copy key and value from line */
  
  for (k=startkey; k<=endkey; k++) key[k-startkey] = line[k];
  key[endkey-startkey+1] = '\0';
  
  for (k=startvalue; k<=endvalue; k++) value[k-startvalue] = line[k];
  value[endvalue-startvalue+1] = '\0';

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_GetNumTokens
static int 
xf_GetNumTokens(const char *line, const char *delim, int *ntok)
{
  int pos, ndelim, i;
  enum xfe_Bool indelim, indelimcur;
  
  (*ntok) = 0;

  /* Handle trivial case */
  if (line == NULL) return xf_OK;

  /* Get number of delim characters */
  ndelim = strlen(delim);
  if (ndelim <= 0) return xf_Error(xf_STRING_ERROR);

  /* Loop through characters of line and count number of tokens */
  pos = 0;
  indelim = xfe_True;
  while ((line[pos] != '\0') && (line[pos] != '\n')){
    indelimcur = xfe_False;
    for (i=0; i<ndelim; i++)
      if (line[pos] == delim[i]){
	indelimcur = xfe_True;
	break;
      }

    if ((!indelimcur) && (indelim)){
      (*ntok)++;
      indelim = xfe_False;
    }
    indelim = indelimcur;
    pos++;
  }

  return xf_OK;

}



/******************************************************************/
//   FUNCTION Definition: xf_DesiredToken
int 
xf_DesiredToken(const char *line, const char *delim, int ntok, char *value)
{
  int pos, itok, start, end;
  int ndelim, i;
  enum xfe_Bool indelim, indelimcur;
  
  /* Get number of delim characters */
  ndelim = strlen(delim);
  if (ndelim <= 0) return xf_Error(xf_STRING_ERROR);

  /* Loop until hit start of token # ntok */
  pos = 0;
  itok = -1;
  indelim = xfe_True;
  while ((line[pos] != '\0') && (line[pos] != '\n')){
    indelimcur = xfe_False;
    for (i=0; i<ndelim; i++)
      if (line[pos] == delim[i]){
	indelimcur = xfe_True;
	break;
      }

    if ((!indelimcur) && (indelim)){
      itok++; // hit start of token
      indelim = xfe_False;
    }
    indelim = indelimcur;
    if (itok == ntok) break;
    pos++;
  }

  if (itok != ntok) return xf_NOT_FOUND;

  start = pos;

  /* Get end of token */
  while ((line[pos] != '\0') && (line[pos] != '\n') && (line[pos] != delim[0])){
    pos++;
  }

  end = pos-1;

  /* copy desired token into value */
  for (pos=start; pos <= end; pos++)
    value[pos-start] = line[pos];
  value[end-start+1] = '\0';

  return xf_OK;

} /* end xf_DesiredToken */



/******************************************************************/
//   FUNCTION Definition: xf_AllocString
int 
xf_AllocString(char **dest, int maxlen, const char *source){
  int ierr, len;

  if (source == NULL){
    (*dest) = NULL;
    return xf_OK;
  }

  len = min(strlen(source)+1, maxlen);

  if (len == 0) return xf_OK;

  ierr = xf_Error(xf_Alloc((void **)dest, len, sizeof(char)));
  if (ierr != xf_OK) return ierr;

  strncpy((*dest), source, len);
  (*dest)[len-1] = '\0';
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ScanInt
int 
xf_ScanInt( const char *line, int n, int *v){
  int i, ierr;
  char value[xf_MAXSTRLEN];

  for (i=0; i<n; i++){
    ierr = xf_Error(xf_DesiredToken(line, " \t", i, value));
    if (ierr != xf_OK) return ierr;

    ierr = sscanf(value, "%d", v+i);
    if (ierr != 1) return xf_Error(xf_STRING_ERROR);
  } // i

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ScanReal
int 
xf_ScanReal( const char *line, int n, real *v){
  int i, ierr;
  char value[xf_MAXSTRLEN];

  for (i=0; i<n; i++){
    ierr = xf_Error(xf_DesiredToken(line, " \t", i, value));
    if (ierr != xf_OK) return ierr;

    ierr = sscanf(value, "%lf", v+i);
    if (ierr != 1) return xf_Error(xf_STRING_ERROR);
  } // i

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ScanXIntReal
static int 
xf_ScanXIntReal( const char *line, int *n, int *vi, real *vr)
{
  int ierr, pos, i;
  int ndelim, iv;
  const char delim[] = " \t";
  char value[xf_MAXSTRLEN];
  enum xfe_Bool indelim, indelimcur;
  
  /* return 0 if null */
  if (line == NULL){
    (*n) = 0;
    return xf_OK;
  }

  /* Error if vi, vr both null or both not null */
  if ( ((vi == NULL) && (vr == NULL)) ||  ((vi != NULL) && (vr != NULL)))
    return xf_Error(xf_INPUT_ERROR);

  /* Get number of delim characters */
  ndelim = strlen(delim);
  if (ndelim <= 0) return xf_Error(xf_STRING_ERROR);

  /* Loop until hit start of token # ntok */
  (*n) = 0;
  pos = 0;
  indelim = xfe_True;
  while ((line[pos] != '\0') && (line[pos] != '\n')){
    indelimcur = xfe_False;
    for (i=0; i<ndelim; i++)
      if (line[pos] == delim[i]){
	indelimcur = xfe_True;
	break;
      }

    if ((!indelimcur) && (indelim)){
      // hit start of token
      indelim = xfe_False;
      iv = 0;
    }
    if ((indelimcur) && (!indelim)){
      // hit end of token
      value[iv] = '\0';
      if (vi != NULL)
	ierr = sscanf(value, "%d", vi+(*n));
      else
	ierr = sscanf(value, "%lf", vr+(*n));
      if (ierr != 1) return xf_Error(xf_STRING_ERROR);
      (*n)++;
    }
    indelim = indelimcur;
    if (!indelim){
      value[iv] = line[pos];
      iv++;
    }
    pos++;
  }

  if (!indelim){
    value[iv] = '\0';
    if (vi != NULL)
      ierr = sscanf(value, "%d", vi+(*n));
    else
      ierr = sscanf(value, "%lf", vr+(*n));
    if (ierr != 1) return xf_Error(xf_STRING_ERROR);
    (*n)++;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ScanXReal
int 
xf_ScanXReal( const char *line, int *n, real *v)
{
  int ierr;

  ierr = xf_Error(xf_ScanXIntReal(line, n, NULL, v));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ScanXRealAlloc
int 
xf_ScanXRealAlloc( const char *line, int *n, real **pv)
{
  int ierr, nn;
  const char delim[] = " \t";
  
  (*n) = 0;

  /* return 0 if null */
  if (line == NULL) return xf_OK;

  /* Get number of tokens */
  ierr = xf_Error(xf_GetNumTokens(line, delim, &nn));
  if (ierr != xf_OK) return ierr;

  if (nn == 0) return xf_OK; // no tokens

  /* Allocate pv */
  ierr = xf_Error(xf_Alloc((void **) pv, nn, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  /* call ScanXReal */
  ierr = xf_Error(xf_ScanXReal(line, n, (*pv)));
  if (ierr != xf_OK) return ierr;

  /* check that nn == n */
  if ((*n) != nn) return xf_Error(xf_STRING_ERROR);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ScanXInt
int 
xf_ScanXInt( const char *line, int *n, int *v)
{
  int ierr;

  ierr = xf_Error(xf_ScanXIntReal(line, n, v, NULL));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ScanXIntAlloc
int 
xf_ScanXIntAlloc( const char *line, int *n, int **pv)
{
  int ierr, nn;
  const char delim[] = " \t";
  
  (*n) = 0;

  /* return 0 if null */
  if (line == NULL) return xf_OK;

  /* Get number of tokens */
  ierr = xf_Error(xf_GetNumTokens(line, delim, &nn));
  if (ierr != xf_OK) return ierr;

  /* Allocate pv */
  ierr = xf_Error(xf_Alloc((void **) pv, nn, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  /* call ScanXReal */
  ierr = xf_Error(xf_ScanXInt(line, n, (*pv)));
  if (ierr != xf_OK) return ierr;

  /* check that nn == n */
  if ((*n) != nn) return xf_Error(xf_STRING_ERROR);
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_ScanXStringAlloc
int 
xf_ScanXStringAlloc( const char *line, int maxlen, int *nString, char ***pStrings)
{
  int ierr, i;
  const char delim[] = " \t";
  
  (*nString) = 0;
  (*pStrings) = NULL;

  /* return 0 if null */
  if (line == NULL) return xf_OK;

  /* Get number of tokens */
  ierr = xf_Error(xf_GetNumTokens(line, delim, nString));
  if (ierr != xf_OK) return ierr;

  /* Allocate pStrings */
  ierr = xf_Error(xf_Alloc2((void ***) pStrings, (*nString), maxlen, sizeof(char)));
  if (ierr != xf_OK) return ierr;

  /* Read tokens */
  for (i=0; i<(*nString); i++){
    ierr = xf_Error(xf_DesiredToken(line, delim, i, (*pStrings)[i]));
    if (ierr != xf_OK) return ierr;
    if (strlen((*pStrings)[i]) > maxlen) return xf_Error(xf_STRING_ERROR);
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_PickRealsFromLine
int 
xf_PickRealsFromLine(const char *line, const char *keys, 
		     const real *data, real *values)
{
  int ierr, ikey, iline, nkey, nline;
  enum xfe_Bool found;
  char **Keys = NULL;
  char **Line = NULL;
  
  // take care of trivial cases right away
  if ((line == NULL) && (keys != NULL)) return xf_NOT_FOUND;
  if  (keys == NULL) return xf_OK;

  // pull off individual Keys
  ierr = xf_Error(xf_ScanXStringAlloc(keys, xf_MAXSTRLEN, &nkey, &Keys));
  if (ierr != xf_OK) return ierr;

  // break up line into individual strings
  ierr = xf_Error(xf_ScanXStringAlloc(line, xf_MAXSTRLEN, &nline, &Line));
  if (ierr != xf_OK) return ierr;
  
  // set keys
  for (ikey=0; ikey<nkey; ikey++){
    found = xfe_False;
    for (iline=0; iline<nline; iline++)
      if (strcmp(Line[iline], Keys[ikey]) == 0){
	values[ikey] = data[iline];
	found = xfe_True;
      }
    if (!found) return xf_NOT_FOUND;
  }
 
  xf_Release2( (void **) Keys);
  xf_Release2( (void **) Line);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SetRealsFromLine
int 
xf_SetRealsFromLine(const char *line, const char *keys, 
		    const real *values, real *data)
{
  int ierr, ikey, iline, nkey, nline;
  enum xfe_Bool found;
  char **Keys = NULL;
  char **Line = NULL;
  
  // take care of trivial cases right away
  if ((line == NULL) && (keys != NULL)) return xf_NOT_FOUND;
  if  (keys == NULL) return xf_OK;

  // pull off individual Keys
  ierr = xf_Error(xf_ScanXStringAlloc(keys, xf_MAXSTRLEN, &nkey, &Keys));
  if (ierr != xf_OK) return ierr;

  // break up line into individual strings
  ierr = xf_Error(xf_ScanXStringAlloc(line, xf_MAXSTRLEN, &nline, &Line));
  if (ierr != xf_OK) return ierr;
  
  // set keys
  for (ikey=0; ikey<nkey; ikey++){
    found = xfe_False;
    for (iline=0; iline<nline; iline++)
      if (strcmp(Line[iline], Keys[ikey]) == 0){
	data[iline] = values[ikey];
	found = xfe_True;
      }
    if (!found) return xf_NOT_FOUND;
  }
 
  xf_Release2( (void **) Keys);
  xf_Release2( (void **) Line);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SetStateFromHeader
int 
xf_SetStateFromHeader(const int *HeadEnum, const real *Data, 
		      int nEnum, int nData, int sr, const int *PIS, 
		      real *U, real *U_UI)
{
  int k, i, he, l;

  if (HeadEnum != NULL){
    // set desired values
    for (k=0; k<nEnum; k++){
      if ((he = HeadEnum[k]) < 0) continue; // this data index is not a state
      i = PIS[he];
      if (i < 0){
	xf_printf("State component # %d not part of current state vector!\n", he);
	return xf_Error(xf_INPUT_ERROR); 
      }
      U[i] = Data[k];
      if (U_UI != NULL) for (l=0; l<sr; l++) U_UI[i*sr+l] = 0.;
    }
  }
  else{
    return xf_Error(xf_INPUT_ERROR);
  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_WriteStringBinary
int 
xf_WriteStringBinary( const char *s, FILE *fid){
  int ierr, len;

  len = ((s == NULL) ? -1 : strlen(s));

  ierr = fwrite(&len, sizeof(int), 1, fid);
  if (ierr!=1) return xf_Error(xf_FILE_WRITE_ERROR);
  if (len >= 0){
    ierr = fwrite(s, sizeof(char), len, fid);
    if (ierr!=len) return xf_Error(xf_FILE_WRITE_ERROR);
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteStringBinaryParallel
int 
xf_WriteStringBinaryParallel( const char *s, FILE *fid){
  int ierr, len;

  len = ((s == NULL) ? -1 : strlen(s));

  ierr = xf_Error(xf_fwrite(&len, sizeof(int), 1, fid));
  if (ierr != xf_OK) return ierr;
  if (len >= 0){
    ierr = xf_Error(xf_fwrite(s, sizeof(char), len, fid));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadStringBinary
int 
xf_ReadStringBinary( FILE *fid, int maxlen, char *ss, char **ps){
  int ierr, len;
  char *s;

  ierr = fread(&len, sizeof(int), 1, fid);
  if (ierr!=1) return xf_Error(xf_FILE_READ_ERROR);

  if (len == -1){ // signifies NULL string
    if (ps != NULL) 
      (*ps) = NULL;
    else
      return xf_Error(xf_STRING_ERROR); 
    // error above because cannot make ss NULL (already preallocated)
    return xf_OK;
  }

  if (len < 0) return xf_Error(xf_OUT_OF_BOUNDS);
  if ((maxlen >= 0) && (len >= maxlen)) return xf_Error(xf_OUT_OF_BOUNDS);

  if (ps != NULL){
    ierr = xf_Error(xf_Alloc((void **) ps, len+1, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    s = (*ps);
  }
  else
    s = ss;

  ierr = fread(s, sizeof(char), len, fid);
  if (ierr!=len) return xf_Error(xf_FILE_READ_ERROR);

  s[len] = '\0'; // NULL terminating character

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadEnumBinary
int 
xf_ReadEnumBinary( FILE *fid, char *EnumName[], int EnumLast, int *val){
  int ierr, len;
  char s[xf_MAXSTRLEN];

  ierr = fread(&len, sizeof(int), 1, fid);
  if (ierr!=1) return xf_Error(xf_FILE_READ_ERROR);

  if ((len < 0)  || (len >= xf_MAXSTRLEN)) return xf_Error(xf_OUT_OF_BOUNDS);
  
  ierr = fread(s, sizeof(char), len, fid);
  if (ierr!=len) return xf_Error(xf_FILE_READ_ERROR);

  s[len] = '\0'; // NULL terminating character

  ierr = xf_Error(xf_Value2Enum(s, EnumName, EnumLast, val));
  if( ierr != xf_OK ) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadStringBinaryParallel
int 
xf_ReadStringBinaryParallel( FILE *fid, int maxlen, char *ss, char **ps){
  int ierr, len;
  char *s;

  // read len
  ierr = xf_Error(xf_fread(fid, sizeof(int), 1, &len));
  if (ierr != xf_OK) return ierr;

  if (len == -1){ // signifies NULL string
    if (ps != NULL) 
      (*ps) = NULL;
    else
      return xf_Error(xf_STRING_ERROR); 
    // error above because cannot make ss NULL (already preallocated)
    return xf_OK;
  }

  if (len < 0) return xf_Error(xf_OUT_OF_BOUNDS);
  if ((maxlen >= 0) && (len >= maxlen)) return xf_Error(xf_OUT_OF_BOUNDS);

  if (ps != NULL){
    ierr = xf_Error(xf_Alloc((void **) ps, len+1, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    s = (*ps);
  }
  else
    s = ss;

  // read s
  ierr = xf_Error(xf_fread(fid, sizeof(char), len, s));
  if (ierr != xf_OK) return ierr;

  s[len] = '\0'; // NULL terminating character

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadEnumBinaryParallel
int 
xf_ReadEnumBinaryParallel( FILE *fid, char *EnumName[], int EnumLast, int *val){
  int ierr, len;
  char s[xf_MAXSTRLEN];
  
  // read len
  ierr = xf_Error(xf_fread(fid, sizeof(int), 1, &len));
  if (ierr != xf_OK) return ierr;

  if ((len < 0)  || (len >= xf_MAXSTRLEN)) return xf_Error(xf_OUT_OF_BOUNDS);
  
  // read s
  ierr = xf_Error(xf_fread(fid, sizeof(char), len, s));
  if (ierr != xf_OK) return ierr;

  s[len] = '\0'; // NULL terminating character

  ierr = xf_Error(xf_Value2Enum(s, EnumName, EnumLast, val));
  if( ierr != xf_OK ) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Prototype: xf_SetFilePositionAfterString_Serial
int
xf_SetFilePositionAfterString_Serial(FILE *fid, enum xfe_Bool Rewind, 
                                     char *string, char *line){
  
  char Line[xf_MAXSTRLEN];
  
  if ((fid == NULL) || (string == NULL))
    return xf_INPUT_ERROR;
  
  if (Rewind)
    rewind(fid);
  
  while (!feof(fid)) {
    if (fgets(Line, xf_MAXSTRLEN, fid) == NULL)
      return xf_NOT_FOUND;
    if (strncmp(Line, string,strlen(string)) == 0) {
      if (line != NULL)
        strcpy(line,Line);
      return xf_OK;
    }
  }
  return xf_NOT_FOUND;
  
}

/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeString
int 
xf_ParallelizeString( char **ps){
  
  int ierr, ival = -1;
  int myRank, nProc;
  char *s;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (nProc <= 1) return xf_OK;

  if (ps == NULL) return xf_Error(xf_INPUT_ERROR);

  s = (*ps);

  if (myRank == 0) ival = ((s == NULL) ? -1 : strlen(s)+1 );
  ierr = xf_Error(xf_MPI_Bcast((void *) &ival, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  if (myRank > 0){
    if (ival > 0){
      ierr = xf_Error(xf_Alloc((void **) ps, ival, sizeof(char)));
      if (ierr != xf_OK) return ierr;
    }
    else (*ps) = NULL;
  }

  if (ival > 0){
    ierr = xf_Error(xf_MPI_Bcast((void *) (*ps), ival*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;

}


#if( UNIT_TEST==1 )
#include "xf_String.test.in"
#endif

