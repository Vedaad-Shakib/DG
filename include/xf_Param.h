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

#ifndef _xf_Param_h
#define _xf_Param_h 1

/*
  FILE:  xf_Param.h

  This file contains the headers for functions dealing with parameters.

*/

/******************************************************************/
//   FUNCTION Prototype: xf_InitKeyValue
extern int 
xf_InitKeyValue( xf_KeyValue *KeyValue);
/*
PURPOSE:

  Initializes a KeyValue structure

INPUTS:

  KeyValue: pointer to a KeyValue structure

OUTPUTS: 

  KeyValue: initialized KeyValue

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_DestroyKeyValue
extern int 
xf_DestroyKeyValue( xf_KeyValue *KeyValue);
/*
PURPOSE:

  Destroys the data within a KeyValue structure (does not release the
  pointer to KeyValue itself).  KeyValue->nKey and nKey0 are set to 0.

INPUTS:

  KeyValue: pointer to a KeyValue structure

OUTPUTS: 

  None

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_CopyKeyValue
extern int 
xf_CopyKeyValue( xf_KeyValue *KeyValue1, xf_KeyValue *KeyValue2);
/*
PURPOSE:

  Copies contents of KeyValue1 to KeyValue2.  KeyValue2 must have been
  at least initialized, as Key and Value get reallocated.

INPUTS:

  KeyValue1 : source key value

OUTPUTS: 

  KeyValue2 : destination key value

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_CreateParam
extern int 
xf_CreateParam( xf_Param **Param);
/*
PURPOSE:

  Creates an Param structure and all of its children.  Memory is
  allocated and initial (zero) values are set.

INPUTS:

  Param : pointer to Param

OUTPUTS: 

  None: Param is allocated

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_DestroyParam
extern int 
xf_DestroyParam( xf_Param *Param);
/*
PURPOSE:

  Destroys a Param structure and all of its children.  Memory is
  released

INPUTS:

  Param : pointer to Param

OUTPUTS: 

  None: Param is destroyed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_CopyParam
extern int 
xf_CopyParam( xf_Param *Param1, xf_Param *Param2);
/*
PURPOSE:

  Copies contents of Param1 to Param2.  Param2 must be initialized, so
  that Param2->KeyValue is valid.

INPUTS:

  Param1 : source param

OUTPUTS: 

  Param2 : destination param

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_AddKeyValue
extern int 
xf_AddKeyValue(xf_KeyValue *KeyValue, const char key[], 
	       const char value[], enum xfe_Bool OverWriteFlag);
/*
PURPOSE:

  Adds (key,value) pair to the KeyValue structure.  Overwrites value
  if key exists (returns xf_OVERWROTE).  Re-allocates if necessary.

INPUTS:

  KeyValue: pointer to a KeyValue structure
  key: input key string
  value: input value string
  OverWriteFlag : if False and key exists in KeyValue, the corresponding
                  value will not be overwritten.

OUTPUTS: 

  KeyValue: modified structure with (key,value) added

RETURN:

  xf_OVERWROTE: An existing key was found in KeyValue; value was overwritten
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Value2Enum
extern int 
xf_Value2Enum(const char value[], char *EName[], int ELast, 
	      int *ival);
/*
PURPOSE:

  Converts a value string to an enumerated type by matching against an
  input enumerated type name string.

INPUTS:

  value: string to convert
  EName: array of enumerated type names
  ELast: length of EName array

OUTPUTS:

  ival: returned enumerated value (if matched)

RETURN:

  xf_OK: value is matched in EName
  xf_NOT_FOUND: if value is not matched

*/

/******************************************************************/
//   FUNCTION Prototype: xf_Line2EnumsAlloc
extern int 
xf_Line2EnumsAlloc(const char *line, char *EName[], int ELast, 
		   int *nval, int **pival);
/*
PURPOSE:

  Converts a line of value strings (space separted) to a vector of
  enumerated types by matching against an input enumerated type name
  string.  Unmatched enumerated types have ther ivals set to -1

INPUTS:

  line : list of string values to convert (space separated)
  EName: array of enumerated type names
  ELast: length of EName array

OUTPUTS:

  nval     : number of values converted
  (*pival) : vector of enumerated values; gets allocated here;
             unmatched strings have this set to -1

RETURN:

  xf_OK: values are all matched from line string

*/



/******************************************************************/
//   FUNCTION Prototype: xf_DumpKeyValue
extern void
xf_DumpKeyValue(xf_KeyValue KeyValue, const char fmt[], FILE *fid);
/*
PURPOSE:

 Writes key value pairs to file stream fid, according to fmt string.

INPUTS:

  KeyValue: structure containing the key,value pairs
  fmt: format string; e.g. "%s = %s\n"
  fid: file stream

OUTPUTS:

  None

RETURN:

  None
*/

/******************************************************************/
//   FUNCTION Prototype: xf_GetKeyValue
extern int 
xf_GetKeyValue(xf_KeyValue KeyValue, const char key[], char value[]);
/*
PURPOSE:

Gets the value associated with key using the KeyValue structure

INPUTS:

  KeyValue: structure containing the key,value pairs
  key: lookup key

OUTPUTS:

  value: value corresponding to key

RETURN:

  xf_OK: if key matched
  xf_NOT_FOUND: if key is not found

*/

/******************************************************************/
//   FUNCTION Prototype: xf_GetKeyValueEnum
extern int 
xf_GetKeyValueEnum(xf_KeyValue KeyValue, const char key[], 
		   char *EnumName[], int EnumLast, int *val);
/*
PURPOSE: 

  Wrapper for xf_GetKeyValue.  After calling xf_GetValue, attempts to
  convert the string to an enumerated type.

INPUTS:

  KeyValue: structure containing the key,value pairs
  key: lookup key
  EnumName: string array of names of enumerated keys
  EnumLast: number of enumerated types

OUTPUTS:

  val: returned enumerated value

RETURN:
  xf_OK: if key matched and is an enumerated type
  Error code if key not found or if could not convert to enumerated type.
*/


/******************************************************************/
//   FUNCTION Prototype: xf_GetKeyValueBool
extern int 
xf_GetKeyValueBool(xf_KeyValue KeyValue, const char key[], enum xfe_Bool *bval);
/*
PURPOSE: 

  Wrapper for xf_GetKeyValue.  After calling xf_GetValue, attempts to
  convert the string to a boolean enumerated type.

INPUTS:

  KeyValue: structure containing the key,value pairs
  key: lookup key

OUTPUTS:

  bval: returned enumerated boolean

RETURN:
  xf_OK: if key matched and is a boolean
  Error code if key not found or if could not convert to a boolean.
*/


/******************************************************************/
//   FUNCTION Prototype: xf_GetKeyValueInt
extern int 
xf_GetKeyValueInt(xf_KeyValue KeyValue, const char key[], int *val);
/*
PURPOSE: 

  Wrapper for xf_GetKeyValue.  After calling xf_GetValue, attempts to
  convert the string to an integer.

INPUTS:

  KeyValue: structure containing the key,value pairs
  key: lookup key

OUTPUTS:

  val: returned integer

RETURN:
  xf_OK: if key matched and is an integer
  Error code if key not found or if could not convert to an integer.
*/


/******************************************************************/
//   FUNCTION Prototype: xf_GetKeyValueReal
extern int 
xf_GetKeyValueReal(xf_KeyValue KeyValue, const char key[], real *val);
/*
PURPOSE: 

  Wrapper for xf_GetKeyValue.  After calling xf_GetValue, attempts to
  convert the string to a real.

INPUTS:

  KeyValue: structure containing the key,value pairs
  key: lookup key

OUTPUTS:

  val: returned real

RETURN:
  xf_OK: if key matched and is a real
  Error code if key not found or if could not convert to an real.
*/


/******************************************************************/
//   FUNCTION Prototype: xf_SetKeyValue
extern int 
xf_SetKeyValue(xf_KeyValue KeyValue, const char key[], const char value[]);
/*
PURPOSE: 

  Looks for key in KeyValue and sets appropriate value

INPUTS:

  KeyValue: structure containing the key,value pairs
  key: lookup key
  value: value string

OUTPUTS:

  None, KeyValue is modified if key is found

RETURN:

  xf_NOT_FOUND if key not found
*/

/******************************************************************/
//   FUNCTION Prototype: xf_SetKeyValueBool
extern int 
xf_SetKeyValueBool(xf_KeyValue KeyValue, const char key[], enum xfe_Bool val);
/*
PURPOSE: 

  Wrapper for xf_SetKeyValue.  Calls xf_SetValue with a %s
  formatted string containing the name string of the boolean val.

INPUTS:

  KeyValue: structure containing the key,value pairs
  key: lookup key
  val: boolean value to set

OUTPUTS:

  None

RETURN:

  xf_NOT_FOUND if key not found
*/

/******************************************************************/
//   FUNCTION Prototype: xf_SetKeyValueInt
extern int 
xf_SetKeyValueInt(xf_KeyValue KeyValue, const char key[], int val);
/*
PURPOSE: 

  Wrapper for xf_SetKeyValue.  Calls xf_SetValue with a %d
  formatted string containing the integer val.

INPUTS:

  KeyValue: structure containing the key,value pairs
  key: lookup key
  val: integer value to set

OUTPUTS:

  None

RETURN:

  xf_NOT_FOUND if key not found
*/

/******************************************************************/
//   FUNCTION Prototype: xf_SetKeyValueReal
extern int 
xf_SetKeyValueReal(xf_KeyValue KeyValue, const char key[], real val);
/*
PURPOSE: 

  Wrapper for xf_SetKeyValue.  Calls xf_SetValue with a %.16E
  formatted string containing the real val.

INPUTS:

  KeyValue: structure containing the key,value pairs
  key: lookup key
  val: real value to set

OUTPUTS:

  None

RETURN:

  xf_NOT_FOUND if key not found
*/

/******************************************************************/
//   FUNCTION Prototype: xf_SetKeyValueVerb
extern int 
xf_SetKeyValueVerb(xf_KeyValue KeyValue, const char key[], enum xfe_Verbosity val);
/*
 PURPOSE: 
 
 Wrapper for xf_SetKeyValue.  Calls xf_SetValue with a %s
 formatted string containing the name string of the Verbosity val.
 
 INPUTS:
 
 KeyValue: structure containing the key,value pairs
 key: lookup key
 val: Verbosity value to set
 
 OUTPUTS:
 
 None
 
 RETURN:
 
 xf_NOT_FOUND if key not found
 */

/******************************************************************/
//   FUNCTION Prototype: xf_AddKeyValueList
extern int 
xf_AddKeyValueList(xf_KeyValue *pKeyValue, char **List, 
		   enum xfe_Bool OverWriteFlag, enum xfe_Bool CheckFlag);

/*
PURPOSE:

  Adds (key,value) pairs from List to the KeyValue structure.
  Overwrites value if key exists (no warning given).

INPUTS:

  pKeyValue: pointer to a KeyValue structure
  List : string array of key values:
         List = { "key1", "value1", "key2", "value2", ... "\0"}
	 Null terminating string is required.
  OverWriteFlag : if False, List keys that exist in *pKeyValue will
                  not have their value overwritten.
  CheckFlag : if True, an error occurs if a key in *pKeyValue does
              not appear in List, or if a key in List is duplicated.

OUTPUTS: 

  pKeyValue: modified structure with (key,value) pairs added

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_MergeKeyValue
extern int
xf_MergeKeyValue(xf_KeyValue *KeyValue1, xf_KeyValue KeyValue2, 
		 int WarnOverwrite);
/*
PURPOSE: 

  Merges key-value pairs in KeyValue2 with key-value pairs in
  KeyValue1.  Modifies KeyValue1.  Warns of overwritten key-value
  pairs if WarnOverwrite is 1.  Note, in case of overwriting, the
  values specified in KeyValue2 are the ones that are retained.  If
  WarnOverWrite is 2, warning are printed for every key in KeyValue2
  that does not exist in KeyValue1.

INPUTS:

  KeyValue1:  first structure containing key-value pairs.  The merged pairs are stored here.
  KeyValue2:  second structure containing key-value pairs.  Not modified.
  WarnOverwrite: if 1, warnings will be printed on every overwrite; if 2, warnings will be printed
                 for every key in KeyValue2 that does not exist in KeyValue1.  If 0, no warnings
                 will be printed.

OUTPUTS:

  KeyValue1:  Contains merged key-value pairs.

RETURN:

  xf_OVERWROTE: An existing key-value in KeyValue1 was overwritten with a value from KeyValue2
  xf_OK: no overwriting or errors.
  Other Error code: if there was a problem with the key value management
*/


/******************************************************************/
//   FUNCTION Prototype: xf_AllocFillRParam
extern int
xf_AllocFillRParam(xf_KeyValue KeyValue, int nRP, char *RPName[],
		   const real *RPDef, real **pRP, int *pnset);
/*
PURPOSE: 

  Allocates and fills real parameter vector (*prP) based on key-values
  in KeyValue.  Optional output (*pnset) is the number of params set
  from KeyValue.

INPUTS:

  KeyValue  : key-value pairs to use
  nRP       : number of real parameters to look for
  RPName    : names of the nRP real parameters of interest
  RPDef     : default values to assign to real parameters, 
              before looking in KeyValue

OUTPUTS:

  (*pRP)    : real parameter vector (is reallocated in this function)
              with values set
  (*pnset)  : number of parameters set from KeyValue

RETURN:

  Error code: if there was a problem with the key value management

*/

/******************************************************************/
//   FUNCTION Prototype: xf_AllocFillIParam
extern int
xf_AllocFillIParam(xf_KeyValue KeyValue, int nIP, char *IPName[],
		   const int *IPDef, int **pIP, int *pnset);
/*
PURPOSE: 

  Allocates and fills integer parameter vector (*pIP) based on key-values
  in KeyValue.  Optional output (*pnset) is the number of params set
  from KeyValue.

INPUTS:

  KeyValue  : key-value pairs to use
  nIP       : number of integer parameters to look for
  IPName    : names of the nIP integer parameters of interest
  IPDef     : default values to assign to integer parameters, 
              before looking in KeyValue

OUTPUTS:

  (*pIP)    : integer parameter vector (is reallocated in this function)
              with values set
  (*pnset)  : number of parameters set from KeyValue

RETURN:

  Error code: if there was a problem with the key value management

*/

/******************************************************************/
//   FUNCTION Prototype: xf_SetDefaults
extern int 
xf_SetDefaults(xf_Param *Param);
/*
PURPOSE: 

  Sets Default values for parameters in Param.  The default values are
  defined in the structure xf_DefaultParamList in xf_ParamDefault.h

INPUTS:

  None

OUTPUTS:

  Param: parameter list with keys initialized to default values

RETURN:

  Error code 
*/



/******************************************************************/
//   FUNCTION Prototype: xf_WriteKeyValueBinary
extern int 
xf_WriteKeyValueBinary( xf_KeyValue KeyValue, FILE *fid);
/*
PURPOSE:

  Writes KeyValue to a binary file

INPUTS:

  KeyValue : KeyValue structure to write
  fid : file to write to

OUTPUTS: 

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ReadKeyValueBinary
extern int 
xf_ReadKeyValueBinary( FILE *fid, xf_KeyValue *KeyValue);
/*
PURPOSE:

  Reads KeyValue from a binary file

INPUTS:

  fid : file from which to read

OUTPUTS: 

  KeyValue : pointer to KeyValue structure to read

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_WriteParamBinary
extern int 
xf_WriteParamBinary( xf_Param *Param, FILE *fid);
/*
PURPOSE:

  Writes Param to a binary file

INPUTS:

  Param : Param structure to write
  fid : file to write to

OUTPUTS: 

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_ReadParamBinary
extern int 
xf_ReadParamBinary( FILE *fid, xf_Param *Param);
/*
PURPOSE:

  Reads Param from a binary file

INPUTS:

  fid : file from which to read

OUTPUTS: 

  Param : Param structure to read

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ParallelizeKeyValue
extern int
xf_ParallelizeKeyValue( xf_KeyValue *KeyValue);
/*
PURPOSE:

  Parallelizes the KeyValue structure.  Each proc gets a copy of the
  KeyValue structure that initially only resides on proc 0.

INPUTS:

  KeyValue : pointer to KeyValue.  Each proc must have this allocated.

OUTPUTS: 

  None: KeyValue is parallelized

RETURN:

  Error Code
*/




#endif // end ifndef _xf_Param_h
