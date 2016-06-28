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

#ifndef _xf_String_h
#define _xf_String_h 1

/*
  FILE:  xf_String.h

  This file contains the headers for functions dealing with strings.

*/

#include <stdio.h>
#include <string.h>

/******************************************************************/
//   FUNCTION Prototype: xf_TrimAndCheckBlank
extern enum xfe_Bool
xf_TrimAndCheckBlank(char **pline, int len);
/*
PURPOSE:

  Trims front whitespace from (*pline) and checks if (*pline) is a
  blank line or a comment line.

INPUTS:

  pline : Pointer to input string.
  len: number of characters in string

OUTPUTS:

  pline : modified line with leading whitespace removed

RETURN:

  True: (*pline) is blank or a comment line
  False: (*pline) is neither a blank nor a comment line
*/


/******************************************************************/
//   FUNCTION Prototype: xf_NotNull
enum xfe_Bool
xf_NotNull( const char *value);
/*
PURPOSE:

  Returns False if value is either "NULL" or "None".

INPUTS:

  value: string to check

OUTPUTS:

  None

RETURN:

  False: string is a null string (as defined above)
  True: string is not a null string
*/



/******************************************************************/
//   FUNCTION Prototype: xf_LineFromFileOrStrings
extern int 
xf_LineFromFileOrStrings( FILE *fid, char **InputStrings, int *piString,
			  char *line0, char **pline);
/*
PURPOSE:

  Reads a line from a file stream (fid), or, if fid is NULL, from
  InputStrings, which is a set of null-terminated strings, with an end
  indicated by a "\0" string.

INPUTS:

  fid : file stream from which to read (can be NULL)
  InputStrings : collection of null-terminated strings.  Read in
                 if fid == NULL
  piString : current position in InputStrings.  Modified (+1) when
             the next string is read.
  line0 : a pointer to memory that can be used to store the line read
          in (only used for readin from a file stream)

OUTPUTS:
  
  line : pointer to line read in

RETURN:

  xf_END_OF_FILE : if end of file is reached
  xf_INPUT_ERROR : fid and InputStrings passed in both as NULL

*/


/******************************************************************/
//   FUNCTION Prototype: xf_LongLineFromFileOrStrings
extern int 
xf_LongLineFromFileOrStrings( FILE *fid, char **InputStrings, int *piString,
			      char *line0, char **pline);
/*
PURPOSE:

  Same as LineFromFileOrStrings, but reads a line of max length
  MAXLONGLINELEN

INPUTS:

  fid : file stream from which to read (can be NULL)
  InputStrings : collection of null-terminated strings.  Read in
                 if fid == NULL
  piString : current position in InputStrings.  Modified (+1) when
             the next string is read.
  line0 : a pointer to memory that can be used to store the line read
          in (only used for readin from a file stream)

OUTPUTS:
  
  line : pointer to line read in

RETURN:

  xf_END_OF_FILE : if end of file is reached
  xf_INPUT_ERROR : fid and InputStrings passed in both as NULL

*/


/******************************************************************/
//   FUNCTION Prototype: xf_RewindFileOrStrings
extern int 
xf_RewindFileOrStrings( FILE *fid, int *piString);
/*
PURPOSE:

  Rewinds a file stream or a collection of null-terminated strings
  (via pointer into position in the strings).

INPUTS:

  fid : file stream from which to read (rewound)
  piString : current position in input strings (set to zero)

OUTPUTS:
 
  (*piString) : set to zero if provided

RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReadKey
extern int
xf_ReadKey( char line[], char delim[], char key[], char value[], int maxlen);
/*
PURPOSE:

  Reads a "key" = "value" pair

INPUTS:

  line : Input string.
  delim: Token delimiter.
  maxlen: maximum length for key and value

OUTPUTS:

  key: key string with leading and trailing whitespace removed
  value: value string, also with leading/trailing whitespace removed

RETURN:

  Error code
  xf_NOT_FOUND if delimeter not found

*/

/******************************************************************/
//   FUNCTION Prototype: xf_DesiredToken
extern int
xf_DesiredToken(const char *line, const char *delim, int ntok, char *value);
/*
PURPOSE:

  Obtains token # ntok from a string which has tokens separated by a
  delimiting character[s].  For example, the following string has 4
  tokens, [0..3], separated by the delimeter " ":

  line = "This is a string\0"
  
  A call with ntok = 0 would return value = "This\0"
  A call with ntok = 2 would return value = "a\0"

  line must terminate with either "\0" or "\n".  value must be large
  enough to handle the desired token (including a terminating null
  character).


INPUTS:

  line : Input string
  delim: Token delimiters.  A string containing more than one delimiter 
         characters is allowed -- in this case, a sequence of characters
	 matching any of the characters in the delimiting strings will be
	 treated as one delimiting sequence.
  ntok : Desired token number (starting at 0)

OUTPUTS:

  value: desired token string, terminated with "\0"

RETURN:

  Error code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_AllocString
extern int 
xf_AllocString(char **dest, int maxlen, const char *source);
/*
PURPOSE:

  Allocates space for string (*dest) and copies contents of source
  into dest.  No more than maxlen characters are allocated/copied, and
  (*dest) is ensured to have a null terminating character.

  source == NULL is allowed, in which case (*dest) is set to NULL

INPUTS:

  (*dest)  : pointer to destination string; will get allocated
  maxlen   : maximum number of characters to allocate/copy
  source   : (optional) source string from which to copy characters

OUTPUTS:

  dest is allocated and filled in

RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_ScanInt
extern int 
xf_ScanInt( const char *line, int n, int *v);
/*
PURPOSE:

  Reads n integers from the string line.  Stores them in v.

INPUTS:

  line : Input string.
  n: number of integers to read

OUTPUTS:

  v: vector of n integers

RETURN:

  Error code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_ScanReal
extern int 
xf_ScanReal( const char *line, int n, real *v);
/*
PURPOSE:

  Reads n reals from the string line.  Stores them in v.

INPUTS:

  line : Input string.
  n: number of reals to read

OUTPUTS:

  v: vector of n reals (must be pre-allocated)

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ScanXInt
extern int 
xf_ScanXInt( const char *line, int *n, int *v);
/*
PURPOSE:

  Reads an a-priori-unknown number of ints from the string line.
  Stores them in v.  The number of ints read is stored in n.

INPUTS:

  line : input string.

OUTPUTS:
 
  (*n): number of ints read
  v: vector of n ints (memory must be preallocated)

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ScanXIntAlloc
extern int 
xf_ScanXIntAlloc( const char *line, int *n, int **pv);
/*
PURPOSE:

  Reads an a-priori-unknown number of ints from the string line.
  Stores them in v.  The number of ints read is stored in n.
  Allocates (*pv).

INPUTS:

  line : input string.

OUTPUTS:
 
  (*n): number of ints read
  v: vector of n ints (memory should NOT be preallocated)

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ScanXReal
extern int 
xf_ScanXReal( const char *line, int *n, real *v);
/*
PURPOSE:

  Reads an a-priori-unknown number of reals from the string line.
  Stores them in v.  The number of reals read is stored in n.

INPUTS:

  line : input string.

OUTPUTS:
 
  (*n): number of reals read
  v: vector of n reals (memory must be preallocated)

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ScanXRealAlloc
extern int 
xf_ScanXRealAlloc( const char *line, int *n, real **pv);
/*
PURPOSE:

  Reads an a-priori-unknown number of reals from the string line.
  Stores them in v.  The number of reals read is stored in n.
  Allocates (*pv).

INPUTS:

  line : input string.

OUTPUTS:
 
  (*n): number of reals read
  v: vector of n reals (memory should NOT be preallocated)

RETURN:

  Error code
*/


/*******************************************************************/
//   FUNCTION Prototype: xf_ScanXStringAlloc
extern int 
xf_ScanXStringAlloc( const char *line, int maxlen, int *nString, 
		     char ***pStrings);
/*
PURPOSE:

  Reads an a-priori-unknown number of strings from the string line.
  Stores them in (*pStrings).  A 2d array is allocated for the
  strings, with each string given a maximum of maxlen characters.  The
  number of strings read is stored in nString.  Allocates (*pStrings).

INPUTS:

  line : input string.

OUTPUTS:
 
  (*nString): number of strings read
  pStrings: Array of n strings (memory should NOT be preallocated)

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_PickRealsFromLine
extern int 
xf_PickRealsFromLine(const char *line, const char *keys, 
		     const real *data, real *values);
/*
PURPOSE:

  Pulls of data into values according to header in line and desired
  keys in keys.

  For example, if  line = "a b c d", 
               and data = [1 2 3 4]
               and keys = "b a d"
	       then values = [2 1 4]
INPUTS:

  line : list of string values to convert (space separated)
  keys : list of string keys (space sparated)
  data : real values corresponding to line

OUTPUTS:

  values : desired real values, as described above

RETURN:  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_SetRealsFromLine
extern int 
xf_SetRealsFromLine(const char *line, const char *keys, 
		    const real *values, real *data);
/*
PURPOSE:

  Writes values into data according to header in line and desired
  keys in keys.

  For example, if  line = "a b c d", 
               and data = [1 2 3 4]
               and keys = "b a d"
	       and values = [17 18 19]
               then on return, data = [18 17 3 19]

INPUTS:

  line : list of string values to convert (space separated)
  keys : list of string keys (space sparated)
  values: real values to set, one for each key
  data : original data vector

OUTPUTS:

  data : modified data vector, as described above

RETURN:  Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_SetStateFromHeader
extern int 
xf_SetStateFromHeader(const int *HeadEnum, const real *Data, 
		      int nEnum, int nData, int sr, const int *PIS, 
		      real *U, real *U_UI);
/*
PURPOSE:

  Sets components of a real vector (the state U) based on a vector of
  enumerated indices (HeadEnum) and the position in the state vector
  of each enumerated index (PIS).  The pseudocode is:

    for k = 1:nEnum,
      U[PIS[HeadEnum[k]]] = Data[k];
    end

INPUTS:

  HeadEnum : enumerated header indices
  Data : real values that are used to set the state
  nEnum : number of enumerated indices in the header
  nData : number of data entries (not used)
  sr : rank of the state vector
  PIS : position in the state vector of each valid enumerated index
  keys : list of string keys (space sparated)
  data : real values corresponding to line

OUTPUTS:

  U : state vector
  U_UI : derivative of state vector is modified, zeroed out,
         for all components of U that are set.

RETURN:  Error code

*/





/******************************************************************/
//   FUNCTION Prototype: xf_WriteStringBinary
extern int 
xf_WriteStringBinary( const char *s, FILE *fid);
/*
PURPOSE:

  Writes a string to a binary file.  The length of the string is
  written as an integer first, followed by the characters.  If the
  string is NULL, -1 is written for the length -- this will be
  understood when the string is read back in by ReadStringBinary.

INPUTS:

  s : string to write
  fid : file pointer to write to

OUTPUTS:
 
  fid : modified file pointer

RETURN:

  Error code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_WriteStringBinaryParallel
extern int 
xf_WriteStringBinaryParallel( const char *s, FILE *fid);
/*
PURPOSE:

  Identical to WriteStringBinary, except that error info is
  broadcasted to all procs upon write.  If called in parallel, needs
  to be called by all procs.
  
INPUTS:

  s : string to write
  fid : file pointer to write to

OUTPUTS:
 
  fid : modified file pointer

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReadStringBinary
extern int 
xf_ReadStringBinary(  FILE *fid, int maxlen, char *ss, char **ps);
/*
PURPOSE:

  Reads a string from a binary file which contains the length of the
  string written as an integer followed by the characters of the
  string.  If ps is not NULL, (*ps) is allocated to len+1 characters
  (to include the NULL terminating character) and the read string is
  stored in (*ps).  If ps is NULL, the read string is stored in ss,
  which must be not NULL and of adequate size to contain the read
  string.

INPUTS:

  fid : file pointer from which to read
  maxlen : maximum length allowed upon read.  Negative means no 
           max enforced

OUTPUTS:

  ss  : string read in here if ps is NULL
  ps  : if not null, allocated and string read in here

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReadEnumBinary
extern int 
xf_ReadEnumBinary( FILE *fid, char *EnumName[], int EnumLast, int *val);
/*
PURPOSE:

  Reads a string from a binary file which contains the length of the
  string written as an integer followed by the characters of the
  string; then converts the string to an enumerated type and stores
  the value in val.

INPUTS:

  fid : file pointer from which to read
  EName: array of enumerated type names
  ELast: length of EName array

OUTPUTS:

  ival: returned enumerated value (if matched)

RETURN:

  xf_NOT_FOUND: if value is not matched
  Other Error code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_ReadStringBinaryParallel
extern int 
xf_ReadStringBinaryParallel(  FILE *fid, int maxlen, char *ss, char **ps);
/*
PURPOSE:

  Identical to ReadStringBinary except that info is broadcasted to all
  procs upon read.

INPUTS:

  fid : file pointer from which to read
  maxlen : maximum length allowed upon read.  Negative means no 
           max enforced

OUTPUTS:

  ss  : string read in here if ps is NULL
  ps  : if not null, allocated and string read in here

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReadEnumBinaryParallel
extern int 
xf_ReadEnumBinaryParallel( FILE *fid, char *EnumName[], int EnumLast, int *val);
/*
PURPOSE:

  Identical to ReadEnumBinary except that info is broadcasted to all
  procs upon read.

INPUTS:

  fid : file pointer from which to read
  EName: array of enumerated type names
  ELast: length of EName array

OUTPUTS:

  ival: returned enumerated value (if matched)

RETURN:

  xf_NOT_FOUND: if value is not matched
  Other Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ParallelizeString
extern int 
xf_ParallelizeString( char **ps);
/*
PURPOSE:

  Parallelizes string (*ps) that initially just sits on the root
  processor.

INPUTS:

  (*ps) : string, allocated on root

OUTPUTS:

  (*ps) : copy (just the right length) on all procs

RETURN:  Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_SetFilePositionAfterString_Serial
extern int
xf_SetFilePositionAfterString_Serial(FILE *fid, enum xfe_Bool Rewind, 
                                     char *string, char *line);
/*
 PURPOSE:
 
 Searches for "string" in ASCII file "fid" and sets the file pointer
 at the next line of the file.  Should be called in serial or only by
 the root processor in parallel.
 
 INPUTS:
 
 fid : file pointer.
 Rewind: if True, search starts from beginning of the file.
 string: string to search for.
 line: if not NULL, the function returns the line that contains
       "string"
 
 OUTPUTS: 
 
 None: pointer "fid" is modified
 
 RETURN:
 Error Code
 */

#endif // end ifndef _xf_String_h
