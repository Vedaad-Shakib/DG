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

#ifndef _xf_EqnSet_h
#define _xf_EqnSet_h 1

/*
  FILE:  xf_EqnSet.h

  This file contains the headers for functions dealing with the EqnSet structure.

*/

/******************************************************************/
//   FUNCTION Prototype: xf_DestroyResTerm
extern int 
xf_DestroyResTerm( xf_ResTerm *ResTerm);
/*
PURPOSE:

  Destroys memory associated with ResTerm.  Does not release
  self-pointer to ResTerm.

INPUTS:

  ResTerm : residual term

OUTPUTS: None

RETURN:  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyResTerms
extern int 
xf_DestroyResTerms( xf_ResTerms *ResTerms);
/*
PURPOSE:

  Destroys memory associated with ResTerms.  Releases self pointer

INPUTS:

  ResTerms : residual terms

OUTPUTS: None

RETURN:  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_CopyResTerm
extern int 
xf_CopyResTerm( xf_ResTerm *ResTerm1, xf_ResTerm *ResTerm2);
/*
PURPOSE:

  Copies contents of ResTerm1 to ResTerm2

INPUTS:

  ResTerm1 : source residual term

OUTPUTS: 

  ResTerm2 : destination residual term

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_AddResTerm
extern int 
xf_AddResTerm( xf_ResTerms *ResTerms, enum xfe_ResTermType Type, char **KeyValueList);
/*
PURPOSE:

  Adds a residual term of type Type to ResTerms.  Key-value pairs in the
  KeyValueList are also added with the residual term.

INPUTS:
 
  ResTerms: set of residual terms
  Type   : type of residual term
  KeyValueList : null-string-terminated list of key-value pairs
                 associated with this term

OUTPUTS: None, residual term is added

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReAllocOutputs
extern int 
xf_ReAllocOutputs( xf_Outputs *Outputs, int nOutput);
/*
PURPOSE:

  Reallocates Outputs->Output to size nOutput.  If increasing # of
  outputs, the new outputs are initialized to default values.  If
  decreasing # of outputs, the last (Outputs->nOutput - nOutput)
  outputs are destroyed.

INPUTS:

  Outputs : structure in which Outputs->Output is reallocated
  nOutput : desired number of outputs

OUTPUTS: 

  None : Outputs->Output is reallocated

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_CreateEqnSet
extern int 
xf_CreateEqnSet( xf_EqnSet **EqnSet);
/*
PURPOSE:

  Creates an EqnSet structure and all of its children.  Memory is
  allocated and initial (zero) values are set.

INPUTS:

  EqnSet : pointer to EqnSet

OUTPUTS: 

  None: EqnSet is allocated

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyEqnSet
extern int 
xf_DestroyEqnSet( xf_EqnSet *EqnSet, enum xfe_Bool DestroySelf);
/*
PURPOSE:

  Destroys a EqnSet structure and all of its children.  Memory is
  released

INPUTS:

  EqnSet : pointer to EqnSet
  DestroySelf : if True, pointer to EqnSet is released

OUTPUTS: 

  None: EqnSet is destroyed

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ReadEqnSetFile
extern int 
xf_ReadEqnSetFile( char *EqnSetFile, char **InputStrings, xf_EqnSet *EqnSet);
/*
PURPOSE:

  Reads an EqnSet from a .eqn ASCII file.  The .eqn file specifies
  parameters, ResTerms, ICs, BCs, and Outputs.  Blocks not present in
  the .eqn file result in NULL being stored for that block.  Omitting
  blocks is useful for restarts from All files when no changes are
  requested in the blocks.  Input can optionally be provided via a set
  of strings, InputStrings, which contain the lines of the .eqn file.

INPUTS:

  EqnSetFile : Name of EqnSet file (or NULL)
  InputStrings: strings containing lines of input file (or NULL)

OUTPUTS: 

  EqnSet : pointer to EqnSet structure where info will be stored.
           Preallocated only in that (*EqnSet) exists.

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_MergeEqnSet
extern int 
xf_MergeEqnSet(xf_EqnSet *EqnSet1, xf_EqnSet *EqnSet2);
/*
PURPOSE:

  EqnSet1 acquires KeyValues, ResTerms, ICs, BCs, and Outputs from
  EqnSet2, for any non-null structs in EqnSet2.

    KeyValues are merged (union
    ResTerms are overwritten (if EqnSet->ResTerms != NULL)
    ICs      are overwritten (if EqnSet->ICs      != NULL)
    BCs      are overwritten (if EqnSet->Bcs      != NULL)
    Outputs  are overwritten (if EqnSet->Outputs  != NULL)

INPUTS:

  EqnSet1 : first equation set
  EqnSet2 : second equation set

OUTPUTS: 

  EqnSet1 : structures updated/overwritten with values from EqnSet2

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_WriteResTermBinary
extern int 
xf_WriteResTermBinary( xf_ResTerm *ResTerm, FILE *fid);
/*
PURPOSE:

  Writes ResTerm to a binary file

INPUTS:

  ResTerm : ResTerm structure to write
  fid : file to write to

OUTPUTS: 

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ReadResTermBinary
extern int 
xf_ReadResTermBinary( FILE *fid, xf_ResTerm *ResTerm);
/*
PURPOSE:

  Reads ResTerm from a binary file

INPUTS:

  fid : file from which to read

OUTPUTS: 

  ResTerm : ResTerm structure to read

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_WriteEqnSetBinary
extern int 
xf_WriteEqnSetBinary( xf_EqnSet *EqnSet, FILE *fid);
/*
PURPOSE:

  Writes EqnSet to a binary file

INPUTS:

  EqnSet : EqnSet structure to write
  fid : file to write to

OUTPUTS: 

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_WriteEqnSetBinarySerial
extern int 
xf_WriteEqnSetBinarySerial( xf_EqnSet *EqnSet, FILE *fid);
/*
PURPOSE:

  Writes EqnSet to a binary file.  Can be called by a single proc,
  even if running in parallel.

INPUTS:

  EqnSet : EqnSet structure to write
  fid : file to write to

OUTPUTS: 

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_WriteEqnSetBinary
extern int 
xf_WriteEqnSetBinary( xf_EqnSet *EqnSet, FILE *fid);
/*
PURPOSE:

  Writes EqnSet to a binary file.  Should be called by all procs if
  running in parallel.

INPUTS:

  EqnSet : EqnSet structure to write
  fid : file to write to

OUTPUTS: 

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_ReadEqnSetBinary
extern int 
xf_ReadEqnSetBinary( FILE *fid, xf_EqnSet *EqnSet);
/*
PURPOSE:

  Reads EqnSet from a binary file

INPUTS:

  fid : file from which to read

OUTPUTS: 

  EqnSet : EqnSet structure to read

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReParallelizeEqnSet
extern int 
xf_ReParallelizeEqnSet( xf_EqnSet *EqnSet);
/*
PURPOSE:

  EqnSet structure is destroyed on all procs except on root.  Root
  EqnSet is parallelized to all procs.

INPUTS:

  EqnSet : EqnSet structure to parallelize (root one is used)

OUTPUTS: 

  EqnSet : EqnSet structure on all procs

RETURN:

  Error Code
*/


#endif // end ifndef _xf_EqnSet_h
