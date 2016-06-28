// This file is included in xf_Mesh.h



/******************************************************************/
//   FUNCTION Prototype: xf_CreateBamgGeomPoints
extern int 
xf_CreateBamgGeomPoints(xf_BamgGeomPoints **pBamgPoints);
/*
PURPOSE:

  Creates a structure for storing a BAMG geometry definition
  (i.e. points).

INPUTS:

  (*pBamgPoints) : pointer to structure to be allocated/initialized

OUTPUTS: 

  None: (*pBamgPoints) is allocated and initialized

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_DestroyBamgGeomPoints
extern int 
xf_DestroyBamgGeomPoints(xf_BamgGeomPoints *BamgPoints);
/*
PURPOSE:

  Destroys BAMG points structure, including self

INPUTS:

  BamgPoints : structure to destroy

OUTPUTS: None

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReadBamgFile
extern int
xf_ReadBamgFile( const char *InputFile, char **InputStrings, xf_Mesh *Mesh);
/*
 PURPOSE:
 
   Reads a Fluent/Gambit InputFile into Mesh.
 
 INPUTS:
 
   InputFile: name of file to read
   InputStrings : optional, instead of file
 
 OUTPUTS: 
 
   Mesh : pointer to Mesh structure; must have been allocated before
          the call to this function.
 
 RETURN:  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_WriteBamgGeometry_Points
extern int 
xf_WriteBamgGeometry_Points( xf_BamgGeomPoints *BP, char *OutputFile );
/*
 PURPOSE:
 
   Writes a BAMG geometry file using BAMG-specific points structure
 
 INPUTS:

   BP  : BAMG-specific structure storing geometry information (points) 
   OutputFile: name of file to write
 
 OUTPUTS:  None

 RETURN:  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_WriteBamgFile
extern int
xf_WriteBamgFile( xf_Mesh *Mesh, enum xfe_Bool AnisotropicFlag,
		  real *NodeMetric, const char *OutputFile);
/*
 PURPOSE:
 
   Writes a BAMG output file, and optionally a node-based metric.
 
 INPUTS:
 
   Mesh : pointer to Mesh structure that is written
   AnisotropicFlag : if True, metric (if provided)is anisotropic
   NodeMetric : node-based metric (optional)
   OutputFile: name of file to write
 
 OUTPUTS:  None

 RETURN:  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_RefineMeshBamg
extern int 
xf_RefineMeshBamg( xf_All *All, enum xfe_Bool AnisotropicFlag,
		   xf_Vector *EM, real *NM, xf_Mesh **pOldMesh);
/*
 PURPOSE:
 
   Refines the mesh in All->Mesh using BAMG according to the metric
   specified in EM or NM.  Data in All->Data is not "transferred" to
   this new mesh in that the linkages likely cease to be valid.
 
 INPUTS:
 
   All : pointer to All structure
   AnisotropicFlag : if True, provided metric is anisotropic
   EM  : element-based metric (optional)
   NM  : node-based metric (if EM == NULL)
   pOldMesh : Optional output -- pointer to previous/old mesh. If this is
              NULL, then the previous mesh is destroyed.
 
 OUTPUTS:  None

 RETURN:  Error Code
*/
