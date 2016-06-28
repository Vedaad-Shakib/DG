//Yu's own unsteady adaptive solver for compressible navier-stokes. 
#include "xf_AllStruct.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_SolverStruct.h"
#include "xf_SolverTools.h"
#include "xf_Data.h"
#include "xf_Param.h"
//#include "xfYu_Residual.h"
#include "xf_Residual.h"
#include "xf_ResidualStab.h"
#include "xf_LinearSolver.h"
#include "xf_Line.h"
#include "xf_Basis.h"
#include "xf_Math.h"
#include "xf_DataMath.h"
#include "xf_Memory.h"
#include "xf_Quad.h"
#include "xf_EqnSetHook.h"
#include "xf_MeshTools.h"
#include "xf_Log.h"
#include "xf_Output.h"
#include "xf_Adapt.h"
#include "xf_AdaptStruct.h"
#include "xf_ErrEst.h"
#include "xf_All.h"
#include "xf_Penalty.h"
#include "xf_MeshMotion.h"
#include "xfYu_Limiter.h"
#include "xfYu_Solver.h"
#include "xf_Partition.h"
#include "xfYu_EntropyBounding.h"
//testing parameter for unsteady adaptation 

static real Yu_time_step_constaint;  //the smallest allowable time step size

static int MaxOrder = 4;
static int MinOrder = 1;
static int SpongeOrder = 0;
static real refinefrac = 0.10;
static real coarsenfrac = 0.10;
static real Yu_adapt_time_size = -1.0;
static real xfa_file_write_inv = -1.0;
static real Yu_load_balance_time_size = 1.e+10;
//log from Jan 31

// For SciTech 2016 to maintain fixed # DOF for comparison of different feature-based indicators
static int MaxDOF = 55000; 
static int FirstAdapt = 1; 
static int InitOrder;
static int ElemType; // 1 for tri, 2 for quad, 3 for tet, 4 for hex - get from Model
static int totalDOF;
static xf_Vector *VOrder;
static int nref; // replace nelemref
static int ncoarse; // replace nelemcoarse
static int ntest = 0;
static int nIsland, nVoid;
static int totDOFchk = 0;
// static int ncomp = 5;

static int
ec_GetDOF(int p, int *tmpnDOF)
{
	// gets number for DOF for one element with given polynomial order p
	if (ElemType == 1) (*tmpnDOF) = (p+1)*(p+2)/2;
	if (ElemType == 2) (*tmpnDOF) = pow(p+1,2);
	if (ElemType == 3) (*tmpnDOF) = (p+1)*(p+2)*(p+3)/6;
	if (ElemType == 4) (*tmpnDOF) = pow(p+1,3);
	//xf_printf("tmpnDOF = %d \n",(*tmpnDOF));
	return xf_OK;
	
}
//only support p-enrichment 
//have to take account for the initialization error
/******************************************************************/
//   FUNCTION Definition: xf_ErrInd2ElemPos
static int
xf_ErrInd2ElemPos(xf_All *All, xf_Vector *SpongeFlag, xf_Vector *ErrIndicator, int ***pElemPos)
{
   /*
    PURPOSE:
    
    Determines the global position numbers of each element when sorted
    according to ErrIndicator.  Handles parallel.
    
    INPUTS:
    
    All : All structure
    ErrIndicator : real-valued vector of abs-value error estimates per elem
    
    OUTPUTS:
    
    ElemPos[egrp][elem] = global position number of egrp,elem
    
    RETURN: Error code
    
    */
   int ierr, egrp, elem, negrp, sr, k, Order;
   int nelemref, **ElemPos = NULL;
   int nelemtot, *nElem = NULL, *RI;
   int nelemtot_glob;
   real **Indicator = NULL;
   xf_Mesh *Mesh;
   
   Mesh = All->Mesh;
   negrp = Mesh->nElemGroup;
   FILE *ind_file;

   if (FirstAdapt == 1)  {ind_file = fopen("ErrIndicator1.dat","w"); xf_printf("File being written\n");}
   
   // used for sorting more than one value per element
   
   //!!!!!!!error indicator might have different ranks
   //here we start from considering only one
   sr = 1; 
   //sr = ErrIndicator->StateRank;
   
   // create indicator vector over all elems (for sorting -> fixed fraction)
   ierr = xf_Error(xf_GetnElem(Mesh, &nElem, &nelemtot));
   if (ierr != xf_OK) return ierr;
   
   ierr = xf_Error(xf_VAlloc2((void ***) &Indicator, negrp, nElem, sr*sizeof(real)));
   if (ierr != xf_OK) return ierr;
   
   // copy over indicator data
   for (egrp=0; egrp<negrp; egrp++)
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
      {
         //find the order of primal state
         //Order = xf_InterpOrder(State, egrp, elem);
         if(SpongeFlag == NULL){
         for (k=0; k<sr; k++){
            //if(Order >= MaxOrder)
            //   Indicator[egrp][sr*elem+k] = 0.0; //if cell reaches the maxorder; have change to the next one in list
            //else
               Indicator[egrp][sr*elem+k] = ErrIndicator->GenArray[egrp].rValue[elem][k];
           	   if (FirstAdapt == 1) fprintf(ind_file, "%d %d %g\n", egrp, elem, Indicator[egrp][sr*elem+k]);
           	   //xf_printf("%.10lf %.10lf \n", ErrIndicator->GenArray[egrp].rValue[elem][k], Indicator[egrp][sr*elem+k]);
           	}
         }
         else
         {
         for (k=0; k<sr; k++)
            if(SpongeFlag->GenArray[egrp].rValue[elem][0] > 0.) // in sponge
               Indicator[egrp][sr*elem+k] = -1.0e+16;
            else
               Indicator[egrp][sr*elem+k] = ErrIndicator->GenArray[egrp].rValue[elem][k];
         }
      }

   // create vector for storing element position post sorting
   ierr = xf_Error(xf_VAlloc2((void ***) pElemPos, negrp, nElem, sr*sizeof(int)));
   if (ierr != xf_OK) return ierr;
   
   // sort indicator (ascending)
   ierr = xf_Error(xf_SortRealParallel(Indicator[0], sr*nelemtot, xfe_False, (*pElemPos)[0]));
   if (ierr != xf_OK) return ierr;
   
   xf_Release( (void *) nElem);
   xf_Release2( (void **) Indicator);

   if (FirstAdapt == 1) fclose(ind_file);
   
   return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_AdaptSmoothRef
static int
xf_AdaptSmoothRef(xf_All *All, xf_Vector *RefIndicator)
{
   /*
    PURPOSE:
    
    Smoothes adaptation to eliminate islands, fill voids, etc.
    
    INPUTS:
    
    All : All structure
    RefIndicator : coming in, elements to be refined have this set to 1;
    all other elements have this set to 0
    
    OUTPUTS:
    
    RefIndicator : modified refinement indicator, based on smoothing operations
    
    RETURN:
    
    Error Code
    */
   int ierr, egrp, elem, face;
   int count, nneigh, nface, negrp;
   int egN, eN, faceN, iiface, hang;
   int *RI;
   real frac;
   const real frac_island = 0.7;
   const real frac_void   = 0.7;
   xf_Mesh *Mesh;
   int p, tmpDOF1, tmpDOF2;
   
   Mesh = All->Mesh;
   
   xf_printf("Smoothing refinement indicator.\n");
   
   negrp = Mesh->nElemGroup;
   nIsland = 0; nVoid = 0;
   
   // remove islands
   for (egrp=0; egrp<negrp; egrp++)
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
         RI = RefIndicator->GenArray[egrp].iValue[elem];
         if (RI[0] == 0) continue;
         nface = Mesh->ElemGroup[egrp].nFace[elem];
         count  = 0;
         nneigh = 0;
         for (face=0; face<nface; face++){
            ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, face, &egN, &eN, &faceN));
            if (ierr != xf_OK) return ierr;
            if (egN >= 0){
               nneigh++;
               // increment count if face is not hanging and other side is not refined
               iiface = Mesh->ElemGroup[egrp].Face[elem][face].Number;
               if ((Mesh->IFace[iiface].HangNumber            == 0) &&
                   (RefIndicator->GenArray[egN].iValue[eN][0] == 0)) count++;
            }
         }
         frac = ((real) count) / ((real) nneigh);
         if (frac > frac_island){
         	//xf_printf("Island removed\n");
	 		p = VOrder->GenArray[egrp].iValue[elem][0]; 
	 		// previously assigned change in order
		    ierr = xf_Error(ec_GetDOF(p+RI[0], &tmpDOF1));
		    if (ierr != xf_OK) return ierr;
		    // propose no change (remove island)
		    ierr = xf_Error(ec_GetDOF(p, &tmpDOF2));
		    if (ierr != xf_OK) return ierr;
		    // check if max DOF would be exceeded - if not proceed
		    if ((totalDOF + tmpDOF2 - tmpDOF1) > MaxDOF) continue;
			else{
	         	nIsland++;
	          	RI[0] = 0;
	          	// update DOF
	          	totalDOF = totalDOF + tmpDOF2 - tmpDOF1;
	        }
          }
      }
   
   
   // fill voids
   for (egrp=0; egrp<negrp; egrp++)
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
         RI = RefIndicator->GenArray[egrp].iValue[elem];
         if (RI[0] > 0) continue;
         nface = Mesh->ElemGroup[egrp].nFace[elem];
         count  = 0;
         nneigh = 0;
         for (face=0; face<nface; face++){
            ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, face, &egN, &eN, &faceN));
            if (ierr != xf_OK) return ierr;
            if (egN >= 0){
               nneigh++;
               // check if on coarse side of a hanging face
               ierr = xf_Error(xf_CheckHangFace(Mesh, egrp, elem, face, &hang, NULL, NULL, NULL));
               if (ierr != xf_OK) return ierr;
               if ((hang) || (RefIndicator->GenArray[egN].iValue[eN][0]==1)) count++;
            }
         }
         if (nneigh == 0) continue;
         frac = ((real) count) / ((real) nneigh);
         if (frac > frac_void){
         	// xf_printf("Void filled\n");
	 		p = VOrder->GenArray[egrp].iValue[elem][0]; 
	 		// previously assigned order
		    ierr = xf_Error(ec_GetDOF(p+RI[0], &tmpDOF1));
		    if (ierr != xf_OK) return ierr;
		    // proposed change (fill void)
		    ierr = xf_Error(ec_GetDOF(p+1, &tmpDOF2));
		    if (ierr != xf_OK) return ierr;
		    // check if max DOF or max order would be exceeded - if not proceed
		    if ((totalDOF + tmpDOF2 - tmpDOF1 > MaxDOF) || (p+1 > MaxOrder)) continue;
			else{
	         	nVoid++;
	         	RI[0] = 1;
	         	// update DOF
	         	totalDOF = totalDOF + tmpDOF2 - tmpDOF1;
	        }
         }
      }
   
   return xf_OK;
}


/**************************************************************************/
// support dynamic p-adaptation
enum xfe_Bool WriteVOrder = xfe_True; //if output the cell-wise order specification 
static int
Yu_AdaptOrderRef(xf_All *All, xf_Vector *SpongeFlag, const char *SavePrefix, xf_Vector *RefIndicator, 
                 xf_Vector *BasisOrderElem)
{
   int ierr, i, j;
   int Order, ref, ChgOrder;
   char OutputFile[xf_MAXSTRLEN];
   enum xfe_Bool Interpolated, ParallelFlag, NeedRefine;
   xf_Vector *V;
   //xf_Vector *VOrder;
   xf_DataSet *DataSet;
   xf_Data *D;
   xf_Mesh *Mesh;
   int p, tmpDOF, egrp, elem, negrp;

   // locate a vector for specifying desired order
   if (FirstAdapt == 1){
	   ierr = xf_Error(xf_FindVector(All, "VOrder", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
	                                 NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, xfe_False,
	                                 xfe_False, NULL, &VOrder, NULL));
	   if (ierr != xf_OK) return ierr;
   }

   Mesh = All->Mesh;
   DataSet = All->DataSet;
   D = DataSet->Head;

   //check if running in parallel mode
   if (Mesh->ParallelInfo != NULL) ParallelFlag = xfe_True;
   else ParallelFlag = xfe_False;

   //pass refinement indicator
   if (ParallelFlag){
       ierr = xf_Error(xf_HaloExchangeVectorBegin(RefIndicator));
       if (ierr != xf_OK) return ierr;
       ierr = xf_Error(xf_HaloExchangeVectorEnd(RefIndicator));
       if (ierr != xf_OK) return ierr;
      // ierr = xf_Error(xf_HaloExchangeVectorBegin(SpongeFlag));
      // if (ierr != xf_OK) return ierr;
      // ierr = xf_Error(xf_HaloExchangeVectorEnd(SpongeFlag));
      // if (ierr != xf_OK) return ierr;
   }

   // project existing vectors based on RefIndicator + basedline order specify
   while (D != NULL) {
      if (D->Type == xfe_Vector){
         V = (xf_Vector *) D->Data;
         //Interpolated = ((V->Basis != NULL) && (V->Order != NULL));

         //!!!!for debug
         Interpolated = xfe_True;

         // on project primal state
         if ((V->Linkage == xfe_LinkageGlobElem) && (Interpolated) &&
             (V->SolverRole == xfe_SolverRolePrimalState)){

            NeedRefine = xfe_False;

            //set VOrder
            for (i=0; i<VOrder->nArray; i++)
               for (j=0; j<VOrder->GenArray[i].n; j++){
                  //current order
                  Order = xf_InterpOrder(V, i, j);

                  // check
                  if (FirstAdapt == 0){
                  	if (VOrder->GenArray[i].iValue[j][0] != Order) xf_printf("MISMATCHING ORDER!!!!\n");
                  }

                  if(SpongeFlag!=NULL && SpongeFlag->GenArray[i].rValue[j][0]>0.)
                  {
                     //in sponge
                     if(Order != SpongeOrder)
                        NeedRefine = xfe_True;

                     VOrder->GenArray[i].iValue[j][0] = SpongeOrder;

                  }
                  
                  // if sponge is not used; or this cell not in sponge
                  if(SpongeFlag==NULL || (SpongeFlag!=NULL && SpongeFlag->GenArray[i].rValue[j][0]<0.)){
                     //refine or coarsen ?
                     //ChgOrder= max(MinOrder,Order + RefIndicator->GenArray[i].iValue[j][0]);
                     //ChgOrder= min(ChgOrder, MaxOrder);
                  	 ChgOrder = Order + RefIndicator->GenArray[i].iValue[j][0];
                     VOrder->GenArray[i].iValue[j][0] = ChgOrder; 

                     if (RefIndicator->GenArray[i].iValue[j][0] != 0) 
                        NeedRefine = xfe_True;
                  }
               }//for;for

            if(ParallelFlag){
               ierr = xf_Error(xf_MPI_Allreduce(&NeedRefine, 1, xfe_SizeInt, xfe_MPI_MAX));
               if (ierr != xf_OK) return ierr;
            }

            if(NeedRefine)//only do projection for cells requesting order change
            {
               //project V in place
               ierr = xf_Error(xf_ProjectVectorInPlace_VOrder(All->Mesh, All->DataSet, V, V->Basis,
                                                              xfe_BasisLast, VOrder));
               if (ierr != xf_OK) return ierr;

               //if requested, write VOrder into file
               if(WriteVOrder){
                  sprintf(OutputFile, "%s_VOrder0.data\0", SavePrefix);
                  ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "VOrder", VOrder, OutputFile));
                  if (ierr != xf_OK) return ierr;
               }
            }//if

         }//if
      }
      D = D->Next;
   }

   // Get totalDOF through different way (to check)
   totDOFchk = 0;
   Mesh = All->Mesh;
   negrp = Mesh->nElemGroup;
   for (egrp=0; egrp<negrp; egrp++){
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
	 		p = VOrder->GenArray[egrp].iValue[elem][0];
			ierr = xf_Error(ec_GetDOF(p, &tmpDOF));
			if (ierr != xf_OK) return ierr;
			totDOFchk = totDOFchk + tmpDOF;
		}
	}


   //destroy VOrder
   // ierr = xf_Error(xf_DestroyVector(VOrder, xfe_True));
   // if (ierr != xf_OK) return ierr;
   // if (ntest == ncomp){
	  //  sprintf(OutputFile, "%s_VOrder1.data\0", "test");
	  //  ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "VOrder", VOrder, OutputFile));
	  //  if (ierr != xf_OK) return ierr;
   // }

   return xf_OK;
}
/*****************************************************************************************/
/*
static int
Yu_AdaptInd2RefInd_OLD(xf_All *All, xf_Vector *Primal_State, xf_Vector *ErrIndicator, xf_Vector *RefIndicator,
                   int *pnelemref, int *pnelemcoarse)
{
   int ierr, egrp, elem, negrp, i;
   int nelemref, nelemcoarse, **ElemPos = NULL, *RI;
   int nelemtot, nelemtot_glob;
   xf_Mesh *Mesh;
   
   Mesh = All->Mesh;
   negrp = Mesh->nElemGroup;
   
   // create indicator vector over all elems (for sorting -> fixed fraction)
   ierr = xf_Error(xf_GetnElem(Mesh, NULL, &nelemtot));
   if (ierr != xf_OK) return ierr;
   
   // Obtain ElemPos[egrp][elem] = global position number
   // key function call to obtain the rank of element for refinement
   ierr = xf_Error(xf_ErrInd2ElemPos(All, ErrIndicator, &ElemPos));
   if (ierr != xf_OK) return ierr;
   
   // nelemtot_glob = global number of elements
   nelemtot_glob = nelemtot;
   ierr = xf_Error(xf_MPI_Allreduce(&nelemtot_glob, 1, xfe_SizeInt, xfe_MPI_SUM));
   if (ierr != xf_OK) return ierr;
   
   // number of elements to refine
   nelemref = nelemtot_glob*refinefrac;
   (*pnelemref) = nelemref;
   
   // flag all elements with ElemPos >= nelement_glob - nelemref;
   for (egrp=0; egrp<negrp; egrp++){
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
         RI = RefIndicator->GenArray[egrp].iValue[elem];
         RI[0] = 0;
         
         if (ElemPos[egrp][elem] >= (nelemtot_glob-nelemref))
            RI[0] = 1;
         
      }//elem
   }//egrp
   
   // number of elements to coarsen
   nelemcoarse = nelemtot_glob*coarsenfrac;
   (*pnelemcoarse) = nelemcoarse;
   
   //if coarsen functionality is evoked
   if(nelemcoarse > 0 )  {
      
      // flag all elements with ElemPos <= nelemcoarse
      for (egrp=0; egrp<negrp; egrp++){
         for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
            RI = RefIndicator->GenArray[egrp].iValue[elem];
            
            if (ElemPos[egrp][elem] <= nelemcoarse)
               RI[0] = -1;
            
         }//elem
      }//egrp
      
   }//if
   
   //enable smooth refinement
   //ierr = xf_Error(xf_AdaptSmoothRef(All, RefIndicator));
   //if (ierr != xf_OK) return ierr;
   
   xf_Release2( (void **) ElemPos);
   
   return xf_OK;
}
 */
/*****************************************************************************************/
static int
Yu_AdaptInd2RefInd(xf_All *All, xf_Vector *SpongeFlag, xf_Vector *ErrIndicator,
                   xf_Vector *RefIndicator, int glob_neleminside, int glob_nelemsponge,
                   int nelemref, int nelemcoarse)
{
   int ierr, egrp, elem, negrp;
   int **ElemPos = NULL, *RI;
   int nelemtot, nelemtot_glob;
   xf_Mesh *Mesh;
   int InitDOF, tmpDOF1, tmpDOF2, target, p;
   int *ascend_elem = NULL;
   int *ascend_egrp = NULL;
   int idx, cstart, c, pc, i, j;
   char OutputFile[xf_MAXSTRLEN];
   
   Mesh = All->Mesh;
   negrp = Mesh->nElemGroup;
   nelemtot_glob = glob_nelemsponge + glob_neleminside;

   nref = 0; ncoarse = 0;

   // if (ntest == ncomp+1){
	  //  sprintf(OutputFile, "%s_VOrder2.data\0", "test");
	  //  ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "VOrder", VOrder, OutputFile));
	  //  if (ierr != xf_OK) return ierr;
   // }

   // Obtain ElemPos[egrp][elem] = global position number
   // key function call to obtain the rank of element for refinement
   xf_printf("before xf_ErrInd2ElemPos\n");
   ierr = xf_Error(xf_ErrInd2ElemPos(All, SpongeFlag, ErrIndicator, &ElemPos));
   if (ierr != xf_OK) return ierr;
   xf_printf("after xf_ErrInd2ElemPos\n");

   // Organize ElemPos into ascend_elem and ascend_egrp, which output elem and egrp, respectively
   // They are indexed such that higher indices mean higher indicator values
   // Also set RefIndicator to zero everywhere
   ierr = xf_Error(xf_Alloc((void **) &ascend_elem, nelemtot_glob, sizeof(int)));
   if (ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_Alloc((void **) &ascend_egrp, nelemtot_glob, sizeof(int)));
   if (ierr != xf_OK) return ierr;
   for (egrp=0; egrp<negrp; egrp++){
   	   for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
   	       idx = ElemPos[egrp][elem];
   	       ascend_elem[idx] = elem; ascend_egrp[idx] = egrp;
   	       if (idx == 0) xf_printf("ElemPos = 0\n");
	       xf_printf("%d %d %d\n",egrp,elem,idx);
   	       // set RefIndicator to zero
   	       RefIndicator->GenArray[egrp].iValue[elem][0] = 0;
   	   }
   }

   xf_printf("after organizing ElemPos\n");   	           
   // if first adaptation, make as many elements max order as possible without exceeding max DOF
   if (FirstAdapt == 1)
   {
	   ierr = xf_Error(ec_GetDOF(InitOrder, &tmpDOF1));
	   if (ierr != xf_OK) return ierr;
	   InitDOF = nelemtot_glob*tmpDOF1; // initial number of total DOF
	   xf_printf("total elem glob = %d, tmpDOF1 = %d, initial order = %d, elem type = %d \n",nelemtot_glob,tmpDOF1, InitOrder,ElemType);
	   xf_printf("initial nDOF = %d.\n",InitDOF);
	   ierr = xf_Error(ec_GetDOF(MaxOrder, &tmpDOF2));
	   if (ierr != xf_OK) return ierr;
	   nelemref = (MaxDOF - InitDOF)/(tmpDOF2-tmpDOF1); 
	   nref = nelemref;

	   totalDOF = tmpDOF2*nref + tmpDOF1*(nelemtot_glob-nref);
	   xf_printf("total nDOF at first adapt = %d.\n",totalDOF);
   	   // flag all elements with ElemPos >= nelement_glob - nelemref;
	   for (egrp=0; egrp<negrp; egrp++){
	      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
	         RI = RefIndicator->GenArray[egrp].iValue[elem];
	         RI[0] = 0;
	         
	         if (ElemPos[egrp][elem] >= (nelemtot_glob-nref))
	            	RI[0] = MaxOrder - InitOrder;	         
	      }//elem
	   }//egrp
   }
   else{
   		cstart = 0;
		for (idx = nelemtot_glob-1; idx >= 0; idx--){
			// Unpack 
			egrp = ascend_egrp[idx]; elem = ascend_elem[idx];
	 		p = VOrder->GenArray[egrp].iValue[elem][0];
	        RI = RefIndicator->GenArray[egrp].iValue[elem];
	 		if (p == MaxOrder) continue; // if already at max order, skip
	 		if (RI[0] != 0) continue; // if already selected for coarsening, skip

	 		// check if increasing order would exceed max DOF
	 		ierr = xf_Error(ec_GetDOF(p, &tmpDOF1));
		    if (ierr != xf_OK) return ierr;
		    ierr = xf_Error(ec_GetDOF(p+1, &tmpDOF2));
		    if (ierr != xf_OK) return ierr;

		    if (totalDOF + tmpDOF2 - tmpDOF1 <= MaxDOF){
		    	RI[0] = 1; 
		    	totalDOF = totalDOF + tmpDOF2 - tmpDOF1;
		    	nref++;
		    }
		    else{
		    	// if increasing order by 1 would exceed max DOF, check if we can coarsen starting from low indicator values
		    	for (c = cstart; c < idx; c++){
		    		i = ascend_egrp[c]; j = ascend_elem[c];
	 				pc = VOrder->GenArray[i].iValue[j][0];
	 				if ((pc > p) && (RefIndicator->GenArray[i].iValue[j][0] == 0)){ 
	 				// this polynomial order is greater than above order, and its hasn't been modified
	 					// refine above
	 					RI[0] = 1; 
	 					// coarsen this
				        RefIndicator->GenArray[i].iValue[j][0] = -1;
				        // update totalDOF from the refinement
				        totalDOF = totalDOF + tmpDOF2 - tmpDOF1;
				        // update totalDOF from the coarsening
	 					ierr = xf_Error(ec_GetDOF(pc, &tmpDOF1));
		    			if (ierr != xf_OK) return ierr;
		    			ierr = xf_Error(ec_GetDOF(pc-1, &tmpDOF2));
		    			if (ierr != xf_OK) return ierr;
		    			totalDOF = totalDOF + tmpDOF2 - tmpDOF1;

		    			nref++; ncoarse++;
				        break;
	 				}

		    	}
		    }//if


	    }//idx

   }//if
   
   // //if coarsen functionality is evoked
   // if(nelemcoarse > 0 )  {
      
   //    // flag all elements with ElemPos <= nelemcoarse
   //    for (egrp=0; egrp<negrp; egrp++){
   //       for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
   //          RI = RefIndicator->GenArray[egrp].iValue[elem];
            
   //          if (ElemPos[egrp][elem] <= nelemcoarse + glob_nelemsponge
   //              && ElemPos[egrp][elem]>=glob_nelemsponge)
   //             RI[0] = -1;
            
   //       }//elem
   //    }//egrp
      
   // }//if
   
   //enable smooth refinement
   // ierr = xf_Error(xf_AdaptSmoothRef(All, RefIndicator));
   // if (ierr != xf_OK) return ierr;
   
   xf_Release2( (void **) ElemPos);
   
   return xf_OK;
}
//******************************************************************************************//
//initial preparation for p-adaptation with or without sponge
static int
Yu_DynPadapt_config(xf_All *All, const int sr, Yu_Sponge *SP, xf_Vector *SpongeFlag, int *pglob_ntotelem,
                    int *pglob_nspgelem, int *pnelemref, int *pnelemcoarse)
{
   int ierr, elem, egrp, negrp;
   int dim, Order, nn, i, k, *Node, flag_tot;
   int glob_ntotelem, ntotelem, glob_nspgelem, nspgelem;
   real spg_pt[3][3], point[3], flag[4], tmp, eps=1.0e-6;
   char multistateflag[30];
   enum xfe_BasisType Basis;
   FILE *spg_file;

   xf_Mesh *Mesh;
   
   Mesh = All->Mesh;
   negrp = Mesh->nElemGroup;
   dim = Mesh->Dim;

   if(SpongeFlag != NULL)
   {
      spg_file = fopen("sponge","r");

      for(i=0; i<3; i++)
         for(k=0; k<2; k++)
            fscanf(spg_file, "%lf", &spg_pt[i][k]);

      //read in target state
      for(i=0; i<sr; i++)
         fscanf(spg_file, "%lf", &(SP->State_Inf[i]));
 
      //determine whether to specify multistate in the sponge 
      fscanf(spg_file, "%s", multistateflag);

      SP->WheMultiStateSponge = xfe_False; 
      if(strcmp(multistateflag, "multistate") == 0)
         SP->WheMultiStateSponge = xfe_True;

      fclose(spg_file);
      
      //shift the x-coordinate of first point toward exterior a bit
      spg_pt[0][0] -= eps;
   
      //save sponge info to Yu_Sponge struct
      SP->point[0*3+0]=spg_pt[0][0]; SP->point[0*3+1]=spg_pt[0][1]; 
      SP->point[1*3+0]=spg_pt[1][0]; SP->point[1*3+1]=spg_pt[1][1]; 
      SP->point[2*3+0]=spg_pt[2][0]; SP->point[2*3+1]=spg_pt[2][1]; 
      SP->point[3*3+0]=spg_pt[0][0]; SP->point[3*3+1]=(spg_pt[1][1] + spg_pt[2][1]) - spg_pt[0][1];
     
      //sponge boundary line
      SP->line[0*3+0]= -(spg_pt[1][1]-spg_pt[0][1]);
      SP->line[0*3+1]= (spg_pt[1][0]-spg_pt[0][0]);
      SP->line[0*3+2]= (spg_pt[1][1]-spg_pt[0][1])*spg_pt[0][0] - (spg_pt[1][0]-spg_pt[0][0])*spg_pt[0][1];
      SP->line[1*3+0]= (spg_pt[1][1]-spg_pt[2][1]);
      SP->line[1*3+1]= -(spg_pt[1][0]-spg_pt[2][0]);
      SP->line[1*3+2]= -(spg_pt[1][1]-spg_pt[2][1])*spg_pt[2][0] + (spg_pt[1][0]-spg_pt[2][0])*spg_pt[2][1];
      tmp = (spg_pt[1][1] + spg_pt[2][1]) - spg_pt[0][1];
      SP->line[2*3+0]= -(tmp-spg_pt[2][1]);
      SP->line[2*3+1]= (spg_pt[0][0]-spg_pt[2][0]);
      SP->line[2*3+2]= (tmp-spg_pt[2][1])*spg_pt[2][0] - (spg_pt[0][0]-spg_pt[2][0])*spg_pt[2][1];
      SP->line[3*3+0]= -1.;
      SP->line[3*3+1]= 0.;
      SP->line[3*3+2]= spg_pt[0][0];
   }

   ntotelem = 0;
   nspgelem = 0;
   for (egrp=0; egrp<negrp; egrp++){
      Basis = Mesh->ElemGroup[egrp].QBasis;
      Order = Mesh->ElemGroup[egrp].QOrder;
      
      ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
      if (ierr != xf_OK) return ierr;
      
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
         Node = Mesh->ElemGroup[egrp].Node[elem];
         
         if(SpongeFlag != NULL) //sponge is applied
         {
         //means in sponge
         SpongeFlag->GenArray[egrp].rValue[elem][0] = -1.0;
         for(i=0; i<nn; i++)
         {
            for(k=0; k<dim; k++)
               point[k] = Mesh->Coord[Node[i]][k];
           
            //dim = 3 assume sponge rotating 360 along x-axis
            if(dim == 3)
               point[1] = sign(point[1]) * sqrt(point[1]*point[1] + point[2]*point[2]);

            // < means in
            for(k=0, flag_tot=0; k<4; k++)
            {   
               flag[k] = (SP->line[k*3+0]) * point[0] + (SP->line[k*3+1]) * point[1] + (SP->line[k*3+2]);
               if(flag[k] > eps)
                  flag_tot += pow(10, k);
        
               //for cylinder case; temporary change here
               //flag[k] = sqrt(point[0]*point[0] + point[1]*point[1])- 10.;
               //if(flag[k] > eps)
               //   flag_tot = 10.;
            }

            //means in sponge
            if(flag_tot > 0)
            {
               if(flag_tot > (int) SpongeFlag->GenArray[egrp].rValue[elem][0])
                  SpongeFlag->GenArray[egrp].rValue[elem][0] = (real) flag_tot;
            }
         
            //if inflow sponge is applied; check in or out
            //....
         
         }//nn
         
         //if(SP->WheMultiStateSponge == xfe_False)
         {
                //only consider the element inside domain/out of sponge
                if(SpongeFlag->GenArray[egrp].rValue[elem][0] == -1.0)
                   ntotelem++;
                else
                   nspgelem++;
   
          }
      /*   else
         {
                //only consider the element 
                if(SpongeFlag->GenArray[egrp].rValue[elem][0] != -1.0)
                {          
                   ntotelem++;
                   SpongeFlag->GenArray[egrp].rValue[elem][0] = -1.0;
                }
                else
                {
                   nspgelem++;
                   SpongeFlag->GenArray[egrp].rValue[elem][0] = 1.0;
                }
         }
      */          


         }
         else
         {
            //sponge is not applied
            ntotelem++;
         }
         
      }//elem
   }//egrp
   
   //sum over all processes
   glob_ntotelem = ntotelem;
   ierr = xf_Error(xf_MPI_Allreduce(&glob_ntotelem, 1, xfe_SizeInt, xfe_MPI_SUM));
   if (ierr != xf_OK) return ierr;
   
   glob_nspgelem = nspgelem;
   ierr = xf_Error(xf_MPI_Allreduce(&glob_nspgelem, 1, xfe_SizeInt, xfe_MPI_SUM));
   if (ierr != xf_OK) return ierr;

   (*pglob_ntotelem) = glob_ntotelem;
   (*pglob_nspgelem) = glob_nspgelem;
   //glob eliminate sponge
   (*pnelemref)    = glob_ntotelem * refinefrac;
   (*pnelemcoarse) = glob_ntotelem * coarsenfrac;
   return xf_OK;
}

//******************************************************************************************//
static int
SpongeSource(int nq, int sr, int dim, real *xglob, real *U, real *S, real dt)
{
   int ierr, i, j;
   real *Un, *Sn;
   real rhoInf, rhouInf, rhovInf, rhoEInf, omg;

   //debug !!!
   //
   rhoInf = 1.0;
   rhouInf = 0.5; rhovInf = 0.5;
   rhoEInf = 1.0/(1.4-1.) + 0.5 * (rhouInf*rhouInf + rhovInf*rhovInf)/rhoInf;
   //damping factor 
   omg = sqrt(1.4) * 5. / 10.;
   for(i=0; i<nq; i++)
   {
      Un = U + i*sr;
      Sn = S + i*sr;

      Sn[0] = (1. - exp(-omg*dt)) * (U[0] - rhoInf);
      Sn[1] = (1. - exp(-omg*dt)) * (U[1] - rhouInf);
      Sn[2] = (1. - exp(-omg*dt)) * (U[2] - rhovInf);
      Sn[3] = (1. - exp(-omg*dt)) * (U[3] - rhoEInf);
   }

   return xf_OK;
}


//******************************************************************************************//
Yu_ResetPointers(xf_All *All, Yu_Model *Model, xf_Vector **State, xf_Vector **AdaptInd, xf_Vector **RefInd)
{
   int ierr;
   xf_Data *Target;

   ierr = xf_FindDataByTitle(All->DataSet, "ElemMaxCharSpeed", xfe_Vector, &Target);
   if(ierr == xf_NOT_FOUND)
   {       
      xf_printf("Cannot find maximum characteristic speed...\n");
      return ierr;
   }
   else    
      Model->MaxCharSpeed = (xf_Vector *) Target->Data;

   ierr = xf_FindDataByTitle(All->DataSet, "ElemMinFaceLen", xfe_Vector, &Target);
   if(ierr == xf_NOT_FOUND)
   {       
      xf_printf("Cannot find minimum length of elements...\n");
      return ierr;
   }
   else    
      Model->MinFaceLen = (xf_Vector *) Target->Data;

   //reset pointer for state
   ierr = xf_Error(xf_FindPrimalState(All->DataSet, 0, &Target, NULL));
   if (ierr != xf_OK) return ierr;
   (*State) = (xf_Vector *)Target->Data;

   //reset pointer for adaptive indicator
   ierr = xf_FindDataByTitle(All->DataSet, "AdaptIndicator", xfe_Vector, &Target);
   if(ierr == xf_NOT_FOUND)
   {           
      xf_printf("Cannot find adaptive indicator...\n");
      return ierr;
   }   
   else        
      (*AdaptInd) = (xf_Vector *) Target->Data;
   
   //reset pointer for refinement indicator
   ierr = xf_FindDataByTitle(All->DataSet, "RefIndicator", xfe_Vector, &Target);
   if(ierr == xf_NOT_FOUND)
   {           
      xf_printf("Cannot find refinement indicator...\n");
      return ierr;
   }   
   else        
      (*RefInd) = (xf_Vector *) Target->Data;
   

   if(Model->EntropyBdFlag){
   ierr = xf_FindDataByTitle(All->DataSet, "MaxPressure", xfe_Vector, &Target);
   if(ierr == xf_NOT_FOUND)
   {           
      xf_printf("Cannot find adaptive indicator...\n");
      return ierr;
   }   
   else        
      Model->EntropyVec = (xf_Vector *) Target->Data;
   }

   Yu_MinFaceLengthScale(All, Model, *State);

   //reset cpu balane indicator
   Model->cputime_indicator = 0.;

   return xf_OK;
}

//******************************************************************************************//
//rule of p-adaptation
int
xfYu_ApplyUnsteadyAdapt(xf_All *All, Yu_Model *Model, Yu_Limiter ** Limiter, const char *SavePrefix,
                        enum xfe_Bool RestartFlag, xf_Vector *U0, xf_TimeHistData *TimeHistData)
{
   int ierr, i, j, nelemcoarse, nelemref, glob_neleminside, glob_nelemsponge;
   int total_num_adapt, iadapt, init_write_offset;
   int total_num_write, adapt_sr, xfa_file_write_index, load_balance_multiple;
   real Time, EndTime, num_write_per_slot, xfa_file_write_next_time;
   char OutputFile[xf_MAXSTRLEN];
   FILE *spg_file;
   enum xfe_Bool load_imbalance;
   xf_Vector *AdaptIndicator, *RefIndicator;
   xf_Vector *SpongeFlag=NULL;
   xf_Data *D;
   xf_ElemSearchStruct *ESS = NULL;
   real *xpoint;

   //limiter and AV is temporally not supported
   if(Model->LimiterFlag || Model->AVmodel)
      return xf_CODE_LOGIC_ERROR;

   //read in parameters from Model
   MaxOrder    = Model->Dyn_p_Adapt_param->MaxOrder;
   MinOrder    = Model->Dyn_p_Adapt_param->MinOrder;
   SpongeOrder = Model->Dyn_p_Adapt_param->spongeOrder;
   refinefrac  = Model->Dyn_p_Adapt_param->refinefrac;
   coarsenfrac = Model->Dyn_p_Adapt_param->coarsenfrac;
   Yu_adapt_time_size = Model->Dyn_p_Adapt_param->Yu_adapt_time_size;
   xfa_file_write_inv = Model->Dyn_p_Adapt_param->xfa_file_write_inv;
   Yu_load_balance_time_size = Model->Dyn_p_Adapt_param->Yu_load_balance_time_size;

   //for load balancing
   Model->Dyn_p_Adapt_param->best_time = 1.e+16;
   Model->Dyn_p_Adapt_param->best_imp  = 1.;   
   Model->Dyn_p_Adapt_param->adapt_index = 0;

   if(Yu_adapt_time_size < 0. || xfa_file_write_inv < 0.){
      xf_printf("P-adaptation initialization error!\n");
      return xf_CODE_LOGIC_ERROR;
   }
   xfa_file_write_index=1;
   xfa_file_write_next_time=xfa_file_write_inv;

   //currently sponge is not considered
   SpongeFlag = NULL;
 /*  
   spg_file = fopen("sponge","r");
   if(spg_file != NULL){
      //if a sponge file exists
      fclose(spg_file);
      //determine the sponge flag
      ierr = xf_Error(xf_FindSimilarVector(All, Model->MinFaceLen, "SpongeFlag", xfe_True, xfe_True, NULL, &SpongeFlag, NULL));
      if (ierr != xf_OK) return ierr;
   
   
      //make a record in Model struct
      Model->Sponge = SpongeFlag;
      ierr = xf_Error(xf_Alloc( (void **) &(Model->SpongeParam), 1, sizeof(Yu_Sponge)));
      if (ierr != xf_OK) return ierr;
   
      //if sponge is defiend, check farpBC for state variables
      //currently only support on far-field BC
    //  for(i=0; i<Model->nBCs; i++)
    //     if(Model->typeBCs[i] == farPBC){
    //        for(j=0; j<5; j++)
    //           Model->SpongeParam->State_Inf[j] = Model->paraBCs[i][j];
    //        break;
    //     }
   }
*/

   // initial polynomial order
   InitOrder = Model->order;
   // basis type
   // 1 for tri, 2 for quad, 3 for tet, 4 for hex
   if(strcmp(Model->basis, "TriHierarch")==0)           		 {ElemType = 1; xf_printf("ElemType = %d - Tri.\n",ElemType);}
         else if(strcmp(Model->basis, "QuadLegendre")==0)    {ElemType = 2; xf_printf("ElemType = %d - Quad.\n",ElemType);}
         else if(strcmp(Model->basis, "TetHierarch")==0) 	 {ElemType = 3; xf_printf("ElemType = %d - Tet.\n",ElemType);}
         else if(strcmp(Model->basis, "HexLegendre")==0) 	 {ElemType = 4; xf_printf("ElemType = %d - Hex.\n",ElemType);}
         else 
         { printf("Basis not properly defined (adapt)\n");   return xf_Error(xf_NOT_SUPPORTED); }
   
   xf_printf("Initial polynomial order is %d.\n",InitOrder);
   xf_printf("Maximum DOF allowed is %d.\n", MaxDOF);

   ierr = xf_Error(Yu_DynPadapt_config(All, Model->nVars, Model->SpongeParam, SpongeFlag, &glob_neleminside, &glob_nelemsponge, 
                                       &nelemref, &nelemcoarse));
   if (ierr != xf_OK) return ierr;
      
   if(SpongeFlag != NULL)
   {
      ierr = xf_Error(xf_HaloExchangeVectorBegin(SpongeFlag));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_HaloExchangeVectorEnd(SpongeFlag));
      if (ierr != xf_OK) return ierr;
   }
   
   xf_printf("inside: %d; sponge: %d %d %d\n", glob_neleminside, glob_nelemsponge, nelemref, nelemcoarse);

   //initialization
   adapt_sr = NUM_VAR_DYM_P;
   ierr = xf_Error(Yu_CreateOrFindAdaptIndicator(All, xfe_False, adapt_sr, &AdaptIndicator));
   if (ierr != xf_OK) return ierr;
                         
   ierr = xf_Error(xf_SetZeroVector(AdaptIndicator));
   if (ierr != xf_OK) return ierr;
   
   //call supporting routine in xfYu_Adaptoin.c for creatinga refinement indicator
   ierr = xf_Error(Yu_CreateRefIndicator(All, AdaptIndicator, &RefIndicator, NULL));
   if (ierr != xf_OK) return ierr;
  
   //initial adaptation to incorporate the order variation in sponge 
   //possible include static order specification for other purposes
   ierr = xf_Error(Yu_AdaptOrderRef(All, SpongeFlag, SavePrefix, RefIndicator, NULL));
   if (ierr != xf_OK) return ierr;

   if (xf_NotNull(SavePrefix)){
      sprintf(OutputFile, "%s_init.xfa\0", SavePrefix);
      ierr = xf_Error(xf_WriteAllBinary(All, OutputFile));
      if (ierr!=xf_OK) return ierr;
   }

   //figure out adaptation configuration
   ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &Time));
   if (ierr != xf_OK) return ierr;
      
    
   ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "EndTime", &EndTime));
   if (ierr != xf_OK) return ierr;

   total_num_adapt = (EndTime - Time)/Yu_adapt_time_size;

   if(total_num_adapt == 0){
      Yu_adapt_time_size = EndTime - Time;
      total_num_adapt = 1;
   }

   load_balance_multiple = Yu_load_balance_time_size/Yu_adapt_time_size;

   ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UnsteadyWriteInterval", 
                                     &total_num_write));
   if (ierr != xf_OK) return ierr;

   num_write_per_slot = (real)total_num_write/(real)total_num_adapt;

   ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "WriteOffset", &init_write_offset));
   if (ierr != xf_OK) return ierr;

   //synchronize output time at beginning
   for(i=0; i<Model->nOutput; i++) 
      Model->Output[i].pretime = Time; 
   
   xf_printf("adaptation calculation starts......\n");
   //solve & adaptation on each time slot
   for(iadapt=0; iadapt<total_num_adapt; iadapt++)
   {
      //forward solve
      ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time+iadapt*Yu_adapt_time_size));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "EndTime", Time+(iadapt+1)*Yu_adapt_time_size));
      if (ierr != xf_OK) return ierr;

      if(num_write_per_slot >= 1.0){
      ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "UnsteadyWriteInterval", (int)num_write_per_slot));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "WriteOffset", init_write_offset+iadapt*(int)num_write_per_slot));
      if (ierr != xf_OK) return ierr;
      }
      else
      {
         if((iadapt+1)%(int)(1./num_write_per_slot) == 0){
            ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "UnsteadyWriteInterval", 1));
            if (ierr != xf_OK) return ierr;
         }
         else
         {
            ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "UnsteadyWriteInterval", 0));
            if (ierr != xf_OK) return ierr;
         }

      ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "WriteOffset", init_write_offset+(int)((real)iadapt*num_write_per_slot)));
      if (ierr != xf_OK) return ierr;
      }

      //adaptation indicator will be cleaned up in the solver routine.
      //call solver
      ierr = xf_Error(xfYu_ApplyTimeScheme(All, Model, NULL, SavePrefix, xfe_False, U0, TimeHistData));
      if (ierr != xf_OK) return ierr;



      //load balance at specific interval; skip this step if finished
      if((iadapt+1)%load_balance_multiple == 0 && fabs(Time+(iadapt+1)*Yu_adapt_time_size - EndTime)>MEPS)
      {
         //generate a work load report
         ierr = xf_Error(xf_CPUWorkReport(All, U0, &load_imbalance));
         if (ierr!=xf_OK) return ierr;

         if(load_imbalance){
            xf_printf("cpu load imbalance; now requesting re-partition...\n");

            //delete some unuseful data 
            DestroyEntropyBoundStruct();

            //before repartition; make sure to synchronize the output data
            if(Model->Output != NULL)
            {
               ierr = xf_Error(Yu_OutputPointValueDump(All, &Model->Output[0]));
               if (ierr!=xf_OK) return ierr;
            }


            ierr = xf_Error(xfYu_CPULoadBalance(All, &U0, 1));
            if (ierr!=xf_OK) return ierr;
    
   /*          
             if (ProcViz) {
                 ierr = xf_Error(xf_FindDataByTitle(All->DataSet, "ProcID", xfe_Vector, &D));
                 if (ierr != xf_OK) return ierr;
                 ierr = xf_Error(xf_DataSetAdd(DataSet, D->Title, xfe_Vector,
                                               xfe_True, D->Data, NULL));
                 if (ierr != xf_OK) return ierr;
             }
   */
           
             //reset pointers for vectors in Model struct
             ierr = xf_Error(Yu_ResetPointers(All, Model, &U0, &AdaptIndicator, &RefIndicator));
             if( ierr != xf_OK) return ierr;
            
             //if output struct is specified
             if(Model->Output != NULL)
             {
                //local coordinates require re-initialization
                ierr = xf_Error(xf_Alloc((void **) &ESS, 1, sizeof(xf_ElemSearchStruct)));
                if (ierr != xf_OK) return ierr;
                ierr = xf_Error(xf_Alloc((void **) &xpoint, 3*(Model->Output->nPoints), sizeof(real)));
                if (ierr != xf_OK) return ierr;

                ierr = xf_Error(xf_BuildElemSearchStructure(All, ESS));
                if (ierr != xf_OK) return ierr;

                for(i=0; i<Model->Output->nPoints; i++){
               
                   xpoint[i*3 + 0] = *(Model->Output->PointData + i*(3 + Model->Output->nVars) + 0);
                   xpoint[i*3 + 1] = *(Model->Output->PointData + i*(3 + Model->Output->nVars) + 1);
                   xpoint[i*3 + 2] = *(Model->Output->PointData + i*(3 + Model->Output->nVars) + 2);
                
                }
                   
                ierr = xf_Error(xf_FindElemUsingSearchStructure(All, Model->Output->nPoints, xpoint, ESS,
                                 Model->Output->egrpLocal, Model->Output->elemLocal, Model->Output->xref));
                if (ierr != xf_OK) return ierr;

                ierr = xf_Error(xf_DestroyElemSearchStructure(ESS));
                if (ierr != xf_OK) return ierr;
                xf_Release((void *) xpoint);

             }

             xf_printf("done.\n");fflush(stdout);
             xf_printf("load balance finished!\n");
          //  break;
         }

      }//load balance
       
       //skip adaptation step if simulation is finished
       if((fabs(Time+(iadapt+1)*Yu_adapt_time_size - EndTime) > MEPS))
       {
       	   ntest++;
       	   //p-adaptation
           //steps: 0-update refinement indicator; 1-refine the mesh using order increase, all in parallel
           //ierr = xf_Error(Yu_AdaptInd2RefInd_OLD(All, U0, AdaptIndicator, RefIndicator, &nelemref, &nelemcoarse));
           //if (ierr != xf_OK) return ierr;
           xf_printf("before Yu_AdaptInd2RefInd\n");
           ierr = xf_Error(Yu_AdaptInd2RefInd(All, SpongeFlag, AdaptIndicator, RefIndicator, glob_neleminside, glob_nelemsponge,
                                              nelemref, nelemcoarse));
           if (ierr != xf_OK) return ierr;
           
           xf_printf("execute p-adaptation......\n");
           xf_printf("%d elements have order reduction.\n", ncoarse);
           xf_printf("%d elements have order increase.\n", nref);
           xf_printf("%d islands removed.\n",nIsland);
           xf_printf("%d voids filled.\n",nVoid);
           xf_printf("Current # DOF is %d.\n",totalDOF);
           
           // check if max DOF exceeded
           if (totalDOF > MaxDOF) {xf_printf("WARNING - MAX DOF EXCEEDED\n"); return xf_Error(xf_NOT_SUPPORTED);}

           ierr = xf_Error(Yu_AdaptOrderRef(All, SpongeFlag, SavePrefix, RefIndicator, NULL));
           if (ierr != xf_OK) return ierr;

           xf_printf("Alternate # DOF is %d.\n",totDOFchk);
           
           if(xfa_file_write_next_time <= Time+(iadapt+1)*Yu_adapt_time_size){
               //write a .xfa file after each adaptation (contains .msh info)
               if (xf_NotNull(SavePrefix)){
                   sprintf(OutputFile, "%s_B%02d.xfa\0", SavePrefix, xfa_file_write_index);
                   ierr = xf_Error(xf_WriteAllBinary(All, OutputFile));
                   if (ierr!=xf_OK) return ierr;
                   
                   xfa_file_write_index++;
                   xfa_file_write_next_time += xfa_file_write_inv;
               }
           }
           FirstAdapt = 0; // no longer first adaptation
       }
       else if (iadapt == total_num_adapt-1)
       {
       	   int egrp, elem, negrp, k, sr = 1;
           xf_Mesh *Mesh;
		   Mesh = All->Mesh;
  	  	   negrp = Mesh->nElemGroup;
  	  	   FILE *ind_file;

   	       ind_file = fopen("ErrIndicator2.dat","w");
           // write indicator values to file - end of simulation
       	   xf_printf("Writing final indicator values to file\n");
       	   for (egrp=0; egrp<negrp; egrp++)
           for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
      	   {
      	   	   for (k=0; k<sr; k++){
	           	   fprintf(ind_file, "%d %d %g\n", egrp, elem, AdaptIndicator->GenArray[egrp].rValue[elem][k]);
           	   }
      	   }
      	   fclose(ind_file);
       }
   }

   //destroy unuseful data
   ierr = xf_Error(xf_DataSetRemove(All->DataSet, "MaxPressure", xfe_False));
   if (ierr!=xf_OK) return ierr;
   ierr = xf_Error(xf_DataSetRemove(All->DataSet, "RefIndicator", xfe_False));
   if (ierr!=xf_OK) return ierr;
   ierr = xf_Error(xf_DataSetRemove(All->DataSet, "ElemMaxCharSpeed", xfe_False));
   if (ierr!=xf_OK) return ierr;

   if (xf_NotNull(SavePrefix)){
      sprintf(OutputFile, "%s_final.xfa\0", SavePrefix);
      ierr = xf_Error(xf_WriteAllBinary(All, OutputFile));
      if (ierr!=xf_OK) return ierr;
   }

   //destry unuseful data 
   //destroy RefIndicator
   //ierr = xf_Error(xf_DestroyVector(RefIndicator, xfe_True));
   //if (ierr != xf_OK) return ierr;

   return xf_OK;
}
