//here implementation of Yu's own adaptation algorithm
//The following routines support the static H-refinement capability
//
//necessary declaration
static int xf_PreAdaptRobust(xf_All *All, xf_KeyValue *pKeyValueOrig);
static int xf_PostAdaptRobust(xf_All *All, xf_KeyValue *KeyValueOrig);
static int xf_FindAdaptIndicator(xf_All *All, enum xfe_Bool WithRefOpt,
                                 xf_Vector **pAdaptIndicator);
static int xf_PreAdaptHang(xf_All **pAll, xf_Vector **pAdaptIndicator,
                           xf_Vector **pRefIndicator, int **OrderToRefine,
                           int nelemref);
static int xf_PostAdaptHang(xf_All *All, xf_Vector *RefIndicator, 
                            int **OrderToRefine);

//allocate or find the adaptation indicator
static int  
Yu_FindAdaptIndicator(xf_All *All, enum xfe_Bool WithRefOpt, int numrank,  
                            xf_Vector **pAdaptIndicator)
{
   int ierr, egrp, nopt;
   int *rvec = NULL;
   enum xfe_ShapeType Shape;
   xf_Data *D;
   xf_Mesh *Mesh;
               
   Mesh = All->Mesh;
                 
   if (WithRefOpt){
   
      ierr = xf_Error(xf_Alloc( (void **) &rvec, Mesh->nElemGroup, sizeof(int)));
      if (ierr != xf_OK) return ierr;
                                
      
      for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
         ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
         if (ierr != xf_OK) return ierr;
         
         ierr = xf_Error(xf_GetNumRefOpt(Shape, &nopt));
         if (ierr != xf_OK) return ierr;
         
         rvec[egrp] = nopt; // includes no refinement as an option
      } // egrp
   }
                   
   ierr = xf_Error(xf_FindVector(All, "AdaptIndicator", xfe_LinkageGlobElem, numrank, NULL, 0, 0, 
                                 NULL, NULL, NULL, NULL, rvec, xfe_SizeReal, xfe_True,  xfe_True, &D,  
                                 pAdaptIndicator, NULL));
   if (ierr != xf_OK) return ierr;
   
   D->ReadWrite = xfe_True; // make indicator writeable
                         
   xf_Release( (void *) rvec);
                           
   return xf_OK;
}

int
Yu_CreateOrFindAdaptIndicator(xf_All *All, enum xfe_Bool WithRefOpt, int numrank,
                        xf_Vector **pAdaptIndicator)
{
   int ierr;
   ierr = xf_Error(Yu_FindAdaptIndicator(All, WithRefOpt, numrank, pAdaptIndicator));
   if (ierr != xf_OK) return ierr;
   
   return xf_OK;
}


int
Yu_CreateRefIndicator(xf_All *All, xf_Vector *AdaptIndicator,
                      xf_Vector **pRefIndicator, int *nelemref)
{
   int ierr, negrp, egrp, elem, *RI, nelemtot, i;
   real tmp;
   xf_Vector *RefIndicator;
   xf_Data *D; 
   xf_Mesh *Mesh;

   Mesh = All->Mesh;
   negrp = Mesh->nElemGroup;

   //create or find indicator vector
   ierr = xf_Error(xf_FindVector(All, "RefIndicator", xfe_LinkageGlobElem,
                                 1, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL,
                                 xfe_SizeInt, xfe_True,  xfe_True, &D,
                                 pRefIndicator, NULL));
   RefIndicator = (*pRefIndicator);

   D->ReadWrite = xfe_True; // make indicator writeable
   
   //for Yu's purpose
   //simply set errIndicator = refIndicator
   ierr = xf_Error(xf_SetZeroVector(RefIndicator));
   if (ierr != xf_OK) return ierr;
   
   for (egrp=0; egrp<negrp; egrp++)
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
         RI = RefIndicator->GenArray[egrp].iValue[elem];
         tmp = 1.e+12 * (AdaptIndicator->GenArray[egrp].rValue[elem][0]);
         RI[0] = (int) tmp;
         printf("%lf\n", RI[0]);
      }

   return xf_OK;
}

int
Yu_MeshAdaptation(xf_All *All, Yu_Model *Model)
{

   int ierr, sr, egrp, elem, n, nn, s;
   int iAdapt, nelemtot, myRank, nProc;
   char filename[xf_MAXSTRLEN], SavePrefix[xf_MAXSTRLEN], OutputFile[xf_MAXSTRLEN];
   xf_KeyValue KeyValueOrig;
   xf_Vector *dt, *Gp, *dU, *R, *EG;
   xf_Vector *AdaptIndicator, *RefIndicator;
   xf_JacobianMatrix *R_U;
   xf_Data *D;
   
   if(Model->Stat_h_Adapt){
      //backup before conduct adaptation
      ierr = xf_Error(xf_PreAdaptRobust(All, &KeyValueOrig));
      if (ierr != xf_OK) return ierr;
      
      xf_printf("Now, Let's start adapting mesh!~\n");
      
      //find adaptation indicator
      ierr = xf_Error(xf_FindAdaptIndicator(All, xfe_False, &AdaptIndicator));
      if (ierr != xf_OK) return ierr;
      AdaptIndicator->SolverRole = xfe_SolverRoleOther;
      
      //perform adaptation (call xf_AdaptAll)
      /* If SavePrefix is None or NULL, will not write anything */
      ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
      if (ierr != xf_OK) return ierr;
      
      if (xf_NotNull(SavePrefix)){
         /* Write out .xfa file before performing adaptation */
         iAdapt = 1;
         sprintf(OutputFile, "%s_A%02d.xfa\0", SavePrefix, iAdapt);
         ierr = xf_Error(xf_WriteAllBinary(All, OutputFile));
         if (ierr!=xf_OK) return ierr;
      }
      
      // adaptation on mesh starts here
      ierr = xf_Error(xf_GetnElem(All->Mesh, NULL, &nelemtot));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_MPI_Allreduce(&nelemtot, 1, xfe_SizeInt, xfe_MPI_SUM));
      if (ierr != xf_OK) return ierr;
      xf_printf(" Current number of elements = %d.\n", nelemtot);
      
      xf_printf(" Adapting the mesh.\n");
      ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
      if (ierr != xf_OK) return ierr;
      
      // adaptation procedure
      /* according to XFlow
       0-Create refinement indicator
       1-Seralize the mesh and data
       2-Calculate the direction of refinement
       3-Clean-up error estimation vectors
       4-Separate h and p refinement indicators
       5-Refine mesh/order using hanging-node/vorder
       6-Parallelize mesh and data
       */
      //step 0
      ierr = xf_Error(Yu_CreateRefIndicator(All, AdaptIndicator,
                                            &RefIndicator, NULL));
      if (ierr != xf_OK) return ierr;
      
      //step 1
      ierr = xf_Error(xf_PreAdaptHang(&All, &AdaptIndicator,
                                      &RefIndicator, NULL, 0));
      if (ierr != xf_OK) return ierr;
      
      //step 2 (have done since assume isotrapic refinement
      
      
      if (myRank == 0){
         //step 3
         ierr = xf_Error(xf_ErrEstOutput(All, NULL, NULL, xfe_False,
                                         xfe_True, NULL, xfe_False, NULL));
         if (ierr != xf_OK) return ierr;
         
         //step 4 (since not consider P-refinement)
         
         //step 5
         ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
         if (ierr != xf_OK) return ierr;
      } //myRank = 0
      xf_MPI_Barrier();
      //6-Parallelize the mesh and data.
      ierr = xf_Error(xf_PostAdaptHang(All, RefIndicator, NULL));
      if (ierr != xf_OK) return ierr;
      
      //information on new meshes
      ierr = xf_Error(xf_GetnElem(All->Mesh, NULL, &nelemtot));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_MPI_Allreduce(&nelemtot, 1, xfe_SizeInt,
                                       xfe_MPI_SUM));
      if (ierr != xf_OK) return ierr;
      xf_printf(" New number of elements = %d.\n", nelemtot);
      
      //write out a new .xfa file after performing adaptation
      if (xf_NotNull(SavePrefix)){
         sprintf(OutputFile, "%s_B%02d.xfa\0", SavePrefix, iAdapt+1);
         ierr = xf_Error(xf_WriteAllBinary(All, OutputFile));
         if (ierr!=xf_OK) return ierr;
      }
      
      //copy back the values to All structure
      ierr = xf_Error(xf_PostAdaptRobust(All, &KeyValueOrig));
      if (ierr != xf_OK) return ierr;
      
   }
   
   
   return xf_OK;
}

/************************************************************************/
//purpose: project solution to P0 on elements which will be adapt
//doing this imposes entropy bounded the elements and robustness enhance
Yu_SolutionProcessBeforeAdaptation(xf_All *All, Yu_Model *Model, xf_Vector *U)
{
   int ierr, i, j, k, iq;
   int elem, egrp, sr, Order, pOrder, QuadOrder;
   int nq, nn, dim, nface, endIndx, bgnIndx;
   int Neigh_egrp, Neigh_elem, Neigh_face;
   real *EU, *QuadStat;
   real *xq, Velem, FaceJ;
   enum xfe_ShapeType Shape;
   enum xfe_Bool QuadChanged, flag;
   enum xfe_BasisType Basis;
   xf_QuadData *QuadDataElem, *QuadDataFace;
   xf_JacobianData *JData;
   xf_BasisData *PhiData;
   xf_Mesh *Mesh;
   xf_Vector *AdaptIndicator;
   real Umean[50], AdaptFlag;
   
   //////////////////////
   Mesh = All->Mesh;
   sr   = Model->nVars;
   dim  = Mesh->Dim;
   Basis = U->Basis[0];
   
   //////////////////////
   QuadDataFace = NULL;
   QuadDataElem = NULL;
   JData        = NULL;
   PhiData      = NULL;

   //find the adaptation indicator
   ierr = xf_Error(Yu_CreateOrFindAdaptIndicator(All, xfe_False, 1, &AdaptIndicator));
   if (ierr != xf_OK) return ierr;
   
   //loopping over all elements
   for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
         
         EU = U->GenArray[egrp].rValue[elem]; //state vector on element
         Order = xf_InterpOrder(U, egrp, elem);
         AdaptFlag = AdaptIndicator->GenArray[egrp].rValue[elem][0];
         
         //quadrature points inside element
         ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, Order, &QuadOrder));
         if (ierr != xf_OK) return ierr;
         
         ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadDataElem, &QuadChanged));
         if (ierr != xf_OK) return ierr;
         
         //element volume Jacobian information
         QuadChanged = xfe_False;
         ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, QuadDataElem->nquad,
                                         QuadDataElem->xquad, xfb_detJ, QuadChanged, &JData));
         if (ierr != xf_OK) return ierr;
         
         //compute the element volume
         Velem = 0.0;
         for(k=0; k<QuadDataElem->nquad; k++)
            Velem += QuadDataElem->wquad[k]*JData->detJ[k*(JData->nq!=1)];
         
         //allocate memory for state
         ierr = xf_Error(xf_Alloc((void **) &QuadStat, QuadDataElem->nquad*sr, sizeof(real)));
         if (ierr != xf_OK) return ierr;
         
         //basis coefficients for all quadrature points 
         ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, QuadDataElem->nquad, 
                                      QuadDataElem->xquad, xfb_Phi, &PhiData));
         if (ierr != xf_OK) return ierr;
                           
         nn = PhiData->nn;

         //evaluate state at quadrature points
         xf_MxM_Set(PhiData->Phi, EU, QuadDataElem->nquad, nn, sr, QuadStat);
         
         //compute mean quantities
         for (j=0; j<sr; j++)
            Umean[j] = 0.;
         for(k=0; k<QuadDataElem->nquad; k++)
            for (j=0; j<sr; j++)
               Umean[j] += QuadStat[k*sr+j] * QuadDataElem->wquad[k]*JData->detJ[k*(JData->nq!=1)];
         for (j=0; j<sr; j++)
            Umean[j] /= Velem;
         
         //if element will be adapted; projection to P0
         if(AdaptFlag > 1.0e-14)
         {
             switch(Basis) {
               case xfe_QuadLegendre:
               case xfe_HexLegendre:
                  for(j=0; j<sr; j++)
                     EU[0*sr+j] = Umean[j];
                  for(i=1; i<nn; i++)
                     for(j=0; j<sr; j++)
                        EU[i*sr+j] *= 0.0;
                  break;
                  
               case xfe_TriHierarch:
               case xfe_TetHierarch:
                  for(i=0; i<dim+1; i++)
                     for(j=0; j<sr; j++)
                        EU[i*sr+j] = Umean[j];
                  for(i=dim+1; i<nn; i++)
                     for(j=0; j<sr; j++)
                        EU[i*sr+j] *= 0.0;
                  break;
                  
               default:
                  return xf_Error(xf_UNKNOWN_BASIS);
                  break;
            }
         } 
         
         //renew all the allocation
         if(QuadDataElem != NULL){
            ierr = xf_Error(xf_DestroyGenericQuadData(QuadDataElem));
            if (ierr != xf_OK) return ierr;
         }
         
         if(QuadDataFace != NULL){
            ierr = xf_Error(xf_DestroyGenericQuadData(QuadDataFace));
            if (ierr != xf_OK) return ierr;
         }
         
         if(JData != NULL){
            ierr = xf_Error(xf_DestroyJacobianData(JData));
            if (ierr != xf_OK) return ierr;
         }

         if (PhiData != NULL){
            ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
            if (ierr != xf_OK) return ierr;
         }
         
         xf_Release( (void *) QuadStat);
         
         QuadDataFace = NULL;
         QuadDataElem = NULL;
         JData        = NULL;
         QuadStat     = NULL;
         PhiData = NULL;
         
      }//elem
   }//egrp
   return xf_OK;
}

//following is the routines for supporting p-enrichment capability 
//this p-enrichment is designed for unsteady problem 

