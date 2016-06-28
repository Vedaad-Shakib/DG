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
 FILE:  xf_Partition.c
 
 This file contains functions for mesh partitioning.
 
 */


#include "xf_AllStruct.h"
#include "xf_MPIStruct.h"
#include "xf_MPI.h"
#include "xf_Memory.h"
#include "xf_MeshTools.h"
#include "xf_Residual.h"
#include "xf_LinearSolverStruct.h"
#include "xf_Line.h"
#include "xf_Param.h"
#include "xf_Solver.h"
#include "xf_Data.h"
#include "xf_Mesh.h"
#include "xf_All.h"
#include "xf_Math.h"
#include <parmetis.h>
#include <mpi.h>

#if PARMETIS_MAJOR_VERSION > 3
typedef idx_t idxtype;
#endif


/******************************************************************/
//   FUNCTION Definition: xf_PartitionMesh
int 
xf_PartitionMesh( xf_Mesh *Mesh, int ***pElem2Proc, int **ElemWeight, 
                 int *ConnectWeight)
{
  
  int ierr, nProc, negrp, i;
  int egrp, elem, face;
  int iadj, nface, egN, eN;
  int wgtflag, numflag, ncon, edgecut, nparts;
  int nelemtot, options[4];
  int *eglob, *nElem, *vwgt, *adjwgt;
  idxtype *vtxdist, *xadj, *adjncy, *part;
  float *tpwgts, *ubvec;
  MPI_Comm comm;
  
  // obtain nProc = # domains
  ierr = xf_Error(xf_MPI_GetRank(NULL, &nProc));
  if (ierr != xf_OK) return ierr;
  
  
  /* Create nelemtot and eglob[egrp] = running count of global elem
   number. Also create nElem = # elems in each elem group. */
  negrp = Mesh->nElemGroup;
  ierr = xf_Error(xf_Alloc( (void **) &eglob, negrp, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &nElem, negrp, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  for (egrp=0, nelemtot=0; egrp<Mesh->nElemGroup; egrp++){
    eglob[egrp] = nelemtot;
    nelemtot += (nElem[egrp] = Mesh->ElemGroup[egrp].nElem);
  }
  
  /*
   vtxdist: This array describes how the vertices of the graph are
   distributed among the processors. (See discussion in
   Section 4.1). Its contents are identical for every
   processor.
   */
  
  ierr = xf_Error(xf_Alloc( (void **) &vtxdist, 2, sizeof(idxtype)));
  if (ierr != xf_OK) return ierr;
  
  vtxdist[0] = 0;
  vtxdist[1] = nelemtot;
  
  /*
   xadj, adjncy: These store the (local) adjacency structure of the
   graph at each processor. (See discussion in Section 4.1).
   */
  ierr = xf_Error(xf_Alloc( (void **) &xadj, nelemtot+1, sizeof(idxtype)) );
  if (ierr != xf_OK) return ierr;
  
  /*
   vwgt, adjwgt These store the weights of the vertices and
   edges. (See discussion in Section 4.1).
   
   wgtflag: This is used to indicate if the graph is
   weighted. wgtflag can take one of four values:
   
   0   No weights (vwgt and adjwgt are both NULL).
   1   Weights on the edges only (vwgt is NULL).
   2   Weights on the vertices only (adjwgt is NULL).
   3   Weights on both the vertices and edges.
   */
  
  vwgt = NULL;
  adjwgt = NULL;
  
  wgtflag = 0;
  if (ElemWeight != NULL){
    ierr = xf_Error(xf_Alloc( (void **) &vwgt, nelemtot, sizeof(idxtype)) );
    if (ierr != xf_OK) return ierr;
    wgtflag += 2;
  }
  
  ierr = xf_Error(xf_Alloc( (void **) &adjncy, 2*Mesh->nIFace, sizeof(idxtype)) );
  if (ierr != xf_OK) return ierr;
  
  if (ConnectWeight != NULL) {
    ierr = xf_Error(xf_Alloc( (void **) &adjwgt, 2*Mesh->nIFace, sizeof(idxtype)) );
    if (ierr != xf_OK) return ierr;
    wgtflag += 1;
  }
  
  xadj[0] = 0;
  iadj = 0;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      nface = Mesh->ElemGroup[egrp].nFace[elem];
      for (face=0; face<nface; face++){
        ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, face, 
                                              &egN, &eN, NULL));
        if (ierr != xf_OK) return ierr;
        if (egN < 0) continue;
        adjncy[iadj] = eglob[egN] + eN;
        if (adjwgt != NULL)
          adjwgt[iadj] = ConnectWeight[Mesh->ElemGroup[egrp].Face[elem][face].Number];
        iadj++;
      }
      xadj[eglob[egrp]+elem+1] = iadj;
      if (vwgt != NULL)
        vwgt[eglob[egrp]+elem] = ElemWeight[egrp][elem];
    } // elem
  } // egrp
  
  
  /*
   numflag: This is used to indicate the numbering scheme that is
   used for the vtxdist, xadj, adjncy, and part
   arrays. numflag can take one of two values: 0 C-style
   numbering that starts from 0.  1 Fortran-style
   numbering that starts from 1.
   */
  numflag = 0;
  
  
  /*
   ncon: This is used to specify the number of weights that each
   vertex has. It is also the number of balance constraints
   that must be satisfied.
   */
  ncon = 1;
  
  
  /*
   nparts: This is used to specify the number of sub-domains that are
   desired. Note that the number of sub- domains is
   independent of the number of processors that call this
   routine.
   */
  nparts = nProc;
  
  
  /*
   tpwgts: An array of size ncon x nparts that is used to specify the
   fraction of vertex weight that should be distributed to
   each sub-domain for each balance constraint. If all of the
   sub-domains are to be of the same size for every vertex
   weight, then each of the ncon x nparts elements should be
   set to a value of 1/nparts. If ncon is greater than one,
   the target sub-domain weights for each sub-domain are
   stored contiguously (similar to the vwgt array). Note that
   the sum of all of the tpwgts for a give vertex weight
   should be one.
   */
  ierr = xf_Error(xf_Alloc( (void *) &tpwgts, ncon*nparts, sizeof(float)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < ncon*nparts; i++) tpwgts[i] = 1.0/nparts;
  
  
  /*
   ubvec: An array of size ncon that is used to specify the imbalance
   tolerance for each vertex weight, with 1 being perfect
   balance and nparts being perfect imbalance. A value of 1.05
   for each of the ncon weights is recommended.
   */
  ierr = xf_Error(xf_Alloc( (void *) &ubvec, ncon, sizeof(float)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < ncon; i++) ubvec[i] = 1.05;
  
  
  /*
   options: This is an array of integers that is used to pass
   additional parameters for the routine. If options[0]=0,
   then the default values are used. If options[0]=1, then
   the remaining two elements of options are interpreted as
   follows: 
   
   options[1]: This specifies the level of information to be
   returned during the execution of the algo- rithm. Timing
   information can be obtained by setting this to
   1. Additional options for this parameter can be obtained
   by looking at the the file defs.h in the ParMETIS- Lib
   directory. The numerical values there should be added to
   obtain the correct value.  The default value is 0.
   
   options[2]: This is the random number seed for the
   routine. The default value is 15.  edgecut Upon
   successful completion, the number of edges that are cut
   by the partitioning is written to this parameter.
   */  
  options[0] = 0; options[1] = 0; options[2] = 0; options[3] = 0;
  
  
  /*
   edgecut: Upon successful completion, the number of edges that are cut by
   the partitioning is written to this parameter.
   */
  edgecut = 0;
  
  
  /*
   part: This is an array of size equal to the number of
   locally-stored vertices. Upon successful completion the
   partition vector of the locally-stored vertices is written
   to this array. (See discussion in Section 4.4).  
   */
  if (sizeof(idxtype) != sizeof(int)) return xf_Error(xf_NOT_SUPPORTED);
  ierr = xf_Error(xf_VAlloc2((void *) pElem2Proc, negrp, nElem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  part = (idxtype *) (*pElem2Proc)[0];
  
  
  /*
   comm: This is a pointer to the MPI communicator of the processes
   that call ParMETIS.
   */
  comm = MPI_COMM_SELF;
  
  
  /*
   From ParMetis 3.1 Manual:
   
   ParMETIS V3 PartKway (idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, 
   idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, 
   int *ncon, int *nparts, float *tpwgts, float *ubvec, 
   int *options, int *edgecut, idxtype *part, MPI Comm *comm)
   
   Description:
   
   This routine is used to compute a k-way partitioning of a graph on p
   processors using the multilevel k-way multi-constraint partitioning
   algorithm.
   
   */
  ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, 
                       &ncon, &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
  
  xf_Release( (void *)   eglob);
  xf_Release( (void *)   nElem);
  xf_Release( (void *)    xadj);
  xf_Release( (void *)  adjncy);
  xf_Release( (void *)    vwgt);
  xf_Release( (void *)  adjwgt);  
  xf_Release( (void *)  tpwgts);
  xf_Release( (void *)   ubvec);
  xf_Release( (void *) vtxdist);
  
  
  return xf_OK;
}

/******************************************************************/
// provide element weight for mesh partition
real
Yu_ElementComputWeight(int Order)
{
    real tmp;
    
    switch (Order) {
        case 0:
            tmp = 1.;
            break;
            
        case 1:
            tmp = 1.9;
            break;
            
        case 2:
            tmp = 5.0;
            break;
            
        case 3:
            tmp = 9.8;
            break;
            
        case 4:
            tmp = 36.6;
            break;
            
        default:
            tmp = 100.;
            break;
    }
    
    return tmp;
}

real
Yu_ElementCommunWeight(int Order)
{
    real tmp;
    
    switch (Order) {
        case 0:
            tmp = 1.;
            break;
            
        case 1:
            tmp = 4.;
            break;
            
        case 2:
            tmp = 9.;
            break;
            
        case 3:
            tmp = 16.;
            break;
            
        case 4:
            tmp = 25.;
            break;
            
        default:
            tmp = 100.;
            break;
    }
    
    return tmp;
}

/******************************************************************/
//   FUNCTION Definition: xf_CPUWorkReport
/*
int
xf_CPUWorkReport(xf_All *All, xf_Vector *State, enum xfe_Bool *load_imbalance)
{
  int ierr, egrp, elem, face, egrpN, elemN, nElemTot, nNonZeros, nDoF, r, rN;
  int nHaloTot, nHaloNonZeros, nHaloDoF, myRank, nProc, perr;
  real tmp, cost_max, cost_min, cost_avg;
  char SavePrefix[xf_MAXSTRLEN], filename[xf_MAXSTRLEN];
  enum xfe_Bool VariableOrder;
  FILE *report;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", 
                                 SavePrefix));
  if (ierr != xf_OK) return ierr;
  
  sprintf(filename,"%s_%dof%d.m", SavePrefix, myRank, nProc);
  
  perr = xf_OK;
  if ((report = fopen(filename, "w")) == NULL) perr = xf_FILE_WRITE_ERROR;
  //Note: using MPI_MAX because error codes are positive
  ierr = xf_Error(xf_MPI_Allreduce(&perr, 1, xfe_SizeInt, xfe_MPI_MAX));
  if (ierr != xf_OK) return ierr;
  if (perr != xf_OK) return perr;
  
  VariableOrder = ((State->nComp != NULL) && (State->vOrder != NULL));
  
  nElemTot = nNonZeros = nDoF = 0;
  for (egrp = 0; egrp < All->Mesh->nElemGroup; egrp++) {
    nElemTot += All->Mesh->ElemGroup[egrp].nElem;
    for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++){
      if (VariableOrder)
        r = State->GenArray[egrp].vr[elem];
      else 
        r = State->GenArray[egrp].r;
      
      nDoF += r;
      nNonZeros += r*r;
      for (face = 0; face < All->Mesh->ElemGroup[egrp].nFace[elem]; face++){
        ierr = xf_Error(xf_NeighborAcrossFace(All->Mesh, egrp, elem, face, 
                                              &egrpN, &elemN, NULL));
        if (ierr != xf_OK) return ierr;
        if (egrpN>=0){
          if (VariableOrder)
            rN = State->GenArray[egrpN].vr[elemN];
          else 
            rN = State->GenArray[egrpN].r;
          nNonZeros += r*rN;
        }
      }
    }
  }
  //get Halo elements work
  if (All->Mesh->ParallelInfo != NULL){
    nHaloDoF = nHaloNonZeros = nHaloTot = 0;
    for (egrp = All->Mesh->nElemGroup; egrp < 2*All->Mesh->nElemGroup; egrp++){
      nHaloTot += All->Mesh->ElemGroup[egrp].nElem;
      for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++){
        if (VariableOrder)
          r = State->GenArray[egrp].vr[elem];
        else 
          r = State->GenArray[egrp].r;
        nHaloDoF += r;
        nHaloNonZeros += r*r;
      }
    }
  }
  //fprintf(report, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
  //fprintf(report, "nElemTot(%d) = %d;\n",myRank+1,nElemTot);
  //fprintf(report, "nDoF(%d) = %d;\n",myRank+1,nDoF);
  //fprintf(report, "nNonZeros(%d) = %d;\n",myRank+1,nNonZeros);
  //fprintf(report, "nHaloTot(%d) = %d;\n",myRank+1,nHaloTot);
  //fprintf(report, "nHaloDof(%d) = %d;\n",myRank+1,nHaloDoF);
  //fprintf(report, "nHaloNonZeros(%d) = %d;\n",myRank+1,nHaloNonZeros);

  //Yu addition cputime
  //fprintf(report, "cpu_cost(%d) = %lf;\n", myRank+1, All->Model->cputime_indicator);
  fprintf(report, "%lf %d %d %d %d", All->Model->cputime_indicator, nElemTot, nDoF, nHaloTot, nHaloDoF);

  fclose(report);

  //blocking for syncronization 
  perr = xf_OK; 
  ierr = xf_Error(xf_MPI_Allreduce(&perr, 1, xfe_SizeInt, xfe_MPI_MAX));
  if (ierr != xf_OK) return ierr;

  (*load_imbalance) = xfe_False; 
  if(myRank==0)
  {
     //use rank 0 for analysize the load balancing data
     cost_max = -1.e+16;
     cost_min = 1.e+16;
     cost_avg = 0.;
     for(r=0; r<nProc; r++)
     {
        sprintf(filename,"%s_%dof%d.m", SavePrefix, r, nProc);
        report = fopen(filename, "r");
        fscanf(report, "%lf", &tmp);
        //printf("cpu time of rank%d: %lf\n", r, tmp);
        if(tmp > cost_max) cost_max = tmp;
        if(tmp < cost_min) cost_min = tmp;
        cost_avg += tmp;
        fclose(report);
     }

     //load imbalance criterion
     //heuritic (1): the max is twice as the average
     //test
     if(cost_max > 2.0 * cost_avg / (real) nProc) (*load_imbalance) = xfe_True; 
  }

  return xf_OK;
}
*/

/******************************************************************/
// create a simplifed copy of the above routine
//   FUNCTION Definition: xf_CPUWorkReport
int
xf_CPUWorkReport(xf_All *All, xf_Vector *State, enum xfe_Bool *load_imbalance)
{
    int ierr, egrp, elem, face, egrpN, elemN, nElemTot, nNonZeros, nDoF, r, rN;
    int nHaloTot, nHaloNonZeros, nHaloDoF, myRank, nProc, perr, Order;
    real tmp, cost_weight, cost_max, cost_avg, cpu_max, cpu_avg;
    char SavePrefix[xf_MAXSTRLEN], filename[xf_MAXSTRLEN];
    enum xfe_Bool VariableOrder;
    FILE *report;
    
    ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix",
                                   SavePrefix));
    if (ierr != xf_OK) return ierr;
    
  //  sprintf(filename,"%s_%dof%d.m", SavePrefix, myRank, nProc);
    
  //  perr = xf_OK;
  //  if ((report = fopen(filename, "w")) == NULL) perr = xf_FILE_WRITE_ERROR;
 
    //Note: using MPI_MAX because error codes are positive
  //  ierr = xf_Error(xf_MPI_Allreduce(&perr, 1, xfe_SizeInt, xfe_MPI_MAX));
  //  if (ierr != xf_OK) return ierr;
  //  if (perr != xf_OK) return perr;
    
    VariableOrder = ((State->nComp != NULL) && (State->vOrder != NULL));
    
    nElemTot = nNonZeros = nDoF = 0;
    cost_weight = 0.;
    for (egrp = 0; egrp < All->Mesh->nElemGroup; egrp++) {
        nElemTot += All->Mesh->ElemGroup[egrp].nElem;
        for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++){
            if (VariableOrder)
                r = State->GenArray[egrp].vr[elem];
            else
                r = State->GenArray[egrp].r;
            
            nDoF += r;
            
            Order = xf_InterpOrder(State, egrp, elem);
            cost_weight += Yu_ElementComputWeight(Order);
         //   nNonZeros += r*r;
         //   for (face = 0; face < All->Mesh->ElemGroup[egrp].nFace[elem]; face++){
         //       ierr = xf_Error(xf_NeighborAcrossFace(All->Mesh, egrp, elem, face,
         //                                             &egrpN, &elemN, NULL));
         //       if (ierr != xf_OK) return ierr;
         //       if (egrpN>=0){
         //           if (VariableOrder)
         //               rN = State->GenArray[egrpN].vr[elemN];
         //           else
         //               rN = State->GenArray[egrpN].r;
         //           nNonZeros += r*rN;
         //       }
         //   }
            
        }
    }
    //get Halo elements work
    /*
    if (All->Mesh->ParallelInfo != NULL){
        nHaloDoF = nHaloNonZeros = nHaloTot = 0;
        for (egrp = All->Mesh->nElemGroup; egrp < 2*All->Mesh->nElemGroup; egrp++){
            nHaloTot += All->Mesh->ElemGroup[egrp].nElem;
            for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++){
                if (VariableOrder)
                    r = State->GenArray[egrp].vr[elem];
                else
                    r = State->GenArray[egrp].r;
                nHaloDoF += r;
                nHaloNonZeros += r*r;
            }
        }
    }
     */
    //fprintf(report, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
    //fprintf(report, "nElemTot(%d) = %d;\n",myRank+1,nElemTot);
    //fprintf(report, "nDoF(%d) = %d;\n",myRank+1,nDoF);
    //fprintf(report, "nNonZeros(%d) = %d;\n",myRank+1,nNonZeros);
    //fprintf(report, "nHaloTot(%d) = %d;\n",myRank+1,nHaloTot);
    //fprintf(report, "nHaloDof(%d) = %d;\n",myRank+1,nHaloDoF);
    //fprintf(report, "nHaloNonZeros(%d) = %d;\n",myRank+1,nHaloNonZeros);
    
    //Yu addition cputime
    //fprintf(report, "cpu_cost(%d) = %lf;\n", myRank+1, All->Model->cputime_indicator);
    //fprintf(report, "%lf %d %d %d %d", All->Model->cputime_indicator, nElemTot, nDoF, nHaloTot, nHaloDoF);
    
    //fclose(report);
    
    //blocking for syncronization
    ierr = xf_Error(xf_MPI_Reduce(&(All->Model->cputime_indicator), &cpu_max, 1, xfe_SizeReal,
                                     xfe_MPI_MAX, 0));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_MPI_Reduce(&(All->Model->cputime_indicator), &cpu_avg, 1, xfe_SizeReal,
                                     xfe_MPI_SUM, 0));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_MPI_Reduce(&cost_weight, &cost_max, 1, xfe_SizeReal,
                                     xfe_MPI_MAX, 0));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_MPI_Reduce(&cost_weight, &cost_avg, 1, xfe_SizeReal,
                                     xfe_MPI_SUM, 0));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_MPI_Reduce(&nDoF, &nNonZeros, 1, xfe_SizeInt,
                                     xfe_MPI_SUM, 0));
    if (ierr != xf_OK) return ierr;
    
    (*load_imbalance) = xfe_False;
    if(myRank==0)
    {
        //write log file
        report = fopen("./cpu_load.log", "a");
        //total DoF
        fprintf(report, "%d  ", nNonZeros);

        xf_printf("total run time = %lf; ", cpu_max);
        fprintf(report, "%lf  ", cpu_max);

        //compute load imbalance
        cpu_avg /= (real)nProc;
        tmp = (real) nProc;
        tmp = (cpu_max - cpu_avg) / cpu_max * tmp / (tmp - 1.);
        xf_printf("total imbalance percentage = %lf; ", tmp);
        fprintf(report, "%lf  ", tmp);
        cpu_avg = tmp;

        cost_avg /= (real)nProc;
        tmp = (real) nProc;
        tmp = (cost_max - cost_avg) / cost_max * tmp / (tmp - 1.);
        xf_printf("load imbalance percentage = %lf\n", tmp);
      
        fprintf(report, "%lf", tmp);

        //load ballancing algorithm
        if(All->Model->Dyn_p_Adapt_param->adapt_index < 10) 
           (*load_imbalance) = xfe_True;   //alway repartition at first ten adaptation
        else
        {
           cpu_max /= (real)nNonZeros;

           if(tmp > 1.2*(All->Model->Dyn_p_Adapt_param->best_imp) &&
              cpu_max > 1.2*(All->Model->Dyn_p_Adapt_param->best_time))
              (*load_imbalance) = xfe_True;
        }

        if((*load_imbalance) == xfe_True) 
           fprintf(report, "  1\n");  //indicate repartition
        else
           fprintf(report, "  0\n");

        //update records;
        if(tmp < All->Model->Dyn_p_Adapt_param->best_imp)
           All->Model->Dyn_p_Adapt_param->best_imp = tmp;
           
        if(cpu_max < All->Model->Dyn_p_Adapt_param->best_time)
           All->Model->Dyn_p_Adapt_param->best_time = cpu_max; 

        All->Model->Dyn_p_Adapt_param->adapt_index++;

        fclose(report);

        //use rank 0 for analysize the load balancing data
        /*
        cost_max = -1.e+16;
        cost_min = 1.e+16;
        cost_avg = 0.;
        for(r=0; r<nProc; r++)
        {
            sprintf(filename,"%s_%dof%d.m", SavePrefix, r, nProc);
            report = fopen(filename, "r");
            fscanf(report, "%lf", &tmp);
            //printf("cpu time of rank%d: %lf\n", r, tmp);
            if(tmp > cost_max) cost_max = tmp;
            if(tmp < cost_min) cost_min = tmp;
            cost_avg += tmp;
            fclose(report);
        }
        
        //load imbalance criterion
        //heuritic (1): the max is twice as the average
        //test
        if(cost_max > 2.0 * cost_avg / (real) nProc) (*load_imbalance) = xfe_True; 
         */
        
    }
   
     //need to broadbast the load-imbalance signal
     ierr = xf_Error(xf_MPI_Bcast((void *) load_imbalance, sizeof(int), 0));
     if (ierr != xf_OK) return ierr;

    return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CPULoadBalance
int
xf_CPULoadBalance(xf_All *All, xf_Vector **pState, int nState)
{
  int ierr, myRank, iLine, ie, egrp, elem, face, **ElemWeight, nelemtot, i;
  int *nElem, Order, Dim, iface, nProc, *ConnectWeight, egrpN, elemN, OrderN;
  enum xfe_PreconditionerType Preconditioner;
  enum xfe_Bool CRequired, SortLines, DeleteNonEssential, ProcViz, *WriteState;
  xf_SolverData *SolverData;
  xf_Vector *R, *C = NULL, *C_Glob = NULL, *State;
  xf_LineSet *LineSet_Glob;
  xf_Mesh *Mesh, *Mesh_Glob;
  xf_DataSet *DataSet_Glob;
  xf_Data *D;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (nState < 1) return xf_Error(xf_INPUT_ERROR);
  
  State = pState[0];
	
	ierr = xf_Error(xf_Alloc((void **)&WriteState, nState, sizeof(int)));
	if (ierr != xf_OK) return ierr;
    
  // locate preconditioner
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Preconditioner", 
                                     xfe_PreconditionerName, 
                                     (int ) xfe_PreconditionerLast, 
                                     (int *) &Preconditioner));
  if (ierr != xf_OK) return ierr;
  // delete non essential data?
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "DD_DeleteNonEssential", 
                                     &DeleteNonEssential));
  if (ierr != xf_OK) return ierr;
  if (DeleteNonEssential){
    ierr = xf_Error(xf_DataSetDeleteNonEssential(All->DataSet));
    if (ierr != xf_OK) return ierr;
  }
  
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "ProcViz", 
                                     &ProcViz));
  if (ierr != xf_OK) return ierr;
  
  // does the preconditioner use lines?
  xf_PreconditionerLineCheck(Preconditioner, &CRequired, &SortLines);
  
  if (CRequired){//we don't need this if no coupling lines are used
    ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_FindSimilarVector(All, State, "Residual", xfe_False, 
                                         xfe_True, NULL, &R, NULL));
    if (ierr != xf_OK) return ierr;
    
    //Calculate a residual so the element coupling and the lines are created
    ierr = xf_Error(xf_CalculateResidual(All, State, R, NULL, SolverData));
    if (ierr != xf_OK) return ierr;
    
    //store pointer to coupling vector of each processor
    C = SolverData->C;
    SolverData->C = NULL;
    
    ierr = xf_Error(xf_DestroySolverData(SolverData));
    if (ierr != xf_OK) return ierr;
  }
  Mesh = All->Mesh;
  
  if (myRank == 0) {
    ierr = xf_Error(xf_CreateMesh(&Mesh_Glob));
    if (ierr != xf_OK) return ierr;
    
    if (CRequired){
          ierr = xf_Error(xf_CreateVector(&C_Glob));
          if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_CreateDataSet(&DataSet_Glob));
    if (ierr != xf_OK) return ierr;
  }
  else
    C_Glob = NULL;
   
  //Unparallelize Mesh
  ierr = xf_Error(xf_UnParallelizeMesh(Mesh, Mesh_Glob));
  if (ierr != xf_OK) return ierr;
  
  if (CRequired){
    //Unparallelize coupling vector
    ierr = xf_Error(xf_UnParallelizeVector(Mesh, C, C_Glob));
    if (ierr != xf_OK) return ierr;
  }
  //making sure all states get unparallelized
  for (i = 0; i < nState; i++){
    ierr = xf_Error(xf_FindPrimalState(All->DataSet, i, &D, NULL));
    if (ierr != xf_OK) return ierr;
		WriteState[i] = D->ReadWrite;
    D->ReadWrite = xfe_True;
  }
  
  //Unparallelize dataset
  ierr = xf_Error(xf_UnParallelizeDataSet(All, All->DataSet, 
                                          DataSet_Glob));
  if (ierr != xf_OK) return ierr;
  
  //destroy old parallel dataset
  ierr = xf_Error(xf_DestroyDataSet(All->DataSet));
  if (ierr != xf_OK) return ierr;
  //destroy old parallel mesh
  ierr = xf_Error(xf_DestroyMesh(All->Mesh));
  if (ierr != xf_OK) return ierr;
  
  ElemWeight = NULL;
  ConnectWeight = NULL;
  
  All->Mesh = NULL;
  All->DataSet = NULL;
  
  if (myRank == 0) {
    All->Mesh = Mesh_Glob;
    All->DataSet = DataSet_Glob;
    Dim = Mesh_Glob->Dim;
    //reset the pointer to the serial mesh
    ierr = xf_Error(xf_FindPrimalState(DataSet_Glob, 0, &D, NULL));
    if (ierr != xf_OK) return ierr;
    State = (xf_Vector *) D->Data;
    
    //Elemental weights
    ierr = xf_Error(xf_GetnElem(Mesh_Glob, &nElem, &nelemtot));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_VAlloc2((void***)&ElemWeight, Mesh_Glob->nElemGroup, nElem, 
                               sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    for (egrp = 0; egrp < Mesh_Glob->nElemGroup; egrp++){
      for (elem = 0; elem < Mesh_Glob->ElemGroup[egrp].nElem; elem++) {
        Order = xf_InterpOrder(State, egrp, elem);
        ElemWeight[egrp][elem] = (int)pow(Order+1,2*Dim);
      }
    }
    
    if (CRequired){
      ierr = xf_Error(xf_CreateLines(All, C_Glob, NULL, &LineSet_Glob));
      if (ierr != xf_OK) return ierr;
      
      //Elemental connectivity weights
      ierr = xf_Error(xf_Alloc((void**)&ConnectWeight, Mesh_Glob->nIFace, 
                               sizeof(int)));
      if (ierr != xf_OK) return ierr;
      for (iface = 0; iface < Mesh_Glob->nIFace; iface++){
        egrp = Mesh_Glob->IFace[iface].ElemGroupL;
        elem = Mesh_Glob->IFace[iface].ElemL;
        egrpN = Mesh_Glob->IFace[iface].ElemGroupR;
        elemN = Mesh_Glob->IFace[iface].ElemR;
        Order = xf_InterpOrder(State, egrp, elem);
        OrderN = xf_InterpOrder(State, egrpN, elemN);
        ConnectWeight[iface] = (int)pow((Order+1),Dim)+(int)pow((OrderN+1),Dim);
      }
      //loop through lines and fill the weights
      for (iLine = 0; iLine < LineSet_Glob->nLine; iLine++) {
        for (ie = 0; ie < LineSet_Glob->Line[iLine].nelem; ie++) {
          elem = LineSet_Glob->Line[iLine].elem[ie];
          egrp = LineSet_Glob->Line[iLine].egrp[ie];
          face = LineSet_Glob->Line[iLine].face[ie];
          Order = xf_InterpOrder(State, egrp, elem);
          if (face>=0){//not the end of the Line
            ierr = xf_Error(xf_NeighborAcrossFace(Mesh_Glob, egrp, elem, face, 
                                                  &egrpN, &elemN, NULL));
            if (ierr != xf_OK) return ierr;
            OrderN = xf_InterpOrder(State, egrpN, elemN);
            //it can only be an internal face
            iface = Mesh_Glob->ElemGroup[egrp].Face[elem][face].Number;
            ConnectWeight[iface] *= max(Mesh_Glob->ElemGroup[egrp].nFace[elem],
                                        Mesh_Glob->ElemGroup[egrpN].nFace[elemN]);
          }
        }
      }
      ierr = xf_Error(xf_DestroyVector(C_Glob, xfe_True));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_DestroyLineSet(LineSet_Glob));
      if (ierr != xf_OK) return ierr;
    }
    else
      ConnectWeight = NULL;
  
    All->Mesh = NULL;
    All->DataSet = NULL;
    xf_Release((void*)nElem);
  }
  //create empty new mesh
  ierr = xf_Error(xf_CreateMesh(&All->Mesh));
  if (ierr != xf_OK) return ierr;
  //create empty new dataset
  ierr = xf_Error(xf_CreateDataSet(&All->DataSet));
  if (ierr != xf_OK) return ierr;

  //Parallelize mesh
  ierr = xf_Error(xf_ParallelizeMesh(Mesh_Glob, All->Mesh, ElemWeight, 
                                     ConnectWeight));
  if (ierr != xf_OK) return ierr;
  
  //Parallelize dataset
  ierr = xf_Error(xf_ParallelizeDataSet(All, DataSet_Glob, All->DataSet));
  if (ierr != xf_OK) return ierr;
  
  if (ProcViz){
    ierr = xf_Error(xf_ProcViz(All));
    if (ierr != xf_OK) return ierr;
  }
  
  //reset the pointers to the state on all procs.
  for (i = 0; i < nState; i++){
    ierr = xf_Error(xf_FindPrimalState(All->DataSet, i, &D, NULL));
    if (ierr != xf_OK) return ierr;
    pState[i] = (xf_Vector *) D->Data;
		D->ReadWrite = WriteState[i];
  }
  //clean-up
  if (myRank == 0) {
    ierr = xf_Error(xf_DestroyMesh(Mesh_Glob));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_DestroyDataSet(DataSet_Glob));
    if (ierr != xf_OK) return ierr;
    
    xf_Release((void*)ConnectWeight);
    xf_Release2((void**)ElemWeight);
  }
  xf_Release((void*)WriteState);
  
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xfYu_CPULoadBalance
//   note: a clean version of the previous function
int
xfYu_CPULoadBalance(xf_All *All, xf_Vector **pState, int nState)
{
  int ierr, myRank, iLine, ie, egrp, elem, face, **ElemWeight, nelemtot, i;
  int *nElem, Order, Dim, iface, nProc, *ConnectWeight, egrpN, elemN, OrderN;
  enum xfe_PreconditionerType Preconditioner;
  enum xfe_Bool CRequired, SortLines, DeleteNonEssential, ProcViz, *WriteState;
  xf_SolverData *SolverData;
  xf_Vector *R, *C = NULL, *C_Glob = NULL, *State;
  xf_LineSet *LineSet_Glob;
  xf_Mesh *Mesh, *Mesh_Glob;
  xf_DataSet *DataSet_Glob;
  xf_Data *D;
 
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (nState < 1) return xf_Error(xf_INPUT_ERROR);
  
  State = pState[0];

  ierr = xf_Error(xf_Alloc((void **)&WriteState, nState, sizeof(int)));
  if (ierr != xf_OK) return ierr;
    
  // locate preconditioner
  //ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Preconditioner", 
  //                                   xfe_PreconditionerName, 
  //                                   (int ) xfe_PreconditionerLast, 
  //                                   (int *) &Preconditioner));
  //if (ierr != xf_OK) return ierr;
  // delete non essential data?
  //ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "DD_DeleteNonEssential", 
  //                                   &DeleteNonEssential));
  //if (ierr != xf_OK) return ierr;
  //if (DeleteNonEssential){
  //  ierr = xf_Error(xf_DataSetDeleteNonEssential(All->DataSet));
  //  if (ierr != xf_OK) return ierr;
  //}
  
  //ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "ProcViz", 
  //                                   &ProcViz));
  //if (ierr != xf_OK) return ierr;
  
  // does the preconditioner use lines?
  //xf_PreconditionerLineCheck(Preconditioner, &CRequired, &SortLines);
  
  /*if (CRequired){//we don't need this if no coupling lines are used
    ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_FindSimilarVector(All, State, "Residual", xfe_False, 
                                         xfe_True, NULL, &R, NULL));
    if (ierr != xf_OK) return ierr;
    
    //Calculate a residual so the element coupling and the lines are created
    ierr = xf_Error(xf_CalculateResidual(All, State, R, NULL, SolverData));
    if (ierr != xf_OK) return ierr;
    
    //store pointer to coupling vector of each processor
    C = SolverData->C;
    SolverData->C = NULL;
    
    ierr = xf_Error(xf_DestroySolverData(SolverData));
    if (ierr != xf_OK) return ierr;
  }
  */
  Mesh = All->Mesh;
  
  if (myRank == 0) {
    ierr = xf_Error(xf_CreateMesh(&Mesh_Glob));
    if (ierr != xf_OK) return ierr;
    
    //if (CRequired){
    //      ierr = xf_Error(xf_CreateVector(&C_Glob));
    //      if (ierr != xf_OK) return ierr;
    //}
    ierr = xf_Error(xf_CreateDataSet(&DataSet_Glob));
    if (ierr != xf_OK) return ierr;
  
  }
  else
    C_Glob = NULL;
   
  //Unparallelize Mesh
  ierr = xf_Error(xf_UnParallelizeMesh(Mesh, Mesh_Glob));
  if (ierr != xf_OK) return ierr;
  
  //if (CRequired){
    //Unparallelize coupling vector
  //  ierr = xf_Error(xf_UnParallelizeVector(Mesh, C, C_Glob));
  //  if (ierr != xf_OK) return ierr;
  //}
  //making sure all states get unparallelized
  for (i = 0; i < nState; i++){
    ierr = xf_Error(xf_FindPrimalState(All->DataSet, i, &D, NULL));
    if (ierr != xf_OK) return ierr;
		
    WriteState[i] = D->ReadWrite;
    D->ReadWrite = xfe_True;
  }
 
  //Unparallelize dataset
  ierr = xf_Error(xf_UnParallelizeDataSet(All, All->DataSet, 
                                          DataSet_Glob));
  if (ierr != xf_OK) return ierr;
  
  //destroy old parallel dataset
  ierr = xf_Error(xf_DestroyDataSet(All->DataSet));
  if (ierr != xf_OK) return ierr;
  //destroy old parallel mesh
  ierr = xf_Error(xf_DestroyMesh(All->Mesh));
  if (ierr != xf_OK) return ierr;
  
  ElemWeight = NULL;
  ConnectWeight = NULL;
  
  All->Mesh = NULL;
  All->DataSet = NULL;

  if (myRank == 0) {
    All->Mesh = Mesh_Glob;
    All->DataSet = DataSet_Glob;
    Dim = Mesh_Glob->Dim;
    //reset the pointer to the serial mesh
    ierr = xf_Error(xf_FindPrimalState(DataSet_Glob, 0, &D, NULL));
    if (ierr != xf_OK) return ierr;
    State = (xf_Vector *) D->Data;
    
    //Elemental weights
    ierr = xf_Error(xf_GetnElem(Mesh_Glob, &nElem, &nelemtot));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_VAlloc2((void***)&ElemWeight, Mesh_Glob->nElemGroup, nElem, 
                               sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    for (egrp = 0; egrp < Mesh_Glob->nElemGroup; egrp++){
      for (elem = 0; elem < Mesh_Glob->ElemGroup[egrp].nElem; elem++) {
        Order = xf_InterpOrder(State, egrp, elem);
        //ElemWeight[egrp][elem] = (int)pow(Order+1,2*Dim);
        ElemWeight[egrp][elem] = (int) Yu_ElementComputWeight(Order);
      }
    }
    
    /*
    if (CRequired){
      ierr = xf_Error(xf_CreateLines(All, C_Glob, NULL, &LineSet_Glob));
      if (ierr != xf_OK) return ierr;
      
      //Elemental connectivity weights
      ierr = xf_Error(xf_Alloc((void**)&ConnectWeight, Mesh_Glob->nIFace, 
                               sizeof(int)));
      if (ierr != xf_OK) return ierr;
      for (iface = 0; iface < Mesh_Glob->nIFace; iface++){
        egrp = Mesh_Glob->IFace[iface].ElemGroupL;
        elem = Mesh_Glob->IFace[iface].ElemL;
        egrpN = Mesh_Glob->IFace[iface].ElemGroupR;
        elemN = Mesh_Glob->IFace[iface].ElemR;
        Order = xf_InterpOrder(State, egrp, elem);
        OrderN = xf_InterpOrder(State, egrpN, elemN);
        ConnectWeight[iface] = (int)pow((Order+1),Dim)+(int)pow((OrderN+1),Dim);
      }
      //loop through lines and fill the weights
      for (iLine = 0; iLine < LineSet_Glob->nLine; iLine++) {
        for (ie = 0; ie < LineSet_Glob->Line[iLine].nelem; ie++) {
          elem = LineSet_Glob->Line[iLine].elem[ie];
          egrp = LineSet_Glob->Line[iLine].egrp[ie];
          face = LineSet_Glob->Line[iLine].face[ie];
          Order = xf_InterpOrder(State, egrp, elem);
          if (face>=0){//not the end of the Line
            ierr = xf_Error(xf_NeighborAcrossFace(Mesh_Glob, egrp, elem, face, 
                                                  &egrpN, &elemN, NULL));
            if (ierr != xf_OK) return ierr;
            OrderN = xf_InterpOrder(State, egrpN, elemN);
            //it can only be an internal face
            iface = Mesh_Glob->ElemGroup[egrp].Face[elem][face].Number;
            ConnectWeight[iface] *= max(Mesh_Glob->ElemGroup[egrp].nFace[elem],
                                        Mesh_Glob->ElemGroup[egrpN].nFace[elemN]);
          }
        }
      }
      ierr = xf_Error(xf_DestroyVector(C_Glob, xfe_True));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_DestroyLineSet(LineSet_Glob));
      if (ierr != xf_OK) return ierr;
    }
    else
      ConnectWeight = NULL;
    */
      
    //Elemental connectivity weights
    ConnectWeight = NULL;
    /*
    ierr = xf_Error(xf_Alloc((void**)&ConnectWeight, Mesh_Glob->nIFace,
                              sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (iface = 0; iface < Mesh_Glob->nIFace; iface++){
        egrp = Mesh_Glob->IFace[iface].ElemGroupL;
        elem = Mesh_Glob->IFace[iface].ElemL;
        egrpN = Mesh_Glob->IFace[iface].ElemGroupR;
        elemN = Mesh_Glob->IFace[iface].ElemR;
        Order = xf_InterpOrder(State, egrp, elem);
        OrderN = xf_InterpOrder(State, egrpN, elemN);
        //ConnectWeight[iface] = (int)pow((Order+1),Dim)+(int)pow((OrderN+1),Dim);
        ConnectWeight[iface] = (int) Yu_ElementCommunWeight(Order) + (int) Yu_ElementCommunWeight(OrderN);
    }
    */
  
    All->Mesh = NULL;
    All->DataSet = NULL;
    xf_Release((void*)nElem);
  }
  
  //create empty new mesh
  ierr = xf_Error(xf_CreateMesh(&All->Mesh));
  if (ierr != xf_OK) return ierr;
  //create empty new dataset
  ierr = xf_Error(xf_CreateDataSet(&All->DataSet));
  if (ierr != xf_OK) return ierr;
  
  //Parallelize mesh
  ierr = xf_Error(xf_ParallelizeMesh(Mesh_Glob, All->Mesh, ElemWeight, 
                                     ConnectWeight));
  if (ierr != xf_OK) return ierr;
  
  //Parallelize dataset
  ierr = xf_Error(xf_ParallelizeDataSet(All, DataSet_Glob, All->DataSet));
  if (ierr != xf_OK) return ierr;
  
  //if (ProcViz){
  //  ierr = xf_Error(xf_ProcViz(All));
  //  if (ierr != xf_OK) return ierr;
  //}
  
  //reset the pointers to the state on all procs.
  for (i = 0; i < nState; i++){
    ierr = xf_Error(xf_FindPrimalState(All->DataSet, i, &D, NULL));
    if (ierr != xf_OK) return ierr;
    pState[i] = (xf_Vector *) D->Data;
		D->ReadWrite = WriteState[i];
  }
  //clean-up
  if (myRank == 0) {
    ierr = xf_Error(xf_DestroyMesh(Mesh_Glob));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_DestroyDataSet(DataSet_Glob));
    if (ierr != xf_OK) return ierr;
    
    xf_Release((void*)ConnectWeight);
    xf_Release2((void**)ElemWeight);
  }
  xf_Release((void*)WriteState);
  
  
  return xf_OK;
}