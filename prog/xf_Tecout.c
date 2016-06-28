//code in developing
//Prequisit: the input cut plane must stay in mesh domain 
//
#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_All.h"
#include "xf_Mesh.h"
#include "xf_EqnSet.h"
#include "xf_MeshTools.h"
#include "xf_Basis.h"
#include "xf_Math.h"
#include "xf_Data.h"
#include "xf_Memory.h"
#include "xf_Quad.h"
#include "xfYu_Model.h"
//#include "reaction_routine.h"
#include <stdlib.h>

//**************************header for input parameter***********************************//
enum TecOut_InputParams{
    TecOut_xfaFileName,
    TecOut_dataFilePrefix,
    TecOut_dataFileBatchBegin,
    TecOut_dataFileBatchInv,
    TecOut_dataFileBatchEnd,
    TecOut_WhtRectDomainOutput,
/*    TecOut_LineDataOrginX,
    TecOut_LineDataOrginY,
    TecOut_LineDataLengthX,
    TecOut_LineDataLengthY,
    TecOut_OutputMeshSize,
 */
    TecOut_WhtCutLineOutput,
    TecOut_WhtPointOutput,
    TecOut_WhtIntergalQuanOutput,
    TecOut_WhtToleVacuum,
    Tecout_Last
};
static char *TecOut_ParamsName[Tecout_Last] = {
    "xfaFileName",
    /* .xfa file for the case required to deal with */
    
    "dataFilePrefix",
    /* the prefix file for the related .data file*/
    
    "dataFileBatchBegin",
    /* the start index of data batch*/
    
    "dataFileBatchInv",
    /* the index interval of data batch*/
    
    "dataFileBatchEnd",
    /* the final index of data batch*/
    
    "WhtRectDomainOutput",
    /* whether do rectangle domain output*/
    
    "WhtCutLineOutput",
    /* whether do line cut output*/
    
    "WhtPointOutput",
    /* whether do point-wise output*/
    
    "WhtIntergalQuanOutput",
    /* whether do integral quantity output*/
    
    "WhtToleVacuum"
    /* whether tolerate empty region in domain output*/
    
};

//***************************************************************************************//
//headers for different output (No. and names of variables);
//need user to specify
static int  DomainNumVars;
static char DomainNameVars[20][xf_MAXSTRLEN];
static int
DomainOutputHeader(FILE *fdomout)
{
    int k, ierr;
    char dummy[xf_MAXSTRLEN];
    
    //user specification
    DomainNumVars = 1;
    strcpy(DomainNameVars[0],"density");
//    strcpy(DomainNameVars[1],"x-vel");
//    strcpy(DomainNameVars[2],"y-vel");
//    strcpy(DomainNameVars[3],"pressure");
    
    
    //Dump data according to Tecplot block format
    sprintf(dummy,"TITLE = DomainOutput.dat\n");
    ierr = xf_fwrite(dummy,sizeof(char),strlen(dummy),fdomout);
    if (ierr!=xf_OK) return ierr;
    
    //Varible name
    sprintf(dummy,"VARIABLES = \"X\", \"Y\",");
    ierr = xf_fwrite(dummy,sizeof(char),strlen(dummy),fdomout);
    if (ierr!=xf_OK) return ierr;
    
//    if(dim == 3){
//        sprintf(dummy," \"Z\",");
//        ierr = xf_fwrite(dummy,sizeof(char),strlen(dummy),fdomout);
//        if (ierr!=xf_OK) return ierr;
//    }
//    if (ierr != xf_OK) return ierr;
    
    //Dump user specified file header
    for(k=0; k<DomainNumVars; k++){
        sprintf(dummy, "\"%s\",", DomainNameVars[k]);
        ierr = xf_fwrite(dummy,sizeof(char),strlen(dummy),fdomout);
        if (ierr!=xf_OK) return ierr;
    }
    
    return xf_OK;
}
static int LineNumVars;
static char LineNameVars[20][xf_MAXSTRLEN];
static int
LineOutputHeader()
{
    return xf_OK;
}
static int PointNumVars;
static char PointNameVars[20][xf_MAXSTRLEN];
static int
PointOutputHeader()
{
    return xf_OK;
}
static int IntegralNumVars;
static char IntegralNameVars[20][xf_MAXSTRLEN];
static int
IntegralOuputHeader()
{
    return xf_OK;
}
//***************************************************************************************//
//define the variables for different outputs; all the quantities are point-wise.
static int
DomainOutputDataDump(Yu_Model *Model, const real *Uin, real *Uout)
{
    return xf_OK;
}
static int
LineOutputDataDump(Yu_Model *Model, const real *Uin, real *Uout)
{
    return xf_OK;
}
static int
PointOutputDataDump(Yu_Model *Model, const real *Uin, real *Uout)
{
    return xf_OK;
}
static int
IntegralOutputDataDump(Yu_Model *Model, const real *Uin, real *Uout)
{
    return xf_OK;
}
//***************************************************************************************/

//linked list data structure for cutplane
typedef struct tagcutelem
{
   int egrp;
   int elem;
   real x[4];
   real refx[4];
   //point to next
   struct tagcutelem *next;
}cutelem;
//declaration of static function in this .c file
static int
WhtEdgeIntersectRect(const real *x0, const real *x1, const real xorg, const real yorg, 
                     const real xlen, const real ylen, enum xfe_Bool *Intersect);

static int
WhtPointInPolygon(const int npt, const real *xnode, const real x,  const real y, enum xfe_Bool *Inside);

//need input here from user
static int
DumpOutputVariable(Yu_Model *Model, const real *UG, real *Uo, const real Maxp, const real gamma)
{ 
   int i, j, k, Nspe, dim;
   real d, u[3], p, T, meanW, Yn;
   real rate[11], Y[11];
/*
   dim = Model->dim;
   Nspe = sr - dim - 2 + 1;
   //for test
   //for(i=0; i<Nout; i++)
   //{
   Yn = 0.0;
   meanW = 0.0;
   for(i=0; i<Nspe-1; i++)
   {
      Y[i] = UG[i+dim+2]/UG[0];
      Yn += UG[i+dim+2]/UG[0];
      meanW += UG[i+dim+2]/UG[0]/Model->moleW[i]; 
   }
   Yn = 1. - Yn;
   meanW += Yn/Model->moleW[Nspe-1];
   Y[Nspe-1] = Yn;

   //Temperature
   p = 101325.0;
   T = p/UG[0]/(8.314*meanW);

   CK_Reaction_Rate(p, T, Y, rate);
*/
      Uo[0] = UG[0];
//      Uo[1] = UG[1]/UG[0];
//      Uo[2] = UG[2]/UG[0];
//      Uo[3] = UG[5]/UG[0];  
 
//Uo[0] = 0.0;
   return xf_OK;
}

//*****************************************************************************************
static int
ParseInputFile(const char ArguFile[], char xfaFileName[], char dataFilePrefix[], int *dataFileBatchBegin, int *dataFileBatchInv, int *dataFileBatchEnd, enum xfe_Bool *WhtRectDomainOutput, real *DomainOutputParams, enum xfe_Bool *WhtCutLineOutput, real *LineOutputParams, enum xfe_Bool *WhtPointOutput, real *PointOutputParams, enum xfe_Bool *WhtIntergalQuanOutput, enum xfe_Bool *WhtToleVacuum)
{
    int i, j, flag, ierr;
    FILE *fp;
    char dummy[xf_MAXSTRLEN], *pdummy;
    char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
    
    ierr = xf_Error(xf_fopen(ArguFile, "r",&fp));
    if (ierr != xf_OK) return ierr;
    
    /***********~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~************/
    do{
        /* read line of job file */
        if (fgets(dummy, xf_MAXLINELEN, fp) == NULL)
            continue;
        pdummy = dummy;
        
        /* if blank line or comment, continue to next line */
        if (xf_TrimAndCheckBlank(&pdummy, xf_MAXLINELEN)) continue;
        
        /* read key and value from line */
        ierr = xf_Error(xf_ReadKey(pdummy, "=", key, value, xf_MAXSTRLEN));
        if (ierr != xf_OK) return ierr;
        
        /* match the input keyvalue */
        for(i=0; i<Tecout_Last; i++)
            if(strcmp(TecOut_ParamsName[i], key)==0)
                flag = i;
        
        switch (flag) {
            case TecOut_xfaFileName:
                strcpy(xfaFileName, value);
                break;
            case TecOut_dataFilePrefix:
                strcpy(dataFilePrefix, value);
                break;
            case TecOut_dataFileBatchBegin:
                ierr = sscanf(value, "%d", dataFileBatchBegin);
                if (ierr != 1) return xf_STRING_ERROR;
                break;
            case TecOut_dataFileBatchInv:
                ierr = sscanf(value, "%d", dataFileBatchInv);
                if (ierr != 1) return xf_STRING_ERROR;
                break;
            case TecOut_dataFileBatchEnd:
                ierr = sscanf(value, "%d", dataFileBatchEnd);
                if (ierr != 1) return xf_STRING_ERROR;
                break;
            case TecOut_WhtRectDomainOutput:
                if(strcmp("False", value)==0){
                    (*WhtRectDomainOutput) = xfe_False;
                }
                else {
                    (*WhtRectDomainOutput) = xfe_True;
                    //load parameters for domain output
                    for(j=0; j<5; j++)
                        fscanf(fp, "%lf", &DomainOutputParams[j]);
                }
                break;
            case TecOut_WhtCutLineOutput:
                if(strcmp("False", value)==0){
                    (*WhtCutLineOutput) = xfe_False;
                }
                else {
                    (*WhtCutLineOutput) = xfe_True;
                    //load parameters for line output
                    for(j=0; j<5; j++)
                        fscanf(fp, "%lf", &LineOutputParams[j]);
                }
                break;
            case TecOut_WhtPointOutput:
                if(strcmp("False", value)==0){
                    (*WhtPointOutput) = xfe_False;
                }
                else {
                    (*WhtPointOutput) = xfe_True;
                    //load parameters for point output
                    for(j=0; j<2; j++)
                        fscanf(fp, "%lf", &PointOutputParams[j]);
                }
                break;
            case TecOut_WhtIntergalQuanOutput:
                if(strcmp("False", value)==0){
                    (*WhtIntergalQuanOutput) = xfe_False;
                }
                else {
                    (*WhtIntergalQuanOutput) = xfe_True;
                }
                break;
            case TecOut_WhtToleVacuum:
                if(strcmp("False", value)==0){
                    (*WhtToleVacuum) = xfe_False;
                }
                else {
                    (*WhtToleVacuum) = xfe_True;
                }
                break;
            default:
                xf_printf("Input File Error!\n");
                return xf_NOT_SUPPORTED;
                break;
        }
        
    }while(feof(fp) == 0);
    

    /***********~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~************/
    ierr = xf_Error(xf_fclose(fp));
    if (ierr != xf_OK) return ierr;

    return xf_OK;
    
}
/*********************************************************************************/
static int
FoundIntersectElement(xf_All *All, const real xorg, const real yorg, const real xlen,
      const real ylen, const real h, real *Uout, real *refxy, int *Ulap, const int Nx, 
      const int Ny)
{
    int i, j, k, ierr, dim;
    int Nv, egrp, elem;
    int nx0, nxn, ny0, nyn,  n0, n1, n01[4];
    int iedge, nedge, nnode, *Node;
    real dtmp, xmin, xmax, ymin, ymax;
    real *xref, *xglob, *x1, *x0, x[3];
    enum xfe_Bool Intersect, Inside, ErrGotRefXY;
    xf_Mesh *Mesh;

    Mesh = All->Mesh;
    dim  = Mesh->Dim;
    Nv   = dim + DomainNumVars; 

    //loop over the element to find elements with intersection to rectangle
    for(egrp=0;egrp<All->Mesh->nElemGroup;egrp++){
        for(elem=0;elem<All->Mesh->ElemGroup[egrp].nElem;elem++){
            
            Node  = Mesh->ElemGroup[egrp].Node[elem];
            nnode = Mesh->ElemGroup[egrp].nNode;
            
            Intersect = xfe_False;

            ierr = xf_Error(xf_Alloc((void **) &xref, dim*nnode, sizeof(real)));
            if (ierr != xf_OK) return ierr;
            
            for(k=0; k<nnode; k++){
                xglob = Mesh->Coord[Node[k]];
                if(xglob[0]>= xorg && xglob[0]<=(xorg+xlen) && xglob[1]>=yorg && xglob[1]<=(yorg+ylen))
                { Intersect = xfe_True;
                    break; }
            }//k
            
            //order point sequance
            ierr = xf_Error(xf_Basis2nFace(Mesh->ElemGroup[egrp].QBasis, &nedge));
            if (ierr != xf_OK) return ierr;
            
            for(iedge=0; iedge<nedge; iedge++){
                ierr = xf_Error(xf_Q1NodesOnFace(Mesh->ElemGroup[egrp].QBasis, 1, iedge, &k, n01));
                if (ierr != xf_OK) return ierr;
                
                n0 = n01[0]; n1 = n01[1];
                x0 = Mesh->Coord[Node[n0]];
                x1 = Mesh->Coord[Node[n1]];
                
                xref[2*iedge]   = x0[0];
                xref[2*iedge+1] = x0[1];

            }
            
            if(Intersect == xfe_False){
                
                for(iedge=0; iedge<nedge; iedge++){
                    ierr = xf_Error(xf_Q1NodesOnFace(Mesh->ElemGroup[egrp].QBasis, 1, iedge, &k, n01));
                    if (ierr != xf_OK) return ierr;
                    //node of the edge (line segment)
                    n0 = n01[0]; n1 = n01[1];
                    //x0 = xref + 2*n0; x1 = xref + 2*n1;
                    x0 = Mesh->Coord[Node[n0]];
                    x1 = Mesh->Coord[Node[n1]];
                    
                    //Whether edge intersect rectangle
                    ierr = xf_Error(WhtEdgeIntersectRect(x0, x1, xorg, yorg, xlen, ylen, &Intersect));
                    if (ierr != xf_OK) return ierr;
                    
                    if(Intersect)
                        break;
                }
            }//if...
            
            //Finish check whether intersect. If yes, dump data to rectangle mesh; no, continue;
            if(Intersect)
            {
                //first, find the area of bounding box for this element
                xmin = xref[0];  xmax = xmin; ymin = xref[1]; ymax = ymin;
                for(k=1; k<nnode; k++){
                   xglob = xref + 2*k;
                   if(xglob[0]>xmax) xmax = xglob[0];
                   if(xglob[0]<xmin) xmin = xglob[0];
                   if(xglob[1]>ymax) ymax = xglob[1];
                   if(xglob[1]<ymin) ymin = xglob[1];
                }//k
            
                //second, pinning the bounding box inside the output mesh
                dtmp = (xmin - xorg)/h; nx0 = (long)dtmp - 1; nx0 = (nx0 > 0) ? nx0 : 0;
                dtmp = (xmax - xorg)/h; nxn = (long)dtmp + 1; nxn = (nxn > Nx) ? Nx : nxn;
                dtmp = (ymin - yorg)/h; ny0 = (long)dtmp - 1; ny0 = (ny0 > 0) ? ny0 : 0;
                dtmp = (ymax - yorg)/h; nyn = (long)dtmp + 1; nyn = (nyn > Ny) ? Ny : nyn;
            
                //third, check if the targeting point in the element and then try to get
                //corresponding coordinates in reference element
                for(i=nx0; i<nxn; i++)
                    for(j=ny0; j<nyn; j++)
                    {
                        x[0] = Uout[i*Ny*Nv + j*Nv + 0];
                        x[1] = Uout[i*Ny*Nv + j*Nv + 1];
                    
                        //whether this point is really inside the element
                        ierr = xf_Error(WhtPointInPolygon(nnode, xref, x[0], x[1], &Inside));
                        if (ierr != xf_OK) return ierr;
                    
                    
                        if(Inside)
                        {
                            k = Ulap[i*Ny*9 + j*9 + 0];
                            Ulap[i*Ny*9 + j*9 + 0]++;
                        
                            Ulap[i*Ny*9 + j*9 + 2*k + 1] = egrp;
                            Ulap[i*Ny*9 + j*9 + 2*k + 2] = elem;
                            x0 = refxy + i*Ny*8 + j*8 + 2*k + 0;
                       
                            //check the coordinates in reference element
                            ierr = xf_Error(xf_Glob2RefElem(Mesh, egrp, elem, x, 1.e+3*MEPS, xfe_False, x0, &ErrGotRefXY));
                            if (ierr != xf_OK) return ierr;
                        
                            if(!ErrGotRefXY)
                            { xf_printf("Not succeed to find reference coordinates!~\n");
                                return xf_Error(xf_NOT_SUPPORTED); }
                        }
                        else
                            continue;
                    }
            }
            else
                xf_Release((void*) xref);
        
        }//elem
    }//egrp
            
    return xf_OK;
}
/*********************************************************************************/
int
main(int argc, char * argv[])
{
    int ierr, i, j, k, l, nn, nz;
    int bat_bgn, bat_end, bat_jmp, Nout, Ncut, dim, sr;
    int Nx, Ny, Nv, egrp, elem, times;
    int nnode, *Node;
    char InputFile[xf_MAXSTRLEN], DataFile[xf_MAXSTRLEN];
    enum xfe_Bool WhtRectDomainOutput, WhtCutLineOutput, WhtPointOutput, WhtIntergalQuanOutput, WhtToleVacuum;
    real DomainOutputParams[20], LineOutputParams[20], PointOutputParams[20];
    real xlen, ylen, xorg, yorg, h;
    real lx1[2], lx2[2], lr, px[2];
    FILE *fdomout;
    xf_Data *D, *N, *Aux1, *Aux2;
    xf_Vector *U, *V_Aux1, *V_Aux2; 
    
    char dummy[xf_MAXSTRLEN], *pdummy, Name[50][xf_MAXSTRLEN];
    real dtmp, fac, dtmp1, dtmp2;
    real xmin, xmax, ymin, ymax, cutpln[3];
    int nx0, nxn, ny0, nyn, n0, n1, n01[4];
    int *Ulap, iedge, nedge, nlist;
    real *Uout, *UG, *UGm, *UGlin, *Uo, x[3], *EU;
    real *xref, *xglob, *x1, *x0, *xlin, *refxy;
    FILE *fp, *flin;
    enum xfe_Bool Intersect, Inside, CUT, Cutplane;
    enum xfe_Bool ErrGotRefXY, IntegralQuan;
    xf_All *All;
    Yu_Model Model;
    xf_Mesh *Mesh;
    xf_BasisData *PhiData = NULL;
    
    xf_printf("=== lydg_TecOut, postprocessing data for tecplot visualization===\n");
    xf_printf("=== Authored by Yu Lv, at Stanford===\n");
    xf_printf("=== email: ylv@stanford.edu===\n");
        
    //check number of arguments
    if (argc != 2) {
        xf_printf("Please specify correct arguments for this postprocessing!\n");
        xf_printf("<input.pst> with all regularization processing\n");
        xf_printf("Ouput files have the default name according to the types\n");
        xf_printf("\n");
        return xf_Error(xf_INPUT_ERROR);
    }
    
    //Load input file
    ierr = xf_Error(ParseInputFile(argv[1], InputFile, DataFile, &bat_bgn, &bat_jmp, &bat_end, 
                    &WhtRectDomainOutput, DomainOutputParams, &WhtCutLineOutput, LineOutputParams,
                    &WhtPointOutput, PointOutputParams, &WhtIntergalQuanOutput, &WhtToleVacuum));
    if (ierr != xf_OK) return ierr;
    
    
    //routine initialization
    //1.
    ierr = xf_Error(PullinModel(&Model));
    if (ierr != xf_OK) return ierr;
    
    //2.
    ierr = xf_Error(xf_CreateAll(&All, xfe_False));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReadAllBinary(InputFile, All));
    if (ierr!=xf_OK) return ierr;
    
    //3.
    Mesh = All->Mesh;
    dim  = All->Mesh->Dim;
    sr   = Model.nVars;
    
    
    //domain based output
    if(WhtRectDomainOutput){
        xorg = DomainOutputParams[0];
        yorg = DomainOutputParams[1];
        xlen = DomainOutputParams[2];
        ylen = DomainOutputParams[3];
        h    = DomainOutputParams[4];
        
        //header
        fdomout = fopen("DomainOutput.dat","w");
        ierr = xf_Error(DomainOutputHeader(fdomout));
        if (ierr != xf_OK) return ierr;
        
        //make mesh according to specification
        dtmp = xlen/h; Nx = (long) dtmp + 1;
        dtmp = ylen/h; Ny = (long) dtmp + 1;
        Nv   = dim + DomainNumVars;
        ierr = xf_Error(xf_Alloc( (void **) &Uout, Nx*Ny*Nv, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_Alloc( (void **) &refxy, Nx*Ny*8, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_Alloc( (void **) &Ulap, Nx*Ny*9, sizeof(int)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_Alloc( (void **) &UG, Model.nVars, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_Alloc( (void **) &UGm, DomainNumVars, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_Alloc( (void **) &Uo, DomainNumVars, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        
        //get x y coordinates
        for(i=0; i<Nx; i++)
            for(j=0; j<Ny; j++)
            {
                dtmp = (real)i;
                Uout[i*Ny*Nv + j*Nv + 0] = xorg + dtmp*h;
                dtmp = (real)j;
                Uout[i*Ny*Nv + j*Nv + 1] = yorg + dtmp*h;
                
                //record egrp, elem and overlap (occur at interface) correspond to the data point
                Ulap[i*Ny*9 + j*9 + 0] = 0;
            }
        
        //found elements that are inside or intersects the output domain
        ierr = xf_Error(FoundIntersectElement(All, xorg, yorg, xlen, ylen, h, Uout, refxy, Ulap, Nx, Ny));
        if (ierr != xf_OK) return ierr;
       
        //check search algorithm: all the point in output region have to be found.
        for(i=0; i<Nx; i++)
            for(j=0; j<Ny; j++)
                if(Ulap[i*Ny*9 + j*9 + 0] == 0 && !WhtToleVacuum){
                    
                    xf_printf("Output node is not found at (%lf,%lf); defect in algorithm~\n", Uout[i*Ny*Nv + j*Nv + 0], Uout[i*Ny*Nv + j*Nv + 1]);
                    return xf_Error(xf_NOT_SUPPORTED);
                }
            
        
        //load data file one-by-one
        for(nz = bat_bgn; nz<=bat_end; nz+=bat_jmp){
            //get data file name
            memset(dummy,0,sizeof(dummy));
            sprintf(dummy, "%s%d.data", DataFile, nz);
            
            printf("%s\n", dummy);
            
            ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, dummy, All->DataSet));
            if (ierr!=xf_OK) return ierr;
            
            //Get element based data
            ierr = xf_FindPrimalState(All->DataSet, 0, &D, NULL);
            if (ierr!= xf_OK) return ierr;
            
            U = (xf_Vector *) D->Data;
            
            //read auxiliary vector; if not found, still continue;
            ierr = xf_FindDataByTitle(All->DataSet, "MaxPressure", xfe_Vector, &Aux1);
            if(ierr == xf_NOT_FOUND)
            {
                xf_printf("requested auxiliary data is not found; still continue\n");
                Aux1 = NULL;
            }
            else
                V_Aux1 = (xf_Vector *) Aux1->Data;
            
            ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &Aux2);
            if(ierr == xf_NOT_FOUND)
            {
                xf_printf("requested auxiliary data is not found; still continue\n");
                Aux2= NULL;
            }
            else
                V_Aux2 = (xf_Vector *) Aux2->Data;
            
            for(i=0; i<Nx; i++)
                for(j=0; j<Ny; j++)
                {
                    times = Ulap[i*Ny*9 + j*9 + 0];
                    //logistic check
                    if(times > 4)
                    {
                        xf_printf("number of overlapping %d\n", times);
                        return xf_Error(xf_NOT_SUPPORTED);
                    }
                    if(times == 0 && WhtToleVacuum)
                    {
                       for(k=dim; k<Nv; k++)
                          Uout[i*Ny*Nv + j*Nv + k] = 1.0e+10;
                       
                       continue;
                    }
                    
                    
                    for(l=0; l<sr; l++)
                        UGm[l] = 0.0;

                    
                    for(k=0; k<times; k++)
                    {
                        egrp  = Ulap[i*Ny*9 + j*9 + 2*k + 1];
                        elem  = Ulap[i*Ny*9 + j*9 + 2*k + 2];
                        
                        // the auxiliary vector data at first
                        dtmp1 = 0.;
                        if(Aux1 != NULL)
                            dtmp1 += V_Aux1->GenArray[egrp].rValue[elem][0];
                        else
                            dtmp1 = 0.;
                        
                        dtmp2 = 0.;
                        if(Aux2 != NULL)
                            dtmp2 += V_Aux2->GenArray[egrp].rValue[elem][0];
                        else
                            dtmp2 = 0.;

                        
                        x0 = refxy + i*Ny*8 + j*8 + 2*k + 0;
                        
                        //get value of basis function for current point
                        ierr = xf_Error(xf_EvalBasis(Mesh->ElemGroup[egrp].QBasis, Model.order, xfe_True,
                                                     1, x0, xfb_Phi, &PhiData));
                        if (ierr != xf_OK) return ierr;
                        
                        nn = PhiData->nn;
                        
                        EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]
                        
                        xf_MxM_Set(PhiData->Phi, EU, 1, nn, sr, UG);
                        
                        ierr = xf_Error(DumpOutputVariable(&Model, UG, Uo, dtmp1, dtmp2));
                        if (ierr != xf_OK) return ierr;
                        
                        for(l=0; l<DomainNumVars; l++)
                            UGm[l] += Uo[l];
                        
                    }//k
                    
                    for(l=0; l<DomainNumVars; l++)
                        Uo[l] = UGm[l] / (real) times;
                    
                    //put in the storage for output
                    for(k=dim; k<Nv; k++)
                        Uout[i*Ny*Nv + j*Nv + k] = Uo[k-dim];
                    
                }//(i, j)
            
            sprintf(dummy, "\nZONE i=%d, j=%d, f=point\n", Ny, Nx);
            ierr = xf_fwrite(dummy,sizeof(char),strlen(dummy),fdomout);
            if (ierr!=xf_OK) return ierr;
            
            //Begin to dump data
            for(i=0; i<Nx; i++)
                for(j=0; j<Ny; j++){
                    for(k=0; k<Nv; k++)
                    {
                        sprintf(dummy, "%lf ", Uout[i*Ny*Nv + j*Nv + k]);
                        ierr = xf_fwrite(dummy,sizeof(char),strlen(dummy),fdomout);
                        if (ierr != xf_OK) return ierr;
                    }//(i, j, k)
                    
                    sprintf(dummy,"\n");
                    ierr = xf_fwrite(dummy,sizeof(char),strlen(dummy),fdomout);
                    if (ierr != xf_OK) return ierr;
                }
            
            //Destroy all the data in All->DataSet but not free it
            D = All->DataSet->Head;
            while(D != NULL){
                N = D->Next;
                ierr = xf_Error(xf_DestroyDataInSet(All->DataSet, D));
                if (ierr != xf_OK) return ierr;
                D = N;
            }
            
        }//nz (number of zones)
        
        //release
        fclose(fdomout);
        xf_Release((void*) Uout);
        xf_Release((void*) Ulap);
        xf_Release((void*) UG);
        xf_Release((void*) UGm);
        xf_Release((void*) Uo);
        xf_Release((void*) refxy);
    }
    
    if(WhtCutLineOutput){
        lx1[0] = LineOutputParams[0];
        lx1[1] = LineOutputParams[1];
        lx2[0] = LineOutputParams[2];
        lx2[1] = LineOutputParams[3];
        lr     = LineOutputParams[4];  //level of refinement: i.e. if 10 the intersective line-segment will have 10 sub-segment.
        
        //header
        ierr = xf_Error(LineOutputHeader());
        if (ierr != xf_OK) return ierr;
    }
    
    if(WhtPointOutput){
        //now only put one point
        px[0] = PointOutputParams[0];
        px[1] = PointOutputParams[1];
        
        //header
        ierr = xf_Error(PointOutputHeader());
        if (ierr != xf_OK) return ierr;
    }
    
    if(WhtIntergalQuanOutput){
        //header
        ierr = xf_Error(IntegralOuputHeader());
        if (ierr != xf_OK) return ierr;
    }
    
    

    ierr = xf_Error(DestroyModel(&Model));
    if (ierr!=xf_OK) return ierr;
 
    ierr = xf_Error(xf_DestroyAll(All));
    if (ierr!=xf_OK) return ierr;
    
    xf_printf("xf_TecOut finished - Outputted\n");

    return xf_OK;
}
static int
WhtPointInPolygon(const int nvert, const real *xnode, const real testx,  const real testy, enum xfe_Bool *Inside)
{
    int i, j;
    real vertx[10], verty[10];
    (*Inside) = xfe_False;
    
    //first check node
    for (i = 0; i < nvert; i++){
       vertx[i] = xnode[2*i + 0];
       verty[i] = xnode[2*i + 1];
       if(fabs(vertx[i]-testx)<MEPS && fabs(verty[i]-testy)<MEPS)
        {(*Inside) = xfe_True; return xf_OK; }
    }
    
    //second check edge
    for (i = 0, j = nvert-1; i < nvert; j = i++)
        if(fabs((testx - vertx[i]) * (verty[j]-verty[i]) - (vertx[j]-vertx[i]) * (testy-verty[i])) < MEPS)
        {
            //check if in the range
            if(((verty[i]>testy) != (verty[j]>testy)))
                (*Inside) = xfe_True;
            if(((vertx[i]>testx) != (vertx[j]>testx)))
                (*Inside) = xfe_True;
            
            if((*Inside))
                return xf_OK;
        }
    
    for (i = 0, j = nvert-1; i < nvert; j = i++) {
        if ( ((verty[i]>testy) != (verty[j]>testy)) &&
            (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
            (*Inside) = !(*Inside);
    }
    
    return xf_OK;
}
static int
WhtEdgeIntersectRect(const real *x0, const real *x1, const real xorg, const real yorg, 
                     const real xlen, const real ylen, enum xfe_Bool *Intersect)
{
   real xc, yc, tmp, xi, yi;
   xc = xorg + xlen;
   yc = yorg + ylen;

   //remember x0, x1 should not locate inside the rectangle
   //first if both node on the bottom, top, left and right
   if(x0[1]<yorg && x1[1]<yorg){         //bottom
      (*Intersect) = xfe_False;
      return xf_OK;
   }
   else if(x0[1]>yc && x1[1]>yc){        //top
      (*Intersect) = xfe_False;
      return xf_OK;
   }
   else if(x0[0]<xorg && x1[0]<xorg){    //left
      (*Intersect) = xfe_False;
      return xf_OK;
   }
   else if(x0[0]>xc && x1[0]>xc){        //right
      (*Intersect) = xfe_False;
      return xf_OK;
   }

   //second, check edge intersection
   if(fabs(x0[0] - x1[0]) > MEPS)
   {
      tmp = (x1[1] - x0[1])/(x1[0] - x0[0]);
      yi  = tmp*(xorg - x0[0]) + x0[1];
      if(yi >= yorg && yi <= yc){
         (*Intersect) = xfe_True;
         return xf_OK;
      }

      yi = tmp*(xc - x0[0]) + x0[1];
      if(yi >= yorg && yi <= yc){
         (*Intersect) = xfe_True;
         return xf_OK;
      }
   }

   if(fabs(x0[1] - x1[1]) > MEPS)
   {
      tmp = (x1[0] - x0[0])/(x1[1] - x0[1]);
      xi = tmp*(yorg - x0[1]) + x0[0];
      if(xi >= xorg && xi <= xc){
         (*Intersect) = xfe_True;
         return xf_OK;
      }

      xi = tmp*(yc - x0[1]) + x0[0];
      if(xi >= xorg && xi <= xc){
         (*Intersect) = xfe_True;
         return xf_OK;
      }
   }

   (*Intersect) = xfe_False;
   return xf_OK;

   //if correctly implemented; code should not reach here
   xf_printf("Error in checking element intersection with output region!~\n");
   return xf_Error(xf_OUT_OF_BOUNDS); 
  
}
