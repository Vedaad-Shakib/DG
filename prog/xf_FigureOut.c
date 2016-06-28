/*
  FILE:  xf_Plot.c

  This is the xflow plotter.  It plots a solution from a binary .xfa
  with optional supplementary .data files.

*/

#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_All.h"
#include "xf_Mesh.h"
#include "xf_Param.h"
#include "xf_Basis.h"
#include "xf_Data.h"
#include "xf_EqnSet.h"
#include "xf_Math.h"
#include "xf_MeshTools.h"
#include "xf_Quad.h"
#include "xf_SolverStruct.h"
#include "xf_SolverTools.h"
#include "xf_Residual.h"
#include "xf_EqnSetHook.h"
#include "xf_ParamDefault.h"
#include "xf_MeshMotion.h"
#include "xf_MeshMotionGCL.h"
#include "xf_MeshMotionIO.h"
#include "xfYu_Statistics.h"
#include "xf_Arg.h"
#include "xfYu_Model.h"

#include <stdlib.h>
#ifdef __APPLE__

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#else

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#endif



/*------------------*/
/* Enumerated Types */
/*------------------*/


/* Types of display groups (structure defined further down) */
enum xfe_DisplayGroupType{
  xfe_DGCutPlane,
  xfe_DGBoundary,
  xfe_DGIsoSurf,
  xfe_DGLast
};

/* corresponding names */
static char *xfe_DisplayGroupName[xfe_DGLast] = {
  "CutPlane",
  "Boundary",
  "IsoSurf",
};

/* Derived-types for vectors */
enum xfe_VectorDeriveType{
  xfe_VectorDeriveNone,
  xfe_VectorDeriveScalar,
  xfe_VectorDeriveVariableSet,
  xfe_VectorDeriveLast
};

/* corresponding names */
static char *xfe_VectorDeriveName[xfe_VectorDeriveLast] = {
  "None",
  "Scalar",
  "VariableSet",
};


/* Types of rendering color maps*/
enum xfe_RenderMapType{
  xfe_RenderMapStandard,
  xfe_RenderMapInverse,
  xfe_RenderMapBW,
  xfe_RenderMapInverseBW,
  xfe_RenderMapLast
};

/* corresponding names */
static char *xfe_RenderMapName[xfe_RenderMapLast] = {
  "Standard",
  "Inverse",
  "Black-White",
  "Inverse Black-White",
};


/*------------*/
/* Structures */
/*------------*/


/* Each DisplayTri structure is a group of triangles used for
   displaying the solution.  These groups correspond to boundary
   groups or cut-planes.  The triangles are related back to the
   element from which they came.  Multiple triangles can point to the
   same element. */
typedef struct
{
  int Version; // Current version of the display structure

  // Nodes
  int nNode, nNode0, DNode;
  real *xglob;  // spatial coordinates of display nodes [nNode*3]
  real *xref;   // elem ref coordinates of display nodes [nNode*3]
  real *xplane; // cut-plane coordinates of display nodes [nNode*3]

  real *Xglob;  // ref-space coords in case of mesh motion [nNode*3]

  // Triangles
  int nTri, nTri0, DTri;
  int *TriList; // list of triangle nodes [nTri*3]

  // Normals
  int nNormal;
  real *Normals;

  // Active flag
  int *Active; // vector indicatiting whether tri is active [nTri]

  // element info
  int *Tri2Elem; // stores parent egrp, elem [nTri*2]
}
xf_DisplayTri;

/* Display loops are used for element boundaries.  For example, when a
   cut plane intersects an 3D element, the surface boundary of the
   element yields one or more closed loops in the cut plane.  All of
   these loops for one display group are stored according to this
   structure.*/
typedef struct
{
  int nLoop;
  int *nNode; // number of nodes per loop [nLoop]
  int **NodeList; // node lists for all loops [nLoop][nNode[iLoop]]
  int *Active; // vector indicatiting whether loop is active [nLoop]
               // note, Active < 0 implies loop is sequence of edges
               // |Active| == 2 means loop is active
}
xf_DisplayLoop;

/* Display strips are used for domain boundaries.*/
typedef struct
{
  int nStrip;
  int *nNode; // number of nodes per strip [nStrip]
  int **NodeList; // node lists for all strip [nStrip][nNode[iStrip]]
}
xf_DisplayStrip;


/* A display group could arise from a cut-plane or a boundary-face
   group.  It contains pointers to a set of display triangles and a
   set of display loops.  Additionally, this structure stores info on
   the mesh elements that participate in this display group,
   specifically including the set of nodes belonging to each element.
   The node locations themselves are stored in the DTri field. */
typedef struct
{
  enum xfe_Bool Active; // true means this display group should be displayed
  enum xfe_DisplayGroupType Type; // CutPlane, boundary, etc.
  int ibfgrp; // boundary face group # for boundary type
  real isoval; // iso-surface value
  
  xf_DisplayTri   *DTri;   // display triangles
  xf_DisplayLoop  *DLoop;  // display loops of edges
  xf_DisplayStrip *DStrip; // display strips of boundary edges

  int nelem; // number of elements participating in this display group
  int *Egrp; // vector of participating element groups [nelem]
  int *Elem; // vector of participating elements [nelem]
  int *nDNode; // vector of number of display nodes [nelem]
  int **DNode; // vector of display nodes [nelem, nDNode[ielem]]
}
xf_DisplayGroup;

/* The display node quantities for each scalar are stored according to
   this structure.*/
typedef struct
{
  int Version;  // display structure version for which the data is valid
  int nNode;
  real *Value;  // scalar values [nNode]
  real *ValueX; // ref scalar values in case of mesh motion [nNode]
  real *Color;  // RGB color values [3*nNode]
}
xf_ScalarData;

/* This structure is a wrapper for the ScalarData: each DisplayScalar
   consists of the nodal data along with an identifying scalar name
   and min/max values for plotting. */
typedef struct
{
  char Name[xf_MAXSTRLEN];
  
  xf_ScalarData *Data; // scalar plotting data for each display group
  real Vmin, Vmax; // min and max values, for plotting
}
xf_DisplayScalar;


/* A DisplayVector is a collection of display scalars.  Typically,
   each such vector originates from a vector in the data. */
typedef struct
{
  char DataName[xf_MAXSTRLEN];
  int SetIndex; // index in vectorset, or -1 if not part of a vectorset
  xf_Vector *Vector; // Associated All Vector structure, if any
  enum xfe_VectorDeriveType DeriveType; // how is this DV derived from Vector?
  int nScalar;
  int CurrentDScalar; // current scalar being displayed
  xf_DisplayScalar *Scalar;
}
xf_DisplayVector;


/* This structure stores information on each display window. */
typedef struct
{
  // Assigned by glut
  int Handle;

  // In pixel coords
  int PixTopLeft[2];
  int PixWidth[2];

  // In geometry coords, if applicable
  real Center[2];
  real Width[2];
}
xf_WindowType;


/* The state of the mouse (buttons pressed, coordinates, etc.) */
typedef struct
{
  enum xfe_Bool LeftButtonDown;
  enum xfe_Bool MiddleButtonDown;
  enum xfe_Bool RightButtonDown;
  int Modifier; // shift, ctl, or alt key pressed on last event
  int x, y; // from last event
}
xf_MouseStateType;

/* For setting the 3D view position */
typedef struct
{
  int ViewAxis; // [0,1,2] = [x,y,z]; +1 = elev axis, +2 = azimuth axis
  real Distance; // from origin
  real Azimuth; // angle (deg)
  real Elevation; // angle (deg)
  real Twist; // angle (deg)
  real Origin[3]; // origin x,y,z about which to rotate
}
xf_CameraType;

/* For storing contour segments */
typedef struct
{
  int nSeg;   // number of segments
  int nSeg0;  // number of segments allocated
  int dSeg;   // delta for reallocation
  real *xseg; // [x0,y0,z0, x1,y1,z1], 6*nSeg
  real *gxseg; // global space version of xseg
}
xf_ContoursType;

/*-----------*/
/*  Defines  */
/*-----------*/

#define xfp_CBAR_SIZE 32
#define xfp_MAX_NEG_POINT 100

/*------------------*/
/* Global variables */
/*------------------*/

xf_All *All;  
// pointer to All structure
real ModelBBox[6]; 
// Bounding box for the model [xmin, xmax, ... ]
real ModelSize; 
// Characteristic (e.g. max) dimension
real ModelOffset; 
// Maximum of (distance to origin in each dim)
real CutPlane[4] = {0.0, -1.0, 0.0, 0.0};
// Cut plane position: [a,b,c,d] in a*x + b*y + c*z + d = 0
xf_DisplayVector *DVector = NULL; 
// Array of registered display vectors
int nDVector = 0; 
// Number of display vectors
int CurrentDVector = 0; 
// Current display vector
enum xfe_Bool MeshOn = xfe_True;
// True if mesh is to be shown
float MeshLineWidth = 2.0;
// Mesh line width
enum xfe_Bool SubMeshOn = xfe_False;
// True if sub-mesh is to be shown
enum xfe_Bool RenderFlag = xfe_False; 
// True if rendering a solution
enum xfe_Bool RenderTemporarilyOff = xfe_False; 
// True if rendering temporarily turned off (e.g. for mesh motion)
enum xfe_Bool LightingFlag = xfe_False; 
// True if using lighting (3d)
enum xfe_Bool MeshFillFlag = xfe_False; 
// True if rendering interior of mesh with a solid color (e.g. white)
enum xfe_Bool DrawingHills = xfe_False; 
// True if plotting hills
enum xfe_Bool UseLogScale = xfe_False; 
// True if plotting log(scalar) instead of scalar



xf_MouseStateType MouseState; 
// State of the mouse 
xf_WindowType WinParent; 
// Parent window
xf_WindowType Win2D;  
// 2D window (sub-window in WinParent)
xf_WindowType Win3D;  
// 3D window (sub-window in WinParent)
xf_WindowType WinAux; 
// Auxiliary window (sub-window in WinParent)
xf_WindowType WinBorder; 
// Border window (sub-window in WinParent)
xf_WindowType *WinBig, *WinSmall, *WinTemp; 
// Pointers to the big/small display windows
xf_CameraType Camera; 
// Camera information for the 3d window
enum xfe_Bool AxesFlag = xfe_False; 
// Flag indicating axes on or off
int UpAxis2D = 1; 
// [0,1,2] = [x,y,z]: up direction axis for viewing the 2D cut plane

xf_DisplayGroup *DGroup; 
// vector of display groups [nDGroup]
int nDGroup; 
// Number of display groups (total)
int nActiveDGroup;
// Number of active display groups
int *ActiveDGroup; 
// Indices of active groups [nActiveDGroup]
int iDGroupCutPlane = -1; 
// Index of cut plane display group
int iDGroupIsoSurf = -1; 
// Index of iso-surface display group
int *iDGroupBFG = NULL; 
// Indices of BFG display groups
int **Elem2Refine = NULL; 
// Refinement level for each element Elem2Refine[egrp][elem]

enum xfe_Bool OutputsOn = xfe_False;
// if True, output points will be drawn
int *OutputIsActive = NULL;
// ActiveOutputs[iOutput] == 1 means output is on
int nOutputPoints = 0;
// Number of output points being drawn
real *xOutputPoints = NULL;
// Coordinates of output points
real OutputPointSize = 5.0;
// output point size in pixels

enum xfe_Bool ElemMarkerOn = xfe_False;
// if True, element id marker is on
real xElemMarker[3];
// Coordinates of element marker (for identifying elements)
real ElemMarkerSize = 0.1;
// size of element marker

int ColorMapSize = 256;
// For color indexing (appears not supported in freeglut)
real CBarColor[3*(xfp_CBAR_SIZE+1)];
// Color-bar colors
float xf_FG[3] = {0.0, 0.0, 0.0}; // black
// primary foreground color
float xf_BG[3] = {1.0, 1.0, 1.0}; // white
// primary background color
const char StateFile[] = "plot.state";
// state file name

enum xfe_Bool HaveEqnSet = xfe_False;
// True if we have an equation set

enum xfe_Bool ContoursOn = xfe_False;
// True if contours are on
int nContours = 10;
// Number of contous
xf_ContoursType Contours; 
// structure for storing contour segments

enum xfe_Bool LineProbeOn = xfe_False;
// True if line probe is on
enum xfe_Bool AcquiringLineProbe = xfe_False;
// True if acquiring a line probe (two mouse clicks)
int nLineProbeAcquired = 0;
// Number of points acquired
real LineProbeXY[4] = {0., 0., 0., 1.};

enum xfe_Bool NegVolOn = xfe_False;
// if True, negative volume points will be drawn
int nNegVolPoints = 0;
// number of negative volume points
real xNegVolPoints[3*xfp_MAX_NEG_POINT];
// locations of tested negative volumes

enum xfe_RenderMapType RenderMap = xfe_RenderMapStandard;
// rendering map type (for color)

enum xfe_Bool BoundaryOutline = xfe_False;
// If true, boundary outline is rendered (2d)

GLfloat light_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
// Ambient light value

enum xfe_Bool CountBoxOn = xfe_False;
// if True, elements are counted within a certain box
real CountBox[4] = {0., 0., 0., 0.};

real Time;
// simulation time
enum xfe_Bool MeshMotionInPhysical = xfe_True;
// if True, will plot mesh in physical space for mesh motion cases

/* For interpolation with mesh motion */
int MDInq = 0;
real *MDIxglob = NULL;
xf_MotionData *MDI = NULL;


/* Menu variables */
enum xfe_PlotMenuOptions{
  xfe_PlotMenu_None,        // options need to start at 1
  xfe_PlotMenu_Sub_LoadVector,
  xfe_PlotMenu_Sub_ChangeScalar,
  xfe_PlotMenu_Sub_Mesh,
  xfe_PlotMenu_Sub_Boundaries,
  xfe_PlotMenu_Sub_ColorScheme,
  xfe_PlotMenu_Sub_CutPlane,
  xfe_PlotMenu_Sub_LineProbe,
  xfe_PlotMenu_Sub_Contours,
  xfe_PlotMenu_Sub_IsoSurface,
  xfe_PlotMenu_Sub_Lighting,
  xfe_PlotMenu_Sub_ViewFrom,
  xfe_PlotMenu_Sub_Hardcopy,
  xfe_PlotMenu_Sub_ColorLimits,
  xfe_PlotMenu_DeriveVector,
  xfe_PlotMenu_ToggleRender,
  xfe_PlotMenu_SwapWindows,
  xfe_PlotMenu_Help,
  xfe_PlotMenu_Quit,
  xfe_PlotMenu_Max          // not actually a menu option  
};

int TopMenu = 1; // top menu index
int SubMenuIndex[xfe_PlotMenu_Max] = {0};
enum xfe_Bool nModifyMenu = 0; // number of menus that need mods
int ModifyMenuList[xfe_PlotMenu_Max] = {0};


/* Funtion prototypes for menu processing */
static void xf_ModifySubMenu_LoadVector();
static void xf_ModifySubMenu_ChangeScalar();
static void xf_ModifySubMenu_Boundaries();
static void xf_ModifySubMenu_Mesh();
static void xf_ModifySubMenu_ColorScheme();
static void xf_ModifySubMenu_CutPlane();
static void xf_ModifySubMenu_LineProbe();
static void xf_ModifySubMenu_Contours();
static void xf_ModifySubMenu_IsoSurface();
static void xf_ModifySubMenu_Lighting();
static void xf_ProcessMenu_LoadVector(int k);
static void xf_ProcessMenu_ChangeScalar(int k);
static void xf_ProcessMenu_Boundaries(int k);
static void xf_ProcessMenu_Mesh(int k);
static void xf_ProcessMenu_ColorScheme(int k);
static void xf_ProcessMenu_CutPlane(int k);
static void xf_ProcessMenu_LineProbe(int k);
static void xf_ProcessMenu_Contours(int k);
static void xf_ProcessMenu_IsoSurface(int k);
static void xf_ProcessMenu_Lighting(int k);
static void xf_ProcessMenu_ViewFrom(int k);
static void xf_ProcessMenu_Hardcopy(int k);
static void xf_ProcessMenu_ColorLimits(int k);


/******************************************************************/
//   FUNCTION Definition: xf_SetColorMap
static void 
xf_SetColorMap()
{
  int i;
  real s;
  GLfloat r, g, b;

  
  for (i=0; i<ColorMapSize; i++){ 
    s = ((real) i) / ((real) ColorMapSize - 1.0);
    if (s < 0.25){
      r = 0.0;
      g = 4.0*s;
      b = 1.0;
    }
    else if (s < 0.5){
      r = 0.0;
      g = 1.0;
      b = 2.0 - 4.0*s;
    }
    else if (s < 0.75){
      r = 4.0*s - 2.0;
      g = 1.0;
      b = 0.0;
    }
    else{
      r = 1.0;
      g = 4.0 - 4.0*s;
      b = 0.0;
    }
    xf_printf("s = %.10E, r = %.10E, g = %.10E, b = %.10E\n", s, r,g,b);
    glutSetColor(i, r, g, b);
  } // i

}




/******************************************************************/
//   FUNCTION Definition: xf_Scalar2RGB
static void 
xf_Scalar2RGB(real s, real *rgb)
{
  switch (RenderMap){
  case xfe_RenderMapBW:
    rgb[0] = rgb[1] = rgb[2] = s;
    break;
  case xfe_RenderMapInverseBW:
    rgb[0] = rgb[1] = rgb[2] = 1.-s;
    break;
  case xfe_RenderMapInverse:
    s = 1.0-s;
  case xfe_RenderMapStandard:
    if (s < 0.25){
      rgb[0] = 0.0;
      rgb[1] = 4.0*s;
      rgb[2] = 1.0;
    }
    else if (s < 0.5){
      rgb[0] = 0.0;
      rgb[1] = 1.0;
      rgb[2] = 2.0 - 4.0*s;
    }
    else if (s < 0.75){
      rgb[0] = 4.0*s - 2.0;
      rgb[1] = 1.0;
      rgb[2] = 0.0;
    }
    else{
      rgb[0] = 1.0;
      rgb[1] = 4.0 - 4.0*s;
      rgb[2] = 0.0;
    }
    break;
  default:
    rgb[0] = rgb[1] = rgb[2] = 0.;
    xf_printf("Unrecognized scalar map\n");
    break;
  }
}



/******************************************************************/
//   FUNCTION Definition: xf_FillColorBar
static void 
xf_FillColorBar()
{
  int k;
  // Calculate colors for  Color Bar
  for (k=0; k<=xfp_CBAR_SIZE; k++)
    xf_Scalar2RGB(((real) k) / ((real) xfp_CBAR_SIZE), CBarColor + 3*k);
}


/******************************************************************/
//   FUNCTION Definition: xf_InitGraphics
void
xf_InitGraphics(void)
{
  glClearColor(1.0, 1.0, 1.0, 0.0);
  glClearDepth(1.0);
  glDepthFunc(GL_LESS); 
  glShadeModel(GL_SMOOTH);

  glEnable(GL_DEPTH_TEST); // enables z-buffer
  glEnable(GL_POLYGON_OFFSET_LINE);
  glEnable(GL_LINE_SMOOTH);
}

/******************************************************************/
//   FUNCTION Definition: xf_DrawText
void
xf_DrawText(float x, float y, char *s, int maxlen)
{
  int k;
  glPushMatrix();
  glColor3f(0.0, 0.0, 0.0);
  glRasterPos2f(x, y);
  for (k=0; k<min(maxlen, strlen(s)); k++)
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, s[k]);
  glPopMatrix();
}



/******************************************************************/
//   FUNCTION Definition: xf_DrawDGroup
void 
xf_DrawDGroup(int iDGroup)
{
  int k, i, iDNode;
  int CurrentDScalar;
  real *Color;
  real *x, *c, *n;
  real c0[3] = {1.0, 1.0, 1.0};
  xf_DisplayGroup *DG;
  xf_DisplayTri *DTri;
  xf_DisplayLoop *DLoop;
  xf_DisplayStrip *DStrip;
  
  DG = DGroup + iDGroup;
  DTri = DG->DTri;
  DLoop = DG->DLoop;
  DStrip = DG->DStrip;

  if (MeshOn){  // Draw element mesh
    glLineWidth(MeshLineWidth);
    glColor3f(0.1, 0.1, 0.1);
    
    for (i=0; i<DLoop->nLoop; i++){
      if ((DLoop->Active[i] == 1) || (DLoop->Active[i] == -1))
	continue;
      if (DLoop->Active[i] > 0){
	glBegin(GL_LINE_LOOP);
	for (k=0; k<DLoop->nNode[i]; k++){
	  x = DTri->xglob + 3*DLoop->NodeList[i][k];
	  glVertex3f(x[0], x[1], x[2]);
	}
	glEnd();
      }
      else{ // negative implies no loop, but edge list instead
	glBegin(GL_LINES);
	for (k=0; k<DLoop->nNode[i]; k++){
	  x = DTri->xglob + 3*DLoop->NodeList[i][k];
	  glVertex3f(x[0], x[1], x[2]);
	}
	glEnd();
      }
    }
  }

  if (SubMeshOn){  // Draw subelement mesh
    glLineWidth(MeshLineWidth);
    glColor3f(0.75, 0.75, 0.75);
    
    for (k=0; k<DTri->nTri; k++){
      glBegin(GL_LINE_LOOP);
      x = DTri->xglob + 3*DTri->TriList[3*k+0];
      glVertex3f(x[0], x[1], x[2]);
      x = DTri->xglob + 3*DTri->TriList[3*k+1];
      glVertex3f(x[0], x[1], x[2]);
      x = DTri->xglob + 3*DTri->TriList[3*k+2];
      glVertex3f(x[0], x[1], x[2]);
      glEnd();
    }
  }

  if (BoundaryOutline){  // Draw boundary strips
    glLineWidth(2.0);
    glColor3f(0.1, 0.1, 0.1);
    
    for (i=0; i<DStrip->nStrip; i++){
      glBegin(GL_LINE_STRIP);
      for (k=0; k<DStrip->nNode[i]; k++){
	x = DTri->xglob + 3*DStrip->NodeList[i][k];
	glVertex3f(x[0], x[1], x[2]);
      }
      glEnd();
    }
  }

  if (MeshFillFlag && (!RenderFlag)){
    // Fill submesh elements with white
    c = c0;
    glBegin(GL_TRIANGLES);
    for (k=0; k<DTri->nTri; k++){

      if (DTri->Active[k] == 0) continue;

      if (LightingFlag){
	n = DTri->Normals + 3*k;
	glNormal3f(n[0], n[1], n[2]);
      }

      iDNode = DTri->TriList[3*k+0];
      x = DTri->xglob + 3*iDNode;
      glColor3f(c[0], c[1], c[2]);
      glVertex3f(x[0], x[1], x[2]);

      iDNode = DTri->TriList[3*k+1];
      x = DTri->xglob + 3*iDNode;
      glColor3f(c[0], c[1], c[2]);
      glVertex3f(x[0], x[1], x[2]);

      iDNode = DTri->TriList[3*k+2];
      x = DTri->xglob + 3*iDNode;
      glColor3f(c[0], c[1], c[2]);
      glVertex3f(x[0], x[1], x[2]);
    }
  
    glEnd();

  }

  if (RenderFlag){
    // Render scalar on triangles

    CurrentDScalar = DVector[CurrentDVector].CurrentDScalar;
    Color = DVector[CurrentDVector].Scalar[CurrentDScalar].Data[iDGroup].Color;
    glBegin(GL_TRIANGLES);
    for (k=0; k<DTri->nTri; k++){

      if (DTri->Active[k] == 0) continue;

      if (LightingFlag){
	n = DTri->Normals + 3*k;
	glNormal3f(n[0], n[1], n[2]);
      }

      iDNode = DTri->TriList[3*k+0];
      x = DTri->xglob + 3*iDNode;
      c = Color + 3*iDNode;
      glColor3f(c[0], c[1], c[2]);
      //glIndexd(c[0]);
      glVertex3f(x[0], x[1], x[2]);

      iDNode = DTri->TriList[3*k+1];
      x = DTri->xglob + 3*iDNode;
      c = Color + 3*iDNode;
      glColor3f(c[0], c[1], c[2]);
      //glIndexd(c[0]);
      glVertex3f(x[0], x[1], x[2]);

      iDNode = DTri->TriList[3*k+2];
      x = DTri->xglob + 3*iDNode;
      c = Color + 3*iDNode;
      glColor3f(c[0], c[1], c[2]);
      //glIndexd(c[0]);
      glVertex3f(x[0], x[1], x[2]);
    }
  
    glEnd();

  }

  glutSwapBuffers();

}



/******************************************************************/
//   FUNCTION Definition: xf_DrawContours2D
void 
xf_DrawContours2D()
{
  int k, i;
  real *x;
  
  if (Contours.nSeg <= 0) return;
  
  glLineWidth(1.0);
  glColor3f(0.1, 0.1, 0.1);

  glBegin(GL_LINES);

  for (i=0; i<Contours.nSeg; i++){
    x = Contours.xseg+6*i;
    glVertex3f(x[0], x[1], x[2]);
    glVertex3f(x[3], x[4], x[5]);
  }
  glEnd();

}


/******************************************************************/
//   FUNCTION Definition: xf_DrawContours3D
void 
xf_DrawContours3D()
{
  int k, i;
  real *x;
  
  if (Contours.nSeg <= 0) return;
  
  glLineWidth(1.0);
  glColor3f(0.1, 0.1, 0.1);

  glBegin(GL_LINES);

  for (i=0; i<Contours.nSeg; i++){
    x = Contours.gxseg+6*i;
    glVertex3f(x[0], x[1], x[2]);
    glVertex3f(x[3], x[4], x[5]);
  }
  glEnd();

}



/******************************************************************/
//   FUNCTION Definition: xf_Display2D
void 
xf_Display2D(void)
{
  int k, i, iDNode;
  int CurrentDScalar;
  real *N, *x, up[3], upn, NN;
  float ac[3], axislen;
  xf_DisplayTri *DTri;

  // Need to manipulate the ModelView matrix to move our model around.
  glMatrixMode(GL_MODELVIEW);

  // Reset to 0,0,0; no rotation, no scaling.
  glLoadIdentity(); 

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      
  // turn off color mappings
  glPixelTransferf(GL_MAP_COLOR, GL_FALSE);
  //glDisable(GL_TEXTURE_2D);

  // Draw in up coordinate axis
  if (AxesFlag){
    N = CutPlane;
    for (k=0; k<3; k++) up[k] = 0.;
    up[UpAxis2D] = 1.0;
    for (k=0, upn=0.; k<3; k++) upn += up[k]*N[k];
    axislen = 2.0*ModelSize*sqrt(1.0-upn*upn);

    for (k=0; k<3; k++) ac[k] = 0.;
    ac[UpAxis2D] = 1.0;

    glLineWidth(2.0);
    glBegin(GL_LINES);    
    glColor3f(ac[0], ac[1], ac[2]);
    glVertex3f(0.0, 0.0, 0.1);
    glVertex3f(0.0, axislen, 0.1);
    glEnd();
  }


  // draw contours if on
  if (ContoursOn) xf_DrawContours2D();

  // draw cut plane if on
  if ((iDGroupCutPlane >= 0) && (DGroup[iDGroupCutPlane].Active)){
    if (All->Mesh->Dim == 3){
      DTri = DGroup[iDGroupCutPlane].DTri;
      swap(DTri->xglob, DTri->xplane, x);
      xf_DrawDGroup(iDGroupCutPlane);
      swap(DTri->xglob, DTri->xplane, x);
    }
    else
      xf_DrawDGroup(iDGroupCutPlane);
  }

  // Draw Line Probe if on
  if (LineProbeOn){

    glLineWidth(3.0);
    glBegin(GL_LINES);    
    glColor3f(0.2, 0.2, 0.2);
    glVertex3f(LineProbeXY[0], LineProbeXY[1], 0.1);
    glVertex3f(LineProbeXY[2], LineProbeXY[3], 0.1);
    glEnd();
  }

  // Draw Count Box if on
  if (CountBoxOn){

    glLineWidth(4.0);
    glBegin(GL_LINE_LOOP);    
    glColor3f(0.2, 0.2, 0.2);
    glVertex3f(CountBox[0], CountBox[2], 0.1);
    glVertex3f(CountBox[1], CountBox[2], 0.1);
    glVertex3f(CountBox[1], CountBox[3], 0.1);
    glVertex3f(CountBox[0], CountBox[3], 0.1);
    glEnd();
  }




  glutSwapBuffers();

}


/******************************************************************/
//   FUNCTION Definition: xf_Display3D
void 
xf_Display3D(void)
{
  int i, j;
  int itet[6] = {0,1,2,3,0,1};
  real dist, a;
  float axislen;
  //GLfloat light_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
  //GLfloat light_ambient[] = { 0.8, 0.8, -0.7, 1.0 };
  real xtet[4][3];

  // turn off color mappings
  glPixelTransferf(GL_MAP_COLOR, GL_FALSE);
  //glDisable(GL_TEXTURE_2D);
  //glEnable(GL_TEXTURE_2D);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_POLYGON_OFFSET_LINE);
  glEnable(GL_LINE_SMOOTH);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  // Need to manipulate the ModelView matrix to move our model around.
  glMatrixMode(GL_MODELVIEW);

  // Reset to 0,0,0; no rotation, no scaling. 
  glLoadIdentity();
    
  dist = Camera.Distance;

  glTranslated(0.0, 0.0, -Camera.Distance);
  glRotated(-Camera.Elevation, 1.0, 0.0, 0.0);
  glRotated(-Camera.Azimuth, 0.0, 1.0, 0.0);


  // camera origin translation
  glTranslated(-Camera.Origin[0], -Camera.Origin[1], -Camera.Origin[2]);


  if (Camera.ViewAxis == 0){
    // looking from positive x
    glRotated(-90, 0.0, 1.0, 0.0);
    glRotated(-90, 1.0, 0.0, 0.0);
  }
  else if (Camera.ViewAxis == 1){
    // looking from positive y
    glRotated(90, 1.0, 0.0, 0.0);
    glRotated(90, 0.0, 1.0, 0.0);
  }
  
  // Lighting
  if (LightingFlag){
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    //glLightfv(GL_LIGHT0, GL_POSITION, light_ambient);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
  }
  else{
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
  }


  if (iDGroupBFG != NULL){
    for (i=0; i<All->Mesh->nBFaceGroup; i++)
      if ((iDGroupBFG[i] >= 0) && (DGroup[iDGroupBFG[i]].Active))
	xf_DrawDGroup(iDGroupBFG[i]);
  }

  // Draw in 3 coordinate axes
  if (AxesFlag){
    axislen = 2.0*ModelSize;
    glLineWidth(2.0);
    glBegin(GL_LINES);
    
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(axislen, 0.0, 0.0);
    
    glColor3f(0.0, 1.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, axislen, 0.0);
    
    glColor3f(0.0, 0.0, 1.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, axislen);
    
    glEnd();
  }

  // Draw in points
  if (nOutputPoints > 0){
    a = ModelSize*.01;
    glColor3f(1.0, 0.0, 0.0);
    
    for (i=0; i<nOutputPoints; i++){
      glBegin(GL_TRIANGLE_STRIP);
      for (j=0; j<3; j++) xtet[0][j] = xOutputPoints[3*i+j];
      for (j=0; j<3; j++) xtet[1][j] = xtet[0][j];
      xtet[1][2] -= 1.5*a; xtet[1][0] -= .866*a; xtet[1][1] -= 0.5*a;
      for (j=0; j<3; j++) xtet[2][j] = xtet[0][j];
      xtet[2][2] -= 1.5*a; xtet[2][0] += .866*a; xtet[2][1] -= 0.5*a;
      for (j=0; j<3; j++) xtet[3][j] = xtet[0][j];
      xtet[3][2] -= 1.5*a; xtet[3][1] += a;
      for (j=0; j<6; j++)
	glVertex3f(xtet[itet[j]][0], xtet[itet[j]][1], xtet[itet[j]][2]);
      glEnd();
    }
  }
  
  // draw element marker
  if (ElemMarkerOn){
    a = ElemMarkerSize;
    glColor3f(0.0, 0.0, 1.0);
    
    glBegin(GL_TRIANGLE_STRIP);
    for (j=0; j<3; j++) xtet[0][j] = xElemMarker[j];
    for (j=0; j<3; j++) xtet[1][j] = xtet[0][j];
    xtet[1][2] -= 1.5*a; xtet[1][0] -= .866*a; xtet[1][1] -= 0.5*a;
    for (j=0; j<3; j++) xtet[2][j] = xtet[0][j];
    xtet[2][2] -= 1.5*a; xtet[2][0] += .866*a; xtet[2][1] -= 0.5*a;
    for (j=0; j<3; j++) xtet[3][j] = xtet[0][j];
    xtet[3][2] -= 1.5*a; xtet[3][1] += a;
    for (j=0; j<6; j++)
      glVertex3f(xtet[itet[j]][0], xtet[itet[j]][1], xtet[itet[j]][2]);
    glEnd();
  }

  // Draw in negative volume points points
  if (nNegVolPoints > 0){
    a = ModelSize*.001;
    glColor3f(1.0, 0.0, 0.0);
    
    for (i=0; i<nNegVolPoints; i++){
      glBegin(GL_TRIANGLE_STRIP);
      for (j=0; j<3; j++) xtet[0][j] = xNegVolPoints[3*i+j];
      for (j=0; j<3; j++) xtet[1][j] = xtet[0][j];
      xtet[1][2] -= 1.5*a; xtet[1][0] -= .866*a; xtet[1][1] -= 0.5*a;
      for (j=0; j<3; j++) xtet[2][j] = xtet[0][j];
      xtet[2][2] -= 1.5*a; xtet[2][0] += .866*a; xtet[2][1] -= 0.5*a;
      for (j=0; j<3; j++) xtet[3][j] = xtet[0][j];
      xtet[3][2] -= 1.5*a; xtet[3][1] += a;
      for (j=0; j<6; j++)
	glVertex3f(xtet[itet[j]][0], xtet[itet[j]][1], xtet[itet[j]][2]);
      glEnd();
    }
  }

  // draw contours if on
  if ((ContoursOn) && (All->Mesh->Dim == 3)) xf_DrawContours3D();

  // draw cut plane if on
  if ((iDGroupCutPlane >= 0) && (DGroup[iDGroupCutPlane].Active))
    xf_DrawDGroup(iDGroupCutPlane);

  // draw isosurface if on
  if ((iDGroupIsoSurf >= 0) && (DGroup[iDGroupIsoSurf].Active))
    xf_DrawDGroup(iDGroupIsoSurf);


  glutSwapBuffers();

}


/******************************************************************/
//   FUNCTION Definition: xf_DisplayAux
void 
xf_DisplayAux(void)
{
  int k, CurrentDScalar;
  char title[80], s0[80], s1[80];
  float ytop, dy;
  float x0, x1, y0, y1;
  float dcbar = 0.1;
  xf_DisplayVector *DV;
  xf_DisplayScalar *DS;
  real *c;
  
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity();
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // turn off color mappings
  //glPixelTransferf(GL_MAP_COLOR, GL_TRUE);

  ytop = WinAux.Center[1] + WinAux.Width[1]*0.5;
  
  // Border
  glLineWidth(6.0);
  x0 = WinAux.Center[0] - WinAux.Width[0]*0.5;
  x1 = WinAux.Center[0] + WinAux.Width[0]*0.5;
  y0 = WinAux.Center[1] - WinAux.Width[1]*0.5;
  y1 = WinAux.Center[1] + WinAux.Width[1]*0.5;
  glBegin(GL_LINE_STRIP);
  glColor3f(0.5, 0.5, 0.5);
  glVertex3f(x0, y1, 0.0);
  glVertex3f(x0, y0, 0.0);
  glVertex3f(x1, y0, 0.0);
  glEnd();
    
  xf_DrawText(0.05, ytop-.07, "XFlow plotter", 14);
  xf_DrawText(0.05, ytop-.15, "Left mouse button = menu", 24);
  
  ytop -= 0.3;

  if (RenderFlag){
    glColor3f(1.0, 1.0, 1.0);
    DV = DVector + CurrentDVector;
    DS = DV->Scalar + DV->CurrentDScalar;
    strncpy(s0, DV->DataName, 12); s0[12] = '\0';
    if (DV->SetIndex >= 0){
      strncpy(s1, DV->DataName, 8); s1[8] = '\0';
      sprintf(s0, "%s(%d)", s1, DV->SetIndex);
    }
    strncpy(s1, DS->Name, 12); s1[12] = '\0';
    sprintf(title, "[%s] %s", s0, s1);
    xf_DrawText(0.05, ytop-0.05, title, 30);
    sprintf(title, "%.4g", DS->Vmax);
    xf_DrawText(0.21, ytop-dcbar-.01, title, 11);
    sprintf(title, "%.4g", DS->Vmin);
    xf_DrawText(0.21, ytop-dcbar-1.0, title, 11);

    dy = 1.0/(xfp_CBAR_SIZE);
    
    // Colorbar
    glBegin(GL_QUADS);
    for (k=0; k<xfp_CBAR_SIZE; k++){
      y0 = ytop-dcbar- k   *dy;
      y1 = ytop-dcbar-(k+1)*dy;
      
      c = CBarColor + 3*(xfp_CBAR_SIZE-k);
      glColor3f(c[0], c[1], c[2]);
      glVertex3f(0.05, y0, 0.0);
      glVertex3f(0.20, y0, 0.0);
      
      c = CBarColor + 3*(xfp_CBAR_SIZE-k-1);
      glColor3f(c[0], c[1], c[2]);
      glVertex3f(0.20, y1, 0.0);
      glVertex3f(0.05, y1, 0.0);
    }
    glEnd();
  }
 

  glutSwapBuffers();

}


/******************************************************************/
//   FUNCTION Definition: xf_DisplayBorder
void 
xf_DisplayBorder(void)
{
  float x0, x1, y0, y1;
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // turn off color mappings
  glPixelTransferf(GL_MAP_COLOR, GL_TRUE);

  // border
  x0 = WinBorder.Center[0] - WinBorder.Width[0]*0.5;
  x1 = WinBorder.Center[0] + WinBorder.Width[0]*0.5;
  y0 = WinBorder.Center[1] - WinBorder.Width[1]*0.5;
  y1 = WinBorder.Center[1] + WinBorder.Width[1]*0.5;
  glBegin(GL_QUADS);
  glColor3f(.5, .5, .5);
  glVertex3f(x0, y0, 0.0);
  glVertex3f(x1, y0, 0.0);
  glVertex3f(x1, y1, 0.0);
  glVertex3f(x0, y1, 0.0);
  glEnd();
 
  glutSwapBuffers();
}



/******************************************************************/
//   FUNCTION Definition: xf_GrabPixels
static int
xf_GrabPixels(enum xfe_Bool inColor, unsigned int width, unsigned int height,
	      GLvoid **pbuffer)
{
  // Taken from code by Mark J. Kilgard, 1996
  GLvoid *buffer;
  GLint swapbytes, lsbfirst, rowlength;
  GLint skiprows, skippixels, alignment;
  GLenum format;
  unsigned int size;

  if (inColor) {
    format = GL_RGB;
    size = width * height * 3;
  } else {
    format = GL_RED;
    size = width * height * 1;
  }

  buffer = (GLvoid *) malloc(size*sizeof(float));
  if (buffer == NULL) xf_Error(xf_MEMORY_ERROR);

  /* Save current modes. */
  glGetIntegerv(GL_UNPACK_SWAP_BYTES, &swapbytes);
  glGetIntegerv(GL_UNPACK_LSB_FIRST, &lsbfirst);
  glGetIntegerv(GL_UNPACK_ROW_LENGTH, &rowlength);
  glGetIntegerv(GL_UNPACK_SKIP_ROWS, &skiprows);
  glGetIntegerv(GL_UNPACK_SKIP_PIXELS, &skippixels);
  glGetIntegerv(GL_UNPACK_ALIGNMENT, &alignment);
  /* Little endian machines (DEC Alpha for example) could
     benefit from setting GL_UNPACK_LSB_FIRST to GL_TRUE
     instead of GL_FALSE, but this would require changing the
     generated bitmaps too. */
  glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  /* Actually read the pixels. */
  glReadPixels(0, 0, width, height, format,
	       GL_FLOAT, (GLvoid *) buffer);


  /* Restore saved modes. */
  glPixelStorei(GL_UNPACK_SWAP_BYTES, swapbytes);
  glPixelStorei(GL_UNPACK_LSB_FIRST, lsbfirst);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, rowlength);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, skiprows);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, skippixels);
  glPixelStorei(GL_UNPACK_ALIGNMENT, alignment);

  (*pbuffer) = buffer;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_HardcopyWindowEPS
static int
xf_HardcopyWindowEPS(xf_WindowType Win, char *filename, enum xfe_Bool inColor)
{
  // Taken from code by Mark J. Kilgard, 1996
  FILE *fp;
  int ierr, width, height, orig;
  int components, pos, i;
  char cp;
  float *floatpix, f;
  GLvoid *pixels;

  width  = Win.PixWidth[0];
  height = Win.PixWidth[1];

  // set current window
  orig = glutGetWindow();
  glutSetWindow(Win.Handle);
  
  ierr = xf_Error(xf_GrabPixels(inColor, width, height, &pixels));
  if (ierr != xf_OK) return ierr;
   
  if (orig > 0) glutSetWindow(orig);
  
  if (inColor)
    components = 3;     /* Red, green, blue. */
  else
    components = 1;     /* Luminance. */

  if ((fp = fopen(filename, "w")) == NULL) return xf_Error(xf_FILE_WRITE_ERROR);
  fprintf(fp, "%%!PS-Adobe-2.0 EPSF-1.2\n");
  fprintf(fp, "%%%%Creator: OpenGL pixmap render output\n");
  fprintf(fp, "%%%%BoundingBox: 0 0 %d %d\n", width, height);
  fprintf(fp, "%%%%EndComments\n");
  fprintf(fp, "gsave\n");
  fprintf(fp, "/bwproc {\n");
  fprintf(fp, "    rgbproc\n");
  fprintf(fp, "    dup length 3 idiv string 0 3 0\n");
  fprintf(fp, "    5 -1 roll {\n");
  fprintf(fp, "    add 2 1 roll 1 sub dup 0 eq\n");
  fprintf(fp, "    { pop 3 idiv 3 -1 roll dup 4 -1 roll dup\n");
  fprintf(fp, "        3 1 roll 5 -1 roll put 1 add 3 0 }\n");
  fprintf(fp, "    { 2 1 roll } ifelse\n");
  fprintf(fp, "    } forall\n");
  fprintf(fp, "    pop pop pop\n");
  fprintf(fp, "} def\n");
  fprintf(fp, "systemdict /colorimage known not {\n");
  fprintf(fp, "    /colorimage {\n");
  fprintf(fp, "        pop\n");
  fprintf(fp, "        pop\n");
  fprintf(fp, "        /rgbproc exch def\n");
  fprintf(fp, "        { bwproc } image\n");
  fprintf(fp, "    } def\n");
  fprintf(fp, "} if\n");
  fprintf(fp, "/picstr %d string def\n", width * components);
  fprintf(fp, "%d %d scale\n", width, height);
  fprintf(fp, "%d %d %d\n", width, height, 8);
  fprintf(fp, "[%d 0 0 %d 0 0]\n", width, height);
  fprintf(fp, "{currentfile picstr readhexstring pop}\n");
  fprintf(fp, "false %d\n", components);
  fprintf(fp, "colorimage\n");


  floatpix = (float *) pixels;
  pos = 0;
  for (i = width * height * components; i > 0; i--) {
    f = *floatpix++;
    cp = f*255.0;
    fprintf(fp, "%02hhx", cp);
    if (++pos >= 32) {
      fprintf(fp, "\n");
      pos = 0;
    }
  }
  if (pos)
    fprintf(fp, "\n");

  fprintf(fp, "grestore\n");
  free(pixels);
  fclose(fp);
  return 0;
}


/******************************************************************/
//   FUNCTION Definition: xf_Idle
void 
xf_Idle() 
{
  int k;
  // glutPostRedisplay();
  if (nModifyMenu > 0){
    for (k=0; k<nModifyMenu; k++){
      switch (ModifyMenuList[k]){
      case xfe_PlotMenu_Sub_LoadVector  : xf_ModifySubMenu_LoadVector();   break;
      case xfe_PlotMenu_Sub_ChangeScalar: xf_ModifySubMenu_ChangeScalar(); break;
      case xfe_PlotMenu_Sub_Boundaries  : xf_ModifySubMenu_Boundaries();   break;
      case xfe_PlotMenu_Sub_Mesh        : xf_ModifySubMenu_Mesh();         break;
      case xfe_PlotMenu_Sub_ColorScheme : xf_ModifySubMenu_ColorScheme();  break;
      case xfe_PlotMenu_Sub_CutPlane    : xf_ModifySubMenu_CutPlane();     break;
      case xfe_PlotMenu_Sub_LineProbe   : xf_ModifySubMenu_LineProbe();    break;
      case xfe_PlotMenu_Sub_Contours    : xf_ModifySubMenu_Contours();     break;
      case xfe_PlotMenu_Sub_IsoSurface  : xf_ModifySubMenu_IsoSurface();   break;
      case xfe_PlotMenu_Sub_Lighting    : xf_ModifySubMenu_Lighting();     break;
      default: break;
      }
    } // k
    nModifyMenu = 0;
  }
}

/******************************************************************/
//   FUNCTION Definition: xf_SetChildrenWindowSizes
static void
xf_SetChildrenWindowSizes()
{
  int hh, ww, w, h;
  int bordersize = 4;

  w = WinParent.PixWidth[0];
  h = WinParent.PixWidth[1];
  ww = w/4;
  hh = h/3;
 
  WinBig->PixWidth[0] = w-ww-bordersize;
  WinBig->PixWidth[1] = h;
  WinBig->PixTopLeft[0] = 0;
  WinBig->PixTopLeft[1] = 0;

  WinSmall->PixWidth[0] = ww;
  WinSmall->PixWidth[1] = hh;
  WinSmall->PixTopLeft[0] = w-ww;
  WinSmall->PixTopLeft[1] = h-hh;

  WinAux.PixWidth[0] = ww;
  WinAux.PixWidth[1] = h-hh;
  WinAux.PixTopLeft[0] = w-ww;
  WinAux.PixTopLeft[1] = 0;

  WinBorder.PixWidth[0] = bordersize;
  WinBorder.PixWidth[1] = h;
  WinBorder.PixTopLeft[0] = w-ww-bordersize;
  WinBorder.PixTopLeft[1] = 0;

}

/******************************************************************/
//   FUNCTION Definition: xf_SetProjection2D
static void 
xf_SetProjection2D(xf_WindowType Win) 
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(Win.Center[0]-0.5*Win.Width[0],
	     Win.Center[0]+0.5*Win.Width[0],
	     Win.Center[1]-0.5*Win.Width[1],
	     Win.Center[1]+0.5*Win.Width[1]);

  glMatrixMode(GL_MODELVIEW);
}


/******************************************************************/
//   FUNCTION Definition: xf_SetProjection3D
static void 
xf_SetProjection3D(xf_WindowType Win) 
{
  real dist;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  dist = Camera.Distance;
  
  glOrtho(-0.5*Win.Width[0], 0.5*Win.Width[0],
	  -0.5*Win.Width[1], 0.5*Win.Width[1],
	  -10.0*(ModelSize+dist), 10.0*(ModelSize+dist));

 /*  glFrustum(-0.5*Win.Width[0], 0.5*Win.Width[0],  */
/* 	     -0.5*Win.Width[1], 0.5*Win.Width[1],  */
/* 	     0.1*(ModelSize+dist), 10.0*(ModelSize+dist)); */

  glMatrixMode(GL_MODELVIEW);
}


/******************************************************************/
//   FUNCTION Definition: xf_Resize
void 
xf_Resize(int w, int h) 
{
  int orig;
  real wmax, pw, ph;
  
  /* This function is called for the parent window; children sizes
     need to be set appropriately here. */

  WinParent.PixWidth[0] = w;
  WinParent.PixWidth[1] = h;
  
  xf_SetChildrenWindowSizes();

  orig = glutGetWindow(); // get current window (to set it back later)
  
  // 2D window
  glutSetWindow(Win2D.Handle);
  glViewport(0, 0, Win2D.PixWidth[0], Win2D.PixWidth[1]);
  pw = Win2D.PixWidth[0]; ph = Win2D.PixWidth[1];
  wmax = max(Win2D.Width[0], Win2D.Width[1]);
  if (pw > ph){
    Win2D.Width[1] = (ph*wmax/pw);
    Win2D.Width[0] = wmax;
  }
  else{
    Win2D.Width[0] = (pw*wmax/ph);
    Win2D.Width[1] = wmax;
  }
  glutPositionWindow(Win2D.PixTopLeft[0], Win2D.PixTopLeft[1]);
  glutReshapeWindow(Win2D.PixWidth[0], Win2D.PixWidth[1]);
  xf_SetProjection2D(Win2D);
  glutPostRedisplay();
  
  // 3D window
  glutSetWindow(Win3D.Handle);
  glViewport(0, 0, Win3D.PixWidth[0], Win3D.PixWidth[1]);
  wmax = max(Win3D.Width[0], Win3D.Width[1]);
  pw = Win3D.PixWidth[0]; ph = Win3D.PixWidth[1];
  if (pw > ph){
    Win3D.Width[1] = (ph*wmax/pw);
    Win3D.Width[0] = wmax;
  }
  else{
    Win3D.Width[0] = (pw*wmax/ph);
    Win3D.Width[1] = wmax;
  }
  glutPositionWindow(Win3D.PixTopLeft[0], Win3D.PixTopLeft[1]);
  glutReshapeWindow(Win3D.PixWidth[0], Win3D.PixWidth[1]);
  xf_SetProjection3D(Win3D);
  glutPostRedisplay();

  // Aux window
  glutSetWindow(WinAux.Handle);
  glViewport(0, 0, WinAux.PixWidth[0], WinAux.PixWidth[1]);
  glutPositionWindow(WinAux.PixTopLeft[0], WinAux.PixTopLeft[1]);
  glutReshapeWindow(WinAux.PixWidth[0], WinAux.PixWidth[1]);
  pw = WinAux.PixWidth[0]; ph = WinAux.PixWidth[1];
  WinAux.Width[0] = 1.0;
  WinAux.Width[1] = ph/pw;
  WinAux.Center[0] = WinAux.Width[0]/2.0;
  WinAux.Center[1] = WinAux.Width[1]/2.0;
  xf_SetProjection2D(WinAux);
  glutPostRedisplay();

  // Border window
  glutSetWindow(WinBorder.Handle);
  glViewport(0, 0, WinBorder.PixWidth[0], WinBorder.PixWidth[1]);
  glutPositionWindow(WinBorder.PixTopLeft[0], WinBorder.PixTopLeft[1]);
  glutReshapeWindow(WinBorder.PixWidth[0], WinBorder.PixWidth[1]);
  pw = WinBorder.PixWidth[0]; ph = WinBorder.PixWidth[1];
  WinBorder.Width[0] = 1.0;
  WinBorder.Width[1] = ph/pw;
  WinBorder.Center[0] = WinBorder.Width[0]/2.0;
  WinBorder.Center[1] = WinBorder.Width[1]/2.0;
  xf_SetProjection2D(WinBorder);
  glutPostRedisplay();

  // back to original
  if (orig > 0) glutSetWindow(orig);
 
}

/******************************************************************/
//   FUNCTION Definition: Calculate2DWindowSize
static void
xf_Calculate2DWindowSize()
{
  int i;
  real x, y, minx, miny, maxx, maxy;
  real *xp, pw, ph, wmax;

  // calculate clipping planes for 2D viewport
  if ((iDGroupCutPlane >= 0) && (DGroup[iDGroupCutPlane].Active)
      && (DGroup[iDGroupCutPlane].DTri->nNode > 0)){
    minx = miny =  1e10;
    maxx = maxy = -1e10;
    if (All->Mesh->Dim == 2) 
      xp = DGroup[iDGroupCutPlane].DTri->xglob;
    else
      xp = DGroup[iDGroupCutPlane].DTri->xplane;
    for (i=0; i<DGroup[iDGroupCutPlane].DTri->nNode; i++){
      x = xp[3*i+0];
      y = xp[3*i+1];
      minx = min(minx, x); miny = min(miny, y);
      maxx = max(maxx, x); maxy = max(maxy, y);
    } // i
    
  }
  else{
    minx = 0.0; maxx = 1.0;
    miny = 0.0; maxy = 1.0;
  }

  Win2D.Center[0] = 0.5*(minx + maxx);
  Win2D.Center[1] = 0.5*(miny + maxy);

  Win2D.Width[0] = maxx - minx;
  Win2D.Width[1] = maxy - miny;

  pw = Win2D.PixWidth[0]; ph = Win2D.PixWidth[1];
  wmax = max(Win2D.Width[0], Win2D.Width[1]);
  if (pw > ph){
    Win2D.Width[1] = (ph*wmax/pw);
    Win2D.Width[0] = wmax;
  }
  else{
    Win2D.Width[0] = (pw*wmax/ph);
    Win2D.Width[1] = wmax;
  }

}


/******************************************************************/
//   FUNCTION Definition: xf_CalculateWindowSizes
static int 
xf_CalculateWindowSizes()
{

  // Set suggestion for parent
  WinParent.PixWidth[0] = 1000; WinParent.PixTopLeft[0] = 50;
  WinParent.PixWidth[1] =  750; WinParent.PixTopLeft[1] = 50;

  // set children pixel widths
  if (All->Mesh->Dim == 2){
    WinBig = &Win2D;   WinSmall = &Win3D;
  }
  else{
    WinBig = &Win3D;   WinSmall = &Win2D;
  }
    
  xf_SetChildrenWindowSizes();
 
  xf_Calculate2DWindowSize();

  // 3D window
  Win3D.Center[0] = 0.0;
  Win3D.Center[1] = 0.0;

  Win3D.Width[0] = 2*Camera.Distance;
  Win3D.Width[1] = 2*Camera.Distance*((real) Win3D.PixWidth[1])/((real) Win3D.PixWidth[0]);

  // Aux window
  WinAux.Width[0] = 1.0;
  WinAux.Width[1] = ((real) WinAux.PixWidth[1])/((real) WinAux.PixWidth[0]);
  WinAux.Center[0] = WinAux.Width[0]/2.0;
  WinAux.Center[0] = WinAux.Width[1]/2.0;

  // Border window
  WinBorder.Width[0] = 1.0;
  WinBorder.Width[1] = ((real) WinBorder.PixWidth[1])/((real) WinBorder.PixWidth[0]);
  WinBorder.Center[0] = WinBorder.Width[0]/2.0;
  WinBorder.Center[0] = WinBorder.Width[1]/2.0;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_InitCamera
static void 
xf_InitCamera()
{
  int d;
  Camera.ViewAxis = 2;
  
  Camera.Distance = ModelSize + ModelOffset;
  Camera.Azimuth = 0.0;
  Camera.Elevation = 0.0;
  Camera.Twist = 0.0;
  for (d=0; d<3; d++) Camera.Origin[d] = 0.0;
}


/******************************************************************/
//   FUNCTION Definition: xf_InitContours
static void 
xf_InitContours()
{
  Contours.nSeg  = 0;
  Contours.nSeg0 = 0;
  Contours.dSeg  = 10;
  Contours.xseg  = NULL;
  Contours.gxseg = NULL;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateDTri
static int 
xf_CreateDTri(xf_DisplayTri **pDTri)
{
  int ierr;

  ierr = xf_Error(xf_Alloc((void **) pDTri, 1, sizeof(xf_DisplayTri)));
  if (ierr != xf_OK) return ierr;
  
  (*pDTri)->Version = 0;
  (*pDTri)->nNode  = 0;
  (*pDTri)->nNode0 = 0;
  (*pDTri)->DNode  = 10;
  (*pDTri)->xglob  = NULL;
  (*pDTri)->xref   = NULL;
  (*pDTri)->xplane = NULL;
  (*pDTri)->Xglob  = NULL;
  (*pDTri)->nTri   = 0;
  (*pDTri)->nTri0  = 0;
  (*pDTri)->DTri   = 10;
  (*pDTri)->TriList  = NULL;
  (*pDTri)->nNormal  = 0;
  (*pDTri)->Normals  = NULL;
  (*pDTri)->Active   = NULL;
  (*pDTri)->Tri2Elem = NULL;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CreateDLoop
static int 
xf_CreateDLoop(xf_DisplayLoop **pDLoop)
{
  int ierr;

  ierr = xf_Error(xf_Alloc((void **) pDLoop, 1, sizeof(xf_DisplayLoop)));
  if (ierr != xf_OK) return ierr;
  
  (*pDLoop)->nLoop     = 0;
  (*pDLoop)->nNode     = NULL;
  (*pDLoop)->NodeList  = NULL;
  (*pDLoop)->Active    = NULL;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateDStrip
static int 
xf_CreateDStrip(xf_DisplayStrip **pDStrip)
{
  int ierr;

  ierr = xf_Error(xf_Alloc((void **) pDStrip, 1, sizeof(xf_DisplayStrip)));
  if (ierr != xf_OK) return ierr;
  
  (*pDStrip)->nStrip    = 0;
  (*pDStrip)->nNode     = NULL;
  (*pDStrip)->NodeList  = NULL;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CheckAllocDTri
static int 
xf_CheckAllocDTri(xf_DisplayTri *DTri, int nnode, int ntri,
		  enum xfe_Bool IsCutPlane)
{
  int ierr;
  int nNode0, nTri0;
  
  if ( (DTri->nNode+nnode) > DTri->nNode0){
    nNode0 = max(DTri->nNode+nnode, DTri->nNode0 + DTri->DNode);
    DTri->nNode0 = nNode0;
    ierr = xf_Error(xf_ReAlloc((void **) &DTri->xglob, 3*nNode0, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc((void **) &DTri->xref , 3*nNode0, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    if (IsCutPlane){
      ierr = xf_Error(xf_ReAlloc((void **) &DTri->xplane , 3*nNode0, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
  }

  if ( (DTri->nTri+ntri) > DTri->nTri0){
    nTri0 = max(DTri->nTri+ntri, DTri->nTri0 + DTri->DTri);
    ierr = xf_Error(xf_ReAlloc((void **) &DTri->TriList , 3*nTri0, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc((void **) &DTri->Tri2Elem, 2*nTri0, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc((void **) &DTri->Active, nTri0, sizeof(int)));
    if (ierr != xf_OK) return ierr; 
    DTri->nTri0 = nTri0;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_InterpolateScalar
static int 
xf_InterpolateScalar( const xf_DisplayVector *DV, int egrp, int elem, 
		      enum xfe_BasisType Basis, int Order, int sr, 
		      const real *U, xf_EqnSet *EqnSet,
		      int nq, real *xref,  xf_BasisData **pPhiData,
		      xf_JacobianData **pJData, const int *IParam, 
		      const real *RParam, real *v, real *gv, real *w)
{
  int ierr, i, d, e, de, ed, nn, dim;
  real *point_w, qgU[3*54], gu[9], S;
  real rho, p, Vel[3], Y[50], gamma, *pv;
  enum xfe_Bool MotionOn = xfe_False;
  enum xfe_Bool showp = xfe_False;
  enum xfe_Bool showstate = xfe_False;
  xf_BasisData *PhiData;
  const char *s = NULL;

  dim = All->Mesh->Dim;
  
  ierr = xf_Error(xf_EvalBasis(Basis, Order, xfe_True, nq, xref, 
			       xfb_Phi | xfb_GPhi | xfb_gPhi, pPhiData));
  if (ierr != xf_OK) return ierr;
  
  PhiData = (*pPhiData);
  nn = PhiData->nn;
  MotionOn = ((All->Mesh->Motion != NULL) && (MeshMotionInPhysical));
  
  //Determine if showing orders; if so, do not transform when mesh motion is active
  s = DV->DataName;
  if ((strlen(s) > 5) && (strncmp(s+strlen(s)-5,"Order",5) == 0)){
    showp = xfe_True;
  }
  
  //Determine if showing state; if so, transform when mesh motion is active
  if ((strlen(s) >= 5) && (strncmp(s,"State",5) == 0)){
    showstate = xfe_True;
  }

  // MDI = motion data structure for interpolation
  if (MotionOn){
    if (MDI == NULL){
      ierr = xf_Error(xf_CreateMotionData(All, &MDI));
      if (ierr != xf_OK) return ierr;
    }
    // realloc MDIxglob if necessary
    if (nq > MDInq){
      MDInq = nq;
      ierr = xf_Error(xf_ReAlloc((void **) &MDIxglob, MDInq*dim, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
    // obtain global coords
    ierr = xf_Error(xf_Ref2GlobElem(All->Mesh, egrp, elem, NULL, xfe_True, 
				    nq, xref, MDIxglob));
    if (ierr != xf_OK) return ierr;
    // apply mesh motion map
    ierr = xf_Error(xf_MeshMotionMap( egrp, elem, PhiData, All->Mesh->Motion, nq, dim, 
				      Time, MDIxglob, MDI));
    if (ierr != xf_OK) return ierr;
  }

  // interpolate v
  xf_MxM_Set(PhiData->Phi, U, nq, nn, sr, v);
  
  // if motion on, set v = v/g to obtain the physical state
  if (MotionOn && (!showp) && showstate) xf_ColDiv(v, MDI->gb, nq, sr, 1);

  // interpolate gv, if necessary
  if (gv != NULL){
    /* Compute geometry Jacobian */
    ierr = xf_ElemJacobian(All->Mesh, egrp, elem, nq, xref, xfb_detJ | xfb_iJ, 
			   xfe_True, pJData);
    //if (ierr != xf_OK) return ierr; // ignore zero/negative Jacobian errors
    
    // convert reference basis grads (GPhi) to physical grads, gPhi
    ierr = xf_Error(xf_EvalPhysicalGrad(PhiData, (*pJData)));
    if (ierr != xf_OK) return ierr;

    for (d=0; d<All->Mesh->Dim; d++)
      xf_MxM_Set(PhiData->gPhi+nn*nq*d, U, nq, nn, sr, gv+nq*sr*d);
    
    // Modify state gradient in presence of mesh motion
    if (MotionOn)
      xf_ModMotionPhysGrad(nq, dim, sr, MDI, v, gv);
    
  }
 
  if (DV->DeriveType == xfe_VectorDeriveNone){
    if (DV->nScalar != sr) return xf_Error(xf_NOT_SUPPORTED);
    for (i=0; i<nq*sr; i++) w[i] = v[i];
  }
  else if (DV->DeriveType == xfe_VectorDeriveScalar){
    ierr = xf_Error(Yu_VisualizationOutput(sr, dim, nq, v, gv, w));
    if (ierr != xf_OK) return ierr;
    
  }
  else if (DV->DeriveType == xfe_VectorDeriveVariableSet){
    ierr = xf_Error(xf_EqnSetVariableChange(EqnSet, DV->DataName, IParam, RParam, 
					    nq, v, w, NULL, NULL));
    if (ierr != xf_OK) return ierr;
  }
  else return xf_Error(xf_NOT_SUPPORTED);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CalculateDGroupElems
static int 
xf_CalculateDGroupElems(xf_DisplayGroup *DG)
{
  int ierr, k, i, d, dim, s;
  int egrp, elem, ielem, CurrentDScalar;
  int ibfgrp, ibface, nbface;
  int nnode, prefine, refine;
  int Order;
  enum xfe_Bool PointsChanged, hit;
  enum xfe_ShapeType Shape, pShape;
  int *IParam = NULL;
  real *RParam = NULL;
  real *coord, *xglob, *v, *gv = NULL, *w, *EV, val;
  xf_BFace BFace;
  xf_BasisData *PhiData, *VPhiData;
  xf_JacobianData *JData;
  xf_Vector *V;
  xf_DisplayVector *DV;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;

  Mesh = All->Mesh;
  dim = Mesh->Dim;
  EqnSet = All->EqnSet;


  
  // if calculating an iso-surface, need a scalar
  V = NULL;
  if (DG->Type == xfe_DGIsoSurf){
    if ((nDVector <= 0) || (CurrentDVector >= nDVector))
      return xf_Error(xf_INPUT_ERROR);
    DV = DVector + CurrentDVector;
    V = DV->Vector;
    CurrentDScalar = DV->CurrentDScalar;
    // pull off IParam and RParam if need them
    if ((DV->DeriveType == xfe_VectorDeriveScalar) ||
	(DV->DeriveType == xfe_VectorDeriveVariableSet)){
      ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
      if (ierr != xf_OK) return ierr;
    }
  }


  if (((DG->Type == xfe_DGCutPlane) || (DG->Type == xfe_DGIsoSurf))
      && (dim == 3)){ // cut plane or iso-surface in 3D

    for (egrp=0, DG->nelem=0; egrp < Mesh->nElemGroup; egrp++)
      DG->nelem += Mesh->ElemGroup[egrp].nElem;    

    ierr = xf_Error(xf_ReAlloc((void **) &DG->Egrp, DG->nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_ReAlloc((void **) &DG->Elem, DG->nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    ielem = 0;
    pShape   = -1;
    PhiData  = NULL;
    VPhiData = NULL;
    JData    = NULL;
    coord    = NULL;
    xglob    = NULL;
    v        = NULL;
    for (egrp=0, ielem=0; egrp < Mesh->nElemGroup; egrp++){
      // pull off element shape
      ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
      if (ierr != xf_OK) return ierr;

      prefine = -1;

      for (elem=0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
	
	refine = Elem2Refine[egrp][elem];

	// get refined coords
	if (refine != prefine){
	  ierr = xf_Error(xf_GetRefineCoords(Shape, refine, &nnode, &coord, 
					     NULL, NULL, NULL, NULL));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_ReAlloc((void **) &xglob, 3*nnode, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  if (DG->Type == xfe_DGIsoSurf){
	    ierr = xf_Error(xf_ReAlloc((void **) &v, 5*nnode*max(V->StateRank, DV->nScalar), 
				       sizeof(real)));
	    if (ierr != xf_OK) return ierr;
	    w  = v +   nnode*max(V->StateRank, DV->nScalar);
	    gv = v + 2*nnode*max(V->StateRank, DV->nScalar); 
	  }
	  prefine = refine;
	}

	PointsChanged = ( (refine != prefine) || (Shape != pShape));

	if (DG->Type == xfe_DGIsoSurf){
	  // interpolate distance
	  if (V == NULL) return xf_Error(xf_CODE_LOGIC_ERROR);
	  EV = V->GenArray[egrp].rValue[elem];
	  
	  Order = xf_InterpOrder(V, egrp, elem);
	  
	  ierr = xf_Error(xf_InterpolateScalar(DV, egrp, elem, V->Basis[egrp], Order, 
					       V->StateRank, V->GenArray[egrp].rValue[elem], 
					       EqnSet, nnode, coord, &VPhiData, &JData,
					       IParam, RParam, v, gv, w));
	  if (ierr != xf_OK) return ierr;

	  // pick off the desired scalar
	  for (i=0; i<nnode; i++) 
	    v[i] = w[i*DV->nScalar + CurrentDScalar] - DG->isoval;
	}
	else{
	  ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &PhiData, PointsChanged,
					  nnode, coord, xglob));
	  if (ierr != xf_OK) return ierr;

	}
	
	// loop over sub-node coords and check for a sign change (or 0)
	hit = xfe_False;
	for (k=0; k<nnode; k++){
	  if (DG->Type == xfe_DGIsoSurf){
	    val = v[k];
	  }
	  else{
	    for (d=0, val=CutPlane[3]; d<3; d++) val += CutPlane[d]*xglob[3*k+d];
	  }
	  if (fabs(val) < MEPS){
	    hit = xfe_True;
	    break;
	  }
	  if (k == 0) s = sign(val);
	  else{
	    if (sign(val) != s){
	      hit = xfe_True;
	      break;
	    }
	  }
	} // k
	if (hit){
	  DG->Egrp[ielem] = egrp;
	  DG->Elem[ielem] = elem;
	  ielem++;
	}

      } // elem
      pShape = Shape;
    } // egrp
    
    xf_printf("Number of elements intersecting plane = %d\n", ielem);

    // set true number of intersected elements
    DG->nelem = ielem;

    // Resize vectors in DG
    ierr = xf_Error(xf_ReAlloc((void **) &DG->Egrp, DG->nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_ReAlloc((void **) &DG->Elem, DG->nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;


    /* Destroy Basis Data */
    ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
    if (ierr != xf_OK) return ierr;

    /* Destroy Basis Data for V interpolation*/
    ierr = xf_Error(xf_DestroyBasisData(VPhiData, xfe_True));
    if (ierr != xf_OK) return ierr;

    /* Destroy geometry Jacobian Data */
    ierr = xf_Error(xf_DestroyJacobianData(JData));
    if (ierr != xf_OK) return ierr;

    // release memory
    xf_Release( (void *) coord);
    xf_Release( (void *) xglob);
    xf_Release( (void *) v);
  }
  else if ((DG->Type == xfe_DGCutPlane) && (dim == 2)){ // cut plane in 2D

    for (egrp=0, DG->nelem=0; egrp < Mesh->nElemGroup; egrp++)
      DG->nelem += Mesh->ElemGroup[egrp].nElem;

    ierr = xf_Error(xf_ReAlloc((void **) &DG->Egrp, DG->nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_ReAlloc((void **) &DG->Elem, DG->nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    
    for (egrp=0, ielem=0; egrp < Mesh->nElemGroup; egrp++)
      for (elem=0; elem < Mesh->ElemGroup[egrp].nElem; elem++, ielem++){
	DG->Egrp[ielem] = egrp;
	DG->Elem[ielem] = elem;
      }
	
  }
  else{ // boundary face group
    ibfgrp = DG->ibfgrp;
    nbface = Mesh->BFaceGroup[ibfgrp].nBFace;


    DG->nelem = nbface;
    ierr = xf_Error(xf_ReAlloc((void **) &DG->Egrp, DG->nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_ReAlloc((void **) &DG->Elem, DG->nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    for (ibface=0; ibface<nbface; ibface++){
      BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
      DG->Egrp[ibface] = BFace.ElemGroup;
      DG->Elem[ibface] = BFace.Elem;
    }

  }

  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CutPlaneAxes
static void 
xf_CutPlaneAxes(real *N, real *up, real *right, real *pd)
{
  int k;
  real *n;
  real d, NN, upn;

  // determine plane normal
  n = CutPlane;
  for (k=0, NN=0.; k<3; k++) NN += n[k]*n[k];
  NN = sqrt(NN);
  for (k=0; k<3; k++) N[k] = n[k]/NN;

  // determine up axis
  if ((fabs(N[(UpAxis2D+1)%3]) < MEPS) && (fabs(N[(UpAxis2D+2)%3]) < MEPS))
    UpAxis2D = (UpAxis2D + 1)%3;
  
  // determine up and right vectors
  for (k=0; k<3; k++) up[k] = 0.;
  up[UpAxis2D] = 1.0;
  for (k=0, upn=0.; k<3; k++) upn += up[k]*N[k];
  for (k=0; k<3; k++) up[k] -= upn*N[k];
  for (k=0, NN=0.; k<3; k++) NN += up[k]*up[k]; // normalize up vector
  for (k=0, NN=sqrt(NN); k<3; k++) up[k] = up[k]/NN;

  xf_CrossProduct(up, N, right); // right = up x N

  // determine d = distance to cut plane from origin
  (*pd) = d = CutPlane[3];
  
}

/******************************************************************/
//   FUNCTION Definition: xf_CalculateCutPlaneCoords
static int 
xf_CalculateCutPlaneCoords(xf_DisplayGroup *DG)
{
  int i, k;
  real right[3], up[3], N[3];
  real d, dnu, dnr, xr, xu;
  real *xglob, *xplane;
  xf_DisplayTri *DTri;

  DTri = DG->DTri;

  xf_CutPlaneAxes(N, up, right, &d);

/*   // determine plane normal */
/*   n = CutPlane; */
/*   for (k=0, NN=0.; k<3; k++) NN += n[k]*n[k]; */
/*   NN = sqrt(NN); */
/*   for (k=0; k<3; k++) N[k] = n[k]/NN; */

/*   // determine up axis */
/*   if ((fabs(N[(UpAxis2D+1)%3]) < MEPS) && (fabs(N[(UpAxis2D+2)%3]) < MEPS)) */
/*     UpAxis2D = (UpAxis2D + 1)%3; */
  
/*   // determine up and right vectors */
/*   for (k=0; k<3; k++) up[k] = 0.; */
/*   up[UpAxis2D] = 1.0; */
/*   for (k=0, upn=0.; k<3; k++) upn += up[k]*N[k]; */
/*   for (k=0; k<3; k++) up[k] -= upn*N[k]; */
/*   for (k=0, NN=0.; k<3; k++) NN += up[k]*up[k]; // normalize up vector */
/*   for (k=0, NN=sqrt(NN); k<3; k++) up[k] = up[k]/NN; */

/*   xf_CrossProduct(up, N, right); // right = up x N */

  // determine (d*n).up and (d*n).right
  for (k=0, dnu=0.; k<3; k++) dnu += d*N[k]*up[k];
  for (k=0, dnr=0.; k<3; k++) dnr += d*N[k]*right[k];
  
  // transform coordinates
  for (i=0; i<DTri->nNode; i++){
    xglob = DTri->xglob + 3*i;
    xplane = DTri->xplane + 3*i;
    for (k=0, xr=0.; k<3; k++) xr += xglob[k]*right[k];
    for (k=0, xu=0.; k<3; k++) xu += xglob[k]*up[k];
    xplane[0] = xr + dnr;
    xplane[1] = xu + dnu;
    xplane[2] = 0.0;
  } // i

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: ClipCutPlane
static int
xf_ClipCutPlane()
{
  int ierr, nNode;
  int i, k, sg, nelemtot;
  int *nodeflag, *ev, *nElem, **ElemFlag;
  real x, y, minx, miny, maxx, maxy;
  real *xp;
  xf_DisplayTri *DTri;
  xf_DisplayLoop *DLoop;

  // return if nothing to do
  if ((iDGroupCutPlane < 0) || (!DGroup[iDGroupCutPlane].Active)
      || (DGroup[iDGroupCutPlane].DTri->nNode <= 0)){
    xf_printf("Nothing to clip.\n");
    return xf_OK;
  }

  nNode = DGroup[iDGroupCutPlane].DTri->nNode;

  ierr = xf_Error(xf_Alloc((void **) &nodeflag, nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<nNode; i++) nodeflag[i] = 1;

  if (All->Mesh->Dim == 2) 
    xp = DGroup[iDGroupCutPlane].DTri->xglob;
  else
    xp = DGroup[iDGroupCutPlane].DTri->xplane;

  minx = Win2D.Center[0] - 0.5*Win2D.Width[0];
  maxx = Win2D.Center[0] + 0.5*Win2D.Width[0];
  miny = Win2D.Center[1] - 0.5*Win2D.Width[1];
  maxy = Win2D.Center[1] + 0.5*Win2D.Width[1];

  for (i=0; i<nNode; i++){
    x = xp[3*i+0];
    y = xp[3*i+1];
    if ((x < minx) || (x > maxx) || (y < miny) || (y > maxy))
      nodeflag[i] = 0;
  } // i

  DTri  = DGroup[iDGroupCutPlane].DTri;
  DLoop = DGroup[iDGroupCutPlane].DLoop;

  // flag tris as active or not
  for (i=0; i<DTri->nTri; i++){
    DTri->Active[i] = 1;
    for (k=0; k<3; k++)
      if (nodeflag[DTri->TriList[3*i+k]] == 0) DTri->Active[i] = 0;
  } // i

  // clean up so we don't plot parts of elements
  ierr = xf_Error(xf_GetnElem(All->Mesh, &nElem, &nelemtot));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VAlloc2((void ***) &ElemFlag, All->Mesh->nElemGroup, 
			     nElem, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  for (i=0; i<nelemtot; i++) ElemFlag[0][i] = 1;

  for (i=0; i<DTri->nTri; i++){
    ev = DTri->Tri2Elem + 2*i;
    ElemFlag[ev[0]][ev[1]] = (ElemFlag[ev[0]][ev[1]] && DTri->Active[i]);
  }
  
  for (i=0; i<DTri->nTri; i++){
    ev = DTri->Tri2Elem + 2*i;
    DTri->Active[i] = ElemFlag[ev[0]][ev[1]];
  }

  xf_Release( (void *) nElem);
  xf_Release2( (void **) ElemFlag);

  // flag loops as active or not
  for (i=0; i<DLoop->nLoop; i++){
    sg = (DLoop->Active[i] > 0) ? 1 : -1;
    DLoop->Active[i] = 2*sg;
    for (k=0; k<DLoop->nNode[i]; k++)
      if (nodeflag[DLoop->NodeList[i][k]] == 0) DLoop->Active[i] = sg;
  } // i

  xf_Release( (void *) nodeflag);

  return xf_OK;

}



/******************************************************************/
//   FUNCTION Definition: xf_CalculateDGroup
static int 
xf_CalculateDGroup(xf_DisplayGroup *DG)
{
  int ierr, i, d, dim, n0, CurrentDScalar;
  int egrp, elem, face, nface0, ielem;
  int nelemPrev, nNode0, iNode, iTri;
  int refine, prefine, pface, nnode, nsplit, nbound;
  int nelemtot, iLoop;
  int *vsplit, *vbound, *itemp, **iitemp;
  int Order;
  enum xfe_Bool PointsChanged, doHills;
  enum xfe_ShapeType Shape, pShape;
  int *IParam = NULL;
  real *RParam = NULL;
  real *coord, *xref, *v = NULL, *gv = NULL, *w = NULL, *EV;
  xf_DisplayTri *DTri;
  xf_DisplayLoop *DLoop;
  xf_DisplayStrip *DStrip;
  xf_BasisData *PhiData, *VPhiData;
  xf_JacobianData *JData = NULL;
  xf_Vector *V;
  xf_DisplayVector *DV;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;

  Mesh = All->Mesh;
  dim = Mesh->Dim;
  EqnSet = All->EqnSet;

  // if calculating an iso-surface, need a scalar
  V = NULL;
  if (DG->Type == xfe_DGIsoSurf){
    if (nDVector <= 0) return xf_OK;
    if (CurrentDVector >= nDVector) return xf_OK;
    DV = DVector + CurrentDVector;
    V = DV->Vector;
    CurrentDScalar = DV->CurrentDScalar;
    // pull off IParam and RParam if need them
    if ((DV->DeriveType == xfe_VectorDeriveScalar) ||
	(DV->DeriveType == xfe_VectorDeriveVariableSet)){
      ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // draw hilly surface if certain conditions are met
  doHills = xfe_False;
  if ((dim == 2) && (nDVector > 0) && (CurrentDVector < nDVector) && (DrawingHills)){
    doHills = xfe_True;
    DV = DVector + CurrentDVector;
    V = DV->Vector;
    CurrentDScalar = DV->CurrentDScalar;
    // pull off IParam and RParam if need them
    if ((DV->DeriveType == xfe_VectorDeriveScalar) ||
	(DV->DeriveType == xfe_VectorDeriveVariableSet)){
      ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
      if (ierr != xf_OK) return ierr;
    }
  }

  // allocate and initialize DTri structure if necessary
  if (DG->DTri == NULL){
    ierr = xf_Error(xf_CreateDTri(&DG->DTri));
    if (ierr != xf_OK) return ierr;
  }
  DTri = DG->DTri;

  // allocate and initialize DLoop structure if necessary
  if (DG->DLoop == NULL){
    ierr = xf_Error(xf_CreateDLoop(&DG->DLoop));
    if (ierr != xf_OK) return ierr;
  }
  DLoop = DG->DLoop;

  nelemPrev = DG->nelem; // previous number of participating elems

  // calculate participating elements (Elem2Refine may affect this)
  ierr = xf_Error(xf_CalculateDGroupElems(DG));
  if (ierr != xf_OK) return ierr;

  // reallocate DLoop, nDNode, and DNode if necessary
  DG->DLoop->nLoop = DG->nelem;
  if (DG->nelem > nelemPrev){
    // DLoop
    ierr = xf_Error(xf_ReAlloc((void **) &DG->DLoop->nNode, DG->nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<DG->nelem; i++) DLoop->nNode[i] = 0;
    ierr = xf_Error(xf_ReAlloc((void **) &DG->DLoop->NodeList, DG->nelem, sizeof(int *)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<DG->nelem; i++) DG->DLoop->NodeList[i] = NULL;
    ierr = xf_Error(xf_ReAlloc((void **) &DG->DLoop->Active, DG->nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<DG->nelem; i++) DG->DLoop->Active[i] = 2;
    
    // nDNode
    ierr = xf_Error(xf_ReAlloc((void **) &DG->nDNode, DG->nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<DG->nelem; i++) DG->nDNode[i] = 0;

    // DNode
    ierr = xf_Error(xf_ReAlloc((void **) &DG->DNode, DG->nelem, sizeof(int * )));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<DG->nelem; i++) DG->DNode[i] = NULL;
  }


  // allocate and initialize DStrip structure if necessary
  if (DG->DStrip == NULL){
    ierr = xf_Error(xf_CreateDStrip(&DG->DStrip));
    if (ierr != xf_OK) return ierr;
  }
  DStrip = DG->DStrip;
  DStrip->nStrip = 0;

  // Allocate DStrip memory
  if (dim == 2){

    ierr = xf_Error(xf_ReAlloc( (void **) &DG->DStrip->nNode, 4+2*DG->nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    DG->DStrip->nStrip = 0;
    for (ielem=0; ielem<DG->nelem; ielem++){
      egrp = DG->Egrp[ielem];
      elem = DG->Elem[ielem];
      for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
	if (Mesh->ElemGroup[egrp].Face[elem][face].Group >= 0){
	  DG->DStrip->nNode[DG->DStrip->nStrip++] = Elem2Refine[egrp][elem]+1;
	}
      } // face
    } // ielem

    // reallocate DStrip nNode
    ierr = xf_Error(xf_ReAlloc((void **) &DG->DStrip->nNode, DG->DStrip->nStrip, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    if (DG->DStrip->nStrip > 0){

      // allocate DStrip NodeList
      xf_Release2( (void **) DG->DStrip->NodeList);
      ierr = xf_Error(xf_VAlloc2((void ***) &DG->DStrip->NodeList, DG->DStrip->nStrip, DG->DStrip->nNode, sizeof(int)));
      if (ierr != xf_OK) return ierr;

      // fill in nodes in main loop below
    }

    DG->DStrip->nStrip = 0;
    DStrip = DG->DStrip;
  } // if 2d



  // fill in DTri
  DTri->nNode = 0;
  DTri->nTri  = 0;
  coord    = NULL;
  vsplit   = NULL;
  vbound   = NULL;
  PhiData  = NULL;
  VPhiData = NULL;
  JData    = NULL;
  v       = NULL;
  iLoop   = 0;
  pShape  = -1;
  prefine = -1;
  pface   = -1;
  face    = -1;

  // Xglob needs to be deleted since nodes may change
  if (DTri->Xglob != NULL){
    xf_Release( (void *) DTri->Xglob);
    DTri->Xglob = NULL;
  }

  // loop over all participating elements
  for (ielem=0; ielem<DG->nelem; ielem++){
    egrp = DG->Egrp[ielem];
    elem = DG->Elem[ielem];

    // pull off element shape
    ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
    if (ierr != xf_OK) return ierr;

    refine = Elem2Refine[egrp][elem]; // element refinement level (for plotting)
    if (DG->ibfgrp >= 0) // pull off face number if calculating a bfg
      face = Mesh->BFaceGroup[DG->ibfgrp].BFace[ielem].Face;

    // determine whether reference display coordinates changed
    PointsChanged = ( (refine != prefine) || (Shape != pShape) || (pface != face) ||
		      ((dim == 3) && (DG->ibfgrp < 0)) );
    
    if (PointsChanged){ // calculate new coordinates
      prefine = refine;
      pShape = Shape;
      if (dim == 2){
	ierr = xf_Error(xf_GetRefineCoords(Shape, refine, &nnode, &coord, 
					   &nsplit, &vsplit, &nbound, &vbound));
	if (ierr != xf_OK) return ierr;

	// hills-specific
	if (doHills){
	  ierr = xf_Error(xf_ReAlloc( (void **) &v, 2*nnode*max(V->StateRank, DV->nScalar), 
				      sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  w  = v +   nnode*max(V->StateRank, DV->nScalar);
	}
      }
      else{
	if (DG->Type == xfe_DGBoundary){
	  // boundary group
	  pface = face;
	  ierr = xf_Error(xf_GetRefineCoordsOnFace(Shape, face, refine, &nnode, &coord, 
						   &nsplit, &vsplit, &nbound, &vbound));
	  if (ierr != xf_OK) return ierr;

	}
	else if (DG->Type == xfe_DGCutPlane){
	  // intersect elem with cut plane
	  ierr = xf_Error(xf_IntersectElemWithPlane(Mesh, egrp, elem, refine, CutPlane, NULL,
						    PointsChanged, &PhiData, &nnode, 
						    &coord, &nsplit, &vsplit, &nbound, 
						    &vbound));
	  if (ierr != xf_OK) return ierr;
	  if (nbound > 0) DLoop->Active[iLoop] = -2;
	}
	else if (DG->Type == xfe_DGIsoSurf){
	  // determine iso-surface triangles on elem
	  if (V == NULL) return xf_Error(xf_CODE_LOGIC_ERROR);
	  EV = V->GenArray[egrp].rValue[elem];

	  // get refined coords on elem
	  ierr = xf_Error(xf_GetRefineCoords(Shape, refine, &n0, &coord, 
					     NULL, NULL, NULL, NULL));
	  if (ierr != xf_OK) return ierr;

	  ierr = xf_Error(xf_ReAlloc( (void **) &v, 5*n0*max(V->StateRank, DV->nScalar), 
				      sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  w  = v +   n0*max(V->StateRank, DV->nScalar);
	  gv = v + 2*n0*max(V->StateRank, DV->nScalar); 

	  Order = xf_InterpOrder(V, egrp, elem);

	  ierr = xf_Error(xf_InterpolateScalar(DV, egrp, elem, V->Basis[egrp], Order, 
					       V->StateRank, V->GenArray[egrp].rValue[elem], 
					       EqnSet, n0, coord, &VPhiData, &JData, IParam, 
					       RParam, v, gv, w));
	  if (ierr != xf_OK) return ierr;
	  

	  // pick off the desired scalar
	  for (i=0; i<n0; i++) 
	    v[i] = w[i*DV->nScalar + CurrentDScalar] - DG->isoval;
	  
	  ierr = xf_Error(xf_IntersectElemWithPlane(Mesh, egrp, elem, refine, NULL, v,
						    PointsChanged, &PhiData, &nnode, 
						    &coord, &nsplit, &vsplit, NULL, 
						    NULL));
	  if (ierr != xf_OK) return ierr;
	  nbound = 0;
	}
      } // end else 3D
    }

    // hills-specific
    if (doHills){  
      Order = xf_InterpOrder(V, egrp, elem);
      
      ierr = xf_Error(xf_InterpolateScalar(DV, egrp, elem, V->Basis[egrp], Order, 
					   V->StateRank, V->GenArray[egrp].rValue[elem], 
					   EqnSet, nnode, coord, &VPhiData, &JData, IParam, 
					   RParam, v, NULL, w));
      if (ierr != xf_OK) return ierr;
      // pick off the desired scalar
      for (i=0; i<nnode; i++) v[i] = w[i*DV->nScalar + CurrentDScalar];
    }

    
    ierr = xf_Error(xf_CheckAllocDTri(DTri, nnode, nsplit, (DG->Type == xfe_DGCutPlane)));
    if (ierr != xf_OK) return ierr;

    nNode0 = DTri->nNode; // current number of nodes (before addition of nnode)
    
    ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &PhiData, PointsChanged, 
				    nnode, coord, DTri->xglob+3*nNode0));
    if (ierr != xf_OK) return ierr;

    for (i=0; i<nnode; i++){
      iNode = nNode0 + i;
      for (d=0; d<dim; d++) DTri->xref[3*iNode+d] = coord[dim*i+d];
      if (dim == 2) DTri->xref[3*iNode+2] = 0.0;
    }

    
    if (dim == 2){
      for (i=nnode-1; i>=0; i--){
	DTri->xglob[3*nNode0 + 3*i+2] = ((doHills) ? v[i] : 0.0);
	for (d=dim-1; d>=0; d--) 
	  DTri->xglob[3*nNode0 + 3*i+d] = DTri->xglob[3*nNode0 + 2*i+d];
      }
    }
     
	
    for (i=0; i<nsplit; i++){
      iTri = DTri->nTri+i;
      for (d=0; d<3; d++) DTri->TriList[3*iTri+d] = nNode0 + vsplit[3*i+d];
      DTri->Tri2Elem[2*DTri->nTri + 2*i + 0] = egrp;
      DTri->Tri2Elem[2*DTri->nTri + 2*i + 1] = elem;
      DTri->Active[DTri->nTri+i] = 1;
    }
    
    // Set DNode
    if (nnode != DG->nDNode[ielem]){
      DG->nDNode[ielem] = nnode;
      ierr = xf_Error(xf_ReAlloc((void **) &DG->DNode[ielem], nnode, sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
    for (i=0; i<nnode; i++) DG->DNode[ielem][i] = DTri->nNode + i;

    // Add loop of nodes for Mesh element plotting
    if (nbound != DLoop->nNode[iLoop]){
      DLoop->nNode[iLoop] = nbound;
      ierr = xf_Error(xf_ReAlloc((void **) &DLoop->NodeList[iLoop], nbound, sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
    for (i=0; i<nbound; i++) DLoop->NodeList[iLoop][i] = DTri->nNode + vbound[i];
    if (nbound > 0) iLoop++;

    // Add to DStrip for boundary plotting
    if ((dim == 2) && (nbound > 0)){

      // number of original faces in Shape
      ierr = xf_Error(xf_Shape2nFace(Shape, &nface0));
      if (ierr != xf_OK) return ierr;

      for (face=0; face<nface0; face++){
	if (Mesh->ElemGroup[egrp].Face[elem][face].Group >= 0){
	  for (i=0; i<(nbound/nface0+1); i++) 
	    DStrip->NodeList[DStrip->nStrip][i] = DTri->nNode + vbound[(nbound/nface0*face+i)%nbound];
	  DStrip->nStrip++;
	}
      } // face
    }
    
    DTri->nNode += nnode;
    DTri->nTri  += nsplit;
  } // ielem


  xf_printf("Display stats: nTri = %d, nNode = %d\n", DTri->nTri, DTri->nNode);

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Data for V interpolation */
  ierr = xf_Error(xf_DestroyBasisData(VPhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;

  // Release memory
  xf_Release((void *) coord);
  xf_Release((void *) vsplit);
  xf_Release((void *) vbound);
  xf_Release((void *) v);
  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);

  
  /* Create cut-plane coords if this is a cut plane in 3D*/
  if ((DG->Type == xfe_DGCutPlane) && (dim == 3)){
    ierr = xf_Error(xf_CalculateCutPlaneCoords(DG));
    if (ierr != xf_OK) return ierr;
  }

  
  /* Increment Version*/
  DTri->Version++;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MoveDGroup
//This subroutine is for moving mesh
static int 
xf_MoveDGroup(xf_DisplayGroup *DG)
{
  int ierr, nNode, i;
  xf_DisplayTri *DTri = NULL;
  xf_MeshMotion *Motion;
  xf_MotionData MData;

  DTri = DG->DTri;

  if (DTri == NULL){ xf_printf("Warning, DTri not allocated. Not moving.\n"); return xf_OK;}
  if ((Motion=All->Mesh->Motion) == NULL){
    xf_printf("All->Mesh->Motion does not exist. Not moving.\n");
    return xf_OK;
  }

  nNode = DTri->nNode;

  // allocate Xglob if necessary (first call), and copy over xglob
  if (DTri->Xglob == NULL){
    ierr = xf_Error(xf_Alloc((void **) &DTri->Xglob, 3*nNode, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<3*nNode; i++) DTri->Xglob[i] = DTri->xglob[i];
  }

  // specify what we want from motion
  xf_InitMotionData(&MData);
  MData.npoint = nNode;
  MData.dim    = 3;
  MData.x      = DTri->xglob;

  ierr = xf_Error(xf_MeshMotionMap(-1, -1, NULL, Motion, nNode, 3, Time, DTri->Xglob, &MData));
  if (ierr != xf_OK){
    xf_printf("Error during mesh motion.  Not moving mesh.\n");
    for (i=0; i<3*nNode; i++) DTri->xglob[i] = DTri->Xglob[i]; // revert back
  }


  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_SetActiveDGroups
static void
xf_SetActiveDGroups()
{
  int i;

  nActiveDGroup = 0;

  for (i=0; i<nDGroup; i++)
    if (DGroup[i].Active)
      ActiveDGroup[nActiveDGroup++] = i;
}

/******************************************************************/
//   FUNCTION Definition: xf_CalculateActiveDGroups
static int 
xf_CalculateActiveDGroups()
{
  int ierr, i;

  xf_SetActiveDGroups();

  xf_printf("nActiveDGroup = %d\n", nActiveDGroup);

  for (i=0; i<nActiveDGroup; i++){
    ierr = xf_Error(xf_CalculateDGroup(DGroup + ActiveDGroup[i]));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MovePlotMesh
static int 
xf_MovePlotMesh(enum xfe_Bool RequestTime)
{
  int ierr, i;
  char buf[xf_MAXSTRLEN];

  // return immediately if mesh motion is not on
  if (All->Mesh->Motion == NULL) return xf_OK;

  xf_SetActiveDGroups();

  if (RequestTime){
    xf_printf("Enter Time:\n");
    if ((fgets(buf, xf_MAXSTRLEN, stdin) == NULL) ||
	(sscanf(buf, "%lf", &Time) != 1)){
      xf_printf("Not understood.\n");
      return xf_OK;
    }
  }

  // move display groups
  for (i=0; i<nActiveDGroup; i++){
    ierr = xf_Error(xf_MoveDGroup(DGroup + ActiveDGroup[i]));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_InitElem2Refine
static int 
xf_InitElem2Refine()
{
  int ierr, refine;
  int egrp, elem;
  int *nElem;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  if (Elem2Refine != NULL) return xf_Error(xf_INPUT_ERROR);

  ierr = xf_Error(xf_GetnElem(Mesh, &nElem, NULL));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_VAlloc2((void ***) &Elem2Refine, Mesh->nElemGroup, 
			     nElem, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    refine = max(Mesh->ElemGroup[egrp].QOrder, 1);
    //testing
    refine += 2;
    
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
      Elem2Refine[egrp][elem] = refine;
  }

  xf_Release( (void  *) nElem);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CalculateModelSize
static int 
xf_CalculateModelSize()
{
  int ierr, dim, d, i;
  real *x, offset;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  dim = Mesh->Dim;

  for (d=0; d<dim; d++){
    ModelBBox[2*d+0] =  1e10;
    ModelBBox[2*d+1] = -1e10;
  }
  
  for (i=0; i<Mesh->nNode; i++){
    x = Mesh->Coord[i];
    for (d=0; d<dim; d++){
      ModelBBox[2*d+0] = min(ModelBBox[2*d+0], x[d]);
      ModelBBox[2*d+1] = max(ModelBBox[2*d+1], x[d]);
    }
  }

  ModelSize = 0.;
  for (d=0; d<dim; d++) 
    ModelSize = max(ModelSize, ModelBBox[2*d+1]-ModelBBox[2*d+0]);
  
  
  ModelOffset = 0.;
  for (d=0; d<dim; d++){
    offset = max(ModelBBox[2*d+0], -ModelBBox[2*d+1]);
    ModelOffset = max(offset, ModelOffset);
  }

  // set CutPlane using bounding box
  d = 1;
  CutPlane[3] = 0.5*(ModelBBox[2*d+1]+ModelBBox[2*d+0]);

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_InitDGroup
static void
xf_InitDGroup(xf_DisplayGroup *DG, enum xfe_Bool Active,
	      enum xfe_DisplayGroupType Type, int ibfgrp)
{
  DG->Active = Active;
  DG->Type   = Type;
  DG->ibfgrp = ibfgrp;
  DG->isoval = 0.0;
  DG->DTri = NULL;
  DG->DLoop = NULL;
  DG->DStrip = NULL;

  DG->nelem = 0;
  DG->Egrp  = NULL;
  DG->Elem  = NULL;
  DG->nDNode = NULL; 
  DG->DNode  = NULL;
}


/******************************************************************/
//   FUNCTION Definition: xf_InitDGroups
static int 
xf_InitDGroups()
{
  int ierr, nbfgrp, dim, i;
  int iDGroup;
  int *itemp;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  dim = Mesh->Dim;
  nbfgrp = Mesh->nBFaceGroup;
  
  if (DGroup != NULL) return xf_Error(xf_INPUT_ERROR);

  // all bface groups (if dim == 3) + isosurf (if dim == 3) + a cut plane
  nDGroup = ((dim == 3) ? nbfgrp+2: 1);
  
  ierr = xf_Error(xf_Alloc((void **) &DGroup, nDGroup, sizeof(xf_DisplayGroup)));
  if (ierr != xf_OK) return ierr;

  iDGroup = 0;

  // CutPlane
  iDGroupCutPlane = iDGroup;
  xf_InitDGroup(DGroup + iDGroup, xfe_True, xfe_DGCutPlane, -1);
  iDGroup++;

  
  if (dim == 3){ // 3D-specific groups
  
    // allocate iDGroupBFG if dim == 3
    ierr = xf_Error(xf_Alloc((void **) &iDGroupBFG, nbfgrp, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<nbfgrp; i++){
      iDGroupBFG[i] = iDGroup;
      xf_InitDGroup(DGroup + iDGroup, xfe_True, xfe_DGBoundary, i);
      // inactivate if zero measure (based on title)
      if (strncmp(Mesh->BFaceGroup[i].Title, "ZeroMeasure", 11) == 0) 
	DGroup[iDGroup].Active = xfe_False;
      iDGroup++;
    }
    
    // IsoSurf
    iDGroupIsoSurf = iDGroup;
    xf_InitDGroup(DGroup + iDGroup, xfe_False, xfe_DGIsoSurf, -1);
    iDGroup++;

  }
  
  if (iDGroup != nDGroup) return xf_Error(xf_CODE_LOGIC_ERROR);
  
  // allocate memory for active group vector
  ierr = xf_Error(xf_Alloc((void **) &ActiveDGroup, nDGroup, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // set active groups
  xf_SetActiveDGroups();

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CheckAllocScalars
static int 
xf_CheckAllocScalars(xf_DisplayTri *DTri, xf_DisplayVector *DV, int iDGroup)
{
  int ierr, k;
  xf_DisplayScalar *DS;

  for (k=0; k<DV->nScalar; k++){
    DS = DV->Scalar + k;
    if (DS->Data[iDGroup].nNode < DTri->nNode){
      DS->Data[iDGroup].nNode = DTri->nNode;
      ierr = xf_Error(xf_ReAlloc((void **) &DS->Data[iDGroup].Value, 
				 DS->Data[iDGroup].nNode, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReAlloc((void **) &DS->Data[iDGroup].Color, 
				 3*DS->Data[iDGroup].nNode, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SetScalarColor
static int 
xf_SetScalarColor(xf_ScalarData *SD, real vmin, real vmax)
{
  int ierr, i;
  real val, col;
  real MinDiff;

  if (vmax < vmin) vmin = vmax;

  MinDiff = 1.0e-12; // minimum indistinguishable difference for plotting
  if ((vmax-vmin) < MinDiff){
    vmin -= MinDiff*0.5;
    vmax += MinDiff*0.5;
  }

  for (i=0; i<SD->nNode; i++){
    val = SD->Value[i];
    if (UseLogScale)
      col = (log(val)-log(vmin))/(log(vmax)-log(vmin));
    else
      col = (val-vmin)/(vmax-vmin);

    col = max(0.0, col);
    col = min(1.0, col);
    xf_Scalar2RGB(col, SD->Color+3*i);
  } // i

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_FillDScalars
static int 
xf_FillDScalars( xf_DisplayVector *DV, int iDGroup)
{
  int ierr, d, i, k, dim, sr;
  int egrp, elem, iDNode, nDNode, pnDNode;
  int ielem, Order;
  enum xfe_BasisType Basis;
  enum xfe_Bool Interpolated;
  int *IParam = NULL;
  real *RParam = NULL;
  real *xref, *v, *gv = NULL, *w;
  xf_DisplayTri *DTri;
  xf_DisplayGroup *DG;
  xf_BasisData *PhiData;
  xf_JacobianData *JData;
  xf_Vector *V;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;

  Mesh = All->Mesh;
  EqnSet  = All->EqnSet;
  dim = Mesh->Dim;

  V = DV->Vector;
  
  sr = V->StateRank;

  Interpolated = ((V->Basis != NULL) && (V->Order != NULL));
  
  DG = DGroup + iDGroup;
  DTri = DG->DTri;

  if ((DV->DeriveType == xfe_VectorDeriveScalar) ||
      (DV->DeriveType == xfe_VectorDeriveVariableSet)){
    ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
    if (ierr != xf_OK) return ierr;
  }

  // Make sure Scalars are up to date in size
  ierr = xf_Error(xf_CheckAllocScalars(DTri, DV, iDGroup));
  if (ierr != xf_OK) return ierr;

  // Set Values
  if (Interpolated){

    PhiData  = NULL;
    JData    = NULL;
    xref     = NULL;
    v        = NULL;
    w        = NULL;
    pnDNode  = -1;
    if (Mesh->nElemGroup != V->nArray) return xf_Error(xf_OUT_OF_BOUNDS);

    // check that mesh and vector are compatible
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
      if (Mesh->ElemGroup[egrp].nElem != V->GenArray[egrp].n) 
	return xf_Error(xf_OUT_OF_BOUNDS);
      if (V->GenArray[egrp].Size != xfe_SizeReal) 
	return xf_Error(xf_OUT_OF_BOUNDS);
    }
      
    // loop over all participating elements
    for (ielem=0; ielem<DG->nelem; ielem++){
      egrp = DG->Egrp[ielem];
      elem = DG->Elem[ielem];
    
      Basis = V->Basis[egrp];
      Order = xf_InterpOrder(V, egrp, elem);

      nDNode = DG->nDNode[ielem];
      ierr = xf_Error(xf_ReAlloc( (void **) &xref, dim*nDNode, sizeof(real)));
      if (ierr != xf_OK) return ierr;

      for (i=0; i<nDNode; i++){
	iDNode = DG->DNode[ielem][i];
	for (d=0; d<dim; d++)
	  xref[i*dim+d] = DTri->xref[iDNode*3+d];
      }

      if (nDNode > pnDNode){
	ierr = xf_Error(xf_ReAlloc( (void **) &v, 4*nDNode*max(DV->nScalar, V->StateRank), 
				    sizeof(real)));
	if (ierr != xf_OK) return ierr;
	gv = v + nDNode*max(DV->nScalar, V->StateRank);
	ierr = xf_Error(xf_ReAlloc( (void **) &w, nDNode*max(DV->nScalar, V->StateRank), 
				    sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }

      ierr = xf_Error(xf_InterpolateScalar(DV, egrp, elem, Basis, Order, sr, 
					   V->GenArray[egrp].rValue[elem],
					   EqnSet, nDNode, xref, &PhiData, &JData,
					   IParam, RParam, v, gv, w));
      if (ierr == xf_NON_PHYSICAL){
	xf_printf("Warning, non-physical state in elem=%d,%d. Setting scalars to zero.\n",
		  egrp, elem);
	for (k=0; k<nDNode*DV->nScalar; k++) w[k] = 0.;
      }
      else if (ierr != xf_OK) return ierr;
      
      for (i=0; i<nDNode; i++){
	iDNode = DG->DNode[ielem][i];
	for (k=0; k<DV->nScalar; k++){
	  DV->Scalar[k].Data[iDGroup].Value[iDNode] = w[DV->nScalar*i+k];
	}
      }
	
      pnDNode = nDNode;
    } // ielem
    
    /* Destroy Basis Data */
    ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
    if (ierr != xf_OK) return ierr;

    /* Destroy Basis Data */
    ierr = xf_Error(xf_DestroyJacobianData(JData));
    if (ierr != xf_OK) return ierr;
    
    xf_Release( (void *) v);
    xf_Release( (void *) w);
    xf_Release( (void *) xref);

  } // end if interpolated
  else
    return xf_Error(xf_NOT_SUPPORTED);

  // Set Version of each scalar
  for (k=0; k<DV->nScalar; k++)
    DV->Scalar[k].Data[iDGroup].Version = DTri->Version;

  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_FillDScalarsAll
//Do-all version of xf_FillDScalars for all ActiveDGroups
static int 
xf_FillDScalarsAll(xf_DisplayVector *DV, enum xfe_Bool KeepLimits)
{
  int ierr, i, ii, k;
  int iDGroup;
  real val, vmin, vmax;

  for (i=0; i<nActiveDGroup; i++){
    ierr = xf_Error(xf_FillDScalars(DV, ActiveDGroup[i]));
    if (ierr != xf_OK) return ierr;
  }

  // Find vmin, vmax
  for (k=0; k<DV->nScalar; k++){
    if ((!KeepLimits) || (DV->Scalar[k].Vmax <= DV->Scalar[k].Vmin)){
      vmin =  1e30;
      vmax = -1e30;
      for (ii=0; ii<nActiveDGroup; ii++){
	iDGroup = ActiveDGroup[ii];
	for (i=0; i<DV->Scalar[k].Data[iDGroup].nNode; i++){
	  val = DV->Scalar[k].Data[iDGroup].Value[i];
	  vmin = min(val, vmin);
	  vmax = max(val, vmax);
	} // i
      } // ii
      
      // Set Vmin, Vmax, and normalized Color
      DV->Scalar[k].Vmin = vmin;
      DV->Scalar[k].Vmax = vmax;
    }
    else{
      vmin = DV->Scalar[k].Vmin;
      vmax = DV->Scalar[k].Vmax;
    }
    //here set color for the on-show scalar
    for (ii=0; ii<nActiveDGroup; ii++){
      iDGroup = ActiveDGroup[ii];
      ierr = xf_Error(xf_SetScalarColor(DV->Scalar[k].Data + iDGroup, vmin, vmax));
      if (ierr != xf_OK) return ierr;
    }
  }

  // Calculate colors for Color Bar
  xf_FillColorBar();
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_Plane2GlobCoord
static int 
xf_Plane2GlobCoord(real *xplane, real *xglob)
{  
  int ierr, k;
  real A[9], *right, *up, *N, iA[9];
  real d, b[3];

  if (All->Mesh->Dim == 2){
    for (k=0; k<2; k++) xglob[k] = xplane[k];
    xglob[2] = 0.0;
  }
  else{
    // Solve for xglob from xplane [3x3 system]
    N = A; right = A+3; up = A+6;
    xf_CutPlaneAxes(N, up, right, &d);
    
    ierr = xf_Error(xf_MatDetInv(A, 3, NULL, iA));
    if (ierr != xf_OK) return ierr;
    
    b[0] = -d; // xglob dot N
    for (k=0, b[1]=xplane[0]; k<3; k++) b[1] += d*N[k]*right[k];
    for (k=0, b[2]=xplane[1]; k<3; k++) b[2] += d*N[k]*up[k];
    
    xf_MxV_Set(iA, b, 3, 3, xglob);
  }    

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_AddContourSeg
static int 
xf_AddContourSeg(real *xseg, real *gxseg)
{
  int ierr, iseg, k;
  
  iseg = Contours.nSeg++;

  if (Contours.nSeg >= Contours.nSeg0){
    Contours.nSeg0 = max(Contours.nSeg+1, Contours.nSeg0 + Contours.dSeg);
    Contours.dSeg *= 2;
    ierr = xf_Error(xf_ReAlloc((void **) &Contours.xseg, 6*Contours.nSeg0, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    if (All->Mesh->Dim == 3){
      ierr = xf_Error(xf_ReAlloc((void **) &Contours.gxseg, 6*Contours.nSeg0, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
  }

  for (k=0; k<6; k++) Contours.xseg[6*iseg+k] = xseg[k];
  if (All->Mesh->Dim == 3)
    for (k=0; k<6; k++) Contours.gxseg[6*iseg+k] = gxseg[k];

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FillContours
static int 
xf_FillContours()
{
  int ierr, iTri, i, k, d;
  int i0, i1, j0, j1, nseg;
  int CurrentDScalar;
  int *nvec;
  real *Value, *xplane;
  real Vmax, Vmin, dv, *xvec;
  real v, vmin, vmax, s, v0, v1;
  real xseg[9], N[3], gxseg[9];
  real offset = 1e-2;
  xf_DisplayTri *DTri;

  if (!RenderFlag)  return xf_INPUT_ERROR;

  if ((iDGroupCutPlane < 0) || (!DGroup[iDGroupCutPlane].Active)){
    xf_printf("Cut plane not on to fill contours.\n");
    return xf_INPUT_ERROR;
  }

  nContours = max(nContours, 2);

  Contours.nSeg = 0;

  CurrentDScalar = DVector[CurrentDVector].CurrentDScalar;
  Value = DVector[CurrentDVector].Scalar[CurrentDScalar].Data[iDGroupCutPlane].Value;
  DTri = DGroup[iDGroupCutPlane].DTri;
  Vmax = DVector[CurrentDVector].Scalar[CurrentDScalar].Vmax;
  Vmin = DVector[CurrentDVector].Scalar[CurrentDScalar].Vmin;
  dv = (Vmax - Vmin)/((real) (nContours-1));

  if (dv < MEPS) return xf_OK; // nothing to plot

  // pull off plane normal
  xf_CutPlaneAxes(N, xseg, xseg, &s);

  xplane = ((All->Mesh->Dim == 3) ? DTri->xplane : DTri->xglob);

  for (iTri=0; iTri<DTri->nTri; iTri++){
    
    nvec = DTri->TriList + 3*iTri;
    
    vmin = vmax = Value[nvec[0]];
    for (k=1; k<3; k++){
      vmin = min(vmin, Value[nvec[k]]);
      vmax = max(vmax, Value[nvec[k]]);
    }
    
    i0 = max((int) ((vmin-Vmin)/dv), 0) + 1;
    i1 = min((int) ((vmax-Vmin)/dv), nContours-1);
    
    // continue if no contours
    if (i1 < i0) continue;

    for (i=i0; i<=i1; i++){
      v = Vmin + dv*i;
      nseg = 0;
      for (k=0; k<3; k++){
	j0 = nvec[ k     ];
	j1 = nvec[(k+1)%3];
	v0 = Value[j0];
	v1 = Value[j1];
	if ((v<min(v0,v1)) || (v>max(v0,v1))) continue;
	s = ((v-v0)/(v1-v0));
	for (d=0; d<3; d++)
	  xseg[3*nseg+d] = (1-s)*xplane[3*j0+d] + s*xplane[3*j1+d];
	//	xseg[3*nseg+2] = 1.0; // to make contours visible
	nseg++;
      }
      if (nseg == 2){
	if (All->Mesh->Dim == 3){
	  ierr = xf_Error(xf_Plane2GlobCoord(xseg  , gxseg    ));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_Plane2GlobCoord(xseg+3, gxseg + 3));
	  if (ierr != xf_OK) return ierr;
	  for (k=0; k<6; k++) gxseg[k] -= offset*N[k%3];
	}

	ierr = xf_Error(xf_AddContourSeg(xseg, gxseg));
	if (ierr != xf_OK) return ierr;
      }
    }
  }

  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CreateDScalars
static int 
xf_CreateDScalars(int nScalar, char **Name, xf_DisplayScalar **pScalar)
{
  int ierr, i, j;

  ierr = xf_Error(xf_Alloc((void **) pScalar, nScalar, sizeof(xf_DisplayScalar)));
  if (ierr != xf_OK) return ierr;

  for (i=0; i<nScalar; i++){
    strcpy((*pScalar)[i].Name, "NULL");
    ierr = xf_Error(xf_Alloc((void **) &(*pScalar)[i].Data, nDGroup, 
			     sizeof(xf_ScalarData)));
    if (ierr != xf_OK) return ierr;
    for (j=0; j<nDGroup; j++){
      (*pScalar)[i].Data[j].Version = -1;
      (*pScalar)[i].Data[j].nNode  = 0;
      (*pScalar)[i].Data[j].Value  = NULL;
      (*pScalar)[i].Data[j].ValueX = NULL;
      (*pScalar)[i].Data[j].Color  = NULL;
    }
    (*pScalar)[i].Vmin =  0.0;
    (*pScalar)[i].Vmax = -1.0;

  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_InterpOrIntPrep
static int 
xf_InterpOrIntPrep(xf_Vector *V)
{
  int ierr;
  int i, in, ir;
  xf_GenArray *ga = NULL; 

  // if data is not interpolated set a default p=0 basis + order
  if ((V->Basis == NULL) || (V->Order == NULL)){
    xf_Release((void *) V->Basis);
    xf_Release((void *) V->Order);
    ierr = xf_Error(xf_Alloc((void **) &V->Basis, V->nArray, sizeof(enum xfe_BasisType)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &V->Order, V->nArray, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<V->nArray; i++){
      // irrelevant for p = 0
      V->Basis[i] = (All->Mesh->Dim == 3) ? xfe_TetLagrange : xfe_TriLagrange; 
      V->Order[i] = 0;
    }

  }

  // if data is integer size, make a real copy
  for (i=0; i<V->nArray; i++){
    if (V->GenArray[i].Size == xfe_SizeInt){
      ga = V->GenArray + i;
      ierr = xf_Error(xf_Alloc2((void ***) &(ga->rValue), ga->n, ga->r, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ga->Size = xfe_SizeReal;
      for (in=0; in<ga->n; in++)
	for (ir=0; ir<ga->r; ir++)
	  ga->rValue[in][ir] = (real) ga->iValue[in][ir];
      // no longer need integer part
      xf_Release2( (void **) ga->iValue);
      ga->iValue = NULL;
    }
  } // i

return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_LoadDVector
static int 
xf_LoadDVector(xf_Vector *V, const char *DataName, xf_DisplayVector **pDV)
{
  int ierr, i, k;
  xf_DisplayVector *DV;


  // reallocate DVector
  ierr = xf_Error(xf_ReAlloc((void **) &DVector, nDVector+1, sizeof(xf_DisplayVector)));
  if (ierr != xf_OK) return ierr;

  DV = DVector + nDVector;

  // copy DataName
  if (DataName != NULL){
    strncpy(DV->DataName, DataName, xf_MAXSTRLEN);
    DV->DataName[xf_MAXSTRLEN-1] = '\0';
  }
  else
    strcpy(DV->DataName, "NULL");

  // set default SetIndex
  DV->SetIndex = -1;

  // set pointer to V
  DV->Vector = V;

  // vector is not derived from others
  DV->DeriveType = xfe_VectorDeriveNone;

  // only support certain linkage types
  if (V->Linkage != xfe_LinkageGlobElem) return xf_Error(xf_NOT_SUPPORTED);

  // make sure we have data to plot
  if ((V->nArray <= 0) || (V->GenArray == NULL)){
    xf_printf("LoadDVector: no data found in the vector.\n");
    return xf_NOT_FOUND;
  }

  // Check/prepare integer or non-interpolated vector
  ierr = xf_Error(xf_InterpOrIntPrep(V));
  if (ierr != xf_OK) return ierr;
  
  // Determine number of scalars
  if ((V->Basis != NULL) && (V->Order != NULL)){
    // data is interpolated

    // determine number of scalars
    DV->nScalar = V->StateRank;

    // allocate scalars (do not fill in values)
    ierr = xf_Error(xf_CreateDScalars(DV->nScalar, V->StateName, &DV->Scalar));
    if (ierr != xf_OK) return ierr;
    
    // set names
    if (V->StateName != NULL){
      for (k=0; k<V->StateRank; k++){
	strncpy(DV->Scalar[k].Name, V->StateName[k], xf_MAXSTRLEN);
	DV->Scalar[k].Name[xf_MAXSTRLEN-1] = '\0';
      }
    }

    // set current scalar
    DV->CurrentDScalar = 0;
  }
  else{
    return xf_Error(xf_NOT_SUPPORTED);
  }

  if (pDV != NULL) (*pDV) = DV;

  nDVector++;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_LoadDVector_Scalars
static int 
xf_LoadDVector_Scalars(xf_Vector *V, const char *DataName, xf_DisplayVector **pDV)
{
  int ierr, i, k;
  int nQuant;
  int *IParam;
  real *RParam;
  char **QuantNames;
  xf_DisplayVector *DV;
  xf_EqnSet *EqnSet;

  EqnSet = All->EqnSet;
 
  // determine number of scalars
  //ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
  //if (ierr != xf_OK) return ierr;

  //ierr = xf_EqnSetScalar(EqnSet, NULL, IParam, RParam, 0, NULL, NULL, NULL,
  //			 NULL, &nQuant, &QuantNames, 0.0);
  nQuant = Yu_ScalarLast;
  QuantNames = Yu_ScalarName;
  
  
  //if ((ierr != xf_OK) || (nQuant == 0)){
  if ((nQuant == 0)){
    xf_printf("No EqnSetScalar function or no scalars returned.  Continuing.\n");
    return xf_OK;
  }

  // if data is not interpolated, error
  if ((V->Basis == NULL) || (V->Order == NULL)){
    xf_printf("Data is not interpolated.\n");
    return xf_Error(xf_NOT_SUPPORTED);
   }

  // reallocate DVector
  ierr = xf_Error(xf_ReAlloc((void **) &DVector, nDVector+1, sizeof(xf_DisplayVector)));
  if (ierr != xf_OK) return ierr;

  DV = DVector + nDVector;

  // copy DataName
  if (DataName != NULL){
    strncpy(DV->DataName, DataName, xf_MAXSTRLEN);
    DV->DataName[xf_MAXSTRLEN-1] = '\0';
  } 
  else
    strcpy(DV->DataName, "NULL");

  // set default SetIndex
  DV->SetIndex = -1;

  // set pointer to V
  DV->Vector = V;

  // vector is not derived from others
  DV->DeriveType = xfe_VectorDeriveScalar;

  // only support certain linkage types
  if (V->Linkage != xfe_LinkageGlobElem) return xf_Error(xf_NOT_SUPPORTED);

  // make sure we have data to plot
  if ((V->nArray <= 0) || (V->GenArray == NULL)){
    xf_printf("LoadDVector: no data found in the vector.\n");
    return xf_NOT_FOUND;
  }


  DV->nScalar = nQuant;

  // allocate scalars (do not fill in values)
  ierr = xf_Error(xf_CreateDScalars(DV->nScalar, QuantNames, &DV->Scalar));
  if (ierr != xf_OK) return ierr;
  
  // set names
  for (k=0; k<nQuant; k++){
    strncpy(DV->Scalar[k].Name, QuantNames[k], xf_MAXSTRLEN);
    xf_printf("Creating scalar %d = %s\n", k, QuantNames[k]);
    DV->Scalar[k].Name[xf_MAXSTRLEN-1] = '\0';
  }

  DV->CurrentDScalar = 0;

  if (pDV != NULL) (*pDV) = DV;

  nDVector++;

  //xf_Release( (void *) IParam);
  //xf_Release( (void *) RParam);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_LoadDVector_VariableSet
static int 
xf_LoadDVector_VariableSet(xf_Vector *V, const char *DataName, xf_DisplayVector **pDV)
{
  int ierr, i, k;
  char buf[xf_MAXSTRLEN];
  int nVSet, iVSet;
  int *IParam;
  real *RParam;
  char **VSetNames;
  xf_DisplayVector *DV;
  xf_EqnSet *EqnSet;

  EqnSet = All->EqnSet;
 
  // determine number of scalars
  ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;

  ierr = xf_EqnSetVariableChange(EqnSet, NULL, IParam, RParam, 0, NULL, NULL,
				 &nVSet, &VSetNames);
  if ((ierr != xf_OK) || (nVSet == 0)){
    xf_printf("No EqnSetVariableChange function or no new sets returned.  Continuing.\n");
    return xf_OK;
  }

  // pick index
  xf_printf("\n");
  xf_printf("Available variable sets:\n");
  for (k=0; k<nVSet; k++) xf_printf("%d : %s\n", k, VSetNames[k]);
  xf_printf("Choose variable set:\n");
  iVSet = -1;
  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
      (sscanf(buf, "%d", &k) == 1) && (k >= 0) && (k < nVSet)){
    xf_printf("Chose %s\n", VSetNames[k]);
    iVSet = k;
  }
  if (iVSet < 0){
    xf_printf("Not understood or out of range.\n");
    xf_Release( (void *) IParam);
    xf_Release( (void *) RParam);
    return xf_OK;
  }

  // if data is not interpolated, error
  if ((V->Basis == NULL) || (V->Order == NULL)){
    xf_printf("Data is not interpolated.\n");
    return xf_Error(xf_NOT_SUPPORTED);
  }

  // reallocate DVector
  ierr = xf_Error(xf_ReAlloc((void **) &DVector, nDVector+1, sizeof(xf_DisplayVector)));
  if (ierr != xf_OK) return ierr;

  DV = DVector + nDVector;

  // copy DataName
  strncpy(DV->DataName, VSetNames[iVSet], xf_MAXSTRLEN);
  DV->DataName[xf_MAXSTRLEN-1] = '\0';

  // set default SetIndex
  DV->SetIndex = -1;

  // set pointer to V
  DV->Vector = V;

  // vector is not derived from others
  DV->DeriveType = xfe_VectorDeriveVariableSet;

  // only support certain linkage types
  if (V->Linkage != xfe_LinkageGlobElem) return xf_Error(xf_NOT_SUPPORTED);

  // make sure we have data to plot
  if ((V->nArray <= 0) || (V->GenArray == NULL)){
    xf_printf("LoadDVector: no data found in the vector.\n");
    return xf_NOT_FOUND;
  }


  DV->nScalar = V->StateRank;

  // allocate scalars (do not fill in values)
  ierr = xf_Error(xf_CreateDScalars(DV->nScalar, V->StateName, &DV->Scalar));
  if (ierr != xf_OK) return ierr;
  
  // set names
  if (V->StateName != NULL){
    for (k=0; k<V->StateRank; k++){
      strncpy(DV->Scalar[k].Name, V->StateName[k], xf_MAXSTRLEN);
      DV->Scalar[k].Name[xf_MAXSTRLEN-1] = '\0';
    }
  }

  DV->CurrentDScalar = 0;

  if (pDV != NULL) (*pDV) = DV;

  nDVector++;

  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_BuildEqnSetQuant
static int 
xf_BuildEqnSetQuant(xf_All *All, xf_Vector **pV)
{
  // THIS FUNCTION IS NOT USED
  int ierr, iq, k, dim, nQuant = 0;
  int egrp, elem, nn;
  int pOrder, Order;
  char **QuantNames;
  int *IParam;
  real *RParam, *xn, *EU, *EV;
  real *u;
  xf_Vector *U, *V;
  xf_Data *StateData;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;

  Mesh    = All->Mesh;
  EqnSet  = All->EqnSet;
  dim = Mesh->Dim;

  ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;

  ierr = xf_EqnSetScalar(EqnSet, NULL, IParam, RParam, 0, NULL, NULL, NULL,
			 NULL, &nQuant, &QuantNames, 0.0);
  if ((ierr != xf_OK) || (nQuant == 0)){
    xf_printf("No EqnSetScalar function or no scalars returned.  Continuing.\n");
    return xf_OK;
  }
  
  // Get State, U
  ierr = xf_Error(xf_FindPrimalState(All->DataSet, 0, &StateData, NULL));
  if (ierr != xf_OK){
    xf_printf("Could not find primal state.  Continuing.\n");
    return ierr;
  }

  U = (xf_Vector *) StateData->Data;

  // Create a new vector for the quantities
  ierr = xf_Error(xf_FindVector(All, "EqnSetQuant", xfe_LinkageGlobElem, nQuant, QuantNames,
				0, 0, U->Basis, U->Order, U->nComp, U->vOrder, NULL, xfe_SizeReal, 
				xfe_False, xfe_True, NULL, pV, NULL));
  if (ierr != xf_OK) return ierr;

  V = (*pV);

/*   // copy over names of quantities */
/*   ierr = xf_Error(xf_Alloc2((void ***)&V->StateName, nQuant,  */
/* 			    xf_MAXSTRLEN, sizeof(char))); */
/*   if (ierr != xf_OK) return ierr; */
/*   for (k=0; k<nQuant; k++){ */
/*     xf_printf("k = %d, Name = %s\n", k, QuantNames[k]); */
/*     strcpy(V->StateName[k], QuantNames[k]); */
/*   } */
    
  // calculate quantities
  xn = NULL;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    
    pOrder = -1;

    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

      Order = xf_InterpOrder(U, egrp, elem);

      if (pOrder != Order){
	pOrder = Order;
	ierr = xf_Error(xf_Order2nNode(U->Basis[egrp], Order, &nn));
	if (ierr != xf_OK) return ierr;
	
	ierr = xf_Error(xf_LagrangeNodes(U->Basis[egrp], Order, NULL, NULL, &xn));
	if (ierr != xf_OK) return ierr;
      }

      EU = U->GenArray[egrp].rValue[elem];
      EV = V->GenArray[egrp].rValue[elem];

      ierr = xf_Error(xf_EqnSetScalar(EqnSet, NULL, IParam, RParam, nn, EU, NULL, EV,
				      NULL, NULL, NULL, 0.0));
      if (ierr != xf_OK) return ierr;
     
      for(iq=0; iq<nn; iq++)
      {
         u = EV + 11 * nn;
         u[0] = 1.2;
         u[1] = 1.3;
      }
      
    
    } // elem

  } // egrp


  // Release memory
  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);
  xf_Release( (void *) xn);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_LocateVectors
static int 
xf_LocateVectors(xf_All *All)
{
  int ierr, i;
  xf_Data *D;
  xf_Vector *V;
  xf_VectorSet *VS;
  xf_DataSet *DataSet;
  xf_EqnSet *EqnSet;
  xf_DisplayVector *DV;

  DataSet = All->DataSet;
  EqnSet  = All->EqnSet;
  
  D = DataSet->Head;

  while (D != NULL){
    if (D->Type == xfe_Vector){
      V = (xf_Vector *) D->Data;
      if (V->Linkage == xfe_LinkageGlobElem){
	ierr = xf_Error(xf_LoadDVector(V, D->Title, NULL));
	if (ierr != xf_OK) return ierr;
      }
    }
    else if (D->Type == xfe_VectorSet){
      VS = (xf_VectorSet *) D->Data;
      for (i=0; i<VS->nVector; i++){
	V = VS->Vector + i;
	if (V->Linkage == xfe_LinkageGlobElem){
	  ierr = xf_Error(xf_LoadDVector(V, D->Title, &DV));
	  if (ierr != xf_OK) return ierr;
	  DV->SetIndex = i;
	}
      }
    }
    D = D->Next;
  }

  /* Search for EqnSet scalars */
  HaveEqnSet = xfe_False;
  if (EqnSet->EqnSetLibrary != NULL){
    ierr = xf_LoadEqnSetLibrary(EqnSet->EqnSetLibrary);
    if (ierr != xf_OK){
      xf_printf("Warning, could not load EqnSet library %s. Continuing.\n", 
		EqnSet->EqnSetLibrary);
    }
    else{
      HaveEqnSet = xfe_True;
      ierr = xf_Error(xf_EqnSetRegister(EqnSet));
      if (ierr != xf_OK) return ierr;
      /*     ierr = xf_Error(xf_BuildEqnSetQuant(All, &V)); */
      /*     if (ierr != xf_OK) return ierr; */
      /*     ierr = xf_Error(xf_LoadDVector(V, "EqnSetQuant", NULL)); */
      /*     if (ierr != xf_OK) return ierr; */
    }
  }

  CurrentDVector = 0;

  xf_printf("nDVector = %d\n", nDVector);
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyDTri
static void
xf_DestroyDTri(xf_DisplayTri *DTri)
{
  if (DTri != NULL){
    xf_Release( (void *) DTri->xglob);
    xf_Release( (void *) DTri->xref);
    xf_Release( (void *) DTri->xplane);
    xf_Release( (void *) DTri->Xglob);
    xf_Release( (void *) DTri->TriList);
    xf_Release( (void *) DTri->Normals);
    xf_Release( (void *) DTri->Active);
    xf_Release( (void *) DTri->Tri2Elem);
    xf_Release( (void *) DTri);
  }
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyDLoop
static void
xf_DestroyDLoop(xf_DisplayLoop *DLoop)
{
  int i;
  if (DLoop != NULL){
    for (i=0; i<DLoop->nLoop; i++)
      xf_Release( (void *) DLoop->NodeList[i]);
    xf_Release( (void *) DLoop->nNode);
    xf_Release( (void *) DLoop->NodeList);
    xf_Release( (void *) DLoop->Active);
    xf_Release( (void *) DLoop);
  }
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyDStrip
static void
xf_DestroyDStrip(xf_DisplayStrip *DStrip)
{
  int i;
  if (DStrip != NULL){
    xf_Release2( (void **) DStrip->NodeList);
    xf_Release( (void *) DStrip->nNode);
    xf_Release( (void *) DStrip);
  }
}



/******************************************************************/
//   FUNCTION Definition: xf_DestroyDGroup
static void
xf_DestroyDGroup(xf_DisplayGroup *DG)
{
  int i;
  if (DG != NULL){
    xf_DestroyDTri(DG->DTri);
    xf_DestroyDLoop(DG->DLoop);
    xf_DestroyDStrip(DG->DStrip);
    xf_Release((void *) DG->Egrp);
    xf_Release((void *) DG->Elem);
    xf_Release((void *) DG->nDNode);
    for (i=0; i<DG->nelem; i++)
      xf_Release((void *) DG->DNode[i]);
    xf_Release((void *) DG->DNode);
  }
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyDScalar
static void
xf_DestroyDScalar(xf_DisplayScalar *DS)
{
  int j;
  if (DS != NULL){
    for (j=0; j<nDGroup; j++){
      xf_Release( (void *) DS->Data[j].Value);
      xf_Release( (void *) DS->Data[j].ValueX);
      xf_Release( (void *) DS->Data[j].Color);
    }
    xf_Release( (void *) DS->Data);
  }
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyDVector
static void
xf_DestroyDVector(xf_DisplayVector *DV)
{
  int k;
  if (DV != NULL){
    for (k=0; k<DV->nScalar; k++)
      xf_DestroyDScalar(DV->Scalar+k);
    xf_Release( (void *) DV->Scalar);
  }
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyGlobals
static void
xf_DestroyGlobals()
{
  int egrp, elem, i;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  for (i=0; i<nDGroup; i++) xf_DestroyDGroup(DGroup+i);
  xf_Release( (void *) DGroup);
  xf_Release( (void *) ActiveDGroup);
  xf_Release( (void *) iDGroupBFG);

  xf_Release2( (void **) Elem2Refine);

  for (i=0; i<nDVector; i++)
    xf_DestroyDVector(DVector+i);
  xf_Release( (void *) DVector);

  xf_Release( (void *) OutputIsActive);
  xf_Release( (void *) xOutputPoints);

  xf_Release( (void *) Contours.xseg);
  xf_Release( (void *) Contours.gxseg);

  // destroy windows
  glutDestroyWindow(Win2D.Handle); 
  glutDestroyWindow(Win3D.Handle); 
  glutDestroyWindow(WinAux.Handle); 
  glutDestroyWindow(WinBorder.Handle); 
  glutDestroyWindow(WinParent.Handle); 
}

/******************************************************************/
//   FUNCTION Definition: xf_ChooseScalar
static void 
xf_ChooseScalar () {
  int k;
  char buf[xf_MAXSTRLEN];
  xf_DisplayVector *DV;
  xf_DisplayScalar *DS;

  DV = DVector + CurrentDVector;

  xf_printf("\n");
  xf_printf("Available scalars for current vector [%s]:\n", DV->DataName);
  
  for (k=0; k<DV->nScalar; k++){
    if (DV->CurrentDScalar == k) xf_printf(" * ");
    else xf_printf("   ");
    xf_printf("%d : %s\n", k, DV->Scalar[k].Name);
  }

  xf_printf("Choose scalar for rendering:\n");
  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
      (sscanf(buf, "%d", &k) == 1)){
    if ((k >= 0) && (k < DV->nScalar)){
      DV->CurrentDScalar = k;
      xf_printf("Scalar set to %d [%s]\n", k, DV->Scalar[k].Name);
      if (k <= 9)
	xf_printf("Note, keys (0-9) can be typed directly in the 2D window.\n");
    }
    else
      xf_printf("Out of range.\n");
  }
  else
    xf_printf("Not understood.\n");

}

/******************************************************************/
//   FUNCTION Definition: xf_WriteScalarInfo
static void 
xf_WriteScalarInfo()
{
  char title[80], s0[80], s1[80];
  xf_DisplayVector *DV;
  xf_DisplayScalar *DS; 

  DV = DVector + CurrentDVector;
  DS = DV->Scalar + DV->CurrentDScalar;
  strncpy(s0, DV->DataName, 20); s0[20] = '\0';
  if (DV->SetIndex >= 0){
    strncpy(s1, DV->DataName, 10); s1[10] = '\0';
    sprintf(s0, "%s(%d)", s1, DV->SetIndex);
  }
  strncpy(s1, DS->Name, 20); s1[20] = '\0';
  sprintf(title, "[%s] %s", s0, s1);


  xf_printf("--- Colorbar Info ---\n");
  xf_printf("%s\n", title);
  xf_printf("Min value = %.10E\n", DS->Vmin);
  xf_printf("Max value = %.10E\n", DS->Vmax);
  xf_printf("---------------------\n");
}


/******************************************************************/
//   FUNCTION Definition: xf_GetVectorNames
static void 
xf_GetVectorNames(char ***pNames, int *pn){
  int k, ierr;

  /* Allocate pNames */
  ierr = xf_Error(xf_Alloc2((void ***) pNames, nDVector, xf_MAXSTRLEN, sizeof(char)));
  if (ierr != xf_OK) exit(0);

  /* fill pNames */
  for (k=0; k<nDVector; k++){
    sprintf((*pNames)[k], "%s", DVector[k].DataName);
    if (CurrentDVector == k) strcat((*pNames)[k], " *");
  }

  (*pn) = nDVector;
}


/******************************************************************/
//   FUNCTION Definition: xf_ListVectors
static void 
xf_ListVectors() 
{
  int k, k0, k1, iprev, i;
  
  xf_printf("\n");
  xf_printf("Available vector[sets]:\n");
  
  for (k=0; k<nDVector; k++){
    k0 = k;
    k1 = k;
    // the following assumes vectorset vectors are in order in DVector
    iprev = -1;
    while (((i=DVector[k1].SetIndex) >= 0) && 
	   (k1 < nDVector) && (i > iprev)){
      k1++;
      iprev = i;
    }
    if ((CurrentDVector == k) ||
	((CurrentDVector >= k0) && (CurrentDVector <= k1)))
      xf_printf(" * ");
    else xf_printf("   ");
    xf_printf("%d : %s", k, DVector[k].DataName);
    if (k1>k0) xf_printf("[%d]",k1-k0);
    xf_printf("\n");
    k = max(k, k1-1);
  }

}

/******************************************************************/
//   FUNCTION Definition: xf_ChooseVector
static void 
xf_ChooseVector () {
  int k, j;
  char buf[xf_MAXSTRLEN];

  xf_ListVectors();

  xf_printf("Choose vector for rendering:\n");
  j = 0;
  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
      ((sscanf(buf, "%d %d", &k, &j) == 2) || (sscanf(buf, "%d", &k) == 1))){
    k += j;
    if ((k >= 0) && (k < nDVector)){
      CurrentDVector = k;
      xf_printf("Vector set to %d: %s", k, DVector[k].DataName);
      if (DVector[k].SetIndex >= 0) xf_printf("(%d)", DVector[k].SetIndex);
      xf_printf("\n");
    }
    else
      xf_printf("Out of range.\n");
  }
  else
    xf_printf("Not understood.\n");
}


/******************************************************************/
//   FUNCTION Definition: xf_DeriveVector
static int 
xf_DeriveVector() 
{
  int ierr, k, j, iOrig;
  char buf[xf_MAXSTRLEN];
  char Title[xf_MAXSTRLEN];
  xf_Vector *U;

  //temp change
  /*
  if (!HaveEqnSet){
    xf_printf("Do not have equation set to perform vector derivations.\n");
    return xf_OK;
  }*/

  xf_ListVectors();

  U = NULL;
  xf_printf("Choose vector to derive from:\n");
  j = 0;
  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
      ((sscanf(buf, "%d %d", &k, &j) == 2) || (sscanf(buf, "%d", &k) == 1))){
    k += j;
    if ((k >= 0) && (k < nDVector)){
      xf_printf("Using Vector %s\n", DVector[k].DataName);
      iOrig = k;
      U = DVector[k].Vector;
    }
    else
      xf_printf("Out of range.\n");
  }
  else
    xf_printf("Not understood.\n");

  if (U != NULL){

    xf_printf("\n");
    xf_printf("Available derivations:\n");
    for (k=0; k<xfe_VectorDeriveLast; k++)
      xf_printf("%d : %s\n", k, xfe_VectorDeriveName[k]);

    xf_printf("Choose derivation:\n");
    k = -1;
    if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	((sscanf(buf, "%d", &k) == 1))){      
      if ((k >= 0) && (k < xfe_VectorDeriveLast)){
	xf_printf("Using derivation %s\n", xfe_VectorDeriveName[k]);
      }
      else
	xf_printf("Out of range.\n");
    }
    else
      xf_printf("Not understood.\n");
    
    if (k == xfe_VectorDeriveScalar){
      sprintf(Title, "%s_Scalars", DVector[iOrig].DataName);
      
      ierr = xf_Error(xf_LoadDVector_Scalars(U, Title, NULL));
      if (ierr != xf_OK) return ierr;
      
      CurrentDVector = nDVector-1; // set current vector to one just created
    }
    else if (k == xfe_VectorDeriveVariableSet){
      ierr = xf_Error(xf_LoadDVector_VariableSet(U,  DVector[iOrig].DataName, NULL));
      if (ierr != xf_OK) return ierr;
      
      CurrentDVector = nDVector-1; // set current vector to one just created
    }
    else{
      xf_printf("Not deriving.\n");
    }
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ShowInterpOrderVector
static int 
xf_ShowInterpOrderVector() 
{
  int ierr, k, j, iOrig;
  char buf[xf_MAXSTRLEN];
  char Title[xf_MAXSTRLEN];
  xf_Vector *U, *V;

  xf_ListVectors();

  U = NULL;
  xf_printf("Choose vector for which to display order information:\n");
  j = 0;
  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
      ((sscanf(buf, "%d %d", &k, &j) == 2) || (sscanf(buf, "%d", &k) == 1))){
    k += j;
    if ((k >= 0) && (k < nDVector)){
      xf_printf("Using Vector %s\n", DVector[k].DataName);
      iOrig = k;
      U = DVector[k].Vector;
    }
    else
      xf_printf("Out of range.\n");
  }
  else
    xf_printf("Not understood.\n");

  if (U != NULL){

    if ((U->nComp==NULL) || (U->vOrder==NULL)){
      xf_printf("%s does not have variable orders on an element level.\n", DVector[iOrig].DataName);
      if (U->Order != NULL) xf_printf("Order of first group = %d\n", U->Order[0]);
    }
    else{
    
      sprintf(Title, "%s_Order", DVector[iOrig].DataName);
	
      xf_printf("Orders are stored in %s.\n", Title);

      ierr = xf_Error(xf_BuildVOrder(All, U, &V));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_LoadDVector(V, Title, NULL));
      if (ierr != xf_OK) return ierr;
      
      CurrentDVector = nDVector-1; // set current vector to one just created
    }
  }
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_CompareVectors
static int 
xf_CompareVectors() 
{
  int ierr, k, j, i0, i1;
  int iq, nq, pnq, nScalar;
  int egrp, elem;
  char buf[xf_MAXSTRLEN];
  enum xfe_BasisType Basis0, Basis1;
  enum xfe_Bool QuadChanged;
  int Order0, Order1, QuadOrder;
  int *IParam = NULL;
  real temp, wtemp, err;
  real *RParam = NULL;
  real *xq = NULL;
  real *v0 = NULL, *v1 = NULL;
  real *w0, *w1;
  real *vError = NULL;
  xf_Vector *U0, *U1;
  xf_DisplayVector *DV0;
  xf_DisplayVector *DV1;
  xf_QuadData *QuadData = NULL;
  xf_JacobianData *JData = NULL;
  xf_BasisData  *PhiData0 = NULL,  *PhiData1 = NULL;
  xf_BasisData *vPhiData0 = NULL, *vPhiData1 = NULL;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;

  Mesh    = All->Mesh;
  EqnSet  = All->EqnSet;

  xf_ListVectors();

  U0 = NULL;
  xf_printf("Choose first vector for comparison:\n");
  j = 0;
  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
      ((sscanf(buf, "%d %d", &k, &j) == 2) || (sscanf(buf, "%d", &k) == 1))){
    k += j;
    if ((k >= 0) && (k < nDVector)){
      xf_printf("Using Vector %s\n", DVector[k].DataName);
      i0 = k;
      U0 = DVector[k].Vector;
    }
    else
      xf_printf("Out of range.\n");
  }
  else
    xf_printf("Not understood.\n");

  if (U0 == NULL){
    xf_printf("Vector is null!\n");
    return xf_OK;
  }

  xf_ListVectors();

  U1 = NULL;
  xf_printf("Choose second vector for comparison:\n");
  j = 0;
  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
      ((sscanf(buf, "%d %d", &k, &j) == 2) || (sscanf(buf, "%d", &k) == 1))){
    k += j;
    if ((k >= 0) && (k < nDVector)){
      xf_printf("Using Vector %s\n", DVector[k].DataName);
      i1 = k;
      U1 = DVector[k].Vector;
    }
    else
      xf_printf("Out of range.\n");
  }
  else
    xf_printf("Not understood.\n");

  if (U1 == NULL){
    xf_printf("Vector is null!\n");
    return xf_OK;
  }

  DV0 = DVector + i0;
  DV1 = DVector + i1;

  if (DV0->nScalar != DV1->nScalar){
    xf_printf("Vectors are not compatible in terms of # scalars.\n");
    return xf_OK;
  }

  nScalar = DV0->nScalar;
  ierr = xf_Error(xf_Alloc( (void **) &vError, nScalar, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  for (k=0; k<nScalar; k++) vError[k] = 0.;

  // retrieve function parameters
  ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;

  pnq = -1;

  // loop over element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    // Determine Basis and Order for U and Ue
    Basis0 = U0->Basis[egrp]; 
    Basis1 = U1->Basis[egrp]; 

    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

      Order0 = xf_InterpOrder(U0, egrp, elem);
      Order1 = xf_InterpOrder(U1, egrp, elem);

      // determine required integration order
      ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, 2*max(Order0, Order1), &QuadOrder));
      if (ierr != xf_OK) return ierr;
      
      /* Pull off quad points for the element; will not recalculate if
	 Basis/Order have not changed. */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;

      nq = QuadData->nquad;
      xq = QuadData->xquad;
      

      // compute basis functions
      ierr = xf_Error(xf_EvalBasis(Basis0, Order0, QuadChanged, nq, xq, xfb_Phi, &PhiData0));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_EvalBasis(Basis1, Order1, QuadChanged, nq, xq, xfb_Phi, &PhiData1));
      if (ierr != xf_OK) return ierr;

     
      /* Compute geometry Jacobian; if not constant, compute at quad
         points.  Note if jacobian is constant, only one Jacobian will
         be computed/returned. */
      ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ, QuadChanged, &JData));
      if (ierr != xf_OK) return ierr;

      // re-allocate data if quad points increased
      if (nq > pnq){
	ierr = xf_Error(xf_ReAlloc( (void **) &v0, 2*nq*max(DV0->nScalar, U0->StateRank), 
				    sizeof(real)));
	if (ierr != xf_OK) return ierr;
	w0 = v0 + nq*max(DV0->nScalar, U0->StateRank);
	ierr = xf_Error(xf_ReAlloc( (void **) &v1, 2*nq*max(DV1->nScalar, U1->StateRank), 
				    sizeof(real)));
	if (ierr != xf_OK) return ierr;
	w1 = v1 + nq*max(DV1->nScalar, U1->StateRank);
      }

      // interpolate scalars 0 and 1 at quad points
      ierr = xf_Error(xf_InterpolateScalar(DV0, egrp, elem, Basis0, Order0, U0->StateRank, 
					   U0->GenArray[egrp].rValue[elem],
					   EqnSet, nq, xq, &vPhiData0, NULL,
					   IParam, RParam, v0, NULL, w0));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_InterpolateScalar(DV1, egrp, elem, Basis1, Order1, U1->StateRank, 
					   U1->GenArray[egrp].rValue[elem],
					   EqnSet, nq, xq, &vPhiData1, NULL,
					   IParam, RParam, v1, NULL, w1));
      if (ierr != xf_OK) return ierr;

      // sum (w0-w1)*wq over quad points, add to vError
      for (iq=0; iq<nq; iq++){
	wtemp = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
	for (k=0; k<nScalar; k++){
	  temp = w0[iq*nScalar+k] - w1[iq*nScalar+k];
	  vError[k] += wtemp*temp*temp;
	}
      }

      pnq = nq;
      
    } // elem

  } // egrp

  // total L2 norm
  for (k=0, err=0.; k<nScalar; k++) err += vError[k];
  err = sqrt(err);

  // take sqrt of vError
  for (k=0; k<nScalar; k++) vError[k] = sqrt(vError[k]);
  
  // print out vError
  for (k=0; k<nScalar; k++)
    xf_printf("iScalar = %d: L2 Error = %.12E\n", k, vError[k]);
  xf_printf("Total error = %.10E\n\n", err);
  
  // Destroy QuadData
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData0, xfe_True));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyBasisData(PhiData1, xfe_True));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyBasisData(vPhiData0, xfe_True));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyBasisData(vPhiData1, xfe_True));
  if (ierr != xf_OK) return ierr;

   /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;

  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);

  xf_Release( (void *) v0);
  xf_Release( (void *) v1);
  xf_Release( (void *) vError);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ToggleBFG
static void 
xf_ToggleBFG() {
  int i, k, nbfgrp;
  enum xfe_Bool done;
  char c;
  char buf[xf_MAXSTRLEN];
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  nbfgrp = Mesh->nBFaceGroup;
  if (Mesh->Dim == 2) return; // nothing to do in 2D

  done = xfe_False;
  while (!done){
    xf_printf("\n");
    xf_printf("BFaceGroup list (* = active):\n");
    
    for (i=0; i<nbfgrp; i++){
      c = ((DGroup[iDGroupBFG[i]].Active) ? '*' : ' ');
      xf_printf("%d %c : %s\n", i, c, Mesh->BFaceGroup[i].Title);
    }
    
    xf_printf("Choose BFG # to toggle (-1 = all, enter = done):\n");

    if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	(sscanf(buf, "%d", &k) == 1)){
      if (k == -1){
	for (i=0; i<nbfgrp; i++) 
	  DGroup[iDGroupBFG[i]].Active = !DGroup[iDGroupBFG[i]].Active;
      }
      else if ((k >= 0) && (k < nbfgrp))
	DGroup[iDGroupBFG[k]].Active = !DGroup[iDGroupBFG[k]].Active;
    }
    else done = xfe_True;
  }
}



/******************************************************************/
//   FUNCTION Definition: xf_ToggleOutputs
static int
xf_ToggleOutputs() {
  int ierr, i, k, d, nOutput;
  enum xfe_Bool done;
  char c;
  char buf[xf_MAXSTRLEN];
  xf_EqnSet *EqnSet;
  xf_Output *Output;

  EqnSet = All->EqnSet;
  OutputsOn = !OutputsOn;


  if (!OutputsOn){
    nOutputPoints = 0;
    xf_printf("Turned Outputs Off.\n");
    return xf_OK;
  }

  if ((EqnSet->Outputs == NULL) || ((nOutput = EqnSet->Outputs->nOutput) == 0)){
    xf_printf("No Outputs exist.\n");
    return xf_OK;
  }

  Output = EqnSet->Outputs->Output;

  if (OutputIsActive == NULL){
    ierr = xf_Error(xf_Alloc((void **) &OutputIsActive, nOutput, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<nOutput; i++) OutputIsActive[i] = 0;
    ierr = xf_Error(xf_Alloc((void **) &xOutputPoints, 3*nOutput, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<3*nOutput; i++) xOutputPoints[i] = 0.0;
  }
  

  done = xfe_False;
  while (!done){
    xf_printf("\n");
    xf_printf("Output list (* = active):\n");
    
    for (i=0; i<nOutput; i++){
      c = ((OutputIsActive[i]) ? '*' : ' ');
      xf_printf("%d %c : %s\n", i, c, Output[i].Name);
    }
    
    xf_printf("Choose Output # to toggle (-1 = all, enter = done):\n");

    if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	(sscanf(buf, "%d", &k) == 1)){
      if (k == -1){
	for (i=0; i<nOutput; i++) 
	  OutputIsActive[i] = !OutputIsActive[i];
      }
      else if ((k >= 0) && (k < nOutput))
	OutputIsActive[k] = !OutputIsActive[k];
    }
    else done = xfe_True;
  }

  // Determine active points
  nOutputPoints = 0;
  xf_printf("% (Output index) (location)\n");
  for (i=0; i<nOutput; i++)
    if ((OutputIsActive[i]) && (Output[i].Type == xfe_PointValue)){
      ierr = xf_Error(xf_Ref2GlobElem(All->Mesh, Output[i].egrp, Output[i].elem, 
				      NULL, xfe_True, 1, Output[i].xref, 
				      xOutputPoints+3*nOutputPoints));
      if (ierr != xf_OK) return ierr;
      xf_printf("%d", i);
      for (d=0; d<All->Mesh->Dim; d++) xf_printf(" %.8E", xOutputPoints[3*nOutputPoints+d]);
      xf_printf("\n");
      nOutputPoints++;
    }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_IdentifyElement
static int
xf_IdentifyElement() {
  int ierr, egrp, elem, d;
  char buf[xf_MAXSTRLEN];
  real xref[3] = {0.3, 0.3, 0.3};
  real xref2[3] = {0.2, 0.2, 0.2};
  real xglob[3];

  xf_printf("Which element to identify? [two numbers: egrp elem]:\n");
  ElemMarkerOn = xfe_False;
  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
      (sscanf(buf, "%d %d", &egrp, &elem) == 2)){
    if ((egrp >= 0) && (egrp < All->Mesh->nElemGroup)){
      if ((elem >= 0) && (elem < All->Mesh->ElemGroup[egrp].nElem)){
	ierr = xf_Error(xf_Ref2GlobElem(All->Mesh, egrp, elem, 
					NULL, xfe_True, 1, xref, 
					xElemMarker));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_Ref2GlobElem(All->Mesh, egrp, elem, 
					NULL, xfe_True, 1, xref2, 
					xglob));
	ElemMarkerSize = 0.0;
	for (d=0; d<All->Mesh->Dim; d++)
	  ElemMarkerSize += (xElemMarker[d]-xglob[d])*(xElemMarker[d]-xglob[d]);
	ElemMarkerSize = sqrt(ElemMarkerSize);
        xf_printf("ElemMarkerSize = %.10E\n", ElemMarkerSize);
	ElemMarkerOn = xfe_True;
      }
      else xf_printf("Out of range.\n");
    }
    else xf_printf("Out of range.\n");
    
  }
  else xf_printf("Not recognized.\n");
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ChangeColorLimits
static int 
xf_ChangeColorLimits()
{  
  int ierr, ii, iDGroup;
  char buf[xf_MAXSTRLEN];
  real vmin, vmax;
  xf_DisplayVector *DV;
  xf_DisplayScalar *DS;

  if (!RenderFlag){
    xf_printf("Not currently rendering.\n");
    return xf_OK;
  }

  DV = DVector + CurrentDVector;
  DS = DV->Scalar + DV->CurrentDScalar;

  xf_printf("Set new limits (two reals, space separated) [%.4E %.4E]:\n",
	    DS->Vmin, DS->Vmax);
  
  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
      (sscanf(buf, "%lf %lf", &vmin, &vmax) == 2)){
    DS->Vmin = vmin;
    DS->Vmax = vmax;
    for (ii=0; ii<nActiveDGroup; ii++){
      iDGroup = ActiveDGroup[ii];
      ierr = xf_Error(xf_SetScalarColor(DS->Data + iDGroup, vmin, vmax));
      if (ierr != xf_OK) return ierr;
    } // ii
      
  }
  else{
    xf_printf("Not understood.\n");
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteStateFile
static int 
xf_WriteStateFile()
{  
  int ierr, k;
  xf_DisplayVector *DV;
  xf_DisplayScalar *DS;
  FILE *fid;

  if ((fid = fopen(StateFile, "w")) == NULL) return xf_Error(xf_FILE_WRITE_ERROR);

  // CutPlane
  fprintf(fid, "%% CutPlane\n");
  fprintf(fid, "%.15E %.15E %.15E %.15E\n", 
	    CutPlane[0], CutPlane[1], CutPlane[2], CutPlane[3]);

  // LineProbe on?
  fprintf(fid, "%% Line probe on?\n");
  fprintf(fid, "%d\n", LineProbeOn);

  // LineProbe
  fprintf(fid, "%% LineProbeXY [x1 y1 x2 y2]\n");
  fprintf(fid, "%.15E %.15E %.15E %.15E\n", 
	    LineProbeXY[0], LineProbeXY[1], LineProbeXY[2], LineProbeXY[3]);

  // Camera
  fprintf(fid, "%% Camera [ViewAxis Distance Azimuth Elevation Twist]\n");
  fprintf(fid, "%d %.15E %.15E %.15E %.15E\n", 
	  Camera.ViewAxis, Camera.Distance, Camera.Azimuth, Camera.Elevation, Camera.Twist);
  fprintf(fid, "%% Camera Origin [x y z]\n");
  fprintf(fid, "%.15E %.15E %.15E\n", Camera.Origin[0], Camera.Origin[1], Camera.Origin[2]);
  
  // parent pixel width + height
  fprintf(fid, "%% Parent window [width height] in pixels]\n");
  fprintf(fid, "%d %d\n", WinParent.PixWidth[0], WinParent.PixWidth[1]);

  // 2D up axis
  fprintf(fid, "%% 2D up axis [0=x, 1=y, 2=z]\n");
  fprintf(fid, "%d\n", UpAxis2D);

  // 2D/3D window center and width
  fprintf(fid, "%% 2D window center [x y]\n");
  fprintf(fid, "%.15E %.15E\n", Win2D.Center[0], Win2D.Center[1]);
  fprintf(fid, "%% 2D window width [x y]\n");
  fprintf(fid, "%.15E %.15E\n", Win2D.Width[0], Win2D.Width[1]);
  fprintf(fid, "%% 3D window center [x y]\n");
  fprintf(fid, "%.15E %.15E\n", Win3D.Center[0], Win3D.Center[1]);
  fprintf(fid, "%% 3D window width [x y]\n");
  fprintf(fid, "%.15E %.15E\n", Win3D.Width[0], Win3D.Width[1]);

  // CountBox
  fprintf(fid, "%% CountBox\n");
  fprintf(fid, "%.15E %.15E %.15E %.15E\n", 
	    CountBox[0], CountBox[1], CountBox[2], CountBox[3]);

  // Which window is big
  fprintf(fid, "%% 2D window is the big window?\n");
  fprintf(fid, "%d\n", (WinBig == &Win2D));

  // Rendering
  fprintf(fid, "%% Rendering is on?\n");
  fprintf(fid, "%d\n", RenderFlag);

  // Current display vector
  fprintf(fid, "%% Current display vector\n");
  fprintf(fid, "%d\n", CurrentDVector);

  // Current display scalar
  fprintf(fid, "%% Current display scalar\n");
  fprintf(fid, "%d\n", DVector[CurrentDVector].CurrentDScalar);

  // Scalar limits
  fprintf(fid, "%% Current scalar limits\n");
  DV = DVector + CurrentDVector;
  DS = DV->Scalar + DV->CurrentDScalar;
  fprintf(fid, "%lf %lf\n", DS->Vmin, DS->Vmax);
  
  // Which display groups are active?
  fprintf(fid, "%% Active display groups\n");
  fprintf(fid, "%d\n", nDGroup);
  for (k=0; k<nDGroup; k++)
    fprintf(fid, "%d\n", DGroup[k].Active);

  // Lighting flag
  fprintf(fid, "%% Lighting is on?\n");
  fprintf(fid, "%d\n", LightingFlag);

  fclose(fid);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadStateFile
static int 
xf_ReadStateFile()
{  
  int ierr, k, n;
  enum xfe_Bool Win2DIsBig;
  char line[xf_MAXLINELEN];
  xf_DisplayVector *DV;
  xf_DisplayScalar *DS;
  FILE *fid;

  // no biggie if not found; just return quietly with an error
  if ((fid = fopen(StateFile, "r")) == NULL) return xf_FILE_READ_ERROR;

  // CutPlane
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%lf %lf %lf %lf\n", 
	     CutPlane+0, CutPlane+1, CutPlane+2, CutPlane+3) != 4)
    return xf_Error(xf_FILE_READ_ERROR);

  // LineProbeOn
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%d\n", (int *) &LineProbeOn) != 1)
    return xf_Error(xf_FILE_READ_ERROR);

  // LineProbeXY
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%lf %lf %lf %lf\n", 
	     LineProbeXY+0, LineProbeXY+1, LineProbeXY+2, LineProbeXY+3) != 4)
    return xf_Error(xf_FILE_READ_ERROR);

  // Camera
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%d %lf %lf %lf %lf\n", 
	     &Camera.ViewAxis, &Camera.Distance, &Camera.Azimuth, 
	     &Camera.Elevation, &Camera.Twist) != 5)
    return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%lf %lf %lf\n", 
	     Camera.Origin+0, Camera.Origin+1, Camera.Origin+2) != 3)
    return xf_Error(xf_FILE_READ_ERROR);
  
  // parent pixel width + height
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%d %d\n", WinParent.PixWidth+0, WinParent.PixWidth+1) != 2)
    return xf_Error(xf_FILE_READ_ERROR);

  // 2D up axis
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%d\n", &UpAxis2D) != 1)
    return xf_Error(xf_FILE_READ_ERROR);

  // 2D/3D window center and width
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%lf %lf\n", Win2D.Center+0, Win2D.Center+1) != 2)
    return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%lf %lf\n", Win2D.Width+0, Win2D.Width+1) != 2)
    return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%lf %lf\n", Win3D.Center+0, Win3D.Center+1) != 2)
    return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%lf %lf\n", Win3D.Width+0, Win3D.Width+1) != 2)
    return xf_Error(xf_FILE_READ_ERROR);

  // CountBox
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%lf %lf %lf %lf\n", 
	     CountBox+0, CountBox+1, CountBox+2, CountBox+3) != 4)
    return xf_Error(xf_FILE_READ_ERROR);

  // 2D window is big one?
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%d\n", (int *) &Win2DIsBig) != 1)
    return xf_Error(xf_FILE_READ_ERROR);
  if ((Win2DIsBig) && (WinBig != &Win2D)){
    swap(WinBig, WinSmall, WinTemp);
    xf_SetChildrenWindowSizes();
  }

  // RenderFlag
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%d\n", (int *) &RenderFlag) != 1)
    return xf_Error(xf_FILE_READ_ERROR);

  // CurrentDVector
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%d\n", &CurrentDVector) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (CurrentDVector >= nDVector) CurrentDVector = 0; // bogus vector

  // CurrentDScalar
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%d\n", &DVector[CurrentDVector].CurrentDScalar) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  // Current scalar limits
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  DV = DVector + CurrentDVector;
  DS = DV->Scalar + DV->CurrentDScalar;
  if (sscanf(line, "%lf %lf\n", &DS->Vmin, &DS->Vmax) != 2) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  // Active display groups
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%d\n", (int *) &n) != 1)
    return xf_Error(xf_FILE_READ_ERROR);
  if (n == nDGroup){
    for (k=0; k<nDGroup; k++){
      if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
      if (sscanf(line, "%d\n", (int *) &DGroup[k].Active) != 1)
	return xf_Error(xf_FILE_READ_ERROR);
    } // k
  }

  // LightingFlag
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%d\n", (int *) &LightingFlag) != 1)
    return xf_Error(xf_FILE_READ_ERROR);

  fclose(fid);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FillLighting
static int 
xf_FillLighting()
{
  int ierr, d, i, iTri;
  int iDGroup;
  xf_DisplayTri *DTri;
  real *x0, *x1, *x2, v1[3], v2[3], *n, NN, sg;

  for (i=0; i<nActiveDGroup; i++){
    iDGroup = ActiveDGroup[i];
    DTri = DGroup[iDGroup].DTri;
    if (DTri->nNormal < DTri->nTri){
      DTri->nNormal = DTri->nTri;
      ierr = xf_Error(xf_ReAlloc((void **) &DTri->Normals, 3*DTri->nNormal, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
    for (iTri=0; iTri<DTri->nTri; iTri++){
      x0 = DTri->xglob + 3*DTri->TriList[3*iTri+0];
      x1 = DTri->xglob + 3*DTri->TriList[3*iTri+1];
      x2 = DTri->xglob + 3*DTri->TriList[3*iTri+2];
      for (d=0; d<3; d++) v1[d] = x1[d] - x0[d];
      for (d=0; d<3; d++) v2[d] = x2[d] - x0[d];
      n = DTri->Normals + 3*iTri;
      xf_CrossProduct(v1, v2, n);
      for (d=0, NN=0.; d<3; d++) NN += n[d]*n[d];
      sg = -1.0;
      if ((iDGroup == iDGroupCutPlane) && (All->Mesh->Dim == 2)) sg = 1.0;
      for (d=0, NN=sqrt(NN); d<3; d++) n[d] = sg*n[d]/NN;
    } // iTri
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MakeMovie
static int 
xf_MakeMovie(int dim)
{  
  int ierr;
  int i, iStart, iStep, iEnd;
  char buf[xf_MAXSTRLEN];
  char DataRoot[xf_MAXSTRLEN];
  char GCLRoot[xf_MAXSTRLEN];
  char DataFile[xf_MAXSTRLEN];
  char EpsFile[xf_MAXSTRLEN];
  char Title[xf_MAXSTRLEN];
  char GCLTitle[xf_MAXSTRLEN];
  char TimeHistFile[xf_MAXSTRLEN] = "NULL";
  char *s = NULL;
  enum xfe_Bool found, showp, derivescal, PlotVector, GotTimeHist, UseGCL, GCLRequired;
  xf_DisplayVector *DV;
  xf_DataSet *DataSet;
  xf_Data *D;
  xf_Vector *V, *VOrig, *Vnew = NULL;
  xf_Vector *GCL = NULL, *VGCL = NULL;
  xf_TimeHistData *TimeHistData = NULL;
  xf_GenArray *gA = NULL;
  FILE *fid;

  xf_printf("Movie animation: writes out a sequence of .eps files.\n");
  xf_printf("Each frame consists of a .data file of the form:\n");
  xf_printf("     root#.data,\n");
  xf_printf("where # is a user-specified range of integers.\n");
  xf_printf("If no data is provided, just mesh is shown.\n");
  xf_printf("\n");
  
  PlotVector = xfe_True;
  xf_printf("Enter root name for movie .data files (one string) or <enter> for none:\n");
  if ((fgets(buf, xf_MAXSTRLEN, stdin) == NULL) || 
      (sscanf(buf, "%s", DataRoot) != 1)){
    xf_printf("No .data files.\n");
    PlotVector = xfe_False;
  }
  xf_printf("Enter range for time step indices (3 ints: start step end):\n");
  if ((fgets(buf, xf_MAXSTRLEN, stdin) == NULL) || 
      (sscanf(buf, "%d %d %d", &iStart, &iStep, &iEnd) != 3)){
    xf_printf("Not understood.\n");
    return xf_INPUT_ERROR;
  }
  GotTimeHist = xfe_True;
  xf_printf("Enter name of time history file (e.g. for mesh motion) or <enter> for none:\n");
  if ((fgets(buf, xf_MAXSTRLEN, stdin) == NULL) || 
      (sscanf(buf, "%s", TimeHistFile) != 1)){
    xf_printf("No time history.\n");
    GotTimeHist = xfe_False;
  }
  else{
    ierr = xf_Error(xf_ReadTimeHistData(TimeHistFile, NULL, &TimeHistData));
    if (ierr != xf_OK) return ierr;
  }
  
  // Check to make sure .data files exist before starting movie
  if (PlotVector)
    for (i=iStart; i<=iEnd; i+=iStep){
      sprintf(DataFile, "%s%d.data", DataRoot, i);
      if ((fid = fopen(DataFile, "rb")) == NULL){
	xf_printf("File %s required for animation not found.\n", DataFile);
	return xf_INPUT_ERROR;
      }
      fclose(fid);
    } // i
  

  if (PlotVector){
    // pull off current vector and current scalar
    DV = DVector + CurrentDVector;
    VOrig = DV->Vector; // store original pointer to vector
    
    // Title
    sprintf(Title, "%s", DV->DataName);
    
    // check if showing orders for this vector
    showp = xfe_False;
    s = DV->DataName;
    if ((strlen(s) > 6) && (strncmp(s+strlen(s)-6,"_Order",6) == 0)){
      showp = xfe_True;
      Title[strlen(s)-6] = '\0';
    }
    
    // check if deriving scalars for this vector
    derivescal = xfe_False;
    s = DV->DataName;
    if ((strlen(s) > 8) && (strncmp(s+strlen(s)-8,"_Scalars",8) == 0)){
      derivescal = xfe_True;
      Title[strlen(s)-8] = '\0';
    }
  }  

  // determine if using a Geometric Conservation Law
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  if (ierr != xf_OK) return ierr;
  
  // If GCL was used, prompt user to determine whether it's necessary for animation
  if (UseGCL){
    xf_printf("\nGCL is on, but is it required to transform state for current animation? (y/n):\n");
    
    if ((fgets(buf, xf_MAXSTRLEN, stdin) == NULL) || (strncmp(buf,"y", 1)) == 0){
      GCLRequired = xfe_True;
      xf_printf("Using GCL.\n");
    }
    else{
      xf_printf("Not using GCL.\n");
      GCLRequired = xfe_False;
    }
  }
  
  // If GCL is necessary, prompt user for the title of the GCL vector to use
  if (GCLRequired){
    xf_printf("\nInput title of GCL vector to be used: \n");
    if ((fgets(buf, xf_MAXSTRLEN, stdin) == NULL) || (sscanf(buf, "%s", GCLTitle) != 1)){
    xf_printf("Not understood.\n");
    return xf_INPUT_ERROR;
    }
  }
  
  //reset UseGCL based on whether GCL is on and whether it's actually required
  UseGCL = ((UseGCL) && (GotTimeHist) && (GCLRequired)); 
  
  if (UseGCL){
    xf_printf("Loading GCL vectors from the state .data files.\n");
    ierr = xf_Error(xf_FindMeshMotionGCLVector(All, &GCL));
    if (ierr != xf_OK) UseGCL = xfe_False;
  }

  // Loop over frames
  for (i=iStart; ((iStep > 0) ? (i<=iEnd) : (i>=iEnd)); i+=iStep){
    
    if (GotTimeHist){
      if (i >= TimeHistData->nTime) return xf_Error(xf_OUT_OF_BOUNDS);
      Time = TimeHistData->Time[i];
      ierr = xf_Error(xf_MovePlotMesh(xfe_False));
      if (ierr != xf_OK) exit(ierr);
    }

    if (!PlotVector){
      sprintf(EpsFile, "plot%03d.eps", i);   // eps file i
    }
    else{

      sprintf(DataFile, "%s%d.data", DataRoot, i);  // data file i
      sprintf(EpsFile, "%s%05d.eps", DataRoot, i);  // eps file i

      // create a dataset
      ierr = xf_Error(xf_CreateDataSet(&DataSet));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ReadDataSetBinary(NULL, NULL, DataFile, DataSet));
      if (ierr!=xf_OK) return ierr;
      
      // locate appropriate vector
      D = DataSet->Head;
      found = xfe_False;
      while (D != NULL){
	if (D->Type == xfe_Vector){
	  V = (xf_Vector *) D->Data;
	  if (strcmp(D->Title, Title) == 0){
	    DV->Vector = V;
	    found = xfe_True;
	    break;
	  }
	}
	D = D->Next;
      }
      if (!found) {
        xf_printf("Not found: Title = %s\n", Title);
        return xf_Error(xf_NOT_FOUND); // data of same name not found
      }
      
      Vnew = NULL;
      
      // build order if current vector has _Order in DataName
      if (showp){
	ierr = xf_Error(xf_BuildVOrder(All, DV->Vector, &Vnew));
	if (ierr != xf_OK) return ierr;
	
	// Check/prepare integer or non-interpolated vector
	ierr = xf_Error(xf_InterpOrIntPrep(Vnew));
	if (ierr != xf_OK) return ierr;
	// set for plotting
	DV->Vector = Vnew;
      }
      // point to new scalars if derivescal is True
      if (derivescal){
	// Nothing to do since DV->Vector is already pointing to V!
      }
      
      // find GCL vector in dataset if GCL is on
      if (UseGCL){
	ierr = xf_FindDataByTitle(DataSet, GCLTitle, xfe_Vector, &D);
	if (ierr == xf_NOT_FOUND)
	  xf_printf("GCL vector not found in %s -- not using GCL vector.\n", DataFile);
	else if (ierr == xf_OK){  // GCL found
	  VGCL = (xf_Vector *) D->Data;

	  // use the just-loaded GCL vector
	  swap(GCL->GenArray, VGCL->GenArray, gA);
	}	
      }

      // Check/prepare integer or non-interpolated vector
      ierr = xf_Error(xf_InterpOrIntPrep(DV->Vector));
      if (ierr != xf_OK) return ierr;
      
      // recompute isosurface if present
      if ((iDGroupIsoSurf >= 0) && (DGroup[iDGroupIsoSurf].Active)){
	ierr = xf_Error(xf_CalculateDGroup(DGroup + iDGroupIsoSurf));
	if (ierr != xf_OK) return ierr;
      }

      // draw hills if flag set
      if ((DrawingHills) && (iDGroupCutPlane >= 0) && 
	  (DGroup[iDGroupCutPlane].Active)){
	ierr = xf_Error(xf_CalculateDGroup(DGroup + iDGroupCutPlane));
	if (ierr != xf_OK) exit(ierr);
      }
      
      // fill scalars (keep state limits)
      ierr = xf_Error(xf_FillDScalarsAll(DV, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    
    // re-calculate lighting if on
    if (LightingFlag){
      ierr = xf_Error(xf_FillLighting());
      if (ierr != xf_OK) exit(ierr);
    }
    
    // display contours if on
    if (ContoursOn){
      ierr = xf_Error(xf_FillContours());
      if (ierr != xf_OK) exit(ierr);
    }

    
    // re-display window
    if (dim == 2)
      xf_Display2D();
    else
      xf_Display3D();

    // capture picture
    if (dim == 2){
      ierr = xf_Error(xf_HardcopyWindowEPS(Win2D, EpsFile, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    else{
      ierr = xf_Error(xf_HardcopyWindowEPS(Win3D, EpsFile, xfe_True));
      if (ierr != xf_OK) return ierr;
    }

    if (PlotVector){

      // swap out GCL data before destroying dataset if using GCL
      if (UseGCL) swap(GCL->GenArray, VGCL->GenArray, gA);

      // destroy dataset
      ierr = xf_Error(xf_DestroyDataSet(DataSet));
      if (ierr != xf_OK) return ierr;
      
      // destroy Vnew
      if (Vnew != NULL){
	ierr = xf_Error(xf_DestroyVector(Vnew, xfe_True));
	if (ierr != xf_OK) return ierr;
      }
    }
    
  } // i

  /* Destroy Time history */
  ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
  if (ierr != xf_OK) return ierr;

  // reset vector pointer
  if (PlotVector)
    DV->Vector = VOrig;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_PostRedisplayAll
static void
xf_PostRedisplayAll()
{
  int orig;

  orig = glutGetWindow();

  glutSetWindow(Win2D.Handle);
  glutPostRedisplay();
  glutSetWindow(Win3D.Handle);
  glutPostRedisplay();
  glutSetWindow(WinAux.Handle);
  glutPostRedisplay();
  glutSetWindow(WinBorder.Handle);
  glutPostRedisplay();

  if (orig > 0) glutSetWindow(orig);

}

/******************************************************************/
//   FUNCTION Definition: xf_TurnRenderOn
static void 
xf_TurnRenderOn() {
  
  int ierr, i, orig;
  int CurrentDScalar;
  enum xfe_Bool Redisplay = xfe_False;
  xf_DisplayGroup *DG;

  // initialize to False in case we cannot turn rendering on
  RenderFlag = xfe_False;

  // make sure we have a vector with data
  if (nDVector <= 0){
    xf_printf("No Vectors loaded for rendering.  Use 'v' to load a vector.\n");
    fflush(stdout);
  }
  else if ((CurrentDVector < 0) || (CurrentDVector >= nDVector)){
    xf_printf("Invalid CurrentDVector = %d.\n", CurrentDVector);
    fflush(stdout);
  }
  else if (DVector[CurrentDVector].nScalar <= 0){
    xf_printf("No scalars present in current vector.\n");
    fflush(stdout);
  }
  else if ((DVector[CurrentDVector].CurrentDScalar < 0) ||
	   (DVector[CurrentDVector].CurrentDScalar >= DVector[CurrentDVector].nScalar)){
    xf_printf("Invalid CurrentDScalar = %d.\n", DVector[CurrentDVector].CurrentDScalar);
    fflush(stdout);
  }
  else{
    RenderFlag = xfe_True;
    CurrentDScalar = DVector[CurrentDVector].CurrentDScalar;
    Redisplay = xfe_False;
    for (i=0; i<nActiveDGroup; i++){
      DG = DGroup + ActiveDGroup[i];
      if (DVector[CurrentDVector].Scalar[CurrentDScalar].Data[i].Version != 
	  DG->DTri->Version){
	Redisplay = xfe_True;
	break;
      }
    } // i
    if (Redisplay){
      ierr = xf_Error(xf_FillDScalarsAll(DVector + CurrentDVector, xfe_True));
      if (ierr != xf_OK) exit(ierr);
    }
    xf_WriteScalarInfo();
  }
  
  orig = glutGetWindow();
  glutSetWindow(WinAux.Handle);
  glutPostRedisplay();
  if (orig > 0) glutSetWindow(orig);

}
    

/******************************************************************/
//   FUNCTION Definition: xf_DisplayKeyMenuCommon
static void 
xf_DisplayKeyMenuCommon () {
  xf_printf("\n");
  xf_printf(" +----------------- Common Key Menu ------------------------+\n");
  xf_printf(" | v : load a new vector        a : disp-refine all elem    |\n");
  xf_printf(" | s : change scalar (also 0-9) A : disp-coarsen all elem   |\n");
  xf_printf(" | r : toggle render scalar     p : show interp orders      |\n");
  xf_printf(" | m : toggle mesh              P : <unused> |\n");
  xf_printf(" | M : toggle display mesh      e : disp-refine elem @ x,y  |\n");
  xf_printf(" | c : toggle contours          E : disp-coarsen elem @ x,y |\n");
  xf_printf(" | C : set # of contours        g : toggle log color axis   |\n");
  xf_printf(" | b : toggle BFGs (3D)         G : [unused]                |\n");
  xf_printf(" | f : change color limits      X : toggle axes             |\n");
  xf_printf(" | w : write xflow state        o : render output points    |\n");
  xf_printf(" | i : identify an element      d : derive a vector         |\n");
  xf_printf(" | F : toggle mesh transparent  n : check for neg volumes   |\n");
  xf_printf(" | B : black-white rendering    D : compare vectors L2      |\n");
  xf_printf(" | W : read xflow state         t : set mesh motion time    |\n");
  xf_printf(" |                                                          |\n");
  xf_printf(" | Esc, q: Exit                                             |\n");
  xf_printf(" | Tab : swap 2D and 3D windows                             |\n");
  xf_printf(" |                                                          |\n");
  xf_printf(" | Mouse:                                                   |\n");
  xf_printf(" |  Middle button down + move = scroll                      |\n");
  xf_printf(" |  Shift + Middle button down + move up/down = zoom        |\n");
  xf_printf(" +----------------------------------------------------------+\n");
  xf_printf("\n");
}

/******************************************************************/
//   FUNCTION Definition: xf_DisplayKeyMenu2D
static void 
xf_DisplayKeyMenu2D () {
  xf_printf("\n");
  xf_printf(" +------------------- 2D Window Menu -----------------------+\n");
  xf_printf(" | h : hardcopy plot to .eps    H : hardcopy movie          |\n");
  xf_printf(" | F5: acquire line probe       F6: write line.txt          |\n");
  xf_printf(" | F7: toggle line probe on/off F8: count elems in box      |\n");
  xf_printf(" | b : toggle disp boundaries   l : clip cut plane          |\n");
  xf_printf(" | x,y,z: set up axis in the x, y, or z direction           |\n");
  xf_printf(" | Right button down = print elem info                      |\n");
  xf_printf(" +----------------------------------------------------------+\n");
  xf_printf("\n");
}


/******************************************************************/
//   FUNCTION Definition: xf_DisplayKeyMenu3D
static void 
xf_DisplayKeyMenu3D () {
  xf_printf("\n");
  xf_printf(" +------------------- 3D Window Menu -----------------------+\n");
  xf_printf(" | h : hardcopy plot to .eps    H : hardcopy movie          |\n");
  xf_printf(" | F3: position cut plane       F4: set iso-surface         |\n");
  xf_printf(" | l : toggle lighting                                      |\n");
  xf_printf(" | x,y,z: set view from -x, -y, or -z axis                  |\n");
  xf_printf(" | Ctrl + Middle button down + move = rotate                |\n");
  xf_printf(" +----------------------------------------------------------+\n");
  xf_printf("\n");
}


/*--------------------------*/
/* function declaration for */
/* tecplot data dump        */
/*--------------------------*/
static int
write_string_to_tecplot(FILE *fout, const char *fname, int len)
{
   int i, tmp;
   int *StrIntData;

   for(i=0; i<len; i++){
      tmp = (int)fname[i];
      fwrite(&tmp, sizeof(int), 1, fout);
   }

   //the last digit should be 0
   tmp = 0;
   fwrite(&tmp, sizeof(int), 1, fout);

   return xf_OK;
}
static int
Yu_DumpMeshOnScreenToTecplot(xf_DisplayVector *DV, int iDGroup, FILE *fmsh)
{

   int ierr, d, i, j, k, dim, sr;
   int Nnodes, Nelements;
   real *x;
   xf_DisplayGroup *DG;
   xf_DisplayLoop *DLoop;  //for mesh
   xf_DisplayTri *DTri;    //for element
   
   DG = DGroup + iDGroup;
   DTri = DG->DTri;
   DLoop = DG->DLoop;
   
   //data
   //first find how many nodes and linseg elements
   Nnodes = 0; Nelements = 0;
   for(i=0; i<DLoop->nLoop; i++){
	   if ((DLoop->Active[i] == 1) || (DLoop->Active[i] == -1))
		   continue;
	   
	   if (DLoop->Active[i] > 0)   //there is line loop (i.e. triangle mesh)
		for(k=0; k<DLoop->nNode[i]; k++)
{	
			Nnodes++;
			Nelements++;
}
           else
{
		   for (k=0; k<DLoop->nNode[i]; k++) //there is no line loop (i.e. quad mesh)
		   {
			   Nnodes++;
			   Nelements++;
		   }
		   Nelements -= 2;   
} 
  }
   
   //zone header (in ASCII format)
   fprintf(fmsh, "ZONE DATAPACKING=POINT, NODES=%d, ELEMENTS=%d, ZONETYPE=FELINESEG\n",
           Nnodes, Nelements);
   
		   //write nodes coordinates
    for(i=0; i<DLoop->nLoop; i++){
	    if ((DLoop->Active[i] == 1) || (DLoop->Active[i] == -1))
		continue;
	   
		  //if (DLoop->Active[i] > 0)
		 for (k=0; k<DLoop->nNode[i]; k++)
		 {
			 x = DTri->xglob + 3*DLoop->NodeList[i][k];
			 fprintf(fmsh, "%.12lf %.12lf %.12lf\n", x[0], x[1], x[2]);
		  }  
    }
   
   fprintf(fmsh, "\n");
   
    //write connectivity
    Nnodes = 1;
    //write nodes coordinates
    for(i=0; i<DLoop->nLoop; i++){
	    if ((DLoop->Active[i] == 1) || (DLoop->Active[i] == -1))
            continue;
        
        if (DLoop->Active[i] > 0)  //drawing line loop
{       
         for (k=0; k<DLoop->nNode[i]-1; k++)
        {
            Nnodes++;
            fprintf(fmsh, "%d %d\n", Nnodes-1, Nnodes);
        }
            fprintf(fmsh, "%d %d\n", Nnodes, Nnodes-DLoop->nNode[i]+1);
            Nnodes++;
}
else  //no drawing line loop
{
if(DLoop->nNode[i]%2 == 1)
printf("warning, odd number of lines for GL_Lines\n");

	for (k=0; k<DLoop->nNode[i]/2; k++)
         { 
             Nnodes += 2;
             fprintf(fmsh, "%d %d\n", Nnodes-1, Nnodes-2);
             
          }
}      
    }
		
   return xf_OK;
}
static int
Yu_DumpHeaderToTecplot(xf_DisplayVector *DV, int CurrentDScalar, int iDGroup, FILE *fout)
{
    int ierr, d, i, j, k, dim, sr;
    int egrp, elem, iDNode, nDNode, pnDNode;
    int ielem, Order;
    real *x, *Scalar, *S, xmin[3], xmax[3], Smin, Smax;
    float fCvter;
    double dCvter;
    char buf[xf_MAXSTRLEN];
    enum xfe_BasisType Basis;
    enum xfe_Bool Interpolated;
    xf_DisplayGroup *DG;
    xf_DisplayTri *DTri;    //for element
    xf_DisplayLoop *DLoop;  //for mesh
    xf_DisplayStrip *DStrip;//for boundary outline
    
    DG = DGroup + iDGroup;
    DTri = DG->DTri;
    DLoop = DG->DLoop;
    DStrip = DG->DStrip;
   
    Scalar = DV->Scalar[CurrentDScalar].Data[iDGroup].Value;
   
    //first write zone information in ASCII format
    //fprintf(fout, "ZONE DATAPACKING=POINT, NODES=%d, ELEMENTS=%d, ZONETYPE=FETRIANGLE\n",
    //        3*DTri->nTri, DTri->nTri);

    //zone info in binary format
    //Zone marker 
    fCvter = 299.0;
    fwrite(&fCvter, sizeof(float), 1, fout);
    //Zone number
    sprintf(buf, "ZONE_%d", iDGroup+1);
    write_string_to_tecplot(fout, buf, strlen(buf));
    //ParentZone
    k=-1; fwrite(&k, sizeof(int), 1, fout);
    //StrandID
    k=-2; fwrite(&k, sizeof(int), 1, fout);
    //solution time (not used now)
    dCvter = 0.0;
    fwrite(&dCvter, sizeof(double), 1, fout);
    //-1 value; not useful
    k=-1; fwrite(&k, sizeof(int), 1, fout);
    //ZoneType: 2=FETRIANGLE
    k=2; fwrite(&k, sizeof(int), 1, fout);
    //Datapacking: 0=block
    k=0; fwrite(&k, sizeof(int), 1, fout);
    //Specify Var location
    k=0; fwrite(&k, sizeof(int), 1, fout);
    //raw local 1-to-1
    k=0; fwrite(&k, sizeof(int), 1, fout);
    //misc user-defined face
    //k=0; fwrite(&k, sizeof(int), 1, fout);
    //if FE zone: NumPts
    k=3*DTri->nTri; fwrite(&k, sizeof(int), 1, fout);
    //if FE zone: Numelements
    k=DTri->nTri; fwrite(&k, sizeof(int), 1, fout);
    //CellDim info for future use
    k=0; fwrite(&k, sizeof(int), 1, fout);
    k=0; fwrite(&k, sizeof(int), 1, fout);
    k=0; fwrite(&k, sizeof(int), 1, fout);
    //Auxiliary data name/value pair
    k=0; fwrite(&k, sizeof(int), 1, fout);
    
    return xf_OK;
}
static int
Yu_DumpDataToTecplot(xf_DisplayVector *DV, int CurrentDScalar, int iDGroup, FILE *fout)
{
   int ierr, d, i, j, k, dim, sr;
   int egrp, elem, iDNode, nDNode, pnDNode;
   int ielem, Order;
   real *x, *Scalar, *S, xmin[3], xmax[3], Smin, Smax;
   float fCvter;
   double dCvter;
   char buf[xf_MAXSTRLEN];
   enum xfe_BasisType Basis;
   enum xfe_Bool Interpolated;
   xf_DisplayGroup *DG;
   xf_DisplayTri *DTri;    //for element
   xf_DisplayLoop *DLoop;  //for mesh
   xf_DisplayStrip *DStrip;//for boundary outline
   
   DG = DGroup + iDGroup;
   DTri = DG->DTri;
   DLoop = DG->DLoop;
   DStrip = DG->DStrip;
   
   Scalar = DV->Scalar[CurrentDScalar].Data[iDGroup].Value;
   
   //Header and Data Section separation
   //fCvter=357.0; fwrite(&fCvter, sizeof(float), 1, fout);
   //Zone magic number
   fCvter=299.0; fwrite(&fCvter, sizeof(float), 1, fout);
   //variable precision
   for(i=0; i<4; i++)
   {k=2; fwrite(&k, sizeof(int), 1, fout);}
   //Has passive variables
   k=0; fwrite(&k, sizeof(int), 1, fout);
   //Is variable passive
   k=0; fwrite(&k, sizeof(int), 1, fout);
   //share connectivity
   k=-1; fwrite(&k, sizeof(int), 1, fout);
   
   //first dump min & max for each variable
   //find max-min for each variables
   for(j=0; j<3; j++) {xmin[j] = 1.0e+16; xmax[j] = -1.0e+16;}
   Smin = 1.0e+16; Smax = -1.0e+16;

   for(k=0; k<DTri->nTri; k++){
      
      if(DTri->Active[k] == 0) continue;
      
      for(i=0; i<3; i++){
         iDNode = DTri->TriList[3*k+i];
         x = DTri->xglob + 3*iDNode;
         S = Scalar + iDNode;
         
         //for coordinate
         for(j=0; j<3; j++){
            if(x[j] < xmin[j])
               xmin[j] = x[j];
            
            if(x[j] > xmax[j])
               xmax[j] = x[j];
         }//j
         
         //for scalar
         if(S[0] < Smin)
            Smin = S[0];
         if(S[0] > Smax)
            Smax = S[0];
         
      }
      
   }
   
   //then dump max & min for variables
   for(j=0; j<3; j++){
      dCvter = xmin[j]; fwrite(&dCvter, sizeof(double), 1, fout);
      dCvter = xmax[j]; fwrite(&dCvter, sizeof(double), 1, fout);
   }
   
   dCvter = Smin; fwrite(&dCvter, sizeof(double), 1, fout);
   dCvter = Smax; fwrite(&dCvter, sizeof(double), 1, fout);
   
   //Dump coordinates and variable magnitude
   for(j=0; j<3; j++ ) {
      for(k=0; k<DTri->nTri; k++){
         
         if(DTri->Active[k] == 0) continue;
         
         for(i=0; i<3; i++){
            iDNode = DTri->TriList[3*k+i];
            x = DTri->xglob + 3*iDNode;
            S = Scalar + iDNode;
            
            //for ASCII format
            //coordinates
            //fprintf(fout, "%.12lf %.12lf %.12lf", x[0], x[1], x[2]);
            //values
            //fprintf(fout, "  %.12lf\n", S[0]);
            
            dCvter = x[j]; fwrite(&dCvter, sizeof(double), 1, fout);
            //dCvter = x[1]; fwrite(&dCvter, sizeof(double), 1, fout);
            //dCvter = x[2]; fwrite(&dCvter, sizeof(double), 1, fout);
            //dCvter = S[0]; fwrite(&dCvter, sizeof(double), 1, fout);
         }
      }
   }
   for(k=0; k<DTri->nTri; k++){
      
      if(DTri->Active[k] == 0) continue;
      
      for(i=0; i<3; i++){
         iDNode = DTri->TriList[3*k+i];
         x = DTri->xglob + 3*iDNode;
         S = Scalar + iDNode;
         
         dCvter = S[0]; fwrite(&dCvter, sizeof(double), 1, fout);
         
      }
   }
   
   //for ASCII format
   //fprintf(fout, "\n");
   
   //Dump connectivity (non-conforming)
   //actually arbitrary
   for(k=0; k<DTri->nTri; k++){
      
      if(DTri->Active[k] == 0) continue;
      
      //for ASCII format
      //fprintf(fout, "%d %d %d\n", 3*k+1, 3*k+2, 3*k+3);
      
      i = 3*k+0; fwrite(&i, sizeof(int), 1, fout);
      i = 3*k+1; fwrite(&i, sizeof(int), 1, fout);
      i = 3*k+2; fwrite(&i, sizeof(int), 1, fout);
   }
   
   return xf_OK;
}

/*
char *trimwhitespace(char *str)
{
     char *end;

     //Trim leading space
     while(isspace(*str)) str++;
       
     if(*str == 0)  // All spaces?
     return str;
       
     // Trim trailing space
     end = str + strlen(str) - 1;
     while(end > str && isspace(*end)) end--;
       
     // Write new null terminator
     *(end+1) = 0;
       
     return str;
}
*/ 
static int
Yu_DumpAllDataToTecplot(xf_DisplayVector *DV)
{
    int ierr, i, ii, k;
    int iDGroup;
    int CurrentDScalar;
    char buf[xf_MAXSTRLEN];
    char fname[xf_MAXSTRLEN];
    FILE *fout, *fmsh;
    float fCvter;
    enum xfe_Bool GainFileName;
   
    //ask if user want to visualize mesh file
    xf_printf("Do you want to save current mesh on screen:\n");
    fgets(buf, 80, stdin);

    while(strcmp(buf, "Y\n") != 0 && (strcmp(buf, "N\n") != 0))
    {
       xf_printf("Input is not understandable; Y-yes; N-no!\n");
       xf_printf("Do you want to save current mesh on screen:\n");
       scanf("%s", buf);
    }

    //need output the mesh on screen
    if(strcmp(buf, "Y\n") == 0)
    {
       //output the mesh file (in ASCII format)
		if(MeshOn){
	   
	   fmsh = fopen("ScreenMesh.dat", "w");
       sprintf(buf, "TITLE = MeshOnScreen\n");
       xf_fwrite(buf,sizeof(char),strlen(buf),fmsh);
       sprintf(buf,"VARIABLES = \"X\", \"Y\", \"Z\" \n");
       xf_fwrite(buf,sizeof(char),strlen(buf),fmsh);
	   
	   for(i=0; i<nActiveDGroup; i++){
		   ierr = xf_Error(Yu_DumpMeshOnScreenToTecplot(DV, ActiveDGroup[i], fmsh));
           if(ierr != xf_OK) return ierr;
	   }
	   
	     fclose(fmsh);
       }
	   else
	   {
		   xf_printf("Please put Mesh on screen and try again\n");
	   }
    }
    else if(strcmp(buf, "N\n") == 0)  //no need to output the mesh
    {
       //do nothing
    }
    else
    {
       xf_printf("Logistic error; check xf_Plot(line 5785)\n");
       return xf_OK;
    }
    
    //wait user to give file name
    GainFileName = xfe_True;
    xf_printf("Enter the name of tecplot .dat file:\n");
    if ((fgets(buf, xf_MAXSTRLEN, stdin) == NULL||
         (sscanf(buf, "%s", fname) != 1))){
        xf_printf("No .dat file.\n");
        GainFileName = xfe_False;
    }
    
    if(GainFileName == xfe_True){
        
        fout = fopen(fname,"w");
        //find variable name
        CurrentDScalar = DV->CurrentDScalar;

        //write header section for ASCII format
        /*
        sprintf(buf, "TITLE = %s\n", fname);
        ierr = xf_fwrite(buf,sizeof(char),strlen(buf),fout);
        sprintf(buf,"VARIABLES = \"X\", \"Y\", \"Z\", \"%s\" \n", DV->Scalar[CurrentDScalar].Name);
        ierr = xf_fwrite(buf,sizeof(char),strlen(buf),fout);
        */
     

        //write header section for binary format
        //Magic number & Version number
        sprintf(buf, "#!TDV112");
        fwrite(buf, sizeof(char), 8, fout);
        //Integer value of 1
        k=1; fwrite(&k, sizeof(int), 1, fout);
        //FileType: 0=Full; 1=Grid; 2=Solution
        k=0; fwrite(&k, sizeof(int), 1, fout);
        //TITLE
        write_string_to_tecplot(fout, fname, strlen(fname));
        //number of variables
        k=4; fwrite(&k, sizeof(int), 1, fout);
        //dump variable name
           //X
           sprintf(buf,"X");
           write_string_to_tecplot(fout, buf, strlen(buf));
           //Y
           sprintf(buf,"Y");
           write_string_to_tecplot(fout, buf, strlen(buf));
           //Z
           sprintf(buf,"Z");
           write_string_to_tecplot(fout, buf, strlen(buf));
           //variable name
           sprintf(buf, "%s", DV->Scalar[CurrentDScalar].Name);
           write_string_to_tecplot(fout, buf, strlen(buf)); 
        //after this, see the next function call.


       //output zones info
       for(i=0; i<nActiveDGroup; i++){
          ierr = xf_Error(Yu_DumpHeaderToTecplot(DV, CurrentDScalar, ActiveDGroup[i], fout));
          if(ierr != xf_OK) return ierr;
       }
       
       fCvter=357.0; fwrite(&fCvter, sizeof(float), 1, fout);
       
        xf_printf("here 1; %s\n", DV->Scalar[CurrentDScalar].Name);
       //output data
        for(i=0; i<nActiveDGroup; i++){
            ierr = xf_Error(Yu_DumpDataToTecplot(DV, CurrentDScalar, ActiveDGroup[i], fout));
            if(ierr != xf_OK) return ierr;
        }
       
        fclose(fout);
    }
    else
    {
        xf_printf("Please try data-dump functionality again!\n");
    }
    
    return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_KeyCommon
void 
xf_KeyCommon (unsigned char k, int x, int y) {
  int ierr, i, j, iScalar, orig, dim;
  int egrp, elem, Order;
  int *OrderVec = NULL;
  int CurrentDScalar;
  char buf[xf_MAXSTRLEN];
  enum xfe_Bool Redisplay = xfe_False;
  enum xfe_Bool AskForTime = xfe_False;
  xf_DisplayGroup *DG;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  dim = Mesh->Dim;

  if (k==27 || k=='q') {
    xf_DestroyGlobals();
    xf_printf("xf_Plot finished.\n"); fflush(stdout);
    exit(0);
  }

  if (k=='r') {
    RenderFlag = !RenderFlag;
    if (RenderFlag) xf_TurnRenderOn();
  }

  if (k=='F') {
    MeshFillFlag = !MeshFillFlag;
    orig = glutGetWindow();
    glutSetWindow(WinAux.Handle);
    glutPostRedisplay();
    if (orig > 0) glutSetWindow(orig);
  }
  
  if ((k >= '0') && (k <= '9')){
    iScalar = k-'0';
    if (nDVector <= 0){
      xf_printf("No Vectors loaded for rendering.  Use 'v' to load a vector.\n");
      fflush(stdout);
    }
    else if ((CurrentDVector < 0) || (CurrentDVector >= nDVector)){
      xf_printf("Invalid CurrentDVector = %d.\n", CurrentDVector);
      fflush(stdout);
    }
    else if ( (iScalar = k-'0') >= DVector[CurrentDVector].nScalar){
      xf_printf("Only %d scalars available (index starts at 0).\n",
		DVector[CurrentDVector].nScalar); 
      fflush(stdout);
    }
    else
      DVector[CurrentDVector].CurrentDScalar = iScalar;
    if (ContoursOn){
      ierr = xf_Error(xf_FillContours());
      if (ierr != xf_OK) exit(ierr);
    }
    xf_WriteScalarInfo();
    orig = glutGetWindow();
    glutSetWindow(WinAux.Handle);
    glutPostRedisplay();
    if (orig > 0) glutSetWindow(orig);
  }

  if (k == 's'){
    xf_ChooseScalar();
    xf_WriteScalarInfo();
    orig = glutGetWindow();
    glutSetWindow(WinAux.Handle);
    glutPostRedisplay();
    if (orig > 0) glutSetWindow(orig);
    xf_ModifySubMenu_ChangeScalar();
  }

  if (k == 'v'){
    xf_ChooseVector();
    orig = glutGetWindow();
    glutSetWindow(WinAux.Handle);
    glutPostRedisplay();
    if (orig > 0) glutSetWindow(orig);
    if (RenderFlag) xf_TurnRenderOn();
    ContoursOn = xfe_False;
    xf_ModifySubMenu_LoadVector();
    xf_ModifySubMenu_ChangeScalar();
  }

  if (k == 'd'){
    ierr = xf_Error(xf_DeriveVector());
    if (ierr != xf_OK) exit(ierr);
    orig = glutGetWindow();
    glutSetWindow(WinAux.Handle);
    glutPostRedisplay();
    if (orig > 0) glutSetWindow(orig);
    RenderFlag = xfe_False;
    ContoursOn = xfe_False;
    xf_ModifySubMenu_LoadVector();
    xf_ModifySubMenu_ChangeScalar();
  }

  if (k == 'p'){
    ierr = xf_Error(xf_ShowInterpOrderVector());
    if (ierr != xf_OK) exit(ierr);
    orig = glutGetWindow();
    glutSetWindow(WinAux.Handle);
    glutPostRedisplay();
    if (orig > 0) glutSetWindow(orig);
    RenderFlag = xfe_False;
    ContoursOn = xfe_False;
  }


  if (k == 'D'){
    ierr = xf_Error(xf_CompareVectors());
    if (ierr != xf_OK) exit(ierr);
  }

  if (k == 'm'){
    MeshOn = !MeshOn;
    xf_ModifySubMenu_Mesh();
  }
  if (k == 'M'){
    SubMeshOn = !SubMeshOn;
    xf_ModifySubMenu_Mesh();
  }

  if (k == 'a'){
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
	Elem2Refine[egrp][elem] += 1;
    ierr = xf_Error(xf_CalculateActiveDGroups());
    if (ierr != xf_OK) exit(ierr);
    if (MeshMotionInPhysical){
      ierr = xf_Error(xf_MovePlotMesh(xfe_False));
      if (ierr != xf_OK) exit(ierr);
    }
    if (RenderFlag){
      ierr = xf_Error(xf_FillDScalarsAll(DVector + CurrentDVector, xfe_True));
      if (ierr != xf_OK) exit(ierr);
    }
    if (LightingFlag){
      ierr = xf_Error(xf_FillLighting());
      if (ierr != xf_OK) exit(ierr);
    }
    if (ContoursOn){
      ierr = xf_Error(xf_FillContours());
      if (ierr != xf_OK) exit(ierr);
    }
  }

  if (k == 'A'){
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
	Elem2Refine[egrp][elem] = max(1, Elem2Refine[egrp][elem]-1);
    ierr = xf_Error(xf_CalculateActiveDGroups());
    if (ierr != xf_OK) exit(ierr);
    if (RenderFlag){
      ierr = xf_Error(xf_FillDScalarsAll(DVector + CurrentDVector, xfe_True));
      if (ierr != xf_OK) exit(ierr);
    }
    if (LightingFlag){
      ierr = xf_Error(xf_FillLighting());
      if (ierr != xf_OK) exit(ierr);
    }
    if (ContoursOn){
      ierr = xf_Error(xf_FillContours());
      if (ierr != xf_OK) exit(ierr);
    }
  }

  if (k == '?') xf_DisplayKeyMenuCommon();

  if (k == 9){ // tab
    // swap 2d and 3d windows
    swap(WinBig, WinSmall, WinTemp);
    xf_Resize(WinParent.PixWidth[0], WinParent.PixWidth[1]);
  }

  if (k == 'b'){
    xf_ToggleBFG();

    for (i=0; i<nActiveDGroup; i++){
      if (DGroup[ActiveDGroup[i]].Type != xfe_DGCutPlane){
	ierr = xf_Error(xf_CalculateDGroup(DGroup + ActiveDGroup[i]));
	if (ierr != xf_OK) exit(ierr);
      }
    }
    //ierr = xf_Error(xf_CalculateActiveDGroups());
    //if (ierr != xf_OK) exit(ierr);
    if (RenderFlag){
      ierr = xf_Error(xf_FillDScalarsAll(DVector + CurrentDVector, xfe_True));
      if (ierr != xf_OK) exit(ierr);
    }
    if (dim == 2) BoundaryOutline = !BoundaryOutline;
    xf_ModifySubMenu_Boundaries();
  }

  if (k == 'f'){
    ierr = xf_Error(xf_ChangeColorLimits());
    if (ierr != xf_OK) exit(ierr);
    if (ContoursOn){
      ierr = xf_Error(xf_FillContours());
      if (ierr != xf_OK) exit(ierr);
    }
    xf_PostRedisplayAll();
  }

  if (k == 'g'){
    UseLogScale = !UseLogScale;
    if (UseLogScale)
      xf_printf("Now using log color scale.\n");
    else
      xf_printf("Now using linear color scale.\n");
    // fill scalars (keep state limits)
    ierr = xf_Error(xf_FillDScalarsAll(DVector + CurrentDVector, xfe_True));
    if (ierr != xf_OK) exit(ierr);
    orig = glutGetWindow();
    glutSetWindow(WinAux.Handle);
    glutPostRedisplay();
    if (orig > 0) glutSetWindow(orig);
  }

  if (k == 'X'){
    AxesFlag = !AxesFlag;
    if (AxesFlag) xf_printf("red = x, green = y, blue = z\n");
  }

  if (k == 'w'){
    ierr = xf_Error(xf_WriteStateFile());
    if (ierr != xf_OK)
      xf_printf("Error writing state to file.\n");
    else
      xf_printf("Wrote plotter state to file %s\n", StateFile);
  }

  if (k == 'W'){
    ierr = xf_Error(xf_ReadStateFile());
    if (ierr != xf_OK)
      xf_printf("Error reading state from file.\n");
    else
      xf_printf("Read plotter state file %s\n", StateFile);
    xf_Resize(WinParent.PixWidth[0], WinParent.PixWidth[1]);
    xf_PostRedisplayAll();
  }

  if (k == 'o'){
    ierr = xf_Error(xf_ToggleOutputs());
    if (ierr != xf_OK) xf_printf("Error %d while Toggling Outputs.\n", ierr);
  }

  if (k == 'i'){
    ierr = xf_Error(xf_IdentifyElement());
    if (ierr != xf_OK) xf_printf("Error %d while Identifying Element.\n", ierr);
  }


  if (k == 'c'){
    if (ContoursOn)
      ContoursOn = xfe_False;
    else if (RenderFlag){
      ierr = xf_Error(xf_FillContours());
      if (ierr == xf_OK){
	ContoursOn = xfe_True;
	xf_printf("Contours are on.\n");
      }
      else
	xf_printf("Error %d while filling contours.\n", ierr);
    }
    xf_ModifySubMenu_Contours();
  }

  if (k == 'C'){
    xf_printf("nContours = %d.  Enter new number of contours:\n", nContours);
    if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	(sscanf(buf, "%d", &i) == 1))
      nContours = min(max(i,2), 200);
    else
      xf_printf("Not understood.\n");
    if (ContoursOn){
      ierr = xf_Error(xf_FillContours());
      if (ierr != xf_OK) exit(ierr);
    }
    xf_ModifySubMenu_Contours();
  }

  if (k == 'n'){

    NegVolOn = !NegVolOn;
    
    if (NegVolOn){
      xf_printf("Enter interpolation order for negative volume check:\n");
      if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	  (sscanf(buf, "%d", &i) == 1)){
	Order = max(1,i);
	
	ierr = xf_Error(xf_Alloc( (void **) &OrderVec, Mesh->nElemGroup, sizeof(int)));
	if (ierr != xf_OK) exit(ierr);
	
	for (j=0; j<Mesh->nElemGroup; j++) OrderVec[j] = Order;
	
	xf_printf("Checking volumes ...\n");
	ierr = xf_Error(xf_CheckVolumes(Mesh, OrderVec, xfe_False, xfp_MAX_NEG_POINT, 
					&nNegVolPoints, xNegVolPoints));
	if (ierr != xf_OK) exit(ierr);
	
	xf_printf("nNegVolPoints = %d\n", nNegVolPoints);

	xf_Release( (void *) OrderVec);
      }
      else
	xf_printf("Not understood.\n");
    }

  }

  if (k == 'B'){
    RenderMap = (RenderMap+1)%xfe_RenderMapLast;
    xf_printf("Render map = %s\n", xfe_RenderMapName[RenderMap]);
    if (RenderFlag){
      ierr = xf_Error(xf_FillDScalarsAll(DVector + CurrentDVector, xfe_True));
      if (ierr != xf_OK) exit(ierr);
    }
  }

  if (k == 't'){
    if (MeshMotionInPhysical){
      AskForTime = xfe_False;
      Time = 0.;
      xf_printf("Showing mesh at Time=0\n");
    }
    else{
      AskForTime = xfe_True;
      xf_printf("Will show physical mesh at requested time\n");
    }
    MeshMotionInPhysical = !MeshMotionInPhysical;

    ierr = xf_Error(xf_MovePlotMesh(AskForTime));
    if (ierr != xf_OK) exit(ierr);
    orig = glutGetWindow();
    glutSetWindow(WinAux.Handle);
    glutPostRedisplay();
    if (orig > 0) glutSetWindow(orig);

    if (RenderFlag){
      ierr = xf_Error(xf_FillDScalarsAll(DVector + CurrentDVector, xfe_False));
      if (ierr != xf_OK) exit(ierr);
    }

  }

  //I would like to use "T" capitalized to write tecplot file for ActiveDGroup
  //Not
  if (k == 'T')
  {

     //if not scalar field shown; ignor it
     if(RenderFlag == xfe_False){
     
        xf_printf("No scalar field for data dump!Please turn render on\n");
     }
     else
     {
        ierr = xf_Error(Yu_DumpAllDataToTecplot(DVector + CurrentDVector));
        if (ierr != xf_OK) exit(ierr);
        xf_printf("Tecplot-formatted data dump successfully!\n");
     }

  }



}


/******************************************************************/
//   FUNCTION Definition: xf_Key2D
void 
xf_Key2D (unsigned char k, int x, int y) {
  int ierr;

  xf_KeyCommon(k, x, y);

  if (k == '?') xf_DisplayKeyMenu2D();

  if ((k >= 'x') && (k <= 'z')){
    UpAxis2D = k - 'x';
    if ((iDGroupCutPlane >= 0) && (DGroup[iDGroupCutPlane].Active)){
      ierr = xf_Error(xf_CalculateCutPlaneCoords(DGroup + iDGroupCutPlane));
      if (ierr != xf_OK) exit(ierr);
      xf_Calculate2DWindowSize();
      xf_Resize(WinParent.PixWidth[0], WinParent.PixWidth[1]);
    }
  }
  
  if (k == 'h'){
    ierr = xf_Error(xf_HardcopyWindowEPS(Win2D, "plot.eps", xfe_True));
    if (ierr != xf_OK) exit(ierr);
    xf_printf("Wrote plot.eps.\n");
  }

  if (k == 'H'){
    ierr = xf_MakeMovie(2);
    if (ierr == xf_OK) xf_printf("Wrote animation.\n");
  }

  if (k == 'l'){
    ierr = xf_Error(xf_ClipCutPlane());
    if (ierr != xf_OK) exit(ierr);
  }
  
  xf_PostRedisplayAll();
}

/******************************************************************/
//   FUNCTION Definition: xf_Key3D
void 
xf_Key3D (unsigned char k, int x, int y) {
  int i, ierr;
  char buf[xf_MAXSTRLEN];

  xf_KeyCommon(k, x, y);

  if (k == '?') xf_DisplayKeyMenu3D();

  if ((k >= 'x') && (k <= 'z')){
    Camera.ViewAxis = k-'x';
    for (i=0; i<3; i++) Camera.Origin[i] = 0.;
    Camera.Azimuth = 0.0;
    Camera.Elevation = 0.0;
  }

  if (k == 'h'){
    ierr = xf_Error(xf_HardcopyWindowEPS(Win3D, "plot.eps", xfe_True));
    if (ierr != xf_OK) exit(ierr);
    xf_printf("Wrote plot.eps.\n");
  }

  if (k == 'H'){
    ierr = xf_MakeMovie(3);
    if (ierr == xf_OK) xf_printf("Wrote animation.\n");
  }

  if (k == 'l'){
    LightingFlag = !LightingFlag;
    if (LightingFlag){
      xf_printf("ambient light = [%.3f %.3f %.3f %.3f].  Enter new values:\n", 
		light_ambient[0], light_ambient[1], light_ambient[2], light_ambient[3]);
      if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	  (sscanf(buf, "%f %f %f %f", light_ambient+0, light_ambient+1, 
		  light_ambient+2, light_ambient+3) == 4))
	xf_printf("New ambient light set.\n");
      else
	xf_printf("Keeping original values for ambient light.\n");
      
      ierr = xf_Error(xf_FillLighting());
      if (ierr != xf_OK) exit(ierr);
      xf_printf("Lighting has been turned on.\n");
    }
    if (!LightingFlag){
      xf_printf("Lighting has been turned off.\n");
    }
    xf_ModifySubMenu_Lighting();
  }

  
  xf_PostRedisplayAll();
}

/******************************************************************/
//   FUNCTION Definition: xf_KeyAux
void 
xf_KeyAux (unsigned char k, int x, int y) {

  xf_KeyCommon(k, x, y);
  
  xf_PostRedisplayAll();
}


/******************************************************************/
//   FUNCTION Definition: xf_PositionCutPlane
static int 
xf_PositionCutPlane()
{  
  int ierr, k;
  char buf[xf_MAXSTRLEN];
  real n[3], d, NN;

  xf_printf("Set cut-plane normal (enter = no change) [%.4E %.4E %.4E]:\n",
	    CutPlane[0], CutPlane[1], CutPlane[2]);
  
  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
      (sscanf(buf, "%lf %lf %lf", n+0, n+1, n+2) == 3)){
    for (k=0, NN=0.; k<3; k++) NN += n[k]*n[k];
    NN = sqrt(NN);
    if (NN < MEPS){
      n[0] = 1.0; n[1] = 0.0; n[2] = 0.0;
      xf_printf("Zero-length normal provided.  Setting to [1,0,0].\n");
    }
    else{
      for (k=0; k<3; k++) n[k] /= NN;
    }
    CutPlane[0] = n[0];
    CutPlane[1] = n[1];
    CutPlane[2] = n[2];
  }
  else
    xf_printf("No change to cut plane normal.\n");
  
  xf_printf("Set cut-plane distance (enter = no change) [%.4E]:\n",
	    -CutPlane[3]);
  
  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
      (sscanf(buf, "%lf", &d) == 1)){
    CutPlane[3] = -d;
  }
  else
    xf_printf("No change to cut plane distance.\n");
  
  // re-calculate cut-plane
  if (iDGroupCutPlane >= 0){
    ierr = xf_Error(xf_CalculateDGroup(DGroup + iDGroupCutPlane));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SetIsoSurf
static int 
xf_SetIsoSurf()
{  
  int ierr, k;
  char buf[xf_MAXSTRLEN];
  real val;

  if (iDGroupIsoSurf < 0) return xf_OK; // nothing to do

  DGroup[iDGroupIsoSurf].Active = !DGroup[iDGroupIsoSurf].Active;
  xf_SetActiveDGroups();

  if (!DGroup[iDGroupIsoSurf].Active){ // just turned off
    xf_printf("Iso-surface has been turned off.\n");
    return xf_OK; 
  }

  xf_printf("Set iso-surface value (one real number):\n");
  
  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
      (sscanf(buf, "%lf", &val) == 1)){
    DGroup[iDGroupIsoSurf].isoval = val;
  }
  else
    xf_printf("Not understood.\n");
  

  // re-calculate iso-surface
  ierr = xf_Error(xf_CalculateDGroup(DGroup + iDGroupIsoSurf));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FindCutPlaneElem
static int 
xf_FindCutPlaneElem(real *xplane, enum xfe_Bool Verbose, real *xref,
		    int *pegrp, int *pelem)
{
  int ierr, iTri, k, iv[3], d;
  int egrp = -1, elem = -1;
  enum xfe_Bool inside;
  xf_DisplayTri *DTri;
  real *xg, *x[3], v[2], vp[2], J[4], iJ[4];
  real xi[2], phi[3];

  if (iDGroupCutPlane < 0) return xf_OK; // nothing to do
  if ((DTri = DGroup[iDGroupCutPlane].DTri) == NULL) return xf_OK; // nothing to do
  
  if (All->Mesh->Dim == 3)
    xg = DTri->xplane;
  else
    xg = DTri->xglob;

  // Brute force search
  for (iTri=0; iTri<DTri->nTri; iTri++){
    for (k=0; k<3; k++) iv[k] = DTri->TriList[3*iTri+k];
    for (k=0; k<3; k++) x[k] = xg + 3*iv[k];
    
    // check if xplane is inside triangle
    inside = xfe_True;
    for (k=0; k<3; k++){
      v[0] = x[(k+1)%3][0] - x[k][0];
      v[1] = x[(k+1)%3][1] - x[k][1];
      vp[0] = xplane[0] - x[k][0];
      vp[1] = xplane[1] - x[k][1];
      if ((v[0]*vp[1]-vp[0]*v[1]) < 0){
	inside = xfe_False;
	break;
      }
    } // k
 
    if (inside){
      /*       xf_printf("xplane = %.6E %.6E\n", xplane[0], xplane[1]); */
      /*       xf_printf("x[0] = %.6E %.6E\n", x[0][0], x[0][1]); */
      /*       xf_printf("x[1] = %.6E %.6E\n", x[1][0], x[1][1]); */
      /*       xf_printf("x[2] = %.6E %.6E\n", x[2][0], x[2][1]); */
      egrp = DTri->Tri2Elem[2*iTri+0];
      elem = DTri->Tri2Elem[2*iTri+1];
      break;
    }
    
  } // iTri

  if (pegrp != NULL) (*pegrp) = egrp;
  if (pelem != NULL) (*pelem) = elem;

  if (Verbose){
    if (egrp == -1)
      xf_printf("Element not found.\n");
    else{
      xf_printf(" egrp = %d, elem = %d\n", egrp, elem);
      J[0] = x[1][0]-x[0][0];  J[1] = x[2][0]-x[0][0];
      J[2] = x[1][1]-x[0][1];  J[3] = x[2][1]-x[0][1];
      ierr = xf_Error(xf_MatDetInv(J, 2, NULL, iJ));
      if (ierr != xf_OK) return ierr;
      v[0] = xplane[0] - x[0][0];
      v[1] = xplane[1] - x[0][1];
      xf_MxV_Set(iJ, v, 2, 2, xi);
      phi[0] = 1 - xi[0] - xi[1];
      phi[1] = xi[0];
      phi[2] = xi[1];
      for (d=0; d<All->Mesh->Dim; d++)
	for (k=0, xref[d]=0.; k<3; k++)
	  xref[d] += DTri->xref[3*iv[k]+d]*phi[k];
      if (All->Mesh->Dim == 2) xref[2] = 0.;
      xf_printf(" xref = %.6E %.6E %.6E\n", xref[0], xref[1], xref[2]);
    }
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_Pix2GeomCoord
static int 
xf_Pix2GeomCoord(xf_WindowType Win, int px, int py, real *xgeom)
{  
  real x, y;

  x = ((real) px)/((real) Win.PixWidth[0]);
  y = ((real) py)/((real) Win.PixWidth[1]);
  xgeom[0] = Win.Center[0] + Win.Width[0]*(x-0.5);
  xgeom[1] = Win.Center[1] + Win.Width[1]*(0.5-y);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_PrintElemInfo
static int 
xf_PrintElemInfo2D(int px, int py)
{  
  int ierr, i, k;
  int egrp, elem;
  real xglob[3], xplane[2];
  int *IParam = NULL;
  real *RParam = NULL;
  real xref[3], *v = NULL, *gv = NULL, *w = NULL;
  xf_BasisData *PhiData;
  xf_JacobianData *JData;
  xf_DisplayVector *DV;
  xf_Vector *V;
  xf_EqnSet *EqnSet;

  
  ierr = xf_Error(xf_Pix2GeomCoord(Win2D, px, py, xplane));
  if (ierr != xf_OK) return ierr;

  // corresponding glob coord
  ierr = xf_Error(xf_Plane2GlobCoord(xplane, xglob));
  if (ierr != xf_OK) return ierr;

  xf_printf("Point plane coords: %.10E %.10E %.10E\n", 
	    xplane[0], xplane[1], xplane[2]);
  xf_printf("Point global coords: %.10E %.10E %.10E\n", 
	    xglob[0], xglob[1], xglob[2]);

  // locate element
  ierr = xf_Error(xf_FindCutPlaneElem(xplane, xfe_True, xref, &egrp, &elem));
  if (ierr != xf_OK) return ierr;


  // print out scalar values if rendering
  if ((RenderFlag) && (egrp >= 0) && (elem >= 0)) {
    DV = DVector + CurrentDVector;
    V = DV->Vector;

    EqnSet = All->EqnSet;

    // pull off IParam/RParam if using a derived vector
    if ((DV->DeriveType == xfe_VectorDeriveScalar) ||
	(DV->DeriveType == xfe_VectorDeriveVariableSet)){
      ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
      if (ierr != xf_OK) return ierr;
    }

    // allocate memory
    ierr = xf_Error(xf_Alloc((void **) &v, 4*max(DV->nScalar, V->StateRank), sizeof(real)));
    if (ierr != xf_OK) return ierr;
    gv = v + max(DV->nScalar, V->StateRank);
    ierr = xf_Error(xf_Alloc((void **) &w, max(DV->nScalar, V->StateRank), sizeof(real)));
    if (ierr != xf_OK) return ierr;


    // interpolate V at xref
    PhiData = NULL;
    JData = NULL;

    ierr = xf_Error(xf_InterpolateScalar(DV, egrp, elem, V->Basis[egrp], xf_InterpOrder(V,egrp,elem), 
					 V->StateRank, V->GenArray[egrp].rValue[elem], 
					 EqnSet, 1, xref, &PhiData, &JData, IParam, 
					 RParam, v, gv, w));
    if (ierr != xf_OK) return ierr;
    
    // print out scalars
    for (k=0; k<DV->nScalar; k++)
      xf_printf("%20s = %.15E\n", DV->Scalar[k].Name, w[k]);

    // destroy basisdata
    ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
    if (ierr != xf_OK) return ierr;

    /* Destroy geometry Jacobian Data */
    ierr = xf_Error(xf_DestroyJacobianData(JData));
    if (ierr != xf_OK) return ierr;

    // release memory
    xf_Release( (void *) v);
    xf_Release( (void *) w);
    

  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_IntersectLineProbe
static int 
xf_IntersectLineProbe(int *npoint, int **pnvec, real **pslocvec, int **pevec)
{
  // intersects the line probe with the triangulation
  int ierr, i, j, k;
  int iTri, nedge;
  int elem, egrp;
  int n1, n2;
  int eg[3], el[3];
  int *nvec, *evec, *pos;
  int *node2ehash;
  xf_DisplayTri *DTri;
  real n[2], c, NN;
  real v1[2], v2[2];
  real *x1, *x2, *xg;
  real d1, d2;
  real sloc, sglob, rtemp;
  real *slocvec, *sglobvec;
  real xplane[2];
  xf_EdgeHash *edgehash;


  if ((!LineProbeOn) || (iDGroupCutPlane < 0)  || ((DTri = DGroup[iDGroupCutPlane].DTri) == NULL)){
    xf_printf("No cutplane or associated triangles for line probe intersection.\n");
    return xf_OK;
  }
  
  nedge = 2*DTri->nTri; // over-estimate

  // allocate edge hash: node2ehash + edgehash
  ierr = xf_Error(xf_Alloc( (void **) &node2ehash, DTri->nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<DTri->nNode; i++) node2ehash[i] = -1;
  ierr = xf_Error(xf_Alloc( (void **) &edgehash, nedge, sizeof(xf_EdgeHash)));
  if (ierr != xf_OK) return ierr;

  // over allocate memory
  ierr = xf_Error(xf_Alloc( (void **) &nvec, 2*nedge, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &evec, 2*nedge, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &slocvec, nedge, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &sglobvec, nedge, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  nedge = 0; // counter of number of intersected edges


  // get line into form a*x + b*y + c = 0
  // (a,b) = normal vector components
  // c = normal dot (-x1)
  n[0] = -LineProbeXY[3] + LineProbeXY[1];
  n[1] =  LineProbeXY[2] - LineProbeXY[0];
  NN = sqrt(n[0]*n[0]+n[1]*n[1]);
  n[0] /= NN; n[1] /= NN;

  c = -n[0]*LineProbeXY[0] - n[1]*LineProbeXY[1];

  // global coordinates
  xg = (All->Mesh->Dim == 3) ? DTri->xplane : DTri->xglob;

  // Brute force search
  for (iTri=0; iTri<DTri->nTri; iTri++){
    
    // element info
    egrp = DTri->Tri2Elem[2*iTri+0];
    elem = DTri->Tri2Elem[2*iTri+1];

    // loop over edges
    for (k=0; k<3; k++){
      // node numbers (ordered)
      n1 = DTri->TriList[3*iTri+ k     ];
      n2 = DTri->TriList[3*iTri+(k+1)%3];
      if (n1>n2) swap(n1,n2,i);
      // corresponding plane coordinates
      x1 = xg + 3*n1;
      x2 = xg + 3*n2;
      // distances from line
      d1 = n[0]*x1[0] + n[1]*x1[1] + c;
      d2 = n[0]*x2[0] + n[1]*x2[1] + c;

      if (d1*d2 < 0.){ // intersection occurs
	// svalue of the intersection
	sloc = fabs(d1)/(fabs(d1) + fabs(d2));
	
	// intersection plane coord
	for (i=0; i<2; i++) xplane[i] = (1.-sloc)*x1[i] + sloc*x2[i];

	// check that intersection is between line endpoints
	for (i=0; i<2; i++) v1[i] = xplane[i] - LineProbeXY[i];
	for (i=0; i<2; i++) v2[i] = LineProbeXY[2+i] - xplane[i];
	if ((v1[0]*v2[0] + v1[1]*v2[1]) > 0.){

	  // global s value
	  for (i=0, d1=0.; i<2; i++) d1 += (xplane[i]-LineProbeXY[  i])*(xplane[i]-LineProbeXY[  i]);
	  for (i=0, d2=0.; i<2; i++) d2 += (xplane[i]-LineProbeXY[2+i])*(xplane[i]-LineProbeXY[2+i]);
	  d1 = sqrt(d1);
	  d2 = sqrt(d2);
	  sglob = d1/(d1+d2);

	  // add to hash
	  i = nedge;
	  ierr = xf_Error(xf_EdgeHashAdd(n1, n2, node2ehash, edgehash, &nedge));
	  if (ierr != xf_OK) return ierr;
	  if (i != nedge){ // edge was added on this pass
	    slocvec [i] = sloc;
	    sglobvec[i] = sglob;
	    evec[2*i+0] = egrp;
	    evec[2*i+1] = elem;
	  }
	} // end if valid intersection
      }

    } // end for k over edges
    
  } // iTri

  (*npoint) = nedge;

  // allocate memory
  ierr = xf_Error(xf_Alloc( (void **) pnvec, 2*nedge, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) pevec, 2*nedge, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) pslocvec, nedge, sizeof(real)));
  if (ierr != xf_OK) return ierr;


  ierr = xf_Error(xf_Alloc( (void **) &pos, nedge, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // sort edge intersections according to sglobvec
  ierr = xf_Error(xf_SortRealParallel(sglobvec, nedge, xfe_False, pos));
  if (ierr != xf_OK) return ierr;

  // build nvec and slocvec (list of n1, n2)
  for (i=0; i<nedge; i++){
    j = pos[i];  // j is post-sorting rank of intersection
    (*pnvec)[2*j + 0] = edgehash[i].n0;
    (*pnvec)[2*j + 1] = edgehash[i].n1;
    (*pslocvec)[j] = slocvec[i];
    (*pevec)[2*j+0] = evec[2*i+0];
    (*pevec)[2*j+1] = evec[2*i+1];
  }

  //check that element edges are together
  for (i=0; i<nedge-2; i++){
    // pull off groups and elems of three nodes in a row
    for (j=0; j<3; j++){
      eg[j] = (*pevec)[2*(i+j)+0];
      el[j] = (*pevec)[2*(i+j)+1];
    }
    if ( ((eg[0] == eg[2]) && (el[0] == el[2])) &&
	 ((eg[0] != eg[1]) || (el[0] != el[1])) ){
      // we need to swap 1 and 2
      swap( (*pevec)[2*(i+1)+0], (*pevec)[2*(i+2)+0], k);
      swap( (*pevec)[2*(i+1)+1], (*pevec)[2*(i+2)+1], k);
      swap( (*pnvec)[2*(i+1)+0], (*pnvec)[2*(i+2)+0], k);
      swap( (*pnvec)[2*(i+1)+1], (*pnvec)[2*(i+2)+1], k);
      swap( (*pslocvec)[i+1], (*pslocvec)[i+2], rtemp);
    }
  } // j
  
  // line probe does not include the starts + ends (can add this later)
  // locate end of line
  /*   ierr = xf_Error(xf_FindCutPlaneElem(LineProbeXY+2, &egrp, &elem, xfe_False)); */
  /*   if (ierr != xf_OK) return ierr; */
  /*   if ((egrp >= 0) && ((*npoint) > 0)) // end exists, remove last point */
  /*     (*npoint)--; */


  xf_Release( (void *) node2ehash);
  xf_Release( (void *) edgehash);
  xf_Release( (void *) evec);
  xf_Release( (void *) slocvec);
  xf_Release( (void *) sglobvec);
  xf_Release( (void *) pos);
    
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_WriteLineData
static int 
xf_WriteLineData()
{  
  int ierr, i, k;
  int n1, n2;
  int ipoint, npoint;
  int egrp, elem;
  int *nvec = NULL;
  int *evec = NULL;
  xf_DisplayTri *DTri;
  xf_DisplayVector *DV;
  real *slocvec = NULL;
  real xplane[2], xglob[3];
  real sloc, sglob, d1, d2;
  real *x1, *x2, *xg;
  real val1, val2, val;
  FILE *fp;

  if (!LineProbeOn){
    xf_printf("Line probe is not on.  Not writing data.\n");
    return xf_OK;
  }
  if ((iDGroupCutPlane < 0)  || ((DTri = DGroup[iDGroupCutPlane].DTri) == NULL)){
    xf_printf("No cutplane or associated triangles for line probe intersection.\n");
    return xf_OK;
  }

  // intersect the line probe with the display triangles (result is sorted)
  ierr = xf_Error(xf_IntersectLineProbe(&npoint, &nvec, &slocvec, &evec));
  if (ierr != xf_OK) return ierr;

  // global coordinates
  xg = (All->Mesh->Dim == 3) ? DTri->xplane : DTri->xglob;

  // dump info to file "line.txt"
  if ((fp = fopen("line.txt", "w")) == NULL) return xf_Error(xf_FILE_WRITE_ERROR);

  // current vector and current scalar
  DV = DVector + CurrentDVector;
  
  // print out header
  fprintf(fp, "%% %21s %21s %21s  %21s ", "x", "y", "z", "s");
  if (RenderFlag){    // print out header for scalars
    for (k=0; k<DV->nScalar; k++){
      fprintf(fp, "%21s",  DV->Scalar[k].Name);
    }
  }
  fprintf(fp, "\n");

  // loop over all intersecting points
  for (ipoint = 0; ipoint < npoint; ipoint++){
    n1 = nvec[2*ipoint+0];    // node 1
    n2 = nvec[2*ipoint+1];    // node 2
    sloc = slocvec[ipoint];   // position between n1 and n2, [0,1]
    
    egrp = evec[2*ipoint+0];  // element group
    elem = evec[2*ipoint+1];  // element number
    
    // plane coordinates of the two nodes
    x1 = xg + 3*n1; // node 1
    x2 = xg + 3*n2; // node 2

    // intersection plane coord
    for (i=0; i<2; i++) xplane[i] = (1.-sloc)*x1[i] + sloc*x2[i];
    
    // line global s-value
    for (i=0, d1=0.; i<2; i++) d1 += (xplane[i]-LineProbeXY[  i])*(xplane[i]-LineProbeXY[  i]);
    for (i=0, d2=0.; i<2; i++) d2 += (xplane[i]-LineProbeXY[2+i])*(xplane[i]-LineProbeXY[2+i]);
    d1 = sqrt(d1);
    d2 = sqrt(d2);
    if ((d1+d2) == 0.){
      xf_printf("Warning, zero-length line!\n");
      sglob = -1.0;
    }
    else sglob = d1/(d1+d2);

    // corresponding glob coord
    ierr = xf_Error(xf_Plane2GlobCoord(xplane, xglob));
    if (ierr != xf_OK) return ierr;

    fprintf(fp, "%21.12E %21.12E %21.12E  %21.12E", 
	    xglob[0], xglob[1], xglob[2], sglob);

    // print out any scalars
    if (RenderFlag){
      for (k=0; k<DV->nScalar; k++){
	val1 = DV->Scalar[k].Data[iDGroupCutPlane].Value[n1];
	val2 = DV->Scalar[k].Data[iDGroupCutPlane].Value[n2];
	val = (1.-sloc)*val1 + sloc*val2;
	fprintf(fp, "% 21.12E", val);
      }
    }

    fprintf(fp, "\n");

  } // ipoint


  // close file
  fclose(fp);

  xf_printf("Wrote line.txt with line probe data.\n");

  // Release memory
  xf_Release( (int *) nvec);
  xf_Release( (int *) slocvec);
  xf_Release( (int *) evec);


  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CountElemsInBox
static int 
xf_CountElemsInBox()
{  
  int ierr, count;
  char buf[xf_MAXSTRLEN];
  int egrp, elem;
  enum xfe_ShapeType Shape;
  real xref_tri[2]  = {0.33333, 0.33333};
  real xref_quad[2] = {0.5, 0.5};
  real *xref, xglob[3];
  xf_BasisData *PhiData;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  if (Mesh->Dim != 2){
    xf_printf("Elem Count only supported in 2D meshes.\n");
    return xf_OK;
  }

  xf_printf("[xmin xmax ymin ymax] = %.6E %.6E %.6E %.6E\n",
	    CountBox[0], CountBox[1], CountBox[2], CountBox[3]);
  
  xf_printf("Enter new (or blank to keep default)\n");

  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
      (sscanf(buf, "%lf %lf %lf %lf", CountBox, CountBox+1,
	      CountBox+2, CountBox+3) == 4)){
    xf_printf("Set new box coordinates.\n");
  }
  else
    xf_printf("Not changing box coordinates.\n");

  count = 0;
  PhiData = NULL;

  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
    if (ierr != xf_OK) return ierr;
    
    xref = ((Shape == xfe_Triangle) ? xref_tri : xref_quad);

    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      ierr = xf_Error(xf_Ref2GlobElem(All->Mesh, egrp, elem, 
				      &PhiData, xfe_True, 1, xref, 
				      xglob));
    
      if ((xglob[0] > CountBox[1]) || (xglob[0] < CountBox[0]) ||
	  (xglob[1] > CountBox[3]) || (xglob[1] < CountBox[2]))
	continue;

      count++;
    } // elem
  } // egrp

  xf_printf("# elems in box = %d\n", count);
  
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
 

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SpecialKeyCommon
static void 
xf_SpecialKeyCommon(int k, int x, int y) {
  int ierr;
}

/******************************************************************/
//   FUNCTION Definition: xf_SpecialKey2D
void 
xf_SpecialKey2D(int k, int x, int y) {
  int ierr;

  xf_SpecialKeyCommon(k, x, y);

  if (k == GLUT_KEY_F5){
    // acquire line probe
    if (LineProbeOn) LineProbeOn = xfe_False; // turn off if on
    AcquiringLineProbe = !AcquiringLineProbe;
    nLineProbeAcquired = 0.;
    if (AcquiringLineProbe)
      xf_printf("Select 2 points in the 2D window using the right mouse button.\n");
    else
      xf_printf("No longer acquiring the line probe.\n");
  }

  if (k == GLUT_KEY_F6){
    // write line probe data
    ierr = xf_Error(xf_WriteLineData());
    if (ierr != xf_OK) exit(ierr);
  }

  if (k == GLUT_KEY_F7){
    // toggle line probe
    LineProbeOn = !LineProbeOn;
  }

  if (k == GLUT_KEY_F8){
    CountBoxOn = !CountBoxOn;
    if (CountBoxOn){
      // count number of elements within a box
      ierr = xf_Error(xf_CountElemsInBox());
      if (ierr != xf_OK) exit(ierr);
    }
  }


  xf_PostRedisplayAll();

}


/******************************************************************/
//   FUNCTION Definition: xf_SpecialKey3D
void 
xf_SpecialKey3D(int k, int x, int y) {
  int ierr;

  xf_SpecialKeyCommon(k, x, y);

  if (k == GLUT_KEY_F3){
    // cut plane positioning
    ierr = xf_Error(xf_PositionCutPlane());
    if (ierr != xf_OK) exit(ierr);
    if ((iDGroupCutPlane >= 0) && (DGroup[iDGroupCutPlane].Active)){
      if (RenderFlag){
	ierr = xf_Error(xf_FillDScalarsAll(DVector + CurrentDVector, xfe_True));
	if (ierr != xf_OK) exit(ierr);
      }
      if (LightingFlag){
	ierr = xf_Error(xf_FillLighting());
	if (ierr != xf_OK) exit(ierr);
      }
      ierr = xf_Error(xf_CalculateCutPlaneCoords(DGroup + iDGroupCutPlane));
      if (ierr != xf_OK) exit(ierr);
      xf_Calculate2DWindowSize();
      xf_Resize(WinParent.PixWidth[0], WinParent.PixWidth[1]);
    }    
  }

  if (k == GLUT_KEY_F4){
    // iso-surface setting
    ierr = xf_Error(xf_SetIsoSurf());
    if (ierr != xf_OK) exit(ierr);
    if ((iDGroupIsoSurf >= 0) && (DGroup[iDGroupIsoSurf].Active)){
      if (RenderFlag){
	ierr = xf_Error(xf_FillDScalarsAll(DVector + CurrentDVector, xfe_True));
	if (ierr != xf_OK) exit(ierr);
      }
      if (LightingFlag){
	ierr = xf_Error(xf_FillLighting());
	if (ierr != xf_OK) exit(ierr);
      }
    }  
    xf_ModifySubMenu_IsoSurface();
  }

  if (k == GLUT_KEY_F8){
    // hills plotting
    if (All->Mesh->Dim == 2){
      DrawingHills = !DrawingHills;
      if (DrawingHills) xf_printf("DrawingHills on\n");
      else xf_printf("DrawingHills off\n");
      if ((iDGroupCutPlane >= 0) && (DGroup[iDGroupCutPlane].Active)){
	ierr = xf_Error(xf_CalculateDGroup(DGroup + iDGroupCutPlane));
	if (ierr != xf_OK) exit(ierr);
	if (RenderFlag){
	  ierr = xf_Error(xf_FillDScalarsAll(DVector + CurrentDVector, xfe_False));
	  if (ierr != xf_OK) exit(ierr);
	}
	if (LightingFlag){
	  ierr = xf_Error(xf_FillLighting());
	  if (ierr != xf_OK) exit(ierr);
	}
      }    
    }
  }




  xf_PostRedisplayAll();
}

/******************************************************************/
//   FUNCTION Definition: xf_SpecialKeyAux
void 
xf_SpecialKeyAux(int k, int x, int y) {
  int ierr;

  xf_SpecialKeyCommon(k, x, y);
  xf_PostRedisplayAll();
}


/******************************************************************/
//   FUNCTION Definition: xf_Mouse2D
void 
xf_Mouse2D (int button, int state, int x, int y) {

  int ierr;

  switch (button){
  case GLUT_LEFT_BUTTON:
    MouseState.LeftButtonDown = (state == GLUT_DOWN);
    break;
  case GLUT_MIDDLE_BUTTON:
    MouseState.MiddleButtonDown = (state == GLUT_DOWN);
    break;
  case GLUT_RIGHT_BUTTON:
    MouseState.RightButtonDown = (state == GLUT_DOWN);
    if ((MouseState.RightButtonDown) && (AcquiringLineProbe)){
      ierr = xf_Error(xf_Pix2GeomCoord(Win2D, x, y, LineProbeXY+nLineProbeAcquired*2));
      if (ierr != xf_OK) exit(ierr);
      nLineProbeAcquired++;
      xf_printf("Got point %d of the line probe.\n", nLineProbeAcquired);
      if (nLineProbeAcquired == 2){
	AcquiringLineProbe = xfe_False;
	xf_printf("Finished acquiring the line probe.\n");
	LineProbeOn = xfe_True;
	xf_ModifySubMenu_LineProbe();
      }
    }
    else if (MouseState.RightButtonDown){
      // on right-button down press, print out elem info
      ierr = xf_Error(xf_PrintElemInfo2D(x, y));
      if (ierr != xf_OK) exit(ierr);
    }
    break;
  default:
    break;
  }
  MouseState.x = x;
  MouseState.y = y;
  MouseState.Modifier = glutGetModifiers();

  xf_PostRedisplayAll();
}


/******************************************************************/
//   FUNCTION Definition: xf_Mouse3D
void 
xf_Mouse3D (int button, int state, int x, int y) {

  switch (button){
  case GLUT_LEFT_BUTTON:
    MouseState.LeftButtonDown = (state == GLUT_DOWN);
    break;
  case GLUT_MIDDLE_BUTTON:
    MouseState.MiddleButtonDown = (state == GLUT_DOWN);
    if ((state == GLUT_DOWN) && (RenderFlag)){
      RenderFlag = xfe_False;
      RenderTemporarilyOff = xfe_True;
    }
    else if ((state == GLUT_UP) && (RenderTemporarilyOff)){
      RenderFlag = xfe_True;
      RenderTemporarilyOff = xfe_False;
      glutPostRedisplay();
    }
    break;
  case GLUT_RIGHT_BUTTON:
    MouseState.RightButtonDown = (state == GLUT_DOWN);
    break;
  default:
    break;
  }
  MouseState.x = x;
  MouseState.y = y;
  MouseState.Modifier = glutGetModifiers();
}

/******************************************************************/
//   FUNCTION Definition: xf_PassiveMotion
void 
xf_PassiveMotion ( int x, int y) {

  xf_Idle();

}

/******************************************************************/
//   FUNCTION Definition: xf_Motion2D
void 
xf_Motion2D ( int x, int y) {
  int dx, dy;
  real fac;

  // get deltas in x and y
  dx = x - MouseState.x;
  dy = y - MouseState.y;
 
  if (MouseState.MiddleButtonDown){
    if (MouseState.Modifier & GLUT_ACTIVE_SHIFT){
      // zoom using dy
      fac = 1.0+3.0*dy/Win2D.PixWidth[1];
      fac = min(fac, 1.5);
      fac = max(fac, 0.5);
      Win2D.Width[0] *= fac;
      Win2D.Width[1] *= fac;
    }
    else{
      // move left/right or up/down if middle mouse button is pressed
      Win2D.Center[0] -= (Win2D.Width[0]*dx)/Win2D.PixWidth[0];
      Win2D.Center[1] += (Win2D.Width[1]*dy)/Win2D.PixWidth[1];
    }
    xf_SetProjection2D(Win2D);
    glutPostRedisplay();
  }
  
  MouseState.x = x;
  MouseState.y = y;

}


/******************************************************************/
//   FUNCTION Definition: xf_Motion3D
void 
xf_Motion3D( int x, int y) {
  int dx, dy, ax, ay, az;
  real fac, dox, doy;
  real sa, ca, se, ce;

  // get deltas in x and y
  dx = x - MouseState.x;
  dy = y - MouseState.y;
 
  if (MouseState.MiddleButtonDown){
    if (MouseState.Modifier & GLUT_ACTIVE_SHIFT){
      // zoom using dy
      fac = 1.0+3.0*dy/Win3D.PixWidth[1];
      fac = min(fac, 1.5);
      fac = max(fac, 0.5);
      Camera.Distance *= fac;
      Win3D.Width[0] *= fac;
      Win3D.Width[1] *= fac;
    }
    else if (MouseState.Modifier & GLUT_ACTIVE_CTRL){
      Camera.Azimuth -= (2.0*(real) dx)/Win3D.PixWidth[0]*180.0/PI;
      Camera.Elevation -= (2.0*(real) dy)/Win3D.PixWidth[1]*180.0/PI;
    }
    else{
      az = 2; 
      ax = (az+1)%3;
      ay = (ax+1)%3;

      // move left/right or up/down if middle mouse button is pressed
      ca = cos(Camera.Azimuth*PI/180.0);
      sa = sin(Camera.Azimuth*PI/180.0);
      dox = Win3D.Width[0]*dx/Win3D.PixWidth[0];//*Camera.Distance;
      Camera.Origin[ax] -= dox*ca;
      Camera.Origin[az] += dox*sa;
      ce = cos(Camera.Elevation*PI/180.0);
      se = sin(Camera.Elevation*PI/180.0);
      doy = Win3D.Width[1]*dy/Win3D.PixWidth[1];//*Camera.Distance;
      Camera.Origin[ax] += doy*sa*se;
      Camera.Origin[ay] += doy*ce;
      Camera.Origin[az] += doy*se*ca;
    }
    xf_SetProjection3D(Win3D);
    glutPostRedisplay();
  }

  // this call enables menu regeneration
  if (!MouseState.LeftButtonDown) xf_Idle();
  
  MouseState.x = x;
  MouseState.y = y;
}

/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_MouseWheel2D */
/* void  */
/* xf_MouseWheel2D ( int button, int dir, int x, int y) { */
/*   real fac; */

/*   fac = (dir > 0) ? 0.5 : 1.5; */
  
/*   Win2D.Width[0] *= fac; */
/*   Win2D.Width[1] *= fac; */
/*   xf_SetProjection2D(Win2D); */
/*   glutPostRedisplay(); */
/* } */

/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_MouseWheel3D */
/* void  */
/* xf_MouseWheel3D ( int button, int dir, int x, int y) { */
/*   real fac; */

/*   fac = (dir > 0) ? 0.5 : 1.5; */
  
/*   Win3D.Width[0] *= fac; */
/*   Win3D.Width[1] *= fac; */
/*   xf_SetProjection3D(Win3D); */
/*   glutPostRedisplay(); */
/* } */

/******************************************************************/
//   FUNCTION Definition: xf_InitWindows
static void
xf_InitWindows()
{
  Win2D.Handle = -1;
  Win3D.Handle = -1;
  WinAux.Handle = -1;
  WinBorder.Handle = -1;
}

static void selectMessage(int value){}
static void selectColor(int value){
  xf_printf("in selectColor: value = %d\n", value);
}

/* Function definitions for menu modification and processing */

/******************************************************************/
//   FUNCTION Definition: xf_ModifySubMenu_LoadVector
static void
xf_ModifySubMenu_LoadVector()
{
  int k, n;
  enum xfe_PlotMenuOptions opt;
  char **Names = NULL;

  opt = xfe_PlotMenu_Sub_LoadVector;
  // destroy existing submenu if it exists
  if (SubMenuIndex[opt] > 0) glutDestroyMenu(SubMenuIndex[opt]);
  // create new submenu
  SubMenuIndex[opt] = glutCreateMenu(xf_ProcessMenu_LoadVector);
  
  // add to submenu
  xf_GetVectorNames(&Names, &n);
  for (k=0; k<n; k++) glutAddMenuEntry(Names[k], k);
  if (n==0) glutAddMenuEntry("None", -1);

  xf_Release2((void **) Names);

  // attach submenu to top menu
  glutSetMenu(TopMenu);
  glutChangeToSubMenu(opt, "Change vector", SubMenuIndex[opt]);
  
}

/******************************************************************/
//   FUNCTION Definition: xf_ModifySubMenu_ChangeScalar
static void
xf_ModifySubMenu_ChangeScalar()
{
  int k;
  char s[xf_MAXSTRLEN];
  enum xfe_PlotMenuOptions opt;
  xf_DisplayVector *DV;

  DV = DVector + CurrentDVector;

  opt = xfe_PlotMenu_Sub_ChangeScalar;
  // destroy existing submenu if it exists
  if (SubMenuIndex[opt] > 0) glutDestroyMenu(SubMenuIndex[opt]);
  // create new submenu
  SubMenuIndex[opt] = glutCreateMenu(xf_ProcessMenu_ChangeScalar);
  
  // add to submenu
  if (DV != NULL){
    for (k=0; k<DV->nScalar; k++){
      sprintf(s, "%s", DV->Scalar[k].Name);
      if (DV->CurrentDScalar == k) strcat(s, " *");
      glutAddMenuEntry(s, k);
    }
  }
  else glutAddMenuEntry("None", -1);


  // attach submenu to top menu
  glutSetMenu(TopMenu);
  glutChangeToSubMenu(opt, "Change scalar", SubMenuIndex[opt]);
}

/******************************************************************/
//   FUNCTION Definition: xf_ModifySubMenu_Mesh
static void
xf_ModifySubMenu_Mesh()
{
  char s[xf_MAXSTRLEN];
  enum xfe_PlotMenuOptions opt;

  opt = xfe_PlotMenu_Sub_Mesh;
  // destroy existing submenu if it exists
  if (SubMenuIndex[opt] > 0) glutDestroyMenu(SubMenuIndex[opt]);
  // create new submenu
  SubMenuIndex[opt] = glutCreateMenu(xf_ProcessMenu_Mesh);
  
  // add to submenu
  sprintf(s, "%s Mesh",    (MeshOn   ) ? "Hide" : "Show");
  glutAddMenuEntry(s, 0);
  sprintf(s, "%s SubMesh", (SubMeshOn) ? "Hide" : "Show");
  glutAddMenuEntry(s, 1);
  glutAddMenuEntry("Thicker lines", 2);
  glutAddMenuEntry("Thinner lines", 3);

  glutSetMenu(TopMenu);
  glutChangeToSubMenu(opt, "Mesh", SubMenuIndex[opt]);
}


/******************************************************************/
//   FUNCTION Definition: xf_ModifySubMenu_Boundaries
static void
xf_ModifySubMenu_Boundaries()
{
  int i;
  char c;
  char s[xf_MAXSTRLEN];
  enum xfe_PlotMenuOptions opt;

  opt = xfe_PlotMenu_Sub_Boundaries;
  // destroy existing submenu if it exists
  if (SubMenuIndex[opt] > 0) glutDestroyMenu(SubMenuIndex[opt]);
  // create new submenu
  SubMenuIndex[opt] = glutCreateMenu(xf_ProcessMenu_Boundaries);
  
  // add to submenu
  if (All->Mesh->Dim == 2){
    sprintf(s, "%s boundary outline", (BoundaryOutline) ? "Turn off" : "Turn on");
    glutAddMenuEntry(s, 0);
  }
  else{
    for (i=0; i<All->Mesh->nBFaceGroup; i++){
      c = ((DGroup[iDGroupBFG[i]].Active) ? '*' : ' ');
      sprintf(s, "%s %c\n", All->Mesh->BFaceGroup[i].Title, c);
      glutAddMenuEntry(s, i);
    }
    glutAddMenuEntry("Toggle all", All->Mesh->nBFaceGroup);
  }

  glutSetMenu(TopMenu);
  glutChangeToSubMenu(opt, "Boundaries", SubMenuIndex[opt]);
}

/******************************************************************/
//   FUNCTION Definition: xf_ModifySubMenu_ColorScheme
static void
xf_ModifySubMenu_ColorScheme()
{
  int k;
  char s[xf_MAXSTRLEN];
  enum xfe_PlotMenuOptions opt;

  opt = xfe_PlotMenu_Sub_ColorScheme;
  // destroy existing submenu if it exists
  if (SubMenuIndex[opt] > 0) glutDestroyMenu(SubMenuIndex[opt]);
  // create new submenu
  SubMenuIndex[opt] = glutCreateMenu(xf_ProcessMenu_ColorScheme);
  
  // add to submenu
  for (k=0; k<xfe_RenderMapLast; k++){
    sprintf(s, "%s", xfe_RenderMapName[k]);
    if (RenderMap == k) strcat(s, " *");
    glutAddMenuEntry(s, k);
  }

  glutSetMenu(TopMenu);
  glutChangeToSubMenu(opt, "Color scheme", SubMenuIndex[opt]);
}


/******************************************************************/
//   FUNCTION Definition: xf_ModifySubMenu_CutPlane
static void
xf_ModifySubMenu_CutPlane()
{
  enum xfe_Bool CutPlaneOn;
  char s[xf_MAXSTRLEN];
  enum xfe_PlotMenuOptions opt;

  opt = xfe_PlotMenu_Sub_CutPlane;
  // destroy existing submenu if it exists
  if (SubMenuIndex[opt] > 0) glutDestroyMenu(SubMenuIndex[opt]);
  // create new submenu
  SubMenuIndex[opt] = glutCreateMenu(xf_ProcessMenu_CutPlane);
  
  CutPlaneOn = ((iDGroupCutPlane >= 0) && (DGroup[iDGroupCutPlane].Active));

  // add to submenu
  sprintf(s, "%s cut plane", (CutPlaneOn) ? "Hide" : "Show");
  glutAddMenuEntry(s, 0);
  glutAddMenuEntry("Position cut plane [i]", 1);

  glutSetMenu(TopMenu);
  glutChangeToSubMenu(opt, "Cut plane", SubMenuIndex[opt]);
}


/******************************************************************/
//   FUNCTION Definition: xf_ModifySubMenu_LineProbe
static void
xf_ModifySubMenu_LineProbe()
{
  char s[xf_MAXSTRLEN];
  enum xfe_PlotMenuOptions opt;

  opt = xfe_PlotMenu_Sub_LineProbe;
  // destroy existing submenu if it exists
  if (SubMenuIndex[opt] > 0) glutDestroyMenu(SubMenuIndex[opt]);
  // create new submenu
  SubMenuIndex[opt] = glutCreateMenu(xf_ProcessMenu_LineProbe);

  // add to submenu
  sprintf(s, "%s line probe", (LineProbeOn) ? "Hide" : "Show");
  glutAddMenuEntry(s, 0);
  glutAddMenuEntry("Acquire line probe", 1);

  glutSetMenu(TopMenu);
  glutChangeToSubMenu(opt, "Line probe", SubMenuIndex[opt]);
}

/******************************************************************/
//   FUNCTION Definition: xf_ModifySubMenu_Contours
static void
xf_ModifySubMenu_Contours()
{
  char s[xf_MAXSTRLEN];
  enum xfe_PlotMenuOptions opt;

  opt = xfe_PlotMenu_Sub_Contours;
  // destroy existing submenu if it exists
  if (SubMenuIndex[opt] > 0) glutDestroyMenu(SubMenuIndex[opt]);
  // create new submenu
  SubMenuIndex[opt] = glutCreateMenu(xf_ProcessMenu_Contours);

  // add to submenu
  sprintf(s, "%s contours", (ContoursOn) ? "Hide" : "Show");
  glutAddMenuEntry(s, 0);
  glutAddMenuEntry("Set number of contours [i]", 1);

  glutSetMenu(TopMenu);
  glutChangeToSubMenu(opt, "Contours", SubMenuIndex[opt]);
}

/******************************************************************/
//   FUNCTION Definition: xf_ModifySubMenu_IsoSurface
static void
xf_ModifySubMenu_IsoSurface()
{
  char s[xf_MAXSTRLEN];
  enum xfe_Bool IsoSurfaceOn;
  enum xfe_PlotMenuOptions opt;

  opt = xfe_PlotMenu_Sub_IsoSurface;
  // destroy existing submenu if it exists
  if (SubMenuIndex[opt] > 0) glutDestroyMenu(SubMenuIndex[opt]);
  // create new submenu
  SubMenuIndex[opt] = glutCreateMenu(xf_ProcessMenu_IsoSurface);

  IsoSurfaceOn = ((iDGroupIsoSurf >= 0) && (DGroup[iDGroupIsoSurf].Active));

  // add to submenu
  sprintf(s, "%s iso-surface", (IsoSurfaceOn) ? "Hide" : "Show");
  glutAddMenuEntry(s, 0);
  glutAddMenuEntry("Set iso-surface level [i]", 1);

  glutSetMenu(TopMenu);
  glutChangeToSubMenu(opt, "Iso-surface", SubMenuIndex[opt]);
}

/******************************************************************/
//   FUNCTION Definition: xf_ModifySubMenu_Lighting
static void
xf_ModifySubMenu_Lighting()
{
  char s[xf_MAXSTRLEN];
  enum xfe_PlotMenuOptions opt;

  opt = xfe_PlotMenu_Sub_Lighting;
  // destroy existing submenu if it exists
  if (SubMenuIndex[opt] > 0) glutDestroyMenu(SubMenuIndex[opt]);
  // create new submenu
  SubMenuIndex[opt] = glutCreateMenu(xf_ProcessMenu_Lighting);

  // add to submenu
  sprintf(s, "%s lighting", (LightingFlag) ? "Turn off" : "Turn on");
  glutAddMenuEntry(s, 0);
  glutAddMenuEntry("Set lighting properties [i]", 1);

  glutSetMenu(TopMenu);
  glutChangeToSubMenu(opt, "Lighting", SubMenuIndex[opt]);
}

/******************************************************************/
//   FUNCTION Definition: xf_ModifySubMenu_ViewFrom
static void
xf_ModifySubMenu_ViewFrom()
{
  enum xfe_PlotMenuOptions opt;

  opt = xfe_PlotMenu_Sub_ViewFrom;
  // destroy existing submenu if it exists
  if (SubMenuIndex[opt] > 0) glutDestroyMenu(SubMenuIndex[opt]);
  // create new submenu
  SubMenuIndex[opt] = glutCreateMenu(xf_ProcessMenu_ViewFrom);

  // add to submenu
  glutAddMenuEntry("From x axis", 0);
  glutAddMenuEntry("From y axis", 1);
  glutAddMenuEntry("From z axis", 2);

  glutSetMenu(TopMenu);
  glutChangeToSubMenu(opt, "Set 3D view", SubMenuIndex[opt]);
}

/******************************************************************/
//   FUNCTION Definition: xf_ModifySubMenu_Hardcopy
static void
xf_ModifySubMenu_Hardcopy()
{
  enum xfe_PlotMenuOptions opt;

  opt = xfe_PlotMenu_Sub_Hardcopy;
  // destroy existing submenu if it exists
  if (SubMenuIndex[opt] > 0) glutDestroyMenu(SubMenuIndex[opt]);
  // create new submenu
  SubMenuIndex[opt] = glutCreateMenu(xf_ProcessMenu_Hardcopy);

  // add to submenu
  glutAddMenuEntry("Write plot.eps", 0);
  glutAddMenuEntry("Make movie [i]", 1);

  glutSetMenu(TopMenu);
  glutChangeToSubMenu(opt, "Hardcopy", SubMenuIndex[opt]);
}


/******************************************************************/
//   FUNCTION Definition: xf_ModifySubMenu_ColorLimits
static void
xf_ModifySubMenu_ColorLimits()
{
  enum xfe_PlotMenuOptions opt;

  opt = xfe_PlotMenu_Sub_ColorLimits;
  // destroy existing submenu if it exists
  if (SubMenuIndex[opt] > 0) glutDestroyMenu(SubMenuIndex[opt]);
  // create new submenu
  SubMenuIndex[opt] = glutCreateMenu(xf_ProcessMenu_ColorLimits);

  // add to submenu
  glutAddMenuEntry("Default color limits", 0);
  glutAddMenuEntry("Set color limits [i]", 1);

  glutSetMenu(TopMenu);
  glutChangeToSubMenu(opt, "Color limits", SubMenuIndex[opt]);
}


/******************************************************************/
//   FUNCTION Definition: xf_ProcessMenu_Top
static void 
xf_ProcessMenu_Top(int value){
  int ierr, dim;

  switch (value){
  case xfe_PlotMenu_DeriveVector:
    xf_KeyCommon('d', 0, 0);
    break;
  case xfe_PlotMenu_SwapWindows:
    xf_KeyCommon(9, 0, 0);
    break;
  case xfe_PlotMenu_ToggleRender:
    RenderFlag = !RenderFlag;
    if (RenderFlag) xf_TurnRenderOn();
    xf_PostRedisplayAll();
    break;
  case xfe_PlotMenu_Help:
    xf_DisplayKeyMenuCommon();
    xf_DisplayKeyMenu2D();
    xf_DisplayKeyMenu3D();
    break;
  case xfe_PlotMenu_Quit:
    xf_DestroyGlobals();
    xf_printf("xf_Plot finished.\n"); fflush(stdout);
    exit(0);
    break;
  default:
    xf_printf("Unknown menu option.\n");
    break;
  }

}

/******************************************************************/
//   FUNCTION Definition: xf_ProcessMenu_LoadVector
static void
xf_ProcessMenu_LoadVector(int k)
{
  CurrentDVector = k;
  if (RenderFlag) xf_TurnRenderOn();
  // modify scalar and vector menus
  nModifyMenu = 2;
  ModifyMenuList[0] = xfe_PlotMenu_Sub_LoadVector;
  ModifyMenuList[1] = xfe_PlotMenu_Sub_ChangeScalar;
  
  xf_PostRedisplayAll();
}

/******************************************************************/
//   FUNCTION Definition: xf_ProcessMenu_ChangeScalar
static void
xf_ProcessMenu_ChangeScalar(int k)
{
  xf_DisplayVector *DV;
  DV = DVector + CurrentDVector;
  DV->CurrentDScalar = k;
  xf_WriteScalarInfo();
  nModifyMenu = 1;
  ModifyMenuList[0] = xfe_PlotMenu_Sub_ChangeScalar;

  xf_PostRedisplayAll();
}

/******************************************************************/
//   FUNCTION Definition: xf_ProcessMenu_Mesh
static void
xf_ProcessMenu_Mesh(int k)
{
  if (k==0) MeshOn    = !MeshOn;
  if (k==1) SubMeshOn = !SubMeshOn;
  if (k==2) MeshLineWidth = MeshLineWidth + 0.5;
  if (k==3) MeshLineWidth = max(0., MeshLineWidth-0.5);
  nModifyMenu = 1;
  ModifyMenuList[0] = xfe_PlotMenu_Sub_Mesh;
  
  xf_PostRedisplayAll();
}

/******************************************************************/
//   FUNCTION Definition: xf_ProcessMenu_Boundaries
static void
xf_ProcessMenu_Boundaries(int k)
{
  int ierr;
  int i, i0, i1;

  if (All->Mesh->Dim == 2){
    BoundaryOutline = !BoundaryOutline;
  }
  else{
    if (k>=All->Mesh->nBFaceGroup){
      i0 = 0; i1 = All->Mesh->nBFaceGroup-1;
    }
    else{
      i0 = k; i1 = k;
    }
    for (i=i0; i<=i1; i++){
      DGroup[iDGroupBFG[i]].Active = !DGroup[iDGroupBFG[i]].Active;
      if (DGroup[iDGroupBFG[i]].Active){
	ierr = xf_Error(xf_CalculateDGroup(DGroup + iDGroupBFG[i]));
	if (ierr != xf_OK) exit(ierr);
	if (RenderFlag){
	  ierr = xf_Error(xf_FillDScalarsAll(DVector + CurrentDVector, xfe_True));
	  if (ierr != xf_OK) exit(ierr);
	}
      }
    } // i
  }
  nModifyMenu = 1;
  ModifyMenuList[0] = xfe_PlotMenu_Sub_Boundaries;
  xf_PostRedisplayAll();
}

/******************************************************************/
//   FUNCTION Definition: xf_ProcessMenu_ColorScheme
static void
xf_ProcessMenu_ColorScheme(int k)
{
  int ierr;
  RenderMap = k;
  if (RenderFlag){
    ierr = xf_Error(xf_FillDScalarsAll(DVector + CurrentDVector, xfe_True));
    if (ierr != xf_OK) exit(ierr);
  }
  nModifyMenu = 1;
  ModifyMenuList[0] = xfe_PlotMenu_Sub_ColorScheme;

}

/******************************************************************/
//   FUNCTION Definition: xf_ProcessMenu_CutPlane
static void
xf_ProcessMenu_CutPlane(int k)
{
  switch (k){
  case 0:
    if (iDGroupCutPlane >= 0) 
      DGroup[iDGroupCutPlane].Active = !DGroup[iDGroupCutPlane].Active;
    break;
  case 1:
    xf_SpecialKey3D(GLUT_KEY_F3, 0, 0);
    break;
  default:
    xf_printf("Unknown menu item\n");
    break;
  }
  nModifyMenu = 1;
  ModifyMenuList[0] = xfe_PlotMenu_Sub_CutPlane;
  xf_PostRedisplayAll();
}


/******************************************************************/
//   FUNCTION Definition: xf_ProcessMenu_LineProbe
static void
xf_ProcessMenu_LineProbe(int k)
{
  switch (k){
  case 0:
    LineProbeOn = !LineProbeOn;
    break;
  case 1:
    xf_SpecialKey2D(GLUT_KEY_F5, 0, 0);
    break;
  default:
    xf_printf("Unknown menu item\n");
    break;
  }
  nModifyMenu = 1;
  ModifyMenuList[0] = xfe_PlotMenu_Sub_LineProbe;
  xf_PostRedisplayAll();
}


/******************************************************************/
//   FUNCTION Definition: xf_ProcessMenu_Contours
static void
xf_ProcessMenu_Contours(int k)
{
  switch (k){
  case 0:
    xf_KeyCommon('c', 0, 0);
    break;
  case 1:
    xf_KeyCommon('C', 0, 0);
    break;
  default:
    xf_printf("Unknown menu item\n");
    break;
  }
  nModifyMenu = 1;
  ModifyMenuList[0] = xfe_PlotMenu_Sub_Contours;
  xf_PostRedisplayAll();
}


/******************************************************************/
//   FUNCTION Definition: xf_ProcessMenu_IsoSurface
static void
xf_ProcessMenu_IsoSurface(int k)
{
  if (iDGroupIsoSurf >= 0){
    switch (k){
    case 0:
      DGroup[iDGroupIsoSurf].Active = !DGroup[iDGroupIsoSurf].Active;
      if ((DGroup[iDGroupIsoSurf].Active) && (DGroup[iDGroupIsoSurf].DTri == NULL)){
	xf_printf("No iso-surface exists; create one by setting a value.\n");
	DGroup[iDGroupIsoSurf].Active = xfe_False; // does not exist, so don't turn on
      }
      xf_SetActiveDGroups();
      break;
    case 1:
      DGroup[iDGroupIsoSurf].Active = xfe_False;
      xf_SpecialKey3D(GLUT_KEY_F4, 0, 0);
      break;
    default:
      xf_printf("Unknown menu item\n");
      break;
    }
    nModifyMenu = 1;
    ModifyMenuList[0] = xfe_PlotMenu_Sub_IsoSurface;
    xf_PostRedisplayAll();
  }
}

/******************************************************************/
//   FUNCTION Definition: xf_ProcessMenu_Lighting
static void
xf_ProcessMenu_Lighting(int k)
{
  int ierr;
  switch (k){
  case 0:
    LightingFlag = !LightingFlag;
    if (LightingFlag){
      ierr = xf_Error(xf_FillLighting());
      if (ierr != xf_OK) exit(ierr);
      xf_printf("Lighting has been turned on.\n");
    }
    else xf_printf("Lighting has been turned off.\n");
    break;
  case 1:
    LightingFlag = xfe_False;
    xf_Key3D('l', 0, 0);
    break;
  default:
    xf_printf("Unknown menu item\n");
    break;
  }
  nModifyMenu = 1;
  ModifyMenuList[0] = xfe_PlotMenu_Sub_Lighting;
  xf_PostRedisplayAll();
}

/******************************************************************/
//   FUNCTION Definition: xf_ProcessMenu_ViewFrom
static void
xf_ProcessMenu_ViewFrom(int k)
{
  xf_Key3D('x'+k, 0, 0);
}


/******************************************************************/
//   FUNCTION Definition: xf_ProcessMenu_Hardcopy
static void
xf_ProcessMenu_Hardcopy(int k)
{
  int ierr, dim;

  switch (k){
  case 0:
    ierr = xf_Error(xf_HardcopyWindowEPS((*WinBig), "plot.eps", xfe_True));
    if (ierr != xf_OK) exit(ierr);
    xf_printf("Wrote plot.eps of big window.\n");
    break;
  case 1:
    dim = (WinBig == &Win3D) ? 3 : 2;
    ierr = xf_MakeMovie(dim);
    if (ierr == xf_OK) xf_printf("Wrote animation of big window.\n");
    break;
  }

}

/******************************************************************/
//   FUNCTION Definition: xf_ProcessMenu_ColorLimits
static void
xf_ProcessMenu_ColorLimits(int k)
{
  int ierr, dim;

  switch (k){
  case 0:
    ierr = xf_Error(xf_FillDScalarsAll(DVector + CurrentDVector, xfe_False));
    if (ierr != xf_OK) exit(ierr);
    if (ContoursOn){
      ierr = xf_Error(xf_FillContours());
      if (ierr != xf_OK) exit(ierr);
    }
    break;
  case 1:
    ierr = xf_Error(xf_ChangeColorLimits());
    if (ierr != xf_OK) exit(ierr);
    if (ContoursOn){
      ierr = xf_Error(xf_FillContours());
      if (ierr != xf_OK) exit(ierr);
    }
    break;
  }

  xf_PostRedisplayAll();

}



/******************************************************************/
//   FUNCTION Definition: CreatePlotMenu
static int
xf_CreatePlotMenu(void)
{
  int k, orig;
  enum xfe_PlotMenuOptions opt;

  // make top-level menu
  TopMenu = glutCreateMenu(xf_ProcessMenu_Top); 
  
  // Modify/Attach sub-menus
  for (k=1; k<xfe_PlotMenu_Max; k++){
    SubMenuIndex[k] = 0;  // zero out submenu index
    glutAddMenuEntry("unused", k); // dummy menu entries
  }
  xf_ModifySubMenu_LoadVector();
  xf_ModifySubMenu_ChangeScalar();
  xf_ModifySubMenu_Mesh();
  xf_ModifySubMenu_Boundaries();
  xf_ModifySubMenu_ColorScheme();
  xf_ModifySubMenu_CutPlane();
  xf_ModifySubMenu_LineProbe();
  xf_ModifySubMenu_Contours();
  xf_ModifySubMenu_IsoSurface();
  xf_ModifySubMenu_Lighting();
  xf_ModifySubMenu_ViewFrom();
  xf_ModifySubMenu_Hardcopy();
  xf_ModifySubMenu_ColorLimits();

  // Add entries
  glutSetMenu(TopMenu);
  opt = xfe_PlotMenu_DeriveVector; glutChangeToMenuEntry(opt, "Derive a vector [i]", opt);
  opt = xfe_PlotMenu_SwapWindows;  glutChangeToMenuEntry(opt, "Swap big/small windows", opt);
  opt = xfe_PlotMenu_ToggleRender; glutChangeToMenuEntry(opt, "Toggle rendering", opt);
  opt = xfe_PlotMenu_Help;         glutChangeToMenuEntry(opt, "Help", opt);
  opt = xfe_PlotMenu_Quit;         glutChangeToMenuEntry(opt, "Quit", opt);
  
  // Attach same menu to all windows
  orig = glutGetWindow();
  glutSetWindow(Win2D.Handle);
  glutAttachMenu(GLUT_LEFT_BUTTON);
  glutSetWindow(Win3D.Handle);
  glutAttachMenu(GLUT_LEFT_BUTTON);
  glutSetWindow(WinAux.Handle);
  glutAttachMenu(GLUT_LEFT_BUTTON);
  glutSetWindow(WinBorder.Handle);
  glutAttachMenu(GLUT_LEFT_BUTTON);
  glutSetWindow(orig);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Visualize
static int 
xf_Visualize(int argc, char *argv[])
{
  int ierr, k;
  xf_Data *StateData;

  // Initialize windows
  xf_InitWindows();

  // calculate default refinement level for each element
  ierr = xf_Error(xf_InitElem2Refine());
  if (ierr != xf_OK) return ierr;

  // calculate model size
  ierr = xf_Error(xf_CalculateModelSize());
  if (ierr != xf_OK) return ierr;

  // initialize display groups
  ierr = xf_Error(xf_InitDGroups());
  if (ierr != xf_OK) return ierr;

  // calculate active display groups
  ierr = xf_Error(xf_CalculateActiveDGroups());
  if (ierr != xf_OK) return ierr;

  // Set 3d camera default
  xf_InitCamera();
  
  // set contours default;
  xf_InitContours();


  // Set suggestion for parent and calculate children window sizes
  ierr = xf_Error(xf_CalculateWindowSizes());
  if (ierr != xf_OK) return ierr;

  // Set default mouse state
  MouseState.LeftButtonDown = xfe_False;
  MouseState.MiddleButtonDown = xfe_False;
  MouseState.RightButtonDown = xfe_False;
  MouseState.Modifier = 0;
  MouseState.x = MouseState.y = 0;

  // set default mesh state
  MeshOn    = xfe_False;
  SubMeshOn = xfe_False;

  // load vectors available for display
  ierr = xf_Error(xf_LocateVectors(All));
  if (ierr != xf_OK) return ierr;

  // try to derive a default set of eqnset scalars
  if (HaveEqnSet){
    for (k=0; k<nDVector; k++){
      if (strcmp(DVector[k].DataName, "State") == 0){
	ierr = xf_Error(xf_LoadDVector_Scalars(DVector[k].Vector, 
					       "State_Scalars", NULL));
	if (ierr != xf_OK) return ierr;
	break;
      }
    } // k
  }

  // Calculate colors for Color Bar
  xf_FillColorBar();

  // Read .state file if available
  ierr = xf_ReadStateFile();
  if (ierr == xf_OK) xf_printf("Read plotter state file: %s\n", StateFile); 

  // re-calculate active display groups
  ierr = xf_Error(xf_CalculateActiveDGroups());
  if (ierr != xf_OK) return ierr;

  //////for testing
  //
  
  if (nDVector <= 0){
    xf_printf("No Vectors loaded for rendering.  Use 'v' to load a vector.\n");
    fflush(stdout);
  }
  else if ((CurrentDVector < 0) || (CurrentDVector >= nDVector)){
    xf_printf("Invalid CurrentDVector = %d.\n", CurrentDVector);
    fflush(stdout);
  }
  else if (DVector[CurrentDVector].nScalar <= 0){
    xf_printf("No scalars present in current vector.\n");
    fflush(stdout);
  }
  else if ((DVector[CurrentDVector].CurrentDScalar < 0) ||
	   (DVector[CurrentDVector].CurrentDScalar >= DVector[CurrentDVector].nScalar)){
    xf_printf("Invalid CurrentDScalar = %d.\n", DVector[CurrentDVector].CurrentDScalar);
    fflush(stdout);
  }
  else{
     xf_printf("Render is turned on.\n");
     ierr = xf_Error(xf_FillDScalarsAll(DVector + CurrentDVector, xfe_True));
     if (ierr != xf_OK) exit(ierr);
  }
  ierr = xf_Error(Yu_DumpAllDataToTecplot(DVector + CurrentDVector));
  if (ierr != xf_OK) exit(ierr);
  xf_printf("Tecplot-formatted data dump successfully!\n");
  return xf_OK;
  
  ///////////////////

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  
  // Parent Window
  glutInitWindowSize(WinParent.PixWidth[0], WinParent.PixWidth[1]);
  glutInitWindowPosition(WinParent.PixTopLeft[0], WinParent.PixTopLeft[1]);
  WinParent.Handle = glutCreateWindow("XFlow Plotter");

  glutDisplayFunc(xf_Idle);

  //glutKeyboardFunc(xf_KeyCommon);
  glutReshapeFunc(xf_Resize);


  // 2D Window
  Win2D.Handle = glutCreateSubWindow(WinParent.Handle, Win2D.PixTopLeft[0], 
				     Win2D.PixTopLeft[1], Win2D.PixWidth[0], 
				     Win2D.PixWidth[1]);
  glutKeyboardFunc(xf_Key2D);
  glutSpecialFunc(xf_SpecialKey2D);
  glutMouseFunc(xf_Mouse2D);
  glutMotionFunc(xf_Motion2D);
  glutDisplayFunc(xf_Display2D);
  glutPassiveMotionFunc(xf_PassiveMotion);
  //glutMouseWheelFunc(xf_MouseWheel2D);
  //glutIdleFunc(xf_Idle);
  xf_InitGraphics();

  // 3D Window
  Win3D.Handle = glutCreateSubWindow(WinParent.Handle, Win3D.PixTopLeft[0], 
				     Win3D.PixTopLeft[1], Win3D.PixWidth[0], 
				     Win3D.PixWidth[1]);
  glutKeyboardFunc(xf_Key3D);
  glutSpecialFunc(xf_SpecialKey3D);
  glutMouseFunc(xf_Mouse3D);
  glutDisplayFunc(xf_Display3D);
  glutMotionFunc(xf_Motion3D);
  glutPassiveMotionFunc(xf_PassiveMotion);
  //glutMouseWheelFunc(xf_MouseWheel3D);
  //glutIdleFunc(xf_Idle);
  xf_InitGraphics();


  // Aux Window
  WinAux.Handle = glutCreateSubWindow(WinParent.Handle, WinAux.PixTopLeft[0], 
				      WinAux.PixTopLeft[1], WinAux.PixWidth[0], 
				      WinAux.PixWidth[1]);
  glutKeyboardFunc(xf_KeyAux);
  glutSpecialFunc(xf_SpecialKeyAux);
  glutDisplayFunc(xf_DisplayAux);
  glutPassiveMotionFunc(xf_PassiveMotion);
  //glutIdleFunc(xf_Idle);
  xf_InitGraphics();

  // Border Window
  WinBorder.Handle = glutCreateSubWindow(WinParent.Handle, WinBorder.PixTopLeft[0], 
				      WinBorder.PixTopLeft[1], WinBorder.PixWidth[0], 
				      WinBorder.PixWidth[1]);
  glutDisplayFunc(xf_DisplayBorder);
  glutPassiveMotionFunc(xf_PassiveMotion);
  //glutIdleFunc(xf_Idle);
  xf_InitGraphics();

 
  // Set Current Window
  glutSetWindow(Win2D.Handle);
  
  // Create menu
  ierr = xf_Error(xf_CreatePlotMenu());
  if (ierr != xf_OK) return ierr;

  // Turn rendering on if flag is true
  if (RenderFlag) xf_TurnRenderOn();
  
  // Calculate lighting if on
  if (LightingFlag){
    ierr = xf_Error(xf_FillLighting());
    if (ierr != xf_OK) exit(ierr);
  }

  // turn over execution to GLUT event handler
  glutMainLoop();
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadGCLFromDataFile
static int
xf_ReadGCLFromDataFile(char* DataFile)
{
  int ierr;
  xf_DataSet *DataSet;
  xf_Data *D;
  xf_GenArray *gA = NULL;
  xf_Vector *GCL = NULL, *VGCL = NULL;
   
  //Create a dataset
  ierr = xf_Error(xf_CreateDataSet(&DataSet));
  if (ierr != xf_OK) return ierr;
    
  ierr = xf_Error(xf_ReadDataSetBinary(NULL, NULL, DataFile, DataSet));
  if (ierr!=xf_OK) return ierr;

  //Find and use GCLVectorMain if it exists in .data file; otherwise, use maximum-i GCL_i vector
  ierr = xf_FindDataByTitle(DataSet, "GCLVectorMain", xfe_Vector, &D);
  if (ierr == xf_NOT_FOUND)
    ierr = xf_FindDataByTitle(DataSet, "GCL_2", xfe_Vector, &D);
  if (ierr == xf_NOT_FOUND)
    ierr = xf_FindDataByTitle(DataSet, "GCL_1", xfe_Vector, &D);
  if (ierr == xf_NOT_FOUND)
      ierr = xf_FindDataByTitle(DataSet, "GCL_0", xfe_Vector, &D);
  if (ierr == xf_NOT_FOUND)
    xf_printf("No GCL vectors found in data file -- plotted data may be unreliable.\n");
     
  //If GCL is found in .data file, replace GCL in All with this new GCL data
  if (ierr == xf_OK){
    xf_printf("\n Note: Using %s as the default GCL vector.\n \n", D->Title);
    VGCL = (xf_Vector *) D->Data;
      
    //Find GCL in the All struct
    ierr = xf_Error(xf_FindMeshMotionGCLVector(All, &GCL));
    if (ierr != xf_OK) return ierr;

    //Swap just-loaded VGCL contents into All struct GCL (by swapping pointers)
    swap(GCL->GenArray, VGCL->GenArray, gA);
    
    /*Destroy dataset. Note that if loaded GCL was of different size than GCL in All, 
    not all created memory is destroyed here (may want to change) */
    ierr = xf_Error(xf_DestroyDataSet(DataSet));
    if (ierr != xf_OK) return ierr;
  } 

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Convert1Dto2D
static int 
xf_Convert1Dto2D(xf_All *All){

  int ierr;
  int i, j, k, l;
  int sr, r1, r2, n2, nlevy;
  int egrp, nelem, elem;
  int nNode0;
  int *nFaceNew;
  enum xfe_Bool Interpolated;
  real yval, yratio = 0.1; // ratio for extent in y-direction
  int *vr = NULL;
  real **rValue, **rtemp;
  xf_Face NullFace = {xf_NULLFACE, 0};
  xf_ElemGroup *EG;
  xf_DataSet *DataSet;
  xf_Data *D;
  xf_Vector *V;
  xf_GenArray *ga;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  if (All->Mesh->Dim != 1 ) return xf_Error(xf_INPUT_ERROR);

  // calculate model size
  ierr = xf_Error(xf_CalculateModelSize());
  if (ierr != xf_OK) return ierr;
  yval = yratio*(ModelBBox[1] - ModelBBox[0]);

  xf_printf("yval = %.10E\n", yval);

  // Dim
  Mesh->Dim = 2;

  // Coord and nNode
  ierr = xf_Error(xf_ReAllocCopy2((void ***) &Mesh->Coord, Mesh->nNode, 1, 
				  2*Mesh->nNode, 2, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  nNode0 = Mesh->nNode;
  for (i=0; i<nNode0; i++){
    Mesh->Coord[Mesh->nNode+i][0] = Mesh->Coord[i][0];
    Mesh->Coord[i][1] = 0.;
    Mesh->Coord[Mesh->nNode+i][1] = yval;
  }
  Mesh->nNode *= 2;

  // IFace -- not changing yet
  // BFaceGroup -- not changing yet
  
  // ElemGroup
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    EG = Mesh->ElemGroup + egrp;
    if ((EG->CutFlag) || (EG->QBasis != xfe_SegLagrange) || (EG->QOrder != 1) ||
	(EG->nNode != 2))
      return xf_Error(xf_NOT_SUPPORTED);
    
    nelem = EG->nElem;
    EG->QBasis = xfe_QuadLagrange;
    EG->nNode = 4;

    // nFace
    ierr = xf_Error(xf_Alloc((void **) &nFaceNew, nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (elem=0; elem<nelem; elem++) nFaceNew[elem] = 4;

    // Face
    ierr = xf_Error(xf_VReAllocCopy2((void ***) &Mesh->ElemGroup[egrp].Face, nelem, 
				     Mesh->ElemGroup[egrp].nFace, nelem, nFaceNew, 
				     sizeof(xf_Face)));
    if (ierr!=xf_OK) return ierr;
    xf_Release( (void *) Mesh->ElemGroup[egrp].nFace);
    Mesh->ElemGroup[egrp].nFace = nFaceNew;
    for (elem=0; elem<nelem; elem++){
      Mesh->ElemGroup[egrp].Face[elem][3] = Mesh->ElemGroup[egrp].Face[elem][0];
      Mesh->ElemGroup[egrp].Face[elem][0] = NullFace;
      Mesh->ElemGroup[egrp].Face[elem][2] = NullFace;
    }

    // Node
    ierr = xf_Error(xf_ReAllocCopy2((void ***) &Mesh->ElemGroup[egrp].Node, nelem, 2, 
				    nelem, 4, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
    for (elem=0; elem<nelem; elem++){
      Mesh->ElemGroup[egrp].Node[elem][2] = Mesh->ElemGroup[egrp].Node[elem][0]+nNode0;
      Mesh->ElemGroup[egrp].Node[elem][3] = Mesh->ElemGroup[egrp].Node[elem][1]+nNode0;
    }
  }

  // Data
  DataSet = All->DataSet;
  D = DataSet->Head;
  while (D != NULL){
    if (D->Type == xfe_Vector){
      V = (xf_Vector *) D->Data;
      Interpolated = ((V->Basis != NULL) && (V->Order != NULL));
      if ((V->Linkage == xfe_LinkageGlobElem) && (Interpolated)){
	sr = V->StateRank;

	// loop over each array
	for (i=0; i<V->nArray; i++){
	  
	  // change basis to Quad
	  if (V->Basis[i] == xfe_SegLagrange)
	    V->Basis[i] = xfe_QuadLagrange;
	  else if (V->Basis[i] == xfe_SegLagrangeGauss)
	    V->Basis[i] = xfe_QuadLagrangeGauss;

	  // current number of unknowns
	  r1 = V->GenArray[i].r;

	  // desired number of unknowns
	  ierr = xf_Error(xf_Order2nNode(V->Basis[i], V->Order[i], &n2));
	  if (ierr != xf_OK) return ierr;
	  r2 = n2*sr;
	  
	  if ((r2 % r1) != 0) return xf_Error(xf_INPUT_ERROR);

	  if (V->GenArray[i].vr == NULL){
	    // allocate an rValue with r2
	    ierr = xf_Error(xf_Alloc2((void ***) &rValue, V->GenArray[i].n, r2, sizeof(real)));
	    if (ierr != xf_OK) return ierr;
	    vr = NULL;
	  }
	  else{
	    // variable order
	    ierr = xf_Error(xf_Alloc((void **) &vr, V->GenArray[i].n, sizeof(int)));
	    if (ierr != xf_OK) return ierr;

	    for (j=0; j<V->GenArray[i].n; j++){
	      r1 = V->GenArray[i].r;
	      
	      // desired number of unknowns
	      ierr = xf_Error(xf_Order2nNode(V->Basis[i], V->vOrder[i][j], &n2));
	      if (ierr != xf_OK) return ierr;

	      vr[j] = n2*sr;
	    } // j

	    // allocate an rValue with vr
	    ierr = xf_Error(xf_VAlloc2((void ***) &rValue, V->GenArray[i].n, vr, sizeof(real)));
	    if (ierr != xf_OK) return ierr;

	  }
	  
	  // fill rValue
	  ga = V->GenArray+i;
	  for (j=0; j<ga->n; j++){
	    r1 = ((ga->vr==NULL) ? ga->r : ga->vr[j]);
	    r2 = ((ga->vr==NULL) ? r2    : vr[j]);
	    nlevy = r2/r1; // number of levels in y
	    for (l=0; l<nlevy; l++)
	      for (k=0; k<r1; k++)
		rValue[j][l*r1 + k] = ga->rValue[j][k];
	  }
	  // swap rValue with V->GenArray[i].rValue
	  swap(rValue, V->GenArray[i].rValue, rtemp);
    
	  // release rValue
	  xf_Release2( (void **) rValue);

	  // release vr
	  xf_Release( (void *) vr);

	  
	  // set V->GenArray[i].r = r2;
	  V->GenArray[i].r = r2;
	  
	} // i

      }
    }
    D = D->Next;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int len, ierr;
  int i, nData;
  char *ArgIn[] = {"xfa", "NULL", ".xfa file to read",
		   "data", "NULL", "alternate .data file(s) to read",
		   "eqn", "NULL", "alternate .eqn file to use",
		   "mm", "NULL", "alternate .mm (mesh motion) file",
		   "\0"};

  char InputFile[xf_MAXSTRLEN];
  char DataFile[xf_MAXSTRLEN];
  char EqnSetFile[xf_MAXSTRLEN];
  char MeshMotionFile[xf_MAXSTRLEN];
  char *pext, **DataFiles=NULL;
  xf_DataSet *DataSet = NULL;
  xf_Data *GammaDat;
  xf_EqnSet *EqnSet;
  xf_KeyValue KeyValueArg;
  enum xfe_Bool UseGCL = xfe_False;

  xf_printf("\n");
  xf_printf("=== xf_Plot ===\n");
  xf_printf("\n");
  xf_printf(" ** Type \"?\" in plotting window for menu of available options **\n");
  xf_printf("\n");

      
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValueArg));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValueArg);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);

  /* Get InputFile */
  ierr = xf_GetKeyValue(KeyValueArg, "xfa", InputFile);
  if (ierr != xf_OK) return ierr;

  // handle special usage shortcut of just a .xfa file
  if (!xf_NotNull(InputFile)){
    if (argc == 2){
      strcpy(InputFile, argv[1]);
      len = strlen(InputFile);
      pext = InputFile + len - 4; // pointer to extension
      if ((len < 4) || ((strncmp(pext, ".xfa", 4) != 0) &&
			(strncmp(pext, ".gri", 4) != 0))){
	xf_printf("Error, InputFile requires .xfa or .gri extension.\n");
	return xf_Error(xf_INPUT_ERROR);
      }
    }
    else{
      xf_printf("Quick usage (xfa/gri file only):\n");
      xf_printf("  xf_Plot <xfafile>\n");
      xf_printf("or general usage (run with no args to see all options):\n");
      xf_printf("  xf_Plot -xfa <xfafile> -data <datafiles> ... \n");
      xf_printf("\n");
      return xf_OK;
    }
  }

  /* Create .xfa structure */
  ierr = xf_Error(xf_CreateAll(&All, xfe_False));
  if (ierr != xf_OK) return ierr;

  /* Read .xfa or .gri file */
  len = strlen(InputFile);
  if (len < 4) return xf_Error(xf_INPUT_ERROR);
  pext = InputFile + len - 4; // pointer to extension
  if (strncmp(pext, ".xfa", 4) == 0){
    ierr = xf_Error(xf_ReadAllBinary(InputFile, All));
    if (ierr!=xf_OK) return ierr;
  }
  else if (strncmp(pext, ".gri", 4) == 0){
    ierr = xf_Error(xf_ReadGriFile(InputFile, NULL, All->Mesh));
    if (ierr!=xf_OK) return ierr;

    // set default parameters upon reading .gri
    ierr = xf_AddKeyValueList(&All->Param->KeyValue, xf_DefaultParamList,
                              xfe_True, xfe_True);
    if ((ierr != xf_OK) && (ierr != xf_NOT_FOUND)) return xf_Error(ierr);

  }


  /* Get DataFile */
  ierr = xf_GetKeyValue(KeyValueArg, "data", DataFile);
  if (ierr != xf_OK) return ierr;

  /* Read DataFile if specified */
  if (xf_NotNull(DataFile)){
    
    xf_printf("Reading alternate DataSet file(s): %s\n", DataFile);

    ierr = xf_Error(xf_ScanXStringAlloc(DataFile, xf_MAXSTRLEN, &nData, &DataFiles));
    if (ierr != xf_OK) return ierr;
    
    for (i=0; i<nData; i++){
      
      ierr = xf_Error(xf_CreateDataSet(&DataSet));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ReadDataSetBinary(NULL, NULL, DataFiles[i], DataSet));
      if (ierr!=xf_OK) return ierr;
      
      ierr = xf_Error(xf_DataSetMerge(DataSet, All->DataSet));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_DestroyDataSet(DataSet));
      if (ierr != xf_OK) return ierr;

    }

    xf_Release2( (void **) DataFiles);

  }

  //Yu: have data already; let's do gamma vector
 /* Yu_Model Model;
  ierr = xf_Error(PullinModel(&Model));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(Yu_GammaVectorCreate(All, &Model, NULL));
  if(ierr != xf_OK) return ierr;
*/
  
  
  /* Convert 1D mesh+data to 2D for purpose of plotting */
  if (All->Mesh->Dim == 1){
    ierr = xf_Error(xf_Convert1Dto2D(All));
    if (ierr != xf_OK) return ierr;
  }

  /* Get EqnSetFile */
  ierr = xf_GetKeyValue(KeyValueArg, "eqn", EqnSetFile);
  if (ierr != xf_OK) return ierr;
  
  /* Read EqnSetFile if specified */
  if (xf_NotNull(EqnSetFile)){

    xf_printf("Reading alternate EqnSet file: %s\n", EqnSetFile);

    ierr = xf_Error(xf_CreateEqnSet(&EqnSet));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ReadEqnSetFile(EqnSetFile, NULL, EqnSet));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_DestroyEqnSet(All->EqnSet, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    All->EqnSet = EqnSet;
    All->EqnSet->Dim = All->Mesh->Dim;
    
    // set default params in All
    ierr = xf_Error(xf_AddKeyValueList(&All->Param->KeyValue, xf_DefaultParamList,
				       xfe_True, xfe_True));
    if (ierr != xf_OK) return ierr;


  }

  /* Get MeshMotionFile */
  ierr = xf_GetKeyValue(KeyValueArg, "mm", MeshMotionFile);
  if (ierr != xf_OK) return ierr;
  
  /* Read MeshMotion if specified */
  if (xf_NotNull(MeshMotionFile)){

    xf_printf("Reading mesh motion file: %s\n", MeshMotionFile);

    // first destroy any current mesh motion
    ierr = xf_Error(xf_DestroyMeshMotion(All->Mesh->Motion));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_CreateMeshMotion(&All->Mesh->Motion));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ReadMeshMotionFile(MeshMotionFile, NULL, All->Mesh->Motion));
    if (ierr != xf_OK) return ierr;

  }
  
  // Determine if using a Geometric Conservation Law
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  if (ierr == xf_NOT_FOUND){
    xf_printf("\"UseGCL\" not found in the parameter list. Maybe old .xfa\n");
  }
  else if (ierr != xf_OK) return ierr;
    
  // If GCL on, fill values from .data file
  if((UseGCL) && xf_NotNull(DataFile)){
    ierr = xf_Error(xf_ReadGCLFromDataFile(DataFile));
    if (ierr != xf_OK) return ierr;
  }

  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValueArg));
  if (ierr!=xf_OK) return ierr;

  /* Visualize solution */
  ierr = xf_Error(xf_Visualize(argc, argv));
  if (ierr != xf_OK) return ierr;

  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;

  xf_printf("xf_Plot finished.\n");

  return xf_OK;
}


