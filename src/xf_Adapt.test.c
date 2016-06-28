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


#include "xf_Unit.h"
#include "xf_All.h"
#include "xf_Data.h"
#include "xf_SolverTools.h"
#include <stdlib.h>

// Functions for setting up and running a case within a unit test
#include "xf_UnitRun.c"


// auxiliary function for some tests below
static int 
xf_MakeDefaultTimeHistData( real StartTime, real EndTime, int nTimeStep,
                            xf_TimeHistData **pTimeHistData)
{
  int ierr, nTime, i;
  enum xfe_TimeSchemeType TimeScheme = xfe_TimeSchemeDG1;
  real TimeStep;
  xf_TimeHistData *TimeHistData;
  
  // allocate structure for TimeHistData
  ierr = xf_Error(xf_CreateTimeHistData(pTimeHistData));
  if (ierr != xf_OK) return ierr;
  TimeHistData = (*pTimeHistData);
   
  // uniform spacing means constant time step
  TimeHistData->ConstTimeStep = xfe_True;
  
  // Calculate time step
  TimeStep = (EndTime - StartTime)/nTimeStep;
    
  // For FE in time, we count time slabs
  nTime = nTimeStep;

  // set nTime, allocate Time, TimeStep, TimeScheme
  TimeHistData->nTime = nTime;
  ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->Time, nTime, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->TimeStep, nTime, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->TimeScheme, nTime,
                             sizeof(enum xfe_TimeSchemeType)));
  if (ierr != xf_OK) return ierr;
  
  // fill Time, TimeStep, TimeScheme
  for (i=0; i<nTime; i++) TimeHistData->Time[i]       = StartTime + i*TimeStep;
  for (i=0; i<nTime; i++) TimeHistData->TimeStep[i]   = TimeStep;
  for (i=0; i<nTime; i++) TimeHistData->TimeScheme[i] = TimeScheme;
  
  return xf_OK;
}


TEST_xf_AdaptTimeHistData()
{
  /*
    Tests adaptation of time history
  */  
  int ierr, i;
  int RefIndicatorTime[10] = {0};
  xf_TimeHistData *TimeHistData;
  xf_TimeHistData *OldTimeHistData;

  // get a default time history
  ierr = xf_Error(xf_MakeDefaultTimeHistData(1.0, 3.0, 4, &TimeHistData));
  if (ierr != xf_OK) return ierr;

  // uniform refinement
  xf_printf("\n*** Refinement ***\n");
  for (i=0; i<4; i++) RefIndicatorTime[i] = 1;
  ierr = xf_Error(xf_AdaptTimeHistData(TimeHistData, RefIndicatorTime, &OldTimeHistData));
  if (ierr != xf_OK) return ierr;
  xf_AssertEqual(TimeHistData->nTime, 8);
  xf_AssertWithin(TimeHistData->Time[1], 1.25, UTOL0);
  xf_AssertEqual(OldTimeHistData->nTime, 4);

  // destroy old time history
  ierr = xf_Error(xf_DestroyTimeHistData(OldTimeHistData));
  if (ierr != xf_OK) return ierr;

  // uniform coarsening
  xf_printf("\n*** Coarsening ***\n");
  for (i=0; i<8; i++) RefIndicatorTime[i] = -1;
  ierr = xf_Error(xf_AdaptTimeHistData(TimeHistData, RefIndicatorTime, &OldTimeHistData));
  if (ierr != xf_OK) return ierr;
  xf_AssertEqual(TimeHistData->nTime, 4);
  xf_AssertWithin(TimeHistData->Time[1], 1.5, UTOL0);
  xf_AssertEqual(OldTimeHistData->nTime, 8);


  // one slab coarsening
  xf_printf("\n*** One slab coarsening ***\n");
  for (i=0; i<4; i++) RefIndicatorTime[i] = 0; 
  RefIndicatorTime[0] = -1;
  ierr = xf_Error(xf_AdaptTimeHistData(TimeHistData, RefIndicatorTime, &OldTimeHistData));
  if (ierr != xf_OK) return ierr;
  xf_AssertEqual(TimeHistData->nTime, 4); // rounding takes place (3.5->4)
  xf_AssertWithin(TimeHistData->Time[0], 1.0, UTOL0);
  xf_AssertWithin(TimeHistData->Time[1], 1.5+0.5*3/8, UTOL0);
  xf_AssertWithin(TimeHistData->Time[2], 2.0+0.5*2/8, UTOL0);
  xf_AssertWithin(TimeHistData->Time[3], 2.5+0.5*1/8, UTOL0);
  xf_AssertWithin(TimeHistData->TimeStep[0], 0.5+0.5*3/8, UTOL0);
  xf_AssertWithin(TimeHistData->TimeStep[1], 0.5*7/8, UTOL0);
  xf_AssertWithin(TimeHistData->TimeStep[2], 0.5*7/8, UTOL0);
  xf_AssertWithin(TimeHistData->TimeStep[3], 0.5*7/8, UTOL0);
  xf_AssertEqual(OldTimeHistData->nTime, 4);
  xf_AssertWithin(OldTimeHistData->Time[0], 1.0, UTOL0);
  xf_AssertWithin(OldTimeHistData->Time[2], 2.0, UTOL0);


  // destroy time histories
  ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyTimeHistData(OldTimeHistData));
  if (ierr != xf_OK) return ierr;

   // get a default time history
  ierr = xf_Error(xf_MakeDefaultTimeHistData(1.0, 3.0, 4, &TimeHistData));
  if (ierr != xf_OK) return ierr;

  // one coarsening and one refining
  xf_printf("\n*** One slab coarsening, one refining ***\n");
  for (i=0; i<4; i++) RefIndicatorTime[i] = 0; 
  RefIndicatorTime[1] = -1;
  RefIndicatorTime[2] =  1;
  ierr = xf_Error(xf_AdaptTimeHistData(TimeHistData, RefIndicatorTime, &OldTimeHistData));
  if (ierr != xf_OK) return ierr;
  xf_AssertEqual(TimeHistData->nTime, 5); // rounding takes place (4.5->5)
  xf_AssertWithin(TimeHistData->Time[0], 1.0, UTOL0);
  xf_AssertWithin(TimeHistData->Time[1], 1.0+0.5*9/10, UTOL0);
  xf_AssertWithin(TimeHistData->Time[2], 2.0+0.5*3/20, UTOL0);
  xf_AssertWithin(TimeHistData->Time[3], 2.0+0.5*12/20, UTOL0);
  xf_AssertWithin(TimeHistData->Time[4], 2.5+0.5*1/10, UTOL0);
  xf_AssertWithin(TimeHistData->TimeStep[0], 0.5*9/10, UTOL0);
  xf_AssertWithin(TimeHistData->TimeStep[1], 0.5*1/10+0.5+0.5*3/20, UTOL0);
  xf_AssertWithin(TimeHistData->TimeStep[2], 0.5*9/20, UTOL0);
  xf_AssertWithin(TimeHistData->TimeStep[3], 0.5*8/20+0.5*1/10, UTOL0);
  xf_AssertWithin(TimeHistData->TimeStep[4], 0.5*9/10, UTOL0);
  xf_AssertEqual(OldTimeHistData->nTime, 4);
  xf_AssertWithin(OldTimeHistData->Time[0], 1.0, UTOL0);
  xf_AssertWithin(OldTimeHistData->Time[2], 2.0, UTOL0);


  // destroy time histories
  ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyTimeHistData(OldTimeHistData));
  if (ierr != xf_OK) return ierr;

   // get a default time history
  ierr = xf_Error(xf_MakeDefaultTimeHistData(1.0, 3.0, 4, &TimeHistData));
  if (ierr != xf_OK) return ierr;

  // two coarsening, adjacent
  xf_printf("\n*** Two coarsening, adjacent ***\n");
  for (i=0; i<4; i++) RefIndicatorTime[i] = 0; 
  RefIndicatorTime[1] = -1;
  RefIndicatorTime[2] = -1;
  ierr = xf_Error(xf_AdaptTimeHistData(TimeHistData, RefIndicatorTime, &OldTimeHistData));
  if (ierr != xf_OK) return ierr;
  xf_AssertEqual(TimeHistData->nTime, 3); // rounding takes place (4.5->5)
  xf_AssertWithin(TimeHistData->Time[0], 1.0, UTOL0);
  xf_AssertWithin(TimeHistData->Time[1], 1.5, UTOL0);
  xf_AssertWithin(TimeHistData->Time[2], 2.5, UTOL0);
  xf_AssertWithin(TimeHistData->TimeStep[0], 0.5, UTOL0);
  xf_AssertWithin(TimeHistData->TimeStep[1], 1.0, UTOL0);
  xf_AssertWithin(TimeHistData->TimeStep[2], 0.5, UTOL0);


  // destroy time histories
  ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyTimeHistData(OldTimeHistData));
  if (ierr != xf_OK) return ierr;


  return xf_OK;  
}



TEST_xf_AdaptTimeHistDataRand()
{
  /*
    Tests adaptation of a pseudo-random time history
  */  
  int ierr, i, j;
  int N;
  int RefIndicatorTime[100] = {0};
  real StartTime, EndTime, r;
  real *Time, *TimeStep;
  xf_TimeHistData *TimeHistData;
  xf_TimeHistData *OldTimeHistData;

  StartTime = 1.0;
  EndTime   = 20.0;
  N         = 100;

  // get a default time history
  ierr = xf_Error(xf_MakeDefaultTimeHistData(StartTime, EndTime, N, &TimeHistData));
  if (ierr != xf_OK) return ierr;

  // initialize pseudo-random seed
  srand(17);

  // make time slabs pseduo-random
  xf_printf("\n*** Pseudo-random test ***\n");
  Time     = TimeHistData->Time;
  TimeStep = TimeHistData->TimeStep;
  for (i=1; i<N; i++){
    r = .01 + .09*((real) rand())/((real) RAND_MAX); // in [0.01,.1]
    Time[i] = Time[i-1] + (EndTime-Time[i-1])*r;
    //xf_printf("i=%d, Time = %.10E\n", i, Time[i]);
    TimeStep[i-1] = Time[i] - Time[i-1];
  }
  TimeStep[N-1] = EndTime - Time[N-1];

  // make a pseduo-random refinement indicator
  for (i=0; i<N; i++){
    j = ((19*i + 3*i*i + (1717%(i+1))) % 17);
    if (j<5) RefIndicatorTime[i] = -1;
    else if (j<13) RefIndicatorTime[i] = 0;
    else RefIndicatorTime[i] = 1;
    //xf_printf("RefIndicator[%d] = %d\n", i, RefIndicatorTime[i]);
  }
  
  // try adapting the time history
  ierr = xf_Error(xf_AdaptTimeHistData(TimeHistData, RefIndicatorTime, &OldTimeHistData));
  if (ierr != xf_OK) return ierr;

  // destroy time histories
  ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyTimeHistData(OldTimeHistData));
  if (ierr != xf_OK) return ierr;


  return xf_OK;  
}





TEST_xf_SplitTimeHistData()
{
  /*
    Tests splitting of time history
  */  
  int ierr, i;
  int RefIndicatorTime[10] = {0};
  xf_TimeHistData *TimeHistData;

  // get a default time history
  ierr = xf_Error(xf_MakeDefaultTimeHistData(1.0, 3.5, 5, &TimeHistData));
  if (ierr != xf_OK) return ierr;

  // split slab 1
  ierr = xf_Error(xf_SplitTimeHistData(TimeHistData, NULL, 1));
  if (ierr != xf_OK) return ierr;
  xf_AssertEqual(TimeHistData->nTime, 6);
  xf_AssertWithin(TimeHistData->Time[2], 1.75, UTOL0);
  xf_AssertWithin(TimeHistData->TimeStep[2], 0.25, UTOL0);

  // use a refinement indicator
  for (i=0; i<6; i++) RefIndicatorTime[i] = 0;
  RefIndicatorTime[2] = 1;
  RefIndicatorTime[3] = 1;
  ierr = xf_Error(xf_SplitTimeHistData(TimeHistData, RefIndicatorTime, 0));
  if (ierr != xf_OK) return ierr;
  xf_AssertEqual(TimeHistData->nTime, 8);
  xf_AssertWithin(TimeHistData->Time[3], 1.875, UTOL0);
  xf_AssertWithin(TimeHistData->Time[5], 2.25, UTOL0);

  // destroy time history
  ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
  if (ierr != xf_OK) return ierr;


  return xf_OK;  
}
