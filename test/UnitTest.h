/*******************************************************************************                              
 **                                                                                                           
 ** "UnitTest.h": Headers and wrappers for unit testing
 ** Author: Vedaad Shakib
 **                                                        
 *******************************************************************************/

#ifndef _UnitTest_h
#define _UnitTest_h 1

#define UTOL0 1e-15
#define UTOL1 1e-14
#define UTOL2 1e-13
#define UTOL3 1e-12
#define UTOL4 1e-11
#define UTOL5 1e-10

using namespace std;

::testing::AssertionResult AssertRealVectorWithin(double *expected, double *Actual, int size, double UTOL);
::testing::AssertionResult AssertIntVectorEqual(int *expected, int *Actual, int size);
::testing::AssertionResult AssertEqual(long expected, long actual);
::testing::AssertionResult AssertWithin(double expected, double actual, double UTOL);

::testing::AssertionResult AssertEqual(long expected,
				       long actual) {
    if (expected != actual) {
	return ::testing::AssertionFailure() << "expected value: " << expected
	                                     << "actual: " << actual;
    }

    return ::testing::AssertionSuccess();
}

::testing::AssertionResult AssertWithin(double expected,
					double actual,
					double UTOL) {
    if (fabs(expected-actual) > UTOL) {
	return ::testing::AssertionFailure() << "expected value: " << expected
					     << "; actual: " << actual;
    }

    return ::testing::AssertionSuccess();
}

::testing::AssertionResult AssertIntVectorEqual(int *expected,
						int *actual,
                                                int size) {
    for (int i = 0; i < size; i++){
	if (expected[i] != actual[i]) {
	    return ::testing::AssertionFailure() << "array[" << i
						 << "] (" << actual[i] << ") != expected[" << i
						 << "] (" << expected[i] << ")";
	}
    }

    return ::testing::AssertionSuccess();
}

::testing::AssertionResult AssertRealVectorWithin(double *expected,
						  double *actual,
						  int size,
						  double UTOL) {
    for (int i = 0; i < size; i++){
	if (fabs(expected[i]-actual[i]) > UTOL) {
	    return ::testing::AssertionFailure() << "array[" << i
						 << "] (" << actual[i] << ") != expected[" << i
						 << "] (" << expected[i] << ")";
	}
    }
    
    return ::testing::AssertionSuccess();
}



#endif // end ifndef _xf_Unit_h

