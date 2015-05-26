/* $Id: randomF77.c 2861 2002-11-22 13:24:41Z hothorn $
*
*  wrapper for calling R's random number generator from
*  the original FORTRAN code
*
*/

#include <R.h>

void F77_SUB(rndstart)(void) { GetRNGstate(); }
void F77_SUB(rndend)(void) { PutRNGstate(); }
double F77_SUB(unifrnd)(void) { return unif_rand(); }
