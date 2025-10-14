#include <string.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

/* 
 Calculates the moving-average of a matrix, column-wise.
 */

void mexFunction ( int nlhs, mxArray * plhs [], int rrhs, const mxArray * prhs [] ) {
  double *data, *cdata;
  double span;
  
  double *sdata, *csdata;
  double sdatum;
  
  int width, hwidth;
  int nsamp, ncol;
  int sindex, cindex;
  
  
  /* Checks the inputs. */
  if ( rrhs < 2 )
      mexErrMsgTxt ( "Not enough input arguments." );
  if ( rrhs > 2 )
      mexErrMsgTxt ( "Too many input arguments." );
  
  /* Gets the dimension of the problem. */
  nsamp = mxGetM ( prhs [0] );
  ncol  = mxGetN ( prhs [0] );
  
  if ( !mxIsScalar ( prhs [1] ) )
      mexErrMsgTxt ( "Second argument must be a scalar." );
  
  
  /* Sets the output. */
  plhs [0] = mxCreateDoubleMatrix ( nsamp, ncol, mxREAL );
  
  /* Gets the pointers to the real data. */
  sdata = mxGetData ( plhs [0] );
  data  = mxGetData ( prhs [0] );
  span  = mxGetScalar ( prhs [1] );
  
  if ( span > nsamp )
      mexErrMsgTxt ( "The span to average cannot be greater than the data." );
  
  
  /* Forces the averaging width to be simmetric. */
  width  = (int) span - 1 + (int) span % 2;
  hwidth = ( width - 1 ) / 2;
  
  
  /* Iterates through columns. */
  for ( cindex = 0; cindex < ncol; cindex ++ ) {
    
    /* Gets the pointer to the current column. */
    cdata  = data  + cindex * nsamp;
    csdata = sdata + cindex * nsamp;
    
    /* Initializes the buffer. */
    sdatum = 0;
    
    /* Loads the buffer. */
    for ( sindex = 0; sindex < width; sindex ++ ) {
      
      /* Stores the current value in the buffer. */
      sdatum += cdata [ sindex ];
      
      /* Uses the odd parts of the buffer for the rise. */
      if ( sindex % 2 == 0 ) {
        
        /* Stores the rising value. */
        csdata [ sindex / 2 ] = sdatum / ( sindex + 1 );
      }
    }
    
    /* Goes through each valid sample. */
    for ( sindex = 0; sindex < nsamp - width; sindex ++ ) {
      
      /* Stores the current value of the buffer */
      csdata [ sindex + hwidth ] = sdatum / width;
      
      /* Updates the buffer */
      sdatum += ( cdata [ sindex + width ] - cdata [ sindex ] );
    }
    
    /* Initializes the buffer. */
    sdatum = 0;
    
    /* Calculates the unrise. */
    for ( sindex = 0; sindex < width; sindex ++ ) {
      
      /* Stores the current value in the buffer. */
      sdatum += cdata [ nsamp - sindex - 1 ];
      
      /* Uses the odd parts of the buffer for the unrise. */
      if ( sindex % 2 == 0 ) {
        
        /* Stores the rising value. */
        csdata [ nsamp - sindex / 2 - 1 ] = sdatum / ( sindex + 1 );
      }
    }
  }
}
