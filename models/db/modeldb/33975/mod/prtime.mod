: $Header: /home/cvsroot/localstep/rings/prtime.mod,v 1.1.1.1 2004/07/04 13:05:55 hines Exp $
COMMENT
from Bill Lytton's Misc. routines:
prtime() // gives date/time
ENDCOMMENT
                           
NEURON {
    SUFFIX nothing
}

VERBATIM
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include <stdio.h>
ENDVERBATIM

:* PROCEDURE prtime ()
FUNCTION prtime () {
VERBATIM
  double prt;
  static double PRTIME;
  prt = (clock()-PRTIME)/CLOCKS_PER_SEC;
  // UINT_MAX for 32 bit machine -- see 'man clock'
  if (prt<0) prt += UINT_MAX/CLOCKS_PER_SEC; 
  PRTIME=clock();
  _lprtime = prt;
ENDVERBATIM
}
