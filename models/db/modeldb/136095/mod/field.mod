: $Id: field.mod,v 1.5 2009/12/24 02:25:27 samn Exp $

VERBATIM
#include "misc.h"
#define NSW 5  // store into this many ww field vecs
static void*   ww[NSW];
static double* wwo[NSW];
static int ixp1;
ENDVERBATIM

NEURON {
  POINT_PROCESS FIELD
  POINTER p1
  RANGE fldID
}

PARAMETER {
}

ASSIGNED {
  nsw
  wwsz
  p1 
  fldID
}

CONSTRUCTOR {
VERBATIM
{
  fldID = ifarg(2)?(int)*getarg(2):1; // assign a field ID so can record > 1 field
  if (fldID>=NSW) { printf("FIELD CONSTRUCTOR WARN: can only store %d ww vecs\n",NSW); hxe();}
}
ENDVERBATIM
}

INITIAL {
  VERBATIM
  ixp1=0;  
  ENDVERBATIM
}

:** initwrec(veclist) sets up recording of sim field potential
PROCEDURE initwrec () {
  VERBATIM 
  {int i, k, num, cap;  Object* ob;
    ob =   *hoc_objgetarg(1); // list of vectors
    num = ivoc_list_count(ob);
    if (num>NSW) { printf("FIELD initwrec() WARN: can only store %d ww vecs\n",NSW); hxe();}
    nsw=(double)num;
    for (k=0;k<num;k++) {
      cap = list_vector_px2(ob, k, &wwo[k], &ww[k]);
      if (k==0) wwsz=cap; else if (wwsz!=cap) {
        printf("FIELD initwrec ERR w-vecs size err: %d,%d,%d",k,wwsz,cap); hxe(); }
    }
  }
  ENDVERBATIM
}

BREAKPOINT {
  VERBATIM
  int k;
  if (t>wwo[0][ixp1]) { // time has moved on
    ixp1++;
    wwo[0][ixp1]=t;
  }
  if (ixp1>=wwsz) { // grow
    if (wwsz==0) wwsz=1e4; else wwsz*=2;
    for (k=0;k<(int)nsw;k++) wwo[k]=vector_newsize(ww[k],wwsz);
  }
  wwo[(int)fldID][ixp1]+=(v-p1);
  ENDVERBATIM  
}

PROCEDURE global_fini () {
  VERBATIM
  int k;
  for (k=0;k<(int)nsw;k++) vector_resize(ww[k], ixp1);
  ENDVERBATIM
}
