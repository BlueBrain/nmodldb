: $Id: vecst.mod,v 1.178 2005/03/07 22:26:47 billl Exp $
 
:* COMMENT
COMMENT
randwd   randomly chooses n bits to set to 1
hamming  v.hamming(v1) is hamming distance between 2 vecs
flipbits v.flipbits(scratch,num) flips num rand chosen bits
flipbalbits v.flipbalbits(scratch,num) balanced flipping
vpr      v.vpr prints out vector as 1 (x[i]>0) or 0 (x[i]<=0)
thresh   turns analog vec to BVBASE,1 vec separating at thresh (scalar or vec)
triplet  return location of a triplet in a vector
onoff    turns state vec on or off depending on values in other vecs
bpeval   used by backprop: vo=outp*(1-outp)*del
w        like where but sets chosen nums // modified from .wh in 1.23
xing     determines where vector crosses in a positive direction (destructive)
xzero    count 0 crossings by putting in a 1 or -1
indset(ind,val) sets spots indexed by ind to val
ismono([arg]) checks if a vector is monotonically increaseing (-1 decreasing)
count(num) count the number of nums in vec
scr.fewind(ind,veclist) // uses ind as index into other vecs
ind.findx(vecAlist,vecBlist) // uses ind as index into vecA's to load in vecB's
ind.sindx(vecAlist,vecBlist) // replace ind elements in vecA's with vecB's values
ind.sindv(vecAlist,valvec) // replace ind elements in vecA's with vals
scr.nind(ind,vec1,vec2[,vec3,vec4]) // uses ind to elim elems in other vecs
ind.keyind(key,vec1,vec2[,vec3,vec4]) // pick out bzw. values from vectors
ind.slct(key,args,veclist) // pick out bzw. values from vectors
ind.slor(key,args,veclist) // do OR rather than AND function of slct
vdest.intrp(flag) // interpolate numbers replacing numbers given as flag
v.insct(v1,v2)   // return v1 intersect v2
vdest.cull(vsrc,key)   // remove values found in key
vdest.redundout(vsrc[,INDFLAG])  // remove repeat values, with flag return indices
v.cvlv(v1,v2)   // convolve v1 with v2
fac      not vec related - returns factorial
logfac   not vec related - returns log factorial
v2d      copy into a double block
d2v      copy out of a double block NB: same as vec.from_double()
smgs     rewrite of sumgauss form ivovect.cpp
iwr      write integers
ird      read integers
info     give pointer addresses and sizes
lcat     concatentate all vectors from a list
Non-vector routines
vseed    set some C level randomizer seeds
isojt    compare 2 objects to see if they're of the same type
eqojt    compare 2 object pointers to see if they point to same thing
ENDCOMMENT

NEURON {
  SUFFIX nothing
  GLOBAL BVBASE, RES, VECST_INSTALLED  : bit vector base number (typically 0 or -1)
}

PARAMETER {
  BVBASE = -1.
  VECST_INSTALLED=0
  : misc
  ERR=-1.3479e121

  : 0 args
  ALL=-1.3479e120
  NEG=-1.3478e120
  POS=-1.3477e120
  CHK=-1.3476e120
  NOZ=-1.3475e120
  : 1 arg
  GTH=-1.3474e120
  GTE=-1.3473e120
  LTH=-1.3472e120
  LTE=-1.3471e120
  EQU=-1.3470e120
  EQV=-1.3469e120 : value equal to same row value in parallel vector
  EQW=-1.3468e120 : value found in other vector
  NEQ=-1.3467e120
  SEQ=-1.3466e120
  RXP=-1.3465e120
  : 2 args
  IBE=-1.3464e120
  EBI=-1.3463e120
  IBI=-1.3462e120
  EBE=-1.3461e120
}

ASSIGNED { RES }

VERBATIM
#include <stdlib.h>
#include <math.h>
#include <values.h> /* contains MAXLONG */
#include <sys/time.h> 
extern double* hoc_pgetarg();
extern double hoc_call_func(Symbol*, int narg);
extern FILE* hoc_obj_file_arg(int narg);
extern Object** hoc_objgetarg();
extern void vector_resize();
extern int vector_instance_px();
extern void* vector_arg();
extern double* vector_vec();
extern double hoc_epsilon;
extern void set_seed();
extern int ivoc_list_count(Object*);
extern Object* ivoc_list_item(Object*, int);
static int list_vector_px();
static int list_vector_px2();
static int list_vector_resize();

typedef struct BVEC {
 int size;
 int bufsize;
 short *x;
 Object* o;
} bvec;
ENDVERBATIM
 
:* v1.slope(num) does a linear regression to find the slope, assuming num=timestep of vector

VERBATIM
static double slope(void* vv) {
	int i, n;
	double *x, *y;
        double timestep, sigxy, sigx, sigy, sigx2;
	/* how to get the instance data */
	n = vector_instance_px(vv, &y);

        if(ifarg(1)) { 
          timestep = *getarg(1); 
        } else { printf("You must supply a timestep\n"); return 0; }

        sigxy= sigx= sigy= sigx2=0; // initialize these

        x = (double *) malloc(sizeof(double)*n);
        for(i=0; i<n; i++) {
          x[i] = timestep*i;
          sigxy += x[i] * y[i];
          sigx  += x[i];
          sigy  += y[i];
          sigx2 += x[i]*x[i];
        }
        return (n*sigxy - sigx*sigy)/(n*sigx2 - sigx*sigx);
}
ENDVERBATIM
 
:* v1.vslope(v2) does a linear regression, using v2 as the x-coords

VERBATIM
static double vslope(void* vv) {
	int i, n;
	double *x, *y;
        double timestep, sigxy, sigx, sigy, sigx2;
	/* how to get the instance data */
	n = vector_instance_px(vv, &y);

        if(ifarg(1)) {
          if(vector_arg_px(1, &x) != n ) {
            hoc_execerror("Vector size doesn't match.", 0); 
          }
          sigxy= sigx= sigy= sigx2=0; // initialize these

          for(i=0; i<n; i++) {
            sigxy += x[i] * y[i];
            sigx  += x[i];
            sigy  += y[i];
            sigx2 += x[i]*x[i];
          }
        }         
        return (n*sigxy - sigx*sigy)/(n*sigx2 - sigx*sigx);
}
ENDVERBATIM
 
:* v1.stats(num) does a linear regression, assuming num=timestep of vector

VERBATIM
static double stats(void* vv) {
	int i, n;
	double *x, *y;
        double timestep, sigxy, sigx, sigy, sigx2, sigy2;
        double r, m, b;
	/* how to get the instance data */
	n = vector_instance_px(vv, &y);

        if(ifarg(1)) { 
          timestep = *getarg(1); 
        } else { printf("You must supply a timestep\n"); return 0; }

        sigxy= sigx= sigy= sigx2=sigy2= 0; // initialize these

        x = (double *) malloc(sizeof(double)*n);
        for(i=0; i<n; i++) {
          x[i] = timestep*i;
          sigxy += x[i] * y[i];
          sigx  += x[i];
          sigy  += y[i];
          sigx2 += x[i]*x[i];
          sigy2 += y[i]*y[i];
        }
        m = (n*sigxy - sigx*sigy)/(n*sigx2 - sigx*sigx);
        b = (sigy*sigx2 - sigx*sigxy)/(n*sigx2 - sigx*sigx);
        r = (n*sigxy - sigx*sigy)/(sqrt(n*sigx2-sigx*sigx) * sqrt(n*sigy2-sigy*sigy));

        printf("Examined %d data points\n", n);
        printf("slope     = %f\n", m);
        printf("intercept = %f\n", b);
        printf("R         = %f\n", r);
        printf("R-squared = %f\n", r*r);
        return 1;
}
ENDVERBATIM
 
:* v1.vstats(v2) does a linear regression, using v2 as the x-coords

VERBATIM
static double vstats(void* vv) {
	int i, n;
	double *x, *y;
        double timestep, sigxy, sigx, sigy, sigx2, sigy2;
        double r, m, b;
	/* how to get the instance data */
	n = vector_instance_px(vv, &y);

        if(ifarg(1)) {
          if(vector_arg_px(1, &x) != n ) {
            hoc_execerror("Vector size doesn't match.", 0); 
          }
          sigxy= sigx= sigy= sigx2=sigy2=0; // initialize these

          for(i=0; i<n; i++) {
            sigxy += x[i] * y[i];
            sigx  += x[i];
            sigy  += y[i];
            sigx2 += x[i]*x[i];
            sigy2 += y[i]*y[i];
          }
          m = (n*sigxy - sigx*sigy)/(n*sigx2 - sigx*sigx);
          b = (sigy*sigx2 - sigx*sigxy)/(n*sigx2 - sigx*sigx);
          r = (n*sigxy - sigx*sigy)/(sqrt(n*sigx2-sigx*sigx) * sqrt(n*sigy2-sigy*sigy));

          printf("Examined %d data points\n", n);
          printf("slope     = %f\n", m);
          printf("intercept = %f\n", b);
          printf("R         = %f\n", r);
          printf("R-squared = %f\n", r*r);
          return 1;
        } else {
          printf("You must supply an x vector!\n");
          return 0;
        }
}
ENDVERBATIM
 
:* v1.randwd(num[,v2]) will randomly flip num bits from BVBASE to 1
: does v1.fill(BVBASE); optionally fill v2 with the indices
VERBATIM
static double randwd(void* vv) {
	int i, ii, jj, nx, ny, flip, flag;
	double* x, *y;
	/* how to get the instance data */
	nx = vector_instance_px(vv, &x);
        flip = (int) *getarg(1);
        if (ifarg(2)) { /* write a diff vector to z */
          flag = 1; ny = vector_arg_px(2, &y);
          if (ny!=flip) { hoc_execerror("Opt vector must be size for # of flips", 0); }
        } else { flag = 0; }
        if (flip>=nx) { hoc_execerror("# of flips exceeds (or ==) vector size", 0); }
	for (i=0; i < nx; i++) { x[i] = BVBASE; }
	for (i=0,jj=0; i < flip; i++) { /* flip these bits */
	  ii = (int) ((nx+1)*drand48());
	  if (x[ii]==BVBASE) {
	    x[ii] = 1.; 
            if (flag) { y[jj] = ii; jj++; }
	  } else {
	    i--;
	  }
	}
	return flip;
}
ENDVERBATIM
 
:* v1.hamming(v2[,v3]) compares v1 and v2 for matches, v3 gives diff vector
VERBATIM
static double hamming(void* vv) {
	int i, nx, ny, nz, prflag;
	double* x, *y, *z,sum;
	sum = 0.;
	nx = vector_instance_px(vv, &x);
	ny = vector_arg_px(1, &y);
        if (ifarg(2)) { /* write a diff vector to z */
          prflag = 1; nz = vector_arg_px(2, &z);
        } else { prflag = 0; }
	if (nx!=ny || (prflag && nx!=nz)) {
	  hoc_execerror("Vectors must be same size", 0);
	}
	for (i=0; i < nx; ++i) {
	  if (x[i] != y[i]) { sum++; 
            if (prflag) { z[i] = 1.; }
          } else if (prflag) { z[i] = 0.; }
        }
	return sum;
}
ENDVERBATIM
 
:* v1.info() gives addresses and sizes
VERBATIM
static double info(void* vv) {
  int nx,bsz; double* x;
  nx = vector_instance_px(vv, &x);
  bsz=vector_buffer_size(vv);
  printf("Obj*%x Dbl*%x Size: %d Bufsize: %d\n",vv,x,nx,bsz);
}
ENDVERBATIM
 
:* v1.indset(ind,x[,y]) sets indexed values to x and other values to optional y
VERBATIM
static double indset(void* vv) {
	int i, nx, ny, nz, flag;
	double* x, *y, *z, val, val2 ;
	nx = vector_instance_px(vv, &x);
	ny = vector_arg_px(1, &y);
        if (hoc_is_object_arg(2)) { 
          flag=1;
          nz = vector_arg_px(2, &z); 
          if (ny!=nz) { hoc_execerror("v.indset: Vector sizes don't match.", 0); }
        } else { flag=0; val=*getarg(2); }
        if (ifarg(3)) { 
          val2 = *getarg(3); 
          for (i=0; i<nx; i++) { x[i]=val2; }
        }
	for (i=0; i<ny; i++) {
	  if (y[i] > nx) { hoc_execerror("v.indset: Index exceeds vector size", 0); }
          if (flag) x[(int)y[i]]=z[i]; else x[(int)y[i]]=val;
        }
	return i;
}
ENDVERBATIM
 
VERBATIM
/* Maintain parallel int vector to avoid slowness of repeated casts */
static int scrsz=0;
static int *scr;
static double dcr[100]; // scratch area for doubles
ENDVERBATIM

:* tmp.fewind(ind,veclist)
: picks out numbers from multiple vectors using index ind
VERBATIM
static double fewind (void* vv) {
  int i, j, nx, ni, nv[10], num;
  Object* ob;
  double *x, *ind, *vvo[10];
  nx = vector_instance_px(vv, &x);
  ni = vector_arg_px(1, &ind);
  ob = *hoc_objgetarg(2);
  num = ivoc_list_count(ob);
  if (num>10) hoc_execerror("ERR: fewind can only handle 10 vectors", 0);
  for (i=0;i<num;i++) { 
    nv[i] = list_vector_px(ob, i, &vvo[i]);
    if (nx!=nv[i]) { printf("fewind ERR %d %d %d\n",i,nx,nv[i]);
      hoc_execerror("Vectors must all be same size: ", 0); }
  }
  if (ni>scrsz) { 
    if (scrsz>0) { free(scr); scr=(int *)NULL; }
    scrsz=ni+10000;
    scr=(int *)ecalloc(scrsz, sizeof(int));
  }
  for (i=0;i<ni;i++) { scr[i]=(int)ind[i]; /* copy into integer array */
    if (scr[i]>=nx || scr[i]<0) { printf("fewind ERR1 %d %d\n",scr[i],nx);
    hoc_execerror("Index vector out-of-bounds", 0); }
  }
  for (j=0;j<num;j++) {
    for (i=0;i<ni;i++) x[i]=vvo[j][scr[i]];    
    for (i=0;i<ni;i++) vvo[j][i]=x[i];   
    list_vector_resize(ob, j, ni);
  }
  return ni;
}
ENDVERBATIM


:* ind.findx(vecAlist,vecBlist) 
: uses ind as index into vecA's to load in vecB's (for select); a nondestructive fewind()
VERBATIM
static double findx (void* vv) {
  int i, j, ni, nx, av[11], bv[11], num;
  Object *ob1, *ob2;
  double *ind, *avo[11], *bvo[11];
  ni = vector_instance_px(vv, &ind);
  ob1 = *hoc_objgetarg(1);
  ob2 = *hoc_objgetarg(2);
  num = ivoc_list_count(ob1);
  i = ivoc_list_count(ob2);
  if (i!=num) hoc_execerror("findx ****ERRA****: lists have different counts", 0);
  if (num>11) hoc_execerror("findx ****ERRB****: can only handle 11 vectors", 0);
  for (i=0;i<num;i++) { 
    av[i]=list_vector_px(ob1, i, &avo[i]); // source vectors
    if (av[0]!=av[i]) { printf("findx ****ERRC**** %d %d %d\n",i,av[0],av[i]);
      hoc_execerror("Src vectors must all be same size: ", 0); }
  }
  nx=av[0]; // size of source vecs
  for (i=0;i<num;i++) { 
    bv[i]=list_vector_px2(ob2, i, &bvo[i], &vv); /* dest vectors */
    if (vector_buffer_size(vv)<ni) { 
      printf("findx ****ERRD**** arg#%d need:%d sz:%d\n",num+i+1,ni,vector_buffer_size(vv));
      hoc_execerror("Destination vector with insufficient size: ", 0); 
    } else {
      vector_resize(vv, ni);
    }
  }
  if (ni>scrsz) { 
    if (scrsz>0) { free(scr); scr=(int *)NULL; }
    scrsz=ni+10000;
    scr=(int *)ecalloc(scrsz, sizeof(int));
  }
  for (i=0;i<ni;i++) { scr[i]=(int)ind[i]; // copy into integer array
    if (scr[i]>=nx || scr[i]<0) { printf("findx ****ERRE**** **** IND:%d SZ:%d\n",scr[i],nx);
      hoc_execerror("Index vector out-of-bounds", 0); }
  }
  for (j=0,i=0;j<num;j++) for (i=0;i<ni;i++) bvo[j][i]=avo[j][scr[i]];    
  return ni;
}
ENDVERBATIM

:* ind.sindx(vecAlist,vecBlist)
: uses ind as index into vecA's to replace with elements from vecB's
VERBATIM
static double sindx (void* vv) {
  int i, j, ni, nx, av[11], bv[11], num;
  Object *ob1, *ob2;
  double *ind, *avo[11], *bvo[11];
  ni = vector_instance_px(vv, &ind);
  ob1 = *hoc_objgetarg(1);
  ob2 = *hoc_objgetarg(2);
  num = ivoc_list_count(ob1);
  i = ivoc_list_count(ob2);
  if (num!=i) hoc_execerror("sindx ****ERRA****: two vec lists have different counts", 0);
  if (num>11) hoc_execerror("sindx ****ERRB****: can only handle 11 vectors", 0);
  for (i=0;i<num;i++) { 
    av[i]=list_vector_px(ob1, i, &avo[i]); // dest vectors
    if (av[0]!=av[i]) { printf("sindx ****ERRC**** %d %d %d\n",i,av[0],av[i]);
      hoc_execerror("Dest. vectors must all be same size: ", 0); }
  }
  nx=av[0]; // size of dest vecs
  for (i=0;i<num;i++) { 
    bv[i]=list_vector_px(ob2, i, &bvo[i]); /* source vectors */
  if (bv[i]!=ni) { 
    printf("sindx ****ERRD**** arg#%d does note match ind length %d vs %d\n",num+i+1,ni,bv[i]);
    hoc_execerror("Source vector with insufficient size: ", 0); }
  }
  if (ni>scrsz) { 
    if (scrsz>0) { free(scr); scr=(int *)NULL; }
    scrsz=ni+10000;
    scr=(int *)ecalloc(scrsz, sizeof(int));
  }
  for (i=0;i<ni;i++) { scr[i]=(int)ind[i]; // copy into integer array
    if (scr[i]>=nx || scr[i]<0) { 
      printf("sindx ****ERRE**** IND:%d SZ:%d\n",scr[i],nx);
      hoc_execerror("Index vector out-of-bounds", 0); }
  }
  for (j=0,i=0;j<num;j++) for (i=0;i<ni;i++) avo[j][scr[i]]=bvo[j][i];
  return ni;
}

ENDVERBATIM
:* ind.sindx(vecAlist,valvec)
: uses ind as index into vecA's to replace with values in valvec
VERBATIM
static double sindv(void* vv) {
  int i, j, ni, nx, av[11], bv, num;
  Object* ob;
  double *ind, *avo[11], *bvo;
  ni = vector_instance_px(vv, &ind);
  ob = *hoc_objgetarg(1);
  bv=vector_arg_px(2, &bvo); // source vector
  num = ivoc_list_count(ob);
  if (num>11) hoc_execerror("sindv ****ERRA****: can only handle 11 vectors", 0);
  for (i=0;i<num;i++) { 
    av[i]=list_vector_px(ob, i, &avo[i]); // dest vectors
    if (av[0]!=av[i]) { printf("sindv ****ERRC**** %d %d %d\n",i,av[0],av[i]);
      hoc_execerror("Dest. vectors must all be same size: ", 0); }
  }
  nx=av[0]; // size of source vecs
  if (bv!=num) { 
    printf("sindv ****ERRD**** Vector arg does note match list count %d vs %d\n",num,bv);
    hoc_execerror("Source vector is wrong size: ", 0); }
  if (ni>scrsz) { 
    if (scrsz>0) { free(scr); scr=(int *)NULL; }
    scrsz=ni+10000;
    scr=(int *)ecalloc(scrsz, sizeof(int));
  }
  for (i=0;i<ni;i++) { scr[i]=(int)ind[i]; // copy into integer array
    if (scr[i]>=nx || scr[i]<0) { 
      printf("sindv ****ERRE**** IND:%d SZ:%d\n",scr[i],nx);
      hoc_execerror("Index vector out-of-bounds", 0); }
  }
  for (j=0,i=0;j<num;j++) for (i=0;i<ni;i++) avo[j][scr[i]]=bvo[j];
  return ni;
}
ENDVERBATIM

:* ind.slct(key,args,veclist)
: picks out indices of numbers in key from multiple vectors
VERBATIM
static double slct(void* vv) {
  int i, j, k, m, n, p, ni, nk, na, nv[10], num, fl, lc;
  Object* lob;
  double *ind, *key, *arg, *vvo[10];

  ni = vector_instance_px(vv, &ind); // vv is ind
  nk = vector_arg_px(1, &key);
  na = vector_arg_px(2, &arg);
  lob = *hoc_objgetarg(3);
  num = ivoc_list_count(lob);
  if (num>10) hoc_execerror("ERR: vecst::slct can only handle 10 vectors", 0);

  for (i=0,j=0,k=0;i<num;i++,j++) { 
    nv[i] = list_vector_px(lob, i, &vvo[i]);
    if (ni!=nv[i] && (key[j]!=EQW || k!=1)) {
      printf("vecst::slct ERR %d %d %d %d %d\n",i,j,k,ni,nv[i]);
      hoc_execerror("index and searched vectors must all be same size: ", 0); 
    }
    if (key[j]==EQW || key[j]==EQV) if (k==0){j--;k++;} else k=0; // EQW,EQV take 2 vector args 
  }
  if (2*nk!=na) { printf("vecst::slct ERR3 %d %d\n",nk,na); 
    hoc_execerror("Arg vector must be double key length",0); }
  for (i=0,n=0;i<nk;i++) if (key[i]==EQV || key[i]==EQW) n++; // special cases take 2 vec args
  if (nk+n!=num) { 
    printf("vecst::slct ERR2 %d(keys)+%d(EQV/W)!=%d(vecs)\n",nk,n,num); 
    hoc_execerror("Key length must be number of vecs + num of EQV/W",0); }
  for (j=0,k=0,m=0;j<ni;j++) { // j steps through elements of vectors
    for (i=0,m=0,n=0,fl=1;i<num;i++,n++,m+=2) { // i steps thru key, m thru args
      if (key[n]==ALL) continue; // OK - do nothing
      if (key[n]==NOZ)        { if (vvo[i][j]==0.) {fl=0; break;} else continue; 
      } else if (key[n]==POS) { if (vvo[i][j]<=0.) {fl=0; break;} else continue; 
      } else if (key[n]==NEG) { if (vvo[i][j]>=0.) {fl=0; break;} else continue; 
      } else if (key[n]==GTH) { if (vvo[i][j]<=arg[m]) {fl=0; break;} else continue; 
      } else if (key[n]==GTE) { if (vvo[i][j]< arg[m]) {fl=0; break;} else continue; 
      } else if (key[n]==LTH) { if (vvo[i][j]>=arg[m]) {fl=0; break;} else continue; 
      } else if (key[n]==LTE) { if (vvo[i][j]> arg[m]) {fl=0; break;} else continue; 
      } else if (key[n]==EQU) { if (vvo[i][j]!=arg[m]) {fl=0; break;} else continue; 
      } else if (key[n]==EQV) { if (vvo[i][j]!=vvo[i+1][j]) { 
          fl=0; break;} else { i++; continue; }
      } else if (key[n]==EQW) {  // check value against values in following vec
        fl=0; // assume it's not going to match
        for (p=0;p<nv[i+1];p++) if (vvo[i][j]==vvo[i+1][p]) {fl=1; break;}
        if (fl==0) break; else { i++; continue; }
      } else if (key[n]==NEQ) { if (vvo[i][j]==arg[m]) {fl=0; break;} else continue; 
      } else if (key[n]==IBE) { if ((vvo[i][j]< arg[m])||(vvo[i][j]>=arg[m+1])) {
          fl=0; break; } else continue;  // IBE="[)" include-bracket-exclude
      } else if (key[n]==EBI) { if ((vvo[i][j]<=arg[m])||(vvo[i][j]> arg[m+1])) { 
          fl=0; break; } else continue;   // "(]" : exclude-bracket-include
      } else if (key[n]==IBI) { if ((vvo[i][j]< arg[m])||(vvo[i][j]> arg[m+1])) {
          fl=0; break; } else continue;   // "[]" : include-bracket-include
      } else if (key[n]==EBE) { if ((vvo[i][j]<=arg[m])||(vvo[i][j]>=arg[m+1])) {
          fl=0; break; } else continue;   // "()" : exclude-bracket-exclude
      } else {printf("vecst::slct ERR4 %g\n",key[n]); hoc_execerror("Unknown key",0);}
    }
    if (fl) ind[k++]=j; // all equal
  }
  vector_resize(vv, k);
  return k;
}
ENDVERBATIM

:* ind.slor(key,args,vec1,vec2[,vec3,vec4,...])
: picks out indices of numbers in key from multiple vectors
VERBATIM
static double slor(void* vv) {
  int i, j, k, m, n, p, ni, nk, na, nv[10], num, fl;
  Object* lob;
  double *ind, *key, *arg, *vvo[10];

  ni = vector_instance_px(vv, &ind); // vv is ind
  nk = vector_arg_px(1, &key);
  na = vector_arg_px(2, &arg);
  lob = *hoc_objgetarg(3);
  num = ivoc_list_count(lob);
  if (num>10) hoc_execerror("ERR: vecst::slor can only handle 10 vectors", 0);

  for (i=0,j=0,k=0;i<num;i++,j++) { 
    nv[i] = list_vector_px(lob, i, &vvo[i]);
    if (ni!=nv[i] && (key[j]!=EQW || k!=1)) {
      printf("vecst::slct ERR %d %d %d %d %d\n",i,j,k,ni,nv[i]);
      hoc_execerror("index and searched vectors must all be same size: ", 0); 
    }
    if (key[j]==EQW || key[j]==EQV) if (k==0){j--;k++;} else k=0; // k counts 2 vecs
  }
  if (2*nk!=na) { printf("vecst::slor ERR3 %d %d\n",nk,na); 
    hoc_execerror("Arg vector must be double key length",0); }
  for (i=0,n=0;i<nk;i++) if (key[i]==EQV || key[i]==EQW) n++; // special case takes 2 vec args
  if (nk+n!=num) { 
    printf("vecst::slor ERR2 %d(keys)+%d(EQV)!=%d(vecs)\n",nk,n,num); 
    hoc_execerror("Key length must be number of vecs + num of EQV",0); }
  for (j=0,k=0,m=0;j<ni;j++) { // j steps through elements of vectors
    for (i=0,m=0,n=0,fl=0;i<num;i++,n++,m+=2) { // i steps thru key, m thru args
      if (key[n]==ALL) {fl=1; break;} // OK - do nothing
      if (key[n]==NOZ)        { if (vvo[i][j]==0.) continue; else {fl=1; break;} 
      } else if (key[n]==POS) { if (vvo[i][j]<=0.) continue; else {fl=1; break;} 
      } else if (key[n]==NEG) { if (vvo[i][j]>=0.) continue; else {fl=1; break;} 
      } else if (key[n]==GTH) { if (vvo[i][j]<=arg[m]) continue; else {fl=1; break;} 
      } else if (key[n]==GTE) { if (vvo[i][j]< arg[m]) continue; else {fl=1; break;} 
      } else if (key[n]==LTH) { if (vvo[i][j]>=arg[m]) continue; else {fl=1; break;} 
      } else if (key[n]==LTE) { if (vvo[i][j]> arg[m]) continue; else {fl=1; break;} 
      } else if (key[n]==EQU) { if (vvo[i][j]!=arg[m]) continue; else {fl=1; break;} 
      } else if (key[n]==EQV) { if (vvo[i][j]!=vvo[i+1][j]) continue; else { 
          i++; fl=1; break; }
      } else if (key[n]==EQW) {  // check value against values in following vec
        fl=0; // assume it's not going to match
        for (p=0;p<nv[i+1];p++) if (vvo[i][j]==vvo[i+1][p]) {fl=1; break;}
        if (fl==1) break; else { i++; continue; }
      } else if (key[n]==NEQ) { if (vvo[i][j]==arg[m]) continue; else {fl=1; break;} 
      } else if (key[n]==IBE) { if ((vvo[i][j]< arg[m])||(vvo[i][j]>=arg[m+1])) {
          continue; } else {fl=1; break;}  // IBE="[)" include-bracket-exclude
      } else if (key[n]==EBI) { if ((vvo[i][j]<=arg[m])||(vvo[i][j]> arg[m+1])) { 
          continue; } else {fl=1; break;}   // "(]" : exclude-bracket-include
      } else if (key[n]==IBI) { if ((vvo[i][j]< arg[m])||(vvo[i][j]> arg[m+1])) {
          continue; } else {fl=1; break;}   // "[]" : include-bracket-include
      } else if (key[n]==EBE) { if ((vvo[i][j]<=arg[m])||(vvo[i][j]>=arg[m+1])) {
          continue; } else {fl=1; break;}   // "()" : exclude-bracket-exclude
      } else {printf("vecst::slor ERR4 %g\n",key[n]); hoc_execerror("Unknown key",0);}
    }
    if (fl) ind[k++]=j; // all equal
  }
  vector_resize(vv, k);
  return k;
}
ENDVERBATIM

:* v.iwr()
VERBATIM
static double iwr(void* vv) {
  int i, j, nx;
  double *x;
  FILE* f, *hoc_obj_file_arg();
  f = hoc_obj_file_arg(1);
  nx = vector_instance_px(vv, &x);
  if (nx>scrsz) { 
    if (scrsz>0) { free(scr); scr=(int *)NULL; }
    scrsz=nx+10000;
    scr=(int *)ecalloc(scrsz, sizeof(int));
  }
  for (i=0;i<nx;i++) scr[i]=(int)x[i]; /* copy into integer array */
  fwrite(&nx,sizeof(int),1,f);  // write out the size
  fwrite(scr,sizeof(int),nx,f);
  return nx;
}
ENDVERBATIM

:* v.ird()
VERBATIM
static double ird(void* vv) {
  int i, j, nx, n;
  double *x;
  FILE* f, *hoc_obj_file_arg();
  f = hoc_obj_file_arg(1);
  nx = vector_instance_px(vv, &x);
  fread(&n,sizeof(int),1,f);  // size
  if (n>scrsz) { 
    if (scrsz>0) { free(scr); scr=(int *)NULL; }
    scrsz=n+10000;
    scr=(int *)ecalloc(scrsz, sizeof(int));
  }
  if (n!=nx) { 
    nx=vector_buffer_size(vv);
    if (n<=nx) {
      vector_resize(vv, n); nx=n; 
    } else {
      printf("%d > %d :: ",n,nx);
      hoc_execerror("Vector max capacity too small for ird ", 0);
    }
  }
  fread(scr,sizeof(int),n,f);
  for (i=0;i<nx;i++) x[i]=(double)scr[i];
  return n;
}
ENDVERBATIM

:* ind.insct(v1,v2)
: return v1 intersect v2
VERBATIM
static double insct(void* vv) {
	int i, j, k, nx, nv1, nv2, maxsz;
	double *x, *v1, *v2;
	nx = vector_instance_px(vv, &x);
        maxsz=vector_buffer_size(vv);
        vector_resize(vv, maxsz);
	nv1 = vector_arg_px(1, &v1);
	nv2 = vector_arg_px(2, &v2);
        for (i=0,k=0;i<nv1;i++) for (j=0;j<nv2;j++) if (v1[i]==v2[j]) {
          if (k<maxsz) { x[k++]=v1[i]; } else {k++;}}  /* v1[i] found in both vectors */
        if (k>maxsz) { 
          printf("\tinsct WARNING: ran out of room: %d<%d\n",maxsz,k);
        } else { vector_resize(vv, k); }
	return (double)k;
}
ENDVERBATIM

:* vec.cull(src,key)
: remove numbers in vec that are found in the key
VERBATIM
static double cull(void* vv) {
	int i, j, k, nx, nv1, nv2, maxsz, flag;
	double *x, *v1, *v2;
	nx = vector_instance_px(vv, &x);
        maxsz=vector_buffer_size(vv);
        vector_resize(vv, maxsz);
	nv1 = vector_arg_px(1, &v1);
	nv2 = vector_arg_px(2, &v2);
        for (i=0,k=0;i<nv1;i++) {
          flag=1;
          for (j=0;j<nv2;j++) if (v1[i]==v2[j]) flag=0;
          if (flag) {if (k<maxsz) { x[k++]=v1[i]; } else { k++; }}
        }
        if (k>maxsz) { 
          printf("\tcull WARNING: ran out of room: %d<%d\n",maxsz,k);
        } else { vector_resize(vv, k); }
	return (double)k;
}
ENDVERBATIM

:* dest.redundout(src[,INDFLAG])
: flag redundant numbers; must sort src first
VERBATIM
static double redundout(void* vv) {
	int i, j, nx, nv1, maxsz, indflag;
	double *x, *v1, val;
        if (ifarg(2)) indflag=1; else indflag=0; 
	nx = vector_instance_px(vv, &x);
        maxsz=vector_buffer_size(vv);
        vector_resize(vv, maxsz);
	nv1 = vector_arg_px(1, &v1);
        val=v1[0]; x[0]=(indflag?0:val);
        for (j=1,i=1;i<nv1&&j<maxsz;i++) if (v1[i]!=val) { val=v1[i]; x[j++]=(indflag?i:val); }
        if (j==maxsz) { 
          printf("\tredundout WARNING: ran out of room: %d<needed\n",maxsz);
        } else { vector_resize(vv, j); }
	return (double)j;
}
ENDVERBATIM

:* vdest.cvlv(vsrc,vfilt)
: convolution
VERBATIM
static double cvlv(void* vv) {
  int i, j, k, nx, nsrc, nfilt;
  double *x, *src, *filt, sum, lpad, rpad;
  nx = vector_instance_px(vv, &x);
  nsrc = vector_arg_px(1, &src);
  nfilt = vector_arg_px(2, &filt);
  if (nx!=nsrc) { hoc_execerror("Vectors not same size: ", 0); }
  if (nfilt>nsrc) { hoc_execerror("Filter bigger than source ", 0); }
  for (i=0;i<nx;i++) {
    x[i]=0.0;
    for (j=0,k=i-(int)(nfilt/2);j<nfilt;j++,k++) {
      if (k>0 && k<nsrc-1) x[i]+=filt[j]*src[k];
    }
  }
}
ENDVERBATIM

:* vdest.intrp(flag)
: interpolate numbers replacing numbers given as flag
VERBATIM
static double intrp(void* vv) {
  int i, la, lb, nx;
  double *x, fl, a, b;
  nx = vector_instance_px(vv, &x);
  fl = *getarg(1);
  i=0; a=x[0]; la=0;
  if (a==fl) a=0;
  while (i<nx-1) {
    for (i=la+1;x[i]==fl && i<nx-1; i++) ; /* find the next one */
    b=x[i]; lb=i;
    for (i=la+1; i<lb; i++) x[i]= a + (b-a)/(lb-la)*(i-la);
    a=b; la=lb;
  }
}
ENDVERBATIM

:* tmp.nind(ind,vec1,vec2[,vec3,vec4,...])
: picks out numbers not in ind from multiple vectors
: ind must be sorted
VERBATIM
static double nind(void* vv) {
	int i, j, k, m, nx, ni, nv[10], num, c, last;
	double *x, *ind, *vvo[10];
	nx = vector_instance_px(vv, &x);
        for (i=0;ifarg(i);i++);
	if (i>9) hoc_execerror("ERR: nind can only handle 9 vectors", 0);
	num = i-2; /* number of vectors to be picked apart */
        for (i=0;i<num;i++) { 
          nv[i] = vector_arg_px(i+2, &vvo[i]);
          if (nx!=nv[i]) { printf("nind ERR %d %d %d\n",i,nx,nv[i]);
            hoc_execerror("Vectors must all be same size: ", 0); }
        }
        ni = vector_arg_px(1, &ind);
        c = nx-ni; /* the elems indexed are to be eliminated */
        if (ni>scrsz) { 
          if (scrsz>0) { free(scr); scr=(int *)NULL; }
          scrsz=ni+10000;
          scr=(int *)ecalloc(scrsz, sizeof(int));
        }
        for (i=0,last=-1;i<ni;i++) { 
          scr[i]=(int)ind[i]; /* copy into integer array */
          if (scr[i]<0 || scr[i]>=nx) hoc_execerror("nind(): Index out of bounds", 0);
          if (scr[i]<=last) hoc_execerror("nind(): indices should mono increase", 0);
          last=scr[i];
        }
        for (j=0;j<num;j++) { /* each output vec */
          for (i=0,last=-1,m=0;i<ni;i++) { /* the indices of ind */
            for (k=last+1;k<scr[i];k++) { x[m++]=vvo[j][k]; }
            last=scr[i];
          }
          for (k=last+1;k<nx;k++,m++) { x[m]=vvo[j][k]; }
          for (i=0;i<c;i++) vvo[j][i]=x[i];   
          vv=vector_arg(j+2); vector_resize(vv, c);
        }
	return c;
}
ENDVERBATIM

:* ind.keyind(key,vec1,vec2[,vec3,vec4,...])
: picks out indices of numbers in key from multiple vectors
VERBATIM
static double keyind(void* vv) {
  int i, j, k, ni, nk, nv[10], num;
  double *ind, *key, *vvo[10];
  ni = vector_instance_px(vv, &ind); // vv is ind
  for (i=0;ifarg(i);i++); i--; // drop back by one to get numarg()
  if (i>10) hoc_execerror("ERR: keyind can only handle 9 vectors", 0);
  num = i-1; /* number of vectors to be picked apart */
  for (i=0;i<num;i++) { 
    nv[i] = vector_arg_px(i+2, &vvo[i]);
    if (ni!=nv[i]) { printf("keyind ERR %d %d %d\n",i,ni,nv[i]);
    hoc_execerror("Non-key vectors must be same size: ", 0); }
  }
  nk = vector_arg_px(1, &key);
  if (nk!=num) { printf("keyind ERR2 %d %d\n",nk,num); 
    hoc_execerror("Key length must be number of vecs",0); }
  k=0;
  for (j=0;j<ni;j++) { // j steps through elements of vectors
    for (i=0;i<nk;i++) { // i steps through the key
      if (key[i]==ALL) continue; // OK - do nothing
      if (key[i]==NOZ)        { if (vvo[i][j]==0.) break; else continue; 
      } else if (key[i]==POS) { if (vvo[i][j]<=0.) break; else continue; 
      } else if (key[i]==NEG) { if (vvo[i][j]>=0.) break; else continue; 
      } else if (key[i]!=vvo[i][j]) break; // default
    }
    if (i==nk) ind[k++]=j; // all equal
  }
  vector_resize(vv, k);
  return k;
}
ENDVERBATIM

:* v1.flipbits(scratch,num) flips num bits
: uses scratch vector of same size as v1 to make sure doesn't flip same bit twice
VERBATIM
static double flipbits(void* vv) {	
	int i, nx, ny, flip, ii;
	double* x, *y;

	nx = vector_instance_px(vv, &x);
	ny = vector_arg_px(1, &y);
        flip = (int)*getarg(2);
	if (nx != ny) {
	  hoc_execerror("Scratch vector must be same size", 0);
	}
	for (i=0; i<nx; i++) { y[i]=x[i]; } /* copy */
	for (i=0; i < flip; i++) { /* flip these bits */
	  ii = (int) ((nx+1)*drand48());
	  if (x[ii]==y[ii]) { /* hasn't been touched */
	    x[ii]=((x[ii]==1.)?BVBASE:1.);
	  } else {
	    i--; /* do it again */
	  }
	}
	return flip;
}
ENDVERBATIM
 
:* v1.flipbalbits(scratch,num) flips num bits making sure to balance every 1
: flip with a 0 flip to preserve initial power
: uses scratch vector of same size as v1 to make sure doesn't flip same bit twice
VERBATIM
static double flipbalbits(void* vv) {	
	int i, nx, ny, flip, ii, next;
	double* x, *y;

	nx = vector_instance_px(vv, &x);
	ny = vector_arg_px(1, &y);
        flip = (int)*getarg(2);
	if (nx != ny) {
	  hoc_execerror("Scratch vector must be same size", 0);
	}
	for (i=0; i<nx; i++) { y[i]=x[i]; } /* copy */
        next = 1; /* start with 1 */
	for (i=0; i < flip;) { /* flip these bits */
	  ii = (int) ((nx+1)*drand48());
	  if (x[ii]==y[ii] && y[ii]==next) { /* hasn't been touched */
	    next=x[ii]=((x[ii]==1.)?BVBASE:1.);
            i++;
	  }
	}
	return flip;
}
ENDVERBATIM
 
:* v1.vpr() prints out neatly
VERBATIM
static double vpr(void* vv) {
  int i, nx;
  double* x;
  FILE* f;
  nx = vector_instance_px(vv, &x);
  if (ifarg(1)) { 
    f = hoc_obj_file_arg(1);
    for (i=0; i<nx; i++) {
      if (x[i]>BVBASE) { fprintf(f,"%d",1); 
      } else { fprintf(f,"%d",0); }
    }
    fprintf(f,"\n");
  } else {
    for (i=0; i<nx; i++) {
      if (x[i]>BVBASE) { printf("%d",1); 
      } else { printf("%d",0); }
    }
    printf("\n");
  }
  return 1.;
}
ENDVERBATIM
 
:* v1.thresh() threshold above and below thresh
VERBATIM
static double thresh(void* vv) {
  int i, nx, ny, cnt;
  double *x, *y, th;
  nx = vector_instance_px(vv, &x);
  cnt=0;
  if (hoc_is_object_arg(1)) { 
    ny = vector_arg_px(1, &y); th=0;
    if (nx!=ny) { hoc_execerror("Vector sizes don't match in thresh.", 0); }
    for (i=0; i<nx; i++) { if (x[i]>=y[i]) { x[i]= 1.; cnt++;} else { x[i]= BVBASE; } }
  } else { th = *getarg(1);
    for (i=0; i<nx; i++) { if (x[i] >= th) { x[i]= 1.; cnt++;} else { x[i]= BVBASE; } }
  }
  return cnt;
}
ENDVERBATIM
 
:* v1.triplet() return location of a triplet
VERBATIM
static double triplet(void* vv) {
  int i, nx;
  double *x, *y, a, b;
  nx = vector_instance_px(vv, &x);
  a = *getarg(1); b = *getarg(2);
  for (i=0; i<nx; i+=3) if (x[i]==a&&x[i+1]==b) break;
  if (i<nx) return (double)i; else return -1.;
}
ENDVERBATIM

:* state.onoff(volts,OBon,thresh,dur,refr)
:  looks at volts vector to decide if have rearch threshold thresh
:  OBon takes account of burst dur and refractory period
VERBATIM
static double onoff(void* vv) {
  int i, j, n, nv, non, nt, nd, nr, num;
  double *st, *vol, *obon, *thr, *dur, *refr;
  n = vector_instance_px(vv, &st);
  nv   = vector_arg_px(1, &vol);
  non  = vector_arg_px(2, &obon);
  nt   = vector_arg_px(3, &thr);
  nd   = vector_arg_px(4, &dur);
  nr   = vector_arg_px(5, &refr);
  if (n!=nv||n!=non||n!=nt||n!=nd||n!=nr) {
    hoc_execerror("v.onoff: vectors not all same size", 0); }
  for (i=0,num=0;i<n;i++) {
    obon[i]--;
    if (obon[i]>0.) { st[i]=1.; continue; } /* cell must fire */
    if (vol[i]>=thr[i] && obon[i]<= -refr[i]) { /* past refractory period */
        st[i]=1.; obon[i]=dur[i]; num++;
    } else { st[i]= BVBASE; }
  }
  return num;
}
ENDVERBATIM
 
:* vo.bpeval(outp,del)
:  service routine for back-prop: vo=outp*(1-outp)*del
VERBATIM
static double bpeval(void* vv) {
  int i, n, no, nd, flag=0;
  double add,div;
  double *vo, *outp, *del;
  n = vector_instance_px(vv, &vo);
  no   = vector_arg_px(1, &outp);
  nd   = vector_arg_px(2, &del);
  if (ifarg(3) && ifarg(4)) { add= *getarg(3); div= *getarg(4); flag=1;}
  if (n!=no||n!=nd) hoc_execerror("v.bpeval: vectors not all same size", 0);
  if (flag) {
    for (i=0;i<n;i++) vo[i]=((outp[i]+add)/div)*(1.-1.*((outp[i]+add)/div))*del[i];
  } else {
    for (i=0;i<n;i++) vo[i]=outp[i]*(1.-1.*outp[i])*del[i];
  }
}
ENDVERBATIM
 
:* v1.sedit() string edit
VERBATIM
static double sedit(void* vv) {
  int i, n, ni, f=0;
  double *x, *ind, th, val;
  Symbol* s; char *op;
  op = gargstr(1);
  n = vector_instance_px(vv, &x);
  sprintf(op,"hello world");
  return (double)n;
}
ENDVERBATIM

:* v1.w() a .where that sets elements in source vector
VERBATIM
static double w(void* vv) {
  int i, n, ni, f=0;
  double *x, *ind, th, val;
  Symbol* s; char *op;
  if (! ifarg(1)) { 
    printf("v1.w('op',thresh[,val,v2])\n"); 
    printf("  a .where that sets elements in v1 to val (default 0), if v2 => only look at these elements\n");
    printf("  'op'=='function name' is a .apply targeted by v2 called as func(x[i],thresh,val)\n");
    return -1.;
  }
  op = gargstr(1);
  n = vector_instance_px(vv, &x);
  th = *getarg(2);
  if (ifarg(3)) { val = *getarg(3); } else { val = 0.0; }
  if (ifarg(4)) {ni = vector_arg_px(4, &ind); f=1;} // just look at the spots indexed
  if (!strcmp(op,"==")) { 
    if (f==1) {for (i=0; i<ni;i++) {if (x[(int)ind[i]]==th) { x[(int)ind[i]] = val;}}
    } else {for (i=0; i<n; i++) {if (x[i]==th) { x[i] = val;}}}
  } else if (!strcmp(op,"!=")) {
    if (f==1) {for (i=0; i<ni;i++) {if (x[(int)ind[i]]!=th) { x[(int)ind[i]] = val;}}
    } else {for (i=0; i<n; i++) {if (x[i]!=th) { x[i] = val;}}}
  } else if (!strcmp(op,">")) {
    if (f==1) {for (i=0; i<ni;i++) {if (x[(int)ind[i]]>th) { x[(int)ind[i]] = val;}}
    } else {for (i=0; i<n; i++) {if (x[i]>th) { x[i] = val;}}}
  } else if (!strcmp(op,"<")) {
    if (f==1) {for (i=0; i<ni;i++) {if (x[(int)ind[i]]<th) { x[(int)ind[i]] = val;}}
    } else {for (i=0; i<n; i++) {if (x[i]<th) { x[i] = val;}}}
  } else if (!strcmp(op,">=")) {
    if (f==1) {for (i=0; i<ni;i++) {if (x[(int)ind[i]]>=th) { x[(int)ind[i]] = val;}}
    } else {for (i=0; i<n; i++) {if (x[i]>=th) { x[i] = val;}}}
  } else if (!strcmp(op,"<=")) {
    if (f==1) {for (i=0; i<ni;i++) {if (x[(int)ind[i]]<=th) { x[(int)ind[i]] = val;}}
    } else {for (i=0; i<n; i++) {if (x[i]<=th) { x[i] = val;}}}
  } else if ((s=hoc_lookup(op))) { // same as .apply but only does indexed ones
    if (f==1) {for (i=0; i<ni;i++) {
        hoc_pushx(x[(int)ind[i]]); hoc_pushx(th); hoc_pushx(val);
      x[(int)ind[i]]=hoc_call_func(s, 3);}
    } else {for (i=0; i<n;i++) {hoc_pushx(x[i]); hoc_pushx(th); hoc_pushx(val);
        x[i]=hoc_call_func(s, 3);}}
  }
  return (double)i;
}
ENDVERBATIM
  
:* dest.xing(src,tvec,thresh) 
: a .where that looks for threshold crossings and then doesn't pick another till
: comes back down again; places values from tvec in dest; interpolates
VERBATIM
static double xing(void* vv) {
  int i, j, nsrc, ndest, ntvec, f, maxsz;
  double *src, *dest, *tvec, th;
  ndest = vector_instance_px(vv, &dest);
  nsrc = vector_arg_px(1, &src);
  ntvec = vector_arg_px(2, &tvec);
  th = *getarg(3);
  maxsz=vector_buffer_size(vv);
  vector_resize(vv, maxsz);
  if (nsrc!=ntvec) hoc_execerror("v.xing: vectors not all same size", 0);
  for (i=0,f=0,j=0; i<nsrc; i++) {
    if (src[i]>th) { /* ? passing thresh */
      if (f==0) { 
        if (j>=maxsz) {
          printf("(%d) :: ",maxsz);
          hoc_execerror("Dest vec too small in xing ", 0);
        }
        if (i>0) { /* don't record if first value is above thresh */
          dest[j++] = tvec[i-1] + (tvec[i]-tvec[i-1])*(th-src[i-1])/(src[i]-src[i-1]);
        }
        f=1; 
      }
    } else {       /* below thresh */
      if (f==1) { f=0; } /* just passed going down */
    }
  }
  vector_resize(vv, j);
  return (double)i;
}
ENDVERBATIM

:* v1.xzero() looks for zero crossings 
VERBATIM
static double xzero(void* vv) {
  int i, n, nv, up;
  double *x, *vc, th, cnt=0.;
  n = vector_instance_px(vv, &x);
  nv = vector_arg_px(1, &vc);
  if (ifarg(2)) { th = *getarg(2); } else { th=0.0;}
  if (vc[0]<th) up=0; else up=1;  /* F or T */
  if (nv!=n) hoc_execerror("Vector size doesn't match.", 0);
  for (i=0; i<nv; i++) {
    x[i]=0.;
    if (up) { /* look for passing down */
      if (vc[i]<th) { x[i]=-1; up=0; cnt++; }
    } else if (vc[i]>th) { up=x[i]=1; cnt++; }
  }
  return cnt;
}
ENDVERBATIM
  
:* v1.sw(FROM,TO) switchs all FROMs to TO
VERBATIM
static double sw(void* vv) {
  int i, n;
  double *x, fr, to;
  n = vector_instance_px(vv, &x);
  fr = *getarg(1);
  to = *getarg(2);
  for (i=0; i<n; i++) {
    if (x[i]==fr) { x[i] = to;}
  }
  return n;
}
ENDVERBATIM
 
:* v.b2v(bytevec) copies from vector to bytevec
VERBATIM
static double b2v(void* vv) {
  int i, n, num;
  double *x; bvec* to; Object *ob;
  n = vector_instance_px(vv, &x);
  ob = *(hoc_objgetarg(1));
  to = (bvec*)ob->u.this_pointer; // doesn't check that this is actually a bvec
  if (to->size!=n) { hoc_execerror("Vector and bytevec sizes don't match.", 0); }
  for (i=0; i<n; i++) x[i] = (double)to->x[i];
  return n;
}
ENDVERBATIM

:* v.v2d(&x) copies from vector to double area -- a seg error waiting to happen
VERBATIM
static double v2d(void* vv) {
  int i, n, num;
  double *x, *to;
  n = vector_instance_px(vv, &x);
  to = hoc_pgetarg(1);
  if (ifarg(2)) { num = *getarg(2); } else { num=-1;}
  if (num>-1 && num!=n) { hoc_execerror("Vector size doesn't match.", 0); }
  for (i=0; i<n; i++) {to[i] = x[i];}
  return n;
}
ENDVERBATIM
 
:* v.d2v(&x) copies from double area to vector -- a seg error waiting to happen
VERBATIM
static double d2v(void* vv) {
  int i, n, num;
  double *x, *fr;
  n = vector_instance_px(vv, &x);
  fr = hoc_pgetarg(1);
  if (ifarg(2)) { num = *getarg(2); } else { num=-1;}
  if (num>-1 && num!=n) { hoc_execerror("Vector size doesn't match.", 0); }
  for (i=0; i<n; i++) {x[i] = fr[i];}
  return n;
}
ENDVERBATIM
 
:* v.lcat(LIST)
VERBATIM
static double lcat(void* vv) {
  int i, j, k, n, lc, cap, maxsz;
  Object *ob1, *ob2;
  double *x, *fr; void *vw;
  n = vector_instance_px(vv, &x);
  vector_resize(vv,maxsz=vector_buffer_size(vv)); // open it up fully
  ob1 = *hoc_objgetarg(1);
  lc = ivoc_list_count(ob1);
  for (i=0,j=0;i<lc && j<maxsz;i++) {
    cap = list_vector_px2(ob1, i, &fr, &vw);
    for (k=0;k<cap && j<maxsz;k++,j++) x[j]=fr[k];
  }
  if (i<lc || k<cap) printf("vecst lcat WARN: not all vecs copied\n");
  vector_resize(vv,j);
  return (double)j;
}
ENDVERBATIM
 
VERBATIM
//* list_vector_px(LIST,ITEM#,DOUBLE PTR ADDRESS) 
// modeled on vector_arg_px() picks up a vec from a list
static int list_vector_px(Object *ob, int i, double** px) {
  Object* obv;
  int sz;
  obv = ivoc_list_item(ob, i);
  sz = vector_capacity(obv->u.this_pointer);
  *px = vector_vec(obv->u.this_pointer);
  return sz;
}

//* list_vector_px(LIST,ITEM#,DOUBLE PTR ADDRESS,VEC POINTER ADDRESS) 
//  returns the vector pointer as well as the double pointer
static int list_vector_px2 (Object *ob, int i, double** px, void** vv) {
  Object* obv;
  int sz;
  obv = ivoc_list_item(ob, i);
  sz = vector_capacity(obv->u.this_pointer);
  *px = vector_vec(obv->u.this_pointer);
  *vv = (void*) obv->u.this_pointer;
  return sz;
}

//* list_vector_resize(LIST,ITEM#,NEW SIZE)
static int list_vector_resize (Object *ob, int i, int sz) {
  Object* obv;
  int maxsz;
  obv = ivoc_list_item(ob, i);
  maxsz = vector_buffer_size(obv->u.this_pointer);
  if (sz>maxsz) {
    printf("max:%d request:%d ",maxsz,sz);
    hoc_execerror("Can't grow vector in list_vector_resize ", 0);
    return -1;
  }
  vector_resize(obv->u.this_pointer,sz);
  return sz;
}
ENDVERBATIM

:* v1.ismono([arg]) asks whether is monotonically increasing, with arg==-1 - decreasing
:  with arg==0 - all same
VERBATIM
static double ismono1(double *x, double n, int flag) {
  int i; double last;
  last=x[0]; 
  if (flag==1) {
    for (i=1; i<n && x[i]>=last; i++) last=x[i];
  } else if (flag==-1) {
    for (i=1; i<n && x[i]<=last; i++) last=x[i];
  } else if (flag==0) {
    for (i=1; i<n && x[i]==last; i++) ;
  }
  if (i==n) return 1.; else return 0.;
}

static double ismono(void* vv) {
  int i, n, flag;
  double *x,last;
  n = vector_instance_px(vv, &x);
  if (ifarg(1)) { flag = (int)*getarg(1); } else { flag = 1; }
  return ismono1(x,n,flag);
}
ENDVERBATIM
 
:* v1.count(num) returns number of instances of num
VERBATIM
static double count(void* vv) {
  int i, n, cnt=0;
  double *x,num;
  n = vector_instance_px(vv, &x);
  num = *getarg(1);
  for (i=0; i<n; i++) if (x[i]==num) cnt++;
  return cnt;
}
ENDVERBATIM

:* fac (n) 
: from numerical recipes p.214
FUNCTION fac (n) {
VERBATIM {
    static int ntop=4;
    static double a[101]={1.,1.,2.,6.,24.};
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
      -1.231739516,0.120858003e-2,-0.536382e-5};
    int j,n;
    n = (int)_ln;
    if (n<0) { hoc_execerror("No negative numbers ", 0); }
    if (n>100) { /* gamma function */
      double x,tmp,ser;
      x = _ln;
      tmp=x+5.5;
      tmp -= (x+0.5)*log(tmp);
      ser=1.0;
      for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
      }
      return exp(-tmp+log(2.50662827465*ser));
    } else {
      while (ntop<n) {
        j=ntop++;
        a[ntop]=a[j]*ntop;
      }
    return a[n];
    }
}
ENDVERBATIM
}
 
:* logfac (n)
: from numerical recipes p.214
FUNCTION logfac (n) {
VERBATIM {
    static int ntop=4;
    static double a[101]={1.,1.,2.,6.,24.};
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
      -1.231739516,0.120858003e-2,-0.536382e-5};
    int j,n;
    n = (int)_ln;
    if (n<0) { hoc_execerror("No negative numbers ", 0); }
    if (n>100) { /* gamma function */
      double x,tmp,ser;
      x = _ln;
      tmp=x+5.5;
      tmp -= (x+0.5)*log(tmp);
      ser=1.0;
      for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
      }
      return (-tmp+log(2.50662827465*ser));
    } else {
      while (ntop<n) {
        j=ntop++;
        a[ntop]=a[j]*ntop;
      }
    return log(a[n]);
    }
}
ENDVERBATIM
}

PROCEDURE Expo (x) {
  TABLE RES FROM -20 TO 0 WITH 5000
  RES = exp(x)
}

FUNCTION AAA (x) {
  Expo(x)
  AAA = RES
}

:* dest.smgs(src,low,high,step,var)
:  rewrite on v.sumgauss() in nrn5.3::ivoc/ivocvect.cpp:1078
VERBATIM
static double smgs (void* vv) {	
  int i, j, nx, xv, nsum, points, maxsz;
  double *x, *sum;
  double  low , high , step , var , svar , scale , arg;

  nsum = vector_instance_px(vv, &sum);
  nx = vector_arg_px(1,&x);
  low = *getarg(2);
  high = *getarg(3);
  step = *getarg(4);
  var = *getarg(5);

  points = (int)((high-low)/step+hoc_epsilon);
  if (nsum!=points) { 
    maxsz=vector_buffer_size(vv);
    if (points<=maxsz) {
      vector_resize(vv, nsum); nsum=points; 
    } else {
      printf("%d > %d :: ",points,maxsz);
      hoc_execerror("Vector max capacity too small in smgs ", 0);
    }
  }

  svar = -2.*var*var/step/step;
  scale = 1./sqrt(2.*M_PI)/var;

  for (j=0; j<points;j++) sum[j] = 0.;
  for (i=0;i<nx;i++) {
    xv = (int)((x[i]-low)/step + 0.5);
    for (j=xv; j<points && (arg=(j-xv)*(j-xv)/svar)>-20;j++) {
      Expo(arg);
      sum[j] += RES;
    }
    for (j=xv-1; j>=0 && (arg=(j-xv)*(j-xv)/svar)>-20;j--) {
      Expo(arg);
      sum[j] += RES;
    }
  }
  for (j=0; j<points;j++) sum[j] *= scale;
  return svar;
}
ENDVERBATIM

:* PROCEDURE install_vecst()
PROCEDURE install_vecst () {
  VECST_INSTALLED=1
VERBATIM
  install_vector_method("slope", slope);
  install_vector_method("vslope", vslope);
  install_vector_method("stats", stats);
  install_vector_method("vstats", vstats);
  install_vector_method("indset", indset);
  install_vector_method("randwd", randwd);
  install_vector_method("hamming", hamming);
  install_vector_method("flipbits", flipbits);
  install_vector_method("flipbalbits", flipbalbits);
  install_vector_method("vpr", vpr);
  install_vector_method("thresh", thresh);
  install_vector_method("triplet", triplet);
  install_vector_method("onoff", onoff);
  install_vector_method("bpeval", bpeval);
  install_vector_method("w", w);
  install_vector_method("sedit", sedit);
  install_vector_method("xing", xing);
  install_vector_method("cvlv", cvlv);
  install_vector_method("intrp", intrp);
  install_vector_method("xzero", xzero);
  install_vector_method("sw", sw);
  install_vector_method("ismono", ismono);
  install_vector_method("count", count);
  install_vector_method("fewind", fewind);
  install_vector_method("findx", findx);
  install_vector_method("sindx", sindx);
  install_vector_method("sindv", sindv);
  install_vector_method("nind", nind);
  install_vector_method("keyind", keyind);
  install_vector_method("slct", slct);
  install_vector_method("slor", slor);
  install_vector_method("insct", insct);
  install_vector_method("cull", cull);
  install_vector_method("redundout", redundout);
  install_vector_method("d2v", d2v);
  install_vector_method("v2d", v2d);
  install_vector_method("b2v", b2v);
  install_vector_method("iwr", iwr);
  install_vector_method("ird", ird);
  install_vector_method("smgs", smgs);
  install_vector_method("info", info);
  install_vector_method("lcat", lcat);
ENDVERBATIM
}

: unable to get the drand here to recognize the same fseed used in rand
PROCEDURE vseed (seed) {
  VERBATIM
  srand48((unsigned)_lseed);
  set_seed(_lseed);
  ENDVERBATIM
}

FUNCTION isojt () {
  VERBATIM
  Object *ob1, *ob2;
  ob1 = *hoc_objgetarg(1); ob2 = *hoc_objgetarg(2);
  if (!ob1) if (!ob2) return 1; else return 0;
  if (!ob2 || ob1->template != ob2->template) {
    return 0;
  }
  return 1;
  ENDVERBATIM
}

FUNCTION eqojt () {
  VERBATIM
  Object *ob1, *ob2;
  ob1 = *hoc_objgetarg(1); ob2 = *hoc_objgetarg(2);
  if (ob1 && ob2 && ob1==ob2) {
    return 1;
  }
  return 0;
  ENDVERBATIM
}
