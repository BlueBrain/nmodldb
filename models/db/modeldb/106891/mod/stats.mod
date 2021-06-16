: $Id: stats.mod,v 1.74 2007/12/01 00:24:14 billl Exp $
 
:* COMMENT
COMMENT
randwd   randomly chooses n bits to set to 1
hamming  v.hamming(v1) is hamming distance between 2 vecs
flipbits v.flipbits(scratch,num) flips num rand chosen bits
flipbalbits v.flipbalbits(scratch,num) balanced flipping
vpr      v.vpr prints out vector as 1 (x[i]>0) or 0 (x[i]<=0)
fac      not vec related - returns factorial
logfac   not vec related - returns log factorial
vseed    set some C level randomizer seeds
slope(num) does a linear regression to find the slope, assuming num=timestep of vector
vslope(v2) does a linear regression to find the slope, assuming num=timestep of vector
stats(num,[out]) does a linear regression, assuming num=timestep of vector
vstats(v2,[out]) does a linear regression, using v2 as the x-coords
setrnd(v,flag) does set rand using 1:rand, 2:drand48
v.hash(list)  // make a hash out values in vecs in list
v.unnan([nan_value,][inf_value])  // remove nan's and inf's from a vector
ENDCOMMENT

NEURON {
  SUFFIX stats
  GLOBAL  INSTALLED,seed,kmeasure,verbose
}

PARAMETER {
  : BVBASE = 0. : defined in vecst.mod
  INSTALLED=0
  kmeasure=0
  verbose=0
}

ASSIGNED { seed }

VERBATIM
#include <stdlib.h>
#include <math.h>
#include <limits.h> /* contains LONG_MAX */
#include <float.h>
#include <sys/time.h> 
extern double BVBASE;
extern double* hoc_pgetarg();
Symbol *hoc_get_symbol();
char ** hoc_pgargstr();
extern Point_process* ob2pntproc(Object*);
extern double hoc_call_func(Symbol*, int narg);
extern FILE* hoc_obj_file_arg(int narg);
extern Object** hoc_objgetarg();
extern void vector_resize();
extern double *vector_newsize();
extern int vector_instance_px();
extern void* vector_arg();
extern double* vector_vec();
extern double hoc_epsilon;
extern void mcell_ran4_init(unsigned int *idum);
extern double mcell_ran4(unsigned int* idum,double* ran_vec,unsigned int n,double range);
extern void set_seed();
extern unsigned int *scrset();
extern int ivoc_list_count(Object*);
extern Object* ivoc_list_item(Object*, int);
extern int list_vector_px2();
extern int hoc_is_double_arg(int narg);
extern Objectdata *hoc_objectdata;
extern int openvec(int, double **);
extern char* hoc_object_name(Object*);
extern int nrn_mlh_gsort();
extern int cmpdfn();
int list_vector_px();
int list_vector_resize();
static void hxe() { hoc_execerror("",0); }
unsigned int valseed;

typedef struct BVEC {
 int size;
 int bufsize;
 short *x;
 Object* o;
} bvec;

union dblint {
  int i[2];
  double d;
};

#define VRRY 50
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
        free(x);
        return (n*sigxy - sigx*sigy)/(n*sigx2 - sigx*sigx);
}
ENDVERBATIM
 
:* v1.vslope(v2) does a linear regression, using v2 as the x-coords
VERBATIM
static double vslope (void* vv) {
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

VERBATIM
//computes mean,max squared error of data points
//off a line model with m=slope , b=y_intercept
//x is independent variable
//y is dependent variable
//n is # of data points
//meansqerr is output
//maxsqerr is output
double getsqerr(double* x,double* y,double m,double b,int n,double* meansqerr,double* maxsqerr){
  int i; double val;
  if(!n){
    return -1.0;
  }
  val=0.0;
  *meansqerr=0.0;
  *maxsqerr=0.0;
  for(i=0;i<n;i++){
    val = y[i] - (m*x[i]+b);
    val = val*val;
    if(val>*maxsqerr) *maxsqerr = val;
    *meansqerr += val;
  }
  *meansqerr=*meansqerr/(double)n;
  return *meansqerr;
}
ENDVERBATIM
 
:* v1.stats(num) does a linear regression, assuming num=timestep of vector

VERBATIM
static double stats(void* vv) {
	int i, n;
	double *x, *y, *out;
        double timestep, sigxy, sigx, sigy, sigx2, sigy2;
        double r, m, b, dmeansqerr,dmaxsqerr;
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
        getsqerr(x,y,m,b,n,&dmeansqerr,&dmaxsqerr); //mean,max squared error
        if(ifarg(2)){ //save results to output
          out=vector_newsize(vector_arg(2),5);
          out[0]=m; out[1]=b; out[2]=r; out[3]=dmeansqerr; out[4]=dmaxsqerr;
        } else {
          printf("Examined %d data points\n", n);
          printf("slope     = %f\n", m);
          printf("intercept = %f\n", b);
          printf("R         = %f\n", r);
          printf("R-squared = %f\n", r*r);
          printf("MeanSQErr = %f\n",dmeansqerr);
          printf("MaxSQErr  = %f\n",dmaxsqerr);
        }
        free(x);
        return 1;
}

/* v1.pcorrels2(v2) does a Pearson correlation*/
static double pcorrels2(double *x, double* y, int n) {
  int i;
  double sigxy, sigx, sigy, sigx2, sigy2, dn;
  sigxy=sigx=sigy=sigx2,sigy2=0; // initialize these
  dn=(double)n;
  for(i=0; i<n; i++) {
    sigxy += x[i] * y[i];
    sigx  += x[i];
    sigy  += y[i];
    sigx2 += x[i]*x[i];
    sigy2 += y[i]*y[i];
  }
  sigxy -= (sigx * sigy) / n;
  sigx2 -= (sigx * sigx) / n;
  sigy2 -= (sigy * sigy) / n;
  if(sigx2 <= 0) return 0;
  if(sigy2 <= 0) return 0;
  sigxy = sigxy / sqrt(sigx2*sigy2);
  return sigxy;
}

static double pcorrel(void* vv) {
 int i, n;
  double *x, *y;
  n = vector_instance_px(vv, &x);
  if ((i=vector_arg_px(1, &y)) != n ) {printf("pcorrelsERRA: %d %d\n",n,i); hxe();}
  return pcorrels2(x,y,n);
}
ENDVERBATIM
 
:* vec.vnan() will reset nans and infs to selected values -- default 0,0
VERBATIM
static double unnan(void *vv) {
  int i,nx,cnt; double newnan,newinf;
  union dblint xx;
  double *x;
  newnan=newinf=0;
  nx = vector_instance_px(vv, &x);
  if (ifarg(1)) newinf=newnan=*getarg(1);
  if (ifarg(2)) newinf=*getarg(2);
  for (i=0,cnt=0;i<nx;i++) { 
    xx.d=x[i];
    if (xx.i[0]==0x0 && xx.i[1]==0xfff80000) {x[i]=newnan; cnt++;}
    if (xx.i[0]==0x0 && xx.i[1]==0x7ff00000) {x[i]=newinf; cnt++;}
  }
  return (double)cnt;
}
ENDVERBATIM

:* v1.vstats(v2) does a linear regression, using v2 as the x-coords
VERBATIM
static double vstats(void* vv) {
	int i, n;
	double *x, *y, *out;
        double sigxy, sigx, sigy, sigx2, sigy2;
        double r, m, b, dmeansqerr,dmaxsqerr;
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
          getsqerr(x,y,m,b,n,&dmeansqerr,&dmaxsqerr);//mean,max squared error
          if(ifarg(2)){ //save results to output
            out=vector_newsize(vector_arg(2),5);
            out[0]=m; out[1]=b; out[2]=r; out[3]=dmeansqerr; out[4]=dmaxsqerr;
          } else {
            printf("Examined %d data points\n", n);
            printf("slope     = %f\n", m);
            printf("intercept = %f\n", b);
            printf("R         = %f\n", r);
            printf("R-squared = %f\n", r*r);
            printf("MeanSQErr = %f\n",dmeansqerr);
            printf("MaxSQErr  = %f\n",dmaxsqerr);
          }
          return 1;
        } else {
          printf("You must supply an x vector\n");
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
 
:* v1.hash(veclist)
VERBATIM
static double hash (void* vv) {
  int i, j, nx, nv[VRRY], num, vfl; 
  union dblint xx;
  Object* ob;
  double *x, *vvo[VRRY], big, y, prod;
  nx = vector_instance_px(vv, &x);
  if (ifarg(1)) {
    vfl=0;
    ob=*hoc_objgetarg(1); 
    num = ivoc_list_count(ob);
    if (num>VRRY) {printf("vecst:hash ERR: can only handle %d vecs: %d\n",VRRY,num); hxe();}
    for (i=0;i<num;i++) { nv[i] = list_vector_px(ob, i, &vvo[i]);
      if (nx!=nv[i]) { printf("vecst:hash ERR %d %d %d\n",i,nx,nv[i]);hxe();}}
  } else {
    vfl=1; num=nx; nx=1; // outer loop will go only once
  }
  big=pow(DBL_MAX,1./(double)num); // base biggest # on nth root of num of values being used
  for (i=0;i<nx;i++) {
    for (j=0,prod=1;j<num;j++) {
      if (vfl) {  xx.d=x[j];       // break the double up into 2 unsigned ints
      } else   {  xx.d=vvo[j][i]; }
      if (xx.i[0]==0) { xx.i[0]=xx.i[1]; xx.i[0]<<=4; } // high order bits may be 0
      if (xx.i[1]==0) { xx.i[1]=xx.i[0]; xx.i[1]<<=4; } // low order bits unlikely 0
      mcell_ran4_init(&xx.i[1]);
      mcell_ran4(&xx.i[0], &y, 1, big); // generate a pseudorand number based on these
      prod*=y;  // keep multiplying these out
    }
    if (! vfl) x[i]=prod; else return prod; // just return the 1 value
  }
  return (double)nx;
}
ENDVERBATIM

:* v1.setrnd(flag) performs setrand()
: flag: 1 rand(); 2 drand48(); 3 scop_random(); 4 mcell_ran4(); 5 integers via mcell_ran4()
: v1.setrand(4 or 5[,MAX_VAL DEFAULT=1, SEED])
VERBATIM
static double setrnd (void* vv) {
  int flag, i; unsigned int nx;
  double *x, y, n;
  unsigned long value;
  value=1;
  nx = vector_instance_px(vv, &x);
  flag = (int) *getarg(1);
  if (flag==1) {
    for (i=0; i < nx; i++) x[i] = (double)rand()/RAND_MAX; 
  }  else if (flag==2) {
    for (i=0; i < nx; i++) x[i] = drand48(); 
  } else if (flag==3) { // scop_random()'s cheap and dirty rand
    unsigned long a = 2147437301, c = 453816981, m = ~0;
    value = (unsigned long) seed;
    for (i=0; i < nx; i++) {
      value = a * value + c;
      x[i] = (fabs((double) value / (double) m));
    }
    seed=(double)value;
  } else if (flag==4) { // mcell_ran4
    n=ifarg(2)?(*getarg(2)):1.;
    if (ifarg(3)) valseed=(unsigned int)(*getarg(3));
    mcell_ran4(&valseed, x, nx, n);
    return (double)valseed;
  } else if (flag==5) { // nx integers from 0 to n
    n=ifarg(2)?(*getarg(2)):1.;
    if (ifarg(3)) valseed=(unsigned int)(*getarg(3));
    mcell_ran4(&valseed, x, nx, n);
    for (i=0;i<nx;i++) x[i]=floor(x[i]);
    return (double)valseed;
  }
  return (double)nx;
}
ENDVERBATIM
 
:* v1.hamming(v2[,v3]) compares v1 and v2 for matches, v3 gives diff vector
VERBATIM
static double hamming (void* vv) {
  int i, nx, ny, nz, prflag;
  double* x, *y, *z,sum;
  sum = 0.;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  if (ifarg(2)) { // write a diff vector to z
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

:* v1.distance(v2) euclidean distance
VERBATIM
static double distance (void* vv) {
  int i, nx, ny;
  double* x, *y, sum;
  sum = 0.;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  if (nx!=ny) {printf("Vectors must be same size %d %d\n",nx,ny); hxe();}
  for (i=0; i<nx; i++) sum+=(x[i]-y[i])*(x[i]-y[i]); 
  return sqrt(sum);
}
ENDVERBATIM

:* v1.ndprd(v2) normalized dot prod distance (cos of angle)
VERBATIM
static double ndprd (void* vv) {
  int i, nx, ny;
  double* x, *y, sum, sumx, sumy;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  if (nx!=ny) {printf("Vectors must be same size %d %d\n",nx,ny); hxe();}
  for (i=0, sum=0., sumx=0., sumy=0.; i<nx; i++) {
    sum+=x[i]*y[i]; sumx+=x[i]*x[i]; sumy+=y[i]*y[i]; 
  }
  if (ifarg(2)) { return sum/sqrt(sumx)/sqrt(sumy);                   // cos of angle
  } else {        return acos(sum/sqrt(sumx)/sqrt(sumy))*180./M_PI; } // angle in degrees
}
ENDVERBATIM

:* v1.flipbits(scratch,num) flips num bits
: uses scratch vector of same size as v1 to make sure doesn't flip same bit twice
VERBATIM
static double flipbits (void* vv) {	
  int i, j, nx, ny, flip, ii;
  double *x, *y;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  flip = (int)*getarg(2);
  if (ny<nx) {hoc_execerror("flipbits:Scratch vector must adequate size", 0);}
  mcell_ran4(&valseed, y, (unsigned int)ny, (double)nx); // indices to flip
  for (i=0,j=0; i<flip && j<ny; j++) { // flip these bits
    ii=(int)y[j];
    if        (x[ii]==BVBASE) { x[ii]= 1e9; i++;  // mark location 
    } else if (x[ii]==1)      { x[ii]=-1e9; i++; }
  }
  j=i;
  for (i=0; i<nx; i++) if (x[i]==1e9) x[i]=1; else if (x[i]==-1e9) x[i]=BVBASE;
  return (double)j;
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
 
:* v1.vpr([BASE]) prints out neatly in binary -- optional arg allows base 10,16,64
: generally prints out 1 char per entry
VERBATIM
static double vpr (void* vv) {
  int i, nx, flag;
  double* x; char c;
  FILE* f;
  nx = vector_instance_px(vv, &x);
  if (ifarg(1)) { 
    if (hoc_is_double_arg(1)) {
      flag=(int) *getarg(1);
      if (flag==0) {
        for (i=0; i<nx; i++) printf("%s%d%s",x[i]>=10?"|":"",(int)x[i],x[i]>=10?"|":"");
      } else if (flag==-1) { 
        for (i=0; i<nx; i++) printf("%04.2f|",x[i]);
      } else if (flag==10) { // decimal
        for (i=0; i<nx; i++) if (x[i]>=10) printf("+"); else printf("%d",(int)x[i]);
      } else if (flag==16) { // hex
        for (i=0; i<nx; i++) if (x[i]>=16) printf("+"); else printf("%x",(int)x[i]);
      } else if (flag==64) { // base 64
        for (i=0; i<nx; i++) {
          if (x[i]>=64) {         printf("+"); 
          } else if (x[i]<16) {  printf("%x",(int)x[i]);    // 0-f   0-15
          } else if (x[i]<36) {  printf("%c",(int)x[i]+87); // g-z  16-35
          } else if (x[i]<62) {  printf("%c",(int)x[i]+29); // A-Z  36-61
          } else if (x[i]<63) { printf("@");                // @    62
          } else if (x[i]<64) { printf("=");                // =    63   
          } else printf("ERROR");
        }
      }
      if (!ifarg(2)) printf("\n"); else printf(" ");
    } else {
      f = hoc_obj_file_arg(1);
      for (i=0; i<nx; i++) {
        if (x[i]>BVBASE) { fprintf(f,"%d",1); 
        } else { fprintf(f,"%d",0); }
      }
      fprintf(f,"\n");
    }
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

:* v1.bin(targ,invl[min,max]) place counts for each interval
:* v1.bin(list,invl[min,max]) using NQS("count","index") for easy sorting
: like .hist() but doesn't throw away values <min or >max
: note that optional max denotes the start of an interval not end
VERBATIM
static double bin (void* vv) {	
  int i, j, nx, maxsz, lfl;
  double* x, *y, *ix, invl, min, max, maxf, jj;
  Object* ob;
  void* voi[2];

  min=0; max=1e9; maxf=-1e9;
  nx = vector_instance_px(vv, &x);
  ob =   *hoc_objgetarg(1);
  if (strncmp(hoc_object_name(ob),"Vector",6)==0) { lfl=0; // lfl is list flag
    if ((maxsz=openvec(1, &y))==0) hxe();
    voi[0]=vector_arg(1);
  } else { // list of 2
    lfl=1;
    maxsz = list_vector_px3(ob, 0, &y, &voi[0]);
    if (maxsz!=(i=list_vector_px3(ob,1,&ix,&voi[1]))){printf("binERRA %d %d\n",maxsz,i); hxe();}
  }
  invl = *getarg(2);
  if (ifarg(4)) {min=*getarg(3); max=*getarg(4);
  } else if (ifarg(3)) max=*getarg(3);
  for (j=0; j<maxsz; j++) y[j]=0.;
  for (i=0; i<nx; i++) {
    if (x[i]>=max) {        y[(j=(int)(jj=(max-min)/invl))]++;
    } else if (x[i]<=min) { y[(j=0)]++; jj=0.;
    } else {
      if (x[i]>maxf) maxf=x[i]; // max found
      j=(int)(jj=(x[i]-min)/invl);
      if (j>maxsz-1) {printf("bin ERRB OOB: %d>%d\n",j,maxsz-1); hxe();}
      y[j]++;
    }
    if (lfl) ix[j]=jj+min;
  }
  maxsz=(max==1e9)?(int)(maxf/invl+1):(int)((max-min)/invl+1);
  vector_resize(voi[0], maxsz);
  if (lfl) vector_resize(voi[1], maxsz);
  return (double)maxsz;
}
ENDVERBATIM

:* PROCEDURE install_stats()
PROCEDURE install () {
  if (INSTALLED==1) {
    printf("$Id: stats.mod,v 1.74 2007/12/01 00:24:14 billl Exp $")
  } else {
  INSTALLED=1
VERBATIM
  valseed=1;
  install_vector_method("slope", slope);
  install_vector_method("vslope", vslope);
  install_vector_method("stats", stats);
  install_vector_method("pcorrel", pcorrel);
  install_vector_method("vstats", vstats);
  install_vector_method("randwd", randwd);
  install_vector_method("hamming", hamming);
  install_vector_method("flipbits", flipbits);
  install_vector_method("flipbalbits", flipbalbits);
  install_vector_method("vpr", vpr);
  install_vector_method("bin", bin);
  install_vector_method("setrnd", setrnd);
  install_vector_method("distance", distance);
  install_vector_method("ndprd", ndprd);
  install_vector_method("hash", hash);
  install_vector_method("unnan", unnan);
ENDVERBATIM
  }
}

PROCEDURE prhash (x) {
VERBATIM {
  union dblint xx;
  xx.d=_lx;
  printf("%08x%08x\n",xx.i[0],xx.i[1]);
}
ENDVERBATIM
}

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

FUNCTION getseed () {
  VERBATIM
  seed=(double)valseed;
  return seed;
  ENDVERBATIM
}

: unable to get the drand here to recognize the same fseed used in rand
FUNCTION vseed () {
  VERBATIM
#ifdef WIN32
  if (ifarg(1)) seed=*getarg(1); else {
    printf("TIME ACCESS NOT PRESENT IN WINDOWS\n");
    hxe();
  }
  srand48((unsigned)seed);
  set_seed(seed);
  return seed;
#else
  struct  timeval tp;
  struct  timezone tzp;
  if (ifarg(1)) seed=*getarg(1); else {
    gettimeofday(&tp,&tzp);
    seed=tp.tv_usec;
  }
  srand48((unsigned)seed);
  set_seed(seed);
  srandom(seed);
  valseed=(unsigned int)seed;
  return seed;
#endif
  ENDVERBATIM
}

: mc4seed() set valseed for mccell_ran4() as series of factors so as not to lose low order 
: digits before casting double to unsigned int
FUNCTION mc4seed () {
  VERBATIM
  int i;
  valseed=(unsigned int)(*getarg(1));
  for (i=2;ifarg(i);i++) {
    valseed*=(unsigned int)(*getarg(i));
  }
  return valseed;
  ENDVERBATIM
}

: from Numerical Recipes in C
FUNCTION gammln (xx) {
  VERBATIM {
    double x,tmp,ser;
    static double cof[6]={76.18009173,-86.50532033,24.01409822,-1.231739516,0.120858003e-2,-0.536382e-5};
    int j;
    x=_lxx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.0;
    for (j=0;j<=5;j++) {
      x += 1.0;
      ser += cof[j]/x;
    }
    return -tmp+log(2.50662827465*ser);
  }
  ENDVERBATIM
}

FUNCTION betai(a,b,x) {
VERBATIM {
  double bt;
  double gammln(),betacf();

  if (_lx < 0.0 || _lx > 1.0) {printf("Bad x in routine BETAI\n"); hxe();}
  if (_lx == 0.0 || _lx == 1.0) bt=0.0;
  else
  bt=exp(gammln(_la+_lb)-gammln(_la)-gammln(_lb)+_la*log(_lx)+_lb*log(1.0-_lx));
  if (_lx < (_la+1.0)/(_la+_lb+2.0))
  return bt*betacf(_la,_lb,_lx)/_la;
  else
  return 1.0-bt*betacf(_lb,_la,1.0-_lx)/_lb;
 }
ENDVERBATIM
}

VERBATIM
#define ITMAX 100
#define EPS 3.0e-7
ENDVERBATIM

FUNCTION betacf(a,b,x) {
VERBATIM {
  double qap,qam,qab,em,tem,d;
  double bz,bm=1.0,bp,bpp;
  double az=1.0,am=1.0,ap,app,aold;
  int m;
  void nrerror();

  qab=_la+_lb;
  qap=_la+1.0;
  qam=_la-1.0;
  bz=1.0-qab*_lx/qap;
  for (m=1;m<=ITMAX;m++) {
    em=(double) m;
    tem=em+em;
    d=em*(_lb-em)*_lx/((qam+tem)*(_la+tem));
    ap=az+d*am;
    bp=bz+d*bm;
    d = -(_la+em)*(qab+em)*_lx/((qap+tem)*(_la+tem));
    app=ap+d*az;
    bpp=bp+d*bz;
    aold=az;
    am=ap/bpp;
    bm=bp/bpp;
    az=app/bpp;
    bz=1.0;
    if (fabs(az-aold) < (EPS*fabs(az))) return az;
  }
  printf("a or b too big, or ITMAX too small in BETACF"); return -1.;
}
ENDVERBATIM
}

FUNCTION symval() {
VERBATIM {
  Symbol *sym;
  sym = hoc_get_symbol(* hoc_pgargstr(1));
  // should do type check eg sym->type == VAR
#if defined(t)
  return *(hoc_objectdata[sym->u.oboff]._pval);
#else
  return *(hoc_objectdata[sym->u.oboff].pval);
#endif
 }
ENDVERBATIM
}
