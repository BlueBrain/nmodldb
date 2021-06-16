/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init_
#define _nrn_initial _nrn_initial_
#define nrn_cur _nrn_cur_
#define _nrn_current _nrn_current_
#define nrn_jacob _nrn_jacob_
#define nrn_state _nrn_state_
#define _net_receive _net_receive_ 
#define install_vecst install_vecst_ 
#define vseed vseed_ 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_fac(void);
 static void _hoc_install_vecst(void);
 static void _hoc_logfac(void);
 static void _hoc_vseed(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_vecst", _hoc_setdata,
 "fac", _hoc_fac,
 "install_vecst", _hoc_install_vecst,
 "logfac", _hoc_logfac,
 "vseed", _hoc_vseed,
 0, 0
};
 extern double fac( double );
 extern double logfac( double );
 /* declare global and static user variables */
 double BVBASE = -1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 0,0
};
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "BVBASE", &BVBASE,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"vecst",
 0,
 0,
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 0, _prop);
 	/*initialize range parameters*/
 	_prop->param = _p;
 	_prop->param_size = 0;
 
}
 static void _initlists();
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _vecst_reg() {
	int _vectorized = 0;
  _initlists();
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 vecst /Users/awile/workbench/modldb/nmodldb/models/db/modeldb/7399/mod/vecst.mod\n");
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int install_vecst();
static int vseed(double);
 
/*VERBATIM*/
#include <stdlib.h>
#include <math.h>
/* actually doesn't contain MAXLONG #include <float.h> /* contains MAXLONG */
#include <limits.h> // includes LONG_MAX instead of MAXLONG which is obsolete
extern double* hoc_pgetarg();
extern double hoc_call_func(Symbol*, int narg);
extern FILE* hoc_obj_file_arg(int narg);
extern void vector_resize();
extern int vector_instance_px();
extern void* vector_arg();

/* some machines do not have drand48 and srand48 so use the implementation
at the end of this file */
extern double my_drand48();
extern void my_srand48();
#undef drand48
#undef srand48
#define drand48 my_drand48
#define srand48 my_srand48
 
/*VERBATIM*/
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

#if defined(__MRC__)
/*
apparently needed by something called from sns.inc which is not in
the mac library distributed with the NEURON installation
*/
unsigned int __uitrunc(x) double x; { return (unsigned)(int)(x);}
#endif

 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
/* Maintain parallel int vector to avoid slowness of repeated casts */
static int scrsz=0;
static int *scr;
 
/*VERBATIM*/
static double fewind(void* vv) {
	int i, j, nx, ni, nv[10], num;
	double *x, *ind, *vvo[10];
	nx = vector_instance_px(vv, &x);
        for (i=0;ifarg(i);i++);
	if (i>9) hoc_execerror("ERR: fewind can only handle 9 vectors", 0);
	num = i-2; /* number of vectors to be picked apart */
        for (i=0;i<num;i++) { 
          nv[i] = vector_arg_px(i+2, &vvo[i]);
          if (nx!=nv[i]) { printf("fewind ERR %d %d %d\n",i,nx,nv[i]);
            hoc_execerror("Vectors must all be same size: ", 0); }
        }
        ni = vector_arg_px(1, &ind);
        if (ni>scrsz) { 
          if (scrsz>0) { free(scr); scr=(int *)NULL; }
          scrsz=ni+10000;
          scr=(int *)ecalloc(scrsz, sizeof(int));
        }
        for (i=0;i<ni;i++) scr[i]=(int)ind[i]; /* copy into integer array */
        for (j=0;j<num;j++) {
          for (i=0;i<ni;i++) x[i]=vvo[j][scr[i]];    
          for (i=0;i<ni;i++) vvo[j][i]=x[i];   
          vv=vector_arg(j+2); vector_resize(vv, ni);
        }
	return ni;
}
 
/*VERBATIM*/
static double insct(void* vv) {
	int i, j, k, nx, nv1, nv2;
	double *x, *v1, *v2;
	nx = vector_instance_px(vv, &x);
	nv1 = vector_arg_px(1, &v1);
	nv2 = vector_arg_px(2, &v2);
        for (i=0,k=0;i<nv1&&k<nx;i++) for (j=0;j<nv2&&k<nx;j++) if (v1[i]==v2[j]) x[k++]=v1[i];
        if (k==nx) { 
          printf("\tinsct WARNING: ran out of room: %d\n",k);
          for (;i<nv1;i++,j=0) for (;j<nv2;j++) if (v1[i]==v2[j]) k++;
        } else { vector_resize(vv, k); } /* can't resize to make bigger */
	return (double)k;
}
 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
static double triplet(void* vv) {
  int i, nx;
  double *x, *y, a, b;
  nx = vector_instance_px(vv, &x);
  a = *getarg(1); b = *getarg(2);
  for (i=0; i<nx; i+=3) if (x[i]==a&&x[i+1]==b) break;
  if (i<nx) return (double)i; else return -1.;
}
 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
static double sedit(void* vv) {
  int i, n, ni, f=0;
  double *x, *ind, th, val;
  Symbol* s; char *op;
  op = gargstr(1);
  n = vector_instance_px(vv, &x);
#if defined(MSWIN)
/* not available prior to version 5.3 and not used in this simulation */
#else
  sprintf(op,"hello world");
#endif
  return (double)n;
}
 
/*VERBATIM*/
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
 
/*VERBATIM*/
static double xing(void* vv) {
  int i, n, f=0;
  double *x, th;
  n = vector_instance_px(vv, &x);
  th = *getarg(1);
  for (i=0; i<n; i++) {
    if (x[i]>th) { /* ? passing thresh */
      if (f==0) { x[i]=1.;  f=1; } else { x[i]=0.; } 
    } else {       /* below thresh */
      x[i]=0.;
      if (f==1) { f=0; } /* just passed going down */
    }
  }
  return (double)i;
}
 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
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
 
/*VERBATIM*/
static double ismono(void* vv) {
  int i, n, flag;
  double *x,last;
  n = vector_instance_px(vv, &x);
  if (ifarg(1)) { flag = (int)*getarg(1); } else { flag = 1; }
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
 
/*VERBATIM*/
static double count(void* vv) {
  int i, n, cnt=0;
  double *x,num;
  n = vector_instance_px(vv, &x);
  num = *getarg(1);
  for (i=0; i<n; i++) if (x[i]==num) cnt++;
  return cnt;
}
 
double fac (  double _ln ) {
   double _lfac;
 
/*VERBATIM*/
{
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
 
return _lfac;
 }
 
static void _hoc_fac(void) {
  double _r;
   _r =  fac (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double logfac (  double _ln ) {
   double _llogfac;
 
/*VERBATIM*/
{
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
 
return _llogfac;
 }
 
static void _hoc_logfac(void) {
  double _r;
   _r =  logfac (  *getarg(1) );
 hoc_retpushx(_r);
}
 
/*VERBATIM*/
static double frd(void* vv) {	
	int i, nx, num, size;
        char *xf;
        FILE* f;
	double* x;

	nx = vector_instance_px(vv, &x);
#if defined(MSWIN)
/* not available prior to version 5.3 and not used in this simulation */
#else
        f = hoc_obj_file_arg(1);
#endif
        num = (int)*getarg(2);
        size = (int)*getarg(3);
        xf = (char *)malloc(num * size);
        if (num!=nx) { hoc_execerror("Correct vector size ", 0); }
        fread(xf, num, size, f);
        for (i=0;i<num;++i) x[i] = (double)xf[i];
	return;
}
 
static int  install_vecst (  ) {
   
/*VERBATIM*/
  install_vector_method("slope", slope);
  install_vector_method("vslope", vslope);
  install_vector_method("stats", stats);
  install_vector_method("vstats", vstats);
  install_vector_method("indset", indset);
  install_vector_method("randwd", randwd);
  install_vector_method("hamming", hamming);
  install_vector_method("flipbits", flipbits);
  install_vector_method("frd", frd);
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
  install_vector_method("nind", nind);
  install_vector_method("insct", insct);
  install_vector_method("d2v", d2v);
  install_vector_method("v2d", v2d);
  return 0; }
 
static void _hoc_install_vecst(void) {
  double _r;
   _r = 1.;
 install_vecst (  );
 hoc_retpushx(_r);
}
 
static int  vseed (  double _lseed ) {
   
/*VERBATIM*/
  srand48((unsigned)_lseed);
  return 0; }
 
static void _hoc_vseed(void) {
  double _r;
   _r = 1.;
 vseed (  *getarg(1) );
 hoc_retpushx(_r);
}
 
/*VERBATIM*/
/* http://www.mit.edu/afs/athena/activity/c/cgs/src/math/drand48/ */
/*
 Michael Hines removed  all code not used by srand48 and drand48,
 the code handling non-floating point processor machines, and the
 pdp-11 fragment. Global names have my_ prefix added.
*/


/*	@(#)drand48.c	2.2	*/
/*LINTLIBRARY*/
/*
 *	drand48, etc. pseudo-random number generator
 *	This implementation assumes unsigned short integers of at least
 *	16 bits, long integers of at least 32 bits, and ignores
 *	overflows on adding or multiplying two unsigned integers.
 *	Two's-complement representation is assumed in a few places.
 *	Some extra masking is done if unsigneds are exactly 16 bits
 *	or longs are exactly 32 bits, but so what?
 *	An assembly-language implementation would run significantly faster.
 */
#define N	16
#define MASK	((unsigned)(1 << (N - 1)) + (1 << (N - 1)) - 1)
#define LOW(x)	((unsigned)(x) & MASK)
#define HIGH(x)	LOW((x) >> N)
#define MUL(x, y, z)	{ long l = (long)(x) * (long)(y); \
		(z)[0] = LOW(l); (z)[1] = HIGH(l); }
#define CARRY(x, y)	((long)(x) + (long)(y) > MASK)
#define ADDEQU(x, y, z)	(z = CARRY(x, (y)), x = LOW(x + (y)))
#define X0	0x330E
#define X1	0xABCD
#define X2	0x1234
#define A0	0xE66D
#define A1	0xDEEC
#define A2	0x5
#define C	0xB
#define SET3(x, x0, x1, x2)	((x)[0] = (x0), (x)[1] = (x1), (x)[2] = (x2))
#define SEED(x0, x1, x2) (SET3(x, x0, x1, x2), SET3(a, A0, A1, A2), c = C)

static unsigned x[3] = { X0, X1, X2 }, a[3] = { A0, A1, A2 }, c = C;
static unsigned short lastx[3];
static void next();

double
my_drand48()
{
	static double two16m = 1.0 / (1L << N);

	next();
	return (two16m * (two16m * (two16m * x[0] + x[1]) + x[2]));
}

static void
next()
{
	unsigned p[2], q[2], r[2], carry0, carry1;

	MUL(a[0], x[0], p);
	ADDEQU(p[0], c, carry0);
	ADDEQU(p[1], carry0, carry1);
	MUL(a[0], x[1], q);
	ADDEQU(p[1], q[0], carry0);
	MUL(a[1], x[0], r);
	x[2] = LOW(carry0 + carry1 + CARRY(p[1], r[0]) + q[1] + r[1] +
		a[0] * x[2] + a[1] * x[1] + a[2] * x[0]);
	x[1] = LOW(p[1] + r[0]);
	x[0] = LOW(p[0]);
}

void
my_srand48(seedval)
long seedval;
{
	SEED(X0, LOW(seedval), HIGH(seedval));
}

#if 0
#ifdef DRIVER
/*
	This should print the sequences of integers in Tables 2
		and 1 of the TM:
	1623, 3442, 1447, 1829, 1305, ...
	657EB7255101, D72A0C966378, 5A743C062A23, ...
 */
#include <stdio.h>

main()
{
	int i;

	for (i = 0; i < 80; i++) {
		printf("%4d ", (int)(4096 * my_drand48()));
		printf("%.4X%.4X%.4X\n", x[2], x[1], x[0]);
	}
}
#endif
#endif

static void initmodel() {
  int _i; double _save;_ninits++;
{

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}
