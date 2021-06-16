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
#define install_matrix install_matrix_ 
 
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
 static void _hoc_install_matrix(void);
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
 "setdata_matrix", _hoc_setdata,
 "install_matrix", _hoc_install_matrix,
 0, 0
};
 /* declare global and static user variables */
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
"matrix",
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

 void _matrix_reg() {
	int _vectorized = 0;
  _initlists();
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 matrix /Users/awile/workbench/modldb/nmodldb/models/db/modeldb/7399/mod/matrix.mod\n");
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int install_matrix();
 
/*VERBATIM*/
static double outprod(void* vv) {
  int i, j, nx, ny, nz;
  double *x, *y, *z;
  /* this will be the outer product */
  nx = vector_instance_px(vv, &x);
	
  /* these are the two vectors that make it up */
  ny = vector_arg_px(1, &y); // will be number of columns
  nz = vector_arg_px(2, &z); // will be number of rows
  if (nx != ny*nz) {
    hoc_execerror("Vector size mismatch", 0);
  }
  for (i=0;i<ny;i++) {
    for (j=0;j<nz;j++) {
      x[i*nz+j] = y[i]*z[j];
    }
  }
  return nx;
}
 
/*VERBATIM*/
static double mmult(void* vv) {
  int i, j, nx, ny, nz;
  double *x, *y, *z;
  /* x will be the product of matrix y and vec z */
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  nz = vector_arg_px(2, &z);
  if (ny != nx*nz) {
    hoc_execerror("Vector size mismatch", 0);
  }
  for (i=0;i<nx;i++) {
    x[i] = 0.;
    for (j=0;j<nz;j++) {
      x[i] += y[i*nz+j]*z[j];
    }
  }
  return nx;
}
 
/*VERBATIM*/
static double spltp(void* vv) {
  int ii, jj, nstpr, nstpo, nw, npr, npo, flag, cnt;
  double *stpr, *stpo, *w, *pr, *po;
  extern double hoc_call_func(Symbol*, int narg);

  char func[4] = "ltp";
  Symbol* s = hoc_lookup(func);
  if (! s) { hoc_execerror("Can't find ltp() func", 0); }
  nstpo = vector_instance_px(vv, &stpo);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  nw = vector_arg_px(3, &w);
  nstpr = vector_arg_px(4, &stpr);
  for (ii=0,jj=0,cnt=0;ii<nstpo;ii++) {
    if (stpo[ii]==1.0) { /* connections to these will be changed */ 
      for (;po[jj]<ii;jj++) ; /* move forward till find a po */
      for (;po[jj]==ii;jj++) { /* move through these po's */
	if (stpr[(int)pr[jj]]==1.) { /*  did the presyn spike? */
	  cnt++; hoc_pushx(1.0);
	} else { 
	  cnt--; hoc_pushx(-1.0);
	}
        hoc_pushx(w[jj]);
        w[jj]=hoc_call_func(s, 2);
      }
    }
  }
  return cnt;
}
 
/*VERBATIM*/
/* Maintain a parallel vector of ints to avoid the slowness of repeated casts in spmult */
static int *pr_int;
static int *po_int;
static int cpfl=0;
 
/*VERBATIM*/
static double mkspcp(void* vv) {
  int j, nw, npr, npo;
  double *w, *pr, *po;
  if (! ifarg(1)) { 
    cpfl=0; 
    if (po_int!=NULL) free(po_int); 
    if (pr_int!=NULL) free(pr_int);
    po_int=(int *)NULL; pr_int=(int *)NULL; 
    return 0;
  }
  nw = vector_instance_px(vv, &w);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  pr_int=(int *)ecalloc(nw, sizeof(int));
  po_int=(int *)ecalloc(nw, sizeof(int));
  for (j=0;j<nw;j++) {
    po_int[j]=(int)po[j];
    pr_int[j]=(int)pr[j];
  }
  cpfl=nw;
  return cpfl;
}
 
/*VERBATIM*/
static double chkspcp(void* vv) {
  int j, nw, npr, npo, flag;
  double *w, *pr, *po;
  nw = vector_instance_px(vv, &w);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  flag=1;
  if (po_int==NULL || pr_int==NULL) { cpfl=0; return 0; }
  if (cpfl!=nw) { flag=0;
  } else for (j=0;j<nw;j++) {
    if (po_int[j]!=(int)po[j] || pr_int[j]!=(int)pr[j]) {flag=0; continue;}
  }
  if (flag==0) {
    cpfl=0; free(po_int); free(pr_int); 
    po_int=(int *)NULL; pr_int=(int *)NULL; 
  }
  return flag;
}
 
/*VERBATIM*/
static double spmult(void* vv) {
  int i, j, nx, ny, nw, npr, npo, flag;
  double *x, *y, *w, *pr, *po, xx;
  ny = vector_instance_px(vv, &y);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  nw = vector_arg_px(3, &w);
  nx = vector_arg_px(4, &x);
  if (ifarg(5)) {flag=1;} else {flag=0;}
  if (nw!=npr || nw!=npo) {
    hoc_execerror("Sparse mat must have 3 identical size vecs for pre/post/wt", 0); 
  }
  if (flag==0) for (i=0;i<ny;i++) y[i] = 0.; // clear dest vec
  if (cpfl==0) {
    for (j=0;j<nw;j++) y[(int)po[j]] += (x[(int)pr[j]]*w[j]);
  } else if (cpfl!=nw) { hoc_execerror("cpfl!=nw in spmult", 0); } else {
    for (j=0;j<nw;j++) if (x[pr_int[j]]!=0) { y[po_int[j]] += ((x[pr_int[j]])*w[j]); }
  }
  return nx;
}
 
/*VERBATIM*/
static double spget(void* vv) {
  int j, nw, npr, npo;
  double *w, *pr, *po, row, col;
  nw = vector_instance_px(vv, &w);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  row = *getarg(3);
  col = *getarg(4);
  for (j=0;j<nw;j++) if (row==po[j]&&col==pr[j]) break;
  if (j==nw) return 0.; else return w[j];
}
 
/*VERBATIM*/
static double transpose(void* vv) {
  int i, j, nx, ny, rows, cols;
  double *x, *y;
  /* x will be the transpose of matrix y */
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  rows = (int)*getarg(2);
  cols = (int)*getarg(3);
  if (ny != nx) {
    hoc_execerror("Vector size mismatch", 0);
  }
  for (i=0;i<rows;i++) {
    for (j=0;j<cols;j++) {
      x[j*rows+i] = y[i*cols+j];
    }
  }
  return nx;
}
 
/*VERBATIM*/
static double mprintf(void* vv) {
  int i, j, nx, rows, cols;
  double *x;
  /* x will be printed out */
  nx = vector_instance_px(vv, &x);
  rows = (int)*getarg(1);
  cols = (int)*getarg(2);
  if (nx != rows*cols) {
    hoc_execerror("Vector size mismatch", 0);
  }
  for (i=0;i<rows;i++) {
    for (j=0;j<cols;j++) {
      printf("%g\t",x[i*cols+j]);
    }
    printf("\n");
  }
  return nx;
}
 
/*VERBATIM*/
static double mget(void* vv) {
  int i, j, nx, rows, cols;
  double *x;
  nx = vector_instance_px(vv, &x);
  i = (int)*getarg(1);
  j = (int)*getarg(2);
  cols = (int)*getarg(3);
  if (i*cols+j >= nx) {
    hoc_execerror("Indices out of bounds", 0);
  }
  return x[i*cols+j];
}
 
/*VERBATIM*/
static double mrow(void* vv) {
  int i, j, nx, ny, rows, cols;
  double *x, *y;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  i = (int)*getarg(2);
  cols = (int)*getarg(3);
  if (cols!=nx || i>=ny/cols) {
    hoc_execerror("Indices out of bounds", 0);
  }
  for (j=0;j<nx;j++) { x[j] = y[i*cols+j]; }
  return nx;
}
 
/*VERBATIM*/
static double mcol(void* vv) {
  int i, j, nx, ny, rows, cols;
  double *x, *y;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  j = (int)*getarg(2);
  cols = (int)*getarg(3);
  if (cols!=ny/nx || j>=cols) {
    hoc_execerror("Indices out of bounds", 0);
  }
  for (i=0;i<nx;i++) { x[i] = y[i*cols+j]; }
  return nx;
}
 
/*VERBATIM*/
static double mset(void* vv) {
  int i, j, nx, rows, cols;
  double *x, val;
  nx = vector_instance_px(vv, &x);
  i = (int)*getarg(1);
  j = (int)*getarg(2);
  cols = (int)*getarg(3);
  val = *getarg(4);
  if (i*cols+j >= nx) {
    hoc_execerror("Indices out of bounds", 0);
  }
  return (x[i*cols+j]=val);
}
 
static int  install_matrix (  ) {
   
/*VERBATIM*/
  /* the list of additional methods */
  install_vector_method("outprod", outprod);
  install_vector_method("mmult", mmult);
  install_vector_method("spmult", spmult);
  install_vector_method("spget", spget);
  install_vector_method("mkspcp", mkspcp);
  install_vector_method("chkspcp", chkspcp);
  install_vector_method("spltp", spltp);
  install_vector_method("transpose", transpose);
  install_vector_method("mprintf", mprintf);
  install_vector_method("mget", mget);
  install_vector_method("mset", mset);
  install_vector_method("mrow", mrow);
  install_vector_method("mcol", mcol);
  return 0; }
 
static void _hoc_install_matrix(void) {
  double _r;
   _r = 1.;
 install_matrix (  );
 hoc_retpushx(_r);
}

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
