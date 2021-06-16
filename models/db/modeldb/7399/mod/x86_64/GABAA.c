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
 
#define nrn_init _nrn_init__GABAA
#define _nrn_initial _nrn_initial__GABAA
#define nrn_cur _nrn_cur__GABAA
#define _nrn_current _nrn_current__GABAA
#define nrn_jacob _nrn_jacob__GABAA
#define nrn_state _nrn_state__GABAA
#define _net_receive _net_receive__GABAA 
#define initq initq__GABAA 
#define init_arrays init_arrays__GABAA 
#define prq prq__GABAA 
#define popqh2 popqh2__GABAA 
#define popqh1 popqh1__GABAA 
#define release release__GABAA 
#define setlink setlink__GABAA 
 
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
#define lpst _p[0]
#define maxsyn _p[1]
#define nsyn _p[2]
#define queu _p[3]
#define qlen _p[4]
#define begsyn _p[5]
#define endsyn _p[6]
#define newspk _p[7]
#define qptr _p[8]
#define i _p[9]
#define g _p[10]
#define R _p[11]
#define Ron _p[12]
#define Roff _p[13]
#define synon _p[14]
#define drive _p[15]
#define qq10 _p[16]
#define edt _p[17]
#define edb _p[18]
#define edc _p[19]
#define _g _p[20]
#define _nd_area  *_ppvar[0]._pval
 
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
 extern double celsius;
 /* declaration of user functions */
 static double _hoc_check();
 static double _hoc_code();
 static double _hoc_delay();
 static double _hoc_exptable();
 static double _hoc_getR();
 static double _hoc_getdbx();
 static double _hoc_gmax();
 static double _hoc_initq();
 static double _hoc_init_arrays();
 static double _hoc_link();
 static double _hoc_prq();
 static double _hoc_popqh2();
 static double _hoc_popqh1();
 static double _hoc_pre();
 static double _hoc_post();
 static double _hoc_release();
 static double _hoc_setlink();
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "check", _hoc_check,
 "code", _hoc_code,
 "delay", _hoc_delay,
 "exptable", _hoc_exptable,
 "getR", _hoc_getR,
 "getdbx", _hoc_getdbx,
 "gmax", _hoc_gmax,
 "initq", _hoc_initq,
 "init_arrays", _hoc_init_arrays,
 "link", _hoc_link,
 "prq", _hoc_prq,
 "popqh2", _hoc_popqh2,
 "popqh1", _hoc_popqh1,
 "pre", _hoc_pre,
 "post", _hoc_post,
 "release", _hoc_release,
 "setlink", _hoc_setlink,
 0, 0
};
#define _f_exptable _f_exptable_GABAA
#define check check_GABAA
#define code code_GABAA
#define delay delay_GABAA
#define exptable exptable_GABAA
#define getR getR_GABAA
#define getdbx getdbx_GABAA
#define gmax gmax_GABAA
#define link link_GABAA
#define pre pre_GABAA
#define post post_GABAA
 extern double _f_exptable( double );
 extern double check( );
 extern double code( );
 extern double delay( );
 extern double exptable( double );
 extern double getR( double );
 extern double getdbx( double , double );
 extern double gmax( );
 extern double link( double );
 extern double pre( );
 extern double post( );
 /* declare global and static user variables */
#define Alpha Alpha_GABAA
 double Alpha = 0.53;
#define Beta Beta_GABAA
 double Beta = 0.18;
#define CHAINLEN CHAINLEN_GABAA
 double CHAINLEN = 1;
#define Cdur Cdur_GABAA
 double Cdur = 1;
#define DELAY DELAY_GABAA
 double DELAY = 2;
#define Deadtime Deadtime_GABAA
 double Deadtime = 1;
#define Erev Erev_GABAA
 double Erev = -75;
#define GMAX GMAX_GABAA
 double GMAX = 1;
#define Rtau Rtau_GABAA
 double Rtau = 0;
#define Rinf Rinf_GABAA
 double Rinf = 0;
#define SPKSTORE SPKSTORE_GABAA
 double SPKSTORE = 50;
#define exptemp exptemp_GABAA
 double exptemp = 37;
#define q10 q10_GABAA
 double q10 = 1;
#define usetable usetable_GABAA
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_GABAA", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Cdur_GABAA", "ms",
 "Alpha_GABAA", "/ms",
 "Beta_GABAA", "/ms",
 "Erev_GABAA", "mV",
 "Deadtime_GABAA", "ms",
 "GMAX_GABAA", "uS",
 "Rtau_GABAA", "ms",
 "i", "nA",
 "g", "umho",
 0,0
};
 static double delta_t = 1;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "Cdur_GABAA", &Cdur_GABAA,
 "Alpha_GABAA", &Alpha_GABAA,
 "Beta_GABAA", &Beta_GABAA,
 "Erev_GABAA", &Erev_GABAA,
 "Deadtime_GABAA", &Deadtime_GABAA,
 "GMAX_GABAA", &GMAX_GABAA,
 "DELAY_GABAA", &DELAY_GABAA,
 "SPKSTORE_GABAA", &SPKSTORE_GABAA,
 "CHAINLEN_GABAA", &CHAINLEN_GABAA,
 "q10_GABAA", &q10_GABAA,
 "exptemp_GABAA", &exptemp_GABAA,
 "Rinf_GABAA", &Rinf_GABAA,
 "Rtau_GABAA", &Rtau_GABAA,
 "usetable_GABAA", &usetable_GABAA,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
 static void _destructor(Prop*);
 static void _constructor(Prop*);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"GABAA",
 "lpst",
 0,
 "maxsyn",
 "nsyn",
 "queu",
 "qlen",
 "begsyn",
 "endsyn",
 "newspk",
 "qptr",
 "i",
 "g",
 "R",
 "Ron",
 "Roff",
 "synon",
 0,
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 21, _prop);
 	/*initialize range parameters*/
 	lpst = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 21;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 2, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 if (!nrn_point_prop_) {_constructor(_prop);}
 
}
 static void _initlists();
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _GABAA_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 	register_destructor(_destructor);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 21, 2);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 GABAA /Users/awile/workbench/modldb/nmodldb/models/db/modeldb/7399/mod/GABAA.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_exptable;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int initq();
static int init_arrays(double);
static int prq();
static int popqh2();
static int popqh1(double);
static int release();
static int setlink();
 static double _n_exptable(double);
 
/*VERBATIM*/

#define PSTCAST ((SynS *)((unsigned long) lpst))
#define QUECAST ((QueU *)((unsigned long) queu))
#define PRECAST ((PreL *)((unsigned long) link))
#define QPRCAST ((QptR *)((unsigned long) qptr))

#define PREINC  50.  /* size increment for growing the presyn array */
#define POSTINC 5.   /* size increment for growing the postsyn array */

typedef struct QueU { /* postsynaptic queue */
  double time;  /* records spike time + delay */
  short index;  /* records index of associated SynS struct */
} QueU;

/* holds all the info for accessing the queu and its tail */
typedef struct QptR {
  struct QueU *qq;     /* will point to the queu */
  double *nspk; /* will point to newspk */
  double *head; /* will point to begsyn */
  double *qln;  /* points to qlen */
  double dead;  /* Cdur + Deadtime */
  double qterm;  /* the last time put on this queue */
  int cpost;  /* user code for info about postsyn cell */
} QptR;

typedef struct SynS { /* postsynaptic structure */
  /* this stuff used presynaptically */
  short index;  /* index of this item */
  double del;   /* param: individualized delay set by user */
  double spkt;  /* state: time when next spike can begin */
  struct QptR *qpt;   /* indirection to get at the queu */

  /* for allowing mult spikes to add up */
  short chainlen;
  short chainptr;

  /* user codes */
  int ucode;  /* user code for user convenience, an arbitrary index */
  /* this stuff used postsynaptically */
  double pgm;   /* percent of gmax - set by user */
  double Rcurr; /* holds current value of R (release status), internal */
  double last;  /* time last syn will end, internal */
  int *pcpre;   /* pointer to access the presyn code */
  struct SynS **back;  /* pointer to presyn pointer to this entry */
} SynS;

/* holds all the info for accessing presyn list */
typedef struct PreL {
  int cpre;   /* user code for info about presyn cell */
  struct SynS **plst;        /* will point to post for postsyn nrn*/
  double *nnpre;     /* number of pointers currently there */
  double *mxpre;     /* number permitted before realloc required */
  double link2;       /* a copy of the link information to confirm correct hook up */
} PreL;

 
static int  init_arrays (  double _lnum ) {
   
/*VERBATIM*/
{
    int ii, msn;
    if (_lnum == maxsyn) {
      printf("Clearing array\n");
      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) { /* nullize pointers */
	*(PSTCAST[ii].back) = (SynS *)NULL; /* remove old pre-post pointer */
      }
      nsyn = 0;
    } else {
      maxsyn = (double)((int)_lnum);
      msn = maxsyn*CHAINLEN;
      qlen = SPKSTORE*msn;

      if (lpst == 0) {
        nsyn = 0.;
        lpst = (double)((unsigned long)hoc_Ecalloc(msn,(sizeof(SynS)))); 
        queu = (double)((unsigned long)hoc_Ecalloc(((int)qlen), (sizeof(QueU)))); 
      } else {
        if (nsyn > maxsyn) {
          printf("Shrinking array\n");
          for (ii=maxsyn;ii<nsyn*CHAINLEN;ii+=CHAINLEN) { /* nullize pointers */
	    /* remove old pre-post pointer */
	    *(PSTCAST[ii].back) = (SynS *)NULL;
          }
          nsyn = maxsyn;
        } else {
          printf("Expanding array %09d\n",QPRCAST->cpost);
        }
        lpst = (double)((unsigned long)hoc_Erealloc((PSTCAST),((int)(msn*sizeof(SynS))))); 
        queu = (double)((unsigned long)hoc_Erealloc((QUECAST),((int)(qlen*sizeof(QueU))))); 
        for (ii=0;ii<(int)nsyn*CHAINLEN;ii+=CHAINLEN) {
          *(PSTCAST[ii].back) = &(PSTCAST[ii]); /* restore old pre pointers */
        }
      }

      /* printf("Initializing arrays: Syn= %d bytes,Queue= %d bytes.\n",
	 (int)msn*sizeof(SynS),(int)qlen, sizeof(QueU)); */
      QPRCAST->qq = QUECAST;
    }
  }
  return 0; }
 
static double _hoc_init_arrays(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 init_arrays (  *getarg(1) );
 return(_r);
}
 
/*VERBATIM*/
static void hshake();
 
static int  setlink (  ) {
   
/*VERBATIM*/
{ 
  int ii, x, new;
  SynS *sns;
  double ptemp;
  PreL *ppsyn;
    
  if (ifarg(2)) { 
    x = (int)*getarg(1);
    ptemp = *getarg(2);		/* pick up pointer */
  } else {
    x = nsyn;
    ptemp = *getarg(1);
  }

  if (x >= maxsyn) {
    init_arrays(maxsyn+POSTINC); /* #DEFINE POSTINC 5 */
  } 
  if (x > nsyn) {
    hoc_execerror("Can't leave empty pointers: see nsyn for current array index.", 0);
  }

  sns = &(PSTCAST[x*(int)CHAINLEN]);	/* postsynaptic entry */
  if ((ppsyn = (PreL *)((unsigned long)ptemp)) == (PreL *)NULL) { /* presyn loc */
    hoc_execerror("Presyn not initialized.", 0); }
  if (ppsyn->link2 != ptemp) {
    hoc_execerror("Invalid (changed) link.", 0); }

  if (x == nsyn) {
    nsyn++;		/* a new entry */
    new = 1;
  } else { 
    /* should generate error if try to change something in middle of chain */
    if (sns->chainlen == -2) {
      hoc_execerror("Internal error: Index used must be multiple of CHAINLEN.", 0);
    }
    new = 0; 
  }

  hshake(sns,ppsyn,new);

  x *= (int)CHAINLEN;
  for (ii=x;ii < x + CHAINLEN;ii++) {
    sns = &(PSTCAST[ii]);	/* access the entry */
    /* initialize postsyn stuff to reasonable values */    
    /* delay, gmax and all codes should be initialized separately */
    sns->del = DELAY;
    sns->ucode = -1;
    sns->chainlen = (ii==x)?CHAINLEN:-2;
    sns->pgm = 1.0;

    /* initialize the presynaptic stuff */
    /* Rcurr and last will be initialized with init */
    sns->index = x;
    sns->qpt = QPRCAST; 
  }
}
  return 0; }
 
static double _hoc_setlink(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 setlink (  );
 return(_r);
}
 
/*VERBATIM*/
static void hshake(ss,pl,flag) 
/* ls will be a pointer to presyn cell's array of pointers */
SynS *ss;
PreL *pl;
int flag; /* flag == 1 if this is a brand new entry */
{
  int ii;
  double *nn;
  double *mx;

  nn = pl->nnpre; /* pointer to presyn npre */
  mx = pl->mxpre; /* pointer to presyn maxpre */

  /* erase presyn pointer if this has been set before */
  if (flag == 0) {		/* an old entry */
    /* fall out of loop if a pointer exists already */
    for (ii=0;ii<(*nn) && (pl->plst[ii])!=ss;ii++); /* no body */
    if (ii < (*nn)) {	/* fell out of loop */
      printf("Maintaining pointers.\n");
    } else {
      flag = 1;			/* now we do have to create a new pointer */
      printf("Erasing pointer (C%09d,Pr%09d,Po%09d).\n",
	     ss->ucode,*(ss->pcpre),QPRCAST->cpost);
      *(ss->back) = (SynS *)NULL; /* remove old pre-post pointer */
    }
  }

  if (flag == 1) {		/* create a new pointer */
    /* if necessary, allocate space for the pointer system */
    if ((*nn) == 0.) { 
      (*mx) = PREINC;  /* #DEFINE PREINC 50 */
      pl->plst = (SynS **)ecalloc(((int)(*mx)), sizeof(SynS *)); 
    } else if ((*nn) == (*mx)) { /* need to create more space */
      (*mx) += PREINC;
      pl->plst = (SynS **)realloc(pl->plst,((int)(*mx)) * sizeof(SynS *)); 
      /* reassign all back pointers */
      for (ii=0;ii<(*nn);ii++) {
	(pl->plst[ii])->back = &(pl->plst[ii]);
      }
    }

    /* the handshake */
    pl->plst[(int)(*nn)] = ss;	/* pre  -> post */
    ss->back = &(pl->plst[(int)(*nn)]); /* post -> pre */
    ss->pcpre = &(pl->cpre); /* post -> precode */
    (*nn)++;
  }
}
 
double code (  ) {
   double _lcode;
 
/*VERBATIM*/
{ int ii;
    if (ifarg(1)) {
      ii = (int)*getarg(1);
      if (ii < 0) { ii = nsyn+ii; }
      if (ii >= nsyn || ii < 0) {hoc_execerror("array index out of bounds", 0);}
      if (ifarg(2)) {  /* look for a second arg to do a set */
	(PSTCAST[ii*(int)CHAINLEN]).ucode = (int)*getarg(2); }
      _lcode = (PSTCAST[ii*(int)CHAINLEN]).ucode;
    } else {
      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) {
	  printf("%d,%09d  ",ii,(PSTCAST[ii]).ucode);
	}
      _lcode = 1.0;
    }}      
 
return _lcode;
 }
 
static double _hoc_code(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r =  code (  );
 return(_r);
}
 
double post (  ) {
   double _lpost;
 
/*VERBATIM*/
{ 
    if (ifarg(2)) { QPRCAST->cpost = (int)*getarg(2); }
    _lpost = QPRCAST->cpost;
  }
 
return _lpost;
 }
 
static double _hoc_post(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r =  post (  );
 return(_r);
}
 
double delay (  ) {
   double _ldelay;
 
/*VERBATIM*/
{ int ii,jj;
    if (ifarg(1)) {
      ii = (int)*getarg(1);
      if (ii < 0) { ii = nsyn+ii; }
      if (ii >= nsyn || ii < 0) {hoc_execerror("array index out of bounds", 0);}
      if (ifarg(2)) {  /* look for a second arg to do a set */
	for (jj=ii*CHAINLEN;jj<(ii+1)*CHAINLEN;jj++) {
 	  (PSTCAST[jj]).del = *getarg(2) + DELAY; }}
      _ldelay = (PSTCAST[ii*(int)CHAINLEN]).del;
    } else {
      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) {
	  printf("%d,%g  ",ii,(PSTCAST[ii]).del);
	}
      _ldelay = 1.0;
    }}
 
return _ldelay;
 }
 
static double _hoc_delay(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r =  delay (  );
 return(_r);
}
 
double gmax (  ) {
   double _lgmax;
 
/*VERBATIM*/
{ int ii,jj;
    if (ifarg(1)) {
      ii = (int)*getarg(1);
      if (ii < 0) { ii = nsyn+ii; }
      if (ii >= nsyn || ii < 0) {hoc_execerror("array index out of bounds", 0);}
      if (ifarg(2)) {  /* look for a second arg to do a set */
 	for (jj=ii*CHAINLEN;jj<(ii+1)*CHAINLEN;jj++) {
 	  (PSTCAST[jj]).pgm = *getarg(2); }}
      _lgmax = (PSTCAST[ii*(int)CHAINLEN]).pgm;
    } else {
      printf("Multiply by %g to get effective gmax.\n",GMAX);
      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) {
	  printf("%d,%g  ",ii,(PSTCAST[ii]).pgm);
	}
      _lgmax = 1.0;
    }}
 
return _lgmax;
 }
 
static double _hoc_gmax(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r =  gmax (  );
 return(_r);
}
 
double pre (  ) {
   double _lpre;
 
/*VERBATIM*/
{ int x,ii;
    x = -2;  /* flag -> -1 then print out pointers */
    if (ifarg(1)) { x = (int)*getarg(1); }
    if (x >= 0) {
      if (x >= nsyn) { hoc_execerror("array index out of bounds", 0);}
      x *= CHAINLEN;
      _lpre = *((PSTCAST[x]).pcpre);
    } else {	  
      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) {
	printf("%1d Del:%-2g Gmax: %-4g Cd:%09d Pre:%09d Post:%09d\n",
	       (PSTCAST[ii]).index,
	       (PSTCAST[ii]).del,
	       (PSTCAST[ii]).pgm,
	       (PSTCAST[ii]).ucode,
	       *((PSTCAST[ii]).pcpre),
	       QPRCAST->cpost);
	if (x == -1 || &(PSTCAST[ii]) != *((PSTCAST[ii]).back)) {
	  printf("\t\t%s: %x -> %x\n",
		 ((&(PSTCAST[ii])==*((PSTCAST[ii]).back))?"OK":"POINTER MISMATCH"),
		 &(PSTCAST[ii]),
		 (PSTCAST[ii]).back);
	}
      }
      _lpre = 1.0;
    }
  } 
 
return _lpre;
 }
 
static double _hoc_pre(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r =  pre (  );
 return(_r);
}
 
double link (  double _lii ) {
   double _llink;
 
/*VERBATIM*/
{ 
    if (_lii >= nsyn) { hoc_execerror("array index out of bounds", 0);}
    _llink = (double)(unsigned long)(PSTCAST[(int)_lii]).pcpre;
  }
 
return _llink;
 }
 
static double _hoc_link(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r =  link (  *getarg(1) );
 return(_r);
}
 
double check (  ) {
   double _lcheck;
 
/*VERBATIM*/
{ 
    int ii;
    for (ii=0; ii<nsyn*CHAINLEN && _lcheck==1.; ii+=CHAINLEN) {
      if (&(PSTCAST[ii]) != *((PSTCAST[ii]).back)) {
	printf("****************************************************************\n");
	printf("ERROR:: Index:%3d,C%09d,Pr%09d,Po%09d, Delay:%6g (%x->%x->%x)\n",
	       (PSTCAST[ii]).index,
	       (PSTCAST[ii]).ucode,
	       *((PSTCAST[ii]).pcpre),
	       QPRCAST->cpost,
	       (PSTCAST[ii]).del,
	       &(PSTCAST[ii]),
	       (PSTCAST[ii]).back,
	       *((PSTCAST[ii]).back));
	printf("****************************************************************\n");
	_lcheck = -1.;
      }
    }
    _lcheck = nsyn;  /* all pointers have to be active */
  }
 
return _lcheck;
 }
 
static double _hoc_check(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r =  check (  );
 return(_r);
}
 
static int  initq (  ) {
   
/*VERBATIM*/
{
  int ii;
  QPRCAST->dead = Cdur + Deadtime;
  begsyn = endsyn = newspk = 0.;
  for (ii=0;ii<qlen;ii++) {
    QUECAST[ii].time = 1.e20;
    QUECAST[ii].index = -1;
  }
  QPRCAST->qterm = -1.e2;
  for (ii=0;ii<nsyn*CHAINLEN;ii++) {
    PSTCAST[ii].Rcurr = 0.;
    PSTCAST[ii].last = -1.e2;
    PSTCAST[ii].spkt = -1.e2;
    if (PSTCAST[ii].chainlen > 1) { PSTCAST[ii].chainptr = -1; }
  }
  }
  return 0; }
 
static double _hoc_initq(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 initq (  );
 return(_r);
}
 
static int  popqh1 (  double _laug ) {
   
/*VERBATIM*/
{  
    if (QUECAST[(int)begsyn].time == 1e20) {
      printf("%x %g %g ",QUECAST,newspk,qlen);
      hoc_execerror("Error: queue exhausted.\n",0);
    } else {  /* augment the time by Cdur */
      QUECAST[(int)begsyn].time += _laug;
      begsyn++;
      if (begsyn == qlen) { begsyn = 0.; }
    }
  }
  return 0; }
 
static double _hoc_popqh1(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 popqh1 (  *getarg(1) );
 return(_r);
}
 
static int  popqh2 (  ) {
   
/*VERBATIM*/
{  
    QUECAST[(int)endsyn].time = 1.e20;  /* clear the entry */
    endsyn++;
    if (endsyn == qlen) { endsyn = 0.; }
  }
  return 0; }
 
static double _hoc_popqh2(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 popqh2 (  );
 return(_r);
}
 
double getdbx (  double _lc , double _lx ) {
   double _lgetdbx;
 
/*VERBATIM*/
{ 
  switch ((int)_lc) {
  case 1:
    _lgetdbx = (PSTCAST[(int)_lx]).last;
    break;
  case 2:
    _lgetdbx = (PSTCAST[(int)_lx]).Rcurr;
    break;
  case 3:
    _lgetdbx = (PSTCAST[(int)_lx]).spkt;
    break;
  case 4:
    _lgetdbx = (QUECAST[(int)_lx]).time;
    break;
  case 5:
    _lgetdbx = (QUECAST[(int)_lx]).index;
    break;
  case 6:
    _lgetdbx = QPRCAST->qterm;
    break;
  case 7:
    _lgetdbx = QPRCAST->dead;
    break;
  case 8:
    _lgetdbx = (PSTCAST[(int)_lx]).chainlen;
    break;
  case 9:
    _lgetdbx = (PSTCAST[(int)_lx]).chainptr;
    break;
  default:
    hoc_execerror("UNAVAILABLE IN DBX",0);
    _lgetdbx = -1;
    break;
  }}
 
return _lgetdbx;
 }
 
static double _hoc_getdbx(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r =  getdbx (  *getarg(1) , *getarg(2) );
 return(_r);
}
 
static int  prq (  ) {
   
/*VERBATIM*/
{ int ii;
    printf("new:%g beg:%g end:%g /%g\n",newspk,begsyn,endsyn,qlen);
    for (ii=endsyn;ii!=newspk;ii=((ii==qlen-1)?0:ii+1)) {
        printf("%d %8g%8d\n",ii,(QUECAST[ii]).time,(QUECAST[ii]).index);
    }
  }
  return 0; }
 
static double _hoc_prq(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 prq (  );
 return(_r);
}
 
static int  release (  ) {
   if ( nsyn > 0.0 ) {
     
/*VERBATIM*/

  int who;
  QueU *pqueu;
  SynS *ppst;

  pqueu = (QueU *)((unsigned long) queu);

  while (t >= pqueu[(int)begsyn].time) { /*  somebody spiked delay time ago */
    ppst = (SynS *)((unsigned long) lpst);
    who = pqueu[(int)begsyn].index; /* who spiked? */
    /* calculate the decay that occurred since last activity ended */
    ppst[who].Rcurr *= exptable(-Beta*(t-ppst[who].last));
    /* transfer the value from Roff to Ron */
    Ron += ppst[who].Rcurr;
    Roff -= ppst[who].Rcurr;
    synon += ppst[who].pgm;	/*  weight syn by percent gmax */
    ppst[who].last = t + Cdur;   /* time when syn will turn off */
    popqh1(Cdur);		/* next (also add Cdur to value on the queu) */
  }

  while (t >= pqueu[(int)endsyn].time) { /*  somebody needs to be turned off */
    ppst = (SynS *)((unsigned long) lpst); /* will do this assign twice in rare cases */
    who = pqueu[(int)endsyn].index;   /* who spiked? */
    /* solve Rcurr forward in time till end of syn activity */
    ppst[who].Rcurr = ppst[who].pgm*Rinf*(1.-edc) + ppst[who].Rcurr*edc;
    Ron -= ppst[who].Rcurr;
    Roff += ppst[who].Rcurr;  /* transfer from on to off */
    synon -= ppst[who].pgm;   /* remove this percent of gmax */
    if (synon<1e-11 && synon>-1e-11) { synon=0.; }
    if (synon==0. || Ron < 0.) { Ron = 0.; }
    popqh2();  /* next */
  }

  /*  update R */
  if (synon > 0) {		/*  one or more synapses turned on? */
    Ron = synon*drive + Ron * edt; /*  drive R toward Rinf */
  } 
  Roff *= edb;			/*  Roff always decays toward 0 */
  return 0;
 }
    return 0; }
 
static double _hoc_release(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 release (  );
 return(_r);
}
 
double getR (  double _lx ) {
   double _lgetR;
 
/*VERBATIM*/
{
    SynS pst;
    double rnow;
    double end, rinf;

    pst = (PSTCAST[(int)_lx]);
    end = pst.last;
    rinf = pst.pgm * Rinf;

    if (end < 0.) {		/* not yet turned on */
      rnow = 0.;
    } else if (t >= end - dt/2) {	/* decay */
      rnow = pst.Rcurr * exptable(-Beta*(t-end));
    } else {			/* turning on */
      rnow = rinf + (pst.Rcurr - rinf)*exptable((t-(end-Cdur))/(-Rtau));
    }
    if (pst.pgm != 0.) {
      _lgetR = rnow/pst.pgm;  /* scale it to 1 so it looks like a state variable */
    } else {
      _lgetR = 0.;
    }
  }
 
return _lgetR;
 }
 
static double _hoc_getR(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r =  getR (  *getarg(1) );
 return(_r);
}
 static double _mfac_exptable, _tmin_exptable;
 static void _check_exptable();
 static void _check_exptable() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_exptable =  - 10.0 ;
   _tmax =  0.0 ;
   _dx = (_tmax - _tmin_exptable)/1000.; _mfac_exptable = 1./_dx;
   for (_i=0, _x=_tmin_exptable; _i < 1001; _x += _dx, _i++) {
    _t_exptable[_i] = _f_exptable(_x);
   }
  }
 }

 double exptable(double _lx){ _check_exptable();
 return _n_exptable(_lx);
 }

 static double _n_exptable(double _lx){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 return _f_exptable(_lx); 
}
 _xi = _mfac_exptable * (_lx - _tmin_exptable);
 if (isnan(_xi)) {
  return _xi; }
 if (_xi <= 0.) {
 return _t_exptable[0];
 }
 if (_xi >= 1000.) {
 return _t_exptable[1000];
 }
 _i = (int) _xi;
 return _t_exptable[_i] + (_xi - (double)_i)*(_t_exptable[_i+1] - _t_exptable[_i]);
 }

 
double _f_exptable (  double _lx ) {
   double _lexptable;
 if ( ( _lx > - 10.0 )  && ( _lx < 0.0 ) ) {
     _lexptable = exp ( _lx ) ;
     }
   else {
     _lexptable = 0. ;
     }
   
return _lexptable;
 }
 
static double _hoc_exptable(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
  _r =  exptable (  *getarg(1) );
 return(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("GABAA", "cannot be used with CVODE"); return 0;}
 
static void _constructor(Prop* _prop) {
	_p = _prop->param; _ppvar = _prop->dparam;
{
 {
   
/*VERBATIM*/
{
    qptr = (double)((unsigned long)hoc_Ecalloc(1, (sizeof(QptR)))); 
    /* qptr allows presynaptic cell to easily manipulate the post queu */
    QPRCAST->qln = &qlen;
    QPRCAST->nspk = &newspk;
    QPRCAST->head = &begsyn;
    nsyn = maxsyn = qlen = lpst = queu = 0.;
    if (ifarg(2)) { 
      QPRCAST->cpost = (int)*getarg(2);
    } else { QPRCAST->cpost = -1; }
    if (ifarg(3)) { init_arrays(*getarg(3)); }
  }
 }
 
}
}
 
static void _destructor(Prop* _prop) {
	_p = _prop->param; _ppvar = _prop->dparam;
{
 {
   
/*VERBATIM*/
{ int ii;
      free(QUECAST);  
      if (lpst != 0) {
        for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) { /* nullize pointers */
	  *(PSTCAST[ii].back) = (SynS *)NULL; /* remove old pre-post pointer */
        }      
        free(PSTCAST);
        free(QPRCAST);  
      }
      nsyn = maxsyn = qlen = 0.;
      lpst = queu = qptr = 0.; 
  }
 }
 
}
}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
 {
   Rinf = Alpha / ( Alpha + Beta ) ;
   qq10 = pow( q10 , ( ( celsius - exptemp ) / 10. ) ) ;
   Rtau = 1.0 / ( Alpha + Beta ) / qq10 ;
   edt = exp ( - dt / Rtau ) ;
   edb = exp ( - Beta * dt ) ;
   edc = exp ( - Cdur / Rtau ) ;
   drive = Rinf * ( 1. - edt ) ;
   synon = 0.0 ;
   R = 0.0 ;
   Ron = 0.0 ;
   Roff = 0.0 ;
   if ( nsyn > 0.0 ) {
     initq ( _threadargs_ ) ;
     }
   else {
     
/*VERBATIM*/
    printf("WARNING nsyn==0 ");
 }
   }
  _sav_indep = t; t = _save;

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

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   R = Ron + Roff ;
   g = GMAX * R ;
   i = g * ( v - Erev ) ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
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
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

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
 { error =  release();
 if(error){fprintf(stderr,"at line 131 in file sns.inc:\n    SOLVE release\n"); nrn_complain(_p); abort_run(error);}
 }}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_exptable = makevector(1001*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/awile/workbench/modldb/nmodldb/models/db/modeldb/7399/mod/GABAA.mod";
static const char* nmodl_file_text = 
  "NEURON {  POINT_PROCESS GABAA }\n"
  ":   GABA SYNAPSE (GABA-A receptors)\n"
  ":\n"
  ":   Parameters estimated from whole cell recordings of GABA-A synaptic currents \n"
  ":   from dentate granule cells: Otis TS and Mody I (1992) Modulation of decay \n"
  ":   kinetics and frequency of GABA_A receptor-mediated spontaneous inhibitory\n"
  ":   postsynaptic currents in hippocampal neurons.  Neurosci. 49: 13-32.\n"
  ":\n"
  ":   GABALOW was created to allow the use two types of GABA-A currents in the\n"
  ":   same simulation.\n"
  "PARAMETER {\n"
  "  Cdur	= 1.0	(ms)		: transmitter duration (rising phase)\n"
  "  Alpha	= 0.53	(/ms mM)	: forward (binding) rate\n"
  "  Beta	= 0.18	(/ms)		: backward (unbinding) rate\n"
  "  Erev	= -75	(mV)		: reversal potential\n"
  "  Deadtime = 1.	(ms)		: mimimum time between release events\n"
  "  GMAX = 1	(uS)		: maximum conductance\n"
  "  DELAY = 2\n"
  "}\n"
  ":::INCLUDE \"sns.inc\"\n"
  ":::realpath /Users/awile/workbench/modldb/nmodldb/models/db/modeldb/7399/mod/sns.inc\n"
  ": $Id: sns.inc,v 1.29 1997/03/25 00:05:44 billl Exp $\n"
  "COMMENT\n"
  "USAGE: for most receptors\n"
  "*****************************************************************************\n"
  "    NEURON {\n"
  "      POINT_PROCESS NAME\n"
  "    }\n"
  "\n"
  "    PARAMETER {\n"
  "      Cdur	= 1.08	(ms)		: transmitter duration (rising phase)\n"
  "      Alpha	= 1	(/ms mM)	: forward (binding) rate\n"
  "      Beta	= 0.02	(/ms)		: backward (unbinding) rate\n"
  "      Erev	= -80	(mV)		: reversal potential\n"
  "      Deadtime = 1	(ms)		: mimimum time between release events\n"
  "      GMAX = 1		(umho)		: reference conductance\n"
  "      DELAY = 0         (ms)\n"
  "    }\n"
  "\n"
  "    \n"
  "    INCLUDE \"sns.inc\"\n"
  "*****************************************************************************\n"
  "\n"
  "USAGE: for NMDA receptor\n"
  "*****************************************************************************\n"
  "    NEURON{ POINT_PROCESS NMDA\n"
  "      RANGE B \n"
  "    }\n"
  "\n"
  "    PARAMETER {\n"
  "      mg        = 1.    (mM)     : external magnesium concentration\n"
  "      Cdur	= 1.	(ms)	 : transmitter duration (rising phase)\n"
  "      Alpha	= 4.	(/ms mM) : forward (binding) rate\n"
  "      Beta	= 0.0067 (/ms)	 : backward (unbinding) rate 1/150\n"
  "      Erev	= 0.	(mV)	 : reversal potential\n"
  "      Deadtime = 1	(ms)	 : mimimum time between release events\n"
  "      GMAX     = 1      (umho)   : reference conductance\n"
  "      DELAY    = 0               : axonal delay\n"
  "    }\n"
  "\n"
  "    ASSIGNED { B }\n"
  "\n"
  "    INCLUDE \"sns.inc\"\n"
  "    : EXTRA BREAKPOINT MUST BE BELOW THE INCLUDE\n"
  "    BREAKPOINT {\n"
  "      rates(v)\n"
  "      g = g * B : but don't really need to readjust conductance\n"
  "      i = i * B : i = g*(v - Erev)\n"
  "    }\n"
  "\n"
  "    PROCEDURE rates(v(mV)) {\n"
  "      TABLE B\n"
  "      DEPEND mg\n"
  "      FROM -100 TO 80 WITH 180\n"
  "      B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))\n"
  "    }\n"
  "*****************************************************************************\n"
  "\n"
  "see end for implementation comments\n"
  "ENDCOMMENT\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "  NONSPECIFIC_CURRENT i\n"
  "  RANGE R, g\n"
  "  RANGE Ron, Roff, synon  : accessible for debugging\n"
  "  GLOBAL GMAX, DELAY, Cdur, Alpha, Beta, Erev, Deadtime, Rinf, Rtau, q10, exptemp\n"
  "}\n"
  " \n"
  ":::INCLUDE \"snsarr.inc\"  : array management routines\n"
  ":::realpath /Users/awile/workbench/modldb/nmodldb/models/db/modeldb/7399/mod/snsarr.inc\n"
  ": $Id: snsarr.inc,v 1.63 1998/10/07 15:06:12 billl Exp $\n"
  "\n"
  "COMMENT\n"
  "manage arrays used in sns.inc\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "  RANGE maxsyn, nsyn, qlen                   : scalars\n"
  "  RANGE lpst, queu                            : SynS bzw. QueU structures\n"
  "  RANGE begsyn, endsyn, newspk\n"
  "  GLOBAL SPKSTORE, CHAINLEN\n"
  "  RANGE qptr                                 : QptR structure\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "  SPKSTORE = 50                  : max number of spikes expected during max delay\n"
  "  CHAINLEN = 1                  : how many extra entries in a chain\n"
  "  lpst = 0                      : pointer to postsynaptic array\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "  maxsyn        : max and counter for this array\n"
  "  nsyn		\n"
  "		\n"
  "  queu          : the postsyn queu\n"
  "  qlen          : its length (maxsyn*SPKSTORE)\n"
  "  begsyn        : index into queu that tells time to start\n"
  "  endsyn        : index into queu that tells time to end\n"
  "  newspk        : index into queu to place spk time + delay (tail)\n"
  "\n"
  "  qptr          : will point to the queu, its tail and its length\n"
  "}\n"
  "\n"
  ":::INCLUDE \"snshead.inc\"\n"
  ":::realpath /Users/awile/workbench/modldb/nmodldb/models/db/modeldb/7399/mod/snshead.inc\n"
  ": $Id: snshead.inc,v 1.22 1998/01/07 20:22:30 billl Exp $\n"
  "\n"
  "VERBATIM\n"
  "\n"
  "#define PSTCAST ((SynS *)((unsigned long) lpst))\n"
  "#define QUECAST ((QueU *)((unsigned long) queu))\n"
  "#define PRECAST ((PreL *)((unsigned long) link))\n"
  "#define QPRCAST ((QptR *)((unsigned long) qptr))\n"
  "\n"
  "#define PREINC  50.  /* size increment for growing the presyn array */\n"
  "#define POSTINC 5.   /* size increment for growing the postsyn array */\n"
  "\n"
  "typedef struct QueU { /* postsynaptic queue */\n"
  "  double time;  /* records spike time + delay */\n"
  "  short index;  /* records index of associated SynS struct */\n"
  "} QueU;\n"
  "\n"
  "/* holds all the info for accessing the queu and its tail */\n"
  "typedef struct QptR {\n"
  "  struct QueU *qq;     /* will point to the queu */\n"
  "  double *nspk; /* will point to newspk */\n"
  "  double *head; /* will point to begsyn */\n"
  "  double *qln;  /* points to qlen */\n"
  "  double dead;  /* Cdur + Deadtime */\n"
  "  double qterm;  /* the last time put on this queue */\n"
  "  int cpost;  /* user code for info about postsyn cell */\n"
  "} QptR;\n"
  "\n"
  "typedef struct SynS { /* postsynaptic structure */\n"
  "  /* this stuff used presynaptically */\n"
  "  short index;  /* index of this item */\n"
  "  double del;   /* param: individualized delay set by user */\n"
  "  double spkt;  /* state: time when next spike can begin */\n"
  "  struct QptR *qpt;   /* indirection to get at the queu */\n"
  "\n"
  "  /* for allowing mult spikes to add up */\n"
  "  short chainlen;\n"
  "  short chainptr;\n"
  "\n"
  "  /* user codes */\n"
  "  int ucode;  /* user code for user convenience, an arbitrary index */\n"
  "  /* this stuff used postsynaptically */\n"
  "  double pgm;   /* percent of gmax - set by user */\n"
  "  double Rcurr; /* holds current value of R (release status), internal */\n"
  "  double last;  /* time last syn will end, internal */\n"
  "  int *pcpre;   /* pointer to access the presyn code */\n"
  "  struct SynS **back;  /* pointer to presyn pointer to this entry */\n"
  "} SynS;\n"
  "\n"
  "/* holds all the info for accessing presyn list */\n"
  "typedef struct PreL {\n"
  "  int cpre;   /* user code for info about presyn cell */\n"
  "  struct SynS **plst;        /* will point to post for postsyn nrn*/\n"
  "  double *nnpre;     /* number of pointers currently there */\n"
  "  double *mxpre;     /* number permitted before realloc required */\n"
  "  double link2;       /* a copy of the link information to confirm correct hook up */\n"
  "} PreL;\n"
  "\n"
  "ENDVERBATIM\n"
  ":::end INCLUDE snshead.inc\n"
  ":  QptR qptr;  /* governs access to the postsyn queue */\n"
  ":  PreL prel;  /* governs access to the presyn list */\n"
  "\n"
  ": takes 1-3 arguments: location, post code, maxsyn\n"
  "CONSTRUCTOR {\n"
  "  VERBATIM {\n"
  "    qptr = (double)((unsigned long)hoc_Ecalloc(1, (sizeof(QptR)))); \n"
  "    /* qptr allows presynaptic cell to easily manipulate the post queu */\n"
  "    QPRCAST->qln = &qlen;\n"
  "    QPRCAST->nspk = &newspk;\n"
  "    QPRCAST->head = &begsyn;\n"
  "    nsyn = maxsyn = qlen = lpst = queu = 0.;\n"
  "    if (ifarg(2)) { \n"
  "      QPRCAST->cpost = (int)*getarg(2);\n"
  "    } else { QPRCAST->cpost = -1; }\n"
  "    if (ifarg(3)) { init_arrays(*getarg(3)); }\n"
  "  }\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  "DESTRUCTOR {\n"
  "  VERBATIM { int ii;\n"
  "      free(QUECAST);  \n"
  "      if (lpst != 0) {\n"
  "        for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) { /* nullize pointers */\n"
  "	  *(PSTCAST[ii].back) = (SynS *)NULL; /* remove old pre-post pointer */\n"
  "        }      \n"
  "        free(PSTCAST);\n"
  "        free(QPRCAST);  \n"
  "      }\n"
  "      nsyn = maxsyn = qlen = 0.;\n"
  "      lpst = queu = qptr = 0.; \n"
  "  }\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  "PROCEDURE init_arrays(num) {\n"
  "  VERBATIM {\n"
  "    int ii, msn;\n"
  "    if (_lnum == maxsyn) {\n"
  "      printf(\"Clearing array\\n\");\n"
  "      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) { /* nullize pointers */\n"
  "	*(PSTCAST[ii].back) = (SynS *)NULL; /* remove old pre-post pointer */\n"
  "      }\n"
  "      nsyn = 0;\n"
  "    } else {\n"
  "      maxsyn = (double)((int)_lnum);\n"
  "      msn = maxsyn*CHAINLEN;\n"
  "      qlen = SPKSTORE*msn;\n"
  "\n"
  "      if (lpst == 0) {\n"
  "        nsyn = 0.;\n"
  "        lpst = (double)((unsigned long)hoc_Ecalloc(msn,(sizeof(SynS)))); \n"
  "        queu = (double)((unsigned long)hoc_Ecalloc(((int)qlen), (sizeof(QueU)))); \n"
  "      } else {\n"
  "        if (nsyn > maxsyn) {\n"
  "          printf(\"Shrinking array\\n\");\n"
  "          for (ii=maxsyn;ii<nsyn*CHAINLEN;ii+=CHAINLEN) { /* nullize pointers */\n"
  "	    /* remove old pre-post pointer */\n"
  "	    *(PSTCAST[ii].back) = (SynS *)NULL;\n"
  "          }\n"
  "          nsyn = maxsyn;\n"
  "        } else {\n"
  "          printf(\"Expanding array %09d\\n\",QPRCAST->cpost);\n"
  "        }\n"
  "        lpst = (double)((unsigned long)hoc_Erealloc((PSTCAST),((int)(msn*sizeof(SynS))))); \n"
  "        queu = (double)((unsigned long)hoc_Erealloc((QUECAST),((int)(qlen*sizeof(QueU))))); \n"
  "        for (ii=0;ii<(int)nsyn*CHAINLEN;ii+=CHAINLEN) {\n"
  "          *(PSTCAST[ii].back) = &(PSTCAST[ii]); /* restore old pre pointers */\n"
  "        }\n"
  "      }\n"
  "\n"
  "      /* printf(\"Initializing arrays: Syn= %d bytes,Queue= %d bytes.\\n\",\n"
  "	 (int)msn*sizeof(SynS),(int)qlen, sizeof(QueU)); */\n"
  "      QPRCAST->qq = QUECAST;\n"
  "    }\n"
  "  }\n"
  " ENDVERBATIM\n"
  "\n"
  "}\n"
  "\n"
  "VERBATIM\n"
  "static void hshake();\n"
  "ENDVERBATIM\n"
  ": 2 arguments - index, presyn link\n"
  "PROCEDURE setlink() {\n"
  "VERBATIM { \n"
  "  int ii, x, new;\n"
  "  SynS *sns;\n"
  "  double ptemp;\n"
  "  PreL *ppsyn;\n"
  "    \n"
  "  if (ifarg(2)) { \n"
  "    x = (int)*getarg(1);\n"
  "    ptemp = *getarg(2);		/* pick up pointer */\n"
  "  } else {\n"
  "    x = nsyn;\n"
  "    ptemp = *getarg(1);\n"
  "  }\n"
  "\n"
  "  if (x >= maxsyn) {\n"
  "    init_arrays(maxsyn+POSTINC); /* #DEFINE POSTINC 5 */\n"
  "  } \n"
  "  if (x > nsyn) {\n"
  "    hoc_execerror(\"Can't leave empty pointers: see nsyn for current array index.\", 0);\n"
  "  }\n"
  "\n"
  "  sns = &(PSTCAST[x*(int)CHAINLEN]);	/* postsynaptic entry */\n"
  "  if ((ppsyn = (PreL *)((unsigned long)ptemp)) == (PreL *)NULL) { /* presyn loc */\n"
  "    hoc_execerror(\"Presyn not initialized.\", 0); }\n"
  "  if (ppsyn->link2 != ptemp) {\n"
  "    hoc_execerror(\"Invalid (changed) link.\", 0); }\n"
  "\n"
  "  if (x == nsyn) {\n"
  "    nsyn++;		/* a new entry */\n"
  "    new = 1;\n"
  "  } else { \n"
  "    /* should generate error if try to change something in middle of chain */\n"
  "    if (sns->chainlen == -2) {\n"
  "      hoc_execerror(\"Internal error: Index used must be multiple of CHAINLEN.\", 0);\n"
  "    }\n"
  "    new = 0; \n"
  "  }\n"
  "\n"
  "  hshake(sns,ppsyn,new);\n"
  "\n"
  "  x *= (int)CHAINLEN;\n"
  "  for (ii=x;ii < x + CHAINLEN;ii++) {\n"
  "    sns = &(PSTCAST[ii]);	/* access the entry */\n"
  "    /* initialize postsyn stuff to reasonable values */    \n"
  "    /* delay, gmax and all codes should be initialized separately */\n"
  "    sns->del = DELAY;\n"
  "    sns->ucode = -1;\n"
  "    sns->chainlen = (ii==x)?CHAINLEN:-2;\n"
  "    sns->pgm = 1.0;\n"
  "\n"
  "    /* initialize the presynaptic stuff */\n"
  "    /* Rcurr and last will be initialized with init */\n"
  "    sns->index = x;\n"
  "    sns->qpt = QPRCAST; \n"
  "  }\n"
  "}\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  ": manipulate the presynaptic list remotely\n"
  "VERBATIM\n"
  "static void hshake(ss,pl,flag) \n"
  "/* ls will be a pointer to presyn cell's array of pointers */\n"
  "SynS *ss;\n"
  "PreL *pl;\n"
  "int flag; /* flag == 1 if this is a brand new entry */\n"
  "{\n"
  "  int ii;\n"
  "  double *nn;\n"
  "  double *mx;\n"
  "\n"
  "  nn = pl->nnpre; /* pointer to presyn npre */\n"
  "  mx = pl->mxpre; /* pointer to presyn maxpre */\n"
  "\n"
  "  /* erase presyn pointer if this has been set before */\n"
  "  if (flag == 0) {		/* an old entry */\n"
  "    /* fall out of loop if a pointer exists already */\n"
  "    for (ii=0;ii<(*nn) && (pl->plst[ii])!=ss;ii++); /* no body */\n"
  "    if (ii < (*nn)) {	/* fell out of loop */\n"
  "      printf(\"Maintaining pointers.\\n\");\n"
  "    } else {\n"
  "      flag = 1;			/* now we do have to create a new pointer */\n"
  "      printf(\"Erasing pointer (C%09d,Pr%09d,Po%09d).\\n\",\n"
  "	     ss->ucode,*(ss->pcpre),QPRCAST->cpost);\n"
  "      *(ss->back) = (SynS *)NULL; /* remove old pre-post pointer */\n"
  "    }\n"
  "  }\n"
  "\n"
  "  if (flag == 1) {		/* create a new pointer */\n"
  "    /* if necessary, allocate space for the pointer system */\n"
  "    if ((*nn) == 0.) { \n"
  "      (*mx) = PREINC;  /* #DEFINE PREINC 50 */\n"
  "      pl->plst = (SynS **)ecalloc(((int)(*mx)), sizeof(SynS *)); \n"
  "    } else if ((*nn) == (*mx)) { /* need to create more space */\n"
  "      (*mx) += PREINC;\n"
  "      pl->plst = (SynS **)realloc(pl->plst,((int)(*mx)) * sizeof(SynS *)); \n"
  "      /* reassign all back pointers */\n"
  "      for (ii=0;ii<(*nn);ii++) {\n"
  "	(pl->plst[ii])->back = &(pl->plst[ii]);\n"
  "      }\n"
  "    }\n"
  "\n"
  "    /* the handshake */\n"
  "    pl->plst[(int)(*nn)] = ss;	/* pre  -> post */\n"
  "    ss->back = &(pl->plst[(int)(*nn)]); /* post -> pre */\n"
  "    ss->pcpre = &(pl->cpre); /* post -> precode */\n"
  "    (*nn)++;\n"
  "  }\n"
  "}\n"
  "ENDVERBATIM\n"
  "\n"
  ": 1 or 2 args, get bzw. set user code\n"
  "FUNCTION code() {\n"
  "  VERBATIM { int ii;\n"
  "    if (ifarg(1)) {\n"
  "      ii = (int)*getarg(1);\n"
  "      if (ii < 0) { ii = nsyn+ii; }\n"
  "      if (ii >= nsyn || ii < 0) {hoc_execerror(\"array index out of bounds\", 0);}\n"
  "      if (ifarg(2)) {  /* look for a second arg to do a set */\n"
  "	(PSTCAST[ii*(int)CHAINLEN]).ucode = (int)*getarg(2); }\n"
  "      _lcode = (PSTCAST[ii*(int)CHAINLEN]).ucode;\n"
  "    } else {\n"
  "      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) {\n"
  "	  printf(\"%d,%09d  \",ii,(PSTCAST[ii]).ucode);\n"
  "	}\n"
  "      _lcode = 1.0;\n"
  "    }}      \n"
  "  ENDVERBATIM\n"
  "}\n"
  " \n"
  ": 2 args => set postsyn code\n"
  "FUNCTION post() { \n"
  "  VERBATIM { \n"
  "    if (ifarg(2)) { QPRCAST->cpost = (int)*getarg(2); }\n"
  "    _lpost = QPRCAST->cpost;\n"
  "  }\n"
  "  ENDVERBATIM\n"
  "}\n"
  " \n"
  ": 1 or 2 args, get bzw. set delay\n"
  "FUNCTION delay() {\n"
  "  VERBATIM { int ii,jj;\n"
  "    if (ifarg(1)) {\n"
  "      ii = (int)*getarg(1);\n"
  "      if (ii < 0) { ii = nsyn+ii; }\n"
  "      if (ii >= nsyn || ii < 0) {hoc_execerror(\"array index out of bounds\", 0);}\n"
  "      if (ifarg(2)) {  /* look for a second arg to do a set */\n"
  "	for (jj=ii*CHAINLEN;jj<(ii+1)*CHAINLEN;jj++) {\n"
  " 	  (PSTCAST[jj]).del = *getarg(2) + DELAY; }}\n"
  "      _ldelay = (PSTCAST[ii*(int)CHAINLEN]).del;\n"
  "    } else {\n"
  "      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) {\n"
  "	  printf(\"%d,%g  \",ii,(PSTCAST[ii]).del);\n"
  "	}\n"
  "      _ldelay = 1.0;\n"
  "    }}\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ": 1 or 2 args, get bzw. set individual percent gmax\n"
  ": user expects to see a gmax in uS, but what is stored is a \n"
  ": fraction of the total gmax (a global variable)\n"
  ": therefore can turn off all GABAA with gmax_GABAA = 0\n"
  "FUNCTION gmax() { \n"
  "  VERBATIM { int ii,jj;\n"
  "    if (ifarg(1)) {\n"
  "      ii = (int)*getarg(1);\n"
  "      if (ii < 0) { ii = nsyn+ii; }\n"
  "      if (ii >= nsyn || ii < 0) {hoc_execerror(\"array index out of bounds\", 0);}\n"
  "      if (ifarg(2)) {  /* look for a second arg to do a set */\n"
  " 	for (jj=ii*CHAINLEN;jj<(ii+1)*CHAINLEN;jj++) {\n"
  " 	  (PSTCAST[jj]).pgm = *getarg(2); }}\n"
  "      _lgmax = (PSTCAST[ii*(int)CHAINLEN]).pgm;\n"
  "    } else {\n"
  "      printf(\"Multiply by %g to get effective gmax.\\n\",GMAX);\n"
  "      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) {\n"
  "	  printf(\"%d,%g  \",ii,(PSTCAST[ii]).pgm);\n"
  "	}\n"
  "      _lgmax = 1.0;\n"
  "    }}\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ": called with 0,1 or 2 args\n"
  ": 0 => print out info about all the presyns\n"
  ": 1 (index) => use as index to return this pre code\n"
  ":              -1 => print out pointer structure\n"
  ": 2 (index,code )=> set precode for this index\n"
  "FUNCTION pre() {\n"
  "  VERBATIM { int x,ii;\n"
  "    x = -2;  /* flag -> -1 then print out pointers */\n"
  "    if (ifarg(1)) { x = (int)*getarg(1); }\n"
  "    if (x >= 0) {\n"
  "      if (x >= nsyn) { hoc_execerror(\"array index out of bounds\", 0);}\n"
  "      x *= CHAINLEN;\n"
  "      _lpre = *((PSTCAST[x]).pcpre);\n"
  "    } else {	  \n"
  "      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) {\n"
  "	printf(\"%1d Del:%-2g Gmax: %-4g Cd:%09d Pre:%09d Post:%09d\\n\",\n"
  "	       (PSTCAST[ii]).index,\n"
  "	       (PSTCAST[ii]).del,\n"
  "	       (PSTCAST[ii]).pgm,\n"
  "	       (PSTCAST[ii]).ucode,\n"
  "	       *((PSTCAST[ii]).pcpre),\n"
  "	       QPRCAST->cpost);\n"
  "	if (x == -1 || &(PSTCAST[ii]) != *((PSTCAST[ii]).back)) {\n"
  "	  printf(\"\\t\\t%s: %x -> %x\\n\",\n"
  "		 ((&(PSTCAST[ii])==*((PSTCAST[ii]).back))?\"OK\":\"POINTER MISMATCH\"),\n"
  "		 &(PSTCAST[ii]),\n"
  "		 (PSTCAST[ii]).back);\n"
  "	}\n"
  "      }\n"
  "      _lpre = 1.0;\n"
  "    }\n"
  "  } \n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ": return link for matching with presyn\n"
  "FUNCTION link(ii) { \n"
  "  VERBATIM { \n"
  "    if (_lii >= nsyn) { hoc_execerror(\"array index out of bounds\", 0);}\n"
  "    _llink = (double)(unsigned long)(PSTCAST[(int)_lii]).pcpre;\n"
  "  }\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  "FUNCTION check() {\n"
  "  VERBATIM { \n"
  "    int ii;\n"
  "    for (ii=0; ii<nsyn*CHAINLEN && _lcheck==1.; ii+=CHAINLEN) {\n"
  "      if (&(PSTCAST[ii]) != *((PSTCAST[ii]).back)) {\n"
  "	printf(\"****************************************************************\\n\");\n"
  "	printf(\"ERROR:: Index:%3d,C%09d,Pr%09d,Po%09d, Delay:%6g (%x->%x->%x)\\n\",\n"
  "	       (PSTCAST[ii]).index,\n"
  "	       (PSTCAST[ii]).ucode,\n"
  "	       *((PSTCAST[ii]).pcpre),\n"
  "	       QPRCAST->cpost,\n"
  "	       (PSTCAST[ii]).del,\n"
  "	       &(PSTCAST[ii]),\n"
  "	       (PSTCAST[ii]).back,\n"
  "	       *((PSTCAST[ii]).back));\n"
  "	printf(\"****************************************************************\\n\");\n"
  "	_lcheck = -1.;\n"
  "      }\n"
  "    }\n"
  "    _lcheck = nsyn;  /* all pointers have to be active */\n"
  "  }\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  "PROCEDURE initq() { \n"
  "  VERBATIM {\n"
  "  int ii;\n"
  "  QPRCAST->dead = Cdur + Deadtime;\n"
  "  begsyn = endsyn = newspk = 0.;\n"
  "  for (ii=0;ii<qlen;ii++) {\n"
  "    QUECAST[ii].time = 1.e20;\n"
  "    QUECAST[ii].index = -1;\n"
  "  }\n"
  "  QPRCAST->qterm = -1.e2;\n"
  "  for (ii=0;ii<nsyn*CHAINLEN;ii++) {\n"
  "    PSTCAST[ii].Rcurr = 0.;\n"
  "    PSTCAST[ii].last = -1.e2;\n"
  "    PSTCAST[ii].spkt = -1.e2;\n"
  "    if (PSTCAST[ii].chainlen > 1) { PSTCAST[ii].chainptr = -1; }\n"
  "  }\n"
  "  }\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ": begsyn is the top of the line when the spike is first utilized to generate synapse\n"
  ": updates the value in the queue so endsyn can detect the end of the syn pulse\n"
  ": doesn't return value (since in a struct)\n"
  "PROCEDURE popqh1(aug) { \n"
  "  VERBATIM {  \n"
  "    if (QUECAST[(int)begsyn].time == 1e20) {\n"
  "      printf(\"%x %g %g \",QUECAST,newspk,qlen);\n"
  "      hoc_execerror(\"Error: queue exhausted.\\n\",0);\n"
  "    } else {  /* augment the time by Cdur */\n"
  "      QUECAST[(int)begsyn].time += _laug;\n"
  "      begsyn++;\n"
  "      if (begsyn == qlen) { begsyn = 0.; }\n"
  "    }\n"
  "  }\n"
  "  ENDVERBATIM\n"
  "}\n"
  " \n"
  ": endsyn is queried for the termination time of the synaptic square wave\n"
  "PROCEDURE popqh2() { \n"
  "  VERBATIM {  \n"
  "    QUECAST[(int)endsyn].time = 1.e20;  /* clear the entry */\n"
  "    endsyn++;\n"
  "    if (endsyn == qlen) { endsyn = 0.; }\n"
  "  }\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ": DEBUGGING ROUTINES\n"
  "\n"
  "FUNCTION getdbx(c,x) {\n"
  "VERBATIM { \n"
  "  switch ((int)_lc) {\n"
  "  case 1:\n"
  "    _lgetdbx = (PSTCAST[(int)_lx]).last;\n"
  "    break;\n"
  "  case 2:\n"
  "    _lgetdbx = (PSTCAST[(int)_lx]).Rcurr;\n"
  "    break;\n"
  "  case 3:\n"
  "    _lgetdbx = (PSTCAST[(int)_lx]).spkt;\n"
  "    break;\n"
  "  case 4:\n"
  "    _lgetdbx = (QUECAST[(int)_lx]).time;\n"
  "    break;\n"
  "  case 5:\n"
  "    _lgetdbx = (QUECAST[(int)_lx]).index;\n"
  "    break;\n"
  "  case 6:\n"
  "    _lgetdbx = QPRCAST->qterm;\n"
  "    break;\n"
  "  case 7:\n"
  "    _lgetdbx = QPRCAST->dead;\n"
  "    break;\n"
  "  case 8:\n"
  "    _lgetdbx = (PSTCAST[(int)_lx]).chainlen;\n"
  "    break;\n"
  "  case 9:\n"
  "    _lgetdbx = (PSTCAST[(int)_lx]).chainptr;\n"
  "    break;\n"
  "  default:\n"
  "    hoc_execerror(\"UNAVAILABLE IN DBX\",0);\n"
  "    _lgetdbx = -1;\n"
  "    break;\n"
  "  }}\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  "PROCEDURE prq() {\n"
  "  VERBATIM { int ii;\n"
  "    printf(\"new:%g beg:%g end:%g /%g\\n\",newspk,begsyn,endsyn,qlen);\n"
  "    for (ii=endsyn;ii!=newspk;ii=((ii==qlen-1)?0:ii+1)) {\n"
  "        printf(\"%d %8g%8d\\n\",ii,(QUECAST[ii]).time,(QUECAST[ii]).index);\n"
  "    }\n"
  "  }\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  "\n"
  ":::end INCLUDE snsarr.inc\n"
  "\n"
  "UNITS {\n"
  "  (nA) = (nanoamp)\n"
  "  (mV) = (millivolt)\n"
  "  (umho) = (micromho)\n"
  "  (mM) = (milli/liter)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "  celsius\n"
  "  q10 = 1.\n"
  "  exptemp = 37.\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "  v		(mV)		: postsynaptic voltage\n"
  "  i 		(nA)		: current = g*(v - Erev)\n"
  "  g 		(umho)		: conductance\n"
  "  R				: fraction of open channels, Ron + Roff\n"
  "  Ron                           : activation state while syn's turned on\n"
  "  Roff                          : activation state for decaying syns\n"
  "  Rinf				: steady state channels open\n"
  "  Rtau		(ms)		: time constant of channel binding\n"
  "  synon                         : number of syns turned on at a time\n"
  "  drive                         : drive for ODE toward Rinf\n"
  "  qq10                          : rate speed-up due to Q10\n"
  "  edt                           : rate factor for Ron\n"
  "  edb                           : decay factor for Roff\n"
  "  edc                           : rate factor for increasing Rcurr\n"
  "  dt\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "  : initialize GLOBAL parameters, in case user changes one\n"
  "  : this repeats unnecessarily for each instance\n"
  "  Rinf = Alpha / (Alpha + Beta)\n"
  "  qq10 = q10^((celsius-exptemp)/10.)\n"
  "  Rtau = 1 / (Alpha + Beta) / qq10\n"
  "  edt = exp(-dt/Rtau)\n"
  "  edb = exp(-Beta * dt)\n"
  "  edc = exp(-Cdur/Rtau)\n"
  "  drive = Rinf * (1. - edt)\n"
  "\n"
  "  : initialize RANGE parameters\n"
  "  synon = 0\n"
  "  R = 0\n"
  "  Ron = 0\n"
  "  Roff = 0\n"
  "\n"
  "  : do not initialize QUEU if it hasn't been allocated by init_arrays()\n"
  "  if (nsyn > 0) {\n"
  "    initq()   \n"
  "  } else {\n"
  "    VERBATIM\n"
  "    printf(\"WARNING nsyn==0 \");\n"
  "    ENDVERBATIM\n"
  "  }\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    SOLVE release\n"
  "    R = Ron + Roff\n"
  "    g = GMAX * R\n"
  "    i = g*(v - Erev)\n"
  "}\n"
  "\n"
  "PROCEDURE release() { \n"
  "  if (nsyn>0) { : do not try accessing a queue that hasn't been allocated\n"
  "  VERBATIM \n"
  "  int who;\n"
  "  QueU *pqueu;\n"
  "  SynS *ppst;\n"
  "\n"
  "  pqueu = (QueU *)((unsigned long) queu);\n"
  "\n"
  "  while (t >= pqueu[(int)begsyn].time) { /*  somebody spiked delay time ago */\n"
  "    ppst = (SynS *)((unsigned long) lpst);\n"
  "    who = pqueu[(int)begsyn].index; /* who spiked? */\n"
  "    /* calculate the decay that occurred since last activity ended */\n"
  "    ppst[who].Rcurr *= exptable(-Beta*(t-ppst[who].last));\n"
  "    /* transfer the value from Roff to Ron */\n"
  "    Ron += ppst[who].Rcurr;\n"
  "    Roff -= ppst[who].Rcurr;\n"
  "    synon += ppst[who].pgm;	/*  weight syn by percent gmax */\n"
  "    ppst[who].last = t + Cdur;   /* time when syn will turn off */\n"
  "    popqh1(Cdur);		/* next (also add Cdur to value on the queu) */\n"
  "  }\n"
  "\n"
  "  while (t >= pqueu[(int)endsyn].time) { /*  somebody needs to be turned off */\n"
  "    ppst = (SynS *)((unsigned long) lpst); /* will do this assign twice in rare cases */\n"
  "    who = pqueu[(int)endsyn].index;   /* who spiked? */\n"
  "    /* solve Rcurr forward in time till end of syn activity */\n"
  "    ppst[who].Rcurr = ppst[who].pgm*Rinf*(1.-edc) + ppst[who].Rcurr*edc;\n"
  "    Ron -= ppst[who].Rcurr;\n"
  "    Roff += ppst[who].Rcurr;  /* transfer from on to off */\n"
  "    synon -= ppst[who].pgm;   /* remove this percent of gmax */\n"
  "    if (synon<1e-11 && synon>-1e-11) { synon=0.; }\n"
  "    if (synon==0. || Ron < 0.) { Ron = 0.; }\n"
  "    popqh2();  /* next */\n"
  "  }\n"
  "\n"
  "  /*  update R */\n"
  "  if (synon > 0) {		/*  one or more synapses turned on? */\n"
  "    Ron = synon*drive + Ron * edt; /*  drive R toward Rinf */\n"
  "  } \n"
  "  Roff *= edb;			/*  Roff always decays toward 0 */\n"
  "  return 0;\n"
  "  ENDVERBATIM\n"
  "  }\n"
  "}\n"
  "\n"
  "FUNCTION getR(x) {\n"
  "  VERBATIM {\n"
  "    SynS pst;\n"
  "    double rnow;\n"
  "    double end, rinf;\n"
  "\n"
  "    pst = (PSTCAST[(int)_lx]);\n"
  "    end = pst.last;\n"
  "    rinf = pst.pgm * Rinf;\n"
  "\n"
  "    if (end < 0.) {		/* not yet turned on */\n"
  "      rnow = 0.;\n"
  "    } else if (t >= end - dt/2) {	/* decay */\n"
  "      rnow = pst.Rcurr * exptable(-Beta*(t-end));\n"
  "    } else {			/* turning on */\n"
  "      rnow = rinf + (pst.Rcurr - rinf)*exptable((t-(end-Cdur))/(-Rtau));\n"
  "    }\n"
  "    if (pst.pgm != 0.) {\n"
  "      _lgetR = rnow/pst.pgm;  /* scale it to 1 so it looks like a state variable */\n"
  "    } else {\n"
  "      _lgetR = 0.;\n"
  "    }\n"
  "  }\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ": only gets called for negative numbers\n"
  "FUNCTION exptable(x) { \n"
  "  TABLE  FROM -10 TO 0 WITH 1000\n"
  "  \n"
  "  if ((x > -10) && (x < 0)) {\n"
  "    exptable = exp(x)\n"
  "  } else {\n"
  "    exptable = 0.\n"
  "  }\n"
  "}\n"
  "\n"
  ":::end INCLUDE sns.inc\n"
  ":** AMPA\n"
  ;
#endif
