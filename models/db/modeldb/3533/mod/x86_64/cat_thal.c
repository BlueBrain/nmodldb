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
 
#define nrn_init _nrn_init__cat
#define _nrn_initial _nrn_initial__cat
#define nrn_cur _nrn_cur__cat
#define _nrn_current _nrn_current__cat
#define nrn_jacob _nrn_jacob__cat
#define nrn_state _nrn_state__cat
#define _net_receive _net_receive__cat 
#define states states__cat 
#define trates trates__cat 
 
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
#define gbar _p[0]
#define m _p[1]
#define h _p[2]
#define cai _p[3]
#define cao _p[4]
#define ica _p[5]
#define Dm _p[6]
#define Dh _p[7]
#define _g _p[8]
#define _ion_cai	*_ppvar[0]._pval
#define _ion_cao	*_ppvar[1]._pval
#define _ion_ica	*_ppvar[2]._pval
#define _ion_dicadv	*_ppvar[3]._pval
 
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
 static void _hoc_alph(void);
 static void _hoc_alpm(void);
 static void _hoc_beth(void);
 static void _hoc_betm(void);
 static void _hoc_trates(void);
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
 "setdata_cat", _hoc_setdata,
 "alph_cat", _hoc_alph,
 "alpm_cat", _hoc_alpm,
 "beth_cat", _hoc_beth,
 "betm_cat", _hoc_betm,
 "trates_cat", _hoc_trates,
 0, 0
};
#define alph alph_cat
#define alpm alpm_cat
#define beth beth_cat
#define betm betm_cat
 extern double alph( double );
 extern double alpm( double );
 extern double beth( double );
 extern double betm( double );
 /* declare global and static user variables */
#define a0h a0h_cat
 double a0h = 0.01;
#define a0m a0m_cat
 double a0m = 0.1;
#define carev carev_cat
 double carev = 0;
#define gmh gmh_cat
 double gmh = 0.5;
#define gmm gmm_cat
 double gmm = 0.5;
#define htau htau_cat
 double htau = 0;
#define hinf hinf_cat
 double hinf = 0;
#define mtau mtau_cat
 double mtau = 0;
#define minf minf_cat
 double minf = 0;
#define q10 q10_cat
 double q10 = 3;
#define vhalfh vhalfh_cat
 double vhalfh = 10;
#define vhalfm vhalfm_cat
 double vhalfm = 55;
#define vrest vrest_cat
 double vrest = -75;
#define zetah zetah_cat
 double zetah = 0.1;
#define zetam zetam_cat
 double zetam = 0.1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "minf_cat", "mV",
 "mtau_cat", "ms",
 "hinf_cat", "mV",
 "htau_cat", "ms",
 "carev_cat", "mV",
 "gbar_cat", "mho/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "vrest_cat", &vrest_cat,
 "a0m_cat", &a0m_cat,
 "vhalfm_cat", &vhalfm_cat,
 "zetam_cat", &zetam_cat,
 "gmm_cat", &gmm_cat,
 "a0h_cat", &a0h_cat,
 "vhalfh_cat", &vhalfh_cat,
 "zetah_cat", &zetah_cat,
 "gmh_cat", &gmh_cat,
 "q10_cat", &q10_cat,
 "minf_cat", &minf_cat,
 "mtau_cat", &mtau_cat,
 "hinf_cat", &hinf_cat,
 "htau_cat", &htau_cat,
 "carev_cat", &carev_cat,
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
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"cat",
 "gbar_cat",
 0,
 0,
 "m_cat",
 "h_cat",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 9, _prop);
 	/*initialize range parameters*/
 	gbar = 0.007;
 	_prop->param = _p;
 	_prop->param_size = 9;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* cao */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _cat_thal_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 9, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 cat /Users/awile/workbench/modldb/nmodldb/models/db/modeldb/3533/mod/cat_thal.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 
#define FARADAY _nrnunit_FARADAY[_nrnunit_use_legacy_]
static double _nrnunit_FARADAY[2] = {0x1.78e555060882cp+16, 96485.3}; /* 96485.3321233100141 */
 
#define R _nrnunit_R[_nrnunit_use_legacy_]
static double _nrnunit_R[2] = {0x1.0a1013e8990bep+3, 8.3145}; /* 8.3144626181532395 */
static int _reset;
static char *modelname = "cat";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int trates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   trates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 trates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   trates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0 ) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0 ) ) ) / htau ) - h) ;
   }
  return 0;
}
 
static int  trates (  double _lv ) {
   double _lqt ;
 _lqt = pow( q10 , ( ( celsius - 22.0 ) / 10.0 ) ) ;
   minf = 1.0 / ( 1.0 + exp ( - ( _lv - ( vrest + 52.9 ) ) / 13.6 ) ) ;
   mtau = betm ( _threadargscomma_ _lv ) / ( _lqt * a0m * ( 1.0 + alpm ( _threadargscomma_ _lv ) ) ) ;
   hinf = 1.0 / ( 1.0 + exp ( ( _lv - ( vrest + 4.5 ) ) / 11.8 ) ) ;
   htau = beth ( _threadargscomma_ _lv ) / ( _lqt * a0h * ( 1.0 + alph ( _threadargscomma_ _lv ) ) ) ;
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
   _r = 1.;
 trates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double alpm (  double _lv ) {
   double _lalpm;
 _lalpm = exp ( zetam * ( _lv - vrest - vhalfm ) ) ;
   
return _lalpm;
 }
 
static void _hoc_alpm(void) {
  double _r;
   _r =  alpm (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double betm (  double _lv ) {
   double _lbetm;
 _lbetm = exp ( zetam * gmm * ( _lv - vrest - vhalfm ) ) ;
   
return _lbetm;
 }
 
static void _hoc_betm(void) {
  double _r;
   _r =  betm (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double alph (  double _lv ) {
   double _lalph;
 _lalph = exp ( zetah * ( _lv - vrest - vhalfh ) ) ;
   
return _lalph;
 }
 
static void _hoc_alph(void) {
  double _r;
   _r =  alph (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double beth (  double _lv ) {
   double _lbeth;
 _lbeth = exp ( zetah * gmh * ( _lv - vrest - vhalfh ) ) ;
   
return _lbeth;
 }
 
static void _hoc_beth(void) {
  double _r;
   _r =  beth (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
 {
   trates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
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
  cai = _ion_cai;
  cao = _ion_cao;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   carev = ( 1e3 ) * ( R * ( celsius + 273.15 ) ) / ( 2.0 * FARADAY ) * log ( cao / cai ) ;
   ica = gbar * pow( m , 3.0 ) * h * ( v - carev ) ;
   }
 _current += ica;

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
  cai = _ion_cai;
  cao = _ion_cao;
 _g = _nrn_current(_v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
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
  cai = _ion_cai;
  cao = _ion_cao;
 { error =  states();
 if(error){fprintf(stderr,"at line 59 in file cat_thal.mod:\n        SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/awile/workbench/modldb/nmodldb/models/db/modeldb/3533/mod/cat_thal.mod";
static const char* nmodl_file_text = 
  "TITLE cat\n"
  ": low-threshold Ca current for TC neurons from Williams and Stuart (2000).\n"
  ": M.Migliore Jan. 2002\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX cat\n"
  "	USEION ca READ cai,cao WRITE ica\n"
  "        RANGE gbar,cai\n"
  "	GLOBAL minf, hinf, mtau, htau, carev\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "	(pS) = (picosiemens)\n"
  "	(um) = (micron)\n"
  "	(molar) = (1/liter)\n"
  "	(mM) = (millimolar)\n"
  "	FARADAY = (faraday) (coulomb)\n"
  "	R = (k-mole) (joule/degC)\n"
  "} \n"
  "\n"
  "PARAMETER {\n"
  "	gbar = 0.007   	(mho/cm2)	\n"
  "								\n"
  "	celsius = 22 (degC)\n"
  "	v 		(mV)\n"
  "	vrest=-75\n"
  "	a0m=0.1\n"
  "	vhalfm=55\n"
  "	zetam=0.1\n"
  "	gmm=0.5\n"
  "\n"
  "	a0h=0.01\n"
  "	vhalfh=10\n"
  "	zetah=0.1\n"
  "	gmh=0.5\n"
  "	\n"
  "	q10=3\n"
  "	cai 	= .000050 (mM)		: initial [Ca]i = 50 nM\n"
  "	cao 	= 2	(mM)		: [Ca]o = 2 mM\n"
  "\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "	ica 		(mA/cm2)\n"
  "	minf (mV)	mtau (ms)	 	\n"
  "	hinf (mV)	htau (ms)	 	\n"
  "	carev	(mV)\n"
  "}\n"
  " \n"
  "\n"
  "\n"
  "STATE { m h}\n"
  "\n"
  "BREAKPOINT {\n"
  "        SOLVE states METHOD cnexp\n"
  "	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)\n"
  "	ica = gbar*m^3*h* (v - carev)\n"
  "} \n"
  "\n"
  "INITIAL {\n"
  "	trates(v)\n"
  "	m=minf  \n"
  "	h=hinf  \n"
  "}\n"
  "\n"
  "DERIVATIVE states {   \n"
  "        trates(v)      \n"
  "        m' = (minf-m)/mtau\n"
  "        h' = (hinf-h)/htau\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE trates(v) {  \n"
  "	LOCAL qt\n"
  "        qt=q10^((celsius-22)/10)\n"
  "        minf = 1/(1 + exp(-(v-(vrest+52.9))/13.6))\n"
  "	mtau = betm(v)/(qt*a0m*(1+alpm(v)))\n"
  "\n"
  "        hinf = 1/(1 + exp((v-(vrest+4.5))/11.8))\n"
  "	htau = beth(v)/(qt*a0h*(1+alph(v)))\n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION alpm(v(mV)) {\n"
  "  alpm = exp(zetam*(v-vrest-vhalfm)) \n"
  "}\n"
  "\n"
  "FUNCTION betm(v(mV)) {\n"
  "  betm = exp(zetam*gmm*(v-vrest-vhalfm)) \n"
  "}\n"
  "\n"
  "FUNCTION alph(v(mV)) {\n"
  "  alph = exp(zetah*(v-vrest-vhalfh)) \n"
  "}\n"
  "\n"
  "FUNCTION beth(v(mV)) {\n"
  "  beth = exp(zetah*gmh*(v-vrest-vhalfh)) \n"
  "}\n"
  ;
#endif
