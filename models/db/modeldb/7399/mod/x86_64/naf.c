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
 
#define nrn_init _nrn_init__naf
#define _nrn_initial _nrn_initial__naf
#define nrn_cur _nrn_cur__naf
#define _nrn_current _nrn_current__naf
#define nrn_jacob _nrn_jacob__naf
#define nrn_state _nrn_state__naf
#define _net_receive _net_receive__naf 
#define _f_mh _f_mh__naf 
#define iassign iassign__naf 
#define mh mh__naf 
#define states states__naf 
 
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
#define gmax _p[0]
#define i _p[1]
#define g _p[2]
#define m _p[3]
#define h _p[4]
#define ina _p[5]
#define Dm _p[6]
#define Dh _p[7]
#define _g _p[8]
#define _ion_ina	*_ppvar[0]._pval
#define _ion_dinadv	*_ppvar[1]._pval
 
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
 static void _hoc_FRT(void);
 static void _hoc_alpha(void);
 static void _hoc_beta(void);
 static void _hoc_ghkca(void);
 static void _hoc_iassign(void);
 static void _hoc_mh(void);
 static void _hoc_states(void);
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
 "setdata_naf", _hoc_setdata,
 "FRT_naf", _hoc_FRT,
 "alpha_naf", _hoc_alpha,
 "beta_naf", _hoc_beta,
 "ghkca_naf", _hoc_ghkca,
 "iassign_naf", _hoc_iassign,
 "mh_naf", _hoc_mh,
 "states_naf", _hoc_states,
 0, 0
};
#define FRT FRT_naf
#define alpha alpha_naf
#define beta beta_naf
#define ghkca ghkca_naf
 extern double FRT( double );
 extern double alpha( double , double );
 extern double beta( double , double );
 extern double ghkca( double );
 /* declare global and static user variables */
#define Add Add_naf
 double Add[2];
#define Inf Inf_naf
 double Inf[2];
#define Mult Mult_naf
 double Mult[2];
#define Tau Tau_naf
 double Tau[2];
#define cai cai_naf
 double cai = 0;
#define cao cao_naf
 double cao = 0;
#define erev erev_naf
 double erev = 60;
#define hexp hexp_naf
 double hexp = 1;
#define hq10 hq10_naf
 double hq10 = 3;
#define htemp htemp_naf
 double htemp = 37;
#define hbasetau hbasetau_naf
 double hbasetau = 0.25;
#define hvhalf hvhalf_naf
 double hvhalf = -39;
#define hbaserate hbaserate_naf
 double hbaserate = 0.095;
#define hgamma hgamma_naf
 double hgamma = 0.3;
#define hvalence hvalence_naf
 double hvalence = -6;
#define mexp mexp_naf
 double mexp = 3;
#define mq10 mq10_naf
 double mq10 = 3;
#define mtemp mtemp_naf
 double mtemp = 37;
#define mbasetau mbasetau_naf
 double mbasetau = 0.02;
#define mvhalf mvhalf_naf
 double mvhalf = -33.5;
#define mbaserate mbaserate_naf
 double mbaserate = 4.5;
#define mgamma mgamma_naf
 double mgamma = 0.5;
#define mvalence mvalence_naf
 double mvalence = 4;
#define usetable usetable_naf
 double usetable = 1;
#define vmin vmin_naf
 double vmin = -100;
#define vmax vmax_naf
 double vmax = 100;
#define vrest vrest_naf
 double vrest = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_naf", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "erev_naf", "mV",
 "cao_naf", "mM",
 "cai_naf", "mM",
 "vmax_naf", "mV",
 "vmin_naf", "mV",
 "gmax_naf", "umho",
 "i_naf", "mA/cm^2",
 "g_naf", "mho/cm^2",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "erev_naf", &erev_naf,
 "vrest_naf", &vrest_naf,
 "mvalence_naf", &mvalence_naf,
 "mgamma_naf", &mgamma_naf,
 "mbaserate_naf", &mbaserate_naf,
 "mvhalf_naf", &mvhalf_naf,
 "mbasetau_naf", &mbasetau_naf,
 "mtemp_naf", &mtemp_naf,
 "mq10_naf", &mq10_naf,
 "mexp_naf", &mexp_naf,
 "hvalence_naf", &hvalence_naf,
 "hgamma_naf", &hgamma_naf,
 "hbaserate_naf", &hbaserate_naf,
 "hvhalf_naf", &hvhalf_naf,
 "hbasetau_naf", &hbasetau_naf,
 "htemp_naf", &htemp_naf,
 "hq10_naf", &hq10_naf,
 "hexp_naf", &hexp_naf,
 "cao_naf", &cao_naf,
 "cai_naf", &cai_naf,
 "vmax_naf", &vmax_naf,
 "vmin_naf", &vmin_naf,
 "usetable_naf", &usetable_naf,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 "Inf_naf", Inf_naf, 2,
 "Tau_naf", Tau_naf, 2,
 "Mult_naf", Mult_naf, 2,
 "Add_naf", Add_naf, 2,
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"naf",
 "gmax_naf",
 0,
 "i_naf",
 "g_naf",
 0,
 "m_naf",
 "h_naf",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 9, _prop);
 	/*initialize range parameters*/
 	gmax = 0.03;
 	_prop->param = _p;
 	_prop->param_size = 9;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 2, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 	_ppvar[0]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[1]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _naf_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 9, 2);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 naf /Users/awile/workbench/modldb/nmodldb/models/db/modeldb/7399/mod/naf.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96489.0;
 static double R = 8.31441;
 static double *_t_Add[2];
 static double *_t_Mult[2];
static int _reset;
static char *modelname = "Borg-Graham Channel Model";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_mh(double);
static int iassign();
static int mh(double);
static int states();
 static void _n_mh(double);
 
static int  states (  ) {
   mh ( _threadargscomma_ v * 1.0 ) ;
   m = m * Mult [ 0 ] + Add [ 0 ] ;
   h = h * Mult [ 1 ] + Add [ 1 ] ;
   
/*VERBATIM*/
	return 0;
  return 0; }
 
static void _hoc_states(void) {
  double _r;
   _r = 1.;
 states (  );
 hoc_retpushx(_r);
}
 static double _mfac_mh, _tmin_mh;
 static void _check_mh();
 static void _check_mh() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_dt;
  static double _sav_hbaserate;
  static double _sav_hbasetau;
  static double _sav_hexp;
  static double _sav_hgamma;
  static double _sav_htemp;
  static double _sav_hvalence;
  static double _sav_hvhalf;
  static double _sav_mbaserate;
  static double _sav_mbasetau;
  static double _sav_mexp;
  static double _sav_mgamma;
  static double _sav_mtemp;
  static double _sav_mvalence;
  static double _sav_mvhalf;
  static double _sav_celsius;
  static double _sav_mq10;
  static double _sav_hq10;
  static double _sav_vrest;
  static double _sav_vmin;
  static double _sav_vmax;
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_sav_hbaserate != hbaserate) { _maktable = 1;}
  if (_sav_hbasetau != hbasetau) { _maktable = 1;}
  if (_sav_hexp != hexp) { _maktable = 1;}
  if (_sav_hgamma != hgamma) { _maktable = 1;}
  if (_sav_htemp != htemp) { _maktable = 1;}
  if (_sav_hvalence != hvalence) { _maktable = 1;}
  if (_sav_hvhalf != hvhalf) { _maktable = 1;}
  if (_sav_mbaserate != mbaserate) { _maktable = 1;}
  if (_sav_mbasetau != mbasetau) { _maktable = 1;}
  if (_sav_mexp != mexp) { _maktable = 1;}
  if (_sav_mgamma != mgamma) { _maktable = 1;}
  if (_sav_mtemp != mtemp) { _maktable = 1;}
  if (_sav_mvalence != mvalence) { _maktable = 1;}
  if (_sav_mvhalf != mvhalf) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_sav_mq10 != mq10) { _maktable = 1;}
  if (_sav_hq10 != hq10) { _maktable = 1;}
  if (_sav_vrest != vrest) { _maktable = 1;}
  if (_sav_vmin != vmin) { _maktable = 1;}
  if (_sav_vmax != vmax) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_mh =  vmin ;
   _tmax =  vmax ;
   _dx = (_tmax - _tmin_mh)/200.; _mfac_mh = 1./_dx;
   for (_i=0, _x=_tmin_mh; _i < 201; _x += _dx, _i++) {
    _f_mh(_x);
    for (_j = 0; _j < 2; _j++) { _t_Add[_j][_i] = Add[_j];
}    for (_j = 0; _j < 2; _j++) { _t_Mult[_j][_i] = Mult[_j];
}   }
   _sav_dt = dt;
   _sav_hbaserate = hbaserate;
   _sav_hbasetau = hbasetau;
   _sav_hexp = hexp;
   _sav_hgamma = hgamma;
   _sav_htemp = htemp;
   _sav_hvalence = hvalence;
   _sav_hvhalf = hvhalf;
   _sav_mbaserate = mbaserate;
   _sav_mbasetau = mbasetau;
   _sav_mexp = mexp;
   _sav_mgamma = mgamma;
   _sav_mtemp = mtemp;
   _sav_mvalence = mvalence;
   _sav_mvhalf = mvhalf;
   _sav_celsius = celsius;
   _sav_mq10 = mq10;
   _sav_hq10 = hq10;
   _sav_vrest = vrest;
   _sav_vmin = vmin;
   _sav_vmax = vmax;
  }
 }

 static int mh(double _lv){ _check_mh();
 _n_mh(_lv);
 return 0;
 }

 static void _n_mh(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_mh(_lv); return; 
}
 _xi = _mfac_mh * (_lv - _tmin_mh);
 if (isnan(_xi)) {
  for (_j = 0; _j < 2; _j++) { Add[_j] = _xi;
}  for (_j = 0; _j < 2; _j++) { Mult[_j] = _xi;
}  return;
 }
 if (_xi <= 0.) {
 for (_j = 0; _j < 2; _j++) { Add[_j] = _t_Add[_j][0];
} for (_j = 0; _j < 2; _j++) { Mult[_j] = _t_Mult[_j][0];
} return; }
 if (_xi >= 200.) {
 for (_j = 0; _j < 2; _j++) { Add[_j] = _t_Add[_j][200];
} for (_j = 0; _j < 2; _j++) { Mult[_j] = _t_Mult[_j][200];
} return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 for (_j = 0; _j < 2; _j++) {double *_t = _t_Add[_j]; Add[_j] = _t[_i] + _theta*(_t[_i+1] - _t[_i]);}
 for (_j = 0; _j < 2; _j++) {double *_t = _t_Mult[_j]; Mult[_j] = _t[_i] + _theta*(_t[_i+1] - _t[_i]);}
 }

 
static int  _f_mh (  double _lv ) {
   double _la , _lb , _lj , _lmqq10 , _lhqq10 ;
 _lmqq10 = pow( mq10 , ( ( celsius - mtemp ) / 10. ) ) ;
   _lhqq10 = pow( hq10 , ( ( celsius - htemp ) / 10. ) ) ;
   {int  _lj ;for ( _lj = 0 ; _lj <= 1 ; _lj ++ ) {
     _la = alpha ( _threadargscomma_ _lv , ((double) _lj ) ) ;
     _lb = beta ( _threadargscomma_ _lv , ((double) _lj ) ) ;
     Inf [ _lj ] = _la / ( _la + _lb ) ;
     
/*VERBATIM*/
		switch (_lj) {
			case 0:
		/* Make sure Tau is not less than the base Tau */
				if ((Tau[_lj] = 1 / (_la + _lb)) < mbasetau) {
					Tau[_lj] = mbasetau;
				}
				Tau[_lj] = Tau[_lj] / _lmqq10;
				break;
			case 1:
				if ((Tau[_lj] = 1 / (_la + _lb)) < hbasetau) {
					Tau[_lj] = hbasetau;
				}
				Tau[_lj] = Tau[_lj] / _lhqq10;
				if (hexp==0) {
					Tau[_lj] = 1.; }
				break;
		}

 Mult [ _lj ] = exp ( - dt / Tau [ _lj ] ) ;
     Add [ _lj ] = Inf [ _lj ] * ( 1. - exp ( - dt / Tau [ _lj ] ) ) ;
     } }
    return 0; }
 
static void _hoc_mh(void) {
  double _r;
    _r = 1.;
 mh (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double alpha (  double _lv , double _lj ) {
   double _lalpha;
 if ( _lj  == 1.0 ) {
     if ( hexp  == 0.0 ) {
       _lalpha = 1.0 ;
       }
     else {
       _lalpha = hbaserate * exp ( ( _lv - ( hvhalf + vrest ) ) * hvalence * hgamma * FRT ( _threadargscomma_ htemp ) ) ;
       }
     }
   else {
     _lalpha = mbaserate * exp ( ( _lv - ( mvhalf + vrest ) ) * mvalence * mgamma * FRT ( _threadargscomma_ mtemp ) ) ;
     }
   
return _lalpha;
 }
 
static void _hoc_alpha(void) {
  double _r;
   _r =  alpha (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double beta (  double _lv , double _lj ) {
   double _lbeta;
 if ( _lj  == 1.0 ) {
     if ( hexp  == 0.0 ) {
       _lbeta = 1.0 ;
       }
     else {
       _lbeta = hbaserate * exp ( ( - _lv + ( hvhalf + vrest ) ) * hvalence * ( 1.0 - hgamma ) * FRT ( _threadargscomma_ htemp ) ) ;
       }
     }
   else {
     _lbeta = mbaserate * exp ( ( - _lv + ( mvhalf + vrest ) ) * mvalence * ( 1.0 - mgamma ) * FRT ( _threadargscomma_ mtemp ) ) ;
     }
   
return _lbeta;
 }
 
static void _hoc_beta(void) {
  double _r;
   _r =  beta (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double FRT (  double _ltemperature ) {
   double _lFRT;
 _lFRT = FARADAY * 0.001 / R / ( _ltemperature + 273.15 ) ;
   
return _lFRT;
 }
 
static void _hoc_FRT(void) {
  double _r;
   _r =  FRT (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double ghkca (  double _lv ) {
   double _lghkca;
 double _lnu , _lefun ;
 _lnu = _lv * 2.0 * FRT ( _threadargscomma_ celsius ) ;
   if ( fabs ( _lnu ) < 1.e-6 ) {
     _lefun = 1. - _lnu / 2. ;
     }
   else {
     _lefun = _lnu / ( exp ( _lnu ) - 1. ) ;
     }
   _lghkca = - FARADAY * 2.e-3 * _lefun * ( cao - cai * exp ( _lnu ) ) ;
   
return _lghkca;
 }
 
static void _hoc_ghkca(void) {
  double _r;
   _r =  ghkca (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int  iassign (  ) {
   i = g * ( v - erev ) ;
   ina = i ;
    return 0; }
 
static void _hoc_iassign(void) {
  double _r;
   _r = 1.;
 iassign (  );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("naf", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
 {
   mh ( _threadargscomma_ v ) ;
   if ( usetable  == 0.0 ) {
     m = Inf [ 0 ] ;
     h = Inf [ 1 ] ;
     }
   else {
     m = Add [ 0 ] / ( 1.0 - Mult [ 0 ] ) ;
     h = Add [ 1 ] / ( 1.0 - Mult [ 1 ] ) ;
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
   double _lhexp_val , _lindex , _lmexp_val ;
 _lhexp_val = 1.0 ;
   _lmexp_val = 1.0 ;
   if ( hexp > 0.0 ) {
     {int  _lindex ;for ( _lindex = 1 ; _lindex <= ((int) hexp ) ; _lindex ++ ) {
       _lhexp_val = h * _lhexp_val ;
       } }
     }
   if ( mexp > 0.0 ) {
     {int  _lindex ;for ( _lindex = 1 ; _lindex <= ((int) mexp ) ; _lindex ++ ) {
       _lmexp_val = m * _lmexp_val ;
       } }
     }
   g = gmax * _lmexp_val * _lhexp_val ;
   iassign ( _threadargs_ ) ;
   }
 _current += ina;

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
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
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
 { error =  states();
 if(error){fprintf(stderr,"at line 121 in file bg.inc:\n\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
  for (_i=0; _i < 2; _i++) {  _t_Add[_i] = makevector(201*sizeof(double)); }
  for (_i=0; _i < 2; _i++) {  _t_Mult[_i] = makevector(201*sizeof(double)); }
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/awile/workbench/modldb/nmodldb/models/db/modeldb/7399/mod/naf.mod";
static const char* nmodl_file_text = 
  "NEURON { SUFFIX naf  }\n"
  "  \n"
  "NEURON {  USEION na WRITE ina }\n"
  "\n"
  "ASSIGNED { ina }\n"
  "\n"
  "PARAMETER {\n"
  "	erev 		= 60       (mV)\n"
  "	gmax 		= 0.030     (umho)\n"
  "        vrest           = 0     : added in with bills change to borg\n"
  "\n"
  "	mvalence 	=  4.\n"
  "	mgamma 		=  0.5\n"
  "	mbaserate 	=  4.5\n"
  "	mvhalf 		=  -33.5\n"
  "	mbasetau 	=  0.02\n"
  "	mtemp 		=  37\n"
  "	mq10		=  3.\n"
  "	mexp 		=  3\n"
  "\n"
  "	hvalence 	= -6\n"
  "	hgamma		=  0.3\n"
  "	hbaserate 	=  0.095\n"
  "	hvhalf 		=  -39\n"
  "	hbasetau 	=  0.25\n"
  "	htemp 		=  37\n"
  "	hq10        =  3.\n"
  "	hexp 		=  1\n"
  "\n"
  "	cao                	 (mM)\n"
  "	cai                  (mM)\n"
  "\n"
  "	celsius			     (degC)\n"
  "	dt 				     (ms)\n"
  "	v 			         (mV)\n"
  "\n"
  "	vmax 		=  100     (mV)\n"
  "	vmin 		= -100   (mV)\n"
  "\n"
  "} : end PARAMETER\n"
  "\n"
  ":::INCLUDE \"bg.inc\"\n"
  ":::realpath /Users/awile/workbench/modldb/nmodldb/models/db/modeldb/7399/mod/bg.inc\n"
  ": $Id: bg.inc,v 1.3 1996/04/05 23:20:18 billl Exp $\n"
  "TITLE Borg-Graham Channel Model\n"
  "\n"
  "COMMENT\n"
  "\n"
  "Modeling the somatic electrical response of hippocampal pyramidal neurons, \n"
  "MS thesis, MIT, May 1987.\n"
  "\n"
  "Each channel has activation and inactivation particles as in the original\n"
  "Hodgkin Huxley formulation.  The activation particle mm and inactivation\n"
  "particle hh go from on to off states according to kinetic variables alpha\n"
  "and beta which are voltage dependent.  The form of the alpha and beta\n"
  "functions were dissimilar in the HH study.  The BG formulae are:\n"
  "	alpha = base_rate * Exp[(v - v_half)*valence*gamma*F/RT]\n"
  "	beta = base_rate * Exp[(-v + v_half)*valence*(1-gamma)*F/RT]\n"
  "where,\n"
  "	baserate : no affect on Inf.  Lowering this increases the maximum\n"
  "		    value of Tau\n"
  "	basetau : (in msec) minimum Tau value.\n"
  "	chanexp : number for exponentiating the state variable; e.g.\n"
  "		  original HH Na channel use m^3, note that chanexp = 0\n"
  "		  will turn off this state variable\n"
  "	erev : reversal potential for the channel\n"
  "	gamma : (between 0 and 1) does not affect the Inf but makes the\n"
  "		Tau more asymetric with increasing deviation from 0.5\n"
  "	celsius : temperature at which experiment was done (Tau will\n"
  "		      will be adjusted using a q10 of 3.0)\n"
  "	valence (z) : determines the steepness of the Inf sigmoid.  Higher\n"
  "		  valence gives steeper sigmoid.\n"
  "	vhalf : (a voltage) determines the voltage at which the value\n"
  "		 of the sigmoid function for Inf is 1/2\n"
  "	vmin, vmax : limits for construction of the table.  Generally,\n"
  "		     these should be set to the limits over which either\n"
  "		     of the 2 state variables are varying.\n"
  "        vrest : (a voltage) voltage shift for vhalf\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "	RANGE gmax, g, i\n"
  "	GLOBAL erev, Inf, Tau, Mult, Add, vmin, vmax, vrest\n"
  "} : end NEURON\n"
  "\n"
  "CONSTANT {\n"
  "	  FARADAY = 96489.0	: Faraday's constant\n"
  "	  R= 8.31441		: Gas constant\n"
  "\n"
  "} : end CONSTANT\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "} : end UNITS\n"
  "\n"
  "\n"
  "COMMENT\n"
  "** Parameter values should come from files specific to particular channels\n"
  "PARAMETER {\n"
  "	erev 		= 0    (mV)\n"
  "	gmax 		= 0    (mho/cm^2)\n"
  "        vrest           = 0    (mV)\n"
  "\n"
  "	mvalence 	= 0\n"
  "	mgamma 		= 0\n"
  "	mbaserate 	= 0\n"
  "	mvhalf 		= 0\n"
  "	mbasetau 	= 0\n"
  "	mtemp 		= 0\n"
  "	mq10		= 3\n"
  "	mexp 		= 0\n"
  "\n"
  "	hvalence 	= 0\n"
  "	hgamma		= 0\n"
  "	hbaserate 	= 0\n"
  "	hvhalf 		= 0\n"
  "	hbasetau 	= 0\n"
  "	htemp 		= 0\n"
  "	hq10		= 3\n"
  "	hexp 		= 0\n"
  "\n"
  "	cao                (mM)\n"
  "	cai                (mM)\n"
  "\n"
  "	celsius			   (degC)\n"
  "	dt 				   (ms)\n"
  "	v 			       (mV)\n"
  "\n"
  "	vmax 		= 100  (mV)\n"
  "	vmin 		= -100 (mV)\n"
  "} : end PARAMETER\n"
  "ENDCOMMENT\n"
  "\n"
  "ASSIGNED {\n"
  "	i (mA/cm^2)		\n"
  "	g (mho/cm^2)\n"
  "	Inf[2]		: 0 = m and 1 = h\n"
  "	Tau[2]		: 0 = m and 1 = h\n"
  "	Mult[2]		: 0 = m and 1 = h\n"
  "	Add[2]		: 0 = m and 1 = h\n"
  "} : end ASSIGNED \n"
  "\n"
  "STATE { m h }\n"
  "\n"
  "INITIAL { \n"
  " 	mh(v)\n"
  "	if (usetable==0) {\n"
  " 	  m = Inf[0] h = Inf[1]\n"
  "	} else {\n"
  " 	  m = Add[0]/(1-Mult[0]) h = Add[1]/(1-Mult[1]) \n"
  "	}\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "\n"
  "	LOCAL hexp_val, index, mexp_val\n"
  "\n"
  "	SOLVE states\n"
  "\n"
  "	hexp_val = 1\n"
  "	mexp_val = 1\n"
  "\n"
  "	: Determining h's exponent value\n"
  "	if (hexp > 0) {\n"
  "		FROM index=1 TO hexp {\n"
  "			hexp_val = h * hexp_val\n"
  "		}\n"
  "	}\n"
  "\n"
  "	: Determining m's exponent value\n"
  "	if (mexp > 0) {\n"
  "		FROM index = 1 TO mexp {\n"
  "			mexp_val = m * mexp_val\n"
  "		}\n"
  "	}\n"
  "\n"
  "	:			       mexp			    hexp\n"
  "	: Note that mexp_val is now = m      and hexp_val is now = h \n"
  "	g = gmax * mexp_val * hexp_val\n"
  "	iassign()\n"
  "} : end BREAKPOINT\n"
  "\n"
  ": ASSIGNMENT PROCEDURES\n"
  ": Must be overwritten by user routines in parameters.multi\n"
  ": PROCEDURE iassign () { i = g*(v-erev) ina=i }\n"
  ": PROCEDURE iassign () { i = g*ghkca(v) ica=i }\n"
  "\n"
  ":-------------------------------------------------------------------\n"
  ": I suppose we have 2 choices, to use the DERIVATIVE function or\n"
  ": to explicitly state m+ and h+.  If you were to use the DERIVATIVE\n"
  ": function, then you will do as follows:\n"
  ": DERIVATIVE deriv {\n"
  ":	m' = (-m + minf) / mtau\n"
  ":	h' = (-h + hinf) / htau\n"
  ": }\n"
  ": Else, since m' = (m+ - m) / dt, setting the 2 equations together,\n"
  ": we can solve for m+ and eventually get :\n"
  ": 	m+ = (m * mtau + dt * minf) / (mtau + dt)\n"
  ": and same for h+:\n"
  ": 	h+ = (h * htau + dt * hinf) / (htau + dt)\n"
  ": This is the one we will use, so ...\n"
  "PROCEDURE states() {\n"
  "\n"
  "	: Setup the mh table values\n"
  "\n"
  "	mh (v*1(/mV))\n"
  "	m = m * Mult[0] + Add[0]\n"
  "	h = h * Mult[1] + Add[1]\n"
  "\n"
  "	VERBATIM\n"
  "	return 0;\n"
  "	ENDVERBATIM	\n"
  "}\n"
  "\n"
  ":-------------------------------------------------------------------\n"
  ": NOTE : 0 = m and 1 = h\n"
  "PROCEDURE mh (v) {\n"
  "	LOCAL a, b, j, mqq10, hqq10\n"
  "	TABLE Add, Mult DEPEND dt, hbaserate, hbasetau, hexp, hgamma, htemp, hvalence, hvhalf, mbaserate, mbasetau, mexp, mgamma, mtemp, mvalence, mvhalf, celsius, mq10, hq10, vrest, vmin, vmax  FROM vmin TO vmax WITH 200\n"
  "\n"
  "	mqq10 = mq10^((celsius-mtemp)/10.)	\n"
  "	hqq10 = hq10^((celsius-htemp)/10.)	\n"
  "\n"
  "	: Calculater Inf and Tau values for h and m\n"
  "	FROM j = 0 TO 1 {\n"
  "		a = alpha (v, j)\n"
  "		b = beta (v, j)\n"
  "\n"
  "		Inf[j] = a / (a + b)\n"
  "\n"
  "		VERBATIM\n"
  "		switch (_lj) {\n"
  "			case 0:\n"
  "		/* Make sure Tau is not less than the base Tau */\n"
  "				if ((Tau[_lj] = 1 / (_la + _lb)) < mbasetau) {\n"
  "					Tau[_lj] = mbasetau;\n"
  "				}\n"
  "				Tau[_lj] = Tau[_lj] / _lmqq10;\n"
  "				break;\n"
  "			case 1:\n"
  "				if ((Tau[_lj] = 1 / (_la + _lb)) < hbasetau) {\n"
  "					Tau[_lj] = hbasetau;\n"
  "				}\n"
  "				Tau[_lj] = Tau[_lj] / _lhqq10;\n"
  "				if (hexp==0) {\n"
  "					Tau[_lj] = 1.; }\n"
  "				break;\n"
  "		}\n"
  "\n"
  "		ENDVERBATIM\n"
  "		Mult[j] = exp(-dt/Tau[j])\n"
  "		Add[j]  = Inf[j]*(1. - exp(-dt/Tau[j]))\n"
  "	}\n"
  "} : end PROCEDURE mh (v)\n"
  "\n"
  ":-------------------------------------------------------------------\n"
  "FUNCTION alpha(v,j) {\n"
  "	if (j == 1) {\n"
  "	   if (hexp==0) {\n"
  "	     alpha = 1\n"
  "	   } else {\n"
  "             alpha = hbaserate * exp((v - (hvhalf+vrest)) * hvalence * hgamma * FRT(htemp)) }\n"
  "	} else {\n"
  "          alpha = mbaserate * exp((v - (mvhalf+vrest)) * mvalence * mgamma * FRT(mtemp))\n"
  "	}\n"
  "} : end FUNCTION alpha (v,j)\n"
  "\n"
  ":-------------------------------------------------------------------\n"
  "FUNCTION beta (v,j) {\n"
  "	if (j == 1) {\n"
  "	   if (hexp==0) {\n"
  "                beta = 1\n"
  "	   } else {\n"
  "		beta = hbaserate * exp((-v + (hvhalf+vrest)) * hvalence * (1 - hgamma) * FRT(htemp)) }\n"
  "	} else {\n"
  "		beta = mbaserate * exp((-v + (mvhalf+vrest)) * mvalence * (1 - mgamma) * FRT(mtemp))\n"
  "	}\n"
  "} : end FUNCTION beta (v,j)\n"
  "\n"
  ":-------------------------------------------------------------------\n"
  "FUNCTION FRT(temperature) {\n"
  "	FRT = FARADAY * 0.001 / R / (temperature + 273.15)\n"
  "} : end FUNCTION FRT (temperature)\n"
  "\n"
  ":-------------------------------------------------------------------\n"
  " FUNCTION ghkca (v) { : Goldman-Hodgkin-Katz eqn\n"
  "       LOCAL nu, efun\n"
  "\n"
  "       nu = v*2*FRT(celsius)\n"
  "       if(fabs(nu) < 1.e-6) {\n"
  "               efun = 1.- nu/2.\n"
  "       } else {\n"
  "               efun = nu/(exp(nu)-1.) }\n"
  "       ghkca = -FARADAY*2.e-3*efun*(cao - cai*exp(nu))\n"
  " } : end FUNCTION ghkca()\n"
  ":::end INCLUDE bg.inc\n"
  "\n"
  "PROCEDURE iassign () { i = g*(v-erev) ina=i }\n"
  ":**kdr\n"
  ":________________________________________________________________\n"
  ": Borg-Graham version of delayed rectifier\n"
  "\n"
  ;
#endif
