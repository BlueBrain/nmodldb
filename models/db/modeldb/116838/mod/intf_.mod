: $Id: intf.mod,v 1.789 2008/12/20 16:04:59 billl Exp $

:* main COMMENT
COMMENT
artificial cell incorporating 4 input weights with different time constants and signs
typically a fast AMPA, slow NMDA, fast GABAA, slow GABAB
features:
  1. Mg dependence for NMDA activation
  2. "G-protein" cooperativity for GABAB activation
  3. depolarization blockade
  4. AHP affects both Vm and refractory period  (adaptation)
  5. decrementing excitatory and/or inhibitory activity post spk (another adaptation)
since artificial cells only do calculations when they receive events, a set of vec
  pointers are maintained to allow state var information storage when event arrives
  (see initrec() and record())
ENDCOMMENT

:* main VERBATIM block
VERBATIM
#if !defined(t)
#define _pval pval
#endif
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h> /* contains LONG_MAX */

extern void* vector_arg(int);
extern double hoc_call_func(Symbol*, int);
extern double* hoc_pgetarg(int);
extern FILE* hoc_obj_file_arg(int);
extern int list_vector_px(Object*, int, double**);
extern int list_vector_px2 (Object*, int, double**, void**);
extern int list_vector_resize(Object*, int, int);//* list_vector_resize(LIST,ITEM#,NEW SIZE)
extern Object** hoc_objgetarg(int);
extern int ivoc_list_count(Object*);
extern Object* ivoc_list_item(Object*, int);
extern Symbol *hoc_get_symbol(char *);
extern Symbol *hoc_lookup(const char*);
extern Point_process* ob2pntproc(Object*);
extern double mcell_ran4(unsigned int*,double*,int,double);
extern int hoc_is_double_arg(int narg);
static void hxe() { hoc_execerror("",0); }
static int ctt(unsigned int, char**);
static int setdvi2(double*,double*,int,int);
extern int stoprun;
extern double hoc_epsilon;
extern short *nrn_artcell_qindex_;
extern double nrn_event_queue_stats(double*);
extern void clear_event_queue();
extern double *vector_newsize(void*, int);
extern Objectdata *hoc_objectdata;
extern int nrn_mlh_gsort(double*,int*,int,int());
extern int cmpdfn(double,double);
extern unsigned int scrsz;
extern unsigned int *scr;
extern unsigned int *scrset(int);

typedef struct{
  int isz;
  double** pv;
  int* plen;
} ListVec;

extern void FreeListVec(ListVec** pp);
extern ListVec* AllocListVec(Object* p);

#define PI 3.141592653589793115997963468544
#define nil 0
#define CTYPp 20 // CTYPp>CTYPi from labels.hoc
#define SOP (((id0*) _p_sop)->vp)
#define IDP (*((id0**) &(_p_sop)))
#define NSW 20  // just store voltages
#define NSV 8  // 7 state variables (+ 1 for time)
#define FOFFSET 100 // flag offset for net_receive()
#define WRNUM 5  // a single INTF can store into this many ww field vecs
#define DELM(X,Y) (*(pg->delm+(X)*CTYPi+(Y)))
#define DELD(X,Y) (*(pg->deld+(X)*CTYPi+(Y)))
#define DVG(X,Y) ((int)*(pg->dvg+(X)*CTYPi+(Y)))
#define WMAT(X,Y,Z) (*(pg->wmat+(X)*CTYPi*STYPi+(Y)*STYPi+(Z)))
#define WD0(X,Y,Z)  (*(pg->wd0 +(X)*CTYPi*STYPi+(Y)*STYPi+(Z)))
#define NUMC(X) (*(pg->numc+(X)))
#define HVAL(X) (*(hoc_objectdata[(hoc_get_symbol((X)))->u.oboff]._pval))
#define HPTR(X) (hoc_objectdata[(hoc_get_symbol((X)))->u.oboff]._pval)

typedef struct VPT {
 unsigned int  id;
 unsigned int  size;
 unsigned int  p;
 void*    vv[NSV];
 double* vvo[NSV];
} vpt;

typedef struct POSTGRP { // postsynaptic group
  double *dvg; double *delm; double *deld; double *ix; double *ixe; double *wmat; double *wd0;
  double *numc;
  struct POSTGRP *next;
} postgrp;

typedef struct ID0 {
  vpt     *vp;
  postgrp *pg;
  float    wscale[WRNUM];
  Point_process **dvi; // each cell has a divergence list
  double *del;        // each syn has its own intrinsic delay
  unsigned char *sprob;    // each syn has a firing probability 0-255->0-1
  unsigned int dvt;
  unsigned int  id;
  unsigned int col;
  unsigned int rvb;
  unsigned int rvi;
  unsigned int spkcnt;
  unsigned int blkcnt;
  int rve;
  char   wreci[WRNUM]; // since use -1 as a flag
  char   errflag;
  // type -> vbr MUST REMAIN unbroked BLOCK -- see flag()
  // when adding flags also augment iflags, iflnum
  // only use first 3 letters with flag() -- see iflags
  unsigned char     type;  // | 
  unsigned char     inhib; // | 
  unsigned char     record;// |
  unsigned char     wrec;  // |
  unsigned char     jttr;  // |
  unsigned char     input; // |
  unsigned char     vinflg;// |
  unsigned char     invl0; // |
  unsigned char     jcn;   // |
  unsigned char     dead;  // |
  unsigned char     vbr;   // |
           char     dbx;   // |
           char     flag;  // |
  // end BLOCK
} id0;

// globals -- range vars must be malloc'ed in the CONSTRUCTOR
static vpt *vp; // vp, pg, ip are used as temporary pointers
static id0 *ip, *qp, *rp;
static postgrp *pg;
static Object *ce, *CTYP;
static Point_process *pmt, *tpnt;
static char *name;
// iflags string use to find flags -- note that only 1st 3 chars are used to identify
static char iflags[100]="typ inh rec wre jtt inp vin inv jcn dea vbr dbx fla"; 
static char iflnum=13, iflneg=11, errflag;      // turn on after generating an error message
static int  vspn;
static double *isp, *vsp, *wsp, *jsp, *invlp;
static double *lop(Object *ob, unsigned int i); // accessed by all INTF
static double *jrid, *jrtv;
static void   *jridv, *jrtvv;
static unsigned int jtpt,jtmax,jrmax;
static unsigned long int jri;
static unsigned long int spktot, eventtot;
static double vii[NSV];   // temp storage
static unsigned int wwpt,wwsz,wwaz,cesz; // pointer, size for shared ww vectors
static unsigned int sead, spikes[CTYPp], blockcnt[CTYPp]; // 'sead' vs global 'seed'/ used elsewhere
static unsigned int AMo[CTYPp],NMo[CTYPp],GAo[CTYPp],GBo[CTYPp]; // count overages for types
static char* CNAME[CTYPp]; // 20 should be > CTYPi
static int cty[CTYPp], process;
static int CTYN, CTYPi, STYPi, dscrsz; // from labels.hoc
static double qlimit, *dscr;
FILE *wf1, *wf2, *tf;
void*    ww[NSW];
double* wwo[NSW];
static int AM=0, NM=1, GA=2, GB=3, SU=3, IN=4, DP=2; // from labels.hoc
static double wts[10],hsh[10];  // for jitcons to use as a junk pointer
ENDVERBATIM

:* NEURON, PARAMETER, ASSIGNED blocks
NEURON {
  ARTIFICIAL_CELL INTF
  RANGE VAM, VNM, VGA, VGB, AHP           :::: cell state variables
  RANGE Vm                                :::: derived var
  : parameters
  RANGE tauAM, tauNM, tauGA, tauGB     :::: synaptic params
  RANGE tauahp, ahpwt                  :::: intrinsic params
  RANGE tauRR , RRWght                 :::: relative refrac. period tau, wght of Vblock-VTH for refrac
  RANGE RMP,VTH,Vblock,VTHC,VTHR       :::: Vblock for depol blockade
  RANGE nbur,tbur,refrac               :::: burst size, interval; refrac period and extender
  RANGE invl,oinvl,WINV,invlt           :::: interval bursting params
  RANGE VGBdel,tGB,VGBa,Vbrefrac,rebound,rebob,offsetGB   :::: GABAB and rebound
  RANGE STDAM, STDNM, STDGA             :::: specific amounts of STD for each type of synapse
  GLOBAL EAM, ENM, EGA, EGB, mg         :::: "reverse potential" distance from rest
  GLOBAL tauGBGP,wGBGP,GPkd,Gn          :::: GABAB G-protein dependence (cooperativity)
  GLOBAL spkht, wwwid,wwht              :::: display: spike height, width/ht for pop spikes
  GLOBAL stopoq                         :::: flags: stop if q is empty, use STD
  : other stuff
  POINTER sop                          :::: Structure pointer for other range vars
  RANGE  xloc,yloc,zloc
  RANGE  t0,tg,twg,tGB,refractory,trrs :::: t0,tg,tGB save times for analytic calc
  RANGE  cbur,WEX                      :::: burst statevar
  GLOBAL vdt,nxt,RES,ESIN,Bb,Psk      :::: table look up values for exp,sin
  GLOBAL prnum, nsw, rebeg             :::: for debugging moves
  GLOBAL subsvint, jrsvn, jrsvd, jrtime, jrtm :::: output params
  GLOBAL DEAD_DIV, seedstep            :::: dead cells on div list?
  GLOBAL seaddvioff                    :::: seed offset for dvi/del 
  GLOBAL WVAR,DELMIN                   :::: weight variance must be > 0 <= 1.0
  GLOBAL tmax,installed,verbose        :::: simplest output
  GLOBAL pathbeg,pathend,PATHMEASURE,pathidtarg,pathtytarg,seadsetting,pathlen
}

PARAMETER {
  tauAM = 10 (ms)
  tauNM = 300 (ms)
  tauGA = 10 (ms)
  tauGB = 300 (ms)
  tauGBGP = 50 (ms) : drop off for burst effect on GABAB (G-protein)
  invl =  100 (ms)
  WINV =  0
  wGBGP = 1 (ms) : augmentation of G-protein with a spike
  GPkd  = 100    : maintain between 50 and 500 since table only up to 10 spikes (with wGBGP=1)
  ahpwt = 0
  tauahp= 10 (ms)
  tauRR = 6 (ms)
  refrac = 5 (ms)
  Vbrefrac = 20 (ms)
  RRWght = 0.75
  wwwid = 10
  wwht = 10
  VTH = -45      : fixed spike threshold
  VTHC = -45
  VTHR = -45
  Vblock = -20   : level of depolarization blockade
  vdt = 0.1      : time step for saving state var
  mg = 1         : for NMDA Mg dep.
  sop=0
  nbur=1
  tbur=2
  VGBdel=0
  rebound=0.01 : the is param for one GB mech
  offsetGB=0
  RMP=-65
  EAM = 65
  ENM = 90
  EGA = -15
  EGB = -30
  spkht = 50
  prnum = -1
  nsw=0
  rebeg=0
  subsvint=0
  jrsvn=1e4 jrsvd=1e4 jrtime=-1 jrtm=-1
  seedstep=44340
  seaddvioff=9102098713763e-134
  DEAD_DIV=1
  WVAR=0.2
  stopoq=0
  PATHMEASURE=0
  verbose=0
  seadsetting=0
  pathidtarg=-1
  DELMIN=1e-5 : min delay to bother using queue -- otherwise considered simultaneous
  STDAM=0
  STDNM=0
  STDGA=0
}

ASSIGNED {
  Vm VAM VNM VGA VGB AHP
  VGBa t0 tGB tg twg refractory nxt xloc yloc zloc trrs
  WEX RES ESIN Gn Bb Psk cbur invlt oinvl rebob tmax
  installed
  pathbeg pathend pathtytarg pathlen
}

:* CONSTRUCTOR, DESTRUCTOR, INITIAL
:** CONSTRUCT: create a structure to save the identity of this unit and char integer flags
CONSTRUCTOR {
  VERBATIM 
  { int lid,lty,lin,lco,i; unsigned int sz;
    if (ifarg(1)) { lid=(int) *getarg(1); } else { lid= UINT_MAX; }
    if (ifarg(2)) { lty=(int) *getarg(2); } else { lty= -1; }
    if (ifarg(3)) { lin=(int) *getarg(3); } else { lin= -1; }
    if (ifarg(4)) { lco=(int) *getarg(4); } else { lco= -1; }
    _p_sop = (void*)ecalloc(1, sizeof(id0)); // important that calloc sets all flags etc to 0
    ip = IDP;
    ip->id=lid; ip->type=lty; ip->inhib=lin; ip->col=lco; 
    ip->pg=0x0; ip->dvi=0x0; ip->sprob=0x0;
    ip->dead = ip->invl0 = ip->record = ip->jttr = ip->input = 0; // all flags off
    ip->dvt = ip->vbr = ip->wrec = ip->jcn = 0;
    for (i=0;i<WRNUM;i++) {ip->wreci[i]=-1; ip->wscale[i]=-1.0;}
    ip->rve=-1;
    pathbeg=-1;
    jrmax=0;
    process=(int)getpid();
    CNAME[SU]="SU"; CNAME[DP]="DP"; CNAME[IN]="IN";
    if (installed==2.0) { // jitcondiv was previously run
      sz=ivoc_list_count(ce);
      printf("\t**** WARNING new INTF created: may want to rerun jitcondiv ****\n");
    } else installed=1.0; // set or reset it
  }
  ENDVERBATIM
}

DESTRUCTOR {
  VERBATIM { 
  free(IDP);
  }
  ENDVERBATIM
}

:** INITIAL
INITIAL { LOCAL id
  reset() 
  tg = 0
  twg = 0
  tmax=0
  pathend=-1
  pathlen=0
  VERBATIM
  { int i,ix;
  ip=IDP;
  _lid=(double)ip->id;
  ip->spkcnt=0;
  ip->blkcnt=0;
  ip->errflag=0;
  for (i=0;i<CTYN;i++){ix=cty[i]; blockcnt[ix]=spikes[ix]=AMo[ix]=NMo[ix]=GAo[ix]=GBo[ix]=0;}
  }
  ENDVERBATIM
  jrsvn=jrsvd jrtime=jrtm
  : init with vinset(0) if will turn on via a NetCon with w5=1
  if (vinflag()) { randspk() net_send(nxt,2)}
  if (recflag()) { recini() } : recini() resets for recording, cf recinit()
  if (pathbeg==id) { 
    stoprun=0 
    net_send(0,2) 
  } : send at time 0
  rebeg=0 : will reset this to restart storage for rec,wrec
}

PROCEDURE reset () {
  Vm = RMP
  VAM = 0
  VNM = 0
  VGA = 0
  VGB = 0
  VGBa = 0
  offsetGB=0
  AHP=0
  rebob=-1e9
  invlt = -1
  t0 = t
  tGB = t
  refractory = 0 : 1 means cell is absolute refractory
  VTHC=VTH :set current threshold to absolute threshold value
}

VERBATIM
unsigned int GetDVIDSeedVal(unsigned int id) {
  double x[2];
  if (seadsetting==1) { 
    sead=((unsigned int)ip->id+seaddvioff)*1e6;
  } else { 
    if (seadsetting==2) printf("Warning: GetDVIDSeedVal called with wt rand turned off\n");
    x[0]=(double)id; x[1]=seaddvioff;
    sead=hashseed2(2,&x);
  }
  return sead;
}
ENDVERBATIM

: seed for divergence and delays -- not yet used
FUNCTION DVIDSeed(){
  VERBATIM
  return (double)GetDVIDSeedVal(IDP->id);
  ENDVERBATIM
}

:* NET_RECEIVE
NET_RECEIVE (wAM,wNM,wGA,wGB,wflg) { LOCAL tmp,jcn,id
  INITIAL { wNM=wNM wGA=wGA wGB=wGB wflg=0}
  : intra-burst, generate next spike as needed
VERBATIM
  int prty,poty,prin,ii,sy; double STDf; //@

ENDVERBATIM
  tmax=t
  VERBATIM
  if (stopoq && !qsz()) stoprun=1;
  ip=IDP; pg=ip->pg;
  if (ip->dead) return; // this cell has died
  _ljcn=ip->jcn; _lid=ip->id;
  tpnt = _pnt; // this pnt
  if (PATHMEASURE) { // do all code for this
    if (_lflag==2 || _lflag<0) { // on the callback -- distribute to divergence list
      double idty; int i;
      if (_lflag==2) ip->flag=-1; 
      idty=(double)(FOFFSET+ip->id)+1e-2*(double)ip->type+1e-3*(double)ip->inhib+1e-4;
      for (i=0;i<ip->dvt && !stoprun;i++) if (ip->sprob[i]) {
        (*pnt_receive[ip->dvi[i]->_prop->_type])(ip->dvi[i], wts, idty); 
        _p=_pnt->_prop->param; _ppvar=_pnt->_prop->dparam; ip=IDP; // restore pointers each time
      }
      return;  // else see if destination has been reached
    } else if (_lflag!=2 && (pathtytarg==(double)ip->type || pathidtarg==(double)ip->id)) {
      if (pathend==(double)ip->id) return; // means that coming back here again
      ip->flag=(unsigned char)floor(t)+1; // type-target or id-target
      pathend=(double)ip->id; 
      pathlen=tmax+1; // tmax gives pathlength
      stoprun=1.; 
      return;
      // deadends:visited || no output  ||stopped
    } else if (ip->flag   || ip->dvt==0 || stoprun) {
      return; // inhib cell is a deadend; don't revisit anyone
    } else if (ip->inhib) {
      if (!ip->flag) ip->flag=(unsigned char)floor(t)+1;
    } else { // first callback will be from the stim
      ip->flag=(unsigned char)floor(t)+1;
#if defined(t)
      net_send((void**)0x0, wts,tpnt,t+1.,-1.); // the callback call
#else
      net_send((void**)0x0, wts,tpnt,1.,-1.); // the callback call
#endif
      return;
    }
  }
  if (_lflag<0) { callback(_lflag); return; }
  eventtot+=1;
  ENDVERBATIM
VERBATIM
  if (ip->dbx>2) 
ENDVERBATIM
{ 
    pid() 
    printf("DB0: flag=%g Vm=%g",flag,VAM+VNM+VGA+VGB+RMP+AHP)
    if (flag==0) { printf(" (%g %g %g %g %g)",wAM,wNM,wGA,wGB,wflg) }
    printf("\n")
  }
: causes of spiking: between VTH and Vblock, random from vsp (flag 2), within burst
:** JITcon code
  if (flag>=FOFFSET) { : jitcon -- set up weights on the fly
    VERBATIM {
      // find type of presyn
      poty=(int)ip->type;
      prty=(int)(1e2*(_lflag-floor(_lflag)));
      prin=(int)(1e3*(_lflag-floor(_lflag)-prty*1e-2)); // stuffed into this flag
      sy=prin?GA:AM;
      STDf=_args[0]; // save value
      for (ii=0;ii<=4;ii++) _args[ii]=0.; // may be unnecessary to clear
      if (seadsetting==2) { // fixed weights
        for (ii=sy;ii<sy+2;ii++) _args[ii]=WMAT(prty,poty,ii)*WD0(prty,poty,ii);
      } else { 
        if (seadsetting==1) {
          sead=(unsigned int)(floor(_lflag)*ip->id*seedstep); // all integers
        } else { // hash on presynaptic id+FOFFSET,poid,seedstep
          hsh[0]=floor(_lflag); hsh[1]=(double)ip->id; hsh[2]=seedstep;
          sead=hashseed2(3,&hsh); // hsh[] is just scratch pad
        }
        mcell_ran4(&sead, &_args[sy], 2, 1.);
        for (ii=sy;ii<sy+2;ii++) { // scale appropriately; 
          _args[ii]=2*WVAR*(_args[ii]+0.5/WVAR-0.5)*WMAT(prty,poty,ii)*WD0(prty,poty,ii);
        }
      }
    }
    ENDVERBATIM
VERBATIM
    if (ip->dbx>2) 
ENDVERBATIM
{ 
      pid() 
      printf("DF: flag=%g Vm=%g",flag,VAM+VNM+VGA+VGB+RMP+AHP)
      printf(" (%g %g %g %g %g)",wAM,wNM,wGA,wGB,wflg)
      printf("\n")
    }
:** mid-burst
  } else if (flag==4) { 
    cbur=cbur-1  : count down the spikes
    if (cbur>0) { 
      net_send(tbur,4) 
    } else { 
      net_send(refrac, 3)
    }
    tmp=t
VERBATIM
    if (ip->jttr) 
ENDVERBATIM
{ tmp= t+jttr()/10 } 
    if (jcn) { jitcon(tmp) } else { net_event(tmp) }
VERBATIM
    spikes[ip->type]++; //@

ENDVERBATIM
VERBATIM
    if (ip->dbx>0) 
ENDVERBATIM
{ pid() printf("DBA: mid-burst event at %g, %g\n",tmp,cbur) } 
VERBATIM
    if (ip->record) 
ENDVERBATIM
{ recspk(tmp) } 
VERBATIM
    if (ip->wrec) 
ENDVERBATIM
{ wrecord(t) } 
VERBATIM
    return; //@ done

ENDVERBATIM
    : start reading random spike times (or burst times) from vsp vector pointer
    : this is signaled externally from a netstim with wflg=1, will turn off on next stim 
    : (NB wflg used in completely different context for GABAB) 
    : this is bad -- should use a special netcon that just handles signals
  } else if (flag==0 && wGB==0 && wflg==1) {
VERBATIM
    ip->input=1; //@

ENDVERBATIM
    wflg=2 : set flag to turn off next time an external event comes from here
    randspk() 
    net_send(nxt,2)
VERBATIM
    return; //@ done

ENDVERBATIM
  } else if (flag==0 && wGB==0 && wflg==2) { : flag to stop random spikes
VERBATIM
    ip->input=0; //@ inputs that are read from a vector of times -- see randspk()

ENDVERBATIM
    wflg=1  : flag to turn on next time
VERBATIM
    return; //@ done

ENDVERBATIM
  }
  : update state variables
VERBATIM
  if (ip->record) 
ENDVERBATIM
{ record() } 
VERBATIM
  if (ip->wrec) 
ENDVERBATIM
{ wrecord(1e9) } 
:** update state variables: VAM, VNM, VGA, VGB
  if (VAM>hoc_epsilon)  { VAM = VAM*EXP(-(t - t0)/tauAM) } else { VAM=0 } :AMPA
  if (VNM>hoc_epsilon)  { VNM = VNM*EXP(-(t - t0)/tauNM) } else { VNM=0 } :NMDA
  if (VGA< -hoc_epsilon){ VGA = VGA*EXP(-(t - t0)/tauGA) } else { VGA=0 } :GABAA    
  if(refractory==0){:once refractory period over, VTHC falls back towards VTH
    if(VTHC>VTH) { VTHC = VTH + (VTHR-VTH)*EXP(-(t-trrs)/tauRR) }
  }
  if (VGBdel>0) {
    VGB = esinr(t-tGB) :VGB has to update each t but calc based on triggering tGB and val VGBa
  } else {
    if (VGB< -hoc_epsilon){ 
      VGB = VGB*EXP(-(t - t0)/tauGB) } else { VGB=0 }
  }      
  if (AHP< -hoc_epsilon){ AHP = AHP*EXP(-(t-t0)/tauahp) } else { AHP=0 } : adaptation
  t0 = t : finished using t0
  Vm = VAM+VNM+VGA+VGB+AHP : membrane deviation from rest
  if (Vm> -RMP) {Vm= -RMP}: 65 mV above rest
  if (Vm<  RMP) {Vm= RMP} : 65 mV below rest
:*** only add weights if an external excitation
  if (flag==0 || flag>=FOFFSET) { 
    : AMPA Erev=0 (0-RMP==65 mV above rest)
    if (wAM>0) {
      if (rebob!=1e9 && rebob!=-1e9) {
VERBATIM
        cbur=floor(rebound*rebob/EGB); //@

ENDVERBATIM
        net_send(tbur,4) 
        rebob=1e9
      }
      if (STDAM==0) { VAM = VAM + wAM*(1-Vm/EAM)
      } else        { VAM = VAM + (1-STDAM*STDf)*wAM*(1-Vm/EAM) }
      if (VAM>EAM) { 
VERBATIM
        AMo[ip->type]++; //@

ENDVERBATIM
      } else if (VAM<0) { VAM=0 }
    }
    : NMDA; Mg effect based on total activation in rates()
    if (wNM>0 && VNM<ENM) { 
      rates(RMP+Vm)
      if (STDNM==0) { VNM = VNM + wNM*Bb*(1-Vm/ENM) 
      } else        { VNM = VNM + (1-STDNM*STDf)*wNM*Bb*(1-Vm/ENM) }
      if (VNM>ENM) { 
VERBATIM
        NMo[ip->type]++; //@

ENDVERBATIM
      } else if (VNM<0) { VNM=0 }
    }
    : GABAA and GABAB: note that all wts are positive
    if (wGA>0 && VGA>EGA) { : the neg here gives the inhibition
      if (STDGA==0) {  VGA = VGA - wGA*(1-Vm/EGA) 
      } else {         VGA = VGA - (1-STDGA*STDf)*wGA*(1-Vm/EGA) }
      if (VGA<EGA) { 
VERBATIM
        GAo[ip->type]++; //@

ENDVERBATIM
VERBATIM
  if (ip->dbx>2) 
ENDVERBATIM
{ 
    pid() printf("DB0A: flag=%g Vm=%g",flag,VAM+VNM+VGA+VGB+RMP+AHP)
    if (flag==0) { printf(" (%g %g %g %g %g %g)",wGA,EGA,VGA,Vm,AHP,STDf) }  
VERBATIM
    printf("\nAA:%d:%d\n",GAo[ip->type],ip->type); //@ 

ENDVERBATIM
    printf("\n")
  }
      } else if (VGA>0) { VGA=0 } : if want reversal of VGA need to also edit above
    }
    if (wGB>1e-6) {
      if (VGBdel>0) { net_send(VGBdel,5)  : delayed effect
      } else { : handle wGB immediately
        : wflg will augment each time there is spike coming in through this line
        wflg=wflg*EXP(-(t-tGB)/tauGBGP)+wGBGP 
        coop(wflg)               : cooperativity -- need mult presyn spikes to activate
        VGB = VGB - wGB*(1-Vm/EGB)*Gn
        if (VGB<EGB) { 
VERBATIM
          GBo[ip->type]++; //@

ENDVERBATIM
        }
        if (VGB<rebob && rebob!=1e9 && rebob!=-1e9) { rebob=VGB }
      }
    }
:*** modulated interval firing; cf invlfire.mod
VERBATIM
    if (ip->invl0) 
ENDVERBATIM
{ 
      Vm = RMP+VAM+VNM+VGA+VGB+AHP
      if (Vm>0)   {Vm= 0 }
      if (Vm<-90) {Vm=-90}
      if (invlt==-1) { : activate for first time
        if (Vm>RMP) {
          oinvl=invl
          invlt=t
          net_send(invl,1) 
        }
      } else {
        tmp=shift(Vm)
        if (tmp!=0)  {
          net_move(tmp) 
          if (id()<prnum) {
            pid() printf("**** MOVE t=%g to %g Vm=%g %g,%g\n",t,tmp,Vm,invlt,oinvl) }
        }
      }      
    }
  } else if (flag==1) { : modulated interval firing; cf invlfire.mod
    : Vm=RMP+VAM+VNM+VGA+VGB+AHP
    if (WINV<0) { 
      if (jcn) { jitcon(t) } else { net_event(t) } : bypass activation calculation
VERBATIM
      spikes[ip->type]++; //@

ENDVERBATIM
VERBATIM
      if (ip->dbx>0) 
ENDVERBATIM
{pid() printf("DBC: interval event\n")}  
VERBATIM
      if (ip->record) 
ENDVERBATIM
{ recspk(t) } 
VERBATIM
      if (ip->wrec) 
ENDVERBATIM
{ wrecord(t) } 
    } else {
      tmp = WINV*(1-Vm/EAM)
      VAM = VAM + tmp :: activate interval depolarization
    }
    oinvl=invl
    invlt=t
    net_send(invl,1) 
:** set GABAB weight after a delay
  } else if (flag==5) { : flag==5 to set GABAB weight after a delay
    offsetGB = VGB : current position
    : wflg overloaded; tauGBGP for GB cooperativity
    : wflg will augment each time there is spike coming in through this line
    wflg=wflg*EXP(-(t-tGB)/tauGBGP)+wGBGP 
    coop(wflg)               : cooperativity -- need mult presyn spikes to activate
    : calculate separately based on VGBa and tGB
    if (VGB>EGB) { 
      tmp = wGB*(1-Vm/EGB)*Gn 
      VGB = VGB - tmp
    }
    VGBa= VGB
    tGB=t : restart for VGB
:** flag==2 -- read off external vec (vsp) for next random spike time or single from shock()
  } else if (flag==2) { 
    if (flag==2) { : else for @VERBATIM
VERBATIM
      if (ip->dbx>1) 
ENDVERBATIM
{pid() printf("DBBa: randspk called: %g,%g\n",WEX,nxt)} 
      if (WEX<0) { 
        if (jcn) { jitcon(t) } else { net_event(t) } : bypass activation calculation
VERBATIM
        spikes[ip->type]++; //@

ENDVERBATIM
VERBATIM
        if (ip->dbx>0) 
ENDVERBATIM
{pid() printf("DBB: randspk event\n")} 
        if (WEX<-1 && WEX!=-1e9) { cbur=-WEX  net_send(tbur,4) }
VERBATIM
        if (ip->record) 
ENDVERBATIM
{ recspk(t) } 
VERBATIM
        if (ip->wrec) 
ENDVERBATIM
{ wrecord(t) } 
      } else if (WEX>0) {
        tmp = WEX*(1-Vm/EAM)
        VAM = VAM + tmp
      }
      if (WEX!=-1e9) { : code for single shock
        randspk() : will set WEX for next time
        if (nxt>0) { net_send(nxt,2) }
      }
    }
  } else if (flag==3) { 
    refractory = 0 :end of absolute refractory period    
    trrs = t : save time of start of relative refractory period
VERBATIM
    return; //@ done

ENDVERBATIM
  }
:** check for Vm>VTH -> fire
  Vm = VAM+VNM+VGA+VGB+RMP+AHP : WARNING -- Vm defined differently than above
  if (Vm>0)   {Vm= 0 }
  if (Vm<-90) {Vm=-90}
  if (refractory==0 && Vm>VTHC) {
VERBATIM
    if (!ip->vbr && Vm>Vblock) {//@ do nothing

ENDVERBATIM
VERBATIM
      ip->blkcnt++; blockcnt[ip->type]++; return; }//@

ENDVERBATIM
    AHP = AHP - ahpwt
    tmp=t
    : note that jtt indicates jitter while jit indicates 'just-in-time'
VERBATIM
    if (ip->jttr) 
ENDVERBATIM
{ tmp= t+jttr() }  
    if (jcn) { jitcon(tmp) } else { net_event(tmp) } 
VERBATIM
    spikes[ip->type]++; //@

ENDVERBATIM
VERBATIM
    if (ip->dbx>0) 
ENDVERBATIM
{pid() printf("DBD: %g>VTH(%g) event at %g (STDf=%g)\n",Vm,VTHC,tmp,STDf)} 
VERBATIM
    if (ip->record) 
ENDVERBATIM
{ recspk(tmp) } 
VERBATIM
    if (ip->wrec) 
ENDVERBATIM
{ wrecord(tmp) } 
    VTHC=VTH+RRWght*(Vblock-VTH):increase threshold for relative refrac. period
    VTHR=VTHC :starting thresh value for relative refrac period, keep track of it
    if (nbur>1) { 
      cbur=nbur-1 net_send(tbur,4) 
VERBATIM
      return; //@ done

ENDVERBATIM
    } else if (rebob==1e9) { rebob=0 }
    refractory = 1
VERBATIM
    if (ip->vbr && Vm>Vblock) 
ENDVERBATIM
{ 
      net_send(Vbrefrac,3) 
VERBATIM
      if (ip->dbx>0) 
ENDVERBATIM
{pid() printf("DBE: %g %g\n",Vbrefrac,Vm)} 
VERBATIM
      return; //@ done

ENDVERBATIM
    }
    net_send(refrac, 3)
  }
}

:* ancillary functions
:** jitcon() creates divergence and delays from rand seed
: jcn flags:
: 0 NetCons                            jcn==0
: 3 Jitcon without jitevent            jcn==3 -- eliminated after v669
: 2 Jitcon with callback               jcn==2 -- NOT DEBUGGED
: 1 Jitcon with callback with pointers jcn==1
PROCEDURE jitcon (tm) {
  VERBATIM {
  double mindel, randel, idty, *x; int prty, poty, i, j, k, dv; 
  Point_process *pnt; void* voi;
  // qsz = nrn_event_queue_stats(stt);
  // if (qsz>=qlimit) { printf("qlimit %g exceeded at t=%g\n",qlimit,t); qlimit*=2; }
  ip=IDP;
  if (!pg) {printf("No network defined -- must run jitcondiv() and pgset\n"); hxe();}
  ip->spkcnt++; // jitcon() called from NET_RECEIVE which sets ip
  if (jrmax && jri<jrmax) { 
    jrid[jri]=(double)ip->id; jrtv[jri]=_ltm;
    if (jri>jrsvn) jrsv();
  }
  jri++;  // keep track of number of spikes
  if (jrtm>0) {
    if (t>jrtime) {
      jrtime+=jrtm;
      spkstats2(1.);
    }
  } else if (jrsvd>0 && jri>jrsvn) { 
    jrsvn+=jrsvd; printf("t=%.02f %d ",t,jri);
    spkstats2(1.);
  }
  prty=(int)ip->type;
#if defined(t)
  if (ip->jcn==1) if (ip->dvt>0) net_send((void**)0x0, wts,tpnt,t+ip->del[0],-1.); // first callback
#else
  if (ip->jcn==1) if (ip->dvt>0) net_send((void**)0x0, wts,tpnt,ip->del[0],-1.); // first callback
#endif
  }   
  ENDVERBATIM  
}

: call spkstat from hoc to set global tf if desired
PROCEDURE spkstats () {
VERBATIM {
  if (ifarg(1)) tf=hoc_obj_file_arg(1); else tf=0x0;
}
ENDVERBATIM
}

: flag 1 means print it to a file, 2 means to both places
PROCEDURE spkstats2 (flag) {
VERBATIM {
  int i, ix, flag;
  flag=(int)(_lflag+1e-6);
  if (tf) fprintf(tf,"t=%.02f;%d ",t,jri); else printf("t=%.02f;%d ",t,jri); 
  for (i=0;i<CTYN;i++) {
    ix=cty[i];
    spktot+=spikes[ix];
    if (tf) {
      fprintf(tf,"%s:%d/%d:%d;%d;%d;%d ",CNAME[i],spikes[ix],\
              blockcnt[ix],AMo[ix],NMo[ix],GAo[ix],GBo[ix]);
    } else {
      printf("%s:%d/%d:%d;%d;%d;%d ",CNAME[i],spikes[ix],blockcnt[ix],\
             AMo[ix],NMo[ix],GAo[ix],GBo[ix]);
    }
    blockcnt[ix]=spikes[ix]=0;
    AMo[ix]=NMo[ix]=GAo[ix]=GBo[ix]=0;
  }
  if (tf && flag==2) {  fprintf(tf,"\nt=%g tot_spks: %ld; tot_events: %ld\n",t,spktot,eventtot); 
  } else if (flag==2) {  printf("\ntotal spikes: %ld; total events: %ld\n",spktot,eventtot); 
  } else if (tf) fprintf(tf,"\n"); else printf("\n");
}
ENDVERBATIM
}

PROCEDURE oobpr () {
VERBATIM {
  int i,ix;
  for (i=0;i<CTYN;i++){ 
    ix=cty[i];
    printf("%d:%d/%d:%d;%d;%d;%d ",ix,spikes[ix],blockcnt[ix],AMo[ix],NMo[ix],GAo[ix],GBo[ix]);
  }
  printf("\n");
}
ENDVERBATIM
}

PROCEDURE callback (fl) {
  VERBATIM {
  int i; double idty, del0, ddel; id0 *jp; Point_process *upnt; // these must be local
  i=(unsigned int)((-_lfl)-1); // -1,-2,-3 -> 0,1,2
  jp=IDP; upnt=tpnt; del0=jp->del[i]; ddel=0.;
  idty=(double)(FOFFSET+jp->id)+1e-2*(double)jp->type+1e-3*(double)jp->inhib+1e-4;
  while (ddel<=DELMIN) { // check if this del is worth waiting, else just send now
    if (Vblock<VTHC) { 
      wts[0]=0; // send [0,1] for STD
    } else { // STDf=(1-STD)
      wts[0]=(VTHC-VTH)/(Vblock-VTH); // just send [0,1] for STD
    }
    if (jp->sprob[i]) (*pnt_receive[jp->dvi[i]->_prop->_type])(jp->dvi[i], wts, idty); 
    _p=upnt->_prop->param; _ppvar=upnt->_prop->dparam; // restore pointers
    i++;
    if (i>=jp->dvt) return; // ran out
    ddel=jp->del[i]-del0;   // delays are relative to event; use difference in delays
  }
  // skip over pruned outputs and dead cells:
  while (i<jp->dvt && (!jp->sprob[i] || (*(id0**)&(jp->dvi[i]->_prop->dparam[2]))->dead)) i++;
  if (i<jp->dvt) {
    ddel= jp->del[i] - del0;;
#if defined(t)
    net_send((void**)0x0, wts,upnt,t+ddel,(double) -(i+1)); // next callback
#else
    net_send((void**)0x0, wts,upnt,ddel,(double) -(i+1)); // next callback
#endif
  }
  } 
  ENDVERBATIM
}

: DEAD_DIV not checked in mkdvi()
: mkdvi() create the connectivity vectors for a random network
PROCEDURE mkdvi () {
VERBATIM {
  int i,j,k,prty,poty,dv,dvt,dvii; double *x, *db, *dbs; 
  Object *lb;  Point_process *pnnt, **da, **das;
  ip=IDP; ip->pg=pg; // this should be called right after jitcondiv()
  if (ip->dead) return;
  prty=ip->type;
  sead=GetDVIDSeedVal(ip->id);//seed for divergence and delays
  for (i=0,k=0,dvt=0;i<CTYN;i++) { // dvt gives total divergence
    poty=cty[i];
    dvt+=DVG(prty,poty);
  }
  da =(Point_process **)malloc(dvt*sizeof(Point_process *));
  das=(Point_process **)malloc(dvt*sizeof(Point_process *)); // das,dbs for after sort
  db =(double *)malloc(dvt*sizeof(double)); // delays
  dbs=(double *)malloc(dvt*sizeof(double)); // delays
  for (i=0,k=0,dvii=0;i<CTYN;i++) { // cell types in cty[]
    poty=cty[i];
    dv=DVG(prty,poty);
    if (dv>0) {
      sead+=dv;
      if (dv>dscrsz) {
        printf("B:Divergence exceeds dscrsz: %d>%d for %d->%d\n",dv,dscrsz,prty,poty); hxe(); }
      mcell_ran4(&sead, dscr ,  dv, pg->ixe[poty]-pg->ix[poty]+1);
      for (j=0;j<dv;j++) {
        if (!(lb=ivoc_list_item(ce,(unsigned int)floor(dscr[j]+pg->ix[poty])))) {
          printf("INTF:callback %g exceeds %d for list ce\n",floor(dscr[j]+pg->ix[poty]),cesz); 
          hxe(); }
        pnnt=(Point_process *)lb->u.this_pointer;
        da[j+dvii]=pnnt;
      }
      mcell_ran4(&sead, dscr , dv, 2*DELD(prty,poty));
      for (j=0;j<dv;j++) {
        db[j+dvii]=dscr[j]+DELM(prty,poty)-DELD(prty,poty); // +/- DELD
        if (db[j+dvii]<0) db[j+dvii]=-db[j+dvii];
      }
      dvii+=dv;
    }
  }
  gsort2(db,da,dvt,dbs,das);
  ip->del=dbs;   ip->dvi=das;   ip->dvt=dvt;
  ip->sprob=(unsigned char *)malloc(dvt*sizeof(char *)); // release probability
  for (i=0;i<dvt;i++) ip->sprob[i]=1; // start out with all firing
  free(da); free(db); // keep das,dbs which are assigned to ip->dvi bzw ip->del
  }
ENDVERBATIM
}

:* paths
PROCEDURE patha2b () {
  VERBATIM
  int i; double idty, *x; static Point_process *_pnt; static id0 *ip0;
  pathbeg=*getarg(1); pathidtarg=*getarg(2);
  pathtytarg=-1;  PATHMEASURE=1; pathlen=stopoq=0;
  for (i=0;i<cesz;i++) { lop(ce,i); 
    if ((i==pathbeg || i==pathidtarg) && qp->inhib) {
      pid(); printf("Checking to or from inhib cell\n" ); hxe(); }
    qp->flag=qp->vinflg=0; 
  }
  hoc_call_func(hoc_lookup("finitialize"), 0);
  cvode_fadvance(1000.0); // this call will not return
  ENDVERBATIM
}

:* paths
: pathgrps(vpre,vpos,vout) finds path lengths from pres to posts
FUNCTION pathgrps () {
  VERBATIM
  int i,j,k,na,nb,flag; double idty,*a,*b,*x,sum; static Point_process *_pnt; static id0 *ip0;
  Symbol* s; char **pfl;
  x=0x0;
  s=hoc_lookup("finitialize");
  if (ifarg(2)) {
    na=vector_arg_px(1,&a);
    nb=vector_arg_px(2,&b);
    if (ifarg(3)) x=vector_newsize(vector_arg(3),na*nb);
  } else {
    na=nb=cesz;  // may want to put output into an unsigned char eventually
    if (ifarg(1)) x=vector_newsize(vector_arg(1),na*nb);
  }
  // if (scrsz<cesz) scrset(cesz); 
  pfl = (char **)malloc(cesz * (unsigned)sizeof(char *));
  for (i=0;i<cesz;i++) { lop(ce,i); scr[i]=qp->inhib; pfl[i]=&qp->flag; }
  pathtytarg=-1;  PATHMEASURE=1; pathlen=stopoq=0;
  for (k=0,sum=0;k<na;k++) {
    pathbeg=a[k]; 
    if (scr[(int)pathbeg]) { 
      if (x) for (j=0;j<nb;j++) x[k*nb+j]=0.;
      continue;
    }
    for (j=0;j<nb;j++) { 
      pathidtarg=b[j]; 
      if (scr[(int)pathidtarg]) { if (x) x[k*nb+j]=0.; 
        continue;
      }
      // for (i=0;i<cesz;i++) {lop(ce,i); qp->flag=0;}
      for (i=0;i<cesz;i++) *pfl[i]=0;
      hoc_call_func(s, 0);
      cvode_fadvance(1000.0); // this call will not return
      sum+=pathlen;
      if (x) x[k*nb+j]=pathlen;
    }
  }
  PATHMEASURE=0;
  free(pfl);
  _lpathgrps=sum/na/nb;
  ENDVERBATIM
}

:* intf.getdvi() get divergence
: intf.getdvi(index_vec,delay_vec[,prob_vec,wt1vec,wt2vec]) -- need both wt1vec and wt2vec
: intf.getdvi(getactive.flag,vecs) with flag==1 return types instead of ids
FUNCTION getdvi () {
  VERBATIM 
  {
    int i,j,k,iarg,av1,a2,a3,a4,dvt,getactive=0,idx=0,*pact,prty,poty,sy,ii; 
    double *dbs, *x,*x1,*x2,*x3,*x4,*x5,idty,y[2],flag;
    void* voi, *voi2,*voi3; Point_process **das;
    ip=IDP; ip->pg=pg; // this should be called right after jitcondiv()
    getactive=a2=a3=a4=0;
    if (ip->dead) return;
    dvt=ip->dvt;
    dbs=ip->del;   das=ip->dvi;
    _lgetdvi=(double)dvt; 
    if (!ifarg(1)) return _lgetdvi; // just return the divergence value
    iarg=1;
    if (hoc_is_double_arg(iarg)) {
      av1=2;
      flag=*getarg(iarg++);
      getactive=(int)flag;
      flag-=(double)getactive;
      if (flag!=0) flag=floor(flag*10); // flag is in the decimal place 1.2 has flag of 2
    } else av1=1; // 1st vector arg
    //just get active postsynapses (not dead and non pruned)
    voi=vector_arg(iarg++); 
    if (flag==2.0) { x1=vector_newsize(voi,CTYPi); for (i=0;i<CTYPi;i++) x1[i]=0;
    } else x1=vector_newsize(voi,dvt);
    if (ifarg(iarg)) { voi=vector_arg(iarg++); x2=vector_newsize(voi,dvt);  a2=1; }
    if (ifarg(iarg)) { voi=vector_arg(iarg++); x3=vector_newsize(voi,dvt); a3=1;}
    if (ifarg(iarg)) { // need 2 weight vecs for AM/NM or GA/GB
      voi=vector_arg(iarg++); x4=vector_newsize(voi,dvt); a4=1;
      voi=vector_arg(iarg++); x5=vector_newsize(voi,dvt);
    }
    idty=(double)(FOFFSET+ip->id)+1e-2*(double)ip->type+1e-3*(double)ip->inhib+1e-4;
    prty=ip->type; sy=ip->inhib?GA:AM;
    for (i=0,j=0;i<dvt;i++) {
      qp=*((id0**) &((das[i]->_prop->dparam)[2])); // #define sop	*_ppvar[2].pval
      if (getactive && (qp->dead || ip->sprob[i]==0)) continue;
      if (flag==1.0) { x1[j]=(double)qp->type; 
      } else if (flag==2.0) { x1[qp->type]++; 
      } else x1[j]=(double)qp->id;
      if (a2) x2[j]=dbs[i];
      if (a3) x3[j]=(double)ip->sprob[i];
      if (a4) {
        poty = qp->type;
        if (seadsetting==2) { // no randomization
          for(ii=0;ii<2;ii++) y[ii]=WMAT(prty,poty,sy+ii)*WD0(prty,poty,sy+ii);
        } else {
          if (seadsetting==1) { // old sead setting
            sead=(unsigned int)(FOFFSET+ip->id)*qp->id*seedstep; 
          } else { // hashed sead setting
            hsh[0]=(double)(FOFFSET+ip->id); hsh[1]=(double)(qp->id); hsh[2]=seedstep;
            sead=hashseed2(3,&hsh); 
          }
          mcell_ran4(&sead, y, 2, 1.);
          for(ii=0;ii<2;ii++) {
            y[ii]=2*WVAR*(y[ii]+0.5/WVAR-0.5)*WMAT(prty,poty,sy+ii)*WD0(prty,poty,sy+ii); }
        }
        x4[j]=y[0]; x5[j]=y[1];
      }
      j++;
    }
    if (j!=dvt) for (i=av1;i<iarg;i++) vector_resize(vector_arg(i),j);
    _lgetdvi=(double)j; 
  }
  ENDVERBATIM
}

FUNCTION getconv () {
VERBATIM 
{
  int iarg,i,j,k,dvt,sz,prfl,getactive; double *x,flag;
  void* voi; Point_process **das;
  ip=IDP; ip->pg=pg; // this should be called right after jitcondiv()
  sz=ip->dvt; //  // assume conv similar to div
  getactive=0;
  if (ifarg(iarg=1) && hoc_is_double_arg(iarg)) {
    flag=*getarg(iarg++);
    getactive=(int)flag;
    flag-=(double)getactive;
    if (flag!=0) flag=floor(flag*10); // flag is in the decimal place 1.2 has flag of 2
  }
  if (getactive) printf("WARNING: getactive flag not implemented for getconv()\n");
  if (!ifarg(iarg)) prfl=0; else { prfl=1;
    voi=vector_arg(iarg); 
    if (flag==2.) { x=vector_newsize(voi,CTYPi); for (i=0;i<CTYPi;i++) x[i]=0;
    } else x=vector_newsize(voi,sz); 
  } 
  for (i=0,k=0; i<cesz; i++) {
    lop(ce,i);
    dvt=qp->dvt; das=qp->dvi;
    for (j=0;j<dvt;j++) {
      if (ip==*((id0**) &((das[j]->_prop->dparam)[2]))) {
        if (prfl) {
          if (flag!=2.0 && k>=sz) x=vector_newsize(voi,sz*=2);
          if (flag==1.0) { x[k]=(double)qp->type; 
          } else if (flag==2.0) { x[qp->type]++; 
          } else x[k]=(double)qp->id;
        } 
        k++;
        break;
      }
    }
  }
  if (prfl && flag!=2) vector_resize(voi,k);
  _lgetconv=(double)k;
}
ENDVERBATIM
}

: INTF[0].adjlist(List,[startid,endid,exonly])
: returns adjacency list in first arg
: startid == optional 2nd arg specifies id from which to start
: endid == optional 3rd arg specifies id to end with
: exonly == optional 4th arg specifies to only store excitatory synapse information
FUNCTION adjlist () {
  VERBATIM
  Object* pList = *hoc_objgetarg(1);
  int iListSz=ivoc_list_count(pList),iCell,iStartID=ifarg(2)?*getarg(2):0,\
    iEndID=ifarg(3)?*getarg(3):cesz-1;
  int skipinhib = ifarg(4)?*getarg(4):0, i,j,nv,*pused=(int*)calloc(cesz,sizeof(int)),iSyns=0;
  double **vvo = (double**)malloc(sizeof(double*)*iListSz),\
    *psyns=(double*)calloc(cesz,sizeof(double));
  id0* rp;
  for(iCell=iStartID;iCell<=iEndID;iCell++){
    if(verbose && iCell%1000==0) printf("%d ",iCell);
    lop(ce,iCell);
    if(!qp->dvt || (skipinhib && qp->inhib)){
      list_vector_resize(pList,iCell,0);
      continue;
    }
    iSyns=0;
    for(j=0;j<qp->dvt;j++){      
      rp=*((id0**) &((qp->dvi[j]->_prop->dparam)[2])); // #define sop	*_ppvar[2].pval      
      if(skipinhib && rp->inhib) continue; // if skip inhib cells...
      if(!rp->dead && qp->sprob[j]>0. && !pused[rp->id]){      
        pused[rp->id]=1;
        psyns[iSyns++]=rp->id;
      }
    }
    list_vector_resize(pList, iCell, iSyns);
    list_vector_px(pList, iCell, &vvo[iCell]);
    memcpy(vvo[iCell],psyns,sizeof(double)*iSyns);
    for(j=0;j<iSyns;j++)pused[(int)psyns[j]]=0;
  }
  free(vvo);  free(pused);  free(psyns);
  if (verbose) printf("\n");
  return 1.0;
  ENDVERBATIM
}

FUNCTION rddvi () {
  VERBATIM
  Point_process *pnnt;
  FILE* fp;
  int i, iCell;
  unsigned int iOutID;
  Object* lb;

  fp=hoc_obj_file_arg(1);

  printf("reading: ");

  for(iCell=0;iCell<cesz;iCell++){

    if(iCell%1000==0)printf("%d ",iCell);

    lop(ce,iCell);

    fread(&qp->id,sizeof(unsigned int),1,fp); // read id
    fread(&qp->type,sizeof(unsigned char),1,fp); // read type id
    fread(&qp->col,sizeof(unsigned int),1,fp); // read column id
    fread(&qp->dead,sizeof(unsigned char),1,fp); // read alive/dead status
    fread(&qp->dvt,sizeof(unsigned int),1,fp); // read divergence size

    //free up old pointers
    if(qp->del){
      free(qp->del);
      qp->del=0;
    }
    if(qp->sprob){
      free(qp->sprob);
      qp->sprob=0;
    }
    if(qp->dvi){
      free(qp->dvi);
      qp->dvi=0;
    }

    //if divergence == 0 , continue
    if(!qp->dvt) continue;

    qp->dvi = (Point_process**)malloc(sizeof(Point_process*)*qp->dvt);  
    for(i=0;i<qp->dvt;i++){
      fread(&iOutID,sizeof(unsigned int),1,fp); // id of output cell

      if (!(lb=ivoc_list_item(ce,iOutID))) {
        printf("INTF:callback %g exceeds %d for list ce\n",iOutID,cesz); hxe(); }

      qp->dvi[i]=(Point_process *)lb->u.this_pointer;
    }

    qp->del = (double*)malloc(sizeof(double)*qp->dvt);
    fread(qp->del,sizeof(double),qp->dvt,fp); // read divergence delays
  
    qp->sprob = (unsigned char*)malloc(sizeof(unsigned char)*qp->dvt);
    fread(qp->sprob,sizeof(unsigned char),qp->dvt,fp); // read divergence firing probabilities
  }

  printf("\n");
  
  return 1.0;
  ENDVERBATIM
}

FUNCTION svdvi () {
  VERBATIM
  Point_process *pnnt;
  FILE* fp;
  int i , iCell;
  fp=hoc_obj_file_arg(1);
  printf("writing: ");
  for(iCell=0;iCell<cesz;iCell++){
    if(iCell%1000==0)printf("%d ",iCell);
    lop(ce,iCell);
    fwrite(&qp->id,sizeof(unsigned int),1,fp); // write id
    fwrite(&qp->type,sizeof(unsigned char),1,fp); // write type id
    fwrite(&qp->col,sizeof(unsigned int),1,fp); // write column id
    fwrite(&qp->dead,sizeof(unsigned char),1,fp); // write alive/dead status
    fwrite(&qp->dvt,sizeof(unsigned int),1,fp); // write divergence size
    if(!qp->dvt)continue; //don't write empty pointers if no divergence
    for(i=0;i<qp->dvt;i++){
      pnnt=qp->dvi[i];
      fwrite(&(*(id0**)&(pnnt->_prop->dparam[2]))->id,sizeof(unsigned int),1,fp); // id of output cell
    }
    fwrite(qp->del,sizeof(double),qp->dvt,fp); // write divergence delays
    fwrite(qp->sprob,sizeof(unsigned char),qp->dvt,fp); // write divergence firing probabilities
  }
  printf("\n"); 
  return 1.0;
  ENDVERBATIM
}

: INTF[0].setdvir(wiringlist,delaylist[,flag]) // flag default is 0 to pass to setdvi2()
: INTF[0].setdvir(wiringlist,delaylist,startid,endid)
: INTF[0].setdvir(wiringlist,delaylist,startid,endid,flag)
: INTF[0].setdvir(wiringlist,delaylist,idvec,flag)
: should either use just with flag == 0 to setup all dvi outputs of cells
: or with flag == 1 to incrementally setup outputs from cells and on the last
: set of outputs from a range of cells call with flag == 2 to setup sprob and sort dvi list
: alternatively, can call setdvir with flag == 1, and at end just call INTF.finishdvir to finalize
FUNCTION setdvir () {
  VERBATIM
  ListVec* pListWires,*pListDels;
  int i,dn,flag,dvt,idvfl,iCell,iStartID,iEndID,nidv,end; 
  double *y, *d, *idvec; unsigned char pdead;
  pListWires = AllocListVec(*hoc_objgetarg(1));
  idvfl=flag=0; iStartID=0; iEndID=cesz-1;
  if(!pListWires){printf("setalldvi ERRA: problem initializing wires list arg!\n"); hxe();}
  pListDels = AllocListVec(*hoc_objgetarg(2));
  if(!pListDels){ printf("setalldvi ERRA: problem initializing delays list arg!\n");
    FreeListVec(&pListWires); hxe(); }
  if (ifarg(3) && !ifarg(4)) { 
    flag=(int)*getarg(3); 
  } else if (hoc_is_double_arg(3)) {
    iStartID=(int)*getarg(3);
    iEndID = (int)*getarg(4);
    if(ifarg(5)) flag=(int)*getarg(5);
  } else {
    nidv=vector_arg_px(3, &idvec);
    idvfl=1;
    if (ifarg(4)) flag=(int)*getarg(4);
  }
  end=idvfl?nidv:(iEndID-iStartID+1);
  for (i=0;i<end;i++) {
    if(i%1000==0) printf("%d",i/1000);
    iCell=idvfl?idvec[i]:(iStartID+i);
    lop(ce,iCell);
    if (qp->dead) continue;
    y=pListWires->pv[i]; dvt=pListWires->plen[i];
    if(!dvt) continue; //skip empty div lists
    d=pListDels->pv[i];  dn=pListDels->plen[i];
    if (dn!=dvt) {printf("setdvir() ERR vec sizes for wire,delay list entries not equal %d: %d %d\n",i,dvt,dn); hxe();}
    setdvi2(y,d,dvt,flag);
  }
  FreeListVec(&pListWires);
  FreeListVec(&pListDels);
  return 1.0;
  ENDVERBATIM
}

PROCEDURE clrdvi () {
  VERBATIM
  int i;
  for (i=0;i<cesz;i++) { 
    lop(ce,i);
    free(qp->dvi); free(qp->del); qp->dvt=0;
  }
  ENDVERBATIM
}

: int.setdviv(prevec,postvec,delvec)
FUNCTION setdviv () {
  VERBATIM
  int i,j,k,nprv,dvt; double *prv,*pov,*dlv,x;
  nprv=vector_arg_px(1, &prv);
  i=vector_arg_px(2, &pov);
  j=vector_arg_px(3, &dlv);
  if (nprv!=i || i!=j) {printf("intf:setdviv ERRA: %d %d %d\n",nprv,i,j); hxe();}
  // start by counting the prids so will know the size that we need for realloc()
  if (scrsz<cesz) scrset(cesz); 
  for (i=0;i<cesz;i++) scr[i]=0;
  for (i=0,j=-1;i<nprv;i++) {
    if (j>(int)prv[i]){printf("intf:setdviv ERRA vecs should be sorted by prid vec\n");hxe();}
    j=(int)prv[i];
    scr[j]++;
  }
  for (i=0,x=-1,k=0;i<nprv;i+=dvt) { if(i%1000==0) printf(".");
    if (prv[i]!=x) lop(ce,(unsigned int)(x=prv[i]));
    if (qp->dead) continue;
    dvt=scr[(int)x]; // number of these presyns
    setdvi2(pov+k,dlv+k,dvt,1);
    k+=dvt;
  }
  return (double)k;
  ENDVERBATIM
}

: finishdvi2 () -- finalize dvi , sort dvi , allocate and set sprob
VERBATIM
static int finishdvi2 (struct ID0* p) {
  Point_process **da,**das;
  double *db,*dbs;
  int i, dvt;
  db=p->del;
  da=p->dvi; 
  dvt=p->dvt;
  dbs=(double*)malloc(dvt*sizeof(double)); // sorted delays
  das=(Point_process**)malloc(dvt*sizeof(Point_process*)); // parallel sorted dvi
  gsort2(db,da,dvt,dbs,das);
  p->del=dbs; p->dvi=das; // sorted versions
  free(db); free(da);
  p->sprob=(unsigned char*)realloc((void*)p->sprob,(size_t)dvt*sizeof(char));// release probability
  for (i=0;i<dvt;i++) p->sprob[i]=1; // start out with all firing
}
ENDVERBATIM

: finalize dvi for all cells
PROCEDURE finishdvir () {
  VERBATIM
  int iCell;
  for(iCell=0;iCell<cesz;iCell++){
    lop(ce,iCell);
    finishdvi2(qp);
  }
  ENDVERBATIM
}

: finishdvi() -- finalize dvi , sort dvi, allocate and set sprob, for this single cell
PROCEDURE finishdvi () {
VERBATIM
  finishdvi2(IDP);
ENDVERBATIM
}

PROCEDURE setdvi () {
VERBATIM {
  int i,j,dvt,ddvi,flag; double *d, *y;
  if (! ifarg(1)) {printf("setdvi(v1,v2[,flag]): v1:cell#s; v2:delays\n"); return; }
  ip=IDP; ip->pg=pg; // this should be called right after jitcondiv()
  if (ip->dead) return;
  dvt=vector_arg_px(1, &y);
  i=vector_arg_px(2, &d);
  // flag 1: grow internal vecs; flag 2: grow and do final sort
  if (ifarg(3)) flag=(int)*getarg(3); else flag=0;
  if (i!=dvt || i==0) {printf("setdvi() ERR vec sizes: %d %d\n",dvt,j); hxe();}
  setdvi2(y,d,dvt,flag);
  }
ENDVERBATIM
}

VERBATIM
// setdvi2(divid_vec,del_vec,div_cnt,flag)
// flag 1 means just augment, 0or2: sort by del, 0: clear lists and replace
static int setdvi2 (double *y,double *d,int dvt,int flag) {
  int i,j,ddvi; double *db, *dbs; unsigned char pdead; unsigned int b,e;
  Object *lb; Point_process *pnnt, **da, **das;
  ddvi=(int)DEAD_DIV;
  ip=IDP; 
  if (flag==0) { b=0; e=dvt; // begin to end
    if (ip->dvi) { free(ip->dvi); free(ip->del); ip->dvi=0; ip->del=0; }
  } else {
    if (ip->dvt==0) {
      ip->dvi=(Point_process **)0x0; 
      ip->del=(double*)0x0; // make sure they're null pointers for realloc
    }
    b=ip->dvt; 
    e=ip->dvt+dvt; // dvt is amount to grow
  }
  da=(Point_process **)realloc((void*)ip->dvi,(size_t)(e*sizeof(Point_process *)));
  db=(double*)realloc((void*)ip->del,(size_t)(e*sizeof(double)));
  for (i=b,j=0;j<dvt;j++) { // i thru da[] j thru y, k to append
    // div can grow at lower rate if dead cells are encountered
    if (!(lb=ivoc_list_item(ce,(unsigned int)y[j]))) {
      printf("INTF:callback %g exceeds %d for list ce\n",y[j],cesz); hxe(); }
      pnnt=(Point_process *)lb->u.this_pointer;
      if (ddvi==1 || !(pdead=(*(id0**)&(pnnt->_prop->dparam[2]))->dead)) {
        da[i]=pnnt; db[i]=d[j]; i++;
      }
  }
  if ((dvt=i)<e) { // will need to shrink these arrays
    da=(Point_process **)realloc((void*)da,(size_t)(e*sizeof(Point_process *)));
    db=(double*)realloc((void*)db,(size_t)(e*sizeof(double)));
  }
  ip->dvt=dvt; ip->del=db; ip->dvi=da;
  if (flag!=1) finishdvi2(ip); // do sort
}
ENDVERBATIM

: prune(p[,potype,rand_seed]) // prune synapses with prob p [0,1], ie 0.1 prunes 10% of the divergence
: prune(vec) // fill in the pruning vec with binary values from vec
PROCEDURE prune () {
  VERBATIM 
  {
  id0* ppost; double *x, p; int nx,j,potype;
  ip=IDP; ip->pg=pg;
  if (hoc_is_double_arg(1)) {
    p=*getarg(1);
    if (p<0 || p>1) {printf("INTF:pruneERR0:need # [0,1] to prune [ALL,NONE]: %g\n",p); hxe();}
    if (p==1.) printf("INTFpruneWARNING: pruning 100% of cell %d\n",ip->id);
    if (ip->dvt>dscrsz) {
      printf("INTFpruneB:Div exceeds dscrsz: %d>%d\n",ip->dvt,dscrsz); hxe(); }
    if (p==0.) {
      for (j=0;j<ip->dvt;j++) ip->sprob[j]=1; // unprune completely
      return; // now that unpruning is done, can return
    }
    potype=ifarg(2)?(int)*getarg(2):-1;
    sead=(ifarg(3))?(unsigned int)*getarg(3):GetDVIDSeedVal(ip->id);//seed for divergence and delays
    mcell_ran4(&sead, dscr , ip->dvt, 1.0); // random var (0,1)
    if(potype==-1){ // prune all types of synapses
      for (j=0;j<ip->dvt;j++) if (dscr[j]<p) ip->sprob[j]=0; // prune with prob p
    } else { // only prune synapses with postsynaptic type == potype
      for (j=0;j<ip->dvt;j++){
        ppost=*((id0**) &((ip->dvi[j]->_prop->dparam)[2])); // #define sop *_ppvar[2].pval
        if (ppost->type==potype && dscr[j]<p) ip->sprob[j]=0; // prune with prob p
      }
    }
  } else {
    nx=vector_arg_px(1,&x);
    if (nx!=ip->dvt) {printf("INTF:pruneERRA:Wrong size vector:%d!=%d\n",nx,ip->dvt); hxe();}
    for (j=0;j<ip->dvt;j++) ip->sprob[j]=(unsigned char)x[j];
  }
  }
ENDVERBATIM
}

: turnoff(v1,v2) turn off any connection from a cell in v1 to a cell with number in v2
: a global call that can be called from any INTF
PROCEDURE turnoff () {
  VERBATIM {
  int nx,ny,i,j,k,dvt; double poid,*x,*y; Point_process **das; unsigned char off;
  nx=vector_arg_px(1,&x);
  ny=vector_arg_px(2,&y);
  if (ifarg(3)) off=(unsigned char)*getarg(3); else off=0;
  for (i=0;i<nx;i++) { 
    lop(ce,(unsigned int)x[i]); 
    dvt=qp->dvt; das=qp->dvi;
    for (j=0;j<dvt;j++) {
      ip=*((id0**) &((das[j]->_prop->dparam)[2])); // sop is *_ppvar[2].pval
      poid=(double)ip->id; // postsyn id
      for (k=0;k<ny;k++) {
        if (poid==y[k]) {
          qp->sprob[j]=off; break;
        }
      }
    }
  }
  }
  ENDVERBATIM
}

VERBATIM 
// gsort2() sorts 2 parallel vectors -- delays and Point_process pointers
int gsort2 (double *db, Point_process **da,int dvt,double *dbs, Point_process **das) {
  int i;
  scr=scrset(dvt);
  for (i=0;i<dvt;i++) scr[i]=i;
  nrn_mlh_gsort(db, scr, dvt, cmpdfn);
  for (i=0;i<dvt;i++) {
    dbs[i]=db[scr[i]];
    das[i]=da[scr[i]];
  }
}
ENDVERBATIM

PROCEDURE freedvi () {
  VERBATIM
  { 
    int i, poty; id0 *jp;
    jp=IDP;
    if (jp->dvi) {
      free(jp->dvi);
      free(jp->del);
      jp->dvi=0x0;
      jp->del=0x0;
    }
  }
  ENDVERBATIM
}

: pgset() should only be needed temporarily
PROCEDURE pgset () {
  VERBATIM
  ip->pg=pg; // set the local to the global
  ENDVERBATIM
}

FUNCTION qstats () {
  VERBATIM {
    double stt[3]; int lct,flag; FILE* tfo;
    if (ifarg(1)) {tfo=hoc_obj_file_arg(1); flag=1;} else flag=0;
    lct=cty[IDP->type];
    _lqstats = nrn_event_queue_stats(stt);
    printf("SPIKES: %d (%ld:%ld)\n",IDP->spkcnt,spikes[lct],blockcnt[lct]);
    printf("QUEUE: Inserted %g; removed %g\n",stt[0],stt[2]);
    if (flag) {
      fprintf(tfo,"SPIKES: %d (%ld:%ld);",IDP->spkcnt,spikes[lct],blockcnt[lct]);
      fprintf(tfo,"QUEUE: Inserted %g; removed %g remaining: %g\n",stt[0],stt[2],_lqstats);
    }
  }
  ENDVERBATIM
}

FUNCTION qsz () {
  VERBATIM {
    double stt[3];
    _lqsz = nrn_event_queue_stats(stt);
  }
  ENDVERBATIM
}

PROCEDURE qclr () {
  VERBATIM {
    clear_event_queue();
  }
  ENDVERBATIM
}

: intf.jitcondiv() assigns pointers for hoc symbol storage
PROCEDURE jitcondiv () {
  VERBATIM {
  Symbol *sym; int i,j; unsigned int sz; char *name;
  pg=(postgrp *)malloc(sizeof(postgrp));
  sym = hoc_lookup("ce"); ce = (*(hoc_objectdata[sym->u.oboff].pobj));
  sym = hoc_lookup("CTYP"); CTYP = (*(hoc_objectdata[sym->u.oboff].pobj));
  if (installed==2.0) { // jitcondiv was previously run
    sz=ivoc_list_count(ce);
    if (sz==cesz) printf("\t**** INTF WARNING cesz unchanged: INTF(s) created off-list ****\n");
  } else installed=2.0;
  cesz = ivoc_list_count(ce);
  cty[0]=DP; cty[1]=SU; cty[2]=IN; // set the cell types
  CTYPi=HVAL("CTYPi"); STYPi=HVAL("STYPi"); dscrsz=HVAL("scrsz");
  pg->ix =HPTR("ix"); pg->ixe=HPTR("ixe"); 
  pg->dvg=HPTR("div"); pg->numc=HPTR("numc"); 
  pg->wmat=HPTR("wmat"); pg->wd0=HPTR("wd0");
  pg->delm=HPTR("delm"); pg->deld=HPTR("deld");
  dscr=HPTR("scr");
  if (!ce) {printf("INTF jitcondiv ERRA: ce not found\n"); hxe();}
  if (ivoc_list_count(CTYP)!=CTYPi){
    printf("INTF jitcondiv ERRB: %d %d\n",ivoc_list_count(CTYP),CTYPi); hxe(); }
  for (i=0;i<cesz;i++) { lop(ce,i); qp->pg=pg; } // set all of the pg pointers for now
  // make sure no seg error:
  printf("Checking for possible seg error in double arrays: CTYPi==%d: ",CTYPi);
  // can access arbitrary member dvg[a][b] using (&dvg[a*CTYPi])[b] or dvg+a*CTYPi+b
  printf("%d %d %d ",DVG(CTYPi-1,CTYPi-1),(int)pg->ix[CTYPi-1],(int)pg->ixe[CTYPi-1]);
  printf("%g %g ",WMAT(CTYPi-1,CTYPi-1,STYPi-1),WD0(CTYPi-1,CTYPi-1,STYPi-1));
  printf("%g %g ",DELM(CTYPi-1,CTYPi-1),DELD(CTYPi-1,CTYPi-1));
  printf("%d %g\n",dscrsz,dscr[dscrsz-1]); // scratch area for doubles
  for (i=0,j=0;i<CTYPi;i++) if (ctt(i,&name)!=0) {
    cty[j]=i; CNAME[j]=name;
    j++;
    if (j>=CTYPp) {printf("jitcondiv() INTERRA\n"); hxe();}
  }
  CTYN=j; // number of cell types being used
  for (i=0;i<CTYN;i++) printf("%s(%d)=%g ",CNAME[i],cty[i],NUMC(cty[i]));
  printf("\n%d cell types being used\n",CTYN);
  }
  ENDVERBATIM  
}

: intf.jitrec(vec,tvec)
PROCEDURE jitrec () {
  VERBATIM {
  int i;
  if (! ifarg(2)) { // clear with jitrec() or jitrec(0)
    jrmax=0; jridv=0x0; jrtvv=0x0;
    return;
  }
  i =   vector_arg_px(1, &jrid); // could just set up the pointers once
  jrmax=vector_arg_px(2, &jrtv);
  jridv=vector_arg(1); jrtvv=vector_arg(2);
  jrmax=vector_buffer_size(jridv);
  if (jrmax!=vector_buffer_size(jrtvv)) {
    printf("jitrec() ERRA: not same size: %d %d\n",i,jrmax); jrmax=0; hxe(); }
  jri=0; // needs to be set at beginning of run
  }
  ENDVERBATIM
}

: intf.jrsv()
PROCEDURE jrsv () {
  VERBATIM {
  char fnm[50]; int ty=4;
  sprintf(fnm,"/tmp/aa%d.vec",process);
  unlink(fnm);
  if ( !(wf1 = fopen(fnm,"w"))) { printf("Can't open %s\n",fnm); hxe(); }
  fwrite(&jri,sizeof(int),1,wf1);
  fwrite(&ty,sizeof(int),1,wf1);
  fwrite(jrid,sizeof(double),jri,wf1);
  fwrite(&jri,sizeof(int),1,wf1);
  fwrite(&ty,sizeof(int),1,wf1);
  fwrite(jrtv,sizeof(double),jri,wf1);
  fclose(wf1);
  }
  ENDVERBATIM
}

: intf.scsv()
FUNCTION scsv () {
  VERBATIM {
  int ty=4; int i,j; unsigned int cnt=0;
  name = gargstr(1);
  if ( !(wf1 = fopen(name,"w"))) { printf("Can't open %s\n",name); hxe(); }
  fwrite(&cesz,sizeof(int),1,wf1);
  fwrite(&ty,sizeof(int),1,wf1);
  for (i=0,j=0;i<cesz;i++,j++) { 
    lop(ce,i); 
    if (qp->spkcnt) {
      dscr[j]=(double)(qp->spkcnt); 
      cnt++;
    } else dscr[j]=0.0;
    if (j>=dscrsz) {
      fwrite(dscr,(size_t)sizeof(double),(size_t)dscrsz,wf1);
      fflush(wf1);
      j=0;
    }
  }
  if (j>0) fwrite(dscr,(size_t)sizeof(double),(size_t)j,wf1);
  fclose(wf1);
  _lscsv=(double)cnt;
  }
  ENDVERBATIM
}

: intf.spkcnt(vec[,vec,flag])
: intf.spkcnt(min,max[,vec,flag]) flag=1 means reset all counts to 0
FUNCTION spkcnt () {
  VERBATIM {
  int nx, ny, i,j, ix, c, min, max, flag; unsigned int sum; double *y,*x;
  nx=ny=min=max=flag=0; i=1;
  if (ifarg(i)) {
    if (hoc_is_object_arg(i)) { 
      ny = vector_arg_px(i, &y); i++;
    } else if (ifarg(i+1)) {
      min=(int)*getarg(i); max=(int)*getarg(i+1); i+=2;
    }
  }
  while (ifarg(i)) { // can pick up flag and vector in either order
    if (hoc_is_object_arg(i)) { // output to a vector
      nx = vector_arg_px(i, &x);
    } else flag=(int)*getarg(i);
    i++;
  }
  if (ny) max=ny; else if (max==0) max=cesz; else max+=1; // enter max index wish to graph
  if (nx && nx!=max-min) {
    printf("INTF spkcnt() ERR: Vectors not same size %d %g\n",nx,max-min);hxe();}
  for  (i=min, sum=0;i<max;i++) { 
    if (ny) lop(ce,(int)y[i]); else lop(ce,i);
    if (flag==2) sum+=(c=qp->blkcnt); else sum+=(c=qp->spkcnt);
    if (nx) x[i]=(double)c;
    if (flag==1) qp->spkcnt=qp->blkcnt=0;
  }
  _lspkcnt=(double)sum;
  }
  ENDVERBATIM
}

:** probejcd()
PROCEDURE probejcd () {
  VERBATIM {  int i,a[4];
    for (i=1;i<=3;i++) a[i]=(int)*getarg(i);
    printf("CTYPi: %d, STYPi: %d, ",CTYPi,STYPi);
    // printf("div: %d, ix: %d, ixe: %d, ",DVG(a[1],a[2]),(int)ix[a[1]],(int)ixe[a[1]]);
    printf("wmat: %g, wd0: %g\n",WMAT(a[1],a[2],a[3]),WD0(a[1],a[2],a[3]));
  }
  ENDVERBATIM  
}

:** randspk() sets next to next val in vector, this vector is handled globally
PROCEDURE randspk () {
  VERBATIM 
  ip=IDP;  
  if (ip->rvi > ip->rve) { // pointers go from rvi to rve inclusive
    ip->input=0;           // turn off
    nxt=-1.;
  } else if (t==0) {     // initialization
    nxt=vsp[ip->rvi];
    WEX=wsp[ip->rvi++];
  } else {     // absolute times in vector -> interval
    while ((nxt=vsp[ip->rvi++]-t)<=1e-6) { 
      if (ip->rvi-1 > ip->rve) { printf("randspk() ERRA: "); chk(2.); hxe(); }
    }
    WEX=wsp[ip->rvi-1]; // rvi was incremented
  }
  ENDVERBATIM
  : net_send(nxt,2) : can only be called from INITIAL or NET_RECEIVE blocks
}

:** vers gives version
PROCEDURE vers () {
  printf("$Id: intf.mod,v 1.789 2008/12/20 16:04:59 billl Exp $\n")
}

:** val(t,tstart) fills global vii[] to pass values back to record() (called from record())
VERBATIM
double val (double xx, double ta) { 
  vii[1]=VAM*EXP(-(xx - ta)/tauAM);
  vii[2]=VNM*EXP(-(xx - ta)/tauNM);
  vii[3]=VGA*EXP(-(xx - ta)/tauGA);
  if (VGBdel>0) {
    vii[4]=esinr(xx-tGB);
  } else {
    vii[4]=VGB*EXP(-(xx - ta)/tauGB);
  }  
  vii[5]=AHP*EXP(-(xx - ta)/tauahp);
  vii[6]=vii[1]+vii[2]+vii[3]+vii[4]+vii[5];
  vii[7]=VTH + (VTHR-VTH)*EXP(-(xx-trrs)/tauRR);
}
ENDVERBATIM

:** valps(t,tstart) like val but builds voltages for pop spike
VERBATIM
double valps (double xx, double ta) { 
  vii[1]=VAM*EXP(-(xx - ta)/tauAM);
  vii[2]=VNM*EXP(-(xx - ta)/tauNM);
  vii[3]=VGA*EXP(-(xx - ta)/tauGA);
  vii[4]=esinr(xx-tGB);
  // vii[5]=AHP*EXP(-(xx - ta)/tauahp);
  vii[6]=vii[1]+vii[2]-vii[3];
}
ENDVERBATIM

:** record() stores values since last tg into appropriate vecs
PROCEDURE record () {
  VERBATIM {
  int i,j,k,nz; double ti;
  vp = SOP;
  if (tg>=t) return;
  if (ip->record==1) {
    while ((int)vp->p >= (int)vp->size-(int)((t-tg)/vdt)-10) { 
      vp->size*=2;
      for (k=0;k<NSV;k++) if (vp->vv[k]!=0x0) vp->vvo[k]=vector_newsize(vp->vv[k], vp->size);
      // printf("**** WARNING expanding recording room to %d (type%d id%d at %g)****\n",vp->size,IDP->type,IDP->id,t);
    }
  } else if ((int)vp->p > (int)vp->size-(int)((t-tg)/vdt)) { // shift if record==2
    nz=(int)((t-tg)/vdt);
    for (k=0;k<NSV;k++) if (vp->vv[k]!=0x0) {
      if (nz>vp->size) {pid(); printf("Record WARNING: vec too short: %d %d\n",nz,vp->size);
        vp->p=0;
      } else {
        for (i=nz,j=0; i<vp->size; i++,j++) vp->vvo[k][j]=vp->vvo[k][i];
        vp->p=vp->size-nz;
      }
    }
  }
  for (ti=tg;ti<=t && vp->p < vp->size;ti+=vdt,vp->p++) { 
    val(ti,tg);  
    if (vp->vvo[0]!=0x0) vp->vvo[0][vp->p]=ti;
    for (k=1;k<NSV-1;k++) if (vp->vvo[k]!=0x0) { // not nil pointer
      vp->vvo[k][vp->p]=vii[k]+RMP;
    }
    for (;k<NSV;k++) if (vp->vvo[k]!=0x0) { // not nil pointer
      vp->vvo[k][vp->p]=vii[k]; 
    }
  }
  tg=t;
  }
  ENDVERBATIM
}

:** recspk() records a spike by writing a 10 into the main VM vector
PROCEDURE recspk (x) {
  VERBATIM { int k;
  vp = SOP;
  record();
  if (vp->p > vp->size || vp->vvo[6]==0) return; 
  if (vp->vvo[0]!=0x0) vp->vvo[0][vp->p-1]=_lx;
  vp->vvo[6][vp->p-1]=spkht; // the spike
  tg=_lx;
  }
  ENDVERBATIM
}

:** recclr() clear the vectors pointers
PROCEDURE recclr () {
  VERBATIM 
  {int k;
  if (IDP->record) {
    if (SOP!=nil) {
      vp = SOP;
      vp->size=0; vp->p=0;
      for (k=0;k<NSV;k++) { vp->vv[k]=nil; vp->vvo[k]=nil; }
    } else printf("INTF recclr ERR: nil pointer\n");
  }
  IDP->record=0;
  }
  ENDVERBATIM 
}

:** recfree() free the vpt pointer memory
PROCEDURE recfree () {
  VERBATIM
  if (SOP!=nil) {
    free(SOP);
    SOP=nil;
  } else printf("INTF recfree ERR: nil pointer\n");
  IDP->record=0;
  ENDVERBATIM
}

:** initvspks() sets up vector from which to read random spike times 
: this is a global procedure to set up pieces of a global vector
: all cells share one vector but read from different locations
: (CHANGED from intervals and global proc in v224)
PROCEDURE initvspks () {
  VERBATIM 
  {int max, i,err;
    double last,lstt;
    if (! ifarg(1)) {printf("Return initvspks(ivspks,vspks,wvspks)\n"); return 0.;}
    if (isp!=NULL) clrvspks();
    ip=IDP;  err=0;
    i = vector_arg_px(1, &isp); // could just set up the pointers once
    max=vector_arg_px(2, &vsp);
    if (max!=i) {err=1; printf("initvspks ERR: vecs of different size\n");}
    if (max==0) {err=1; printf("initvspks ERR: vec not initialized\n");}
    max=vector_arg_px(3, &wsp);
    if (max!=i) {err=1; printf("initvspks ERR: 3rd vec is of different size\n");}
    vspn=max;
    if (!ce) {printf("Need global ce for initvspks() since intf.mod501\n"); hxe();}
    for (i=0,last=-1; i<max; ) { // move forward to first
      if (isp[i]!=last) { // new one
        lop(ce,(unsigned int)isp[i]);
        qp->rvb=qp->rvi=i;
        qp->vinflg=1;
        last=isp[i];
        lstt=vsp[i];
        i++;
      }
      for (; i<max && isp[i] == last; i++) { // move forward to last
        if (vsp[i]<=lstt) { err=1; 
          printf("initvspks ERR: nonmonotonic for cell#%d: %g %g\n",qp->id,lstt,vsp[i]); }
          lstt=vsp[i];
      }
      qp->rve=i-1;
      if (subsvint>0) { 
        vsp[qp->rve] = vsp[qp->rvb]+subsvint;
        wsp[qp->rve] = wsp[qp->rvb];
      }
      if (err) { qp->rve=0; hxe(); }
    }
  }
  ENDVERBATIM
}

:** shock() reads random spike times from save db as initvspks() but just sends a single shock
: to each listed cell
: this is a global procedure that calls multiple cells
PROCEDURE shock () {
  VERBATIM 
  {int max, i,err;
    double last, lstt, *isp, *vsp, *wsp;
    if (! ifarg(1)) {printf("Return shock(ivspks,vspks,wvspks)\n"); return 0.;}
    ip=IDP;  err=0;
    i = vector_arg_px(1, &isp); // could just set up the pointers once
    max=vector_arg_px(2, &vsp);
    if (max!=i) {err=1; printf("shock ERR: vecs of different size\n");}
    if (max==0) {err=1; printf("shock ERR: vec not initialized\n");}
    max=vector_arg_px(3, &wsp);
    if (max!=i) {err=1; printf("shock ERR: 3rd vec is of different size\n");}
    vspn=max;
    if (!ce) {printf("Need global ce for shock()\n"); hxe();}
    for (i=0,last=-1; i<max; ) { // move forward to first
      if (isp[i]!=last) { // skip any redund indices
        lop(ce,(unsigned int)isp[i]);
        WEX=-1e9; // code for shock
#if defined(t)
        net_send((void**)0x0, wts,pmt,t+vsp[i],2.0); // 2 is randspk flag
#else
        net_send((void**)0x0, wts,pmt,vsp[i],2.0); // 2 is randspk flag
#endif
        i++;
      }
    }
  }
  ENDVERBATIM
}

PROCEDURE clrvspks () {
 VERBATIM {
 unsigned int i;
 for (i=0; i<cesz; i++) {
   lop(ce,i);
   qp->vinflg=0;
 }   
 }
 ENDVERBATIM
}

: trvsp gets called globally to go through the vector
: first pass (arg 1) it replaces terminal values with 1e9
: second pass (arg 2) it replaces terminal values with first+subsvint
PROCEDURE trvsp ()
{
  VERBATIM 
  int i, flag; 
  double ind, t0;
  ip=IDP;
  flag=(int) *getarg(1);
  if (subsvint==0.) {printf("trvsp"); return(0.);}
  ind=isp[0]; t0=vsp[0];
  if (flag==1) {
    for (i=0; i<vspn; i++) {
      if (isp[i]!=ind) {
        vsp[i-1]=1.e9;
        ind=isp[i];
      }
    }
    vsp[vspn-1]=1.e9;
  } else if (flag==2) {
    for (i=0; i<vspn; i++) {
      if (isp[i]!=ind) {
        vsp[i-1]=t0+subsvint;
        ind=isp[i]; t0=vsp[i];
      }
    }
    vsp[vspn-1]=t0+subsvint;
  } else {printf("trvsp flag %d not recognized\n",flag); hxe();}
  ENDVERBATIM
}

:** initjttr() sets up vector from which to read jitter 
: -- key jtt to avoid confusion with jitcon=='just in time connection'
: this is a global not a range procedure -- just call once
PROCEDURE initjttr () {
  VERBATIM 
  {int max, i, err=0;
    jtpt=0;
    if (! ifarg(1)) {printf("Return initjttr(vec)\n"); return(0.);}
    max=vector_arg_px(1, &jsp);
    if (max==0) {err=1; printf("initjttr ERR: vec not initialized\n");}
    for (i=0; i<max; i++) if (jsp[i]<=0) {err=1;
      printf("initjttr ERR: vec should be >0: %g\n",jsp[i]);}
    if (err) { jsp=nil; jtmax=0.; return(0.); }// hoc_execerror("",0);
    if (max != jtmax) {
      printf("WARNING: resetting jtmax_INTF to %d\n",max); jtmax=max; }
  }
  ENDVERBATIM
}

:* internal routines
VERBATIM
//** lop(LIST,ITEM#) sets qp
// modeled on vector_arg_px(): picks up obj from list and resolves pointers
static double* lop (Object *ob, unsigned int i) {
  Object *lb;
  lb = ivoc_list_item(ob, i);
  if (! lb) { printf("INTF:lop %d exceeds %d for list ce\n",i,cesz); hxe();}
  pmt=ob2pntproc(lb);
  qp=*((id0**) &((pmt->_prop->dparam)[2])); // #define sop *_ppvar[2].pval
  // _hoc_setdata((void*)pmt); // make all the range vars accessible
  _p=pmt->_prop->param; _ppvar=pmt->_prop->dparam;
  return pmt->_prop->param;
}

// use stoppo() as a convenient conditional breakpoint in gdb (gdb watching is too slow)
int stoppo () {
}

//** ctt(ITEM#) find cells that exist by name
static int ctt (unsigned int i, char** name) {
  Object *lb;
  if (NUMC(i)==0) return 0; // none of this cell type
  lb = ivoc_list_item(CTYP, i);
  if (! lb) { printf("INTF:ctt %d exceeds %d for list CTYP\n",i,CTYPi); hxe();}
  {*name=*(lb->u.dataspace->ppstr);}
  return (int)NUMC(i);
}
ENDVERBATIM


PROCEDURE test () {
  VERBATIM
  char *str; int x;
  x=ctt(7,&str); 
  printf("%s (%d)\n",str,x);
  ENDVERBATIM
}

: lof can find object information
PROCEDURE lof () {
VERBATIM {
  Object *ob; int num,i,ii,j,k,si,nx;  double *vvo[7], *par; void *vv[7];
  ob = *(hoc_objgetarg(1));
  si=(int)*getarg(2);
  num = ivoc_list_count(ob);
  if (num!=7) { printf("INTF lof ERR %d>7\n",num); hxe(); }
  for (i=0;i<num;i++) { 
    j = list_vector_px3(ob, i, &vvo[i], &vv[i]);
    if (i==0) nx=j;
    if (j!=nx) { printf("INTF lof ERR %d %d\n",j,nx); hxe(); }
  }
  //  for (i=ix[si],ii=0;i<=ixe[si] && ii<nx;i++,ii++) {
  //   vvo[0][ii]=(double)i;
  //   par=lop(ce,i);
  //   for (j=20,k=1;j<25;j++,k++) { // NB these could move: Vm,VAM,VNM,VGA,VGB
  //     vvo[k][ii]=par[j];
  //   }
  // }
 }
ENDVERBATIM
}

:* initinvl() sets up vector from which to read intervals
: this is a global not a range procedure -- just call once
PROCEDURE initinvl () {
  printf("initinvl() NOT BEING USED\n")
}

: invlflag() used internally; can't set from here; use initinvl() and range invlset()
FUNCTION invlflag () {
  VERBATIM
  ip=IDP;
  if (ip->invl0==1 && invlp==nil) { // err
    printf("INTF invlflag ERR: pointer not initialized\n"); hoc_execerror("",0); 
  }
  _linvlflag= (double)ip->invl0;
  ENDVERBATIM
}

:** shift() returns the appropriate shift
FUNCTION shift (vl) { 
  VERBATIM   
  double expand, tmp, min, max;
//if (invlp==nil) {printf("INTF invlflag ERRa: pointer not initialized\n"); hoc_execerror("",0);}
  if ((t<(invlt-invl)+invl/2) && invlt != -1) { // don't shift if less than halfway through
    _lshift=0.;  // flag for no shift
  } else {
    expand = -(_lvl-(-65))/20; // expand positive if hyperpolarized
    if (expand>1.) expand=1.; if (expand<-1.) expand=-1.;
    if (expand>0.) { // expand interval
      max=1.5*invl;
      tmp=oinvl+0.8*expand*(max-oinvl); // the amount we can add to the invl
    } else {
      min=0.5*invl; 
      tmp=oinvl+0.8*expand*(oinvl-min); // the amount we can reduce current invl
    }
    if (invlt+tmp<t+2) { // getting too near spike time
      _lshift=0.;
    } else {
      oinvl=tmp; // new interval
      _lshift=invlt+oinvl;
    }
  }
  ENDVERBATIM
}

:* recini() called from INITIAL block to set vp->p to zero and open up vectors
PROCEDURE recini () {
  VERBATIM 
  { int k;
  if (SOP==nil) { 
    printf("INTF record ERR: pointer not initialized\n"); hoc_execerror("",0); 
  } else {
    vp = SOP;
    vp->p=0;
    // open up the vector maximally before writing into it; will correct size in fini
    for (k=0;k<NSV;k++) if (vp->vvo[k]!=0) vector_resize(vp->vv[k], vp->size);
  }}
  ENDVERBATIM
}

:** fini() to finish up recording -- should be called from FinishMisc()
PROCEDURE fini () {
  VERBATIM 
  {int k;
  // initialization for next round, this will not be set if job terminates prematurely
  IDP->rvi=IDP->rvb;  // -- see vinset()
  if (IDP->wrec) { wrecord(1e9); }
  if (IDP->record) {
    record(); // finish up
    for (k=0;k<NSV;k++) if (vp->vvo[k]!=0) { // not nil pointer
      vector_resize(vp->vv[k], vp->p);
    }
  }}
  ENDVERBATIM
}

:** chk([flag]) with flag=1 prints out info on the record structure
:                    flag=2 prints out info on the global vectors
PROCEDURE chk (f) {
  VERBATIM 
  {int i,lfg;
  lfg=(int)_lf;
  ip=IDP;
  printf("ID:%d; typ: %d; rec:%d wrec:%d inp:%d jtt:%d invl:%d\n",ip->id,ip->type,ip->record,ip->wrec,ip->input,ip->jttr,ip->invl0);
  if (lfg==1) {
    if (SOP!=nil) {
      vp = SOP;
      printf("p %d size %d tg %g\n",vp->p,vp->size,tg);
      for (i=0;i<NSV;i++) if (vp->vv[i]) printf("%d %x %x;",i,vp->vv[i],vp->vvo[i]);
    } else printf("Recording pointers not initialized");
  }
  if (lfg==2) { 
    printf("Global vectors for input and jitter (jttr): \n");
    if (vsp!=nil) printf("VSP: %x (%d/%d-%d)\n",vsp,ip->rvi,ip->rvb,ip->rve); else printf("no VSP\n");
    if (jsp!=nil) printf("JSP: %x (%d/%d)\n",jsp,jtpt,jtmax); else printf("no JSP\n");
  }
  if (lfg==3) { 
    if (vsp!=nil) { printf("VSP: (%d/%d-%d)\n",ip->rvi,ip->rvb,ip->rve); 
      for (i=ip->rvb;i<=ip->rve;i++) printf("%d:%g  ",i,vsp[i]);
      printf("\n");
    } else printf("no VSP\n");
  }
  if (lfg==4) {  // was used to give invlp[],invlmax
  }
  if (lfg==5) { 
    printf("wwpt %d wwsz %d\n WW vecs: ",wwpt,wwsz);
    printf("wwwid %g wwht %d nsw %g\n WW vecs: ",wwwid,(int)wwht,nsw);
    for (i=0;i<NSW;i++) printf("%d %x %x;",i,ww[i],wwo[i]);
  }}
  ENDVERBATIM
}

:** id() and pid() identify the cell -- printf and function return
FUNCTION pid () {
  VERBATIM 
  printf("INTF%d(%d/%d@%g) ",IDP->id,IDP->type,IDP->col,t);
  _lpid = (double)IDP->id;
  ENDVERBATIM
}

FUNCTION id () {
  VERBATIM
  if (ifarg(1)) IDP->id = (unsigned int) *getarg(1);
  _lid = (double)IDP->id;
  ENDVERBATIM
}

FUNCTION type () {
  VERBATIM
  if (ifarg(1)) IDP->type = (unsigned char) *getarg(1);
  _ltype = (double)IDP->type;
  ENDVERBATIM
}

FUNCTION col () {
  VERBATIM 
  ip = IDP; 
  if (ifarg(1)) ip->col = (unsigned int) *getarg(1);
  _lcol = (double)ip->col;
  ENDVERBATIM
}

FUNCTION dbx () {
  VERBATIM 
  ip = IDP; 
  if (ifarg(1)) ip->dbx = (unsigned char) *getarg(1);
  _ldbx = (double)ip->dbx;
  ENDVERBATIM
}

:** initrec(name,vec) sets up recording of name (see varnum for list) into a vector
PROCEDURE initrec () {
  VERBATIM 
  {int i;
  name = gargstr(1);
  if (SOP==nil) { 
    IDP->record=1;
    SOP = (vpt*)ecalloc(1, sizeof(vpt));
    SOP->size=0;
  }
  if (IDP->record==0) {
    recini();
    IDP->record=1;
  }
  vp = SOP;
  i=(int)varnum();
  if (i==-1) {printf("INTF record ERR %s not recognized\n",name); hoc_execerror("",0); }
  vp->vv[i]=vector_arg(2);
  vector_arg_px(2, &(vp->vvo[i]));
  if (vp->size==0) { vp->size=(unsigned int)vector_buffer_size(vp->vv[i]);
  } else if (vp->size != (unsigned int)vector_buffer_size(vp->vv[i])) {
    printf("INTF initrec ERR vectors not all same size: %d vs %d",vp->size,vector_buffer_size(vp->vv[i]));
    hoc_execerror("", 0); 
  }} 
  ENDVERBATIM
}

:** varnum(statevar_name) returns index number associated with particular variable name
: called by initrec() using global name
FUNCTION varnum () { LOCAL i
  i=-1
  VERBATIM
  if (strcmp(name,"time")==0)      { _li=0.;
  } else if (strcmp(name,"VAM")==0) { _li=1.;
  } else if (strcmp(name,"VNM")==0) { _li=2.;
  } else if (strcmp(name,"VGA")==0) { _li=3.;
  } else if (strcmp(name,"VGB")==0) { _li=4.;
  } else if (strcmp(name,"AHP")==0) { _li=5.;
  } else if (strcmp(name,"V")==0) { _li=6.;
  } else if (strcmp(name,"VM")==0) { _li=6.; // 2 names for V
  } else if (strcmp(name,"VTHC")==0) { _li=7.;
  }
  ENDVERBATIM
  varnum=i
}

:** vecname(INDEX) prints name when given an index
PROCEDURE vecname () {
  VERBATIM
  int i; 
  i = (int)*getarg(1);
  if (i==0)      printf("time\n");
  else if (i==1) printf("VAM\n");
  else if (i==2) printf("VNM\n");
  else if (i==3) printf("VGA\n");
  else if (i==4) printf("VGB\n");
  else if (i==5) printf("AHP\n");
  else if (i==6) printf("V\n");
  else if (i==7) printf("VTHC\n");
  ENDVERBATIM
}

:** initwrec(name,vec) sets up recording of sim field potential
PROCEDURE initwrec () {
  VERBATIM 
  {int i, k, num, cap;  Object* ob;
    ob =   *hoc_objgetarg(1); // list of vectors
    num = ivoc_list_count(ob);
    if (num>NSW) { printf("INTF initwrec() WARN: can only store %d ww vecs\n",NSW); hxe();}
    nsw=(double)num;
    for (k=0;k<num;k++) {
      cap = list_vector_px2(ob, k, &wwo[k], &ww[k]);
      if (k==0) wwsz=cap; else if (wwsz!=cap) {
        printf("INTF initwrec ERR w-vecs size err: %d,%d,%d",k,wwsz,cap); hxe(); }
    }
  }
  ENDVERBATIM
}

: popspk() is paste on gaussian for a pop spk: with vdt=0.1 -20 to 20 is 4 ms
: needs to be above location where is actively accessed
PROCEDURE popspk (x) {
  TABLE Psk DEPEND wwwid,wwht FROM -40 TO 40 WITH 81
  Psk = -wwht*exp(-2.*x*x/wwwid/wwwid)
}

PROCEDURE pskshowtable () {
  VERBATIM 
  int j;
  printf("_tmin_popspk:%g -_tmin_popspk:%g\n",_tmin_popspk,-_tmin_popspk);
  for (j=0;j<=-2*(int)_tmin_popspk+1;j++) printf("%g ",_t_Psk[j]);
  printf("\n");
  ENDVERBATIM 
}

:** wrecord() records voltages onto single global vector
PROCEDURE wrecord (te) {
  VERBATIM 
  {int i,j,k,max,wrp; double ti,scale;
  for (i=0;i<WRNUM && (wrp=(int)IDP->wreci[i])>-1;i++) {
    // wrp: index for multiple field recordings
    scale=(double)IDP->wscale[i];
    if (_lte<1.e9) { // a spike recording
      if (scale>0) {
        max=-(int)_tmin_popspk; // max of table max=-min
        k=(int)floor((_lte-rebeg)/vdt+0.5);
        for (j= -max;j<=max && k+j>0 && k+j<wwsz;j++) {
          wwo[wrp][k+j] += scale*_t_Psk[j+max]; // direct copy from the Psk table
        }
      }
    } else if (twg>=t) { return;
    } else {
      for (ti=twg,k=(int)floor((twg-rebeg)/vdt+0.5);ti<=t && k<wwsz;ti+=vdt,k++) { 
        valps(ti,twg);  // valps() for pop spike calculation
        wwo[wrp][k]+=vii[6];
        if (IDP->dbx==-1) printf("%g:%g ",vii[6],wwo[wrp][k]);
      }
    }
  }
  if (_lte==1.e9) twg=ti;
  }
  ENDVERBATIM
}

: backward compatibility -- note that index was 1-offset; convert to 0 offset here
: wrec() -- return value in wrec0
: wrec(VAL) -- set wrec0
: wrec(VAL,SCALE) -- set wrecIND and scaling for wrecIND
FUNCTION wrec () {
  VERBATIM
  { int k,ix;
  ip=IDP; 
  if (ifarg(1)) {
    ix=(int)*getarg(1);
    if (ix>=1) {
      if (ix-1>=nsw) {
        printf("Attempt to save into ww[%d] but only have %d\n",ix-1,(int)nsw); hxe();}
      ip->wrec=1;
      ip->wreci[0]=(char)ix-1;
      ip->wscale[0]=1.; // default
      if (ifarg(2)) ip->wscale[0]= (float)*getarg(2); 
    } else if (ix<=0) {
      ip->wrec=0;
      for (k=0;k<WRNUM;k++) { ip->wreci[k]=-1; ip->wscale[k]=-1.0; }
    } else {printf("INTF wrec ERR flag(0/1) %d\n",ip->wrec); hxe();
    }
  }
  _lwrec=(double)ip->wrec;
  }
  ENDVERBATIM
}

: wrc() -- return value in wrec0
: wrc(VAL) -- set wrec0
: wrc(IND,VAL) -- set wrec0 and scaling for wrec0
: wrc(IND,VAL,SCALE) -- set wrecIND and scaling for wrecIND
FUNCTION wrc () {
  VERBATIM
  { int i,ix;
  ip=IDP; 
  if (ifarg(1)) {  // 1 or 2 args
    ix=(int)*getarg(1);
    if (ix<0) {
      ip->wrec=0;
      for (i=0;i<WRNUM;i++) { ip->wreci[i]=-1; ip->wscale[i]=-1.0; }
    } else {
      for (i=0;i<WRNUM && ip->wreci[i]!=-1 && ip->wreci[i]!=ix;i++) {};
      if (i==WRNUM) {
        pid(); printf("INFT wrc() ERR: out of wreci pointers (max %d)\n",WRNUM); hxe();}
      if (ix>=nsw) {printf("Attempt to save into ww[%d] but only have %d\n",ix,(int)nsw); hxe();}
      ip->wrec=1; 
      ip->wreci[i]=ix;
      if (ifarg(2)) ip->wscale[i]=(float)*getarg(2); else ip->wscale[i]=1.0;
    }
  } else {
    for (i=0;i<WRNUM;i++) printf("%d:%g ",ip->wreci[i],ip->wscale[i]);
    printf("\n");
  }
  _lwrc=(double)ip->wrec;
  }
  ENDVERBATIM
}

FUNCTION wwszset () {
  VERBATIM
  if (ifarg(1)) wwsz = (unsigned int) *getarg(1);
  _lwwszset=(double)wwsz;
  ENDVERBATIM
}

:** wwfree()
FUNCTION wwfree () {
  VERBATIM
  int k;
  IDP->wrec=0;
  wwsz=0; wwpt=0; nsw=0.;
  for (k=0;k<NSW;k++) { ww[k]=nil; wwo[k]=nil; }
  ENDVERBATIM
}

:** jttr() reads out of a noise vector (call from NET_RECEIVE block)
FUNCTION jttr () {
  VERBATIM 
  if (jtmax>0 && jtpt>=jtmax) {  
    jtpt=0;
    printf("Warning, cycling through jttr vector at t=%g\n",t);
  }
  if (jtmax>0) _ljttr = jsp[jtpt++]; else _ljttr=0;
  ENDVERBATIM
}

:** global_init() initialize globals shared by all INTFs
PROCEDURE global_init () {
  popspk(0) : recreate table if any change in wid or ht
  VERBATIM 
  { int i,j,k; double stt[3];
  if (jridv) { jri=0; vector_resize(jridv, jrmax); vector_resize(jrtvv, jrmax); }
  if (nsw>0. && wwo[0]!=0) { // do just once
    printf("Initializing ww to record for %g (%g)\n",vdt*wwsz,vdt);
    wwpt=0;
    for (k=0;k<(int)nsw;k++) {
      vector_resize(ww[k], wwsz);
      for (j=0;j<wwsz;j++) wwo[k][j]=0.;
    }
  }
  spktot=0;
  jtpt=0;
  eventtot=0;
  errflag=0;
  for (i=0;i<CTYN;i++) blockcnt[cty[i]]=spikes[cty[i]]=0;
  }
  ENDVERBATIM
}

PROCEDURE global_fini () {
  VERBATIM
  int k;
  for (k=0;k<(int)nsw;k++) vector_resize(ww[k], (int)floor(t/vdt+0.5));
  if (jridv && jri<jrmax) {
    vector_resize(jridv, jri); 
    vector_resize(jrtvv, jri);
  }
  ENDVERBATIM
}

:* setting and getting flags: fflag, record,input,jttr
FUNCTION fflag () { fflag=1 }
FUNCTION thrh () { thrh=VTH-RMP }
: reflag() used internally; can't set from here; use recinit()
FUNCTION recflag () { 
  VERBATIM
  _lrecflag= (double)IDP->record;
  ENDVERBATIM
}

: vinflag() used internally; can't set from here; use global initvspks() and range vinset()
FUNCTION vinflag () {
  VERBATIM
  ip=IDP;
  if (ip->vinflg==0 && vsp==nil) { // do nothing
  } else if (ip->vinflg==1 && ip->rve==-1) {
    printf("INTF vinflag ERR: pointer not initialized\n"); hoc_execerror("",0); 
  } else if (ip->rve >= 0) { 
    if (vsp==nil) {
      printf("INTF vinflag ERR1: pointer not initialized\n"); hoc_execerror("",0); 
    }
    ip->rvi=ip->rvb;
    ip->input=1;
  }
  _lvinflag= (double)ip->vinflg;
  ENDVERBATIM
}

:** flag(name,[val,setall]) set or get a flag
:   flag(name,vec) fill vec with flag value from all the cells
: seek names from iflags[] and look at location &ip->type -- beginning of flags
FUNCTION flag () {
  VERBATIM
  char *sf; int ii,i; unsigned char val; double *x; void *voi;
  ip = IDP;
  sf = gargstr(1);
  for (ii=0;ii<iflnum && strncmp(sf, &iflags[ii*4], 3)!=0;ii++) ;
  if (ii==iflnum) {printf("INTF ERR: %s not found as a flag (%s)\n",sf,iflags); hxe();}
  if (ifarg(2)) {
    if (hoc_is_double_arg(2)) { (&ip->type)[ii]=val=(unsigned char)*getarg(2);  
    } else {
      voi=vector_arg(2);  x=vector_newsize(voi,cesz);
      for (i=0;i<cesz;i++) { lop(ce,i); 
        x[i]=((ii>=iflneg)?(double)(char)(&qp->type)[ii]:(double)(&qp->type)[ii]);
      }
    }
  }
  if (ifarg(3)) for (i=0;i<cesz;i++) { 
    lop(ce,i); (&qp->type)[ii]=val; 
  }
  _lflag=(double)(unsigned char)(&ip->type)[ii];
  ENDVERBATIM
}

:** resetall()
PROCEDURE resetall () {
  VERBATIM
  char *sf; int ii,i; unsigned char val;
  for (i=0;i<cesz;i++) { 
    lop(ce,i);
    Vm=RMP; VAM=0; VNM=0; VGA=0; VGB=0; VGBa=0; offsetGB=0; AHP=0; rebob=-1e9; invlt=-1; 
    t0=t; tGB=t; refractory=0;
  }
  ENDVERBATIM
}

:** floc(x,y[,z]) // find a cell by location
FUNCTION floc () {
  VERBATIM
  double x,y,z,r,min,rad, *ix; int ii,i,n,cnt; void* voi;
  cnt=0; n=1000; r=-1;
  ip = IDP;
  x = *getarg(1);
  y = *getarg(2);
  z= ifarg(3)?(*getarg(3)):1e9;
  if (ifarg(5)) {
    r= *getarg(4);
    voi=vector_arg(5);
    ix=vector_newsize(voi,n);
  } 
  for (i=0,min=1e9,ii=-1;i<cesz;i++) { 
    lop(ce,i); 
    rad=(x-xloc)*(x-xloc)+(y-yloc)*(y-yloc)+(z==1e9?0.:((z-zloc)*(z-zloc))); // rad^2
    if (r>0 && rad<r*r) {
      if (cnt>=n) ix=vector_newsize(voi,n*=2);
      ix[cnt]=(double)i;
      cnt++;
    }
    if (rad<min) { min=rad; ii=i; }
  }
  if (r>0) ix=vector_newsize(voi,cnt);
  _lfloc=(double)ii;
  ENDVERBATIM
}

:** invlset([val]) set or get the invl flag
FUNCTION invlset () {
  VERBATIM
  ip=IDP;
  if (ifarg(1)) ip->invl0 = (unsigned char) *getarg(1);
  _linvlset=(double)ip->invl0;
  ENDVERBATIM
}

:** vinset([val]) set or get the input flag (for using shared input from a vector)
FUNCTION vinset () {
  VERBATIM
  ip=IDP;
  if (ifarg(1)) ip->vinflg = (unsigned char) *getarg(1);
  if (ip->vinflg==1) {
    ip->input=1;
    ip->rvi = ip->rvb;
  }
  _lvinset=(double)ip->vinflg;
  ENDVERBATIM
}

:* TABLES
PROCEDURE EXPo (x) {
  TABLE RES FROM -20 TO 0 WITH 5000
  RES = exp(x)
}

FUNCTION EXP (x) {
  EXPo(x)
  EXP = RES
}

FUNCTION esinr (x) {
  ESINo(PI*x/tauGB)
  if        (x<tauGB)   {    esinr= (VGBa-offsetGB)*ESIN +offsetGB
  } else if (x>2*tauGB) {    esinr= 0
  } else {                   esinr= rebound*VGBa*ESIN }
}

PROCEDURE ESINo (x) {
  TABLE ESIN FROM 0 TO 2*PI WITH 3000 : one cycle
  ESIN = sin(x)
}

PROCEDURE rates (vv) {
  TABLE Bb DEPEND mg FROM -100 TO 50 WITH 300
  : from Stevens & Jahr 1990a,b
  Bb = 1 / (1 + exp(0.062 (/mV) * -vv) * (mg / 3.57 (mM)))
}

PROCEDURE coop (x) {
  TABLE Gn DEPEND GPkd FROM 0 TO 10 WITH 100
  : from Destexhe and Sejnowski, PNAS 92:9515 1995
  Gn = (x^4)/(x^4+GPkd) : n=4; kd=100
}

