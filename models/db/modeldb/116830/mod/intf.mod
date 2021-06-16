: $Id: intf.mod,v 1.320 2005/10/05 20:05:22 billl Exp $

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
extern void* vector_arg();
extern FILE* hoc_obj_file_arg(int narg);
static void initmodel();

extern double hoc_epsilon;
#define PI 3.141592653589793115997963468544
#define UINT_MAX	4294967295U  // from /usr/include/limits.h
#define nil 0
#define SOP (((id0*) _p_sop)->vp)
#define IDP (*((id0**) &(_p_sop)))
#define NSW 8  // store output in 8 vecs; switch to 4 if using 4 arg initwrec()
#define WSZ 1000000 // use internal vectors statt external
#define NSV 7  // 6 state variables (+ 1 for time)

typedef struct VPT {
 unsigned int  id;
 unsigned int  size;
 unsigned int  p;
 double* vvo[NSV];
 void*    vv[NSV];
} vpt;

typedef struct ID0 {
 vpt*     vp;
 unsigned int  id;
 unsigned int rvb;
 unsigned int rvi;
 int rve;
 short     type;
 short     record;
 short     wrec;
 short     jitter;
 short     input;
 short     vinflg;
 short     invl0;
} id0;

// globals -- range vars must be malloc'ed in the CONSTRUCTOR
static vpt* vp; // vp and ip are used as temporary pointers
static id0* ip;
static char *name;
static int  errflag;      // turn on after generating an error message
static double *vsp, *wsp, *jsp, *invlp; // accessed by all INTF
static unsigned int jtpt,jitmax;
static double vii[NSV];   // temp storage
static unsigned int wwpt,wwsz,wwaz;   // pointer and size for the shared recording vector
FILE *wf1, *wf2;
void*    ww[NSW];
double* wwo[NSW];
float  wwt[WSZ]; float www[WSZ]; unsigned int wwi[WSZ]; char wws[WSZ];
ENDVERBATIM

:* NEURON, PARAMETER, ASSIGNED blocks
NEURON {
  ARTIFICIAL_CELL INTF
  RANGE VAM, VNM, VGA, VGB, AHP            :::: cell state variables
  RANGE tauAM, tauNM, tauGA, tauGB, tauahp, ahpwt :::: time constants and AHP weight
  RANGE VGBdel,tGB,VGBa,rebound,offsetGB   :::: GABAB and rebound
  RANGE RMP,VTH,Vm,Vblock,refractory       :::: Vblock for depol blockade
  RANGE taum,invl,oinvl,WINV,invlt         :::: interval bursting params
  RANGE t0,tg,tGB,refrac                   :::: t0,tg,tGB save times for analytic calc
  RANGE nbur,tbur,cbur                     :::: burst size, interval and statevar
  POINTER sop                              :::: Structure pointer for other range vars
  GLOBAL AMdec,NMdec,GAdec,GBdec           :::: decrement exc bzw inh activations after a spike
  GLOBAL vdt,next,WEX,mg,RES,ESIN,Bb       :::: table look up values for exp,sin,NMDA-Mg dep
  GLOBAL tauGBGP,wGBGP,GPkd,Gn             :::: GABAB G-protein dependence (cooperativity)
  GLOBAL EAM, ENM, EGA, EGB, spkht         :::: "reverse potential" distance from rest
  GLOBAL prnum                             :::: for debugging moves
}

PARAMETER {
  tauAM = 10 (ms)
  tauNM = 300 (ms)
  tauGA = 10 (ms)
  tauGB = 300 (ms)
  tauGBGP = 50 (ms) : drop off for burst effect on GABAB (G-protein)
  taum =  10 (ms)
  invl =  100 (ms)
  WINV =  0
  wGBGP = 1 (ms) : augmentation of G-protein with a spike
  GPkd  = 100    : maintain between 50 and 500 since table only up to 10 spikes (with wGBGP=1)
  ahpwt = 0
  tauahp= 10 (ms)
  refrac = 5 (ms)
  VTH = -45      : fixed spike threshold
  Vblock = -20   : level of depolarization blockade
  vdt = 0.1      : time step for saving state var
  mg = 1         : for NMDA Mg dep.
  sop=0
  AMdec=1       : default is no fall-off
  NMdec=1
  GAdec=1
  GBdec=1
  nbur=1
  tbur=2
  VGBdel=30
  rebound=0.01 : cannot be set to 0
  offsetGB=0
  RMP=-65
  EAM = 65
  ENM = 90
  EGA = -15
  EGB = -30
  spkht = 50
  prnum = -1
}

ASSIGNED {
  VAM
  VNM
  VGA
  VGB
  VGBa
  AHP
  t0(ms)
  tGB(ms)
  tg(ms)
  refractory
  next
  WEX
  RES
  ESIN
  Gn
  Bb
  cbur
  Vm
  invlt
  oinvl
}

:* CONSTRUCTOR, DESTRUCTOR, INITIAL
: create a structure to save the identity of this unit and short integer flags
CONSTRUCTOR {
  VERBATIM 
  { int lid,lty;
    if (ifarg(2)) { lid=(int) *getarg(2); } else { lid= UINT_MAX; }
    if (ifarg(3)) { lty=(int) *getarg(3); } else { lty= -1; }
    _p_sop = (void*)ecalloc(1, sizeof(id0));
    ip = IDP;
    ip->id=lid; ip->type=lty; 
    ip->invl0 = ip->record = ip->jitter = ip->input = 0; // all flags off
    ip->rve=-1;
  }
  ENDVERBATIM
}

DESTRUCTOR {
  VERBATIM { 
  free(IDP);
  }
  ENDVERBATIM
}

INITIAL {
  VAM = 0
  VNM = 0
  VGA = 0
  VGB = 0
  VGBa = 0
  t0 = t
  tGB = t
  tg = 0
  offsetGB=0
  AHP=0
  invlt = -1
  VERBATIM
  jtpt=0;    // ok to initialize since not altered during init
  errflag=0;
  ENDVERBATIM
  refractory = 0 : 1 means cell is absolute refractory
  : init with vinset(0) if will turn on via a NetCon with w5=1
  if (vinflag()) { randspk() net_send(next,2) }
  if (recflag()) { recini() } : recini() resets for recording, cf recinit()
}

:* NET_RECEIVE
NET_RECEIVE (wAM,wNM,wGA,wGB,wflg) { LOCAL tmp,rflg,wrec,id,jflg,iflg,invlflg,ty
  INITIAL { wNM=wNM wGA=wGA wGB=wGB wflg=0}
  : intra-burst, generate next spike as needed
  VERBATIM
  ip = IDP;  // grab all the flags
  _lrflg=(double)ip->record; _ljflg=(double)ip->jitter; _liflg=(double)ip->input; 
  _linvlflg=(double)ip->invl0; _lwrec=(double)ip->wrec; _lid=(double)ip->id; 
  _lty=(double)ip->type; 
  ENDVERBATIM
  if (flag==4) { : mid-burst
    cbur=cbur-1  : count down the spikes
    if (cbur>0) { net_send(tbur,4) 
    } else { net_send(refrac-AHP/10, 3) }
    if (jflg) { tmp= t+jitter()/10 } else { tmp=t }
    net_event(tmp)
    if (rflg) { recspk(tmp) }
    if (wrec) { wrecord(tmp,-1,0) }
  : start reading random spike times (or burst times) from vsp vector pointer
  : this is signaled externally from a netstim with wflg=1, will turn off on next stim 
  : (NB wflg used in completely different context for GABAB) 
  } else if (flag==0 && wGB==0 && wflg==1) {
    VERBATIM
    ip->input=1;
    ENDVERBATIM
    iflg=1 : two versions of the same flag
    wflg=2 : set flag to turn off next time an external event comes from here
    randspk() 
    net_send(next,2)
  } else if (flag==0 && wGB==0 && wflg==2) { : flag to stop random spikes
    VERBATIM
    ip->input=0;
    ENDVERBATIM
    iflg=0 : two versions of the same flag
    wflg=1  : flag to turn on next time
  } else {  : external input
    if (rflg) { record() }
    : update all statevars
    if (VAM>hoc_epsilon)  { VAM = VAM*EXP(-(t - t0)/tauAM) } else { VAM=0 } :AMPA
    if (VNM>hoc_epsilon)  { VNM = VNM*EXP(-(t - t0)/tauNM) } else { VNM=0 } :NMDA
    if (VGA< -hoc_epsilon){ VGA = VGA*EXP(-(t - t0)/tauGA) } else { VGA=0 } :GABAA    
    VGB = esinr(t-tGB) : VGB has to update each t but calc based on triggering tGB and val VGBa
    if (AHP< -hoc_epsilon){ AHP = AHP*EXP(-(t-t0)/tauahp) } else { AHP=0 } : adaptation
    : for debugging if (VGA<EGA) { pid() printf("CC: %g %g %g %g %g\n",t,VGA,wGA,EGA,Vm) }
    Vm = VAM+VNM+VGA+VGB+AHP : membrane deviation from rest
    if (Vm>100||Vm<-60){ pid() 
      printf("WARN:t=%g Vm=%g (%g,%g,%g,%g,%g)\n",t,Vm,VAM,VNM,VGA,VGB,AHP) }
    if (flag==0) { : only add weights if an external excitation
      : AMPA Erev=0 (0-RMP==65 mV above rest)
      if (wAM>0 && VAM<EAM) { 
        tmp = wAM*(1-Vm/EAM)
        VAM = VAM + tmp
        if (wrec) { wrecord(t,0,tmp) }
      }
      : NMDA; Mg effect based on total activation in rates()
      if (wNM>0 && VNM<ENM) { rates(RMP+Vm)
        tmp = wNM*Bb*(1-Vm/ENM) 
        VNM = VNM + tmp
        if (wrec) { wrecord(t,1,tmp) }
      } 
      if (VNM>1.2*ENM) { pid() : signal if some nasty number creeps in here
        : allow it to creep over by a little which can happend with coincident spikes
        printf("**** ERR: t=%g VNM=%g wNM=%g ENM=%g Vm=%g\n",t,VNM,wNM,ENM,Vm) 
      }
      : GABAA and GABAB: note that all wts are positive
      if (wGA>0 && VGA>EGA) { 
        tmp = wGA*(1-Vm/EGA) 
        VGA = VGA - tmp
        if (wrec) { wrecord(t,2,tmp) }
      }
      if (wGB>0) { net_send(VGBdel,5) } : delayed effect
      if (invlflg) { : repetitive activity weight
        Vm = RMP+VAM+VNM+VGA+VGB+AHP
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
    } else if (flag==5) { : flag==5 to set GABAB weight after a delay
      offsetGB = VGB : current position
      : wflg overloaded; 50 is tauGBGP for GB cooperativity
      : wflg will augment each time there is spike coming in through this line
      wflg=wflg*EXP(-(t-tGB)/tauGBGP)+wGBGP 
      coop(wflg)               : cooperativity -- need mult presyn spikes to activate
      : calculate separately based on VGBa and tGB
      if (VGB>EGB) { 
        tmp = wGB*(1-Vm/EGB)*Gn 
        VGB = VGB - tmp
        if (wrec) { wrecord(t,3,tmp) }
      }
      VGBa= VGB
      tGB=t : restart for VGB
    : flag==2 -- read off external vec (vsp) for next random spike time
    } else if (flag==2) { 
      if (iflg==0) { flag=-1 : turned off in meantime so don't spike below
      } else { 
        randspk() 
        net_send(next,2) 
        if (WEX<0) { 
          net_event(t)   : bypass activation calculation
          if (rflg) { recspk(t) }
          if (wrec) { wrecord(t,-1,0) }
        } else {
          tmp = WEX*(1-Vm/EAM)
          VAM = VAM + tmp
          if (wrec) { wrecord(t,0,tmp) }
        }
      }
    } else if (flag==1) { 
      : Vm=RMP+VAM+VNM+VGA+VGB+AHP
      if (WINV<0) { 
        net_event(t)   : bypass activation calculation
        if (rflg) { recspk(t) }
        if (wrec) { wrecord(t,-1,0) }
      } else {
        tmp = WINV*(1-Vm/EAM)
        VAM = VAM + tmp :: activate interval depolarization
        if (wrec) { wrecord(t,0,tmp) }
      }
      oinvl=invl
      invlt=t
      net_send(invl,1) 
    } else if (flag==3) { 
      refractory = 0 :end of refractory period
    }
    : 3 causes of spiking: between VTH and Vblock, random from vsp (flag 2), within burst (v.s.)
    Vm = VAM+VNM+VGA+VGB+RMP+AHP
    if (refractory==0 && (Vm>VTH && Vm<Vblock)) {
      AHP = AHP - ahpwt
      if (jflg) { tmp=t+jitter() } else { tmp=t }
      net_event(tmp)
      if (rflg) { recspk(tmp) }
      if (wrec) { wrecord(tmp,-2,0) }
      VAM=VAM*AMdec VNM=VNM*NMdec
      VGA=VGA*GAdec VGB=VGB*GBdec
      if (nbur>1) { cbur=nbur-1 net_send(tbur,4) 
      } else { net_send(refrac-AHP/10, 3) } : AHP doubles as a refrac period extender
      refractory = 1
    }
    t0 = t
  }
}

:* ancillary functions
:** randspk() sets next to next val in vector, this vector is handled globally
PROCEDURE randspk () {
  VERBATIM 
  ip=IDP;  
  if (ip->rvi > ip->rve) { ip->input=0;
  } else { 
    WEX=wsp[ip->rvi];
    next=vsp[ip->rvi++]-t; // vector contains absolute times; adjust to make interval
  }
  ENDVERBATIM
  : net_send(next,2) : must be called from appropriate blocks
}

:** val(t,tstart) fills global vii[] to pass values back to record() (called from record())
VERBATIM
double val (double xx, double ta) { 
  vii[1]=VAM*EXP(-(xx - ta)/tauAM);
  vii[2]=VNM*EXP(-(xx - ta)/tauNM);
  vii[3]=VGA*EXP(-(xx - ta)/tauGA);
  vii[4]=esinr(xx-tGB);
  vii[5]=AHP*EXP(-(xx - ta)/tauahp);
  vii[6]=vii[1]+vii[2]+vii[3]+vii[4]+vii[5];
}
ENDVERBATIM

:** record() stores values since last tg into appropriate vecs
PROCEDURE record () {
  VERBATIM {
  int k; double ti;
  vp = SOP;
  if (tg>=t) return;
  if (vp->p >= vp->size) { if (errflag) return; 
    printf("**** WARNING out of recording room for INTF type%d id%d at %g****\n",IDP->type,IDP->id,t);
    printf("**************** WARNING: No further WARNINGS ****************\n");
    errflag=1; return; }
  for (ti=tg;ti<=t && vp->p < vp->size;ti+=vdt,vp->p++) { 
    val(ti,tg);  
    vp->vvo[0][vp->p]=ti;
    for (k=1;k<NSV;k++) if (vp->vvo[k]!=0) { // not nil pointer
      vp->vvo[k][vp->p]=vii[k]+RMP;
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
  if (vp->p >= vp->size || vp->vvo[6]==0) return; 
  vp->vvo[0][vp->p-1]=_lx;
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
: this is a range procedure
: all cells share one vector but read from different locations
: (CHANGED from intervals and global proc in v224)
PROCEDURE initvspks () {
  VERBATIM 
  {int max, i, err=0; 
    double *iv, last;
    if (! ifarg(1)) {printf("Return initvspks(ivspks,vspks,wvspks)\n"); return 0.;}
    ip=IDP;
    i = vector_arg_px(1, &iv);
    max=vector_arg_px(2, &vsp);
    if (max!=i) {err=1; printf("initvspks ERR: vecs of different size\n");}
    if (max==0) {err=1; printf("initvspks ERR: vec not initialized\n");}
    max=vector_arg_px(3, &wsp);
    if (max!=i) {err=1; printf("initvspks ERR: 3rd vec is of different size\n");}
    ip->vinflg=1;
    for (i=0; i<max && (int)iv[i] != ip->id ; i++); // move forward to first
    if (i==max) { 
      printf("initvspks WARN: %d not found in ivspks\n",ip->id); 
      ip->vinflg=0; ip->rve=-1;
      return(0.); 
    }
    ip->rvb=ip->rvi=i;
    last=vsp[i++];
    for (; i<max && (int)iv[i] == ip->id ; i++) { // move forward to last
      if (vsp[i]<=last) { err=1; 
        printf("initvspks ERR: nonmonotonic for cell#%d: %g %g\n",ip->id,last,vsp[i]); }
      last=vsp[i];
    }
    ip->rve=i-1;
    if (err) { ip->rve=0; hoc_execerror("",0); }
  }
  ENDVERBATIM
}

:** initjitter() sets up vector from which to read jitter
: this is a global not a range procedure -- just call once
PROCEDURE initjitter () {
  VERBATIM 
  {int max, i, err=0;
    jtpt=0;
    if (! ifarg(1)) {printf("Return initjitter(vec)\n"); return(0.);}
    max=vector_arg_px(1, &jsp);
    if (max==0) {err=1; printf("initjitter ERR: vec not initialized\n");}
    for (i=0; i<max; i++) if (jsp[i]<=0) {err=1;
      printf("initjitter ERR: vec should be >0: %g\n",jsp[i]);}
    if (err) { jsp=nil; jitmax=0.; return(0.); }// hoc_execerror("",0);
    if (max != jitmax) {
      printf("WARNING: resetting jitmax_INTF to %d\n",max); jitmax=max; }
  }
  ENDVERBATIM
}

:* initinvl() sets up vector from which to read intervals
: this is a global not a range procedure -- just call once
PROCEDURE initinvl () {
  VERBATIM 
  printf("initinvl() NOT BEING USED\n"); return(0.);
  ENDVERBATIM
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
    printf("INTF record ERR: but pointer not initialized\n"); hoc_execerror("",0); 
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
PROCEDURE chk () {
  VERBATIM 
  {int i,lfg;
  ip=IDP;
  printf("ID:%d; typ: %d; rec:%d wrec:%d inp:%d jit:%d invl:%d\n",ip->id,ip->type,ip->record,ip->wrec,ip->input,ip->jitter,ip->invl0);
  if (ifarg(1)) lfg=(int) *getarg(1); else lfg=0;
  if (lfg==1) {
    if (SOP!=nil) {
      vp = SOP;
      printf("p %d size %d tg %g\n",vp->p,vp->size,tg);
      for (i=0;i<NSV;i++) printf("%d %x %x;",i,vp->vv[i],vp->vvo[i]);
    } else printf("Recording pointers not initialized");
  }
  if (lfg==2) { 
    printf("Global vectors for input and jitter: \n");
    if (vsp!=nil) printf("VSP: %x (%d/%d-%d)\n",vsp,ip->rvi,ip->rvb,ip->rve); else printf("no VSP\n");
    if (jsp!=nil) printf("JSP: %x (%d/%d)\n",jsp,jtpt,jitmax); else printf("no JSP\n");
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
    for (i=0;i<NSW;i++) printf("%d %x %x;",i,ww[i],wwo[i]);
  }}
  ENDVERBATIM
}

:** id() and pid() identify the cell -- printf and function return
FUNCTION pid () {
  VERBATIM 
  printf("INTF%d(%d) ",IDP->id,IDP->type);
  _lpid = (double)IDP->id;
  ENDVERBATIM
}

FUNCTION id () {
  VERBATIM 
  _lid = (double)IDP->id;
  ENDVERBATIM
}

FUNCTION type () {
  VERBATIM 
  _ltype = (double)IDP->type;
  ENDVERBATIM
}

:** initrec(name,vec) sets up recording of name (see varnum for list) into a vector
PROCEDURE initrec () {
  VERBATIM 
  {int i; void *vv;
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
  ENDVERBATIM
}

:* rebuild() -- build the vvo vectors from stored wwo information
PROCEDURE rebuild () { LOCAL s0,w0,wwaz,ii,wflg,tmp
  VERBATIM
  int ii,jj; double i0,tstop;
  ip = IDP; vp=SOP; // grab all the flags
  ip->record=1; ip->input=0; i0=wwo[1][0];
  initmodel(); _lwwaz=wwaz; _lii=0.;
  tstop=(ifarg(1))?(*getarg(1)):0.;
  ENDVERBATIM
  WHILE (ii<wwaz) { : no 'for' loop in NMODL
    VERBATIM
    int ii=(int)_lii;
    if (wwo[1][ii]!=i0) {
      printf("ERROR wrong id at %d %g not %g\n",ii,wwo[1][ii],i0); hoc_execerror("", 0);}
    t=wwo[0][ii]; _ls0=wwo[2][ii]; _lw0=wwo[3][ii];
    ENDVERBATIM
    record()
    : update all statevars
    if (VAM>hoc_epsilon)  { VAM = VAM*EXP(-(t - t0)/tauAM) } else { VAM=0 } :AMPA
    if (VNM>hoc_epsilon)  { VNM = VNM*EXP(-(t - t0)/tauNM) } else { VNM=0 } :NMDA
    if (VGA< -hoc_epsilon){ VGA = VGA*EXP(-(t - t0)/tauGA) } else { VGA=0 } :GABAA    
    VGB = esinr(t-tGB) : VGB has to update each t but calc based on triggering tGB and val VGBa
    if (AHP< -hoc_epsilon){ AHP = AHP*EXP(-(t-t0)/tauahp) } else { AHP=0 } : adaptation
    if (s0==0) { VAM = VAM + w0 }
    if (s0==1) { VNM = VNM + w0 }
    if (s0==2) { VGA = VGA - w0 }
    if (s0==3) { 
      offsetGB = VGB : current position
      VGB = VGB - w0
      VGBa= VGB
      tGB=t
    }
    if (s0==-2) { 
      recspk(t) 
      AHP = AHP - ahpwt
      VAM=VAM*AMdec VNM=VNM*NMdec
      VGA=VGA*GAdec VGB=VGB*GBdec
    }
    if (s0==-1) { recspk(t) }
    t0 = t
    ii = ii+1
  }
  if (tstop>0) {
    VERBATIM
    t=tstop;  // not allowed to set t in a mod file
    ENDVERBATIM
    record()
  }
}

:* wrec block
:** initwrec(vec1,vec2,vec3,vec4) sets up recording of external events
PROCEDURE initwrec () {
  VERBATIM 
  {int k;
  if (! ifarg(4)) { // assign 2 files
    wwsz=WSZ;
    wf1 = hoc_obj_file_arg(1);
    wf2 = hoc_obj_file_arg(2);
  } else if (! ifarg(8)) { // assign 4 vectors
    if (NSW!=4) { 
      printf("INTF initwrec ERR w-vecs compiled for 4 args\n");
      hoc_execerror("",0); }
    IDP->wrec=1;
    for (k=0;k<NSW;k++) {
      ww[k]=vector_arg(k+1);
      wwaz=vector_arg_px(k+1, &(wwo[k]));
    }
    if (wwsz==0) wwsz=(unsigned int)vector_buffer_size(ww[0]); 
    for (k=0;k<NSW;k++) if (wwsz!=(unsigned int)vector_buffer_size(ww[k])) {
      printf("INTF initwrec ERR w-vecs size err: %d,%d,%d",k,wwsz,vector_buffer_size(ww[k]));
    }
  } else { // assign 8 vectors
    if (NSW!=8) { 
      printf("INTF initwrec ERR w-vecs compiled for 8 args\n");
      hoc_execerror("",0); }
    IDP->wrec=1;
    for (k=0;k<NSW;k++) {
      ww[k]=vector_arg(k+1);
      wwaz=vector_arg_px(k+1, &(wwo[k]));
    }
    if (wwsz==0) wwsz=(unsigned int)vector_buffer_size(ww[0]); 
    for (k=0;k<NSW;k++) if (wwsz!=(unsigned int)vector_buffer_size(ww[k])) {
      printf("INTF initwrec ERR w-vecs size err: %d,%d,%d",k,wwsz,vector_buffer_size(ww[k]));
    }
  }}
  ENDVERBATIM
}

:** wrecord()
PROCEDURE wrecord (t,s0,w0) {
  VERBATIM {
  int k; double id = (double)IDP->id;
  if (wwpt >= wwsz) { 
    wwpt=0;
    fprintf(wf1,"//b8 %d INTF %g %d\n",WSZ,_lt,ftell(wf2));
    fwrite(&wwt,sizeof(float),WSZ,wf2);  // write out the size
    fwrite(&wwi,sizeof(int),WSZ,wf2);  // write out the size
    fwrite(&wws,sizeof(char),WSZ,wf2);  // write out the size
    fwrite(&www,sizeof(float),WSZ,wf2);  // write out the size
  } 
  // wwo[0][wwpt]=_lt; wwo[1][wwpt]=id; wwo[2][wwpt]=_ls0; wwo[3][wwpt]=_lw0; wwpt++;
  wwt[wwpt]=(float)_lt; 
  wwi[wwpt]=(unsigned int)IDP->id; 
  wws[wwpt]=(char)_ls0; 
  www[wwpt]=(float)_lw0; 
  wwpt++;
  }
  ENDVERBATIM
}

FUNCTION wrec () {
  VERBATIM
  ip=IDP;
  if (ifarg(1)) ip->wrec = (short) *getarg(1);
  _lwrec=(double)ip->wrec;
  ENDVERBATIM
}

FUNCTION wwszset () {
  VERBATIM
  if (ifarg(1)) wwsz = (short) *getarg(1);
  _lwwszset=(double)wwsz;
  ENDVERBATIM
}

:** wwfree()
FUNCTION wwfree () {
  VERBATIM
  int k;
  IDP->wrec=0;
  wwsz=0; wwpt=0;
  for (k=0;k<NSW;k++) { ww[k]=nil; wwo[k]=nil; }
  ENDVERBATIM
}

:* jitter
:** jitter() reads out of a noise vector (call from NET_RECEIVE block)
FUNCTION jitter () {
  if (jitmax>0 && jtpt>=jitmax) {  jtpt=0
    printf("Warning, cycling through jitter vector at t=%g\n",t) }
  if (jitmax>0) {
    VERBATIM 
    _ljitter = jsp[jtpt++];
    ENDVERBATIM
  } else { jitter=0 }
}

:** initialize globals shared by all INTFs
PROCEDURE global_init () {
  VERBATIM
  int k;
  if (wwo[0]!=0) { // do just once
    printf("Initializing global ww vectors\n");
    wwpt=0;
    for (k=0;k<NSW;k++) vector_resize(ww[k], wwsz);
  }
  ENDVERBATIM
}

PROCEDURE global_fini () {
  VERBATIM
  {int k;
  if (IDP->wrec) {
    if (wwo[0]!=0) { 
      for (k=0;k<NSW;k++) vector_resize(ww[k], wwpt);
    } else {
      fprintf(wf1,"//b8 %d INTF %g %d\n",wwpt,t,ftell(wf2));
      fwrite(&wwt,sizeof(float),wwpt,wf2);  // write out the size
      fwrite(&wwi,sizeof(int),wwpt,wf2);  // write out the size
      fwrite(&wws,sizeof(char),wwpt,wf2);  // write out the size
      fwrite(&www,sizeof(float),wwpt,wf2);  // write out the size
      printf("Closing file with wwpt=%d at location %d\n",wwpt,ftell(wf2));
      fclose(wf1); fclose(wf2);
    }
  } else {
  printf("WARNING: global_fini() called from %d:%d with no wrec pointers\n",IDP->type,IDP->id);
  }}
  ENDVERBATIM
}

:** setting and getting flags: fflag, record,input,jitter
FUNCTION fflag () { fflag=1 }
FUNCTION thresh () { thresh=VTH-RMP }

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

: jitset([val]) set or get the jitter flag
FUNCTION jitset () {
  VERBATIM
  ip=IDP;
  if (ifarg(1)) ip->jitter = (short) *getarg(1);
  _ljitset=(double)ip->jitter;
  ENDVERBATIM
}

: invlset([val]) set or get the invl flag
FUNCTION invlset () {
  VERBATIM
  ip=IDP;
  if (ifarg(1)) ip->invl0 = (short) *getarg(1);
  _linvlset=(double)ip->invl0;
  ENDVERBATIM
}

: vinset([val]) set or get the input flag (for using shared input from a vector)
FUNCTION vinset () {
  VERBATIM
  ip=IDP;
  if (ifarg(1)) ip->vinflg = (short) *getarg(1);
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

PROCEDURE rates(vv) {
  TABLE Bb DEPEND mg FROM -100 TO 50 WITH 300
  : from Stevens & Jahr 1990a,b
  Bb = 1 / (1 + exp(0.062 (/mV) * -vv) * (mg / 3.57 (mM)))
}

PROCEDURE coop (x) {
  TABLE Gn DEPEND GPkd FROM 0 TO 10 WITH 100
  : from Destexhe and Sejnowski, PNAS 92:9515 1995
  Gn = (x^4)/(x^4+GPkd) : n=4; kd=100
}
