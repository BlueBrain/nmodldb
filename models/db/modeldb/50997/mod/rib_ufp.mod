COMMENT
Ribbon Synapse with function calls within BREAKPOINT.
Alpha function evaluation for glutamate profile
ENDCOMMENT

DEFINE NumSites 20 : The number of Release Sites per Synapse

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
     POINT_PROCESS ribbon_ca
     EXTERNAL mod_modulator
     POINTER preCA1, preCA2
     RANGE preCAnow
     RANGE CAthresh
     RANGE S, W, K1, powr
     RANGE  Alpha_Max, Alpha_Delay, Alpha_tau, DC_Level, DC_Delay, DC_Off, Ramp_Max
     RANGE  Ramp_Delay, Ramp_Off, Slope_UP, Slope_DOWN
     RANGE R, RM, RA, RA2, RdA, RdA2, RMA, RMA2, RdMA, RdMA2, O, OM
     RANGE g, gmax, erev, rb, rBb, rDb, rDBb, rMb, rMBb, rm, rBm, rB2m, rDBm, rDB2m
     RANGE Rb, RBu, RBd, RBDr, RBb, RB2u
     RANGE RB2d, RB2Dr, RB2o, RB2c, RDBb, RDB2u
     RANGE Rm, RMum, RBm, RBMum, RB2m, RB2Mum
     RANGE RMb, RMBu, RMBb, RMB2u, RMB2o, RMB2c
     RANGE RMB2d, RMB2Dr, RMBd, RMBDr, RMDBb, RMDB2u
     RANGE RDBMum, RDBm, RDB2Mum, RDB2m
     RANGE tau_dAMPA, glumax_dAMPA, e, AbsRefract,  gluConc_dAMPA, pr, StanddAMPA, rate_constant
:  for the NMDA receptor
     RANGE g2, gmaxN, erevN, rb2, mg
     RANGE C0, C1, C2, D, OO, B
     RANGE Rb2, Ru, Rd, Rr, Ro, Rc, alpend
     RANGE tau_NMDA, gluConc_NMDA, glumax_NMDA, StandNMDA
     RANGE Max_RVP,  rate_constantIN
     RANGE Vmax, k, n
     RANGE UFP
     GLOBAL total_count, Orate

: for setting state variable at HOC level
:NMDA
     RANGE xC0, xC1, xC2, xD, xOO
:
     GLOBAL vmin, vmax
     NONSPECIFIC_CURRENT i
}

UNITS {
     (nA) = (nanoamp)
     (mV) = (millivolt)
     (pS) = (picosiemens)
     (umho) = (micromho)
     (mM) = (milli/liter)
     (uM) = (micro/liter)
}

PARAMETER {
     S = 2.84985 : Proportionality Constant rate equation
     W = 0.0498303 : Constant factor in rate equation
     powr = 4 : Power of rate equation
     K1 = 18.3892 : CaX dissociation constant
     CAthresh = .00010001 (mM) : [Ca++] threshold
     preCA1 (mM) : presynaptic [Ca++] for a full RVP
     preCA2 (mM) : presynaptic [Ca++] for a less than full RVP
     preCAnow
     erev = 0 (mV)	: reversal potential of dAMPA
     erevN = 0 (mV)     : reversal potential of NMDA
     gmax = 0 (umho)	: maximal conductance dAMPA
     gmaxN = 0 (umho)   : max conductance NMDA
     StanddAMPA = 0 (mM) : standing level of glutamate seen by dAMPA
     StandNMDA = 0.001 (mM)  : standing level of glutamate seen by NMDA
	xC0 = 0.55004
	xC1 = 0.21319
	xC2 = 0.08263
	xD =  0.10207
	xOO = 0.05207
     mg = 0 (mM)
        AbsRefract = 0 (ms)
        rate_constant (1/s)
        tau_dAMPA=0.055 (ms) 
        tau_NMDA=0.1 (ms)
        glumax_dAMPA=0   (mM) 
        glumax_NMDA=0    (mM)
        e=0     (mV)
        Alpha_Max=0     (umho)
        Alpha_Delay=0 (ms)
        Alpha_tau=.1 (ms)
        DC_Level (mM)
        DC_Delay (ms)
        DC_Off (ms)
        Ramp_Max (mM)
        Ramp_Delay (ms)
        Ramp_Off (ms)
        Slope_UP (ms)
        Slope_DOWN (ms)
        alpend=0
        Max_RVP = 5 : Maximum number of vesicles in RVP for one release site
        rate_constantIN = 0.125 (1/s) : Rate constant of RVP uptake
	  Vmax = 1840.22 : for Hill equation
        k = 86.48823 :        "
        n = 3.30393  :        "

: Rates

     : Partin, Fleck & Mayer (1996) J Nsci 16, 6634-47
     Rb  = 2e-2     (/uM /ms) :binding single ligand
     RBu  = 3e-1    (/ms)     :unbinding from single bound
     RBd  = 1e0     (/ms)     :desensitization from single bound
     RBDr  = 3e-1   (/ms)     :resensitization from single bound
     RBb  =1e-2     (/uM /ms) :binding second ligand
     RB2u  = 1e2    (/ms)     :unbinding from double to single ligand
     RB2d = 8e0     (/ms)     :desensitization from double bound
     RB2Dr = 2e-4   (/ms)     :resensitization from double bound
     RB2o = 3e1     (/ms)     :opening from double bound ligand
     RB2c = 1.5e0   (/ms)     :closing to double bound ligand
     RDBb = 1e-2    (/uM /ms) :desensitized/single bound ligand binding 2nd ligand
     RDB2u = 8.3e-3      (/ms)     :desensitized/double bound ligand unbinding 2nd ligand

: dAMPA modulator kinetics OFF
     Rm  = 0 (/uM /ms)       :binding modulator
     RMum  = 0        (/ms)   :unbinding modulator
     RBm = 0         (/uM /ms)       :single bound ligand binding modulator
     RBMum =0         (/ms)   :single bound ligand unbinding modulator
     RB2m = 0        (/uM /ms)       :double bound ligand binding modulator
     RB2Mum = 0       (/ms)   :double bound ligand unbinding modulator
     RMb = 0  (/uM /ms)      :modulator bound binding single ligand
     RMBu = 0        (/ms)   :modulator bound unbinding single ligand
     RMBb = 0        (/uM /ms)       :modulator bound binding double ligand
     RMB2u = 0        (/ms)   :modulator bound unbinding double ligand
     RMB2o = 0        (/ms)   :opening from modulator plus double bound ligand
     RMB2c = 0     (/ms)   :closing from modulator plus double bound ligand
     RMB2d = 0        (/ms)   :desensitization from double ligand bound with modulator
     RMB2Dr = 0      (/ms)   :resensitization from double ligand bound with modulator
     RMBd = 0 (/ms)   :desensitization from single ligand bound with modulator
     RMBDr = 0       (/ms)   :resensitization to single ligand bound with modulator
     RMDBb = 0       (/uM /ms)       :desensitized/single bound ligand&modulator binding 2nd ligand
     RMDB2u = 0 (/ms)      :desensitized/double bound ligand&modulator unbinding ligand
     RDBm = 0        (/uM /ms)       :desensitized/single bound ligand binding modulator
     RDBMum = 0      (/ms)   :desensitized/single bound ligand unbinding modulator
     RDB2m = 0       (/uM /ms)       :desensitized/double bound ligand binding modulator
     RDB2Mum = 0 (/ms)       :desensitized/double bound ligand unbinding modulator

        : Destexhe, Mainen & Sejnowski, 1996
        Rb2     = 5e-3    (/uM /ms)     : binding
        Ru      = 12.9e-3  (/ms)        : unbinding
        Rd      = 8.4e-3   (/ms)        : desensitization
        Rr      = 6.8e-3   (/ms)        : resensitization
        Ro      = 46.5e-3   (/ms)       : opening
        Rc      = 73.8e-3   (/ms)       : closing

        vmin = -200     (mV)
        vmax = 100      (mV)

}

COMMENT
:    Partin values for modulator
     Rm  = 1e-2 (/uM /ms)       :binding modulator
     RMum  = 5e0        (/ms)   :unbinding modulator
     RBm = 1e-2         (/uM /ms)       :single bound ligand binding modulator
     RBMum =5e0         (/ms)   :single bound ligand unbinding modulator
     RB2m = 1e-2        (/uM /ms)       :double bound ligand binding modulator
     RB2Mum = 5e0       (/ms)   :double bound ligand unbinding modulator
     RMb = 2e-2  (/uM /ms)      :modulator bound binding single ligand
     RMBu = 3e-1        (/ms)   :modulator bound unbinding single ligand
     RMBb = 1e-2        (/uM /ms)       :modulator bound binding double ligand
     RMB2u = 1e2        (/ms)   :modulator bound unbinding double ligand
     RMB2o = 3e1        (/ms)   :opening from modulator plus double bound ligand
     RMB2c = 2.5e-1     (/ms)   :closing from modulator plus double bound ligand
     RMB2d = 8e0        (/ms)   :desensitization from double ligand bound with modulator
     RMB2Dr = 2e-4      (/ms)   :resensitization from double ligand bound with modulator
     RMBd = 1e0 (/ms)   :desensitization from single ligand bound with modulator
     RMBDr = 3e-1       (/ms)   :resensitization to single ligand bound with modulator
     RMDBb = 1e-2       (/uM /ms)       :desensitized/single bound ligand&modulator binding 2nd ligand
     RMDB2u = 8.3e-3 (/ms)      :desensitized/double bound ligand&modulator unbinding ligand
     RDBm = 1e-4        (/uM /ms)       :desensitized/single bound ligand binding modulator
     RDBMum = 5e-2      (/ms)   :desensitized/single bound ligand unbinding modulator
     RDB2m = 1e-4       (/uM /ms)       :desensitized/double bound ligand binding modulator
     RDB2Mum = 5e-2 (/ms)       :desensitized/double bound ligand unbinding modulator
ENDCOMMENT


ASSIGNED {
     gluConc_dAMPA    (mM)
     gluConc_NMDA     (mM)
     v		(mV)	: postsynaptic voltage
     i		(nA)	: current = g*(v - Erev)
     g		(umho)	: conductance dAMPA
     g2         (umho)  : conductance NMDA
     mod_modulator (mM):pointer to concentration of ext. modulator compound (ie cyclothiazide)
     rb		(/ms)	:binding single ligand
     rb2        (/ms)
     rBb		(/ms)	:binding second ligand
     rDBb	(/ms)	:binding second ligand from desensitized state
     rMb		(/ms)	:binding ligand from modulator bound state
     rMBb	(/ms)	:binding second ligand from modulator single ligand bound state
     rm		(/ms)	:binding modulator no ligand bound
     rBm		(/ms)	:binding modulator from single ligand bound state
     rB2m	(/ms)	:binding modulator from double ligand bound state
     rDBm	(/ms)	:binding modulator from single ligand bound desensitized state
     rDB2m	(/ms)	:binding modulator from double ligand bound desensitized state
     dt         (ms)
     last_quanta[NumSites] (ms)   :time since last quantal release
     release_start[NumSites] (ms) :time of last quantal release
     RVP_Size[NumSites] : Current size of RVP
     RVP_out[NumSites] : Number of vesicles released at the current time step
     total_count : count of all vesicular releases
     Orate      (/s)  : overall calculated rate for entire ribbon
     rate[NumSites] :
     rateIN[NumSites] :
     UFP : Current size of ultra-fast pool
}

STATE {
     : Channel states (all fractions) dAMPA
     R		: unbound
     RM		: modulator bound no ligand
     RA		: single ligand bound
     RA2		: double ligand bound
     RdA		: desensitized single bound
     RdA2		: desensitized double bound
     RMA		: modulator bound single ligand bound
     RMA2		: modulator bound double ligand bound
     RdMA		: modulator bound single ligand desensitized
     RdMA2		: modulator bound double ligand desensitized
     O		: open, ligand bound
     OM		: open, modulator  & ligand bound

     : Channel states (all fractions) NMDA
        C0              : unbound
        C1              : single bound
        C2              : double bound
        D               : desensitized
        OO              : open

        B               : fraction free of Mg2+ block
}


INITIAL {
     R = 1
FROM i=0 TO NumSites-1{     
           last_quanta[i] = 40
           release_start[i] = 10000 :initialize to something that will always be larger than any t
RVP_Size[i] = Max_RVP
     }
        rates(v)

        C0 = xC0 
        C1 = xC1
        C2 = xC2
        D = xD
        OO = xOO

total_count = 0
UFP = 0
}

BREAKPOINT {
rates(v)

gluConc_dAMPA = StanddAMPA
gluConc_NMDA = StandNMDA

SOLVE vesicle_release METHOD after_cvode

SOLVE RVPSIZE METHOD after_cvode

SOLVE kstates METHOD sparse

     g = gmax *(O+OM)
     g2 = gmaxN * OO * B
     i = (g * (v - erev) ) + (g2 * (v - erevN) )
}

PROCEDURE vesicle_release(){
Orate = 0
: Calculate size of ultra-fast pool
UFP = 0
FROM i=0 TO NumSites-1
{ 
if (RVP_Size[i] == Max_RVP)
   {UFP = UFP + 1}
}
:
FROM i=0 TO NumSites-1
  {
if (RVP_Size[i] < Max_RVP) {preCAnow = preCA2}
else
{preCAnow = preCA1}
:
if (RVP_Size[i] < 1)
{ rate[i] = S * (W * preCAnow * 1000) / ( ((preCAnow * 1000) / K1) + 1)^powr
}  : factor of 1000 for conversion from NEURON mM to formula uM
else
{
rate_constant =  (Vmax * (1000 * preCAnow)^n) / (k^n + (1000 * preCAnow)^n)

if (RVP_Size[i] == Max_RVP)
 { rate[i] = UFP * rate_constant }
else
 { rate[i] =  RVP_Size[i] * rate_constant }
}
Orate = Orate + rate[i]
:
  last_quanta[i] = last_quanta[i] + dt
: The probability of release is linearly related to rate (rate/(1000/dt)
  if ( ( scop_random() < (rate[i] / (1000/dt) ) ) && (last_quanta[i] >= AbsRefract) )
     {
       release_start[i] = t
       last_quanta[i] = 0
       RVP_out[i] = 1
       total_count = total_count + 1
     }
gluConc_dAMPA= gluConc_dAMPA + (glumax_dAMPA * alpha( (t - release_start[i])/tau_dAMPA ))
gluConc_NMDA = gluConc_NMDA + (glumax_NMDA * alpha( (t - release_start[i])/tau_NMDA ))
  }
}

PROCEDURE RVPSIZE(){
FROM i=0 TO NumSites-1
{
RVP_Size[i] = RVP_Size[i]  - RVP_out[i]
if (RVP_Size[i] < 0) {RVP_Size[i] = 0}
if (RVP_Size[i] > Max_RVP){RVP_Size[i] = Max_RVP}
: refilling of RVP
if (RVP_Size[i] < Max_RVP){
rateIN[i] = rate_constantIN  *  (Max_RVP - RVP_Size[i]) 
 if (scop_random() < rateIN[i] / (1000 / dt))
  { RVP_Size[i] = RVP_Size[i] + 1 }
}
RVP_out[i] = 0
}
}

KINETIC kstates {
     rb = Rb * (1e3) * gluConc_dAMPA		:Initial lower case 'r' indicates
     rBb = RBb * (1e3) * gluConc_dAMPA		:the rates are ligand dependent
     rDBb = RDBb * (1e3) * gluConc_dAMPA		:and are expressed in units of	
     rMb = RMb * (1e3) * gluConc_dAMPA 		:/micromolar/microsec
     rMBb = RMBb * (1e3) * gluConc_dAMPA 
     rm  = Rm * (1e3) * mod_modulator		:Initial lower case 'r' indicates
     rBm = RBm * (1e3) * mod_modulator   	:the rates are dependent on the 
     rB2m = RB2m * (1e3) * mod_modulator	:concentration of the modulator
     rDBm = RDBm * (1e3) * mod_modulator	:(cyclothiozide or aniracetam)
     rDB2m = RDB2m * (1e3) * mod_modulator	:units are /micromolar/microsec

     ~ R <-> RA     (rb,RBu)
     ~ RA <-> RdA   (RBd,RBDr)
     ~ RA <-> RA2   (rBb,RB2u)
     ~ RA2 <-> RdA2 (RB2d,RB2Dr)
     ~ RA2 <-> O    (RB2o,RB2c)
     ~ RdA <-> RdA2 (rDBb,RDB2u)
     ~ R <-> RM     (rm,RMum)
     ~ RA <-> RMA   (rBm,RBMum)
     ~ RA2 <-> RMA2 (rB2m,RB2Mum)
     ~ RM <-> RMA   (rMb,RMBu)
     ~ RMA <-> RMA2 (rMBb,RMB2u)
     ~ RMA2 <-> OM   (RMB2o,RMB2c)
     ~ RMA2 <-> RdMA2 (RMB2d,RMB2Dr)
     ~ RMA <-> RdMA (RMBd,RMBDr)
     ~ RdMA <-> RdMA2 (RMDBb,RMDB2u)
     ~ RdA <-> RdMA (rDBm,RDBMum)
     ~ RdA2 <-> RdMA2 (rDB2m,RDB2Mum)

     CONSERVE R+RM+RA+RA2+RdA+RdA2+RMA+RMA2+RdMA+RdMA2+O+OM = 1


       rb2 = Rb2 * (1e3) * gluConc_NMDA

        ~ C0 <-> C1     (rb2,Ru)
        ~ C1 <-> C2     (rb2,Ru)
        ~ C2 <-> D      (Rd,Rr)
        ~ C2 <-> OO     (Ro,Rc)

        CONSERVE C0+C1+C2+D+OO = 1
}

PROCEDURE rates(v(mV)) {
        TABLE B
        DEPEND mg
        FROM vmin TO vmax WITH 200

        : Stevens & Jahr 1990a,b

        B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}


FUNCTION alpha(x) {
        if (x < 0 || x > 10) {
                alpha = 0
        }else{
                alpha = x * exp(1 - x)
        }
}



