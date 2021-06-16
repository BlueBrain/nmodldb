COMMENT
Ribbon Synapse with inline calculations rather than function calls from BREAKPOINT
Log-Normal evaluation for glutamate profile
ENDCOMMENT

DEFINE NumSites 10 : The number of Release Sites per Synapse

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
     POINT_PROCESS ribbon_tiger
     EXTERNAL mod_modulator
: Two presynaptic compartments are used with identical characteristics to allow the measurement of calcium concentration
: at two depths. Future implementations might utilize a more elegant mechanism.
     POINTER preCA1, preCA2
:  for the AMPA receptor
     RANGE preCAnow
     RANGE R, RM, RA, RA2, RdA, RdA2, RMA, RMA2, RdMA, RdMA2, O, OM
     RANGE g, gmax, erev, rb, rBb, rDb, rDBb, rMb, rMBb, rm, rBm, rB2m, rDBm, rDB2m
     RANGE Rb, RBu, RBd, RBDr, RBb, RB2u
     RANGE RB2d, RB2Dr, RB2o, RB2c, RDBb, RDB2u
     RANGE Rm, RMum, RBm, RBMum, RB2m, RB2Mum
     RANGE RMb, RMBu, RMBb, RMB2u, RMB2o, RMB2c
     RANGE RMB2d, RMB2Dr, RMBd, RMBDr, RMDBb, RMDB2u
     RANGE RDBMum, RDBm, RDB2Mum, RDB2m
     RANGE AbsRefract,  gluConc_dAMPA, pr, StanddAMPA, rate_constant
:  for the NMDA receptor
     RANGE g2, gmaxN, erevN, rb2, mg
     RANGE C0, C1, C2, D, OO, B
     RANGE Rb2, Ru, Rd, Rr, Ro, Rc
     RANGE gluConc_NMDA, StandNMDA
:
     RANGE XC,w,A
     RANGE Max_RVP,  rate_constantIN
     RANGE Vmax, k, n
     RANGE UFP, Orate
     RANGE S,W,K1,powr
     GLOBAL total_count

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
     preCA1 (mM) : presynaptic [Ca++] for a full RVP
     preCA2 (mM) : presynaptic [Ca++] for a less than full RVP
     preCAnow
     erev = 0 (mV)	: reversal potential of dAMPA
     erevN = 0 (mV)     : reversal potential of NMDA
     gmax = 0 (umho)	: maximal conductance dAMPA
     gmaxN = 0 (umho)   : max conductance NMDA
     StanddAMPA = 0 (mM) : standing level of glutamate seen by dAMPA
     StandNMDA = 0.001 (mM)  : standing level of glutamate seen by NMDA
	xC0 = 0.442216   : Allow for the dynamic setting of NMDA rate constants at HOC level
	xC1 = 0.23271    :       "
	xC2 = 0.12248    :       "
	xD =  0.15075    :       "
	xOO = 0.0519     :       "
     mg = 0 (mM)        : [Mg++] for NMDA
        AbsRefract = 0 (ms) : absolute refractory period for vesicular release from a site
        rate_constant (1/s) : Rate constant for exocytosis
        rate_constantIN = 0.25 (1/s) : Rate constant of RVP uptake
        XC : for Log-Normal Eqation
        w  :        "
        A  :        "
        Max_RVP = 5 : Maximum number of vesicles in RVP for one release site
        Vmax = 1840.22 : for Hill equation
        k = 86.48823 :        "
        n = 3.30393  :        "
        S
        W
        K1
        powr
: Rates

: AMPA from Partin, Fleck & Mayer (1996) J Nsci 16, 6634-47
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

: AMPA modulator kinetics OFF
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

: NMDA Hessler NA, Shirke AM, Malinow R (1993)
        Rb2     = 5e-3    (/uM /ms)     : binding
        Ru      = 9.5e-3   (/ms)        : unbinding
        Rd      = 16e-3    (/ms)        : desensitization
        Rr      = 13e-3    (/ms)        : resensitization
        Ro      = 25e-3     (/ms)       : opening
        Rc      = 59e-3     (/ms)       : closing

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
     gluConc_dAMPA    (mM) : glutamate concentration seen by AMPA receptors
     gluConc_NMDA     (mM) : glutamate concentration seen by NMDA receptors
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
     rate[NumSites] (/s) : Current rate of vesicular release
     rateIN[NumSites] (/s) :  Rate of RVP refilling
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
: Set standing levels of glutamate
gluConc_dAMPA = StanddAMPA
gluConc_NMDA = StandNMDA

Orate = 0
: Calculate size of ultra-fast pool
UFP = 0
FROM i=0 TO NumSites-1
{ 
if (RVP_Size[i] == Max_RVP)
   {UFP = UFP + 1}
}
: ------------------------------------------------------------------------------- 
: vesicular release

FROM i=0 TO NumSites-1 
  {
: Determine at which depth [Ca++] is determined
if (RVP_Size[i] < Max_RVP) {preCAnow = preCA2}
else
{preCAnow = preCA1}

: Calculate the instantaneous rate of release for every release site.
: This rate is an average for a stochastic process.
if (RVP_Size[i] < 1)
: When RVP pool is empty use  Rouze & Schwartz (1998) J.Neurosci.18:8614-8624 Fig.7 scaled for a single release site
       { S = 2.84985  W = 0.0498303 K1 = 18.3892 powr=4
         rate[i] = S * (W * preCAnow * 1000) / ( ((preCAnow * 1000) / K1) + 1)^powr }
: factor of 1000 for conversion from NEURON mM to formula uM
else : When RVP is not empty
{
: Determine the rate constant based on the [Ca++]
: factor of 1000 for conversion from NEURON mM to formula uM
: Function of Heidelberger (1994),Fig3a is used, fit by a Hill equation

  rate_constant =  (Vmax * (1000 * preCAnow)^n) / (k^n + (1000 * preCAnow)^n)

  if (RVP_Size[i] == Max_RVP)
   { rate[i] = UFP * rate_constant }
  else
   { rate[i] =  RVP_Size[i] * rate_constant }
}
 
Orate = Orate + rate[i]

: Release is stochastic but dependent on rate of release calculated above
  last_quanta[i] = last_quanta[i] + dt
  if ( ( scop_random() < (rate[i] / (1000/dt) ) ) && (last_quanta[i] >= AbsRefract) )
     {
       release_start[i] = t
       last_quanta[i] = 0
       RVP_out[i] = 1
       total_count = total_count + 1
     }

: A glutamate concentration profile is generated for  AMPA and NMDA receptors based on Mcell simulations
: fit to a Log-Normal function
XC=0.09 w=1.28093 A=0.09262
gluConc_dAMPA= gluConc_dAMPA + transmitter(t - release_start[i])
XC=0.19 w=1.29098 A=0.0256
gluConc_NMDA = gluConc_NMDA + transmitter(t - release_start[i])
  }
: ---------------------------------------------------------------------------------------------------

: Refill of pools
FROM i=0 TO NumSites-1
{
RVP_Size[i] = RVP_Size[i]  - RVP_out[i]
if (RVP_Size[i] < 0) {RVP_Size[i] = 0}
if (RVP_Size[i] > Max_RVP){RVP_Size[i] = Max_RVP}
: refilling of RVP is stochastic, its rate dependent on the size of the RVP
if (RVP_Size[i] < Max_RVP){
rateIN[i] = rate_constantIN  *  (Max_RVP - RVP_Size[i]) 
 if (scop_random() < rateIN[i] / (1000 / dt))
  { RVP_Size[i] = RVP_Size[i] + 1 }
}
RVP_out[i] = 0
}

SOLVE kstates METHOD sparse

     g = gmax *(O+OM) : AMPA conductance
     g2 = gmaxN * OO * B : NMDA conductance
     i = (g * (v - erev) ) + (g2 * (v - erevN) )
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


FUNCTION transmitter(x) {
: Log-Normal function
        if (x <= 0 || x > 10) {
             transmitter = 0
        }else{
             transmitter = A/(sqrt(2*3.1415927)*w*x) * exp(-log(x/XC)^2 / (2*w^2))
        }
}
