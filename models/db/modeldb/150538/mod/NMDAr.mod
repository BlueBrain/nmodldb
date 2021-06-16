: Major modification: add NET_RECEIVE block to be controled by NetCon

TITLE Asymmetric trapping block model of NMDA receptors

COMMENT
------------------------------------------------------------------------ 
-----

     Asymmetric trapping block model of NMDA receptors
     ===============================

     See:
     
     Vargas-Caballero M and Robinson HPC (2004). "Fast and slow voltage-dependent dynamics 
     of magnesium block in the NMDA receptor: the asymmetric trapping block model", J. Neurosci. 24:6171-6180.
     
     10-state gating model:

     Modified from Sobolevsky and Yelshansky, 2000.
     Asymmetric rate constants for Mg bound and unbound states

                  D
                  |
     C  -- C1  -- C2  --  O      Mg-free states
                          |   ---------------------
     CB -- C1B -- C2B -- OB      Mg-bound states
                  |
                  DB

     Voltage dependence of Mg2+ block and unblock reactions 
     from Ascher and Nowak 1988
     
     This version applies to room temperature, 1 mM [Mg2+]_o
------------------------------------------------------------------------ 
-----

------------------------------------------------------------------------ 
-----

M. Vargas Caballero, H.P.C. Robinson and A. Roth, 2005

------------------------------------------------------------------------ 
-----
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
     POINT_PROCESS NMDA_TESTED
     RANGE state_C0, state_C1, state_C2, state_D, state_O, state_B,  
state_DB, state_C2B, state_C1B, state_CB
     RANGE g, nchan, gamma, rb, mg_on, mg_off, C, conc, del, dur
     GLOBAL Erev, mg, kon, koff, kd, kr, beta, alpha
     GLOBAL beta_mg, alpha_mg, koff_mg, rb_mg, kon_mg, kr_mg, kd_mg
     GLOBAL vmin, vmax
     NONSPECIFIC_CURRENT i
}

UNITS {
     (nA) = (nanoamp)
     (mV) = (millivolt)
     (uS) = (microsiemens)
     (umho) = (micromho)
     (mM) = (milli/liter)
     (uM) = (micro/liter)
}

PARAMETER {

     Erev    = -3    (mV) : reversal potential
     nchan  = 10          : number of channels
     gamma = 0.00005 (uS) : single channel conductance
     mg  = 1     (mM)     : external magnesium concentration
     vmin = -180 (mV)
     vmax = 100  (mV)
     normfactor = .06948171

: Rates

     : Lester & Jahr
     kon  = 5e-3    (/uM /ms)    : binding
     koff  = 82e-3  (/ms)      : unbinding
     beta  = 46.5e-3   (/ms)     : opening
     alpha  = 91.6e-3   (/ms)    : closing
     kr  = 1.8e-3   (/ms)        : resensitization
     kd  = 8.4e-3   (/ms)        : desensitization

     kon_mg  = 5e-3    (/uM /ms) : binding + mg
     koff_mg  = 82e-3  (/ms)   : unbinding + mg
     beta_mg  = 41.85e-3   (/ms)  : opening + mg
     alpha_mg  = 229e-3   (/ms) : closing + mg
     kr_mg  = 1.8e-3   (/ms)     : resensitization from the blocked state + mg
     kd_mg  = 8.4e-3   (/ms)     : desensitization from the blocked state + mg
     del = 100 (ms)
     dur = 10 (ms)
     conc = 0       (uM)   : transmitter concentration
}

ASSIGNED {
     v       (mV)        : postsynaptic voltage
     i       (nA)        : current = g*(v - Erev)
     g       (uS)        : conductance
     rb      (/ms)       : binding
     rb_mg   (/ms)   : binding of glutamate to blocked channels
     mg_on   (/mM /ms)       : blocking rate
     mg_off  (/ms)       : unblocking rate
     C       (uM)       : concentration calculated by the mod file
}

STATE {
     : Channel states (all fractions)
     state_C0       : unbound
     state_C1       : single bound
     state_C2       : double bound
     state_D        : desensitized
     state_O        : open
     state_B        : blocked
     state_DB      : desensitised closed
     state_C2B     : double bound closed
     state_C1B     : single bound closed
     state_CB      : unbound closed
}

INITIAL {
SOLVE kstates STEADYSTATE sparse
}

BREAKPOINT {
     transmitter()
     SOLVE kstates METHOD sparse

     g = nchan * state_O * gamma : scale the single channel conductance
     i = g * (v - Erev) : current in nA (mV*uS)
}

KINETIC kstates {
	rates(v)
     rb = kon * C
     rb_mg = kon_mg * C

     ~ state_C0  <-> state_C1   (2*rb,koff)
     ~ state_C1  <-> state_C2   (rb,2*koff)
     ~ state_C2  <-> state_D    (kd,kr)
     ~ state_C2  <-> state_O    (beta,alpha)
     ~ state_O   <-> state_B    (mg_on,mg_off)
     ~ state_C2B <-> state_DB   (kd_mg,kr_mg)
     ~ state_B   <-> state_C2B  (alpha_mg,beta_mg)
     ~ state_C2B <-> state_C1B  (2*koff_mg,rb_mg)
     ~ state_C1B <-> state_CB   (koff_mg,2*rb_mg)

     CONSERVE  
state_C0+state_C1+state_C2+state_D+state_O+state_B+state_DB+state_C2B+state_C1B+state_CB = 1
}

PROCEDURE rates(v(mV)) {

     : from Ascher and Nowak - mg replaced with 1

     mg_on = 610*exp(-v/17)*(1/1000)
     mg_off = 5.4*exp(v/47)
}

: Use procedure only when working in non-stationary conditions,
: when in use, uncomment the call for this
: procedure at the breakpoint

: a brief square pulse may reproduce the synaptic activation
: of the receptor

:///////////////////////////////////////////////////////////////////////////////////////////////////////////////
: spike time is determined by NET_RECEIVE; in this case, w is the value assigned to GLUT concentration
:///////////////////////////////////////////////////////////////////////////////////////////////////////////////

NET_RECEIVE(w (uM)) {

  : state_discontinuity(state_O, state_O + w/normfactor)
           
       del=t
       conc=w

}

PROCEDURE transmitter() {
if ((conc!=0)&&(t>=del)&&(t<=del+dur)) {
   C=conc
} else {
   C=0
   conc=0
}
}

