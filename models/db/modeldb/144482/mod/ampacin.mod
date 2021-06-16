TITLE minimal model of AMPA receptors

COMMENT
-----------------------------------------------------------------------------

        Minimal kinetic model for glutamate AMPA receptors
        ==================================================

  Model of Destexhe, Mainen & Sejnowski, 1994:

        (closed) + T <-> (open)

  The simplest kinetics are considered for the binding of transmitter (T)
  to open postsynaptic receptors.   The corresponding equations are in
  similar form as the Hodgkin-Huxley model:

        dr/dt = alpha * [T] * (1-r) - beta * r

        I = gmax * [open] * (V-Erev)

  where [T] is the transmitter concentration and r is the fraction of 
  receptors in the open form.

  If the time course of transmitter occurs as a pulse of fixed duration,
  then this first-order model can be solved analytically, leading to a very
  fast mechanism for simulating synaptic currents, since no differential
  equation must be solved (see Destexhe, Mainen & Sejnowski, 1994).

-----------------------------------------------------------------------------

  Based on voltage-clamp recordings of AMPA receptor-mediated currents in rat
  hippocampal slices (Xiang et al., J. Neurophysiol. 71: 2552-2556, 1994), this
  model was fit directly to experimental recordings in order to obtain the
  optimal values for the parameters (see Destexhe, Mainen and Sejnowski, 1996).

-----------------------------------------------------------------------------

  This mod file includes a mechanism to describe the time course of transmitter
  on the receptors.  The time course is approximated here as a brief pulse
  triggered when the presynaptic compartment produces an action potential.
  The pointer "pre" represents the voltage of the presynaptic compartment and
  must be connected to the appropriate variable in oc.

-----------------------------------------------------------------------------

  See details in:

  Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  An efficient method for
  computing synaptic conductances based on a kinetic model of receptor binding
  Neural Computation 6: 10-14, 1994.  

  Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  Kinetic models of 
  synaptic transmission.  In: Methods in Neuronal Modeling (2nd edition; 
  edited by Koch, C. and Segev, I.), MIT press, Cambridge, 1998, pp. 1-25.

    (electronic copy available at http://cns.iaf.cnrs-gif.fr)


  Written by Alain Destexhe, Laval University, 1995

Modified by M. Badoual, 2004

-----------------------------------------------------------------------------
ENDCOMMENT



INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
        POINT_PROCESS AMPAKIT
	  RANGE onset,periodpre, periodpost, delta,nbrepre, nbrepost, change, tau0, tau1, g,gmax, e, i,C
	  RANGE Cmax
	  NONSPECIFIC_CURRENT i
      GLOBAL Erev,Cdur
}
UNITS {
    (celsius) = (degC)        
	(nA) = (nanoamp)
    (mV) = (millivolt)
    (umho) = (micromho)
    (mM) = (milli/liter)
}

PARAMETER {
	onset = 100  (ms)
	periodpre = 0 (ms)	:periode
	periodpost=20
	delta= 10 (ms)		:temps entre pre et post
	nbrepre=1			:nbre de repetitions
	nbrepost=1
	tau0 = 1	 (ms)		: 0.34
	tau1 = 3	 (ms)		: 2.0
    Erev = 0    (mV)            : reversal potential
	Cmax = 1	(mM)		: max transmitter concentration
	Cdur	= 1	(ms)		: transmitter duration (rising phase)
    gmax = 0.001 (umho)          : maximum conductance (original value 0.002)
}


ASSIGNED {
        v               (mV)            : postsynaptic voltage
        i               (nA)            : current = g*(v - Erev)
        g               (umho)          : conductance
        C		(mM)		: transmitter concentration
	change
}

LOCAL   a[2]
LOCAL   tpeak
LOCAL   adjust
LOCAL   amp

INITIAL {
	C = 0
}


BREAKPOINT {
    g = cond(t,onset)
	C = trans(t,onset)
	
	if (nbrepre>1) {
	  FROM j=1 TO (nbrepre-1) {
	    g = g+cond(t,onset+j*periodpre)
		C = C+trans(t,onset+j*periodpre)
		
		}
	}
        i = g*(v - Erev)
}

FUNCTION myexp(x) {
	if (x < -100) {
	myexp = 0
	}else{
	myexp = exp(x)
	}
}

FUNCTION cond(x (ms), onset1 (ms)) (umho) {
	tpeak=tau0*tau1*log(tau0/tau1)/(tau0-tau1)
	adjust=1/((1-myexp(-tpeak/tau0))-(1-myexp(-tpeak/tau1)))
	amp=adjust*gmax
	if (x < onset1) {
		cond = 0
	}else{
		a[0]=1-myexp(-(x-onset1)/tau0)
		a[1]=1-myexp(-(x-onset1)/tau1)
		cond = amp*(a[0]-a[1])
	}
}

FUNCTION trans(x (ms), onset1 (ms)) (mM) {
	if ((x>onset1) && (x-onset1<=Cdur)) {
		trans=Cmax
	} else {
		trans=0
	}
}



