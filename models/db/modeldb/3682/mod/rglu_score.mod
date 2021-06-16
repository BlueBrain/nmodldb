COMMENT
$Id: rglu_score.mod,v 1.4 2000/08/14 22:21:27 karchie Exp $

-----------------------------------------------------------------------------
Simple synaptic mechanism derived for first order kinetics of
binding of transmitter to postsynaptic receptors.

A. Destexhe & Z. Mainen, The Salk Institute, March 12, 1993.
-----------------------------------------------------------------------------

Now includes both voltage-dependent (NMDA) and V-independent (AMPA/kainate)
components.  This model also keeps an efficacy score for modeling synaptic
plasticity.

Revised to remove include assert.h (NEURON will supply own version).  2004/07/26

Split from rnmda_score.mod, version 1.2
(in turn split from rnmda.mod, version 1.7)

$Log: rglu_score.mod,v $
Revision 1.4  2000/08/14 22:21:27  karchie
Now using argument of exprand() correctly.

Revision 1.3  2000/06/05 19:34:17  karchie
Corrected declaration of double exprand().

Revision 1.2  2000/06/05 19:32:28  karchie
Now gets exponentially distributed r.v. from SCoPlib instead of random.mod.

Revision 1.1  1999/05/13 20:19:20  karchie
Initial revision


-----------------------------------------------------------------------------


During the arrival of the presynaptic spike (detected by threshold 
crossing), it is assumed that there is a brief pulse (duration=Cdur)
of neurotransmitter C in the synaptic cleft (the maximal concentration
of C is Cmax).  Then, C is assumed to bind to a receptor Rc according 
to the following first-order kinetic scheme:

Rc + C ---(Alpha)--> Ro							(1)
       <--(Beta)--- 

where Rc and Ro are respectively the closed and open form of the 
postsynaptic receptor, Alpha and Beta are the forward and backward
rate constants.  If R represents the fraction of open gates Ro, 
then one can write the following kinetic equation:

dR/dt = Alpha * C * (1-R) - Beta * R					(2)

and the postsynaptic current is given by:

Isyn = gmax * R * (V-Erev)						(3)

where V is the postsynaptic potential, gmax is the maximal conductance 
of the synapse and Erev is the reversal potential.

If C is assumed to occur as a pulse in the synaptic cleft, such as

C     _____ . . . . . . Cmax
      |   |
 _____|   |______ . . . 0 
     t0   t1

then one can solve the kinetic equation exactly, instead of solving
one differential equation for the state variable and for each synapse, 
which would be greatly time consuming...  

Equation (2) can be solved as follows:

1. during the pulse (from t=t0 to t=t1), C = Cmax, which gives:

   R(t-t0) = Rinf + [ R(t0) - Rinf ] * exp (- (t-t0) / Rtau )		(4)

where 
   Rinf = Alpha * Cmax / (Alpha * Cmax + Beta) 
and
   Rtau = 1 / (Alpha * Cmax + Beta)

2. after the pulse (t>t1), C = 0, and one can write:

   R(t-t1) = R(t1) * exp (- Beta * (t-t1) )				(5)

The variables "phase" and "mean_ia_time" control presynaptic spike arrival:
  * phase, if in [0, 1), determines the time of the first of the regular
    arrivals.  If phase < 0, arrivals are a Poisson process.
  * mean_ia_time > 0 is the mean interspike interval.  If arrivals are
    regular, mean_ia_time > Cdur + Deadtime.

ENDCOMMENT

VERBATIM
/* #include <assert.h>  commented out because NEURON has own assert */

static char rcsid[] = "$Id: rglu_score.mod,v 1.4 2000/08/14 22:21:27 karchie Exp $";

extern double exprand(double mean);

ENDVERBATIM


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS rGluSc
	RANGE C, gp_NMDA, lastrelease, score, idx, phase, mean_ia_time, gmax_NMDA, R_NMDA, R0_NMDA, R1_NMDA, gp_AMPA, gmax_AMPA, R_AMPA, R0_AMPA, R1_AMPA
	NONSPECIFIC_CURRENT i
	GLOBAL Cmax, Cdur, Deadtime, score_thresh, score_tau, Alpha_NMDA, Beta_NMDA, Erev_NMDA, Rinf_NMDA, Rtau_NMDA, Alpha_AMPA, Beta_AMPA, Erev_AMPA, Rinf_AMPA, Rtau_AMPA
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	dt (ms)
	: Transmitter release parameters shared by AMPA/NMDA receptors.
	Cmax	= 1	(mM)		: max transmitter concentration
	Cdur	= 1.1	(ms)		: transmitter duration (rising phase)
	Deadtime = 0	(ms)		: minimum time between release events
	phase	= 0			: event time offset / poisson flag
	mean_ia_time = 100 (ms)		: event interarrival time

	: Transmitter binding differs between AMPA and NMDA.
	Alpha_NMDA = 10	(/ms mM)	: forward (binding) rate
	Beta_NMDA = 0.0125 (/ms)	: backward (unbinding) rate
	Alpha_AMPA = 10 (/ms mM)	: forward rate
	Beta_AMPA = 0.5 (/ms)		: backward rate

	: Reversal potential may differ.
	Erev_NMDA	= 0	(mV)		: reversal potential
	Erev_AMPA = 0	(mV)

	: Conductances are specified separately.
	gmax_NMDA	(umho)		: maximum conductance
	gmax_AMPA	(umho)		: maximum conductance

	: Scoring mechanism applies only to NMDA.
	score_thresh = -63	(mV)	: minimum voltage for scoring
        score_tau = 20	(ms)		: time constant for score

	: NMDA kinetics are more complicated.
	eta     = 0.33  (/mM)
	mag     = 1     (mM)
	gamma   = 0.06  (/mV)
}

ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	gp_NMDA		(umho)		: NMDA conductance component
	gp_AMPA		(umho)		: AMPA conductance component

	depol		(mV)		: depolarization = (v - 70)
	C		(mM)		: transmitter concentration

	R_NMDA				: fraction of open NMDA channels
	R0_NMDA				: open channels at start of release
	R1_NMDA				: open channels at end of release
	R_AMPA				: fraction of open AMPA channels
	R0_AMPA				: open channels at start of release
	R1_AMPA				: open channels at end of release

	Rinf_NMDA			: steady state channels open
	Rtau_NMDA	(ms)		: time constant of channel binding
	Rinf_AMPA			: steady state channels open
	Rtau_AMPA	(ms)		: time constant of channel binding

	lastrelease	(ms)		: time of last spike
	ia_time		(ms)		: time between last and next spike
	score				: efficacy score
	score_postsyn_act		: local/recent postsynaptic activity
	idx				: site index (~= synapse location)
}

INITIAL {
	C = 0

	R_NMDA = 0
	R1_NMDA = R_NMDA
	R_AMPA = 0
	R1_AMPA = R_AMPA

	Rinf_NMDA = Cmax*Alpha_NMDA / (Cmax*Alpha_NMDA + Beta_NMDA)
	Rtau_NMDA = 1 / ((Alpha_NMDA * Cmax) + Beta_NMDA)

	Rinf_AMPA = Cmax*Alpha_AMPA / (Cmax*Alpha_AMPA + Beta_AMPA)
	Rtau_AMPA = 1 / ((Alpha_AMPA * Cmax) + Beta_AMPA)

	score = 0
	score_postsyn_act = 0

	: Make sure mean IA time and phase are valid.
	: Valid phase values for regular firing are 0.0 <= phase < 1.0
	: If (phase < 0.0), events are a Poisson process.
VERBATIM
	assert(phase < 1.0);
	if (phase < 0.0) {
		assert(mean_ia_time > 0.0);
	} else {
		assert(mean_ia_time > Cdur + Deadtime);
	}
ENDVERBATIM

	: Set the first spike time.
	if (phase < 0.0) {
		: Events are a Poisson process.
		lastrelease = t - next_ia_time()
		ia_time = next_ia_time()
	} else {
		: Use phase to set first spike time.
		ia_time = mean_ia_time
		lastrelease = t - phase * mean_ia_time
	}	
}

BREAKPOINT {
	SOLVE release
	gp_AMPA = gmax_AMPA * R_AMPA
	gp_NMDA = (gmax_NMDA * R_NMDA)/(1 + eta * mag * exp( - (gamma * v)))
	i = gp_AMPA*(v - Erev_AMPA) + gp_NMDA*(v - Erev_NMDA)

	: Update the score: presyn * postsyn
	depol = v - score_thresh	: depolarization
	if (depol < 0) {
		depol = 0	: don't count hyperpolarization.
	}

	: Postsynaptic activity lags the depolarization.
	score_postsyn_act = score_postsyn_act + dt * (depol - score_postsyn_act) / score_tau

	score = score + dt * gp_NMDA * score_postsyn_act
}

PROCEDURE release() { LOCAL q
	: First determine whether we should start the next release,
	: end the previous release, or neither.
	q = t - lastrelease
	if (q >= ia_time) {			: time for the next release?
		C = Cmax
		R0_NMDA = R_NMDA
		R0_AMPA = R_AMPA
		lastrelease = t
		ia_time = next_ia_time()	: calculate next ia time.
	} else if (q < Cdur) {			: still releasing?
		: do nothing
	} else if (C == Cmax) {			: in dead time after release
		R1_NMDA = R_NMDA
		R1_AMPA = R_AMPA
		C = 0.
	}		
	
	if (C > 0) {				: transmitter being released?
	   R_NMDA = Rinf_NMDA + (R0_NMDA - Rinf_NMDA) * exptable (- (t - lastrelease) / Rtau_NMDA)
	   R_AMPA = Rinf_AMPA + (R0_AMPA - Rinf_AMPA) * exptable (- (t - lastrelease) / Rtau_AMPA)
	} else {				: no release occuring
  	   R_NMDA = R1_NMDA * exptable(-Beta_NMDA * (t - (lastrelease + Cdur)))
  	   R_AMPA = R1_AMPA * exptable(-Beta_AMPA * (t - (lastrelease + Cdur)))
	}

	VERBATIM
	return 0;
	ENDVERBATIM
}


FUNCTION next_ia_time() {
	if (phase < 0.0) {
		: Arrivals are a Poisson process
		next_ia_time = 0
		while (next_ia_time <= Cdur + Deadtime) {
VERBATIM
			_lnext_ia_time =
				exprand(mean_ia_time);
ENDVERBATIM
		}
	} else {
		: Arrivals are regular
		next_ia_time = mean_ia_time
	}
}


FUNCTION exptable(x) { 
	TABLE  FROM -10 TO 10 WITH 2000

	if ((x > -10) && (x < 10)) {
		exptable = exp(x)
	} else {
		exptable = 0.
	}
}
