TITLE model of GABAB receptors

COMMENT
-----------------------------------------------------------------------------

	Kinetic model for GABA-B receptors
	==========================================

	Model of GABAB currents including nonlinear stimulus 
	dependency (fundamental to take into account for GABAB receptors).


	Features:

	  - peak at ~200 ms after burst activation (5@50 Hz); time course fit from experimental IPSPs recorded by J. Schulz
	  - NONLINEAR SUMMATION (psc is much stronger with bursts)
		due to cooperativity of G-protein binding on K+ channels

	Approximations:

	  - single binding site on receptor	
	  - model of alpha G-protein activation (direct) of K+ channel
	  - G-protein dynamics is second-order; simplified as follows:
		- saturating receptor
		- no desensitization
		- Michaelis-Menten of receptor for G-protein production
		- "resting" G-protein is in excess
		- Quasi-stat of intermediate enzymatic forms
	  - binding on K+ channel is fast


	Kinetic Equations of model:

	  dT/dt = -T/tauD -k1 * T * (Bm - B) + k_1 * B 
	  dB/dt = k1 * T * (Bm - B) - (k_1 + k2) * B
	  dR/dt = K1 * T * (1-R) - K2 * R
	  dG/dt = (K3 * R * (1-G) - K4 * G) *f

      
	  R : fraction activated receptor
	  T : transmitter
      B : GABA transporter
	  G : fraction activated G-protein
	  K1,K2,K3,K4 = kinetic rate cst; from Thomson & Destexhe, 1999, Fig. 15 for n=2
      k1,k_1,k2 = kinetic rate cst; from Thomson & Destexhe, 1999
      tauD : decay due to diffusion; from Sanders et al., 2013
      f : factor f to G protein control dynamics
      
  f and K2 adjusted to reach max amplitude ~200 ms after burst start (5@50 Hz)

  n activated G-protein bind to a K+ channel:

	n G + C <-> O		(Alpha,Beta)

  If the binding is fast, the fraction of open channels is given by:

	O = G^n / ( G^n + KD )

  where KD = Beta / Alpha is the dissociation constant

-----------------------------------------------------------------------------

  Also see details in:

  Destexhe, A. and Sejnowski, T.J.	G-protein activation kinetics and
  spill-over of GABA may account for differences between inhibitory responses
  in the hippocampus and thalamus.	Proc. Natl. Acad. Sci. USA	92:
  9515-9519, 1995.

  Thompson, A.M. and Destexhe, A. DUAL INTRACELLULAR RECORDINGS AND COMPUTATIONAL
  MODELS OF SLOW INHIBITORY POSTSYNAPTIC POTENTIALS IN RAT NEOCORTICAL AND HIPPOCAMPAL 
  SLICES. Neuroscience 92: 1193-1215, 1999.
  
  Sanders, H., Berends, M., Major, G., Goldman, M.S. and Lisman, J.E. NMDA and 
  GABAB (KIR) conductances: the "perfect couple" for bistability. J Neurosci 33(2): 424-9, 2013.
  
  Taken from Poirazi, Brannon & Mel. Arithmetic of Subthreshold Synaptic
  Summation in a Model CA1 Pyramidal Cell. Neuron 2003 (Originally written by Alain Destexhe, Laval University, 1995)
  
  Modified by J. Schulz according to Thompson & Destexhe (1999) and Sanders, Berends et al. (2013) 

-----------------------------------------------------------------------------
ENDCOMMENT

NEURON {
	POINT_PROCESS GABABsyn
	RANGE C, R, G, B, g, gmax, tauD
	NONSPECIFIC_CURRENT i
	RANGE vgat,sst,npy,pv,xEff
	RANGE isOn
	GLOBAL K1, K2, K3, K4, KD, k1, k_1, k2, e, Bm
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(uS) = (microsiemens)
}

PARAMETER {

	tauD = 10	(ms)		: decay of transmitter concentration
	K1	= 0.066	(/ms mM)	: forward binding rate to receptor
	K2	= 0.008 (/ms)		: backward (unbinding) rate of receptor
	K3	= 0.27 (/ms)		: rate of G-protein production
	K4	= 0.044 (/ms)		: rate of G-protein decay
	KD	= 0.5				: half maximal coductance at a level of ~0.7 activated G-protein
	n	= 2			: nb of binding sites of G-protein on K+
	e	= -95	(mV)		: reversal potential (E_K)
	gmax		(uS)		: maximum conductance
    f   = 0.1              : factor f controlling the G protein dynamics
	k1	= 30	(/ms mM)	: 30, forward binding rate to transporter
	k_1	= 0.1 (/ms)		: backward (unbinding) rate of transporter
	k2	= 0.02 (/ms)		: clearance of GABA
	Bm = 1 (mM)			: maximum binding capacity of transporter
	vgat=0
	sst=0
	npy=0
	pv=0
	xEff=-1
	isOn=0
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i		(nA)		: current = g*(v - e)
	g		(uS)		: conductance
	Gn
}


STATE {
	C	(mM)		: extracellular transmitter concentration
	R				: fraction of activated receptor
	G				: normalized concentration of activated G-protein
	B	(mM)		: bound GABA transporter
}


INITIAL {
	C = 0
	R = 0
	G = 0
	B = 0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	Gn = G^n
	g = isOn * gmax * Gn / (Gn+KD)
	i = g *(v - e)
}


DERIVATIVE state {

	C' = (-C/tauD -k1 * C * (Bm - B) + k_1 * B) 
	R' = (K1 * C * (1-R) - K2 * R) 
	G' = (K3 * R * (1-G) - K4 * G) * f
	B' = (k1 * C * (Bm - B) - (k_1 + k2) * B) 

}


NET_RECEIVE(weight (mM)) {
	C = C + weight
}
