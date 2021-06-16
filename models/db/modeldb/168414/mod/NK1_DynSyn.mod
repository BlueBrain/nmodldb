TITLE (very)-Empirical model for the effect of NK1 receptor activation

COMMENT

 The effects of the activation of NK1 receptors are not totally
 clear but they seem to include both a membrane conductance change
 and an increase in intra-cellular calcium concentration (assumed
 to result from the release of calcium from intra-cellular buffers)
 The dynamics of this receptor should be used in association with
 dynamics for intra-cellular calcium dynamics. It affects all 
 calcium-dependent currents!

 We have based this model in the paper by Ito et al, 2002, with the
 title: "Substance P mobilizes intracellular calcium and activates
 a nonselective cation conductance in rat spiral ganglion neurons."
 The model includes therefore two mechanisms:
 a) a slow conductance change generating a depolarizing non-specific
    cationic current iNK1R
 b) an increase in cai
 Details are given in the code, in the form of comments.

 The activation of NK1R is subject to short-term dynamics acting on
 SP release & activation:
 "Release of substance P is induced by stressful stimuli, and 
 the magnitude of its release is proportional to the intensity 
 and frequency of stimulation. More potent and more frequent 
 stimuli allow diffusion of substance P farther from the site 
 of release, allowing activation of an approximately 3- to 5-times
 greater number of NK1 receptor-expressing neurons" (Mantyh, 2002)
 We used the short-term plasticity equations from Fuhrmann et al, 
 2002, to model the frequency dependent nature of the NK1R activation.
 Note however that these equations are being used in a more broad
 manner (i.e. the dynamics are intended to encapsulate not only
 synaptic release but other mechanisms, such as SP diffusion, which
 can give rise to facilitation).
 
 Written by Paulo Aguiar, May 2009
 pauloaguiar@fc.up.pt

 My thanks to Ted Carnevale for his fantastic help at
 The NEURON Forum
 "
ENDCOMMENT


NEURON {
	POINT_PROCESS NK1_DynSyn	
	USEION ca WRITE ica	
	RANGE tau_rise, tau_decay
	RANGE U1, tau_rec, tau_fac
	RANGE i, g, e, iNK1R, ica, ca_ratio
	NONSPECIFIC_CURRENT iNK1R
}
    
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}    
    
PARAMETER {
	tau_rise  = 10.0		(ms)  : dual-exponential conductance profile
	tau_decay = 5000.0		(ms)  : IMPORTANT: tau_rise < tau_decay
	U1        = 1.0			(1)   : The parameter U1, tau_rec and tau_fac define _
	tau_rec   = 0.1			(ms)  : the pre-synaptic SP short-term plasticity _
	tau_fac   = 0.1			(ms)  : mechanism (see Fuhrmann et al, 2002)
	e         = 0.0			(mV)  : reversal potential
	ca_ratio = 0.1			(1)   : the increase in cai is assumed here to be proportional to iNK1; this is the constant of proportionality	  
}
    
    
ASSIGNED {
	v		(mV)
	i		(nA)
	g		(umho)
	factor	(1)
	ica		(nA)
	iNK1R	(nA)
}

STATE {
	A
	B
}

INITIAL{
	LOCAL tp
	A = 0
	B = 0
	tp = (tau_rise*tau_decay)/(tau_decay-tau_rise)*log(tau_decay/tau_rise)
	factor = -exp(-tp/tau_rise)+exp(-tp/tau_decay)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B-A
	i = g*(v-e)
	ica = ca_ratio * i
	: the ica current is "silenced" by removing ica from the total current
	iNK1R = i - ica
	: this way we have an increase in cai without a contribution of ica in membrane currents
}

DERIVATIVE state{
	A' = -A/tau_rise
	B' = -B/tau_decay
}

NET_RECEIVE (weight, Pv, P, Use, t0 (ms)){
	INITIAL{
		P=1
		Use=0
		t0=t
	}	

	Use = Use * exp(-(t-t0)/tau_fac)
	Use = Use + U1*(1-Use) 
	P = 1-(1- P) * exp(-(t-t0)/tau_rec)
	Pv= Use * P
	P = P - Use * P
	
	t0=t
	
	A=A + weight*factor*Pv
	B=B + weight*factor*Pv
}

