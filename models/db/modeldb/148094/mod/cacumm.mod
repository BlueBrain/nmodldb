COMMENT
	calcium accumulation into a volume of area*depth next to the
	membrane with a decay (time constant tau) to resting level
	given by the global calcium variable cai0_ca_ion
	Modified to include a resting current (irest) and peak value
	(cmax)
	i is a dummy current needed to force a BREAKPOINT
	Note: irest cancels out resting calcium current contributions
	by mechanisms that write ica.  To initialize irest properly
	use a custom proc init that assigns values to irest as illustrated
	in this excerpt:
	finitialize(Vrest) // use v_init if want to change.
        fcurrent()
	forall if (ismembrane("cacum")) {
	       for(x,0) irest_cacum(x)=ica(x) // (fixed from =-ica minus sign bug)
	}

ENDCOMMENT

NEURON {
	SUFFIX cacum
	USEION ca READ ica WRITE cai
	NONSPECIFIC_CURRENT i
	RANGE depth, tau, cai0, cmax, irest, ca_tmax
}

UNITS {
        (um) = (micron)
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F = (faraday) (coulombs)
}

PARAMETER {
	depth = 0.1 (um)	: assume volume = area*depth
	irest = 0  (mA/cm2)		: to be initialized in hoc	
	tau = 100 (ms) :  ~40 ms phenomenologically fits dendrite shafts, ~350 ms fits spines (murthy et al 2000 PNAS, earlier model default was 100 (ms)
	cai0 = 50e-6 (mM)	: Requires explicit use in INITIAL
			: block for it to take precedence over cai0_ca_ion
			: Do not forget to initialize in hoc if different
			: from this default.
}

ASSIGNED {
        ica      (mA/cm2)
        cmax     (milli/liter)
        ca_tmax  (ms)
        i        (mA/cm2)
}

STATE {
	cai (mM)
}

INITIAL {
	cai = cai0
:	irest = ica : this make simulations depend on end result of prior simulations
	cmax=cai
	ca_tmax=0
}

BREAKPOINT {
	SOLVE integrate METHOD derivimplicit
	if (cai>cmax) {cmax=cai ca_tmax=t}
	if (cai<0) {cai=0}
	i=0
}

DERIVATIVE integrate {
	cai' = -(ica+irest)/depth/F/2 * (1e4) + (cai0 - cai)/tau
}
