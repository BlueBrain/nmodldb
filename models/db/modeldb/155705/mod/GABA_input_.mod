TITLE first_GABA_recept_model (first order kinetics )
COMMENT
This program just explores basics charactieristics of GABAergic synapses including facilitation
according to Synchrony generation in recurrent networks with frequency-dependent synapses
Journal of Neuroscience 20:RC50:1-5, 2000. Originally implemented by Carnevale T.

addapted to be used as GABA_A synapse O. J. Avella
	
ENDCOMMENT

NEURON {
	POINT_PROCESS Gaba_synf
	RANGE e, i, g, area_cell :           => g not in use
	RANGE tau_d, frac_rec
	RANGE tau_1, tau_rec, tau_facil, U, u0
	NONSPECIFIC_CURRENT i
}

UNITS { 
	(nA)   = (nanoamp)
	(mV)   = (millivolt)
	(uS)   = (microsiemens)
	(um)   = (micron)
	}

PARAMETER {	
	:==============================================================================================
	:					unused
	tau_d = 10 (ms)      < 1e-9, 1e9 > : tau_d decay time,
	frac_rec = 0.9 		(1)  <0,1> :frac_rec usage fraction of receptors
	area_cell= 1 		(um2) 	   :cell surface area
	g=1 (1)
	:==============================================================================================
	

	:facilitation (used)
	
	e=-80	(mV)			 : value of boergers and koppel's paper
	tau_1 = 1 (ms) < 1e-9, 1e9 >     : tau_1 was the same for inhibitory and excitatory synapses in the models used by T et al
	tau_rec = 100 (ms) < 1e-9, 1e9 > : tau_rec = 100 ms for inhibitory synapses/800 ms for excitatory,
	tau_facil =0: 1000 (ms) < 0, 1e9 > : tau_facil = 1000 ms for inhibitory synapses/0 ms for excitatory
	u0 = 0 (1) < 0, 1 > 		 : initial value for the "facilitation variable"
	U = 0.04 (1) < 0, 1 > 		 : U = 0.04/0.5 i/e syn; the (1) needed for < 0, 1 > ...
					 :to be effective in limiting the values of U and u0
	:==============================================================================================


}

ASSIGNED {
	v (mV)
	i (nA)
	x
	
}

STATE {
geff (siemens)	
}

INITIAL {	
geff=0	 
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i=geff*(v - e)

}

DERIVATIVE state {
	geff'=-geff/tau_1
}

NET_RECEIVE(w (us),y, z, u,tp(ms)) {:tp time of previous spike
INITIAL {
: these are in NET_RECEIVE to be per-stream
	y = 0
	z = 0
	u = u0 : or u = 0 for no-facilitation
	tp = t
}
	
	: first calculate z at event-    based on prior y and z 
	z = z*exp(-(t - tp)/tau_rec)
	z = z + ( y*(exp(-(t - tp)/tau_1) - exp(-(t - tp)/tau_rec)) / ((tau_1/tau_rec)-1) )
	: now calc y at event-
	y = y*exp(-(t - tp)/tau_1)
	x = 1-y-z
	:(I DON'T UNDERSTAND how Carnevale arrived to the former solution; do you Ronald?, could you explain me? )

	: calc u at event--
	if (tau_facil > 0) {u = u*exp(-(t - tp)/tau_facil)} else {u = U}
	if (tau_facil > 0) {u= u + U*(1-u)}			:updates facil. factor @ event arrival
	geff=geff + w*x:*u					:updates conductance @ event arrival
	y= y + x:*u						:updates number of active "resources @ event arrival"
	tp = t






}
