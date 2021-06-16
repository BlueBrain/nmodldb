COMMENT
alpha function synapse implemented as continuously integrated
kinetic scheme a la Srinivasan and Chiel (Neural Computation) so that
one can give many stimuli which summate.

Onset times are generated from exponentially decay distribution.

Conductance located in state variable G
The amplitude of each individual alpha function is given by stim,
stim * t * exp(-t/tau).
ENDCOMMENT

DEFINE SIZE 1000

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}


NEURON {
	POINT_PROCESS SynAlphaPoisson
	RANGE tau, stim, e, i,onset, offset, mean
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	tau=.25 (ms)
	stim=50 (umho)
	e=0	(mV)
	v	(mV)
	mean  = 1000   (ms)
	onset  =0      (ms)
	offset =1e9    (ms)
}

ASSIGNED {
	index
	i (nA)
	bath (umho)
	k (/ms)
	t_activation (ms)
	interval (ms)
}

STATE {
	A (umho)
	G (umho)
}

INITIAL {
	k = 1/tau
	A = 0
	G = 0
	t_activation=t-1
	
	while( t > t_activation){
	   t_activation=onset-mean+mean*exprand(1)
	}
	net_send(t_activation,1)
}

? current
BREAKPOINT {
	SOLVE state METHOD sparse
	i = G*(v - e)
}

: at each onset time a fixed quantity of material is added to state A
: this material moves through G with the form of an alpha function



? kinetics
KINETIC state {
	~ A <-> G	(k, 0)
	~ G <-> bath	(k, 0)
}

NET_RECEIVE(w){
    if(flag==1){
        A = A + stim
        interval=mean*exprand(1)
        t_activation=t_activation+interval
        if(t_activation<offset){
            net_send(interval,1)
        }
    }
}























