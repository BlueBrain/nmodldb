TITLE stochastic Ih-channels

COMMENT

Author: Stefan Hallermann 

Provides Ih-channel stochastics as described in Kole et al. (2006).

In the initiation the number of channels (N) and the number of open channels (N_open)
is calculated for each segment and stored as a range variable. For each dt the Procedure
noise() is evaluated once. In noise() each open channel has the chance to close according
to the closing rate beta(v) and each close channels has the chance to open according to
the opening rate alpha(v). After this "update" of the number of open channels the
resulting current through the open channels (i) is calculated depending on the local
driving force in the segment.
	
ENDCOMMENT


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 			(mV)
	ehd=-45  		(mV) 		:ih-reversal potential			       
	ghdbar=0.00015 		(S/cm2)		:default Ih conductance; exponential distribution is set in Ri18init.hoc 
	gamma=680e-15		(S)		:single channel cond
	seed
}


NEURON {
	SUFFIX ih
	NONSPECIFIC_CURRENT i
	RANGE ghdbar,N,N_open
}

STATE {
	l
}

ASSIGNED {
	i (mA/cm2)
	dt			(ms)
	area			(um2)
	N			:number of channels
	N_open			:number of open channels
}

INITIAL {								:calculates the number of Ih-channel per segment and the number of open channels at the initial potential
	N=abs(((1e-8*area*ghdbar)/gamma)+0.5)				:area in um2; 1e-8*area in cm2; ghdbar in S/cm2; gamma in S
	N_open=abs(N*alpha(v)/(alpha(v)+beta(v)))
	l=0								:only needed for dummy diff.eq.

	set_seed(seed)
}


BREAKPOINT {
	SOLVE states METHOD cnexp					:only needed to make the proc noise() be evaluated once per dt (breakpoint() is evaluated twice per dt!)
	i = ((N_open*gamma)/(1e-8*area))*(v-ehd)			:cond/cm2 * delta_pot		(cond=N_open*gamma in S)
}


FUNCTION alpha(v(mV)) {
	:opening rate in 1/s
	alpha = 6.43*(v+154.9)/(exp((v+154.9)/11.9)-1)			:parameters are estimated by direct fitting of HH model to activation time constants and voltage activation curve recorded at 34C
}

FUNCTION beta(v(mV)) {
	:closing rate in 1/s
	beta = 193*exp(v/33.1)			
}


DERIVATIVE states {     						
	l' =  l			    					:dummy diff.eq.
	noise()
}

PROCEDURE noise() {
	LOCAL N_close,N_open_merk,N_close_merk,a,b,prob_open,prob_close
	a=alpha(v)
	b=beta(v)
	N_open_merk=N_open
	N_close_merk=N-N_open

	:activation	(all close channels have the chance to open with the probability prob_open depending on dt and the opening rate alpha(v))
	:(the approximation with the 1. and 2. element of exp. infinitive series is only little faster:  prob_open=dt*b/1000)
	prob_open=1-exp(-dt*a/1000)					:(/1000 since dt is in ms and rate beta in 1/s)
	FROM ii=1 TO N_close_merk {
		if (scop_random()<= prob_open)	{			:scop_random uniform between 0 and 1
			N_open=N_open+1
		}
	}

	:deactivation	(all open channels have the chance to close with the probability prob_close depending on dt and the closing rate beta(v))
	prob_close=1-exp(-dt*b/1000)	
	FROM ii=1 TO N_open_merk {
		if (scop_random()<= prob_close)	{		
			N_open=N_open-1
		}
	}
}




