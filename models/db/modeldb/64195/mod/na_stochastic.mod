
COMMENT

Author: Stefan Hallermann 

Provides Na-channel stochastics as described in Kole et al. (2006).

In the initiation the number of channels (N) is calculated for each segment and stored as a
range variable. For each dt the procedure noise() is evaluated once. In noise() each channel
in each state has the chance to move to one of its neighbouring states with the appropriate
probability. After this "update" of the number of channels in each state the resulting current
through the open channels (i) is calculated depending on the local driving force in the
segment. The eight-state reaction model was adopted from Hille (1978) and the resulting
kinetics are identical to Na-kinetics of Mainen and Sejnowski (1996). 

Caution 1: Because of the fast kinetics of Na-channels the iteration time dt should be 1 micro
second.
Caution 2: You should wait at the beginning of the simulation until the Na channels have
reached steady state (i.e. distributed in a voltage dependent manner over the eight states of
the reaction scheme). 

ENDCOMMENT

NEURON {
	SUFFIX na
	USEION na READ ena WRITE ina
	RANGE gbar,m0h0,m1h0,m2h0,m3h0,m0h1,m1h1,m2h1,m3h1
	GLOBAL tha, thi1, thi2, qa, qi, qinf, thinf
	RANGE minf, hinf, mtau, htau
	GLOBAL Ra, Rb, Rd, Rg
	GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

PARAMETER {
	gbar = 1000   	(pS/um2)	: 0.12 mho/cm2
	vshift = -10	(mV)		: voltage shift (affects all)
								
	tha  = -35	(mV)		: v 1/2 for act		(-42)
	qa   = 9	(mV)		: act slope		
	Ra   = 0.182	(/ms)		: open (v)		
	Rb   = 0.124	(/ms)		: close (v)		

	thi1  = -50	(mV)		: v 1/2 for inact 	
	thi2  = -75	(mV)		: v 1/2 for inact 	
	qi   = 5	(mV)	        : inact tau slope
	thinf  = -65	(mV)		: inact inf slope	
	qinf  = 6.2	(mV)		: inact inf slope
	Rg   = 0.0091	(/ms)		: inact (v)	
	Rd   = 0.024	(/ms)		: inact recov (v) 

	temp = 23	(degC)		: original temp 
	q10  = 2.3			: temperature sensitivity

	v 		(mV)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)

	gamma=10e-12		(S)		:single channel cond
	seed
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	gna		(pS/um2)
	ena		(mV)

	minf 		hinf		:not used
	mtau (ms)	htau (ms)	:not used
	tadj				:not used

	dt			(ms)
	area			(um2)
   	N			:number of channels
	m0h0			:inactivated state
	m1h0			:inactivated state
	m2h0			:inactivated state
	m3h0			:inactivated state
	m0h1			:closed state
	m1h1			:closed state
	m2h1			:closed state
	m3h1			:open state
}
 

STATE { m h }		:not used

FUNCTION trap0(v,th,a,q) {
	if (fabs((v-th)/th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}

INITIAL { 
	m = 0
	h = 0
	N=abs(((1e-8*area*1e-4*gbar)/gamma)+0.5)				:area in um2; 1e-8*area in cm2; gbar in pS/um2=S/m2; 1e-4*gbar in S/cm2; gamma in S

	m0h0=0
	m1h0=0
	m2h0=0
	m3h0=0
	m0h1=N		:therefore you should wait at the beginning of the simulation until the Na channels have reached steady state.
	m1h1=0
	m2h1=0
	m3h1=0

	set_seed(seed)
}

BREAKPOINT {
        SOLVE states METHOD cnexp
	ina = tadj*((m3h1*gamma)/(1e-8*area))*(v-ena)		:cond/cm2 * delta_pot		(cond=N_open*gamma_ih in S)
} 


DERIVATIVE states {   
        m' =  m
        h' =  h
	noise()
}

PROCEDURE noise() {
	LOCAL a,b,vm,p,m0h0_merk,m1h0_merk,m2h0_merk,m3h0_merk,m0h1_merk,m1h1_merk,m2h1_merk,m3h1_merk,am,bm,ah,bh

	vm=v+vshift
	a = trap0(vm,tha,Ra,qa)
	b = trap0(-vm,-tha,Rb,qa)
        tadj = q10^((celsius - temp)/10)
	mtau = 1/tadj/(a+b)
	minf = a/(a+b)
	am=tadj*a
	bm=tadj*b


		:"h" inactivation 
	a = trap0(vm,thi1,Rd,qi)
	b = trap0(-vm,-thi2,Rg,qi)
	htau = 1/tadj/(a+b)
	hinf = 1/(1+exp((vm-thinf)/qinf))	:

:	since hinf is not a/(a+b) but the above value,
:	now htau and hinf are used to define again the rates ah and bh
:   	in adittion change nomenclature of "h" inactivation: a to b and b to a
	ah=hinf/htau
	bh=(1-hinf)/htau

	m0h0_merk=m0h0
	m1h0_merk=m1h0
	m2h0_merk=m2h0
	m3h0_merk=m3h0
	m0h1_merk=m0h1
	m1h1_merk=m1h1
	m2h1_merk=m2h1
	m3h1_merk=m3h1

:scop_random gives random number uniform between 0 and 1

	:m0h0
	p=1-exp(-dt*(ah+3*am))				
	FROM ii=1 TO m0h0_merk {
		if (scop_random()<= p)	{					:probability that a channel in the state m0h0 goes to state m0h1 or m1h0
			if (scop_random()<= ah/(ah+3*am))	{		:probability that this channel goes to state m0h1 (via rate ah)

				m0h0=m0h0-1
				m0h1=m0h1+1
			}
			else {							:otherwise this channel goes to m1h0 (via rate 3*am)
				m0h0=m0h0-1
				m1h0=m1h0+1
			}
		}
	}

	:m1h0
	p=1-exp(-dt*(ah+2*am+bm))
	FROM ii=1 TO m1h0_merk {
		if (scop_random()<= p) {	
			if (scop_random()<= ah/(ah+2*am+bm)) {	
				m1h0=m1h0-1
				m1h1=m1h1+1
			}
			else {							:"if you dont like m1h1 (via rate ah), you have to choose again between the two states that are left (m2h0 and m0h0)"
				if (scop_random()<= 2*am/(2*am+bm)) {	
					m1h0=m1h0-1
					m2h0=m2h0+1
				}
				else {
					m1h0=m1h0-1
					m0h0=m0h0+1
				}
			}
		}
	}

	:m2h0
	p=1-exp(-dt*(ah+am+2*bm)) 
	FROM ii=1 TO m2h0_merk {
		if (scop_random()<= p){	
			if (scop_random()<= ah/(ah+am+2*bm))	{
				m2h0=m2h0-1
				m2h1=m2h1+1
			}
			else {
				if (scop_random()<= am/(am+2*bm))	{
					m2h0=m2h0-1
					m3h0=m3h0+1
				}
				else {
					m2h0=m2h0-1
					m1h0=m1h0+1
				}
			}
		}
	}

	:m3h0
	p=1-exp(-dt*(ah+3*bm)) 
	FROM ii=1 TO m3h0_merk {
		if (scop_random()<= p){	
			if (scop_random()<= ah/(ah+3*bm))	{	
				m3h0=m3h0-1
				m3h1=m3h1+1
			}
			else {	
				m3h0=m3h0-1
				m2h0=m2h0+1
			}
		}
	}

:h1------------
	:m0h1
	p=1-exp(-dt*(bh+3*am)) 
	FROM ii=1 TO m0h1_merk {
		if (scop_random()<= p){	
			if (scop_random()<= bh/(bh+3*am))	{	
				m0h1=m0h1-1
				m0h0=m0h0+1
			}
			else {
				m0h1=m0h1-1
				m1h1=m1h1+1
			}
		}
	}

	:m1h1
	p=1-exp(-dt*(bh+2*am+bm)) 
	FROM ii=1 TO m1h1_merk {
		if (scop_random()<= p){	
			if (scop_random()<= bh/(bh+2*am+bm))	{	
				m1h1=m1h1-1
				m1h0=m1h0+1
			}
			else {
				if (scop_random()<= 2*am/(2*am+bm))	{	
					m1h1=m1h1-1
					m2h1=m2h1+1
				}
				else {
					m1h1=m1h1-1
					m0h1=m0h1+1
				}
			}
		}
	}

	:m2h1
	p=1-exp(-dt*(bh+am+2*bm)) 
	FROM ii=1 TO m2h1_merk {
		if (scop_random()<= p){	
			if (scop_random()<= bh/(bh+am+2*bm))	{	
				m2h1=m2h1-1
				m2h0=m2h0+1
			}
			else {
				if (scop_random()<= am/(am+2*bm))	{	
					m2h1=m2h1-1
					m3h1=m3h1+1
				}
				else {
					m2h1=m2h1-1
					m1h1=m1h1+1
				}
			}
		}
	}

	:m3h1
	p=1-exp(-dt*(bh+3*bm)) 
	FROM ii=1 TO m3h1_merk {
		if (scop_random()<= p){	
			if (scop_random()<= bh/(bh+3*bm))	{	
				m3h1=m3h1-1
				m3h0=m3h0+1
			}	
			else {
				m3h1=m3h1-1
				m2h1=m2h1+1
			}
		}
	}
}
