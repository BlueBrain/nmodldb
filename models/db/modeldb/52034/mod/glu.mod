COMMENT
Modified from Zach's code: 6.3.94
Simulates alpha function synapse with synaptic delay.
If pre > PreThresh then a new spike is put on the queue at t+delay.
A delta pulse is added into the alpha function equations when time >
a time on the queue. A maximum of (pspMAX) can be handled on the queue.
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE pspMAX 15				:max # of spikes on queue

NEURON {
	POINT_PROCESS glu
	POINTER pre				:presyn variable to trigger on
	RANGE gmax, i, g, delay
	RANGE delta, flag, erev:, p
	NONSPECIFIC_CURRENT i
	GLOBAL tau,PreThresh
}

ASSIGNED {
	pspSTACK[pspMAX]
	psp0
	psp1
	pspN
	pre
	flag
	delta
:	rn
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	gmax=0.005     (umho)
	erev=0		(mV)
	tau=1.0 			:tau for alpha function
	delay=1.5	(ms)		:synaptic delay
	PreThresh=-20	(mV)		:presyn spike thresh to test for
	v		(mV)
        dt		(ms)
:	p = 0.2
}

STATE {
	X
	Y
}

INITIAL { LOCAL i
	X = 0				:~proportion of bound transmitter/receptors
	Y = 0				:~proportion of channels open (g)
	flag = 0			:set if spike is above thresh
	delta = 0			:set if quantum released
	g = 0
	pspN = 0
	psp0 = 0
	psp1 = 0
:	rn = 1/(2^31)
  FROM i=0 TO pspMAX {
  pspSTACK[i] = 0 
 }
}

ASSIGNED { i (nA)  g (umho)}

BREAKPOINT {			:everything in BREAKPOINT called twice except SOLVE

  SOLVE dstates METHOD cnexp

  g = gmax * Y / (tau * 0.36787944)   : tau * exp(-1)
  i = g*(v - erev)
}

DERIVATIVE dstates {
	CheckThresh()        
	CheckTime()
	X' = delta - X / tau
        Y' = X - Y / tau
}

PROCEDURE CheckThresh() { 

  if (flag) { if (pre < PreThresh) {flag = 0} }
  else {
    if (pre > PreThresh) {
      flag = 1
:      if (random()*rn < p) {
        AddSpike()
:      }
    }
  }
}

PROCEDURE AddSpike() { 

  if (pspN < pspMAX) {
    pspN=pspN+1
    if (psp1 >= pspMAX) { psp1 = 0 }
    pspSTACK[psp1] = t+delay-dt		:time at which next quantum released
    psp1=psp1+1
  }
  else { printf("ERROR: spike queue full - synaptic delay too long?\n") }
}

COMMENT
PROCEDURE CheckTimeo() {		:hacked code to get around pseudo bug
  if ( pspSTACK[psp0] > 0) {		:if delta set for only 1 dt, doesn't work
    if (t-pspSTACK[psp0] > dt) {
  	delta = 0
        pspN=pspN-1
        pspSTACK[psp0] = 0
        psp0 = psp0+1
        if (psp0 >= pspMAX) { psp0 = 0 }
    } else if (t-pspSTACK[psp0] > 0) {
        printf("quantum released\n")
        delta = 1/dt
    }
  }
}
ENDCOMMENT

PROCEDURE CheckTime() {			
delta = 0
  if ( (pspSTACK[psp0] > 0) && (t > pspSTACK[psp0]) ){		
        delta = 1/dt
        pspN=pspN-1
        pspSTACK[psp0] = 0
        psp0 = psp0+1
        if (psp0 >= pspMAX) { psp0 = 0 }
  }
}
