TITLE IzhikevichCurrent

COMMENT

The module implements Izhikevich's model as a  non-specific current. Module is defined as point-process 
to be able send and receive net events.

The module treats Izhikevich's system of equations as
v'=(Iizh(v)+i)/cm
where Iizh(v)=e*v^2+f*v+g-u is a transmembrane current. u is a state variable

Do not forget setup:
 cm in 1uF
 L in 1um
 diam in 1/PI

Here an example, who to use it in hoc file.

objref izh
soma{
	L=1
	diam=1/PI
	nseg=1
	izh = new izhcur(0.5)
	cm=1
}

Parameters are given by Izhikevich matlab code:

       a        b       c      d       I
================================================================================
      0.02      0.2     -65     6      14       % tonic spiking
      0.02      0.25    -65     6       0.5     % phasic spiking
      0.02      0.2     -50     2      15       % tonic bursting
      0.02      0.25    -55     0.05    0.6     % phasic bursting
      0.02      0.2     -55     4      10       % mixed mode
      0.01      0.2     -65     8      30       % spike frequency adaptation
      0.02     -0.1     -55     6       0       % Class 1
      0.2       0.26    -65     0       0       % Class 2
      0.02      0.2     -65     6       7       % spike latency
      0.05      0.26    -60     0       0       % subthreshold oscillations
      0.1       0.26    -60    -1       0       % resonator
      0.02     -0.1     -55     6       0       % integrator
      0.03      0.25    -60     4       0       % rebound spike
      0.03      0.25    -52     0       0       % rebound burst
      0.03      0.25    -60     4       0       % threshold variability
      1         1.5     -60     0     -65       % bistability
      1         0.2     -60   -21       0       % DAP
      0.02      1       -55     4       0       % accomodation
     -0.02     -1       -60     8      80       % inhibition-induced spiking
     -0.026    -1       -45     0      80       % inhibition-induced bursting    

written by Ruben Tikidji-Hamburyan <rth@nisms.krinc.ru>, 2013 - 1015
ENDCOMMENT

NEURON {
	POINT_PROCESS izhcur
	NONSPECIFIC_CURRENT i_izh
	RANGE a,b,c,d,e,f,g,uinit,F,u,I
}

PARAMETER {
	a = 0.01	(1)
	b = 0.2		(1)
	c = -65		(mV)
	d = 2		(1)
	e = 0.04	(1)
	f = 5		(1)
	g = 140		(1)
	uinit = -14	(1)
	F = 1		(1) : speedup / slowdown parameter
	I = 0		(mA/cm2) : doesn't work yet
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

ASSIGNED {
	v (mV)
	i_izh (mA)
}

STATE { u }

BREAKPOINT {
	SOLVE states METHOD cnexp
	i_izh = (-1e-5)*F*(vinfi(v)-u)-I	:minus, because it is inward current	
}

INITIAL {
	u = uinit
	net_send(0,1)					:we have send first event
}

DERIVATIVE states {
	UNITSOFF
	u'= F*a*(b*v-u)
	UNITSON
}

FUNCTION vinfi(v (mV)) {
	UNITSOFF
	vinfi = e*v*v + f*v + g
	UNITSON
}

NET_RECEIVE (w) {
	if (flag == 1) {
		WATCH (v > 30.0) 2
	} else {
		net_event(t)
		v = c
		u = u+d
  }
}
