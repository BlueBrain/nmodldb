TITLE Cerebellum Golgi Cell Model

COMMENT

Author:L. Forti & S. Solinas
Data from: Santoro et al. J Neurosci. 2000
Last revised: April 2006

---
Adapted by Sungho Hong and Claus Lang
Computational Neuroscience Unit, Okinawa Institute of Science and Technology, Japan
Supervision: Erik De Schutter

Correspondence: Sungho Hong (shhong@oist.jp)

September 16, 2017
ENDCOMMENT

NEURON {
        SUFFIX Golgi_hcn1
	NONSPECIFIC_CURRENT ih
	RANGE Q10_diff,Q10_channel,gbar_Q10, ic, fix_celsius
	RANGE o_fast_inf, o_slow_inf, tau_f, tau_s, Erev
	RANGE gbar,r,g, o
}

UNITS {
        (mA) = (milliamp)
	(mV) = (millivolt)
	(S)  = (siemens)
}


PARAMETER {
    fix_celsius = 37 (degC)
	gbar = 5e-5   (S/cm2)
        Erev = -20 (mV)
	Q10_diff	= 1.5
	Q10_channel	= 3

	Ehalf = -72.49 (mV)
	c = 0.11305	(/mV)

	rA = 0.002096 (/mV)
        rB = 0.97596  (1)
        tCf = 0.01371 (1)
        tDf = -3.368  (mV)
	tEf = 2.302585092 (/mV)
	tCs = 0.01451 (1)
        tDs = -4.056  (mV)
	tEs = 2.302585092 (/mV)
}

ASSIGNED {
	ih		(mA/cm2)
        v               (mV)
	g		(S/cm2)
	o_fast_inf
        o_slow_inf
        tau_f           (ms)
	tau_s           (ms)
	gbar_Q10 (mho/cm2)
	Q10 (1)
	ic
	o
    }


INITIAL {
	gbar_Q10 = gbar*(Q10_diff^((fix_celsius-23)/10))
	rate(v)
	o_fast = o_fast_inf
	o_slow = o_slow_inf

}

STATE {	o_fast o_slow }		: fraction of fast and slow open channels


BREAKPOINT {
	SOLVE state METHOD cnexp
	g = gbar_Q10 * (o_fast + o_slow)
        ih = g * (v - Erev)
	ic = ih
	o = o_fast + o_slow
}

DERIVATIVE state {
	rate(v)
	o_fast' = (o_fast_inf - o_fast) / tau_f
	o_slow' = (o_slow_inf - o_slow) / tau_s
}


FUNCTION r(potential (mV))  { 	:fraction of fast component in double exponential
	UNITSOFF
	r =  rA * potential + rB
        UNITSON
}

FUNCTION tau(potential (mV),t1,t2,t3, Q10_channel) (ms) {
	UNITSOFF
	Q10 = Q10_channel^((fix_celsius -23(degC)) / 10(degC))
        tau = exp(((t1 * potential) - t2)*t3) :/Q10 : Dividing by Q10  commented to avoid segmentation violation
	UNITSON
}

FUNCTION o_inf(potential (mV), Ehalf, c)  {
	UNITSOFF
        o_inf = 1 / (1 + exp((potential - Ehalf) * c))
        UNITSON
}

PROCEDURE rate(v (mV)) {
	TABLE o_fast_inf, o_slow_inf, tau_f, tau_s
	DEPEND fix_celsius FROM -100 TO 30 WITH 13000

	: r(v) is the fraction of fast component in double exponential
	o_fast_inf = r(v) * o_inf(v,Ehalf,c)
	o_slow_inf = (1 - r(v)) * o_inf(v,Ehalf,c)

	tau_f =  tau(v,tCf,tDf,tEf,Q10)
	tau_s =  tau(v,tCs,tDs,tEs,Q10)
}
