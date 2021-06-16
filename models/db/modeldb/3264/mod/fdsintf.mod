COMMENT
Implementation of the model of short-term facilitation and depression described in
  Varela, J.A., Sen, K., Gibson, J., Fost, J., Abbott, L.R., and Nelson, S.B.
  A quantitative description of short-term plasticity at excitatory synapses 
  in layer 2/3 of rat primary visual cortex
  Journal of Neuroscience 17:7926-7940, 1997
This is an "integrate and fire" cell that can receive multiple streams of 
synaptic input via NetCon objects.  Each stream keeps track of its own 
weight and activation history.

An input that occurs when the cell is not refractory causes an immediate 
increment of activation by an amount equal to the product of the intrinsic 
synaptic weight and the facilitation and depression terms.  The latter are 
functions of the time that has elapsed since the most recent prior event.
The effect of synaptic input on the facilitation and depression terms is 
computed _after_ the activation is incremented.

An input that occurs when the cell is refractory does nothing to the activation
but does affect facilitation and depression.

The printf() statements are for testing purposes only.
Added FUNCTION yy() to allow smooth plot of "activation"
and FUNCTIONs Fval(), D1val(), and D2val() to allow plots of 
F, D1, and D2 vs. t

20110810 Changed POINT_PROCESS to ARTIFICIAL_CELL, 
to be consistent with best practices in recent versions of NEURON.
ENDCOMMENT


NEURON {
:	POINT_PROCESS FDSIntFire
	ARTIFICIAL_CELL FDSIntFire
	RANGE tau, y, spikedur, refrac, f, tau_F, d1, tau_D1, d2, tau_D2
	: y plays the role of voltage
}

PARAMETER {
	tau = 10 (ms) < 1e-9, 1e9 >
	spikedur = 1 (ms) < 0, 1e9 >
	refrac = 5 (ms) < 0, 1e9 >
	: these values are from Fig.3 in Varela et al. 1997
	: the (1) is needed for the range limits to be effective
	f = 0.917 (1) < 0, 1e9 >	: facilitation
	tau_F = 94 (ms) < 1e-9, 1e9 >
	d1 = 0.416 (1) < 0, 1 >	: fast depression
	tau_D1 = 380 (ms) < 1e-9, 1e9 >
	d2 = 0.975 (1) < 0, 1 >	: slow depression
	tau_D2 = 9200 (ms) < 1e-9, 1e9 >
}

ASSIGNED {
	y
	t0 (ms)
	refractory
}

INITIAL {
	y = 0
	t0 = t
	refractory = 0 : 0-integrates input, 1-spiking, 2-refractory
}

FUNCTION yy() {
	yy = y*exp(-(t - t0)/tau)
}

FUNCTION Fval(F, time) {
	Fval = 1 + (F-1)*exp(-time/tau_F)
}

FUNCTION D1val(D1, time) {
	D1val = 1 - (1-D1)*exp(-time/tau_D1)
}

FUNCTION D2val(D2, time) {
	D2val = 1 - (1-D2)*exp(-time/tau_D2)
}

NET_RECEIVE (w, F, D1, D2, tsyn (ms)) {
INITIAL {
: these are in NET_RECEIVE to be per-stream
	F = 1
	D1 = 1
	D2 = 1
	tsyn = t
: this header will appear once per stream
: printf("t\t t-tsyn\t F\t D1\t D2\t amp\t newF\t newD1\t newD2\n")
}

	if (flag == 0) {
		: calculate present facilitation & depression of the active input
		: from the values they had the last time the input was active
:		F = 1 + (F-1)*exp(-(t - tsyn)/tau_F)
:		D1 = 1 - (1-D1)*exp(-(t - tsyn)/tau_D1)
:		D2 = 1 - (1-D2)*exp(-(t - tsyn)/tau_D2)
		F = Fval(F, t-tsyn)
		D1 = D1val(D1, t-tsyn)
		D2 = D2val(D2, t-tsyn)
: printf("%g\t%g\t%g\t%g\t%g\t%g", t, t-tsyn, F, D1, D2, w*F*D1*D2)
		tsyn = t
	}
	if (refractory == 0) { : inputs affect activation only when excitable
:		y = y*exp(-(t - t0)/tau)
		y = yy()
		y = y + w*F*D1*D2
		t0 = t
		if (y > 1) {
			refractory = 1
			y = 2
			net_send(spikedur, refractory)
			net_event(t)
		}
	} else if (flag == 1) { : don't really need this.
		: this extra event just makes the spike narrow
		refractory = 2
		y = -1
		net_send(refrac, refractory)
	} else if (flag == 2) { : ready to integrate again
		refractory = 0
		y = 0
	}
	if (flag == 0) {
		: add facil and depr produced by the active input
		F = F + f
		D1 = D1 * d1
		D2 = D2 * d2
: printf("\t%g\t%g\t%g\n", F, D1, D2)
	}
}
