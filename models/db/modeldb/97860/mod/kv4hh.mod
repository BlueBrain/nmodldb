: KV4hh.MOD
:
: Josh Held	j-held@northwestern.edu
: 3/2003
:
: Kv4 current (A current) model
:
: Current is defined by
:	i = g * m^3 * h * (v-e)
:
: Data from 9_13_2_1


NEURON {
	SUFFIX kv4hh
	USEION k READ ek WRITE ik
    RANGE minf, tm, hinf, thf, ths, ik, alpha, th
    RANGE gbar
    GLOBAL vhm, vcm
    GLOBAL vhh, vch, p
    GLOBAL vhtm, atm, btm, Ctm, tm0
    GLOBAL vhthf, athf, bthf, Cthf, thf0
    GLOBAL vhths, aths, bths, Cths, ths0
    GLOBAL vha, vca, a0
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	ek		        (mV)

    gbar   = 1 	(mho/cm2)

    vhm     = -59
    vcm     = -21

    vhh     = -86.4
    vch     = 9.6
    p       = 0.02

    vhtm    = -55
    atm     = 24.8
    btm     = 6.4
    Ctm     = 2.75
    tm0     = 1.21

    vhthf   = -70
    athf    = 15
    bthf    = 15
    Cthf    = 20
    thf0    = 20

    vhths   = -70
    aths    = 15
    bths    = 15
    Cths    = 150
    ths0    = 150

    vha     = -30
    vca     = 15
    a0      = 0.36




}

STATE {
	m h1 h2
}

ASSIGNED {
	v		        (mV)
	ik		(mA/cm2)
	minf
	tm		(ms)
	hinf
	thf		(ms)
    ths     (ms)
    alpha
}

BREAKPOINT {
	SOLVE states METHOD cnexp
 	ik = gbar * m^3 * (alpha*h1 + (1-alpha)*h2) * (v-ek)
}


DERIVATIVE states{
	rates(v)
	m' = (minf - m)/tm
    h1' = (hinf - h1)/thf
    h2' = (hinf - h2)/ths
}

UNITSOFF

INITIAL {
	rates(v)
	m = minf
    h1 = hinf
    h2 = hinf
}

PROCEDURE rates(v(mV)) {LOCAL q10
	UNITSOFF
    q10 = 3^((celsius-22)/10)
    minf = 1/(1 + exp((v-vhm)/vcm))^3
    hinf = p + ((1-p)/(1 + exp((v-vhh)/vch)))
    tm = (1/q10)*(tm0 + Ctm/(exp((v-vhtm)/atm) + exp(-(v-vhtm)/btm)))
    thf = (1/q10)*(thf0 + Cthf/(exp((v-vhthf)/athf) + exp(-(v-vhthf)/bthf)))
    ths = (1/q10)*(ths0 + Cths/(exp((v-vhths)/aths) + exp(-(v-vhths)/bths)))

    alpha = a0 + ((1-a0)/(1 + exp((v-vha)/vca)))


}

UNITSON
