TITLE Two-closed one-open state sodium channel model 
: G Baranauskas November 2005, this model was developed by fitting
: voltage clamp data reported in J. Neurosci 26 pp 671-84 2006

NEURON {
     SUFFIX MCna1
     USEION na READ ena WRITE ina
     RANGE gna1bar, gna1, ina   
     GLOBAL cnt1, cnt2, Na_intern, Na_extern 
}

UNITS {
     (mA) = (milliamp)
     (mV) = (millivolt)
}

PARAMETER {

     gna1bar=.120 (mho/cm2) <0,1e9> : though this parameter has the
: same units as conductance because of the Goldman-Hodgkin-Katz
: formalism used here the real meaning of this value is
: P*[Na_extern]*F*F/RT. The conversion factor from this parameter to the
: conductance will depend on voltage and Na_inter and Na_extern defined
: in the BREAKPOINT fomula

     Na_intern = 15 (mM)
     Na_extern = 135 (mM)

}
STATE {
     C1 C2 O I1 I2 IO
}
ASSIGNED {
     v (mV)
     celsius (degC) : 6.3 
     ena (mV)
     cnt1 cnt2
     ina (mA/cm2)
     gna1 (mho/cm2)
}

ASSIGNED {  a1 (/ms)   b1 (/ms)  a2 (/ms)  b2 (/ms) a3 (/ms)  b3 (/ms)}

INITIAL {
	cnt1 = 0
	cnt2 = 0
     C1=1
     rate(v*1(/mV))
     SOLVE states STEADYSTATE sparse
}

BREAKPOINT {
     SOLVE states METHOD sparse
     gna1 = gna1bar*O

     ina = gna1*v*(Na_intern/Na_extern - exp(-v/25.4))/(1-exp(-v/25.4)) 
: the Goldman-Hodgkin-Katz equation was used to determine current 
: amplitude, in this case 15 mM Na_inter and 135 mM Na_extern was used

	cnt1 = cnt1 + 1
}

KINETIC states {
	cnt2 = cnt2 + 1
     rate(v*1(/mV))
:     CONSERVE I1 + I2 + IO + C1 + C2 + O = 1
     ~ I1 <-> I2 (a1, b1) 
     ~ C1 <-> C2 (a1, b1)         
   
     ~ I2 <-> IO (a2, b2)
     ~ C2 <-> O (a2, b2)

     ~ C1 <-> I1 (a3, b3)
     ~ C2 <-> I2 (a3, b3)      
     ~ O <-> IO (a3, b3)   
}

UNITSOFF

PROCEDURE rate(v) {LOCAL q10, q11 
     TABLE a1, a2, a3, b1, b2, b3 DEPEND celsius FROM -100 TO 100 WITH 200
     q10 = (2.8)^((celsius - 13)/10)  : the fit was performed at 13C
     q11 = (2.4)^((celsius - 13)/10) 
 
     a1 = q10*10*exp((v+6)/45) : a1 and b1 correspond to the fast
: transition and a2, b2 - the slow transition that determines kinetics

     a2 = q10*11/(0.4+exp(-(v+6)/12))
     b1 = q10*0.35*exp(-(v+6)/8)
     b2 = q10*0.035/(0.0015 + exp((v+6)/12))

     a3 = q11*2/(2+exp(-(v+6)/12))      : a3 and b3 are inactivation rates
     b3 = q11*0.00005*exp(-(v+6)/13)
}
UNITSON






