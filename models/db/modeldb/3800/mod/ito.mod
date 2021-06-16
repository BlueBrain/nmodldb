TITLE Cardiac Ito  current
: Hodgkin-Huxley type k channel from Courtemanche et al Am J Physiol Am J Physiol 1998 275:H301


NEURON {
	SUFFIX Ito
	USEION k READ ek WRITE ik
	RANGE gto, ik, Tauact, Tauinact
	GLOBAL minf, ninf, mtau, ntau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)
	
}

PARAMETER {
	 gto=0.3304e-3 (S/cm2) <0,1e9>
	Tauact=1 (ms)
	Tauinact=1 (ms)
}
STATE {
	 m n
}

ASSIGNED {
	v (mV)
	celsius (degC) : 37
	ik (mA/cm2)
	minf ninf
	mtau (ms)
        ntau (ms)
	ek (mV)
}
LOCAL k
INITIAL {
	rates(v)
	m = minf
        n = ninf
	
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	ik = gto*m*m*m*n*(v - ek)
}

DERIVATIVE states {	: exact when v held constant
	rates(v)
	m' = (minf - m)/mtau
        n' = (ninf - n)/ntau
}

UNITSOFF
FUNCTION alp(v(mV),i) { LOCAL q10 : order m n
	v = v
	q10 = 2.2^((celsius - 37)/10)
       if (i==0) {
	          alp = q10*0.65/(exp(-(v + 10)/8.5) + exp(-(v - 30)/59))
          } else if (i==1) {
                   alp = q10/(18.53 + exp((v + 113.7)/10.95))
          }
}

FUNCTION bet(v(mV),i) (/ms) { LOCAL q10 : order m n 
	v = v 
	q10 = 2.2^((celsius - 37)/10)
        if (i==0){
	         bet = q10*0.65/(2.5 + exp((v + 82)/17))
        }else if (i==1){
                  bet = q10/(35.56 + exp(-(v + 1.26)/7.44))
        }
}
                
FUNCTION ce(v(mV),i)(/ms) {  :  order m n 
        v = v
        
        if (i==0) {
                ce = 1/(1 + exp(-(v + 20.47)/17.54))
        }else if (i==1){
                ce = 1/(1 + exp((v + 43.1)/5.3))
        }
}


PROCEDURE rates(v) {LOCAL a,b,c :rest = -70
	:TABLE minf, ninf, mtau, ntau DEPEND celsius FROM -100 TO 100 WITH 200
		a = alp(v,0)  b=bet(v,0) c = ce(v,0)
		mtau = 1/(a + b)/3*Tauact
		minf = c
               a = alp(v,1)  b=bet(v,1) c = ce(v,1)
		ntau = 1/(a + b)/3*Tauinact
		ninf = c
}
UNITSON
