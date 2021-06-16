: K-Strom, I_CT,  Stacey, Durand 2000
: eK from Martina

UNITS 
{
        (molar) = (1/liter)
	(mM) = (millimolar)
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
NEURON {
        SUFFIX KCT
	USEION ca READ cai
	USEION k WRITE ik
        RANGE gCTbar, gCT
        GLOBAL cinf, dinf, dtau, ctau
}
 
PARAMETER 
{
        gCTbar = 0.120 (S/cm2)	<0,1e9>
        eK = -95 (mV)
	ctau = 0.55 (ms)
}
 

STATE 
{
        c d
}
 
ASSIGNED 
{
        ik (mA/cm2)
        cai (mM)
        v (mV)
        celsius (degC)
	gCT (S/cm2)
	cinf
	dinf
	dtau (ms)
:	ctau (ms)
}
 

BREAKPOINT 
{
        SOLVE states METHOD cnexp
        gCT = gCTbar*c*c*d
	ik = gCT*(v - eK)
}
 
 
INITIAL 
{
	rates(v)
	c = cinf
	d = dinf
}

DERIVATIVE states 
{  
        rates(v)
 
        c' =  (cinf-c)/ctau
	d' = (dinf-d)/dtau
}
 
LOCAL q10


PROCEDURE rates(v(mV))   :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
{
        LOCAL  alpha, beta, sum, vshift

UNITSOFF
               
        vshift = 40 * log10(cai)
        q10 = 3^((celsius - 6.3)/10)
                :"c" potassium activation system

        alpha = -0.0077 * vtrap(v+vshift+103, -12)
        beta =  1.7 / exp((v+vshift+237)/30)
        sum = alpha + beta
	: see Warman
	:	ctau = 1/sum
        cinf = alpha/sum
	: "d" activation
	alpha = 1/(exp((v+79)/10))
	beta = 4/(exp((v-82)/-27)+1)
	sum = alpha + beta
	dinf = alpha/sum
	dtau = 1/sum 
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
 
UNITSON
