
UNITS 
{
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
NEURON {
        SUFFIX KM
	USEION k READ ek WRITE ik
        RANGE gMbar, gM, timesTau, plusTau
        GLOBAL uinf, utau
}
 
PARAMETER 
{
        gMbar = 0.00034 (S/cm2)	<0,1e9>
        eK = -95 (mV)
	ek (mV)
	timesTau=1
	plusTau=0
}
 

STATE 
{
        u 
}
 
ASSIGNED 
{
        v (mV)
        celsius (degC)
	gM (S/cm2)
        ik (mA/cm2)
        uinf
	utau (ms)
}
 
LOCAL uexp
 
BREAKPOINT 
{
        SOLVE states METHOD cnexp
        gM = gMbar*u*u
	ik = gM*(v - eK)
}
 
 
INITIAL 
{
	rates(v)
	u = uinf
}

DERIVATIVE states 
{  
        rates(v)
        u' =  (uinf-u)/utau
}
 
LOCAL q10


PROCEDURE rates(v(mV))   
                     
{
        LOCAL  alpha, beta, sum

UNITSOFF

        TABLE uinf, utau FROM -150 TO 150 WITH 3000 : Mitti
        
        q10 = 5^((celsius - 23)/10)
               
        alpha = 0.016 / exp((v+52.7)/-23)
        beta =  0.016 / exp((v+52.7)/18.8)
        sum = alpha + beta
	
        uinf = alpha/sum
	utau = timesTau/(sum*q10)+plusTau
}
 

 
UNITSON

