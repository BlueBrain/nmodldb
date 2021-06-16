COMMENT
Na channel from Schweighofer et al 1999.
The referred model is an inferior olive neuron

B. Torben-Nielsen @ HUJI, 21-10-2010
ENDCOMMENT

NEURON {
       SUFFIX iona
       USEION na READ ena WRITE ina
       RANGE gbar,g,i,minf,hinf,tauh,m : now i can access these variables
}

UNITS {
      (S) = (siemens)
      (mS) = (millisiemens)
      (mV) = (millivolt)
      (mA) = (milliamp)
}

PARAMETER {
	  gbar = 70 (mS/cm2)
	  ena = 55 (mV)
}

ASSIGNED {
	 v (mV)
	 ina (mA/cm2)
	 i (mA/cm2)
	 g (mS/cm2)
	 minf 
	 hinf 
	 tauh (ms)
}

STATE {
      h
}

INITIAL {
	rates(v)
	h = hinf
}

BREAKPOINT {
	   rates(v)
	   SOLVE states METHOD cnexp
	   g = gbar *minf*minf*minf*h
	   i = g * (v - ena)*(0.001)
	   ina = i
}

DERIVATIVE states {
	   h' = (hinf -h)/tauh
}

PROCEDURE rates(v (mV)) {
	  LOCAL a_h, b_h, a_m,b_m
	  UNITSOFF
	  if(fabs(v+41.0) < 1e-6) {
	    a_m=(0.1*(v+41.000001)) / ( 1-exp( -(v+41.000001)/10 ) )
	  } else {
	    a_m=(0.1*(v+41)) / ( 1-exp( -(v+41)/10 ) )
	  }
	  b_m=9.0*exp( -(v+66)/20 )
	  minf=a_m/(a_m+b_m)
	  
	  a_h=5.0*exp( -(v+60)/15 )
	  b_h=(v+50) / ( 1-exp( -(v+50)/10 ) )
	  hinf=a_h/(a_h+b_h)
	  tauh=250/( a_h+b_h ) : was 170, 250
	  UNITSON
}