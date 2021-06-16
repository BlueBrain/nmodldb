TITLE Basic Potassium current
 
COMMENT
 from "Gamma Oscillation by Synaptic Inhibition in a Hippocampal Interneuronal Network Model" (Wang and Buzsaki 1996)
 Used in Role of a Striatal Slowly Inactivating Potassion Current in Short-term Facilitation of Corticostriatal Inputs" A computer Simulation Study" (Mahon et al. 2000)
Implemented by Kevin M. Biddell kevin.biddell@gmail.com
7/12/06

NOTE: 1S=1mho Neuron wants the units in mhos not millisiemens, please note the conversion!

Phi =5 and no q10 or temp adjustment according to Bruno Delord 11/13/06

ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
 	SUFFIX Km
	USEION k WRITE ik
	RANGE gkmbar, gkm, ninf, an, Bn
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
  	
	ek	= -90	(mV)
	gkmbar	= 0.006 (mho/cm2) : 6mS
	phi	= 5  < 0, 1e9 > : according to Bruno delord 11/13/06
      Van	= -27 :NOT the original value from wang and Buzsaki
	Kan	= 1
	Vbn	= -37 :NOT the original value from wang and Buzsaki
	Kbn	= 80
}
 
STATE {
        n
}
 
ASSIGNED {
        v  (mV)
	ik (mA/cm2)
	celsius		(degC)
 	ninf
	an
	Bn
        gkm
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gkm = gkmbar*n^4
        ik = gkm*(v - ek)
  
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	n= ninf
}

DERIVATIVE states {  :Computes states variable n 
        rates(v)      :             at the current v and dt.
       
	n'=phi*(an*(1-n)-Bn*n)

}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        
        an = (-0.01*(v-Van)/Kan/(exp(-0.1*(v-Van)/Kan)-1))
        Bn = 0.125*exp(-(v-Vbn)/Kbn)
        ninf = an/(an+Bn)
	      
}
 
UNITSON

