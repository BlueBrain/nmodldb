TITLE A current
 
COMMENT
A-type VGKC. Starting from the formulation of A-type VGKCs channels in De Schutter and Bower,
1994, the kinetic parameters were modified in line with modifications of T-type VGCCs to account for
behaviours at hyperpolarised states. The density was corrected at intermediate states to account
for partial inactivation.

    Current Model Reference: Karima Ait Ouares , Luiza Filipis , Alexandra Tzilivaki , Panayiota Poirazi , Marco Canepari (2018) Two distinct sets of Ca 2+ and K + channels 
are activated at different membrane potential by the climbing fibre synaptic potential in Purkinje neuron dendrites. 

PubMed link: 

Contact: Filipis Luiza (luiza.filipis@univ-grenoble-alpes.fr)
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX KA1
	USEION k WRITE ik
        RANGE  gkbar, gk, minf, hinf, mexp, hexp, ik, vh1, vm1, vm2, t1,t2
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gkbar	= .015 (mho/cm2)
        ek	= -85 (mV)
	mon = 1
	hon = 1
    vh1=85 (mV)
    vm1 (mV)
    vm2 (mV)
    t1=0.3
    t2=50
}
 
STATE {
        m h
}
 
ASSIGNED {
        ik (mA/cm2)
        gk minf hinf mexp hexp 
}
 
BREAKPOINT {
        SOLVE states
        gk = gkbar *m*m*m* m*h 
	ik = gk* (v-ek)
}
 
UNITSOFF
 
INITIAL {
	rates(v,vm1,vm2)
	m = minf
	h = hinf
}

PROCEDURE states() {  :Computes state variables m, h
        rates(v,vm1,vm2)      :             at the current v and dt.
        m = mon * (m + mexp*(minf-m))
        h = hon * (h + hexp*(hinf-h))
}
 
PROCEDURE rates(v, vm1,vm2) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  q10, tinc, alpha, beta, sum
        
        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
                :"m" potassium activation system
        alpha = 1.4/(1+exp((v+vm1)/(-12)))
        beta =  0.49/(1+exp((v+vm2)/4))
        sum = alpha + beta
        minf = (alpha/sum)*(1/(1+exp(-t1*(v+t2))))
        mexp = 1 - exp(tinc*sum)
                :"h" potassium inactivation system
        alpha = 0.0175/(1+exp((v+vh1)/8))
        beta = 1.3/(1+exp((v+13)/(-10)))
        sum = (alpha + beta)
        hinf = alpha/sum
       
        hexp = 1 - exp(tinc*sum)
}

 
UNITSON

