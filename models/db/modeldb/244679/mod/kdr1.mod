TITLE Delayed rectifire
 
COMMENT
HVA-VGKC (HVAK). Starting from the formulation of the “delayed rectifier channel” given in De
Schutter and Bower, 1994, The HVAK kinetic parameters were modified to account for the
behaviours at depolarised states. Specifically, the activation curve was multiplied by a sigmoid
function to track the occurrence of the first Ca 2+ spike. Then the activation time was decreased by
~95% to reproduce the number and shape of the observed Ca 2+ spikes. Notably, this is the only
channel for which the experimental pharmacological block was not available.

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
        SUFFIX Kdr1
	USEION k WRITE ik
        RANGE  gkbar, gk, minf, hinf, mexp, hexp, ik, alpha, beta ,pdr,vsh1,vsh2,hdr,pp
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gkbar	= .6 (mho/cm2)
        ek	= -85 (mV)
        pp=6
        pdr=1
        hdr=10
        vsh1 (mV)
        vsh2 (mV)

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
        gk = gkbar *m*m*h
	ik = gk* (v-ek)
}
 
UNITSOFF
 
INITIAL {
	rates(v,vsh1,vsh2,pdr,hdr,pp)
	m = minf
	h = hinf
}

PROCEDURE states() {  :Computes state variables m,h
        rates(v,vsh1,vsh2,pdr,hdr,pp)      :             at the current v and dt.
        m = m + mexp*(minf-m)
	h = h + hexp*(hinf-h)
}

PROCEDURE rates(v(mV),vsh1 (mV),vsh2 (mV),pdr,hdr,pp) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  q10, tinc, tauh, alpha, beta, gamma, zeta, taum, vshi1,vshi2,corr

        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
                :"m" potassium activation system
                vshi1=vsh1+v
                if(vshi1 == 8) {vshi1 = 8.0001}
                     
        alpha = -0.0047*(vshi1-8)/(exp((vshi1-8)/(-12))-0.9999)
	beta = exp((vshi1+127)/(-30))
     corr=1/(1+exp(-0.4*(v+35)))
	minf = (alpha/(alpha+beta))*corr
    gamma = -0.0047*(vshi1+12)/(exp((vshi1+12)/(-12))-0.9999)
    zeta = exp((vshi1+147)/(-30))
    taum = pdr/(gamma + zeta)
    mexp = 1 - exp(pp*tinc/taum)






                :"h" potassium activation system
                vshi2=vsh2+v
        hinf = 1.0 / (1+exp((vshi2+25)/4))
        if(v<-25) {
	tauh = 1200/hdr
	}else{
	tauh = 10/hdr
	}
	hexp = 1 - exp(tinc/tauh)
       
}



 
UNITSON

