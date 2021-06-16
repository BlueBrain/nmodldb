TITLE klt.mod  The low threshold conductance of cochlear nucleus neurons

COMMENT

NEURON implementation of Jason Rothman's measurements of VCN conductances.

This file implements the low threshold potassium current found in several brainstem
 nuclei of the auditory system, including the spherical and globular bushy cells
  (Manis and Marx, 1991; Rothman and Manis, 2003a,b) and octopus cells (Bal and
  Oertel, 2000) of the ventral cochlear nucleus, principal cells of the medial 
  nucleus of the trapzoid body (Brew and Forsythe, 1995, Wang and Kaczmarek, 
  1997) and neurons of the medial superior olive. The current is likely mediated by 
  heteromultimers of Kv1.1 and Kv1.2 potassium channel subunits. The specific 
  implementation is described in Rothman and Manis, J. Neurophysiol. 2003, in the 
  appendix. Measurements were made from isolated neurons from adult guinea pig, 
  under reasonably stringent voltage clamp conditions. The measured current is 
  sensitive to the mamba snake toxin dendrotoxin-I.


Similar conductrances are found in the homologous neurons of the avian auditory 
system (Reyes and Rubel; Zhang and Trussell; Rathouz and Trussell), and the 
conductance described here, in the absence of more detailed kinetic measurements
, is probably suitable for use in modeling that system.


Original implementation by Paul B. Manis, April (JHU) and Sept, (UNC)1999.

File split implementation, February 28, 2004.

Contact: pmanis@med.unc.edu

Modified and simplified by Christian Roessert in 2010

ENDCOMMENT


UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX klt
        USEION k READ ek WRITE ik
        RANGE gbar, g, ik
        GLOBAL winf, wtau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        ek = -90 (mV)
        gbar = 0 (mho/cm2) <0,1e9>
        tfac = 0.5
}

STATE {
        w
}

ASSIGNED {
    ik (mA/cm2) 
    g (mho/cm2)
    winf 
    wtau (ms)
    }

BREAKPOINT {
        SOLVE states METHOD cnexp 

        g = gbar*(w^4)
        ik = g*(v - ek)
}
 
 
INITIAL {
    rates(v)
    w = winf
}

DERIVATIVE states {  
        rates(v)
        w' =  (winf-w)/wtau
}


UNITSOFF

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

    winf = (1 / (1 + exp(-(v + 48) / 6)))^0.25
    wtau =  tfac*((100 / (6*exp((v+60) / 6) + 16*exp(-(v+60) / 45))) + 1.5)
}

UNITSON
