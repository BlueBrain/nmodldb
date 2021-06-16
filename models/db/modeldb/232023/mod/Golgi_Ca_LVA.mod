TITLE Low threshold calcium current Cerebellum Golgi Cell Model
:
:   Ca++ current responsible for low threshold spikes (LTS)
:   RETICULAR THALAMUS
:   Differential equations
:
:   Model of Huguenard & McCormick, J Neurophysiol 68: 1373-1383, 1992.
:   The kinetics is described by standard equations (NOT GHK)
:   using a m2h format, according to the voltage-clamp data
:   (whole cell patch clamp) of Huguenard & Prince, J Neurosci.
:   12: 3804-3817, 1992.  The model was introduced in Destexhe et al.
:   J. Neurophysiology 72: 803-818, 1994.
:   See http://www.cnl.salk.edu/~alain , http://cns.fmed.ulaval.ca
:
:    - Kinetics adapted to fit the T-channel of reticular neuron
:    - Q10 changed to 5 and 3
:    - Time constant tau_h fitted from experimental data
:    - shift parameter for screening charge
:
:   ACTIVATION FUNCTIONS FROM EXPERIMENTS (NO CORRECTION)
:
:   Reversal potential taken from Nernst Equation
:
:   Written by Alain Destexhe, Salk Institute, Sept 18, 1992
:---
:Adapted by Sungho Hong and Claus Lang
:Computational Neuroscience Unit, Okinawa Institute of Science and Technology, Japan
:Supervision: Erik De Schutter
:
:Correspondence: Sungho Hong (shhong@oist.jp)
:
:September 16, 2017


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
        SUFFIX Golgi_Ca_LVA
        USEION ca2 READ ca2i, ca2o WRITE ica2 VALENCE 2
	RANGE Q10_diff,Q10_channel,gbar_Q10, ic, fix_celsius
        RANGE g, gbar, m_inf, tau_m, h_inf, tau_h, shift
	RANGE ica2, m ,h, ca2rev
	RANGE phi_m, phi_h
}

UNITS {
        (molar) = (1/liter)
        (mV) =  (millivolt)
        (mA) =  (milliamp)
        (mM) =  (millimolar)

        FARADAY = (faraday) (coulomb)
        R = (k-mole) (joule/degC)
}

PARAMETER {
        v               (mV)
        fix_celsius = 37 (degC)
        eca2 (mV)
	gbar  = 2.5e-4 (mho/cm2)
        shift   = 2     (mV)            : screening charge for Ca_o = 2 mM
        ca2i  (mM)           : adjusted for eca=120 mV
        ca2o  (mM)
	Q10_diff	= 1.5
	Q10_channel	= 3
}

STATE {
        m h
}

ASSIGNED {
        ica2     (mA/cm2)
        ca2rev   (mV)
	g        (mho/cm2)
        m_inf
        tau_m   (ms)
        h_inf
        tau_h   (ms)
        phi_m
        phi_h
	gbar_Q10 (mho/cm2)
	ic
}

BREAKPOINT {
        SOLVE ca2state METHOD cnexp
        ca2rev = (1e3) * (R*(fix_celsius+273.15))/(2*FARADAY) * log (ca2o/ca2i)
        g = gbar_Q10 * m*m*h
        ica2 = gbar_Q10 * m*m*h * (v-ca2rev)
	ic = ica2
}

DERIVATIVE ca2state {
        evaluate_fct(v)

        m' = (m_inf - m) / tau_m
        h' = (h_inf - h) / tau_h
}

UNITSOFF
INITIAL {
:
:   Activation functions and kinetics were obtained from
:   Huguenard & Prince, and were at 23-25 deg.
:   Transformation to 36 deg assuming Q10 of 5 and 3 for m and h
:   (as in Coulter et al., J Physiol 414: 587, 1989)
:
	gbar_Q10 = gbar*(Q10_diff^((fix_celsius-23)/10))
        evaluate_fct(v)
        m = m_inf
        h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) {
:
:   Time constants were obtained from J. Huguenard
:
	TABLE m_inf, tau_m, h_inf, tau_h
	DEPEND fix_celsius, shift, phi_m, phi_h FROM -100 TO 30 WITH 13000
        phi_m = 5.0 ^ ((fix_celsius-24)/10)
        phi_h = Q10_channel ^ ((fix_celsius-24)/10)

        m_inf = 1.0 / ( 1 + exp(-(v+shift+50)/7.4) )
        h_inf = 1.0 / ( 1 + exp((v+shift+78)/5.0) )

        tau_m = ( 3 + 1.0 / ( exp((v+shift+25)/10) + exp(-(v+shift+100)/15) ) ) / phi_m
        tau_h = ( 85 + 1.0 / ( exp((v+shift+46)/4) + exp(-(v+shift+405)/50) ) ) / phi_h
}
UNITSON
