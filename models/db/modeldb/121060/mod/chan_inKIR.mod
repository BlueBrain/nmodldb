TITLE inKIR channel
COMMENT
Described by Steephen and Manchanda 2009, which was based on Mermelstein et al. 1998 and Wolf et al. 2005

Steephen, J. E., & Manchanda, R. (2009). Differences in biophysical properties of nucleus accumbens 
medium spiny neurons emerging from inactivation of inward rectifying potassium currents. J Comput Neurosci, 
doi:10.1007/s10827-009-0161-7

Mermelstein, P. G., Song, W. J., Tkatch, T., Yan, Z., & Surmeier, D. J. (1998). Inwardly rectifying potassium
(IRK) currents are correlated with IRK subunit expression in rat nucleus accumbens medium spiny neurons. 
The Journal of Neuroscience, 18, 6650–6661.

Wolf, J. A., Moyer, J. T., Lazarewicz, M. T., Contreras, D., Benoit-Marand, M., O’Donnel, P., et al. (2005). 
NMDA/AMPA ratio impacts state transitions and entrainment to oscillations in a computational model of 
the nucleus accumbens medium spiny projection neuron. The Journal of Neuroscience, 25, 9080–9095.
doi:10.1523/JNEUROSCI.2220-05.2005.

ENDCOMMENT

NEURON {
	SUFFIX inKIR
	USEION  k READ ek WRITE ik
	RANGE  g, ik
	GLOBAL minf, mtau, hinf, htau, eff_hinf, a
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

CONSTANT {	
	Q10 = 3 (1)
}

PARAMETER {
	ek			(mV)
	gmax = 1.4e-4	(mho/cm2)	<0,1e9>
	m_vh = -82	(mV)	: half activation
	m_ve = 13		(mV)	: slope
	a = 0.47 (1)		: 0.47 default, 0.27 for pinKir
:	a = 0.27 (1)		
}

ASSIGNED {
	v	(mV)
	g	(mho/cm2)
	ik	(mA/cm2)
	minf	(1)
	mtau	(ms)
	hinf (1)
	htau (ms)
	eff_hinf (1)
	qt (1)
}

STATE {
	m h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g=gmax*m*(a*h+(1 -a)) 
	ik = g*(v - ek)
}

INITIAL {
	qt = Q10^((celsius-35)/10)
	rates(v)
	m = minf
	h = hinf
}

DERIVATIVE states { 
	rates(v)
	m' = (minf-m)/mtau
	h' = (hinf-h)/htau
}

FUNCTION_TABLE tabmtau(v(mV)) (ms)
FUNCTION_TABLE tabhtau(v(mV)) (ms)
FUNCTION_TABLE tabhinf(v(mV))(1)

: rates() computes rate and other constants at present v
: call once from hoc to initialize inf at resting v

PROCEDURE rates(v(mV)) {
	mtau = tabmtau(v)/qt
:	mtau = tabmtau(v)
	htau = tabhtau(v)/qt
:	htau = tabhtau(v)
	hinf = tabhinf(v)
	eff_hinf = a*hinf+(1 -a)
	minf = 1/(1 + exp((v - m_vh)/m_ve))
}

