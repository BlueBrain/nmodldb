TITLE TRPM8 current for bladder small DRG neuron soma model
: Adapted from Olivares et al., 2015

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)


: based on Voets 2004 with some modifications to resemble data from Malkia 2007
: Written by Orio, P. & Olivares E.  - December 2014 and edited for DRG neuron model 2017
COMMENT 
NONSPECIFIC_CURRENT made ica
DVimf = minf
DV = m
Dvmin = mmin
DVmax = mmax
taudv=taum
gm8=gbar
im8=ica

dE=dE
vhalf= vhalf
C=C
z=z=0.65 ??
em8=0

--Removed--
cam8 removed: calcium conc due to trpm8
tauca removed
d
accel: for accelerating channel rate
n
ENDCOMMENT

NEURON {
    SUFFIX trpm8
	USEION ca READ ica, cai WRITE ica VALENCE 2
    RANGE minf, vhalf
    RANGE p_ca, accel
    RANGE em8, gbar, am8, C, z
	NONSPECIFIC_CURRENT im8
}


UNITS {
    R = (k-mole) (joule/degC)
    (mA) = (milliamp)
    (mV) = (millivolt)
    (mol) = (1)
    (molar) = (1/liter)
    (mM) = (millimolar)
} 

CONSTANT {
    F = 96500        (coulomb)        : moles do not appear in units
}


PARAMETER {
	gbar = 1e-7        (mho/cm2)
    dE	 = 9e3           (joule)
    C	 = 67
    z	 = 0.65

    em8  = 0            (mV)
    mmin = 0		   (mV)
    mmax = 200		   (mV)
    Kca  = 0.0005       (mM)

    p_ca = 0.01    
    taum = 80000       (ms)
}

STATE {
    m  (mV)
}

INITIAL {
	rate(cai)
	m= minf
}

ASSIGNED {
    celsius (degC)
    v       (mV)
    ica     (mA/cm2)
    vhalf   (mV)
    am8
    minf    (mV)
	cai		(mM)
	im8     (mA/cm2)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
	vhalf=(1000)*(C*R*celsius - dE)/(z*F)+m : temperature function + calcium function
	am8=1/(1+exp(-z*F*(v-vhalf)/((1000)*R*(celsius+273.15)))) 
    im8 = (1-p_ca)*gbar*am8*(v-em8) 
	ica = p_ca*gbar*am8*(v-em8) : p_ca used was 0.01 which is the fraction of current by calcium. See Olivares et al., 2015 pg 5 and mod file here:
	: https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=182988&file=/OlivaresEtAl2015/Neuron/trpm8.mod#tabs-2. 
}

DERIVATIVE states {
	rate(cai)
	m' = (minf-m)/taum
}

PROCEDURE rate(ca (mM)){
	minf = mmin+(mmax-mmin)*(ca)/(Kca+ca) : calcium function
}