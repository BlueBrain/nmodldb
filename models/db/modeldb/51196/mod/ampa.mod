TITLE Jonas et al. (1993) kinetic scheme of AMPA receptors 
COMMENT
Author: Elena Saftenku, 2001
ENDCOMMENT
NEURON{
POINT_PROCESS GrC_AMPA
POINTER pNp, pNd, pglu
RANGE   iampa
RANGE gampa, gbarampa, freqdel
RANGE  Erev 
NONSPECIFIC_CURRENT iampa
}

UNITS{
(nA)=(nanoamp)
(mA)=(milliamp)
(mV)=(millivolt)
(umho)=(micromho)
(molar)=(1/liter)
(mM)=(millimolar)
}

PARAMETER {
Erev = 0 (mV)
gbarampa =  1.915e-3 (umho): gpeak=0.45e-3
cam_4_3 = 26.6  (/mM-ms)
cam_3_2 = 13.3 (/mM-ms)
cam_7_6 =2.41 (/mM-ms)
cam_1_5 = 0.109 (/ms)
cam_1_2 = 0.302 (/ms)
cam_2_1 = 4.2 (/ms)
cam_2_3 = 12.5 (/ms)
cam_2_6 = 0.395 (/ms)
cam_3_4 =  6.24 (/ms)
cam_3_7 = 0.513 (/ms)
cam_5_1=0.0334 (/ms)
cam_5_6 = 0.13 (/ms)
cam_6_2 = 0.000546 (/ms)
cam_6_5 = 0.00815 (/ms)
cam_6_7 =0.0571(/ms)
cam_7_3=0.0281(/ms)
freqdel=0.579 (/nA): determines the part of long-term change of AMPAR-mediated synaptic current
}

ASSIGNED{
v (mV)
iampa (nA)
gampa  (umho)
pNp (nA) : reference to protein autophosphorylation Np
pNd (nA) : reference to protein autophosphorylation Nd
pglu (mM): reference to glutamate concentration
celsius (degC)
}

STATE{
St1ampa
St2ampa
St3ampa
St4ampa
St5ampa
St6ampa
St7ampa
}

INITIAL {
St1ampa=0
St2ampa=0
St3ampa=0
St4ampa=1
St5ampa=0
St6ampa=0
St7ampa=0
}

BREAKPOINT
{
SOLVE states METHOD derivimplicit
gampa = gbarampa*St1ampa *(1+freqdel*(pNp-pNd))
iampa =gampa*(v-Erev)
}

DERIVATIVE states {
LOCAL Q10
Q10 = 2^((celsius-20(degC))/10(degC))
  
St1ampa' = Q10*(-(cam_1_5 + cam_1_2)*St1ampa + cam_2_1*St2ampa + 
cam_5_1* St5ampa)
St2ampa' = Q10*( -(cam_2_1 + cam_2_3  + cam_2_6) * St2ampa + cam_1_2*St1ampa+
cam_3_2*pglu*St3ampa + cam_6_2*St6ampa)
St3ampa' = Q10*(-(cam_3_2*pglu+cam_3_4+cam_3_7) * St3ampa + 
cam_4_3 *pglu* St4ampa + cam_7_3*St7ampa + cam_2_3* St2ampa)
St4ampa' = Q10*(-cam_4_3*pglu*St4ampa+cam_3_4*St3ampa)
St5ampa' = Q10*(-(cam_5_1+cam_5_6)*St5ampa+ cam_1_5*St1ampa+cam_6_5*St6ampa)
St6ampa' = Q10*(-(cam_6_2+cam_6_5+cam_6_7)*St6ampa+cam_7_6*pglu*St7ampa
+cam_5_6*St5ampa + cam_2_6*St2ampa)
St7ampa' = Q10*(-(cam_7_6*pglu+cam_7_3)*St7ampa+cam_3_7*St3ampa+
cam_6_7*St6ampa)
}

