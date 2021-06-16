TITLE Contraction
COMMENT
	Sarcomere Dynamics  modelisation from NEGRONI and LASCANO J Mol Cell Cardiol 1996 28, 915
    modified for Neuron by FE GANNIER & CO MALECOT
	francois.gannier@univ-tours.fr (University de TOURS)
ENDCOMMENT
INCLUDE "Unit.inc"
NEURON {
	SUFFIX Cont
	USEION ca READ cai WRITE cai, ica VALENCE 2
	RANGE Force, Fb, Fp, Tr, K, A
	RANGE Qm, Qrel, Qpump, QpumpR, L
	RANGE Qd, Qa, Qb, Qr, Qd1, Qd2, period
}

ASSIGNED {
   ica	 (mA/cm2) 
   Qa    (uM/s)
   Qb    (uM/s)
   Qr    (uM/s)
   Qrel  (uM/s)
   Qpump (uM/s)
   Qd    (uM/s)
   Qd1   (uM/s)
   Qd2   (uM/s)
   Tr	(uM)
   Fb	(mN/mm2)
   Fp	(mN/mm2)
   Force	(mN/mm2)
 }

PARAMETER {
   Y1     = 39      (/uM/s)
   Z1     = 30      (1/s)
   Y2     = 1.3     (1/s)
   Z2     = 1.3     (1/s)
   Y3     = 30      (1/s)
   Z3     = 1560    (/uM/s)
   Y4     = 40      (1/s)
   Yd     = 9       (s/um2)
   Tt     = 70      (uM)
   B      = 1200    (1/s)
   hc     = 0.005   (um)
   La     = 1.17    (um)
   Ra     = 20      (1/um2)
   Kp     = 150     (uM/s)
   Km     = 0.1     (uM)
   Qm     = 1600    (uM/s)
   t1     = 25   	(ms)
   period = 1000   	(ms)
   decal  			(ms)
   
   A	= 1800		(mN/mm2/um/uM)
   K	= 140000	(mN/mm2/um5)
   Lo	= 0.97		(um)
   
   QpumpR = 12.25		(uM/s)
   L	= 1.05 	     (um)
   
}

STATE {
   X	 (um)
   TCa   (uM)
   TCaA  (uM)
   TA    (uM)
   cai   (mM)
}

LOCAL decal
INITIAL {
	VERBATIM
		cai = _ion_cai;
	ENDVERBATIM
	TCa = 0
   TA = 0
   TCaA = 0
   X = L
	decal = 0
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
	Tr = Tt - TCa - TCaA - TA
 	Fb = A * (TCaA+TA) * (L - X)
	Fp = K * (L - Lo)^5
	Force = Fb + Fp
	ica = 0
}

DERIVATIVE state {
   X' = (0.001)*B*(L-X-hc)
   Qd = Y4 * TA
   Qd1 = (1e+6)*(Yd * (X')^2 * TA)
   Qd2 = (1e+6)*(Yd * (X')^2 * TCaA)   
   Qa = Y2 * TCa * exp(-Ra*(L - La)^2) - Z2 * TCaA
   Qb = (1000)*(Y1 * cai * Tr) - (Z1 * TCa)
   Qr = Y3*TCaA - (1000)*(Z3*TA*cai)
   Qpump = ( Kp /(1 + ( Km / ((1000)*cai))^2))

   if (t > decal+period) { 	
		decal = decal + period
	}
   Qrel = (Qm*(((t-(decal))/t1)^4)*exp(4*(1-(t-(decal))/t1))) + QpumpR

   TCa' = (0.001)*(Qb - Qa)
   TCaA' = (0.001)*(Qa - Qr - Qd2)
   TA' = (0.001)*(Qr - Qd - Qd1)
   cai' = (1e-6)*(Qrel - Qpump - Qb + Qr + Qd2)
  
}