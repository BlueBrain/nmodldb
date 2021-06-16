TITLE Kdr channels with DA and truncation/restoration
 
COMMENT
Hodgkin and Huxley potassium channels

Stochastic equations with diffusion aproximation and a Truncation-Restoration procedure
DA equations as in Orio & Soudry (2012) PLoS One, Truncated-Restored algorithm from 
Huang et al. (2013) Phys Rev E 87:012716  DOI: 10.1103/PhysRevE.87.012716

Implemented for Pezo, Soudry and Orio (2014) Front Comp Neurosci 

ENDCOMMENT
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
}
 
NEURON {
        SUFFIX KI
        USEION k READ ek WRITE ik
        RANGE gkbar, NK, scale_a, R, stsum, en, an, bn, N, se
        GLOBAL gu_K
 	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
		se = -1
        gkbar = .004 (S/cm2)	<0,1e9>
        gu_K = 20e-12	(S)
        scale_a = 1.0
}
 
STATE {	
	n0
	n1
	n2
	n3
	n4
}
 
ASSIGNED {
	NK
	area	(micron2)
	v (mV)
	celsius (degC)
	ek (mV)
	dt (ms)
	ik (mA/cm2)
	an	(/ms)
	bn	(/ms)
	stsum
	N
	R[4]	(/ms)
	en[5]
}
 
BREAKPOINT {
    SOLVE states METHOD euler :Restoration procedure inly works with euler's method
	ik = gkbar*n4*(v - ek)   
}
 
 
INITIAL {
	rates(v)
	NK = floor((1e-8)*gkbar*area/gu_K + 0.5) :round to nearest integer
	if (se>=0) {set_seed(se)}	


	N=an/bn
	stsum=(1+N)^4
	n0=1/stsum
	n1=4*N/stsum
	n2=6*N^2/stsum
	n3=4*N^3/stsum
	n4=N^4/stsum
	fijaerror()

}

DERIVATIVE states {  
	rates(v)
	
	n0' = -4*an*n0 + bn*n1 + R[0] + en[0]/dt
	n1' = (-3*an-bn)*n1 + 4*an*n0 + 2*bn*n2 - R[0] + R[1] + en[1]/dt
	n2' = (-2*an-2*bn)*n2 + 3*an*n1 + 3*bn*n3 -R[1] + R[2] + en[2]/dt
	n3' = (-an-3*bn)*n3 + 2*an*n2 + 4*bn*n4 -R[2] + R[3] + en[3]/dt
	n4' = -4*bn*n4 + an*n3 -R[3] + en[4]/dt

	ntrunca()
}
 
PROCEDURE ntrunca() { :Trunca de acuerdo a Huang
	LOCAL N[5], i, j, k, l, nsum, nsumN, ps, aux, aux2, pos, pos2[5], en_aux[5]
	UNITSOFF
	nsum = n0+n1+n2+n3+n4
	N[0]=n0/nsum
	N[1]=n1/nsum
	N[2]=n2/nsum
	N[3]=n3/nsum
	N[4]=n4/nsum

	nsumN = 1
	aux=0
	aux2=0
	l=0
	FROM i=0 TO 4 {
		if (N[i]>1) {
			aux=1
			pos = i
			VERBATIM
			break;
			ENDVERBATIM
			}
		if (N[i]<0) {
			aux=2
			pos2[l] = i
			l=l+1
			}
		}
	if (aux == 0) {
		FROM l=0 TO 4 {en[l]=0}
		}
	if (aux == 1) {
		aux2 = N[pos]
		FROM j=0 TO 4 {
			en[j]=N[j]
			N[j]=0
			}
		en[pos]=aux2-1
		N[pos]=1
	}
	if (aux == 2) {
		FROM n = 0 TO (l-1) {
			ps=pos2[n]
			en_aux[n]=N[ps]
			aux2=aux2+N[ps]
			}
		FROM k = 0 TO 4 {
			en[k]=N[k]*(1-1/(nsumN-aux2))
			N[k]=N[k]/(nsumN-aux2)
			}
		FROM n = 0 TO (l-1) {
			ps=pos2[n]
			en[ps]=en_aux[n]
			N[ps]=0
			}
	}

	:::::::Entrega Truncamiento:::::::
	n0=N[0]
	n1=N[1]
	n2=N[2]
	n3=N[3]
	n4=N[4]
	UNITSON
}

PROCEDURE fijaerror() {
	UNITSOFF
	FROM k=0 TO 4 {en[k]=0}
	UNITSON
}


PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

UNITSOFF
    
    :"n" activation 
    an = scale_a*.01*vtrap(-(v+55),10) 
    bn = scale_a*.125*exp(-(v+65)/80)

	FROM ii=0 TO 3 {R[ii]=normrand(0,1/sqrt(NK*dt))}
	R[0] = R[0]*sqrt(4*an*n0+bn*n1)
	R[1] = R[1]*sqrt(3*an*n1+2*bn*n2)
	R[2] = R[2]*sqrt(2*an*n2+3*bn*n3)
	R[3] = R[3]*sqrt(an*n3+4*bn*n4)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
