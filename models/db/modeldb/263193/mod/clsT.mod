TITLE clsT.mod   squid sodium, potassium, and leak channels

COMMENT
This is a modification of hh.mod, the original Hodgkin-Huxley model provided by NEURON.
Na-K pumps, and the ability to left-shift a portion of the sodium channels, were added to reproduce 'Coupled left-shift of Nav channels: modeling the Na(+)-loading and dysfunctional excitability of damaged axons', by: Boucher PA, Joos B, Morris CE (2012), J Comput Neurosci.
(Note: The final left-shift voltage is set by the parameter vLeftShift, whereas vLS may initially be zero, or ramped etc.).

Temperature effects reported in (https://doi.org/10.1063/1.5040288)
'Cooling reverses pathological bifurcations to spontaneous firing caused by mild traumatic injury', by: B. M. Barlow, B. Joos, A. K. Trinh, and A. Longtin 
are modeled using four so-called "Q-ten's": qPump, qNa, qK, qGate.



The following comments (in quotes) are those of the original hh.mod author:

"This is the original Hodgkin-Huxley treatment for the set of sodium, 
potassium, and leakage channels found in the squid giant axon membrane.
("A quantitative description of membrane current and its application 
conduction and excitation in nerve" J.Physiol. (Lond.) 117:500-544 (1952).)
Membrane voltage is in absolute mV and has been reversed in polarity
from the original HH convention and shifted to reflect a resting potential
of -65 mV.
Remember to set celsius=6.3 (or whatever) in your HOC file.
See squid.hoc for an example of a simulation using this model.
SW Jaslove  6 March, 1992"
ENDCOMMENT

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (mol) = (1)
    (molar) = (mol/liter)
    (mM) = (millimolar)
    (uM) = (micromolar)
    (um) = (micrometer)
    FARADAY = (faraday) (coulombs)
}

: interface
NEURON {
       SUFFIX clsT
       USEION na READ ena, nai, nao WRITE ina, nai, nao
       USEION k READ ek, ki, ko WRITE ik, ki, ko

       NONSPECIFIC_CURRENT il
       ::: vLS = Voltage left-shift VARIABLE (Cathy Morris and Bela Joos)
       RANGE vLeftShift, AC, ACpotassium, timeLS, vLS0, vLS
       RANGE gnabar, gkbar, gna, gk, gnal, gkl, gl, el, INaKmax, ink, Kmko, Kmnai
       RANGE nai0, nao0, ki0, ko0
       RANGE qPump_0, qNa_0, qK_0, qGate_0 
       GLOBAL minf, hinf, ninf, mtau, htau, ntau, mLSinf, hLSinf, nLSinf, mLStau, hLStau, nLStau
  THREADSAFE : assigned GLOBALs will be per thread
}

PARAMETER {
        
        vLeftShift = 0.0 (mV) 
        timeLS = 0.0 (ms)
        vLS0 = 0.0 (mV)
          
        AC = 0.0  <0,1.0> ::: Proportion of affected (left-shifted) SODIUM channels on node
        ACpotassium = 0.0  <0,1.0> ::: Proportion of affected (left-shifted) POTASSIUM channels on node

        INaKmax = 9.09e-2       (mA/cm2) <0,1e6>




        ink0 = 0.010731         (mA/cm2) <0,1e6> 
        Kmnai =    10          (mM)    <0,1e6>
        Kmko =     3.5         (mM)    <0,1e6>



        :::Boucher
        nai0 = 20.0 (mM)
        nao0 =  154.0 (mM)
        ki0 =  150.0 (mM)
        ko0 =  6.0 (mM)


        gnabar = 0.12 (S/cm2) <0,1e9>
        gkbar = 0.036 (S/cm2) <0,1e9>
        gl = 0.0005 (S/cm2) <0,1e9>
        el = -59.9 (mV)
        gnal = 0.00025 (S/cm2) <0,1e9>
        gkl = 0.0001 (S/cm2) <0,1e9>

        :::Temperature Factors (so-called "Q-ten's")
        qPump_0 = 1.9
        qNa_0 = 1.4
        qK_0 = 1.1
        qGate_0 = 3.0


}

STATE {

      m
      h
      n
      mLS
      hLS
      nLS
      nai (mM)
      nao (mM)
      ki (mM)
      ko (mM)
      
}

ASSIGNED {
        v (mV)
        vLS (mV)
        celsius (degC)
        ena (mV)
        ek (mV)

        ink (mA/cm2)

        gna (S/cm2)
        gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf
        hinf
        ninf
        mtau (ms) 
        htau (ms) 
        ntau (ms)

        mLSinf
        hLSinf
        nLSinf
        mLStau (ms) 
        hLStau (ms) 
        nLStau (ms) 

        diam (um)
        L (um)
        
}

INITIAL {
  v = -59.8137886
  vLS = vLS0
  rates(v)
  m = minf
  h = hinf
  n = ninf
  mLS = mLSinf
  hLS = hLSinf
  nLS = nLSinf


  nai = nai0
  nao = nao0
  ki = ki0
  ko = ko0
  ink = ink0
}

: currents
BREAKPOINT {
  LOCAL qPump, qNa, qK

  UNITSOFF
  qPump = qPump_0^((celsius - 20.0)/10.0)
  qNa = qNa_0^((celsius - 20.0)/10.0)
  qK = qK_0^((celsius - 20.0)/10.0)
  UNITSON

  SOLVE states METHOD cnexp 
  ink= qPump*INaKmax/(((1 + (Kmnai/nai))^3)*((1 + Kmko/ko)^2))

  gna = qNa*gnabar*(m*m*m*h*(1.0 - AC) + mLS*mLS*mLS*hLS*AC)
  ina = (gna + gnal)*(v - ena) + 3*ink
  gk = qK*gkbar*(n*n*n*n*(1.0 - ACpotassium) + nLS*nLS*nLS*nLS*ACpotassium )
  ik = (gk + gkl)*(v - ek) - 2*ink 
  il = gl*(v - el)
}



: states
DERIVATIVE states {
        LOCAL nai_prime, ki_prime, qGate
        UNITSOFF
        qGate = qGate_0^((celsius - 20.0)/10.0)
        UNITSON

        rates(v)
        m' =  qGate*(minf-m)/mtau
        h' = qGate*(hinf-h)/htau
        n' = qGate*(ninf-n)/ntau

        mLS' =  qGate*(mLSinf-mLS)/mLStau
        hLS' = qGate*(hLSinf-hLS)/hLStau
        nLS' =  qGate*(nLSinf-nLS)/nLStau


        nai_prime = -(1e4)*4*ina/(FARADAY*diam)
        ki_prime = -(1e4)*4*ik/(FARADAY*diam)
        nai' = nai_prime
        nao' = -nai_prime
        ki' = ki_prime
        ko' = -ki_prime

}




: rates
PROCEDURE rates(v(mV)) {
      :Computes rate and other constants at current v.
      :Call once from HOC to initialize inf at resting v.
      LOCAL  alpha, beta, sum
      :TABLE minf, mtau, hinf, htau, ninf, ntau, mLSinf, mLStau, hLSinf, hLStau DEPEND celsius FROM -100 TO 100 WITH 200

      UNITSOFF
      :"m" sodium activation system
      alpha = .1 * vtrap(-((v)+40),10)
      beta =  4 * exp(-((v)+65)/18)
      sum = alpha + beta
      mtau = 1/(sum)
      minf = alpha/sum
      :"h" sodium inactivation system
      alpha = .07 * exp(-((v)+65)/20)
      beta = 1 / (exp(-((v)+35)/10) + 1)
      sum = alpha + beta
      htau = 1/(sum)
      hinf = alpha/sum
      :"n" potassium activation system
      alpha = .01*vtrap(-(v+55),10) 
      beta = .125*exp(-(v+65)/80)
      sum = alpha + beta
      ntau = 1/(sum)
      ninf = alpha/sum

      if ( t < timeLS) {
               vLS = vLS0
       }else{
               vLS = vLeftShift
       }

      :"mLS" Affected (Left Shifted) sodium activation system
      alpha = .1 * vtrap(-((v+vLS)+40),10)
      beta =  4 * exp(-((v+vLS)+65)/18)
      sum = alpha + beta
      mLStau = 1/(sum)
      mLSinf = alpha/sum
      :"hLS" Affected (Left Shifted) sodium inactivation system
      alpha = .07 * exp(-((v+vLS)+65)/20)
      beta = 1 / (exp(-((v+vLS)+35)/10) + 1)
      sum = alpha + beta
      hLStau = 1/(sum)
      hLSinf = alpha/sum
      :"nLS" Affected (Left Shifted) potassium activation system
      alpha = .01*vtrap(-((v+vLS)+55),10) 
      beta = .125*exp(-((v+vLS)+65)/80)
      sum = alpha + beta
      nLStau = 1/(sum)
      nLSinf = alpha/sum

}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
       if (fabs(x/y) < 1e-6) {
               vtrap = y*(1 - x/y/2)
       }else{
               vtrap = x/(exp(x/y) - 1)
       }
}

UNITSON