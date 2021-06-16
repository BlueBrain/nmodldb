NEURON{
SUFFIX nafast2
USEION na READ ena WRITE ina
 RANGE g, gbar
}

UNITS{
(mV) = (millivolt)
(S) = (siemens)
(mA)=(milliamp)
}

PARAMETER{
gbar=0.035 (S/cm2)
vshift=0 (mV)
}

ASSIGNED{
g (S/cm2)
v (mV)
ena(mV)
ina(mA/cm2)
o1c1 (/ms)
i1c1 (/ms)
o1c2 (/ms)
o1c3 (/ms)
i1c3 (/ms)
c1o1 (/ms)
c2o1 (/ms)
c3o1 (/ms)
c1i1 (/ms)
c3i1 (/ms)
}

STATE {c1 c2 c3 o1 i1 }

BREAKPOINT {SOLVE kin METHOD sparse
g=gbar*o1
ina=g*(v-ena)
}

INITIAL {LOCAL sum
c1=1
c2 = exp(4.4751000 + (0.0128154*v))
c3 = exp(10.1000000 + (0.0988810*v))
o1 = exp(14.4191083 + (0.2904200*v))
i1 = exp(19.1720000 + (0.3424800*v))
sum = c1+c2+c3+o1+i1
c1 = c1/sum
c2 = c2/sum
c3 = c3/sum
o1 = o1/sum
i1 = i1/sum
}

KINETIC kin {
rates(v)
~ c1<->o1 (c1o1, o1c1)
~ c2<->o1 (c2o1, o1c2)
~ c3<->o1 (c3o1, o1c3)
~ c1<->i1 (c1i1, i1c1)
~ c3<->i1 (c3i1, i1c3)
CONSERVE c1 + c2 + c3 + o1 + i1 =1
}

PROCEDURE rates(vm(mV)){
UNITSOFF
o1c1=exp(0.4689458 + (vm*0.0054700))
i1c1=exp(-3.1509423 + (vm*0.0112800))
o1c2=exp(-4.9720042 + (vm*-0.2123123))
o1c3=exp(-2.1595542 + (vm*-0.1994914))
i1c3=exp(-4.5213608 + (vm*-0.0296084))
c1o1=exp(14.8880542 + (vm*0.2958900))
c2o1=exp(4.9720042 + (vm*0.0652923))
c3o1=exp(2.1595542 + (vm*-0.0079524))
c1i1=exp(16.0210577 + (vm*0.3537600))
c3i1=exp(4.5506392 + (vm*0.2139906))
UNITSON
}