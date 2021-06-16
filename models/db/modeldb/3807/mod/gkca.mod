TITLE gkca.mod   
   
UNITS {          
   (uM) = (micro/liter)          
   (mM) = (milli/liter)          
   (mA) = (milliamp)          
   (mV) = (millivolt)   
}          
    
NEURON {          
   SUFFIX gkca         
   USEION ca READ cai      
   USEION k READ ek WRITE ik          
   RANGE gkcabar, ikca
   GLOBAL cinf, cexp
}          

STATE {
c
}          
          
INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }          
          
PARAMETER {          
   celsius = 20 (degC)
   dt (ms)
   v (mV) 
   gkcabar = 0.0008 (mho/cm2)    
   ikca (ma/cm2)
   ek = -68 (mV)           
  }          
          
ASSIGNED {          
   ik  (mA/cm2)          
   cai (mM)          
   cinf cexp
}          

LOCAL tinc, q10, alpha, beta, sum

BREAKPOINT {          
        q10 = 3^((celsius - 20)/10)
        tinc = -dt * q10
        alpha = 0.1 * (cai / 0.01)
        beta =  0.1 
        sum = alpha + beta
        cinf = alpha/sum
        cexp = 1 - exp(tinc*sum)
        c = c + cexp*(cinf-c)
        ikca = gkcabar * c * (v-ek)
        ik = ikca
}

UNITSOFF
 
INITIAL {
        alpha = 0.1 * (cai / 0.01)
        beta =  0.1 
        sum = alpha + beta
        cinf = alpha/sum
        c = cinf     
}
 