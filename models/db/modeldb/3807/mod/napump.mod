NEURON {    
   SUFFIX napump    
   USEION na READ nai WRITE ina    
   USEION k WRITE ik    
   RANGE vmax, khalf,ksteep,ina, ik    
}    
    
UNITS {    
   (uM) = (micro/liter)    
   (mM) = (milli/liter)    
   (mA) = (milliamp)    
}    
    
INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }    
    
PARAMETER {    
   vmax = 0.003
   khalf = 12 (mM)
   ksteep = 1 (mM)     
   celsius (degC)    
}    
    
ASSIGNED {    
   ina (mA/cm2)    
   ik  (mA/cm2)    
   nai (mM)    
}    
    
LOCAL Q, s_celsius    
    
BREAKPOINT {    
   if (s_celsius*1(degC) != celsius) {    
      s_celsius = celsius    
      Q = 3^((celsius - 6.3)/10 (degC))    
   }            
   ina = 3*(vmax/(1+exp((khalf-nai)/ksteep)))*1(ma/cm2)    
   ik = -2/3 * ina   
}    
 