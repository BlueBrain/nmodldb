NEURON {       
   SUFFIX capump       
   USEION ca READ cai WRITE ica       
   RANGE ica  
}       
       
UNITS {       
   (uM) = (micro/liter)       
   (mM) = (milli/liter)       
   (mA) = (milliamp)       
}       
       
INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }       
       
PARAMETER {       
   final_conc = 0.0001 (mM) 
   tau = 150 (ms) 
   celsius (degC)       
}       
       
INITIAL { 
cai = 0.0001 
} 
 
ASSIGNED {       
   ica (mA/cm2)       
   cai (mM)       
}       
       
LOCAL Q, s_celsius       
       
BREAKPOINT {       
   if (s_celsius*1(degC) != celsius) {       
      s_celsius = celsius       
      Q = 3^((celsius - 6.3)/10 (degC))       
   }               
   ica =  (cai-final_conc)/tau   
}       
 