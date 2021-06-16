TITLE caacum.mod: Calcium ion accumulation   
       
NEURON {       
   SUFFIX ca       
   USEION ca READ ica, cai WRITE cai 
}       
       
INDEPENDENT {t FROM 0 TO 1 WITH 10 (ms)}       
       
UNITS {       
   (molar) = (1/liter)       
   (mV) = (millivolt)       
   (um) = (micron)       
   (mM) = (millimolar)       
   (mA) = (milliamp)       
   FARADAY = 96520 (coul)       
   R = 8.3134     (joule/degC)       
}       
       
PARAMETER {       
   celsius=20     (degC)       
   cabath = 1.8   (mM)       
   diam = 1 (um)       
   ica            (mA/cm2)       
}       
       
STATE {       
   cai 
}       
       
       
INITIAL {       
   VERBATIM       
   cai = 0.0001;
   ENDVERBATIM       
}       
       
BREAKPOINT {       
   SOLVE state METHOD derivimplicit
}       
       
DERIVATIVE state {       
   cai' = -ica * 4/(diam*FARADAY) * (1e4)  
}       
          
COMMENT       
This model uses ica but does not WRITE it; thus this model does       
not add anything to the total ionic current.       
       
The initial block works around a difficulty that arises from a STATE in      

this model having the same name as an ion.  (Note: in the cabpump model      

there is no name conflict between the ca[] states and the cai ion       
concentration.) The sequence of events when finitialize is called is       
that the ca_ion's cai is initialized to the global variables       
cai0_ca_ion. Then this model's INITIAL block      
  
is called. By default, cai/cao would be set to the initial state values      

cai0/cao0 implicitly declared in this model and on exit from the intitial    

 
  
block, the ca_ion values would be assigned these local values. We therefore  

 
  
   
    
directly set the local state values to the ca_ion values. See the       
"nocmodl cacum" generated caacum.c file to see the precise sequence on       
the nrn_init() call.       
       
diam is a special range variable in NEURON and refers to the diameter in     

 
microns.  Under scop and hocmodl its default value is specified in the       
PARAMETER block. In NEURON, however, its value is taken from the       
"morphology" mechanism.       
ENDCOMMENT       
  