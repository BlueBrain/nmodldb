# CONNECTION PROPERTIES
#
minPsdSize:               10.0  # Unpot. conn drifts towards this
maxPsdSize:               100.0

minNumCiAmpars:           2.0   # Unpot. conn drifts towards this
minNumCpAmpars:           0.0   # Unpot. conn drifts towards this

potProbK:                 1     # K value for potentiation probability function
potProbHalf:              15    # numStimCycles at which pot prob is 0.5

# UNIT PROPERTIES
#
actFuncK:                 2.0
actThreshold:             0.0

# LAYER PROPERTIES
#
minInhibition:            2.5
maxInhibition:            10
initInhibition:           0
inhibIncr:                0.05
orthogonalPatterns:       true

W:                        5
H:                        5
K:                        .2

hpc.width:                W
hpc.height:               H
hpc.k:                    K     # Prop. units active at equilibrium

acc.width:                W
acc.height:               H
acc.k:                    K

sc0.width:                W
sc0.height:               H
sc0.k:                    K

sc1.width:                W
sc1.height:               H
sc1.k:                    K

# COMMON TRACT PROPERTY VALUES
#
cpAmparRemovalRate01h:       0.1   # constitutive
ciAmparInsertionRate01h:     2.0   # when potentiated (This is linear!)
ciAmparRemovalRate01h:       0.015 # when unpotentiated
psdDecayRate01h:             0.01  # Shrink rate of unpopulated portion of PSD
e3DecayRate01h:              0.03 
maxE3DepotProb01h:           0.05  # 01h prob of E3 depot when e3Level is 1.0
maxPotProb01h:               1.0   # 01h Max probability of potentiation

# Hippocampus tract properties
#
hpc.AcqLearnRate:            0.08
hpc.ConsLearnRate01h:        0.0
hpc.psdDecayRate01h:         psdDecayRate01h
hpc.cpAmparRemovalRate01h:   cpAmparRemovalRate01h
hpc.ciAmparInsertionRate01h: ciAmparInsertionRate01h
hpc.ciAmparRemovalRate01h:   ciAmparRemovalRate01h
hpc.BaseDepotProb01h:        0.002 # spontaneous depotentiation         
hpc.reactE3Level:            1.0
hpc.e3DecayRate01h:          e3DecayRate01h
hpc.maxE3DepotProb01h:       maxE3DepotProb01h
hpc.maxPotProb01h:           maxPotProb01h

# Neocortex tract properties
#
nc.AcqLearnRate:             0.004
nc.ConsLearnRate01h:         0.004
nc.psdDecayRate01h:          psdDecayRate01h
nc.cpAmparRemovalRate01h:    cpAmparRemovalRate01h
nc.ciAmparInsertionRate01h:  ciAmparInsertionRate01h
nc.ciAmparRemovalRate01h:    ciAmparRemovalRate01h
nc.BaseDepotProb01h:         0.0  # spontaneous depotentiation         
nc.reactE3Level:             0.0
nc.e3DecayRate01h:           e3DecayRate01h
nc.maxE3DepotProb01h:        maxE3DepotProb01h
nc.maxPotProb01h:            maxPotProb01h

# SIMULATION CONTROL
#
stopTime:                    60    # days::hours

trainNumStimCycles:          50
consNumStimCycles:            1
reactNumStimCycles:          50
numSettleCycles:             20
timeStepChanges:             { 0::00 1 1::00 24 } # when (dd::mm), step(hh), ...
numBackgroundPatterns:       0

# OUTPUT CONTROL
#
printPatterns:               true
