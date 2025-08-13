# Public API

```@meta
CurrentModule = ChargeFlipPhaser
```

## Input data
### Symmetry and lattice parameters
```@docs
DiffractionData
metric_data_inconsistency
physicalnorm
```
### Reflections
```@docs
add_peak!
```
### Form factors
```@docs
WeightedF0
formfactors_synthetic
```

## Running the program
```@docs
Phaser
ball_autocorr
do_phasing!
```

### GUI monitor
```@docs
PhasingMonitor 
Cut2D 
add_panel!
display
MonitorHooks
```

## Custom phasing algorithms
```@docs
AbstractPhasingAlgorithm
WorkingAmplitudes
```

## Modifying program behavior
```@docs
AbstractHooks
```

## Saving results
```@docs
AbstractSaver
CSVSaver
```

