# RFID_Simulation Model
Simulation Model for RFID Localization System
## Run Simulation
```
SimulationVeify.m
```
* Run SimulationVeify.m to load simulated states, groundtruth states, and plot figures;
* In SimulationVeify.m, change sim = 1 for 1D simulation ; sim = 2 for 2D simulation;

## Generate GroundTruth States
```
groundtruth1D.m
```
```
groundtruth2D.m
```
Parameters in the two scripts can be changed corresponding to motion you want to simulate
Also, ground truth states can also be loaded with .m files

## Phase Concatenation 
```
noisysim.m
```
The simulation model is based on phase unwrapping - concatenating phase in 2\pi periods. Ways of concatenation for different routes can be changed in noisysim.m
