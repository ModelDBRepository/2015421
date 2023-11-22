This code package was developed to generate and analyze simulated neuronal networks in which we simulated the effects of BDNF and glutamate on excitatory and inhibitory neurons and connections. 

Software Requirements: Python3.7; MATLAB2018a or later
Note: The Python scripts rely on Brian2, which can be downloaded at [https://brian2.readthedocs.io/en/stable/introduction/install.html](https://brian2.readthedocs.io/en/stable/introduction/install.html).

It includes, in order of use:
- `bdnfHomeostasisSimulation.py`
- `bdnfHomeostasisMain.m`
  - `convertToMEA.m`
    - `neuronToMEA.m`
  - `bdnfHomeostasisMEA_spikeAnalysis.m`
    - `splitIntoSegments.m`
    - `INSIbyElectrode.m`
    - `fanoFactorFn.m`
  - `bdnfHomeostasisMEA_networkAnalysis.m`
    - `functionalConnectivityMatrix.m`
    - `efficiency_wei.m`
  - `bdnfHomeostasisMEA_stdpAnalysis.m`

The only files that need to be executed by the user are `bdnfHomeostasisSimulation.py` and `bdnfHomeostasisMain.m`. 
`bdnfHomeostasisSimulation.py` runs the simulation and `bdnfHomeostasisMain.m` calls all the subsequent analysis functions.

This code can be accessed at [https://www.seas.upenn.edu/~molneuro/](https://www.seas.upenn.edu/~molneuro/) and on ModelDB.

It was developed by Erin D. Anderson. For questions, contact E. Anderson at anderin@seas.upenn.edu

Citation:

O'Neill KM, Anderson ED, Mukherjee S, Gandu S, McEwan SA, Omelchenko A, Rodriguez AR, Losert W, Meaney DF, Babadi B, Firestein BL. Time-dependent homeostatic mechanisms underlie Brain-Derived Neurotrophic Factor action on neural circuitry. *Comms Bio*, 2023.