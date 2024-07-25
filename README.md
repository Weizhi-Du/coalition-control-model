The optimization of resource distribution is a challenging topic in today's society. To explore this topic, we develop a Coalition Control Model(CCM) based on the Model Predictive Control(MPC) and test it using a fishing model with linear parameters. The fishing model focuses on the problem of distributing fishing fleets in certain regions to maximize fish caught using either exhaustive or heuristic search. Our method introduces a communication mechanism to allow an individual fishing fleet to merge or split, after which new coalitions can be automatically formed. Having stabilized the coalition structure, the system reached an equilibrium state through the Nash-Bargaining process. Our experiments on the hypothetical fishing model demonstrated that the CCM can dynamically distribute limited resources in complex scenarios.


Recommendations for the laboratory environment:
memory of PC: 8GB and above
MATLAB: R2015b and above
(Our experiments were performed on PC with Intel i7-7700CPU, having 3.6-GHz clock speed and 56 GB of memory, and all the methods were implemented and tested using MATLAB R2016a.)

Our code includes all the work for grand coalition, isolated coalition, controlled coalition and accelerated control coalition, including:
1. The folder grand contains the entire contents of the grand coalition, the running script is: grand.m
2. The folder isolated contains the entire contents of the grand coalition, the running script is:  isolated.m
3. The folder coal_control contains the entire contents of the grand coalition, the running script is: coal_control.m
4. The folder ACCM (it means Accelerated Coalition Control Model) contains the entire contents of the accelerated control coalition, the running script is: ACCM.m

In addition, some encapsulated functions such as:
1. Equilibrium function: equilibrium.m
2. Converge mpc: convergempc.m
3. MPC solver: mpc_solver.m
4. Objective function: myfun.m
5. The strategy of coalition: itercol.m (belong to ACCM)

The approximate running times for each time step (each day) of the individual algorithms are as follows (unit: s):

Method	      4x6	4x12	4x18	4x24
isolated      3.0       6.6    9.6     13.1
nocross	      6.1       25.8    156.1   197.3
coal_control  382.7     6391.9  NaN     NaN
accelerate    14.5      46.9    71.2    115.2

PS: The running time is influenced by the convergence of MATLAB's built-in function fmincon when solving and is highly dependent on the input data.
