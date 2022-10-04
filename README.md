# TumorEvasionAntigenLoss
DOI: 10.5281/zenodo.7145286
MATLab code for dynamical modeling of tumor escape via antigen loss.

PassiveCase.m      - Script containing analysis studying passive immune evasion strategies.
OptimalPolicy.m    - Script containing analytic description of the optimal evasio policy
Immunogenicity.m   - Script containing stochastic trajectories and distributional plots of post-escape distribution, immunogenicity, and cumulative mutational burden

EmbededLegend.m               - Function used to embed text into figure plots
pbreakeven.m                  - Function calculating break-even evasion rate given antigen number and recognition rate so that evasion and escape are equally likely
SnDynamics.m                  - Function simulating the dynamics of passive evasion conditioned on a tie
SnDynamicsSaveLastSn.m        - Same as above with recorded final state
StochasticTrajectoriesTiev2.m - Generates stochastic trajectories assuming optimal evasion.
