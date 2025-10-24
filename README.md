# fem
Our code is a two-dimensional finite element code written in Fortran.
The following functionalities can be implemented:
Steady-State Heat Transfer Analysis,
Steady-State Thermal Stress Analysis (only displacement boundary conditions verified),
Verification of thermal stress due to differences in linear expansion coefficient, verification of reference temperature, verification of load boundary conditions,
Transient Solution (Incomplete),
Transient Solution (Completed),
Transient Thermal Stress Analysis (Likely implemented, elastoplasticity not verified),
Cases where temperature conditions are functions of time in transient heat transfer analysis,
Consideration of Freezing Expansion (Elastoplastic behavior became problematic),
Coexistence of Freezing Expansion and Elastoplasticity,
Elastoplasticity originating from thermal stress due to freezing expansion, etc.,
Consideration of Heat Transfer,
Thermal Radiation,
Latent Heat of Fusion / Solar Radiation,
Added a switch to enable/disable elastoplasticity. For the elastic case, modified output to include the required uniaxial compressive strength and required tensile strength to prevent failure,
Creep.
