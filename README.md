Aircraft Panel Generation and Aerodynamic Analysis (APGAA)

APGAA is a MATLAB-based software designed for aerodynamic analysis using the Vortex Lattice Method (VLM) and Doublet Lattice Method (DLM). The solvers in APGAA are derived from the PanelAero software [1], which has been validated against Nastran, ensuring accuracy and reliability. For ease of implementation, APGAA has been validated directly against PanelAero, and several test cases demonstrate that APGAA produces results in perfect agreement with PanelAero.

Key Features:
Aerodynamic Model Generation: APGAA can generate aerodynamic models for the entire aircraft, including both lifting and stabilizing surfaces, as well as control surfaces.
Control Surface Analysis: The software allows for the adjustment of control surface angles to compute maneuver loadings.
Aerodynamic Influence Coefficient (AIC): APGAA calculates the AIC, which can be used for dynamic aeroelastic analysis.

Validation:
APGAA has been rigorously tested and validated against PanelAero, with all test cases showing perfect matches. This ensures that users can trust the accuracy of the aerodynamic analysis performed by the software.

Reference:
[1] Voß, A., “An Implementation of the Vortex Lattice and the Doublet Lattice Method,” Institut für Aeroelastik, Deutsches Zentrum für Luft- und Raumfahrt, Göttingen, Germany, Technical Report DLR-IB-AE-GO-2020-137, October 2020, https://elib.dlr.de/136536/.

How to run:
1) run the 'init.m' file in the main folder to add required directory to MATLAB path.
2) All examples and tutorials are included in the /examples directory. You can browse /examples directory and run the script/live-script files.
3) You can write your old script and input files to run new created model.
