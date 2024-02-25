# mars_ghg_2024

A collection of source code and scripts to reproduce the results of Jorge & Wordsworth, 2024, JGR Planets.

PCM_LBL/ contains PCM_LBL, the line-by-line radiative-convective model used to calculation surface temperature as a function of total pressure, H2 content and solar luminosity. More details on the operation of this code are given in the README file in that directory.

run_PCM_LBL.m is the script used to iterate over greenhouse gases and molar concentrations to calculate surface temperature at a given maximum surface pressure.

plot_TIPS.ipynb is the script used to plot the TIPS ratio and the TIPS temperature approximation.
