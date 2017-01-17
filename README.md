# Tissue-optical-properties

Calculation of optical properties of biological tissues for Monte Carlo simulations.

It contains a MATLAB function for the calculation of the absorption and scattering coefficients, plus a text file with the molar extinction coefficients of HbO2 and HHb according to different wavelengths. 

All the equations for the estimate of the optical properties and different values of the inputs to the MATLAB function are taken from the following paper (also included): 
  Steven L. Jacques, "Optical properties of biological tissues: a review", Phys. Med. Biol. 58 (2013), pp. 37-61

UPDATE: A new version (v.2) is available. It takes into account also CCO (oxidized and reduced), plus additional absorption/extinction coefficients libraries have been added for HbO2, HHb, water, fat and CCO.

UPDATE (Jan. 17th, 2017): A new version (v.3.) is available. The input for HbO2 and HHb have been modified to allow the selection of specific molar conecentrations of each of the chromophores, instead of using total hemoglobin concentration and oxygenation saturation.
