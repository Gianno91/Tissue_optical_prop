%% TISSUE OPTICAL PROPERTIES CALCULATOR:
% Author: Luca Giannoni
% Version: 3
% Date of the current version: December, 14th, 2016
% Biomedical Optics Research Laboratory (BORL),
% Department of Medical Physics and Biomedical Engineering,
% University College London (UCL), London, UK

function [mu_s,mu_a] = Tissue_properties_v3(lambda,g,a,b,C_HbO2,C_HHb,W,F,C_oxCCO,C_redCCO)

% DESCRIPTION: The following function calculates the optical scattering and the total absorption coefficients of a
% specific biological tissue, given the inputs related to such tissue and to the selected wavelenght of light. The
% equations used in the function are mostly based on the following paper: Steven L. Jacques, "Optical properties of
% biological tissues: a review", Phys. Med. Biol. 58 (2013), pp. 37-61.
% OUTPUTS: 
% 1) mu_s = Scattering coefficient of the tissue [cm^-1];
% 2) mu_a = Absorption coefficient of the tissue [cm^-1];
% INPUTS: 
% 1) lambda = Wavelength of light [nm];
% 2) g = Anisotropy;
% 3) a = Tissue-dependent scaling factor for the scattering coefficient equation [cm^-1];
% 4) b = Tissue-dependent scattering power for the scattering coefficient equation;
% 5) C_HbO2 = Concentration of oxyhemoglobin (HbO2) in the tissue [uM];
% 6) C_HHb = Concentration of deoxyhemoglobin (HHb) in the tissue [uM];
% 7) W = Water content in the tissue;
% 8) F = Fat content in the tissue;
% 9) C_oxCCO = Concentration of oxidized cytochrome c-oxidase (oxCCO) in the tissue [uM];
% 10) C_redCCO = Concentration of reduced cytochrome c-oxidase (redCCO) in the tissue [uM];

%% Calculation of the scattering coefficient of the tissue:
mu_s_reduced=a*((lambda/500)^-b); % Reduced scattering coefficient of the tissue [cm^-1];
mu_s=mu_s_reduced/(1-g); % Scattering coefficient of the tissue [cm^-1];

%% Calculation of the total absorption coefficient of the tissue:
% Read the molar absorption coefficients [uM^-1*cm^-1] of HbO2 and HHb from external file:
filename='Molar absorption coeff HbO2.txt'; % Name of the text file containing the molar absorption coefficients of HbO2;
abs_HbO2=dlmread(filename,'\t',2,0); % Read the file containing the molar absorption coefficients of HbO2 into a matrix;
abs_HbO2(:,2)=abs_HbO2(:,2)*10; % Unit conversion of the molar absorption coefficients of HbO2, from [uM^-1*mm^-1] to [uM^-1*cm^-1];
filename='Molar absorption coeff HHb.txt'; % Name of the text file containing the molar absorption coefficients of HHb;
abs_HHb=dlmread(filename,'\t',2,0); % Read the file containing the molar absorption coefficients of HHb into a matrix;
abs_HHb(:,2)=abs_HHb(:,2)*10; % Unit conversion of the molar absorption coefficients of HHb, from [uM^-1*mm^-1] to [uM^-1*cm^-1];

% Calculate the absorption coefficients [cm^-1] of HbO2 and HHb for the given wavelength:
index=abs_HbO2(:,1)==lambda; % Matrix index of HbO2 and HHb molar absorption coefficients corresponding to the given wavelength;
abs_HbO2=abs_HbO2(index,2); % Molar absorption coefficients of HbO2 for the given wavelength [uM^-1*cm^-1];
abs_HHb=abs_HHb(index,2); % Molar absorption coefficients of HHb for the given wavelength [uM^-1*cm^-1];
mu_a_HbO2=C_HbO2*abs_HbO2; % Absorption coefficients of HbO2 for the given wavelength [cm^-1];
mu_a_HHb=C_HHb*abs_HHb; % Absorption coefficients of HHb for the given wavelength [cm^-1];

% Read the absorption coefficients [cm^-1] of water from external file:
filename='Absorption coeff water.txt'; % Name of the text file containing the absorption coefficients of water;
mu_a_H2O=dlmread(filename,'\t',2,0); % Read the file containing the absorption coefficients of water into a matrix;
index=mu_a_H2O(:,1)==lambda; % Matrix index of the absorption coefficients of water corresponding to the given wavelength;
mu_a_H2O=mu_a_H2O(index,2); % Absorption coefficients of water for the given wavelength [cm^-1];

% Read the absorption coefficients [cm^-1] of fat from external file:
filename='Absorption coeff fat.txt'; % Name of the text file containing the absorption coefficients of fat [m^-1];
mu_a_fat=dlmread(filename,'\t',2,0); % Read the file containing the absorption coefficients of fat [m^-1] into a matrix;
mu_a_fat(:,2)=mu_a_fat(:,2)./100; % Unit conversion of the matrix of absorption coefficients of fat, from [m^-1] to [cm^-1];
index=mu_a_fat(:,1)==lambda; % Matrix index of the absorption coefficients of fat corresponding to the given wavelength;
mu_a_fat=mu_a_fat(index,2); % Absorption coefficients of fat for the given wavelength [cm^-1];

% Read the molar extinction coefficients [uM^-1*cm^-1] of oxidized and reduced CCO from external file:
filename='Extinction coeff oxCCO.txt'; % Name of the text file containing the molar extinction coefficients of oxidized CCO;
eps_oxCCO=dlmread(filename,'\t',2,0); % Read the file containing the molar extinction coefficients of oxidized CCO into a matrix;
eps_oxCCO(:,2)=eps_oxCCO(:,2)./1000; % Unit conversion of the molar absorption coefficients of oxidized CCO in [uM^-1*cm^-1];
filename='Extinction coeff redCCO.txt'; % Name of the text file containing the molar extinction coefficients of reduced CCO;
eps_redCCO=dlmread(filename,'\t',2,0); % Read the file containing the molar extinction coefficients of reduced CCO into a matrix;
eps_redCCO(:,2)=eps_redCCO(:,2)./1000; % Unit conversion of the molar absorption coefficients of reduced CCO in [uM^-1*cm^-1];;

% Calculate the absorption coefficients [cm^-1] of of oxidized and reduced CCO for the given wavelength:
index=find(eps_oxCCO(:,1)==lambda); % Matrix index corresponding to the given wavelength band;
eps_oxCCO=eps_oxCCO(index,2); % Molar extinction coefficients of HbO2 for the given wavelength [uM^-1*cm^-1];
eps_redCCO=eps_redCCO(index,2); % Molar extinction coefficients of HHb for the given wavelength [uM^-1*cm^-1];

% Calculate the absorption coefficient [cm^-1] of the tissue for the given wavelength:
mu_a=mu_a_HbO2+mu_a_HHb+W*mu_a_H2O+F*mu_a_fat+...
    log(10)*C_oxCCO*eps_oxCCO+log(10)*C_redCCO*eps_redCCO; % Absorption coefficient of the tissue [cm^-1];

end % [END OF FUNCTION]