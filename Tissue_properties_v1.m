%% TISSUE OPTICAL PROPERTIES CALCULATOR:
% Author: Luca Giannoni
% Version: 1
% Date of the current version: November, 22nd, 2016
% Biomedical Optics Research Laboratory (BORL),
% Department of Medical Physics and Biomedical Engineering,
% University College London (UCL), London, UK

function [mu_s,mu_a] = Tissue_properties_v1(lambda,g,a,b,B,S,C_HGb,...
    W,mu_a_H2O,F,mu_a_fat,M,mu_a_mel,C_bili,eps_bili,C_betaC, eps_betaC)

% DESCRIPTION: The following function calculates the optical scattering and the absorption coefficients of a
% specific biological tissue, given the inputs related to such tissue and to the selected wavelenght of light.
% The equations used in the function are taken from the following paper: Steven L. Jacques, "Optical properties
% of biological tissues: a review", Phys. Med. Biol. 58 (2013), pp. 37-61.
% OUTPUTS: 
% 1) mu_s = Scattering coefficient of the tissue [cm^-1];
% 2) mu_a = Absorption coefficient of the tissue [cm^-1];
% INPUTS: 
% 1) lambda = Wavelength of light [nm];
% 2) g = Anisotropy;
% 3) a = Tissue-dependent scaling factor for the scattering coefficient equation [cm^-1];
% 4) b = Tissue-dependent scattering power for the scattering coefficient equation;
% 5) B = Average blood volume fraction in the tissue;
% 6) S = HGb oxygen saturation of mixed arterio-venous vasculature in the tissue;
% 7) C_HGb = Concentration of hemoglobin in the tissue [M];
% 8) W = Water content in the tissue;
% 9) mu_a_H2O = Absorption coefficient of water at the given wavelength [cm^-1];
% 10) F = Fat content in the tissue;
% 11) mu_a_fat = Absorption coefficient of fat at the given wavelength [cm^-1];
% 12) M = Melanosome volume fraction in the tissue;
% 13) mu_a_mel = Absorption coefficient of melanosome at the given wavelength [cm^-1];
% 14) C_bili = Concentration of bilirubin in the tissue [M];
% 15) eps_bili = Molar extinction coeffcient of bilirubin at the given wavelength [cm^-1*M^-1];
% 16) C_betaC = Concentration of beta-carotene in the tissue [M];
% 17) eps_betaC = Molar extinction coeffcient of beta-carotene at the given wavelength [cm^-1*M^-1];

%% Calculation of the scattering coefficient of the tissue:
mu_s_reduced=a*((lambda/500)^-b); % Reduced scattering coefficient of the tissue [cm^-1];
mu_s=mu_s_reduced/(1-g); % Scattering coefficient of the tissue [cm^-1];

%% Calculation of the absorption coefficient of the tissue:
% Read the molar extinction coefficients [M^-1*cm^-1] of HbO2 and HHb from external file:
filename='Extinction coeff_HbO2_HHb.txt'; % Name of the text file containing the molar extinction coefficients;
Ext_coeff=dlmread(filename,'\t',2,0); % Read the file containing the molar extinction coefficients into a matrix;

% Calculate the absorption coefficients [cm^-1] of HbO2 and HHb for the given wavelength:
Ext_coeff_index=find(Ext_coeff(:,1)==lambda); % Matrix index corresponding to the given wavelength band;
eps_HbO2=Ext_coeff(Ext_coeff_index,2); % Molar extinction coefficients of HbO2 for the given wavelength [M^-1*cm^-1];
eps_HHb=Ext_coeff(Ext_coeff_index,3); % Molar extinction coefficients of HHb for the given wavelength [M^-1*cm^-1];
mu_a_HbO2=log(10)*C_HGb*eps_HbO2; % Absorption coefficients of HbO2 for the given wavelength [cm^-1];
mu_a_HHb=log(10)*C_HGb*eps_HHb; % Absorption coefficients of HHb for the given wavelength [cm^-1];

% Calculate the absorption coefficient [cm^-1] of the tissue for the given wavelength:
mu_a=B*S*mu_a_HbO2+B*(1-S)*mu_a_HHb+W*mu_a_H2O+F*mu_a_fat+...
    M*mu_a_mel+log(10)*C_bili*eps_bili+log(10)*C_betaC*eps_betaC; % Absorption coefficient of the tissue [cm^-1];

end % [END OF FUNCTION]