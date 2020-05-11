function alfa = getAbsorptionFromWavelength(wavelength, ...
                                            chlorophyll, ...
                                            gelbstoff, ...
                                            minerals)
% This function returns the absorption coefficient of seawater for a given
% wavelength. This is carried out using the pure water data from Pope and
% Fry (1995). Their measurements were in cm^-1, and we must scale them to m^-1.
% For chlorophyll we are using Bricaud's dataset (2004). Gelbstoff and
% mineral absorption are modeled as exponentials. Brigaud measured not only
% chlorophyll but also other pigments. These pigments should be considered
% theoretically depending on the type of suspended microalgae. In this work
% only chlorophyll is considered for simplicity. The units are in m^2/g. We
% have to normalize it respect to their specific volumetric concentration,
% usually expressed in g/m^3 of mg/m^3

% pure water absorption
data = readtable('Pope_Fry_measurements.csv');
water = 100*interp1(data.Wavelength, data.Absorption, wavelength);

% chlorophyll absorption
data = readtable('Bricaud_et_al_2004.xlsx');
ChlA = chlorophyll*interp1(data.lambda, data.ChlA, wavelength);
ChlB = chlorophyll*interp1(data.lambda, data.ChlB, wavelength);
% Here we could introduce accurately the absorption spectra of each
% specific pigment and finally carry out a weighed sum to obtain the
% resulting spectrum. Here 50%-50% chlorophylls A and B.
phytoplankton = mean([ChlA; ChlB],1);

% gelbstoff
decaying = exp(-0.0139*(wavelength - 400));

% minerals
inorganic = exp(-0.0069*(wavelength - 400));

% Finally we apply Haltrin's model to include all the contributions
alfa = water + phytoplankton*chlorophyll + ...
    + gelbstoff*decaying + minerals*inorganic;