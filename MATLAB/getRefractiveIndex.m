function n = getRefractiveIndex(wavelength, type, depth)
% This function returns the refractive index according to McNiel or 
% Matthäus approximation. If no depth is selected, it is assumed 15 meters 
% depth. McNeil is assumed by default. Wavelength is in nanometers.

if nargin==1
    type = 'McNeil';
    depth = 15;
elseif nargin == 2
    depth = 15;
end

[T, S] = getParamsFromDepth(depth);

if strcmp(type,'Matthaus')
    L = wavelength/1000;
    n = 1.447824 + 3.011e-4*S - 1.8029e-5*T - 1.6916e-6*T^2 - 0.489*L + ...
        0.728*L^2 - 0.384*L^3 - S*(7.9362e-7*T - 8.06e-9*T^2 + ...
        4.249e-4*L - 5.847e-4*L^2 + 2.812e-4*L^3);
else
    L = wavelength;
    n = 1.3247 - 2.5e-6*T^2 + S*(2e-4 - 8e-7*T) + 3300/L^2 - 3.2e7/L^4;
end