%% Impulse response simulation
% Author: Victor Guerra, PhD
% Date: May 2020

% This script carries out the simulation of an UWOC link using a Modified
% Monte Carlo scheme. It calculates the impulse response at a given point
% as a tuple (x,y,z,theta,phi,t, p). This way, calculating the effective
% impulse response for any receiver is straightforward (Glens, FOV, etc...)

% The way the algorithm operates is easy, and is based on three main steps:
%   - Emit a forced contribution (directed to the receiver)
%   - Emit Nr random contributions following a given radiation pattern
%   - Calculate random impacts with suspended matter

% The algorithm classifies each ray as random or direct, and only direct
% rays will be considered for the final calculation.

%% Simulation parameters
% In general, light energy will scatter indefinitely in an underwater
% scenario. Properly limiting the scenario may help in reducing the
% convergence time of the simulation.

% The limits are expressed in the form (xmin, ymin, zmin, xmax, ymax, zmax)
scenario_limits = [-5, -5, -5, 5, 5, 5]; %1000 m^3 cube

% Wavelength is fundamental in UWOC, since both water absorption and
% scattering are closely related to it.
wavelength = 660; %nm

% Link depth. This is important for parameter estimation.
depth = 15;

% Refractive index of water depends on three parameters: wavelength,
% temperature and salinity. It has been shown that it is closely related to
% depth, since temperature and salinity are.
refractive_index = getRefractiveIndex(wavelength, 'McNeil', depth);

% Regarding scattering, we are using a simplified model based on
% Henyey-Greenstein scattering phase function. This phase function has one
% single parameter, the mean cosine "g".
g = 0.8; % 0 implies isotropic radiation and 1 just forward scattering.

% Chlorophyll concentration
chlorophyll = 0.01; %g/m^3

% Absorption depends on the concentration of several substances, organic
% and inorganic. However, we are considering only chlorophyll.
alfa = getAbsorptionFromWavelength(wavelength, chlorophyll, 0, 0);

% Given alfa and g, it is not straightforward to find a relationship
% between the scattering coefficient beta and the amount of particles.
% Because of that we are approximating beta from a desired expected value
% of extinction c(lambda).
extinction = 0.45; % m^-1

% This is to check if alfa is greater than the expected value of the
% extinction coefficient.
assert(extinction > alfa);

% Thus, beta is the difference between extinction and alfa, and will serve
% to determine the amount of enery a particle absorbs.
beta = extinction - alfa;

% At this point, since we have defined the mean cosine (particle size) and
% chlorophyll concentration, we have implicitly defined beta. However,
% currently we have no clue about the relationship between "g" and particle
% size, and it must be investigated. This is just a reminder about the
% possible misadjustment between the simulated beta and the actual beta.

% Since we are using a Modified Monte Carlo Ray Tracing we will not limit
% the simulation to small-energy rays, we will just discard those rays
% goint outside the cube of interest.

% Transmitter position, must be inside the cube-of-interest
tx_pos = [0,0,-2];

% We define here the Lambertian degeneration factor (m)
lambertian_factor = 1; % Pure Lambertian

% Receiver position, the same constraint as above.
rx_pos = [0,0,4];

% We group parameters to ease the understading of the simulation call
scenario.n = refractive_index;
scenario.alfa = alfa;
scenario.beta = beta;
scenario.extinction = extinction;
scenario.limits = scenario_limits;
scenario.max_hops = 4;
scenario.power_threshold = 1e-14;
scenario.plot = 0;
scenario.info_period = 100; % Display information after 100 iterations

% We have to add the specific simulation parameters
scenario.N = 1000; % Rays on light source
scenario.M = 10;   % Rays on each scatterer

particle.params = g;
particle.type = 'HG';

tx.position = tx_pos;
tx.orientation = [0,0,1];
tx.type = 'Lambertian';
tx.params = lambertian_factor;

rx.position = rx_pos;

%% SIMULATION
impulse_response = monte_carlo(tx, rx, particle, scenario);

%% VISUALIZATION
% At this point we have the aforementioned tuple inside impulse_response.
% Now we have to apply appropriate parameters to the receiver in order to
% orthographically project the result to obtain h(t). It must be taken into
% account that the link range must comply with the far-field approximation
% in order to ensure that dOmega = A/d^2.

% This is a TODO!

optics.type = 'CPC';
optics.params = 30*pi/180; % In this case corresponds to 30 degrees FOV FWHM
optics.area = pi/4*(5e-3)^2; % 5mm diameter receiver
optics.orientation = [0,0,-1]; % pointing vector

[time, h_t] = project_response(impulse_response, scenario, optics);

% Finally we plot the response
scatter(time*1e9, 10*log10(h_t));

[H, BW] = channelParameters(time, h_t);