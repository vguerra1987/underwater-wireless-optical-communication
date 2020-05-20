function [time, h_t] = project_response(impulse_response, scenario, optics)
% This function carries out the projection of the pentadimensional impulse
% response to a time-amplitude vector. Currently we are assuming that the
% optics is always a CPC, but we can extend this funcionality to support
% different types of receiver using inheritance in an OOP scheme.

% Memory allocation
time = zeros(1, length(impulse_response));
h_t = zeros(size(time));

% CPC's cosine of FWHM FOV/2. This is used to check if a ray contributes 
% effectively or not
cos_half_fov = cos(optics.params/2);

% Lens gain (CPC gain)
lens_gain = (optics.n/sin(optics.params/2))^2;

% Counter to iterate on time and h_t
counter = 1;
for ray = impulse_response.storage
    t = ray{1}.distance*scenario.n/3e8;
    impact_cosine = dot(-ray{1}.direction,optics.orientation);
    
    % Entry condition check (ray is within rx's FOV)
    if (cos_half_fov < impact_cosine)
        time(counter) = t;
        h_t(counter) = ray{1}.power*impact_cosine*optics.area*lens_gain;
    end
    
    counter = counter + 1;
end

% Finally, since much rays did not enter, we trim the vectors and sort them
% by time.

% Trimming
indices = h_t > 0;
time = time(indices);
h_t = h_t(indices);

% Sorting
[time, indices] = sort(time);
h_t = h_t(indices);
