function direction = randomDirection(type, params)
% We are using the Gaussian approximation of the lambertian pattern

if strcmp(type, 'Lambertian')
    if params(1) == 1
        theta = asin(rand);
        phi = 2*pi*rand;
    else
        sigma = sqrt(0.9872/(params(1) + 0.4861));
        theta = sigma*randn;
        while(abs(theta) > pi/2)
            theta = sigma*randn;
        end
        phi = pi*randn;
    end
    
elseif strcmp(type, 'HG')
    % TODO
    theta = 2*pi*rand;
    phi = pi*rand;
end

direction(1) = sin(theta)*cos(phi);
direction(2) = sin(theta)*sin(phi);
direction(3) = cos(theta);