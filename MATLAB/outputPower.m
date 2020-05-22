function outputPower = outputPower(type, params, angle)

if strcmp(type,'Lambertian')
    params = params(1);
    outputPower = (params+1)/2/pi*cos(angle).^(params+1);
elseif strcmp(type,'HG')
    outputPower = 1/4/pi*(1 - params^2)/...
        (1 + params^2 - 2*params*cos(angle))^(3/2);
else
    outputPower = 1/4/pi;  
end