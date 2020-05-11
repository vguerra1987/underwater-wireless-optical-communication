function outputPower = outputPower(type, params, angle)

if strcmp(type,'Lambertian')
    params = params(1);
    outputPower = (params+1)/2/pi*cos(angle).^(params+1);
elseif strcmp(type,'HG')
    outputPower = 1/4/pi;
else
    outputPower = 1/4/pi;  
end