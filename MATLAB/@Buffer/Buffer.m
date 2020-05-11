classdef Buffer < handle
% This is a FIFO queue

% The main property is the vector-like storage
properties
    storage = {};
end

% Methods are pop and push
methods
    % Constructor
    function obj = Buffer()
    end
    
    % Data insertion
    function push(obj, data)
        obj.storage{end+1} = data;
    end
    
    % Data extraction
    function aux = pop(obj)
        aux = obj.storage{end};
        obj.storage = obj.storage(1:end-1);
    end
end

end