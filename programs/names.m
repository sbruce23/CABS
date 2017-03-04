function [n] = names(X)
%Function returns the names of the object
%   Input: X is a struct, e.g. object from CABS(), object from setMCMCOptions() etc.
%   Output: n is a cell array of names for X
    if isstruct(X) == 1
        n = fieldnames(X)';
    else
        fprintf(1,'Warning: only struct is supported.\n');
        n = [];
    end
end
