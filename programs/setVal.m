function [X] = setVal(X,varname, val)
% Function sets component of struct from setMCMCOptions() to the value you want 
%   This function would be handy for you to use setMCMCOptions() to set 
%   most of the default values and use setVal() to change couple of them.
%
%   Input: 
%      X:  a struct from setMCMCOptions() 
%      varname:  a character string that defines the variable you want from X
%      val:   value to be set to X  with name of "varname"        
%
%   Output:
%      X:   updated cell array X
%
%     ex:   px =setMCMCOptions();  
%           now, change nbasis to be 10
%           px = setVal(px,'nbasis',10);
%           
% 	showOptionNames
    if isempty(find(strcmp(names(X),varname),1))
        error('Invalid variable name. Type showOptionNames() for more details.');
    else
        if isstruct(X) == 1
            X.(varname) = val;
        else
            error('Only struct is supported.'); 
        end
    end
end
