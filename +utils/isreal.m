function b = isreal(z)
% MYISREAL Checks if a number or array is real.
%    b = isreal(z)
%    This function performs the following computation:
%    b = (~isnan(z))&(~isinf(z))&isreal(z);
b = (~isnan(z))&(~isinf(z))&isreal(z);
end

