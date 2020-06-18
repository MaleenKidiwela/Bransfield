function x = string2vector(string);
% Breaks up a string into a vector

if all(~(string == ']')) && all(~(string == ']')); 
  string = ['[' string ']'];
end

eval(['x = ' string ';'])