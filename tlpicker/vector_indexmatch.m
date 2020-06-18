function index = vector_indexmatch(vec1,vec2)
% Indicies of elements in the first vector equaling those in the second
%
% Usage
%  index = vector_indexmatch(vec1,vec2)
% Inputs
%  vec1 - A vector of numbers (e.g., event ids for shots)
%  vec2 - A second vector of numbers (e.g., event ids for arrivals)
% Outputs
%  index -  For each element in vec2 gives the index of the element in vec1
%           with the same value
%           0   - No match
%           -N  - More than one match 

index = zeros(size(vec2));

for i = 1:length(vec2)
  j = find(vec1 == vec2(i));
  if length(j) == 1;
    index(i) = j;
  elseif length(j) > 1;
    index(i) = - length(j);
  end
end