function tf = vector_match(vec1,vec2)
% T/F output indicate if values in the first vector equal any in the second
%
% Usage
%  tf = index_matched(vec1,vec2)
% Inputs
%  vec1 - A vector of numbers (e.g., event ids for shots)
%  vec2 - A second vector of numbers (e.g., event ids for arrivals)
% Outputs
%  tf -  Logical vector which is true for elements in vec1 that equal at
%  least one element in vec2

tf = false(size(vec1));

for i = 1:length(vec2)
  tf = tf | (vec1 == vec2(i));
end