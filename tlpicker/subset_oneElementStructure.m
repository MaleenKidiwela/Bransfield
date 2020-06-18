function struc = subset_oneElementStructure(struc, n, keep)
% Function to subsets vectors that are fields of a 1 element structure
%
% Usage
%   struc = subset_oneElementStructure(struc, n, keep)
% 
% Inputs
%   STRUC - One element structure
%   N     - Number of elements in vectors that are fields of struc or the
%           name of an a scalar field in struc with this value
%   KEEP  - Logical index vector or regular indicies of elements keep in
%           the vectors of struc
%
% Outputs
%   STRUC - Subset structure

if nargin~=3
  error('subset_oneElementStructure requires 3 input arguments')
end

if islogical(keep)
  keep = find(keep);
else
  if all(keep==0 | keep==1)
    keep = find(keep);
  end
end

if ischar(n)
  nold = getfield(struc, n);
  struc = setfield(struc, n, length(keep));
else
  nold = n;
end

fn = fieldnames(struc);
for i = 1:length(fn)
  data = getfield(struc, fn{i});
  if length(data) == nold
    struc = setfield(struc, fn{i}, data(keep));
  end
end