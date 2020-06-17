function [index,status] = index_eventid(id1,id2)
% Finds the indexes of elements in the the first input those in the second
%
% Usage
%  [index,status] = index_eventid(id1,id2)
% Inputs
%  id1 - A vector of ids (e.g., event ids for shots)
%  id2 - A second vector/matrix of ids (e.g., event ids for arrivals)
% Outputs
%  index_id - Index of id1 value matching each element in id2
%             0 implies no match (e.g., arrival is for unknown shot)
%             -N implies >1 match (e.g., arrival if for >1 shot)
%  status   - 0 All id2 values are matched once
%             1 Some id2 values are unmatched
%             2 Some id2 values are matched more than once

index = zeros(length(id2(:)),1);

for i = 1:length(id2);
  k = find(id2(i)==id1(:));
  status = 0;
  if length(k)==1
    index(i) = k;
  elseif ~~isempty(k);
    index(i) = 0;
    status = max(status,1);
  else
    index(i) = -length(k);
    status = max(status,2);
  end
end

index = reshape(index,size(id2));  

end