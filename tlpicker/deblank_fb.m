function Out = deblank_fb(In)
% Removes leading and trailing blanks from a string or cell array of strings
% 
% Usage 
%   Out = deblank_fb(In);
% Inputs
%   In - Input string or cell array of strings
% Outputs
%   Out - Output string or cell array of strings without leading/trailing
%   blanks

% Deblank end of string(s)
Out = deblank(In);

% Deblank start of string
if iscell(Out)
  for i = 1:length(Out(:))
    if any(Out{i}==' ')
      l = find(find((Out{i}==' '))==(1:sum((Out{i}==' '))),1,'last');
      if ~isempty(l)
        Out(i) = cellstr(Out{i}(l+1:end));
      end
    end
  end
else
  if any(Out==' ')
    l = find(find((Out==' '))==(1:sum((Out==' '))),1,'last');
    if ~isempty(l)
      Out = Out(l+1:end);
    end
  end
end

    
    
