function sortIndex=get_sortIndex(traceMetaData,sortStruc);

fieldNames = fieldnames(traceMetaData);

for i = length(sortStruc):-1:1  
  name = sortStruc(i).name;
  if isempty(name)
    sortStruc(i).index = (1:traceMetaData.ntrace)';
  else
    ascend = true;
    if name(1) == '+'
      name = name(2:end);
    elseif name(1) == '-';
      ascend = false;
      name = name(2:end);
    end 
    k = find(strcmpi(fieldNames,name));
    if isempty(k)
      disp(['Sort name ''' sortStruc(i).name ''' not recognized - Ignored']);
      sortStruc(i).index = (1:traceMetaData.ntrace)';
    else
      eval(['a = traceMetaData.' cell2mat(fieldNames(k)) ';']);
      if length(a)~=traceMetaData.ntrace;
        disp(['Sort name ''' name  ...
          ''' is a valid field but does not have ntrace elements']);
        sortStruc(i).index = (1:traceMetaData.ntrace)';
      else
        for j = length(sortStruc):-1:i+1
          a = a(sortStruc(j).index);
        end
        if ~iscell(a)
          if ascend
            [~,sortStruc(i).index] = sort(a,'ascend');
          else
            [~,sortStruc(i).index] = sort(a,'descend');
          end
        else
          [~,sortStruc(i).index] = sort(a);
          if ~ascend
            sortStruc(i).index = sortStruc(i).index(end:-1:1);
            disp('Warning - due to matlab sort limitations for cell arrays the')
            disp(['descend sort for field ''' name ''' will reverse earlier sorts'])
          end
        end
      end
    end
  end
end

% Apply sorting
sortIndex = (1:traceMetaData.ntrace)';
for i=length(sortStruc):-1:1
  sortIndex = sortIndex(sortStruc(i).index);
end 