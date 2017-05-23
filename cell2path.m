function path = cell2path(stringCell)

path = stringCell{1};
for i = 2:numel(stringCell)
    path = [path,filesep,stringCell{i}]; %#ok<AGROW>
end