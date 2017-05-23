function correl(varargin)

if isdeployed
    echo off
else
    dbstop if error
end

if nargin == 0
    disp('Usage:');
    disp('correl file1 file2 ... [-rowName rowName] [-colName colName] [-headers] [-fisher]');
    return
elseif nargin == 1
    error('not enough inputs');
end
rowNameIdx = find(strcmp('-rowName',varargin));
if isempty(rowNameIdx)
    rowName = '';
else
    rowName = varargin{rowNameIdx + 1};
    varargin(rowNameIdx:rowNameIdx+1) = [];
end
colNameIdx = find(strcmp('-colName',varargin));
if isempty(colNameIdx)
    colName = '';
else
    colName = varargin{colNameIdx + 1};
    varargin(colNameIdx:colNameIdx+1) = [];
end
headers = find(strcmp('-headers',varargin));
if headers
    varargin(headers) = [];
end
fisher = find(strcmp('-fisher',varargin));
if fisher
    varargin(fisher) = [];
end
nrFiles = numel(varargin);
x = [];
for fileNr = 1:nrFiles
    fileContents = dlmread(varargin{fileNr});
    if ~isempty(x) && size(x,1)~=numel(fileContents)
        error('file size mismatch');
    end
    x =cat(2,x,fileContents(:));
end
if headers
    if isempty(colName)
        rowText = '';
    else
        rowText = colName;
    end
    for row = 1:nrFiles-1
        for col = row+1:nrFiles
            rowText = [rowText,' ',num2str(row),'x',num2str(col)];
        end
    end
    disp(rowText);
end
if isempty(rowName)
    rowText = '';
else
    rowText = rowName;
end
for row = 1:nrFiles-1
    for col = row+1:nrFiles
        if fisher
            rowText = [rowText,' ',num2str(atanh(corr(x(:,row),x(:,col))))];
        else
            rowText = [rowText,' ',num2str(corr(x(:,row),x(:,col)))];
        end
    end
end
disp(rowText);
end