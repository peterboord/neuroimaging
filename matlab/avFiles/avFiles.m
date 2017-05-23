function avFiles(varargin)

if isdeployed
    echo off
else
    dbstop if error
end

if nargin == 0
    disp('Usage:');
    disp('avFiles file1 file2 ... outFile');
    return
elseif nargin == 1
    error('not enough inputs');
end
outFile = varargin{end};
varargin(end) = [];
nrFiles = numel(varargin);
system(['cp ',varargin{2},' ',outFile]);

x = [];
for fileNr = 1:nrFiles
    dataIn = dlmread(varargin{fileNr},' ',1,1);
    x = cat(3,x,reshape(dataIn,[size(dataIn),1]));
end
avData = mean(x,3);

outFid = fopen(outFile,'w');
inFid = fopen(varargin{1});
fprintf(outFid,'%s',fgets(inFid));
while ~feof(inFid)
    row = fgets(inFid);
    fprintf(outFid,'%s',row);
end
fclose('all');
