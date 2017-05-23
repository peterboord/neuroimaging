#!/usr/bin/octave -qf
# usage:
# oTan2.m yFile xFile outFile
# outFile = atan2(yFile,xFile)

arg_list = argv();
yFile = arg_list{1};
xFile = arg_list{2};
outFile = arg_list{3};
if nargin < 3
	echo "usage:"
	echo "oTan2.m yFile xFile outFile"
	echo "outFile = atan2(yFile,xFile)"
	exit
end
y = MRIread(yFile);
x = MRIread(xFile);




corrxy = corr(x,y);
nrRows=size(corrxy,1);
nrCols=size(corrxy,2);
if nrRows == 1
	printf(formatString, corrxy);
else
	for rowNr = 1:nrRows
		for colNr = 1:nrCols-1
			printf(formatString,corrxy(rowNr,colNr));
			printf(",");
		end
		printf([formatString,"\n"],corrxy(rowNr,nrCols));
	end
end
