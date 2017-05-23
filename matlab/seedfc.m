function seedfc(imPathIn,imPathOut,xFSL,yFSL,zFSL,sz) % Note yFS,xFS = xFSL,yFSL !!!
%%
% N.B. VERY IMPORTANT!!!!
% Note that FS & FSL have x & y swapped!!!!
%%
if nargin == 0
    disp('seedfc(imPathIn,imPathOut,xFSL,yFSL,zFSL)');
end
if nargin < 6
    sz = 1;
else
    sz = str2double(sz);
end

xFS = str2double(yFSL);
yFS = str2double(xFSL);
zFS = str2double(zFSL);

% tmpS=MRIread(imPathOut);
% disp(tmpS.vol(xFS,yFS,z,1));
% blankS.vol = tmpS.vol; 
% MRIwrite(blankS,imPathOut);

imS = MRIread(imPathIn);
% add 1 for FSL 0-based numbering
% imInS.vol = fc(imInS,[xFS+1,yFS+1,zFS+1],sz);
func2d = reshape(imS.vol,[],imS.nframes)';
imS.vol = reshape(corr(func2d,squeeze(imS.vol(xFS,yFS,zFS,:))),imS.volsize);
% imInS.nframes = (2*sz+1)^3;
% imInS.dim(5) = imInS.nframes;
MRIwrite(imS,imPathOut,'float');
end

function fcOut = fc(imS,coord,sz)
fcOut = zeros([imS.volsize,(2*sz+1)^3]);
func2d = reshape(imS.vol,[],imS.nframes)';
idx = 0;
for i=coord(1)-sz:coord(1)+sz
    for j=coord(2)-sz:coord(2)+sz
        for k=coord(3)-sz:coord(3)+sz
            idx = idx + 1;
            fcOut(:,:,:,idx) = reshape(corr(func2d,squeeze(imS.vol(i,j,k,:))),imS.volsize);
        end
    end
end
end
