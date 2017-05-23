function imgOut = imMean(img, mask, dim, backImg)

switch dim
    case 'x'
        img = shiftdim(img,1);
        mask = shiftdim(mask,1);
    case 'y'
        img = shiftdim(img,2);
        mask = shiftdim(mask,2);
    case 'z'
    otherwise
        error('invalid dimension');
end

if nargin > 3
    imgOut = backImg;
else
    imgOut = max(img(:))*ones(size(img,1),size(img,2));
end
    
for j = 1:size(img,1)
    for k = 1:size(img,2)
        zPos = mask(j,k,:)>0;
        if sum(zPos(:)) > 0
            imgOut(j,k) = mean(img(j,k,zPos),3);
        end
    end
end

switch dim
    case 'y'
        imgOut = shiftdim(imgOut,1);
    otherwise
end

% function img = imMean(imgIn, mask, xIn, yIn)
% 
% img = 0.4*max(imgIn(:))*ones(size(imgIn,1),size(imgIn,2));
% img(xIn,:) = 0;
% img(:,yIn) = 0;
% for x = 1:size(imgIn,1)
%     for y = 1:size(imgIn,2)
%         zPos = mask(x,y,:)>0;
%         if sum(zPos(:)) > 0
%             img(x,y) = mean(imgIn(x,y,zPos),3);
%         end
%     end
% end