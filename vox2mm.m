function mm = vox2mm(vox,hdr,firstIndex)

if firstIndex
    % xfm assumes FSL convention of firstIndex = 0
    vox = vox - 1;
end

xfm = [hdr.hist.srow_x; hdr.hist.srow_y; hdr.hist.srow_z];

vox = cat(2,vox,ones(size(vox,1),1,class(vox)));

mm = (xfm*vox')';