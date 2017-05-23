function vox = mm2vox(mm,hdr,firstIndex)

xfm = [hdr.hist.srow_x; hdr.hist.srow_y; hdr.hist.srow_z; 0 0 0 1];

vox = (xfm\([mm,1]'))';

vox = round(vox(:,1:3));

if firstIndex
    % xfm assumes FSL convention of firstIndex = 0
    vox = vox + 1;
end
