function Afsl_shifted = fsl2afni(fslParams,fslMat)
%(fslMat,sz_x,sz_y,sz_z,vox_sz)
% fslParParams are [rot_x,rot_y,rot_z,dx,dy,dz] as specified in mcflirt par
% file, which has origin at voxel (0,0,0).
% afniMotParams are [roll,pitch,yaw,dS,dL,dP], which has origin at center
% of volume, i.e. voxel [size(volume)/2]
% 1) Get FSL affine xfm, Afsl, using makeFslXfmMatrix.m
% 2) Change origin of xfm, Afsl_shifted = Tafni2fsl*Afsl*Tfsl2afni
% where Tfsl2afni=[size(volume)/2]-[0,0,0]
% and Tafni2fsl=[0,0,0]-[size(volume)/2]
% 3) Change from LPI coordinates (FSL) to RAI coordinates (AFNI), via
% negating first 2 rows then first 2 columns. ref: http://afni.nimh.nih.gov/afni/community/board/read.php?1,65457,65472#msg-65472
% 4) Get rotations and translation from matrix
Rfsl=fslParams(1:3);
Tfsl=fslParams(4:6);
Afsl=makeFslXfmMatrix(Tfsl,Rfsl,[1,1,1],'fsl_0000.mat')
filename=['/projects2/udall/pboord/pic/preproc/nogmwm/fsl_RC4101-1/prefiltered_func_data_mcf.mat/',fslMat];
AfslMat=load(filename,'-ascii')
% mv_x=floor(vox_sz*sz_x/2)
% mv_y=floor(vox_sz*sz_y/2)
% mv_z=floor(vox_sz*sz_z/2)
keyboard
mv_x=-115.641151;
mv_y=77.620110;
mv_z=-60.705669;
mv_x=mv_x+2.859;
mv_y=mv_y+(-40.88);
mv_z=mv_z+2.294;
% mv_x=-mv_x;
% mv_y=-mv_y;
% mv_x=-mv_x;
% mv_y=-mv_y;
% mv_z=-mv_z;
% mv_x=-115.641151-(-0.2241);
% mv_y=77.620110-(-19.4613);
% mv_z=-60.705669-32.0361;
%Afsl=Afsl.*[1 1 1 1;1 1 -1 -1;1 -1 1 -1;0 0 0 1];
%Afsl_shifted=[1 0 0 -mv_x;0 1 0 -mv_y;0 0 1 -mv_z;0 0 0 1]*Afsl*[1 0 0 mv_x;0 1 0 mv_y;0 0 1 mv_z;0 0 0 1];
Afsl=Afsl.*[1 1 -1 -1;1 1 -1 -1;-1 -1 1 1;0 0 0 1];
Afsl_shifted=[1 0 0 mv_x;0 1 0 mv_y;0 0 1 mv_z;0 0 0 1]*Afsl*[1 0 0 -mv_x;0 1 0 -mv_y;0 0 1 -mv_z;0 0 0 1];
disp('dS,dL,dP:');
disp([Afsl_shifted(3,4),Afsl_shifted(1,4),Afsl_shifted(2,4)]);
filename='Afsl_shifted';
M=Afsl_shifted;
% Try to open the file
fid = fopen(filename,'w');
if fid == -1
    error(['Could not open for writing: ' filename]);
end
for k = 1:4
    fprintf(fid,[num2str(M(k,:)) '\n']);
end
fclose(fid);