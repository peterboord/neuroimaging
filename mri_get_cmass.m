function cog = mri_get_cmass(imFile)

imS=MRIread(imFile);
far=imS.vol(:);
nx=imS.volsize(1);
ny=imS.volsize(2);
nz=imS.volsize(3);
nxy=nx*ny;
sum=0;
xx=0;
yy=0;
zz=0;
for kk = 0:nz-1
    koff = kk*nxy;
    for jj = 0:ny-1
        joff = koff + jj * nx;
        for ii = 0:nx-1
            val = abs(far(ii+joff+1));
            sum = sum+val;
            xx = xx + val*ii;
            yy = yy + val*jj;
            zz = zz + val*kk;
        end
    end
end
xx = xx/sum;
yy = yy/sum;
zz = zz/sum;
% MRIwrite swaps x & y, so swap back
cog = [yy,xx,zz];

% for( kk=0 ; kk < nz ; kk++ ){
%      koff = kk * nxy ;
%      for( jj=0 ; jj < ny ; jj++ ){
%        joff = koff + jj * nx ;
%        for( ii=0 ; ii < nx ; ii++ ){
%          val = fabs(far[ii+joff]) ;
%          sum += val ;
%          xx  += val * ii ;
%          yy  += val * jj ;
%          zz  += val * kk ;
%        }
%      }
%    }
%    if( flim != im ) mri_free(flim) ;
% 
%    if( sum > 0.0 ){ xx /= sum ; yy /= sum ; zz /= sum ; }
%    else           { xx = 0.5*(nx-1); yy=0.5*(ny-1); zz=0.5*(nz-1); }
