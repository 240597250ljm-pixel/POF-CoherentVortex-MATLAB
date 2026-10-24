function [Pi,valid,details] = horizontalFlux(U,V,dx,dy,ell_m,opts)
% Horizontal two-component coarse-grained diagnostic -tau_ij S_ij.
% ell is Gaussian FWHM in meters; use 3-sigma symmetric convolution.
% Explicitly NOT a complete 3D turbulent kinetic-energy flux.
% opts.preSmoothSigmaCells defaults to zero per stated coarse-grain method.
% Historical analysis scripts additionally smooth U,V at sigma=0.8 cells;
% if needed set that option EXPLICITLY and disclose the difference.
if nargin<6,opts=struct();end
if ~isfield(opts,'preSmoothSigmaCells'),opts.preSmoothSigmaCells=0;end
assert(isequal(size(U),size(V)) && all(size(U)>=3),'U/V grid mismatch.');
assert(all(isfinite(U(:)))&&all(isfinite(V(:))),'Pass finite QC velocity grid.');
assert(isscalar(ell_m)&&isfinite(ell_m)&&ell_m>0);
assert(all(isfinite([dx dy]))&&dx>0&&dy>0);
U=double(U)-mean(U(:));V=double(V)-mean(V(:));
if opts.preSmoothSigmaCells>0
    warning('pof:historicalPrefilter', ...
       ['An additional prefilter on U/V was requested. ' ...
        'It is NOT specified in manuscript Sec. II E; disclose it separately.']);
    U=pof.gaussianSymmetric(U,opts.preSmoothSigmaCells, ...
                             opts.preSmoothSigmaCells,4);
    V=pof.gaussianSymmetric(V,opts.preSmoothSigmaCells, ...
                             opts.preSmoothSigmaCells,4);
end
sigma_m=ell_m/(2*sqrt(2*log(2)));
sx=sigma_m/dx;sy=sigma_m/dy;
ub=pof.gaussianSymmetric(U,sy,sx);
vb=pof.gaussianSymmetric(V,sy,sx);
txx=pof.gaussianSymmetric(U.*U,sy,sx)-ub.^2;
txy=pof.gaussianSymmetric(U.*V,sy,sx)-ub.*vb;
tyy=pof.gaussianSymmetric(V.*V,sy,sx)-vb.^2;
[ux,uy]=gradient(ub,dx,dy);[vx,vy]=gradient(vb,dx,dy);
Pi=-(txx.*ux + txy.*(uy+vx) + tyy.*vy);
marginY=ceil(3*sy);marginX=ceil(3*sx);
valid=false(size(Pi));
if size(Pi,1)>2*marginY && size(Pi,2)>2*marginX
    valid((1+marginY):(end-marginY),(1+marginX):(end-marginX))=true;
end
Pi(~valid)=NaN;
details=struct('ellFWHM_m',ell_m,'sigma_m',sigma_m,...
   'sigmaGrid',[sy sx],'marginGrid',[marginY marginX],...
   'preSmoothSigmaCells',opts.preSmoothSigmaCells,'validCount',nnz(valid),...
   'isFull3DEnergyFlux',false);
end
