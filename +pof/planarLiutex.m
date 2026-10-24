function [R,Wz,lambdaCi,grad] = planarLiutex(U,V,dx,dy,sigmaCells)
% Horizontal 2x2 planar Liutex (NOT a 3D Liutex vector).
% Rows represent increasing y, columns increasing x; dx/dy in meters.
% R = 2 max(|omega_z|/2 - sqrt[((ux-vy)/2)^2+((uy+vx)/2)^2],0).
% Gaussian smoothing applies ONLY to R, not to U/V. The original
% Fig04/Fig05 production code uses a 4-sigma kernel truncation for R.
if nargin<5, sigmaCells=0.8;end
assert(isequal(size(U),size(V)) && ismatrix(U) && all(size(U)>=3), ...
    'U/V must be equally sized two-dimensional arrays >=3 x 3.');
assert(all(isfinite(U(:)))&&all(isfinite(V(:))), ...
    'Vortex detection requires finite velocities; supply the frozen QC field.');
assert(isscalar(dx)&&dx>0&&isfinite(dx)&&isscalar(dy)&&dy>0&&isfinite(dy));
[ux,uy]=gradient(double(U),dx,dy);
[vx,vy]=gradient(double(V),dx,dy);
Wz=vx-uy;
beta=0.5*Wz;
alpha=0.5*hypot(ux-vy,uy+vx);
R=2*max(abs(beta)-alpha,0);
R(beta.^2<=alpha.^2)=0;
R=pof.gaussianSymmetric(R,sigmaCells,sigmaCells,4);
tr=ux+vy;
detG=ux.*vy-uy.*vx;
lambdaCi=0.5*sqrt(max(4*detG-tr.^2,0));
grad=struct('ux',ux,'uy',uy,'vx',vx,'vy',vy);
end
