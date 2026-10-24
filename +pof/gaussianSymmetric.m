function B = gaussianSymmetric(A,sigmaY,sigmaX,truncationSigma)
% Separable, normalized finite Gaussian with symmetric extension.
% Sigma is in PIV vector-grid cells (not camera pixels).
% truncationSigma = 3 for coarse-grained flux (paper boundary rule),
%                   4 for Liutex R smoothing (original Fig04/05 producer).
% This cut-off is a numerical kernel setting, NOT an additional filtering scale.
if nargin<3 || isempty(sigmaX), sigmaX=sigmaY; end
if nargin<4 || isempty(truncationSigma), truncationSigma=3; end
assert(isscalar(truncationSigma) && isfinite(truncationSigma) && ...
    truncationSigma>0, 'Kernel cut-off must be positive.');
validateattributes(A,{'double','single'},{'2d','real','nonempty'});
assert(isscalar(sigmaX)&&isscalar(sigmaY)&&all(isfinite([sigmaY sigmaX])) && ...
    sigmaX>=0 && sigmaY>=0,'Gaussian sigmas must be finite nonnegative scalars.');
if sigmaX==0 && sigmaY==0, B=double(A);return;end
assert(exist('imfilter','file')==2,'Image Processing Toolbox required (imfilter).');
if sigmaX==0, x=0;kx=1; else
    x=-ceil(truncationSigma*sigmaX):ceil(truncationSigma*sigmaX);kx=exp(-0.5*(x/sigmaX).^2);kx=kx/sum(kx);
end
if sigmaY==0,y=0;ky=1; else
    y=(-ceil(truncationSigma*sigmaY):ceil(truncationSigma*sigmaY))';ky=exp(-0.5*(y/sigmaY).^2);ky=ky/sum(ky);
end
B=imfilter(double(A),ky*kx,'symmetric','conv');
end
