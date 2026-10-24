function out = analyzeFrame(F,opts)
% One-frame complete PAPER-METHOD reference calculation from calibrated
% cropped arrays F.X,F.Y,F.U,F.V,F.dx,F.dy (SI units).
% This is NOT a manuscript figure reproducer; no raw PIV QC is performed.
if nargin<2,opts=struct();end
if ~isfield(opts,'ell_m'),opts.ell_m=20.43e-3;end
if ~isfield(opts,'preSmoothSigmaCells'),opts.preSmoothSigmaCells=0;end
needed={'X','Y','U','V','dx','dy'};
for k=1:numel(needed)
    assert(isfield(F,needed{k}),'F missing %s',needed{k});
end
[R,Wz,lambdaCi]=pof.planarLiutex(F.U,F.V,F.dx,F.dy,0.8);
[T,cores,P80]=pof.detectVortices(R,Wz,F.X,F.Y,F.dx,F.dy);
[Pi,valid,filterInfo]=pof.horizontalFlux(F.U,F.V,F.dx,F.dy, ...
    opts.ell_m,struct('preSmoothSigmaCells',opts.preSmoothSigmaCells));
[r,S3,count]=pof.structureFunctionS3(F.U,F.V,F.dx,F.dy,14);
[near,shell,bg,rho,overlap]=pof.nearFieldMasks(F.X,F.Y, ...
    [T.CentroidX_m T.CentroidY_m],T.Req_m,valid,1.5,2.0);
regional=pof.nearFieldStats(Pi,near,bg);
strong=pof.strongLiutexStats(R,Pi,valid);
out=struct('R',R,'Wz',Wz,'lambdaCi',lambdaCi,'cores',cores,...
    'catalog',T,'thresholdP80',P80,'Pi',Pi,'filterInfo',filterInfo,...
    'valid',valid,'r_m',r,'S3',S3,'S3PairCount',count,...
    'near',near,'transition',shell,'background',bg,'rho',rho,...
    'overlap',overlap,'regional',regional,'strongLiutex',strong,...
    'isCompleteThreeDimensionalEnergyFlux',false);
end
