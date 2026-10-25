# One-frame reproducible MATLAB usage (after original QC)

```matlab
F=pof_read_vec_field(vecPath,114.910004,50,0.05);
U=F.U; V=F.V; X=F.X; Y=F.Y; dx=F.dx; dy=F.dy;
% U,V: 57x57 cropped calibrated m/s velocity arrays on a uniform grid.
% dx,dy in meters, X,Y meshgrid coordinates matching U/V orientation.
frameResult=pof.analyzeFrame(F);  % One-frame reference pipeline
[R,wz,lambdaCi]=pof.planarLiutex(U,V,dx,dy,0.8);
[vortexCatalog,coreMask,p80]=pof.detectVortices(R,wz,X,Y,dx,dy);
[piell,valid,details]=pof.horizontalFlux(U,V,dx,dy,20.43e-3);
[r,s3,count]=pof.structureFunctionS3(U,V,dx,dy,14);
[near,shell,bg,rho,overlap]=pof.nearFieldMasks( ...
    X,Y,[vortexCatalog.CentroidX_m vortexCatalog.CentroidY_m], ...
    vortexCatalog.Req_m,valid,1.5,2);
frameStats=pof.nearFieldStats(piell,near,bg);
strong=pof.strongLiutexStats(R,piell,valid);
% Within-frame two-amplitude matching uses U,V on the same valid domain.
[mn,mb,matchInfo]=pof.matchGradientStrain(U,V,dx,dy,near,bg);
% Run circular MBB on FRAME-BY-FRAME statistics, not on individual pixels.
% Each row represents a frame in acquisition order, NaN for rejected frames.
[medianRatio,CI]=pof.movingBlockCI(activityRatioPerFrame,10,3000,20261022);
```

Default `horizontalFlux`/`analyzeFrame` behavior includes manuscript-specified U/V pre-smoothing (σ=0.8 grid cells), followed by a separate FWHM-ell coarse graining. Set `struct('preSmoothSigmaCells',0)` explicitly to examine the no-pre-smoothing sensitivity case. The main paper's full producer also requires original accepted-frame masks and provenance information; this tutorial does not generate a manuscript-numbered figure.
