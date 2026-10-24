function result = run_method_tests()
% Purely synthetic analytical checks. This is NOT the manuscript regression.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'io')); % for VEC metadata guard only
n=57;dx=3.67712e-3;
[x,y]=meshgrid((0:n-1)*dx,(0:n-1)*dx);
w=2;u=-w*y;v=w*x;
[R,Wz,ci]=pof.planarLiutex(u,v,dx,dx,0);
assert(max(abs(R(:)-2*w))<1e-10,'Rigid rotation Liutex mismatch.');
assert(max(abs(Wz(:)-2*w))<1e-10,'Rigid rotation vorticity mismatch.');
assert(max(abs(ci(:)-w))<1e-10,'Rigid rotation lambda_ci mismatch.');
[Rstr,~,~]=pof.planarLiutex(u+0.75*x,v+0.75*y,dx,dx,0);
assert(max(abs(Rstr(:)-R(:)))<1e-10, ...
 'Planar Liutex must be invariant under isotropic planar expansion.');
[Rs,~,~]=pof.planarLiutex(2*y,zeros(size(y)),dx,dx,0);
assert(max(abs(Rs(:)))<1e-10,'Pure planar shear must have zero Liutex.');
Vort=pof.detectVortices(R,Wz,x,y,dx,dx); %#ok<NASGU> (uniform rotation touches boundary)
% Deliberately nonuniform R weights across an asymmetric connected region:
% final Fig. 5 uses GEOMETRIC centroids, never Liutex-weighted centroids.
Rshape=zeros(n);Wshape=zeros(n);
Rshape(1,5:40)=1; % many low-positive boundary pixels, excluded as one edge component
Rshape(20:21,20:21)=[10 20;12 18]; Wshape(20:21,20:21)=3;
[Tv,maskV,pv]=pof.detectVortices(Rshape,Wshape,x,y,dx,dx,struct('percentile',70)); %#ok<ASGLU>
assert(height(Tv)==1,'Expected one 4-cell synthetic vortex core.');
assert(abs(Tv.CentroidX_m(1)-mean(x(maskV)))<1e-12 && ...
       abs(Tv.CentroidY_m(1)-mean(y(maskV)))<1e-12, ...
 'Final Fig05 geometric centroid mismatch.');
assert(abs(Tv.AbsGamma_m2ps(1)-12*dx^2)<1e-12, ...
 'Vortex circulation area-integration mismatch.');
% Liutex R smoothing must match the original symmetric 4-sigma kernel.
rr=zeros(n);rr(29,29)=1;
g4=pof.gaussianSymmetric(rr,0.8,0.8,4);
g3=pof.gaussianSymmetric(rr,0.8,0.8,3);
assert(g4(29,33)>0 && g3(29,33)==0, ...
 'Liutex smoothing must retain 4-sigma discrete kernel support.');
[Pi,valid,details]=pof.horizontalFlux(zeros(n),zeros(n),dx,dx,20.43e-3);
assert(all(Pi(valid)==0)&&nnz(valid)==1681,'Expected 41x41 main support for ell=20.43mm.');
% Homogeneous 2C coarse-graining must be Galilean invariant and cubic
% under uniform velocity scaling (tau is quadratic, S linear).
U0=.5*x+3*y;V0=-1.7*x+0.4*y;
[P0,v0]=pof.horizontalFlux(U0,V0,dx,dx,20.43e-3);
[Pshift,~]=pof.horizontalFlux(U0+17,V0-11,dx,dx,20.43e-3);
[Pscale,~]=pof.horizontalFlux(2*U0,2*V0,dx,dx,20.43e-3);
assert(max(abs(Pshift(v0)-P0(v0)))<1e-9, 'Flux Galilean invariance failed.');
assert(max(abs(Pscale(v0)-8*P0(v0)))<1e-9, 'Flux cubic scaling failed.');
% NOTE 1681 holds when ceil(3*20.43mm/(2.355*3.67712mm))=8 cells.
assert(details.validCount==1681,'Gaussian boundary-width mismatch.');
[rr,s3,~]=pof.structureFunctionS3(ones(n),zeros(n),dx,dx,14);
assert(numel(rr)==14&&all(s3(isfinite(s3))==0));
[~,sp]=pof.structureFunctionS3(x,zeros(n),dx,dx,14);
[~,sn]=pof.structureFunctionS3(-x,zeros(n),dx,dx,14);
assert(all(abs(sp(isfinite(sp))+sn(isfinite(sn)))<1e-12), ...
 'Longitudinal cubic increment parity failed.');
[near,tr,bg,rho,over]=pof.nearFieldMasks(x,y,[.10 .10],.015,valid); %#ok<ASGLU>
assert(nnz(near&bg)==0 && nnz(near&tr)==0);
assert(nnz((near|tr|bg)~=valid)==0);
stat=pof.nearFieldStats(-ones(n),near,bg);
assert(abs(stat.negativeContributionRatio-1)<1e-10);
assert(abs(stat.identityError)<1e-10);
sv=pof.strongLiutexStats(R,-ones(n),valid);
assert(sv.A90>0 && abs(sv.negativeFraction-sv.A90)<1e-12);
[vs,areas]=pof.voronoiBounded([.25 .25;.75 .25;.25 .75;.75 .75],[0 1 0 1]);
assert(max(abs(areas-.25))<1e-10 && vs.PartitionAreaError<1e-10);
[point,ci,boot]=pof.movingBlockCI((1:50)',10,60,7); %#ok<ASGLU>
assert(isfinite(point)&&numel(ci)==2&&size(boot,1)==60);
header='VARIABLES = "X mm", "Y mm", "U m/s", "V m/s" MicrosecondsPerDeltaT="50"';
m=pof_parse_vec_header(header);
[~,~,U,V,meta]=pof_normalize_vec_units(zeros(2),zeros(2),ones(2),ones(2),m,114.910004,50); %#ok<ASGLU>
assert(all(U(:)==1)&&all(V(:)==1)&&meta.velocityFactor==1);
Tnull=pof.voronoiNulls([.25 .25;.75 .25;.25 .75;.75 .75],.05*ones(4,1),[0 1 0 1],[.1 .9 .1 .9], ...
  struct('nRep',[1 1 1 1],'seed',9));
assert(height(Tnull)==4 && all(Tnull.Success),'Voronoi conditional references failed.');
result=struct('syntheticMathematicalTests','PASSED',...
    'numericalPaperRegression','NOT_TESTED','fullPIVRecompute','NOT_TESTED');
fprintf('PASS: rotation, shear, finite Gaussian support, S3, masks, negative identity,\n');
fprintf('      geometric centroids, R90, bounded Voronoi, block bootstrap, units.\n');
fprintf('NOT TESTED: 2935-frame regression, 81773 catalog, exact paper figures.\n');
end
