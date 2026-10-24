function [near,transition,background,rho,overlap] = nearFieldMasks(X,Y,centersXY,radii_m,valid,nearCut,bgCut)
% rho(x)=min_i distance(x,center_i)/r_i; overlaps count union only once.
if nargin<6,nearCut=1.5;end
if nargin<7,bgCut=2.0;end
assert(isequal(size(X),size(Y),size(valid)),'X/Y/valid sizes mismatch.');
assert(isscalar(nearCut)&&isscalar(bgCut)&&nearCut>0&&bgCut>nearCut);
assert(size(centersXY,2)==2&&numel(radii_m)==size(centersXY,1));
assert(all(isfinite(centersXY(:)))&&all(isfinite(radii_m))&&all(radii_m>0));
rho=inf(size(X));overlap=zeros(size(X));
for k=1:numel(radii_m)
    thisRho=hypot(X-centersXY(k,1),Y-centersXY(k,2))/radii_m(k);
    rho=min(rho,thisRho);
    overlap=overlap+(thisRho<=nearCut);
end
near=valid & rho<=nearCut;
transition=valid & rho>nearCut & rho<=bgCut;
background=valid & rho>bgCut;
end
