function [T,mask,threshold] = detectVortices(R,Wz,X,Y,dx,dy,opts)
% Connected planar-Liutex vortex cores, paper baseline P80 of positive R.
% Core area >=4 cells, holes filled, 8-connectivity, 2-cell exclusion.
% GEOMETRIC (equal-area) component centroid as in final Fig05 regionprops.
% Circulation = integral of signed planar vorticity over the detected core.
if nargin<7,opts=struct();end
if ~isfield(opts,'percentile'),opts.percentile=80;end
if ~isfield(opts,'minCells'),opts.minCells=4;end
if ~isfield(opts,'boundaryCells'),opts.boundaryCells=2;end
assert(isequal(size(R),size(Wz),size(X),size(Y)),'Field size mismatch.');
assert(isfinite(dx)&&isfinite(dy)&&isscalar(dx)&&isscalar(dy)&&dx>0&&dy>0, ...
  'Grid spacings dx and dy must be positive lengths in meters.');
assert(isfinite(dx)&&isfinite(dy)&&isscalar(dx)&&isscalar(dy)&&dx>0&&dy>0, ...
  'Grid spacings dx and dy must be positive lengths in meters.');
assert(all(isfinite(R(:)))&&all(isfinite(Wz(:))) && ...
    all(isfinite(X(:)))&&all(isfinite(Y(:))),'Non-finite detection input.');
assert(opts.percentile>0&&opts.percentile<100&&opts.minCells>=1&&opts.boundaryCells>=0);
assert(exist('imfill','file')==2&&exist('bwconncomp','file')==2, ...
    'Image Processing Toolbox required.');
pos=R(R>0);
mask=false(size(R)); threshold=NaN;
if ~isempty(pos)
    threshold=pof.percentileInclusive(pos,opts.percentile);
    mask=imfill(R>=threshold,'holes');
end
CC=bwconncomp(mask,8);
keep=false(size(mask));
rows=zeros(0,7);[ny,nx]=size(R);
for i=1:CC.NumObjects
    ind=CC.PixelIdxList{i};[rr,cc]=ind2sub(size(R),ind);
    if numel(ind)<opts.minCells || any(rr<=opts.boundaryCells) || ...
       any(cc<=opts.boundaryCells) || any(rr>ny-opts.boundaryCells) || ...
       any(cc>nx-opts.boundaryCells)
        continue;
    end
    keep(ind)=true;
    % Do not intensity-weight the centroid: the manuscript source producer
    % uses regionprops(...,'Centroid') (arithmetic geometric centroid).
    xc=mean(X(ind)); yc=mean(Y(ind));
    area=numel(ind)*dx*dy;
    rows(end+1,:)=[xc yc sqrt(area/pi) sum(Wz(ind))*dx*dy ...
        max(R(ind)) numel(ind) area]; %#ok<AGROW>
end
mask=keep;
T=array2table(rows,'VariableNames', ...
    {'CentroidX_m','CentroidY_m','Req_m','Gamma_m2ps','PeakR_1ps','CellCount','CoreArea_m2'});
% The main manuscript reports |Gamma|. Preserve signed circulation too.
T.AbsGamma_m2ps=abs(T.Gamma_m2ps);
end
