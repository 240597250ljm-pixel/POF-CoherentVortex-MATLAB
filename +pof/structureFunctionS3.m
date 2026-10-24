function [r_m,S3,count] = structureFunctionS3(U,V,dx,dy,nBins)
% Frame-level 3rd-order longitudinal structure function on a 2C grid.
% Unique HALF-PLANE integer offsets; min 2*spacing, max quarter of
% smaller cropped field extent; 14 linearly spaced radial bins by default.
% Average cubic projected increments across all valid pairs and offsets.
% NO inference of an inertial range or inverse energy cascade.
if nargin<5,nBins=14;end
assert(isequal(size(U),size(V))&&all(size(U)>=6),'U/V grid mismatch.');
assert(isscalar(nBins)&&nBins>=1&&nBins==floor(nBins));
assert(dx>0&&dy>0&&all(isfinite([dx dy])));
[ny,nx]=size(U);
minR=2*min(dx,dy);
maxR=.25*min((nx-1)*dx,(ny-1)*dy);
assert(maxR>minR,'Grid too small for requested scale range.');
edges=linspace(minR,maxR,nBins+1);
r_m=(edges(1:end-1)+edges(2:end))'/2;
sums=zeros(nBins,1);count=zeros(nBins,1);
maxmx=min(nx-2,ceil(maxR/dx));maxmy=min(ny-2,ceil(maxR/dy));
for my=0:maxmy
    for mx=-maxmx:maxmx
        if my==0 && mx<=0, continue;end
        r=hypot(mx*dx,my*dy);
        if r<minR || r>maxR,continue;end
        b=find(r>=edges(1:end-1) & r<edges(2:end),1,'first');
        if isempty(b) && r==edges(end),b=nBins;end
        if isempty(b),continue;end
        if mx>=0,xa=1:nx-mx;xb=1+mx:nx;
        else,xa=1-mx:nx;xb=1:nx+mx;end
        ya=1:ny-my;yb=1+my:ny;
        du=double(U(yb,xb))-double(U(ya,xa));
        dv=double(V(yb,xb))-double(V(ya,xa));
        dlong=du*(mx*dx/r)+dv*(my*dy/r);
        values=dlong(isfinite(dlong));
        if ~isempty(values)
            sums(b)=sums(b)+sum(values.^3);
            count(b)=count(b)+numel(values);
        end
    end
end
S3=nan(nBins,1);idx=count>0;S3(idx)=sums(idx)./count(idx);
end
