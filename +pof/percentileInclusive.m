function q = percentileInclusive(x,p)
% Stable linear-interpolation percentile: r=1+(n-1)*p/100.
% This intentionally matches finite_percentile() of the paper source,
% avoiding MATLAB prctile default-algorithm differences by version.
assert(isnumeric(x)&&isnumeric(p)&&all(isfinite(p(:)))&&...
    all(p(:)>=0)&&all(p(:)<=100));
x=sort(double(x(isfinite(x))));
shp=size(p);p=double(p(:));q=nan(size(p));
if isempty(x),q=reshape(q,shp);return;end
if numel(x)==1,q(:)=x;else
 for k=1:numel(p)
    r=1+(numel(x)-1)*p(k)/100; lo=floor(r);hi=ceil(r);
    if lo==hi,q(k)=x(lo);else
       q(k)=x(lo)+(r-lo)*(x(hi)-x(lo));
    end
 end
end
q=reshape(q,shp);
end
