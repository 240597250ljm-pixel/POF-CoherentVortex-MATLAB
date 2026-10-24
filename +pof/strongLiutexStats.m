function s = strongLiutexStats(R,Pi,valid)
% SI S4: R90 positive planar Liutex threshold on the valid 2D area.
% A90 = area(strong)/area(valid), NOT fixed 0.1.
% F90 = negative magnitude in strong region / all negative magnitude.
% Enrichment = F90/A90, defined only for nonempty negative magnitudes.
assert(isequal(size(R),size(Pi),size(valid)));
v=valid & isfinite(R)&isfinite(Pi);
pos=R(v & R>0);
s=struct('R90',NaN,'A90',NaN,'negativeFraction',NaN,...
  'enrichment',NaN,'nValid',nnz(v));
if isempty(pos)||~any(v(:)),return;end
s.R90=pof.percentileInclusive(pos,90);
strong=v & R>=s.R90;
s.A90=nnz(strong)/nnz(v);
w=max(-Pi,0);w(~v)=0;
total=sum(w(:));
if total>0
    s.negativeFraction=sum(w(strong))/total;
    s.enrichment=s.negativeFraction/s.A90;
end
end
