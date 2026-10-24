function [matchedNear,matchedBg,summary] = matchGradientStrain(U,V,dx,dy,near,bg,opts)
% Within-frame 2-scalar matching of |grad u|_F and |sym grad u|_F.
% Six-by-six quantile strata of log magnitudes, greedy without replacement,
% standardized-log 0.60 caliper (R4-C9 historical analysis settings).
% Controls two amplitudes, not causal attribution or the full tensor.
if nargin<7,opts=struct();end
if ~isfield(opts,'nBins'),opts.nBins=6;end
if ~isfield(opts,'caliper'),opts.caliper=0.60;end
if ~isfield(opts,'seed'),opts.seed=20261018;end
assert(isequal(size(U),size(V),size(near),size(bg)));
assert(~any(near(:)&bg(:)));
[ux,uy]=gradient(double(U),dx,dy);
[vx,vy]=gradient(double(V),dx,dy);
G=sqrt(ux.^2+uy.^2+vx.^2+vy.^2);
S=sqrt(ux.^2+vy.^2+0.5*(uy+vx).^2);
pool=(near|bg)&isfinite(G)&isfinite(S);
matchedNear=false(size(near));matchedBg=false(size(bg));
summary=struct('nPairs',0,'coverageNear',0,'coverageBackground',0,...
 'nCommonBins',0,'matchedSMDGradient',NaN,'matchedSMDStrain',NaN);
if nnz(near)<2 || nnz(bg)<2,return;end
g0=max(median(G(pool)),eps)*1e-3;
s0=max(median(S(pool)),eps)*1e-3;
lG=log(G+g0);lS=log(S+s0);
muG=mean(lG(pool));sdG=std(lG(pool));
muS=mean(lS(pool));sdS=std(lS(pool));
if sdG<=eps,sdG=1;end
if sdS<=eps,sdS=1;end
zG=(lG-muG)/sdG;zS=(lS-muS)/sdS;
pc=(1:opts.nBins-1)*100/opts.nBins;
gc=pof.percentileInclusive(lG(pool),pc);sc=pof.percentileInclusive(lS(pool),pc);
gbin=1+sum(bsxfun(@gt,lG(:),gc),2);
sbin=1+sum(bsxfun(@gt,lS(:),sc),2);
gbin=reshape(gbin,size(U));sbin=reshape(sbin,size(U));
old=rng;cleanup=onCleanup(@()rng(old));rng(opts.seed,'twister'); %#ok<NASGU>
for i=1:opts.nBins
 for j=1:opts.nBins
    ni=find(near & pool & gbin==i & sbin==j);
    bi=find(bg & pool & gbin==i & sbin==j);
    if isempty(ni)||isempty(bi),continue;end
    summary.nCommonBins=summary.nCommonBins+1;
    smallIsNear=(numel(ni)<=numel(bi));
    if smallIsNear,small=ni;large=bi;else,small=bi;large=ni;end
    small=small(randperm(numel(small)));
    available=true(numel(large),1);
    for k=1:numel(small)
        candidates=find(available);
        if isempty(candidates),break;end
        distances=(zG(small(k))-zG(large(candidates))).^2 + ...
             (zS(small(k))-zS(large(candidates))).^2;
        [minD,id]=min(distances);
        if isfinite(minD)&&minD<=opts.caliper^2
            partner=large(candidates(id));available(candidates(id))=false;
            if smallIsNear
                matchedNear(small(k))=true;matchedBg(partner)=true;
            else
                matchedBg(small(k))=true;matchedNear(partner)=true;
            end
        end
    end
 end
end
summary.nPairs=nnz(matchedNear);
summary.coverageNear=summary.nPairs/nnz(near);
summary.coverageBackground=nnz(matchedBg)/nnz(bg);
summary.matchedSMDGradient=smd(G(matchedNear),G(matchedBg));
summary.matchedSMDStrain=smd(S(matchedNear),S(matchedBg));
end
function z=smd(a,b)
if numel(a)<2||numel(b)<2,z=NaN;return;end
pooled=sqrt((var(a)+var(b))/2);
if pooled<=eps,z=0;else,z=(mean(a)-mean(b))/pooled;end
end
