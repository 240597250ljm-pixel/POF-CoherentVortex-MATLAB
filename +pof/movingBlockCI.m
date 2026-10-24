function [point,ci,boot] = movingBlockCI(frameValues,blockLength,nBoot,seed)
% Circular moving-block bootstrap of FRAME statistics (never pixels).
% Columns: frame observables; rows: chronological 500-index timeline or
% accepted-frame sequence. Missing frames remain NaN rather than shifted.
% Uses percentile interval on bootstrap median; 10-frame/3000-draw default.
if nargin<2,blockLength=10;end
if nargin<3,nBoot=3000;end
if nargin<4,seed=20261022;end
assert(isnumeric(frameValues)&&ismatrix(frameValues)&&~isempty(frameValues));
assert(blockLength>=1&&blockLength<=size(frameValues,1)&&blockLength==floor(blockLength));
assert(nBoot>=1&&nBoot==floor(nBoot));
% Restore caller RNG state (avoids changing outside analyses).
old=rng;restore=onCleanup(@()rng(old));rng(seed,'twister'); %#ok<NASGU>
[n,p]=size(frameValues);
point=nan(1,p);ci=nan(2,p);boot=nan(nBoot,p);
for j=1:p
    y=frameValues(:,j);yy=y(isfinite(y));if ~isempty(yy),point(j)=median(yy);end
end
nBlocks=ceil(n/blockLength);
for b=1:nBoot
    start=randi(n,nBlocks,1);
    idx=zeros(nBlocks*blockLength,1);
    for k=1:nBlocks
        ii=(k-1)*blockLength+(1:blockLength);
        idx(ii)=mod(start(k)-1+(0:blockLength-1),n)+1;
    end
    idx=idx(1:n);
    for j=1:p
        a=frameValues(idx,j);a=a(isfinite(a));
        if ~isempty(a),boot(b,j)=median(a);end
    end
end
for j=1:p
    a=boot(:,j);a=a(isfinite(a));
    if ~isempty(a),ci(:,j)=pof.percentileInclusive(a,[2.5 97.5])';end
end
end
