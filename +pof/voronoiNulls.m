function T = voronoiNulls(P,radii,domain,detectDomain,opts)
% Count- and radius-conditioned rectangular Voronoi references.
% CSR_FULL: independent uniform points on full domain (12/frame baseline)
% CSR_DETECT: independent uniform points on detection inset (24/frame)
% SSI_CORE: random trial centroids in detection inset constrained by
%     |pi-pj| >= f*(ri+rj); f=0.5 baseline (24), f=0.75 (12).
% Polygons for each realization clipped to FULL observation rectangle.
% Conditional reference distributions, NOT a generative physical model.
if nargin<5,opts=struct();end
% By default, consume the caller's RNG stream; seed ONCE per complete study.
% Optional opts.seed explicitly resets local stream for a standalone call.
% Do not pass the same seed for each frame of a multi-frame null analysis.
if ~isfield(opts,'nRep'),opts.nRep=[12 24 24 12];end
if ~isfield(opts,'maxProposal'),opts.maxProposal=[0 0 2000 5000];end
assert(numel(radii)==size(P,1)&&size(P,2)==2&&size(P,1)>=2);
assert(all(radii(:)>0)&&all(isfinite(radii(:))));
assert(numel(detectDomain)==4 && ...
       detectDomain(1)>=domain(1)&&detectDomain(2)<=domain(2) && ...
       detectDomain(3)>=domain(3)&&detectDomain(4)<=domain(4));
assert(numel(opts.nRep)==4 && numel(opts.maxProposal)==4);
if isfield(opts,'seed') && ~isempty(opts.seed)
    old=rng;closer=onCleanup(@()rng(old));rng(opts.seed,'twister'); %#ok<NASGU>
end
labels={'CSR_FULL','CSR_DETECT','SSI_CORE_f0p5','SSI_CORE_f0p75'};
factors=[0 0 0.5 0.75];domains=[domain;detectDomain;detectDomain;detectDomain];
modelCol=strings(0,1); repCol=zeros(0,1); okCol=false(0,1); trialCol=zeros(0,1);
cvCol=zeros(0,1); smallCol=zeros(0,1); largeCol=zeros(0,1); errCol=zeros(0,1);
for type=1:4
 for rep=1:opts.nRep(type)
    d=domains(type,:);
    if type<=2
        C=[d(1)+(d(2)-d(1))*rand(numel(radii),1), ...
           d(3)+(d(4)-d(3))*rand(numel(radii),1)];
        success=true;trials=numel(radii);
    else
        C=nan(numel(radii),2);success=true;trials=0;
        for i=1:numel(radii)
            accepted=false;
            for attempt=1:opts.maxProposal(type)
                trials=trials+1;
                q=[d(1)+(d(2)-d(1))*rand(), d(3)+(d(4)-d(3))*rand()];
                if i==1 || all(hypot(C(1:i-1,1)-q(1),C(1:i-1,2)-q(2))>= ...
                                    factors(type)*(radii(1:i-1)+radii(i)))
                    C(i,:)=q;accepted=true;break;
                end
            end
            if ~accepted,success=false;break;end
        end
    end
    cv=NaN;small=NaN;large=NaN;partitionErr=NaN;
    if success
        try
            s=pof.voronoiBounded(C,domain);
            cv=s.VoronoiCV;small=s.SmallCellFraction;
            large=s.LargeCellFraction;partitionErr=s.PartitionAreaError;
        catch
            success=false;
        end
    end
    modelCol(end+1,1)=string(labels{type});repCol(end+1,1)=rep;
    okCol(end+1,1)=success;trialCol(end+1,1)=trials;
    cvCol(end+1,1)=cv;smallCol(end+1,1)=small;
    largeCol(end+1,1)=large;errCol(end+1,1)=partitionErr;
 end
end
T=table(modelCol,repCol,okCol,trialCol,cvCol,smallCol,largeCol,errCol, ...
 'VariableNames',{'Model','Replicate','Success','Trials','VoronoiCV', ...
 'SmallCellFraction','LargeCellFraction','PartitionAreaError'});
% The manuscript's finite-window null results require sufficiently complete
% conditional ensembles. Report an error rather than silently averaging an
% unrepresentative subset of only successful inhibition placements.
if isequal(opts.nRep,[12 24 24 12])
    minimum=[12 24 22 10]; % Original final WP6 per-model acceptance gates.
else
    minimum=ceil(0.9*opts.nRep); % For illustrative custom replicate counts.
end
for typ=1:4
    if nnz(T.Success(T.Model==string(labels{typ}))) < minimum(typ)
        error('pof:incompleteNullEnsemble', ...
           'Finite-window null model %s failed the replicate-completeness gate.',labels{typ});
    end
end
end
