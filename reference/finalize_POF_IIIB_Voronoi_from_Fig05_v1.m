function finalize_POF_IIIB_Voronoi_from_Fig05_v1()
% POF SECTION III B: final bounded Voronoi analysis / Fig. 6 / Fig. S5
% Uses the 81,773-vortex catalog exported by the VERIFIED Fig. 5 script,
% which reproduced frozen Fig. 4 on 2935/2935 accepted PIV frames.
% NO DAT reading, NO old Fig10_voronoi_cell_data.csv, NO dialogs.
%
% Outputs: final Fig. 6; separate Fig. S5(a-d) 600 dpi PNG and TIFF;
%          per-frame and per-cell Voronoi CSV; cluster MBB tail CI;
%          endpoint CV MBB contrast; old WP6 comparison audit.
%
% NOTE: Existing WP6 CSR/SSI null values are NOT silently relabeled as
%       results for the slightly different final 81,773-vortex catalog.
%       This script recomputes OBSERVED Voronoi only. Preserve old WP6
%       null sensitivities as historical until consistency is confirmed.

clc; close all;
cfg.fig04CSV = ['E:\DAT\POF_final_figures\Fig04_vortex_statistics\' ...
    'Fig04_final_MBB10\Fig04_final_2935_frame_data.csv'];
cfg.vortexCSV = ['E:\DAT\POF_final_figures\Fig04_vortex_statistics\' ...
    'Fig04_final_MBB10\Fig05_POF_consistent_with_Fig04\Fig05_POF_vortex_samples_2935.csv'];
cfg.wp6Observed = ['E:\DAT\POF_revision_202610\06_voronoi_finite_window\' ...
    'WP6_voronoi_observed_frames.csv'];
cfg.outputDir = 'E:\DAT\POF_final_figures\Fig06_Voronoi_FINAL_81773';
cfg.expectedRa = [4.63e8;9.26e8;1.39e9;1.85e9;2.32e9;2.78e9];
cfg.expectedNFrames = [498;499;499;492;470;477];
cfg.totalFrames = 2935; cfg.totalVortices = 81773;
cfg.gridSize = 63;
cfg.dx_m = 32*114.91e-6;
cfg.cropCells = round(0.05*cfg.gridSize);
cfg.domain_m = [cfg.cropCells*cfg.dx_m, ...
    (cfg.gridSize-1-cfg.cropCells)*cfg.dx_m, ...
    cfg.cropCells*cfg.dx_m, ...
    (cfg.gridSize-1-cfg.cropCells)*cfg.dx_m];
cfg.smallThreshold = 0.5; cfg.largeThreshold = 2;
cfg.bootstrap = 3000; cfg.blockFrames = 10; cfg.seed = 20261022;
cfg.dpi = 600;
cfg.font = 'Times New Roman';
cfg.repRaIdx = [1,3,6];
cfg.cdfX = linspace(0.15,2.8,300);
assert(isfile(cfg.vortexCSV), ['Missing frozen Fig. 5 individual-vortex catalog:\n' ...
    '%s\nRun the previously verified Fig. 5 script (no new PIV work).'],cfg.vortexCSV);
assert(isfile(cfg.fig04CSV),'Missing frozen Fig. 4 frame CSV:\n%s',cfg.fig04CSV);
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end

F=readtable(cfg.fig04CSV);
V=readtable(cfg.vortexCSV);
needF={'Ra','Frame','NVortices','MedianReq_m'};
needV={'Ra','Frame','CentroidX_m','CentroidY_m','Req_m'};
assert(all(ismember(needF,F.Properties.VariableNames)),'Fig. 4 required columns missing.');
assert(all(ismember(needV,V.Properties.VariableNames)),'Fig. 5 required columns missing.');
F=sortrows(F,{'Ra','Frame'});V=sortrows(V,{'Ra','Frame'});
assert(height(F)==cfg.totalFrames,'Frozen Fig. 4 must have 2935 frames.');
assert(height(V)==cfg.totalVortices,'Frozen Fig. 5 must have 81773 vortices.');
assert(isequal(unique(F.Ra),cfg.expectedRa),'Frozen Ra values differ.');
assert(isequal(arrayfun(@(r)nnz(F.Ra==r),cfg.expectedRa),cfg.expectedNFrames), ...
    'Frozen accepted-frame counts differ.');
assert(all(isfinite(V.CentroidX_m)&isfinite(V.CentroidY_m)& ...
    isfinite(V.Req_m)&V.Req_m>0),'Invalid centroid/radius in final Fig. 5.');
keyV=[V.Ra,V.Frame];
[keyFrames,~,group]=unique(keyV,'rows','stable');
assert(size(keyFrames,1)==height(F)&&isequal(keyFrames,[F.Ra,F.Frame]), ...
    'Final Fig. 5 and Fig. 4 frame identities do not match.');
nPer=accumarray(group,1,[height(F),1]);
assert(isequal(nPer,F.NVortices),'Fig. 5 vortex counts differ from frozen Fig. 4.');
medReq=splitapply(@median,V.Req_m,group);
assert(all(abs(medReq-F.MedianReq_m)<1e-7), ...
    'Final Fig. 5 median vortex radius differs from frozen Fig. 4.');

% WP6 EXACT finite-window geometry: rectangular clipping to 5%% cropped PIV.
W=cfg.domain_m;
inside=V.CentroidX_m>=W(1)-1e-8 & V.CentroidX_m<=W(2)+1e-8 & ...
       V.CentroidY_m>=W(3)-1e-8 & V.CentroidY_m<=W(4)+1e-8;
assert(all(inside),'Fig. 5 vortex centroid falls outside cropped Voronoi window.');

vCV=nan(height(F),1);vS=vCV;vL=vCV;vQ95=vCV;vErr=vCV;
area_m2=nan(height(V),1);a=nan(height(V),1);
endRows=[find(diff(group)~=0);height(V)];
startRows=[1;endRows(1:end-1)+1];
fprintf('Reconstructing exact bounded Voronoi cells: %d frames, %d detections.\n', ...
    height(F),height(V));
for k=1:height(F)
    rows=startRows(k):endRows(k);
    P=[V.CentroidX_m(rows),V.CentroidY_m(rows)];
    [areas,~,err]=boundedCells(P,W,false);
    assert(all(isfinite(areas)&areas>0) && err<1e-6, ...
       'Invalid Voronoi partition at Ra=%.0f frame=%d: error %.3g', ...
       F.Ra(k),F.Frame(k),err);
    an=areas/mean(areas);
    area_m2(rows)=areas;a(rows)=an;
    vCV(k)=std(an,0)/mean(an);vS(k)=mean(an<cfg.smallThreshold);
    vL(k)=mean(an>cfg.largeThreshold);vQ95(k)=localPercentile(an,95);
    vErr(k)=err;
    if mod(k,250)==0 || k==height(F)
        fprintf('  %4d/%4d frames, max partition error %.2e\n', ...
            k,height(F),max(vErr(1:k)));
    end
end
R=F(:,{'Ra','Frame','NVortices'});
R.VoronoiCV=vCV;R.SmallCellFraction=vS;R.LargeCellFraction=vL;
R.NormalizedAreaQ95=vQ95;R.PartitionAreaError=vErr;
writetable(R,fullfile(cfg.outputDir,'IIIB_Voronoi_frame_statistics_81773.csv'));
C=V(:,{'Ra','Frame','CentroidX_m','CentroidY_m','Req_m'});
C.VoronoiArea_m2=area_m2;C.NormalizedArea=a;
writetable(C,fullfile(cfg.outputDir,'IIIB_Voronoi_cell_data_81773.csv'));

% WP6 as a REFERENCE VERSION ONLY: does not assert same vortex identities.
if isfile(cfg.wp6Observed)
    W6=readtable(cfg.wp6Observed);
    [hit,idx]=ismember([R.Ra,R.Frame],[W6.Ra,W6.Frame],'rows');
    if all(hit)
        D=table(R.Ra,R.Frame,R.NVortices,W6.NVortices(idx), ...
            R.VoronoiCV,W6.VoronoiCV(idx),R.VoronoiCV-W6.VoronoiCV(idx), ...
            'VariableNames',{'Ra','Frame','FinalN','OldWP6N', ...
            'FinalCV','OldWP6CV','FinalMinusWP6_CV'});
        writetable(D,fullfile(cfg.outputDir,'IIIB_Voronoi_vs_historical_WP6.csv'));
        fprintf('Historical WP6 comparison: %d frames with count differences; max |delta CV| %.4g.\n', ...
           nnz(D.FinalN~=D.OldWP6N),max(abs(D.FinalMinusWP6_CV)));
    else
        warning('Historical WP6 frame coverage differs: audit comparison skipped.');
    end
end

% Frame-level distributions, pooled cell fractions, 10-FRAME CIRCULAR MBB.
rng(cfg.seed,'twister');
B=cfg.bootstrap;
summary=nan(6,13);
for k=1:6
    T=R(R.Ra==cfg.expectedRa(k),:);
    n=T.NVortices;ns=T.SmallCellFraction.*n;nl=T.LargeCellFraction.*n;
    estS=sum(ns)/sum(n);estL=sum(nl)/sum(n);
    bootS=nan(B,1);bootL=nan(B,1);
    for b=1:B
        ii=circularIndices(height(T),cfg.blockFrames);
        bootS(b)=sum(ns(ii))/sum(n(ii));
        bootL(b)=sum(nl(ii))/sum(n(ii));
    end
    summary(k,:)=[cfg.expectedRa(k),height(T),sum(n),median(T.VoronoiCV), ...
        estS,localPercentile(bootS,2.5),localPercentile(bootS,97.5), ...
        estL,localPercentile(bootL,2.5),localPercentile(bootL,97.5), ...
        localPercentile(T.VoronoiCV,25),localPercentile(T.VoronoiCV,75), ...
        mean(T.VoronoiCV)];
end
S=array2table(summary,'VariableNames',{'Ra','NFrames','NVortices', ...
    'FrameMedianCV','Psmall','Psmall_CI95_Lo','Psmall_CI95_Hi', ...
    'Plarge','Plarge_CI95_Lo','Plarge_CI95_Hi','CV_Q25','CV_Q75','FrameMeanCV'});
writetable(S,fullfile(cfg.outputDir,'IIIB_Fig06_pooled_tail_MBB10_summary.csv'));

% Endpoint median CV contrast (independent resampling of two records).
x=R.VoronoiCV(R.Ra==cfg.expectedRa(1));
y=R.VoronoiCV(R.Ra==cfg.expectedRa(6));
delta=median(y)-median(x);bb=nan(B,1);
for b=1:B
    ix=circularIndices(numel(x),cfg.blockFrames);
    iy=circularIndices(numel(y),cfg.blockFrames);
    bb(b)=median(y(iy))-median(x(ix));
end
D=table(median(x),median(y),delta,localPercentile(bb,2.5), ...
    localPercentile(bb,97.5),cfg.blockFrames,B, ...
    'VariableNames',{'T5_MedianCV','T30_MedianCV','DeltaCV', ...
    'DeltaCV_CI95_Lo','DeltaCV_CI95_Hi','BlockFrames','BootstrapReplicates'});
writetable(D,fullfile(cfg.outputDir,'IIIB_CV_endpoint_MBB10.csv'));

% Main Fig. 6 (old script panel definitions, corrected final samples).
plotMainFig06(R,S,cfg);

% Supplementary Fig. S5: representative bounded cells and Poisson delta CDF.
plotSupplementFigS05(R,C,cfg);

fprintf('\nSECTION III B OBSERVATIONAL REPROCESS COMPLETE.\n');
fprintf('Final Fig. 4/5 catalog: %d frames, %d vortex cores.\n',height(F),height(V));
fprintf('T5 to T30 median CV: %.5f -> %.5f, delta %.5f, MBB95 [%.5f, %.5f].\n', ...
    D.T5_MedianCV,D.T30_MedianCV,D.DeltaCV,D.DeltaCV_CI95_Lo,D.DeltaCV_CI95_Hi);
fprintf('IMPORTANT: WP6 CSR/SSI null comparisons require catalog-version review.\n');
fprintf('Results and publication TIFF/PNG: %s\n',cfg.outputDir);
end

function plotMainFig06(R,S,cfg)
colors=[.11 .43 .67; .13 .64 .39];
fig=figure('Color','w','Units','pixels','Position',[100 100 1240 525]);
ax1=axes('Parent',fig,'Position',[0.075 .17 .40 .72]);hold(ax1,'on');
for k=1:6
    x=R.VoronoiCV(R.Ra==cfg.expectedRa(k));
    [y,d]=localDensity(x);
    patch(ax1,[k-0.32*d(:)' k*ones(1,numel(d))], ...
        [y(:)' fliplr(y(:)')],[.42 .67 .83], ...
        'FaceAlpha',.30,'EdgeColor',[.16 .42 .66],'LineWidth',.65);
    low=localPercentile(x,25);m=median(x);hi=localPercentile(x,75);
    rectangle(ax1,'Position',[k-.055 low .11 hi-low], ...
        'FaceColor','w','EdgeColor',[.08 .28 .58],'LineWidth',1.1);
    plot(ax1,[k-.055 k+.055],[m m],'-','Color',[.08 .28 .58],'LineWidth',1.7);
end
plot(ax1,1:6,S.FrameMedianCV,'-d','Color',[.08 .28 .58], ...
    'MarkerFaceColor','w','MarkerSize',5.7,'LineWidth',1.5);
setPaperAxes(ax1,cfg);xlim(ax1,[.5 6.5]);xticks(ax1,1:6);
xticklabels(ax1,compose('%.2f',cfg.expectedRa/1e9));
xlabel(ax1,'$Ra_L$ ($10^9$)','Interpreter','latex');
ylabel(ax1,'Voronoi area $CV_A$','Interpreter','latex');
text(ax1,.02,.98,'(a)','Units','normalized','VerticalAlignment','top', ...
    'FontWeight','bold','FontSize',15,'FontName',cfg.font);
ax2=axes('Parent',fig,'Position',[.565 .17 .38 .72]);hold(ax2,'on');
x=cfg.expectedRa/1e9;dd=.032;
errorbar(ax2,x-dd,S.Psmall,S.Psmall-S.Psmall_CI95_Lo, ...
    S.Psmall_CI95_Hi-S.Psmall,'-o','Color',colors(1,:), ...
    'MarkerFaceColor','w','LineWidth',1.65,'MarkerSize',6,'CapSize',4);
errorbar(ax2,x+dd,S.Plarge,S.Plarge-S.Plarge_CI95_Lo, ...
    S.Plarge_CI95_Hi-S.Plarge,'--s','Color',colors(2,:), ...
    'MarkerFaceColor','w','LineWidth',1.65,'MarkerSize',6,'CapSize',4);
setPaperAxes(ax2,cfg);xlim(ax2,[min(x)-.12 max(x)+.12]);
xticks(ax2,x);xticklabels(ax2,compose('%.2f',x));
xlabel(ax2,'$Ra_L$ ($10^9$)','Interpreter','latex');
ylabel(ax2,'Pooled Voronoi-cell probability','Interpreter','latex');
legend(ax2,{'$P(a<0.5)$','$P(a>2)$'},'Interpreter','latex', ...
    'Location','best','Box','off','FontName',cfg.font);
text(ax2,.02,.98,'(b)','Units','normalized','VerticalAlignment','top', ...
    'FontWeight','bold','FontSize',15,'FontName',cfg.font);
exportgraphics(fig,fullfile(cfg.outputDir,'Fig06_Voronoi_FINAL.png'), ...
    'Resolution',cfg.dpi);
exportgraphics(fig,fullfile(cfg.outputDir,'Fig06_Voronoi_FINAL.tif'), ...
    'Resolution',cfg.dpi);
close(fig);
end

function plotSupplementFigS05(R,C,cfg)
colors=[.05 .35 .65;.05 .60 .68;.28 .68 .48; ...
    .88 .67 .12;.90 .35 .12;.72 .08 .10];
styles={'-','--','-.',':','-','--'};
marks={'o','s','^','d','v','p'};
allA=C.NormalizedArea;
Amax=2.5;Amin=0.2;
for k=1:3
    ra=cfg.expectedRa(cfg.repRaIdx(k));
    T=R(R.Ra==ra,:);
    % A reproducible, near-median representative frame in each condition.
    medN=median(T.NVortices);medCV=median(T.VoronoiCV);
    [~,sel]=min(abs(T.VoronoiCV-medCV)+.004*abs(T.NVortices-medN));
    fr=T.Frame(sel);
    Z=C(C.Ra==ra & C.Frame==fr,:);
    P=[Z.CentroidX_m,Z.CentroidY_m];
    [~,polys,err]=boundedCells(P,cfg.domain_m,true);
    assert(err<1e-6,'Supplementary representative Voronoi cell error.');
    fig=figure('Color','w','Units','pixels','Position',[100 100 670 615]);
    ax=axes('Parent',fig,'Position',[.125 .13 .72 .77]);hold(ax,'on');
    for j=1:height(Z)
        xy=polys{j}*1e3;
        patch(ax,xy(:,1),xy(:,2),Z.NormalizedArea(j), ...
            'EdgeColor',[.20 .20 .20],'LineWidth',.65);
    end
    plot(ax,P(:,1)*1e3,P(:,2)*1e3,'k.','MarkerSize',11);
    xlim(ax,cfg.domain_m(1:2)*1e3);ylim(ax,cfg.domain_m(3:4)*1e3);
    axis(ax,'image');colormap(fig,parula(256));clim(ax,[Amin Amax]);
    setPaperAxes(ax,cfg);grid(ax,'off');
    xlabel(ax,'x (mm)');ylabel(ax,'y (mm)');
    title(ax,sprintf('$Ra_L=%.2f\\times 10^9$ (frame %d)',ra/1e9,fr), ...
        'Interpreter','latex','FontWeight','normal');
    if k==3
        cb=colorbar(ax);cb.FontName=cfg.font;
        ylabel(cb,'$a=A/\langle A\rangle_{\mathrm{frame}}$', ...
            'Interpreter','latex');
    end
    lbl=char('a'+k-1);
    exportgraphics(fig,fullfile(cfg.outputDir, ...
        sprintf('FigS05%s_Voronoi_FINAL.png',lbl)),'Resolution',cfg.dpi);
    exportgraphics(fig,fullfile(cfg.outputDir, ...
        sprintf('FigS05%s_Voronoi_FINAL.tif',lbl)),'Resolution',cfg.dpi);
    close(fig);
end
fig=figure('Color','w','Units','pixels','Position',[100 100 730 580]);
ax=axes('Parent',fig,'Position',[.13 .14 .78 .75]);hold(ax,'on');
x=cfg.cdfX;kshape=3.61;
ref=gammainc(kshape*x,kshape);
for k=1:6
    z=sort(allA(C.Ra==cfg.expectedRa(k)));
    F=histcounts(z,[-inf x inf]); F=cumsum(F(1:end-1))./numel(z);
    plot(ax,x,F-ref,'Color',colors(k,:),'LineStyle',styles{k}, ...
        'LineWidth',1.55,'Marker',marks{k},'MarkerIndices',1:28:numel(x), ...
        'MarkerSize',4.5,'MarkerFaceColor','w');
end
yline(ax,0,'--','Color',[.25 .25 .25],'HandleVisibility','off');
xline(ax,.5,':','Color',[.45 .45 .45],'HandleVisibility','off');
xline(ax,2,':','Color',[.45 .45 .45],'HandleVisibility','off');
setPaperAxes(ax,cfg);xlabel(ax,'Normalized Voronoi area, $a$','Interpreter','latex');
ylabel(ax,'$F_{\mathrm{obs}}(a)-F_{\mathrm{P}}(a)$','Interpreter','latex');
legend(ax,compose('$Ra_L=%.2f\\times 10^9$',cfg.expectedRa/1e9), ...
    'Interpreter','latex','Location','best','Box','off');
xlim(ax,[0 3]);
exportgraphics(fig,fullfile(cfg.outputDir,'FigS05d_Voronoi_deltaCDF_FINAL.png'), ...
    'Resolution',cfg.dpi);
exportgraphics(fig,fullfile(cfg.outputDir,'FigS05d_Voronoi_deltaCDF_FINAL.tif'), ...
    'Resolution',cfg.dpi);
close(fig);
end

function [areas,polys,partitionError]=boundedCells(P,domain,returnPolys)
n=size(P,1);areas=nan(n,1);
if returnPolys,polys=cell(n,1);else,polys={};end
assert(n>=2 && size(unique(P,'rows'),1)==n, ...
    'Voronoi centroids must be finite and distinct.');
rect=[domain(1) domain(3);domain(2) domain(3); ...
      domain(2) domain(4);domain(1) domain(4)];
adj=true(n,n);
if n>=4
    try
        dt=delaunayTriangulation(P); E=edges(dt);
        if ~isempty(E)
            adj=false(n,n);
            for k=1:size(E,1)
                adj(E(k,1),E(k,2))=true;
                adj(E(k,2),E(k,1))=true;
            end
        end
    catch
        adj=true(n,n);
    end
end
for i=1:n
    poly=rect;pi=P(i,:);neighbors=find(adj(i,:));
    neighbors=neighbors(neighbors~=i);
    for j=neighbors
        pj=P(j,:);d=pj-pi;
        c=(pj(1)^2+pj(2)^2-pi(1)^2-pi(2)^2)/2;
        poly=clipHalfplane(poly,d(1),d(2),c);
        if isempty(poly),break;end
    end
    if size(poly,1)<3,areas(i)=0;
    else,areas(i)=polyarea(poly(:,1),poly(:,2));end
    if returnPolys,polys{i}=poly;end
end
windowArea=(domain(2)-domain(1))*(domain(4)-domain(3));
partitionError=abs(sum(areas)-windowArea)/windowArea;
end

function Q=clipHalfplane(P,a,b,c)
if isempty(P),Q=zeros(0,2);return;end
s=a*P(:,1)+b*P(:,2)-c;
Q=zeros(2*size(P,1),2);nQ=0;n=size(P,1);
for i=1:n
    j=mod(i,n)+1;si=s(i);sj=s(j);
    in1=si<=1e-12;in2=sj<=1e-12;
    if in1 && in2
        nQ=nQ+1;Q(nQ,:)=P(j,:);
    elseif in1 && ~in2
        t=si/(si-sj);nQ=nQ+1;Q(nQ,:)=P(i,:)+t*(P(j,:)-P(i,:));
    elseif ~in1 && in2
        t=si/(si-sj);
        nQ=nQ+1;Q(nQ,:)=P(i,:)+t*(P(j,:)-P(i,:));
        nQ=nQ+1;Q(nQ,:)=P(j,:);
    end
end
Q=Q(1:nQ,:);
if size(Q,1)>1
    dd=hypot(diff(Q(:,1)),diff(Q(:,2)));
    Q=Q([true;dd>1e-12],:);
end
end

function idx=circularIndices(n,block)
nBlocks=ceil(n/block); starts=randi(n,nBlocks,1);
idx=zeros(nBlocks*block,1);
for j=1:nBlocks
    rr=(j-1)*block+(1:block);
    idx(rr)=mod(starts(j)-1+(0:block-1),n)+1;
end
idx=idx(1:n);
end

function q=localPercentile(z,p)
z=sort(z(isfinite(z)));n=numel(z);
if n==0,q=nan;return;end
if n==1,q=z;return;end
pos=1+(n-1)*p/100;i=floor(pos);j=ceil(pos);
if i==j,q=z(i);else,q=z(i)+(pos-i)*(z(j)-z(i));end
end

function [y,d]=localDensity(z)
z=z(isfinite(z));
y=linspace(min(z),max(z),140);
if exist('ksdensity','file')==2
    d=ksdensity(z,y);
else
    edges=linspace(min(z),max(z),30);
    h=histcounts(z,edges,'Normalization','pdf');
    c=.5*(edges(1:end-1)+edges(2:end));
    d=interp1(c,h,y,'linear',0);
end
if max(d)>0,d=d/max(d);end
end

function setPaperAxes(ax,cfg)
set(ax,'FontName',cfg.font,'FontSize',11,'LineWidth',.95, ...
    'TickDir','out','Box','on','Layer','top');
grid(ax,'on');ax.GridAlpha=.13;
end
