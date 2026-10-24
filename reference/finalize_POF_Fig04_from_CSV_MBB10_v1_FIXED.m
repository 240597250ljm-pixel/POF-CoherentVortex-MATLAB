function finalize_POF_Fig04_from_CSV_MBB10_v1_FIXED(inputDir)
% FINALIZE_POF_FIG04_FROM_CSV_MBB10_V1_FIXED
% Replot Fig. 4 using existing CSVs only (no DAT/VEC processing).
%
% Required input CSV files:
%   Fig04_frame_quality.csv
%   Fig04_vortex_statistics_frame_data.csv
%
% Final definitions:
%   - 2935 accepted frames, QC: median <= 1 px, P95 <= 4 px,
%     and valid-vector fraction >= 0.90.
%   - Panels (a,b): median of frame-wise vortex-property medians,
%     95% moving-block-bootstrap CI (10-frame blocks, B = 3000).
%   - Panels (c,d): median and interquartile range across frames.
%
% Run without arguments (automatically locate/choose the CSV folder):
%   finalize_POF_Fig04_from_CSV_MBB10_v1_FIXED
% Or supply a folder explicitly:
%   finalize_POF_Fig04_from_CSV_MBB10_v1_FIXED('E:\DAT\POF_final_figures\Fig04_vortex_statistics')

clc;
close all;

%% 1. Locate the two existing CSV files
if nargin < 1 || isempty(inputDir)
    scriptDir = fileparts(mfilename('fullpath'));
    preferredDir = 'E:\DAT\POF_final_figures\Fig04_vortex_statistics';
    folders = {scriptDir, preferredDir, pwd};
    inputDir = '';
    for k = 1:numel(folders)
        testDir = folders{k};
        if isfile(fullfile(testDir,'Fig04_frame_quality.csv')) && ...
           isfile(fullfile(testDir,'Fig04_vortex_statistics_frame_data.csv'))
            inputDir = testDir;
            break;
        end
    end
    if isempty(inputDir)
        inputDir = uigetdir(pwd,'Select folder containing the two Fig04 CSV files');
        if isequal(inputDir,0)
            error('No input folder selected.');
        end
    end
end

inputDir = char(inputDir);
qFile = fullfile(inputDir,'Fig04_frame_quality.csv');
fFile = fullfile(inputDir,'Fig04_vortex_statistics_frame_data.csv');
assert(isfile(qFile),'Cannot find %s',qFile);
assert(isfile(fFile),'Cannot find %s',fFile);

outDir = fullfile(inputDir,'Fig04_final_MBB10');
if ~isfolder(outDir)
    mkdir(outDir);
end

Q = readtable(qFile);
F = readtable(fFile);
requiredQ = {'Ra','Frame','Accepted','MedianDisp_px', ...
             'P95Disp_px','ValidFraction'};
requiredF = {'Ra','Frame','NVortices','MedianAbsGamma_m2ps', ...
             'MedianRpeak_1ps','MedianReq_m'};
assert(all(ismember(requiredQ,Q.Properties.VariableNames)), ...
    'Fig04_frame_quality.csv lacks one or more required columns.');
assert(all(ismember(requiredF,F.Properties.VariableNames)), ...
    'Fig04_vortex_statistics_frame_data.csv lacks one or more required columns.');

%% 2. Apply final QC and match frames by (Ra, Frame)
valid = logical(Q.Accepted) & isfinite(Q.MedianDisp_px) & ...
        isfinite(Q.P95Disp_px) & isfinite(Q.ValidFraction) & ...
        Q.MedianDisp_px <= 1.0 & Q.P95Disp_px <= 4.0 & ...
        Q.ValidFraction >= 0.90;

qKey = [double(Q.Ra(valid)),double(Q.Frame(valid))];
fKey = [double(F.Ra),double(F.Frame)];
assert(size(unique(qKey,'rows'),1) == size(qKey,1), ...
    'Repeated (Ra, Frame) keys in accepted quality records.');
assert(size(unique(fKey,'rows'),1) == size(fKey,1), ...
    'Repeated (Ra, Frame) keys in Fig04 frame records.');

keep = ismember(fKey,qKey,'rows');
D = sortrows(F(keep,:),{'Ra','Frame'});
assert(height(D) == 2935, ...
    'Expected 2935 accepted Fig04 frames; found %d.',height(D));
assert(size(qKey,1) == height(D), ...
    'QC accepted-frame count and matched Fig04 count differ.');

raList = unique(D.Ra);
expectedN = [498;499;499;492;470;477];
assert(numel(raList)==6,'Expected six thermal conditions.');

%% 3. Frame-level summaries and moving-block bootstrap
rng(20261020,'twister');
B = 3000;
blockSize = 10;
rows = cell(6,1);

for i = 1:6
    one = D(D.Ra == raList(i),:);
    assert(height(one) == expectedN(i), ...
        'Unexpected frame count for Ra = %.5g: %d.',raList(i),height(one));

    [gm,gl,gh] = mbb_median(one.MedianAbsGamma_m2ps,B,blockSize);
    [rm,rl,rh] = mbb_median(one.MedianRpeak_1ps,B,blockSize);
    [em,e25,e75] = median_iqr(one.MedianReq_m);
    [nm,n25,n75] = median_iqr(one.NVortices);

    rows{i} = table(raList(i),height(one),sum(one.NVortices), ...
        gm,gl,gh,rm,rl,rh,em,e25,e75,nm,n25,n75, ...
        'VariableNames',{'Ra','NFrames','NVortexDetections', ...
        'Gamma_Median','Gamma_MBB95_Low','Gamma_MBB95_High', ...
        'Rpeak_Median','Rpeak_MBB95_Low','Rpeak_MBB95_High', ...
        'Req_Median_m','Req_Q25_m','Req_Q75_m', ...
        'Count_Median','Count_Q25','Count_Q75'});
end

S = vertcat(rows{:});
writetable(D,fullfile(outDir,'Fig04_final_2935_frame_data.csv'));
writetable(S,fullfile(outDir,'Fig04_final_MBB10_summary.csv'));

fprintf('\nFinal Fig04: %d frames, %d figure-level vortex detections.\n', ...
    height(D),sum(S.NVortexDetections));
disp(S(:,{'Ra','NFrames','Gamma_Median','Rpeak_Median','Count_Median'}));
fprintf('Note: this Fig04 export may differ slightly from the master vortex catalog.\n');

%% 4. Draw the four-panel figure
PUB.blue = [0 0.4470 0.7410];
fig = figure('Color','w','Units','centimeters', ...
    'Position',[2 2 18.2 13.8],'PaperPositionMode','auto');

left=.105; right=.975; bottom=.105; top=.955; hgap=.135; vgap=.145;
aw=(right-left-hgap)/2;
ah=(top-bottom-vgap)/2;
P = [left,         bottom+ah+vgap, aw, ah; ...
     left+aw+hgap, bottom+ah+vgap, aw, ah; ...
     left,         bottom,         aw, ah; ...
     left+aw+hgap, bottom,         aw, ah];
ax = gobjects(4,1);
for j=1:4
    ax(j) = axes('Parent',fig,'Position',P(j,:));
end

plotpanel(ax(1),S.Ra,S.Gamma_Median*1e3, ...
    S.Gamma_MBB95_Low*1e3,S.Gamma_MBB95_High*1e3, ...
    '$|\Gamma|\ (10^{-3}\,\mathrm{m^2\,s^{-1}})$','(a)',PUB);
plotpanel(ax(2),S.Ra,S.Rpeak_Median, ...
    S.Rpeak_MBB95_Low,S.Rpeak_MBB95_High, ...
    '$R_{\mathrm{peak}}\ (\mathrm{s^{-1}})$','(b)',PUB);
plotpanel(ax(3),S.Ra,S.Req_Median_m*1e3, ...
    S.Req_Q25_m*1e3,S.Req_Q75_m*1e3, ...
    '$r_{\mathrm{eq}}\ (\mathrm{mm})$','(c)',PUB);
plotpanel(ax(4),S.Ra,S.Count_Median, ...
    S.Count_Q25,S.Count_Q75, ...
    '$N_{\mathrm{v}}\ (\mathrm{frame^{-1}})$','(d)',PUB);

exportgraphics(fig,fullfile(outDir,'Fig04_vortex_statistics_FINAL_MBB10.png'), ...
    'Resolution',600);
exportgraphics(fig,fullfile(outDir,'Fig04_vortex_statistics_FINAL_MBB10.tif'), ...
    'Resolution',600);

fprintf('All files saved to:\n%s\n',outDir);
end

function [med,lo,hi] = mbb_median(y,nBoot,blockSize)
y = y(:);
assert(all(isfinite(y)), ...
    'Nonfinite frame-level data; do not silently alter time order.');
n = numel(y);
assert(n >= blockSize,'Too few frames for a block bootstrap.');
med = median(y);
numBlocks = ceil(n/blockSize);
boot = zeros(nBoot,1);
for b=1:nBoot
    starts = randi(n-blockSize+1,numBlocks,1);
    idx = zeros(numBlocks*blockSize,1);
    for k=1:numBlocks
        dest = (k-1)*blockSize + (1:blockSize);
        idx(dest) = starts(k) + (0:blockSize-1);
    end
    boot(b) = median(y(idx(1:n)));
end
lo = quantile_linear(boot,2.5);
hi = quantile_linear(boot,97.5);
end

function [med,q25,q75] = median_iqr(y)
y = y(isfinite(y));
assert(~isempty(y),'No finite data for the median/IQR.');
med = median(y);
q25 = quantile_linear(y,25);
q75 = quantile_linear(y,75);
end

function v = quantile_linear(y,p)
y = sort(y(isfinite(y)));
n = numel(y);
assert(n>0,'Cannot compute a quantile of empty data.');
if n==1
    v = y(1);
    return;
end
r=1+(n-1)*(p/100);
a=floor(r);
b=ceil(r);
if a==b
    v=y(a);
else
    v=y(a)+(r-a)*(y(b)-y(a));
end
end

function plotpanel(ax,x,y,lo,hi,yLabel,panel,PUB)
hold(ax,'on');
assert(all((y-lo)>=-1e-12 & (hi-y)>=-1e-12), ...
    'A reported interval does not contain the plotted median.');
errorbar(ax,x,y,max(0,y-lo),max(0,hi-y),'-o', ...
    'Color',PUB.blue,'LineWidth',1.7,'MarkerSize',6.8, ...
    'MarkerFaceColor','w','MarkerEdgeColor',PUB.blue,'CapSize',7);
set(ax,'FontName','Times New Roman','FontSize',10.5, ...
    'LineWidth',.9,'TickDir','out','Box','on','Layer','top', ...
    'XScale','log','XGrid','on','YGrid','on', ...
    'TickLabelInterpreter','latex');
xlim(ax,[4e8 3.1e9]);
set(ax,'XTick',[5e8 1e9 2e9], ...
    'XTickLabel',{'$5\times10^8$','$10^9$','$2\times10^9$'});
xlabel(ax,'$Ra_L$','Interpreter','latex','FontSize',12);
ylabel(ax,yLabel,'Interpreter','latex','FontSize',12);
text(ax,-.19,1.045,panel,'Units','normalized', ...
    'FontName','Times New Roman','FontSize',18,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','top','Clipping','off');
yl=ylim(ax);
span=yl(2)-yl(1);
if isfinite(span) && span>0
    ylim(ax,yl+[-.1 .1]*span);
end
end
