function redraw_POF_Fig03_units_smooth_HQ_v2()
% REDRAW_POF_FIG03_UNITS_SMOOTH_HQ_V2
% POF 1024 Fig.3: six representative planar-Liutex fields with detected
% vortex centroids and equivalent-radius circles, from physical m/s VEC/DAT.
% Uses original 0604 experimental fields (NOT synthetic examples).
%
% No dialogs; no new full-dataset vortex analysis; 600 dpi TIFF + PNG +
% vector PDF. 63x63 raw grid -> 5%-crop 57x57; spatial step 3.67712 mm;
% Liutex from gradients, Gaussian smoothing sigma=0.8 cell; positive P80,
% min 4 cells and 2-cell exclusion margin. Raw VEC already exports m/s:
% NEVER re-apply calibration (114.91 um/px)/(50 us).
%
% DISPLAY-ONLY CHANGE vs v1:
% Once the scientific Liutex field, P80 threshold, vortex positions and
% equivalent radii have been calculated on the ORIGINAL 57x57 grid, only
% the displayed color background is refined by cubic interpolation and
% a small Gaussian blur on the FINE display grid. Nothing is re-detected
% or re-measured on the fine grid. Do not use the interpolated image as
% data for any publication statistic.
%
% REPRESENTATIVE FRAME NOTE:
% The exact six original Fig.3 frame IDs are not present in the available
% legacy metadata. By default, each condition picks a reproducible frame
% closest to median count, gamma and Liutex among the FROZEN 2935 frames.
% To preserve exact previously published physical snapshots, insert their
% known frame IDs into cfg.manualFrames once identified (not mandatory).
% Selected file paths and frame IDs are recorded in output manifest CSV.
close all;clc;
base=fileparts(mfilename('fullpath'));
outDir=fullfile(base,'OUTPUT_PUBLICATION');
if ~exist(outDir,'dir'),mkdir(outDir);end
frameCsv=fullfile(base,'Fig04_2935_UNITS_ONLY_corrected.csv');
assert(isfile(frameCsv),'Missing frozen Fig.4 frame index: %s',frameCsv);
F=sortrows(readtable(frameCsv),{'Ra','Frame'});
assert(height(F)==2935 && sum(F.NVortices)==81773,'Invalid Fig.4 frozen catalog.');

cfg.ra=[4.63e8;9.26e8;1.39e9;1.85e9;2.32e9;2.78e9];
cfg.deltaT=[5 10 15 20 25 30];
cfg.prefix={'02run','03run','04run','05run','06run','07run'};
cfg.dx_m=3.67712e-3;
cfg.cmax_1ps=35/2.2982; % proportional to legacy Fig3 maximum, once corrected

% Display refinement only (original PIV statistics remain unchanged):
cfg.displayRefinement=6;  % six display pixels per measured-grid interval
cfg.displaySigmaFine=1.05; % fine-grid Gaussian sigma in display pixels
cfg.displayInterpolation='cubic';
cfg.circleLineWidth=1.0;
cfg.manualFrames=nan(6,1); % set six integer IDs here to force original frames
cfg.rawRoot='G:\0604';
cfg.fallbackRoot='E:\DAT';

assert(exist('imfilter','file')==2 && exist('bwconncomp','file')==2 ...
    && exist('imfill','file')==2 && exist('imgaussfilt','file')==2, ...
    'Fig.3 requires Image Processing Toolbox (imfilter, bwconncomp, imfill, imgaussfilt).');

f=figure('Color','w','Units','centimeters','Position',[2 2 18.9 12.3]);
t=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
map=blueWhiteRed(256);
manifest=cell(6,7);
for k=1:6
    D=F(F.Ra==cfg.ra(k),:);
    assert(~isempty(D),'Missing Ra in reference catalog.');
    if isfinite(cfg.manualFrames(k))
        idx=find(D.Frame==cfg.manualFrames(k),1);
        assert(~isempty(idx),'Requested manual frame is not in frozen 2935 list.');
    else
        medN=median(D.NVortices);medG=median(D.MedianAbsGamma_m2ps);
        medR=median(D.MedianRpeak_1ps);
        sN=max(std(D.NVortices),1);
        sG=max(std(D.MedianAbsGamma_m2ps),1e-7);
        sR=max(std(D.MedianRpeak_1ps),.1);
        score=abs(D.NVortices-medN)/sN + ...
            abs(D.MedianAbsGamma_m2ps-medG)/sG + ...
            abs(D.MedianRpeak_1ps-medR)/sR;
        [~,idx]=min(score);
    end
    frame=double(D.Frame(idx));
    fpath=locateRaw(cfg,k,frame);
    if isempty(fpath)
        error(['Cannot find original Ra %g, frame %d (prefix %s).\n' ...
            'Expected G:\\0604\\%02d\\... or E:\\DAT\\piv_data_RA%d\\...\n' ...
            'Check cfg.rawRoot, cfg.fallbackRoot; no synthetic replacement is allowed.'],...
            cfg.ra(k),frame,cfg.prefix{k},cfg.deltaT(k),k);
    end
    [U,V]=readOriginalVelocity(fpath);
    assert(isequal(size(U),[63 63]),'Expected original 63x63 input.');
    U=U(4:60,4:60); V=V(4:60,4:60); % 5% border, exactly 57x57
    assert(all(isfinite(U(:)))&&all(isfinite(V(:))), ...
        'Nonfinite velocity samples in selected Fig. 3 frame.');
    [R,centers,radii,threshold,count]=detectLiutex(U,V,cfg.dx_m);
    nref=double(D.NVortices(idx));
    tref=double(D.Threshold_1ps(idx));
    assert(count==nref && abs(threshold-tref)<.03,...
        ['Original field does not reproduce the frozen vortex count/P80.\n' ...
         'Data provenance/orientation/units differ: Ra %g, Frame %d, ' ...
         'count %d vs %d, P80 %.4g vs %.4g.'], ...
        cfg.ra(k),frame,count,nref,threshold,tref);

    ax=nexttile(t,k);hold(ax,'on');
    x=(4:60)*cfg.dx_m*1e3; y=(4:60)*cfg.dx_m*1e3;
    % Display-only refinement: no thresholding or core detection on this image.
    % Do not alter R, centers, radii, threshold, or the manifest.
    [xFine,yFine,Rdisplay]=refineLiutexForDisplay( ...
        x,y,R,cfg.displayRefinement,cfg.displaySigmaFine, ...
        cfg.displayInterpolation);
    imagesc(ax,xFine,yFine,Rdisplay);
    for v=1:numel(radii)
        theta=linspace(0,2*pi,180);
        plot(ax,centers(v,1)+radii(v)*cos(theta), ...
            centers(v,2)+radii(v)*sin(theta), ...
            '-','Color',[.2 .23 .28],'LineWidth',cfg.circleLineWidth);
    end
    axis(ax,'image');set(ax,'YDir','normal');
    xlim(ax,[min(x),max(x)]);ylim(ax,[min(y),max(y)]);
    set(ax,'FontName','Times New Roman','FontSize',9.4,'LineWidth',.95,...
        'TickDir','out','Box','on','Layer','top');
    caxis(ax,[0 cfg.cmax_1ps]);colormap(ax,map);
    if k>3,xlabel(ax,'$x\ (\mathrm{mm})$','Interpreter','latex','FontSize',10.5);end
    if mod(k-1,3)==0,ylabel(ax,'$y\ (\mathrm{mm})$','Interpreter','latex','FontSize',10.5);end
    manifest(k,:)={cfg.ra(k),cfg.deltaT(k),frame,count,threshold,nref,fpath};
end
cb=colorbar;cb.Layout.Tile='east';
cb.Label.String='Planar Liutex R (s^{-1})';
cb.Label.FontName='Times New Roman';cb.Label.FontSize=10.5;
cb.FontName='Times New Roman';cb.FontSize=9.4;
cb.Ticks=0:3:15;
% No panel numbers and no top-row Rayleigh labels as previously requested.
png=fullfile(outDir,'Fig03_POF_corrected_smooth_600dpi.png');
tif=fullfile(outDir,'Fig03_POF_corrected_smooth_600dpi.tif');
pdf=fullfile(outDir,'Fig03_POF_corrected_smooth.pdf');
exportgraphics(f,png,'Resolution',600,'BackgroundColor','white');
exportgraphics(f,tif,'Resolution',600,'BackgroundColor','white');
exportgraphics(f,pdf,'ContentType','vector','BackgroundColor','white');
M=cell2table(manifest,'VariableNames', ...
    {'Ra','TemperatureDelta_K','Frame','DetectedNVortices','P80_1ps',...
     'FrozenNVortices','SourceFile'});
writetable(M,fullfile(outDir,'Fig03_smooth_source_frames_and_QC.csv'));
fprintf('Fig.3 display-smoothed export completed (raw statistics unchanged): %s\n',outDir);
end


function [xq,yq,Zq]=refineLiutexForDisplay(x,y,Z,refineFactor,sigmaFine,method)
% DISPLAY ONLY. Original threshold/core identification uses input Z (57x57).
% The interpolation preserves the original physical coordinate extents.
assert(isvector(x) && isvector(y) && ismatrix(Z),'Invalid display grid');
assert(isequal(size(Z),[numel(y),numel(x)]),'Liutex/grid size mismatch');
assert(refineFactor>=1 && refineFactor==round(refineFactor),'Invalid refineFactor');
assert(sigmaFine>=0,'sigmaFine must be nonnegative');

xq=linspace(x(1),x(end),(numel(x)-1)*refineFactor+1);
yq=linspace(y(1),y(end),(numel(y)-1)*refineFactor+1);
[Xq,Yq]=meshgrid(xq,yq);
Zq=interp2(x,y,double(Z),Xq,Yq,method);
% Protect nonnegativity of planar Liutex against cubic interpolation ringing.
Zq(~isfinite(Zq))=0;
Zq=max(Zq,0);
if sigmaFine>0
    Zq=imgaussfilt(Zq,sigmaFine,'Padding','symmetric');
end
end

function path=locateRaw(c,k,fr)
prefix=c.prefix{k}; short=sprintf('%02d',c.deltaT(k));
candidates={fullfile(c.rawRoot,short), ...
  fullfile(c.rawRoot,num2str(c.deltaT(k))), ...
  c.rawRoot,fullfile(c.fallbackRoot,sprintf('piv_data_RA%d',k))};
path='';
for j=1:numel(candidates)
    d=candidates{j}; if ~isfolder(d),continue;end
    stem=sprintf('%s%06d',prefix,fr);
    patterns={[stem,'.dat'],[stem,'*.vec']};
    for p=1:numel(patterns)
        direct=dir(fullfile(d,patterns{p}));
        if ~isempty(direct),path=fullfile(direct(1).folder,direct(1).name);return;end
        sub=dir(fullfile(d,'**',patterns{p}));
        if ~isempty(sub),path=fullfile(sub(1).folder,sub(1).name);return;end
    end
end
end

function [U,V]=readOriginalVelocity(path)
[~,~,ext]=fileparts(path);
fid=fopen(path,'r');assert(fid>0,'Unable to read %s',path);
clean=onCleanup(@()fclose(fid));
if strcmpi(ext,'.dat')
    fgetl(fid);fgetl(fid);
    raw=textscan(fid,'%f%f%f%f','Delimiter',{' ','\t',','}, ...
        'MultipleDelimsAsOne',true,'CollectOutput',true);
    A=raw{1};
else
    head=fgetl(fid);
    assert(contains(head,'U m/s')&&contains(head,'V m/s'), ...
        'VEC header does not explicitly identify physical m/s components: %s',path);
    raw=textscan(fid,'%f%f%f%f%f','Delimiter',',',...
        'MultipleDelimsAsOne',true,'CollectOutput',true);
    A=raw{1};
end
assert(size(A,1)==3969 && size(A,2)>=4,'Original field is not a 63x63 grid.');
X=reshape(A(:,1),63,63)';Y=reshape(A(:,2),63,63)';
U=reshape(A(:,3),63,63)';V=reshape(A(:,4),63,63)';
if median(diff(X(1,:)))<0, U=fliplr(U);V=fliplr(V);X=fliplr(X);Y=fliplr(Y);end
if median(diff(Y(:,1)))<0,U=flipud(U);V=flipud(V);end
% U,V ALREADY m/s: do not multiply by 2.2982.
end

function [R,xy,radius,thr,nKeep]=detectLiutex(U,V,h)
[ux,uy]=gradient(U,h,h);[vx,vy]=gradient(V,h,h);
om=vx-uy;shear=sqrt((ux-vy).^2+(uy+vx).^2);
R=max(0,abs(om)-shear); R(~isfinite(R))=0;
% Exactly matches the legacy explicit rotation-invariant 2D formula.
sigma=.8;rad=ceil(4*sigma);xx=-rad:rad;
g=exp(-.5*(xx/sigma).^2);g=g/sum(g);
R=imfilter(imfilter(R,g,'symmetric','same'),g','symmetric','same');
positive=R(isfinite(R)&R>0);
assert(~isempty(positive),'No positive planar Liutex in selected frame.');
thr=localPctl(positive,80);
mask=imfill(R>=thr,'holes');
C=bwconncomp(mask,8);xy=zeros(0,2);radius=zeros(0,1);
[ny,nx]=size(R);hmm=h*1000;
for i=1:C.NumObjects
    pix=C.PixelIdxList{i};
    if numel(pix)<4,continue;end
    [row,col]=ind2sub([ny nx],pix);
    if any(row<=2|row>=ny-1|col<=2|col>=nx-1),continue;end
    xy(end+1,:)=[mean(col)+3,mean(row)+3]*hmm; %#ok<AGROW>
    radius(end+1,1)=sqrt(numel(pix)*hmm*hmm/pi); %#ok<AGROW>
end
nKeep=numel(radius);
end

function p=localPctl(x,q)
x=sort(x(isfinite(x)));n=numel(x);assert(n>0);
pos=1+(n-1)*q/100; lo=floor(pos);hi=ceil(pos);
if lo==hi,p=x(lo);else,p=x(lo)+(pos-lo)*(x(hi)-x(lo));end
end

function map=blueWhiteRed(N)
anchors=[0 .87 .92 .98;.22 .70 .82 .92;.45 .98 .97 .95;.70 .96 .62 .51;1 .72 .14 .16];
map=interp1(anchors(:,1),anchors(:,2:end),linspace(0,1,N),'linear');
map=max(0,min(1,map));
end
