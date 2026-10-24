function F = pof_read_vec_field(filename,umPerPixel,pulseUs,cropFraction)
% Strict native Tecplot POINT VEC reader, no legacy second velocity scaling.
% Reads X,Y,U,V,CHC from a VEC with explicit variable-unit header and
% ZONE I=... J=... F=POINT. Does not infer QC from CHC codes.
% The returned physical fields have rows increasing y and columns increasing x.
% Specify calibration explicitly; 0605 pixel-displacement needs same-day value.
if nargin<4,cropFraction=0.05;end
assert(isfile(filename),'VEC does not exist: %s',filename);
assert(cropFraction>=0&&cropFraction<0.45,'Invalid crop fraction.');
meta=pof_parse_vec_header(filename);
fid=fopen(filename,'rt');assert(fid>0,'Could not open VEC.');
closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
zone='';
for k=1:64
  line=fgetl(fid);if ~ischar(line),break;end
  if ~isempty(regexpi(line,'\bZONE\b','once'))
      zone=line;break;
  end
end
assert(~isempty(zone),'VEC ZONE header not found.');
mi=regexp(zone,'\bI\s*=\s*(\d+)','tokens','once');
mj=regexp(zone,'\bJ\s*=\s*(\d+)','tokens','once');
assert(~isempty(mi)&&~isempty(mj),'VEC I,J dims not explicit.');
ni=str2double(mi{1});nj=str2double(mj{1});
assert(ni>=3&&nj>=3&&ni*nj<=1e7,'Suspicious VEC grid dimensions.');
assert(~isempty(regexpi(zone,'F\s*=\s*POINT','once')), ...
    'Only point-wise VEC data supported.');
C=textscan(fid,'%f%f%f%f%f','CollectOutput',true,...
   'Delimiter',{' ',char(9),','},'MultipleDelimsAsOne',true);
raw=C{1};
assert(size(raw,2)==5&&size(raw,1)==ni*nj,...
    'VEC numeric row count mismatch or unsupported columns.');
X=reshape(raw(:,1),ni,nj)';Y=reshape(raw(:,2),ni,nj)';
U=reshape(raw(:,3),ni,nj)';V=reshape(raw(:,4),ni,nj)';
CHC=reshape(raw(:,5),ni,nj)';
if median(diff(X(1,:)))<0
    X=fliplr(X);Y=fliplr(Y);U=fliplr(U);V=fliplr(V);CHC=fliplr(CHC);
end
if median(diff(Y(:,1)))<0
    X=flipud(X);Y=flipud(Y);U=flipud(U);V=flipud(V);CHC=flipud(CHC);
end
[X,Y,U,V,conv]=pof_normalize_vec_units(X,Y,U,V,meta,umPerPixel,pulseUs);
m=round(cropFraction*min(size(U)));
rows=(m+1):(size(U,1)-m);cols=(m+1):(size(U,2)-m);
X=X(rows,cols);Y=Y(rows,cols);
U=U(rows,cols);V=V(rows,cols);CHC=CHC(rows,cols);
assert(all(isfinite(X(:)))&&all(isfinite(Y(:))), ...
    'Missing grid positions.');
assert(all(diff(X(1,:))>0)&&all(diff(Y(:,1))>0), ...
    'Coordinate axes are not strictly increasing.');
dxi=diff(X(1,:));dyi=diff(Y(:,1));
assert(max(abs(dxi-median(dxi)))<=1e-4*median(dxi) && ...
       max(abs(dyi-median(dyi)))<=1e-4*median(dyi), ...
    'Nonuniform vector grid: expected calibrated rectilinear PIV lattice.');
F=struct('X',X,'Y',Y,'U',U,'V',V,'CHC',CHC,...
   'dx',median(diff(X(1,:))),'dy',median(diff(Y(:,1))),...
   'metadata',meta,'conversion',conv,...
   'velocityFieldsArePhysical_mps',true);
end
