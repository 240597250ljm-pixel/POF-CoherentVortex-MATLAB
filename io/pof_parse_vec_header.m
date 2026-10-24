function meta = pof_parse_vec_header(headerOrFile)
% POF_PARSE_VEC_HEADER Read explicit VEC field units and acquisition metadata.
% Accepts a VEC path or a header literal. Never guesses units from magnitudes.
% Unlike v0.3, tolerates metadata spread over several header lines.
if isstring(headerOrFile),headerOrFile=char(headerOrFile);end
assert(ischar(headerOrFile) && ~isempty(headerOrFile), ...
    'Expected nonempty VEC path or literal header.');
if isfile(headerOrFile)
    fid=fopen(headerOrFile,'rt');
    assert(fid>0,'Could not open VEC file.');
    c=onCleanup(@() fclose(fid));
    fragments=cell(0,1);
    for j=1:24
        line=fgetl(fid);
        if ~ischar(line),break;end
        fragments{end+1}=line; %#ok<AGROW>
        if ~isempty(regexpi(line,'(^|\s)ZONE(\s|$)','once')),break;end
    end
    header=strjoin(fragments,' ');
else
    assert(~isempty(regexpi(headerOrFile,'VARIABLES\s*=','once')), ...
        'VEC file not found or literal header lacks VARIABLES.');
    header=headerOrFile;
end
assert(~isempty(regexpi(header,'VARIABLES\s*=','once')), ...
    'VEC header lacks VARIABLES line; no implicit unit inference.');
x=regexp(header,'"X\s+(mm|pixel)"','tokens','once','ignorecase');
y=regexp(header,'"Y\s+(mm|pixel)"','tokens','once','ignorecase');
u=regexp(header,'"U\s+(m/s|pixel)"','tokens','once','ignorecase');
v=regexp(header,'"V\s+(m/s|pixel)"','tokens','once','ignorecase');
assert(~isempty(x)&&~isempty(y)&&~isempty(u)&&~isempty(v), ...
    'Missing explicit X/Y/U/V units. DAT headers without units are invalid.');
assert(strcmpi(x{1},y{1}) && strcmpi(u{1},v{1}), ...
    'X/Y or U/V unit mismatch in VEC header.');
meta=struct();
meta.positionUnit=lower(x{1});
meta.velocityUnit=lower(u{1});
meta.micrometersPerPixel=NaN;
meta.pulseIntervalUs=NaN;
a=regexp(header,'MicrometersPerPixelX\s*=\s*"([0-9.eE+-]+)"', ...
    'tokens','once','ignorecase');
b=regexp(header,'MicrosecondsPerDeltaT\s*=\s*"([0-9.eE+-]+)"', ...
    'tokens','once','ignorecase');
if ~isempty(a),meta.micrometersPerPixel=str2double(a{1});end
if ~isempty(b),meta.pulseIntervalUs=str2double(b{1});end
if ~isnan(meta.micrometersPerPixel)
    assert(isfinite(meta.micrometersPerPixel)&&meta.micrometersPerPixel>0, ...
        'Invalid spatial calibration metadata.');
end
if ~isnan(meta.pulseIntervalUs)
    assert(isfinite(meta.pulseIntervalUs)&&meta.pulseIntervalUs>0, ...
        'Invalid pulse interval metadata.');
end
end
