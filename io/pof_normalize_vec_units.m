function [x_m,y_m,u_mps,v_mps,info] = pof_normalize_vec_units(x,y,u,v,meta,umPerPixel,pulseUs)
% Strict unit boundary for a single VEC/DAT field, WITHOUT PIV reprocessing.
% Metadata MUST come from paired VEC header, not generic DAT numeric values.
% Pass calibration and pulse interval explicitly; do not infer from the day name.
% For native m/s velocity there is no multiplication; pixel means displacement.
assert(isstruct(meta) && isfield(meta,'positionUnit') && isfield(meta,'velocityUnit'), ...
    'Explicit VEC header metadata is required.');
assert(isnumeric(x)&&isnumeric(y)&&isnumeric(u)&&isnumeric(v) && ...
    isequal(size(x),size(y),size(u),size(v)), ...
    'Coordinates and velocity arrays must be numeric and equal sized.');
assert(isnumeric(umPerPixel)&&isscalar(umPerPixel)&&isfinite(umPerPixel)&&umPerPixel>0, ...
    'A verified positive calibration (micrometers/pixel) is required.');
assert(isnumeric(pulseUs)&&isscalar(pulseUs)&&isfinite(pulseUs)&&pulseUs>0, ...
    'A verified positive pulse interval (microseconds) is required.');
if isfield(meta,'micrometersPerPixel') && isfinite(meta.micrometersPerPixel)
    assert(abs(meta.micrometersPerPixel-umPerPixel)<=1e-3*umPerPixel, ...
        'Supplied spatial calibration conflicts with VEC header.');
end
if isfield(meta,'pulseIntervalUs') && isfinite(meta.pulseIntervalUs)
    assert(abs(meta.pulseIntervalUs-pulseUs)<=1e-3*pulseUs, ...
        'Supplied pulse interval conflicts with VEC header.');
end
switch lower(char(meta.positionUnit))
    case 'mm'
        positionFactor=1e-3;
    case 'pixel'
        positionFactor=umPerPixel*1e-6;
    otherwise
        error('Unsupported or unverified coordinate unit.');
end
switch lower(char(meta.velocityUnit))
    case 'm/s'
        velocityFactor=1.0;
    case 'pixel'
        velocityFactor=umPerPixel/pulseUs; % (micrometers/pixel)/(microseconds) -> (m/s)/pixel
    otherwise
        error('Unsupported or unverified velocity unit.');
end
x_m=x*positionFactor;
y_m=y*positionFactor;
u_mps=u*velocityFactor;
v_mps=v*velocityFactor;
info.positionFactor=positionFactor;
info.velocityFactor=velocityFactor;
info.velocitySourceUnit=char(meta.velocityUnit);
info.positionSourceUnit=char(meta.positionUnit);
end
