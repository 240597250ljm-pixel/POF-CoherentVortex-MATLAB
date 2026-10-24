function result = pof_validate_primary_cohort(T)
% Validates the FROZEN final primary cohort, not a historical 2937-frame CSV.
% Needs columns Ra,Frame,NVortices. No PIV vector recalculation.
assert(istable(T),'Expected MATLAB table.');
assert(all(ismember({'Ra','Frame','NVortices'},T.Properties.VariableNames)), ...
    'Required columns are Ra, Frame, NVortices.');
R=[4.63e8;9.26e8;1.39e9;1.85e9;2.32e9;2.78e9];
NF=[498;499;499;492;470;477];
NV=[14840;14434;14182;13480;12676;12161];
assert(height(T)==2935,'Frozen primary cohort must have exactly 2935 frames.');
assert(all(isfinite(T.Ra))&&all(isfinite(T.Frame))&&all(isfinite(T.NVortices)), ...
    'Cohort contains NaN/Inf.');
assert(all(T.Frame==round(T.Frame))&&all(T.NVortices==round(T.NVortices)), ...
    'Non-integer Frame/NVortices.');
assert(height(unique(T(:,{'Ra','Frame'})))==height(T), ...
    'Duplicate Ra+Frame keys.');
for k=1:6
    idx=abs(T.Ra-R(k))<1e-5*R(k);
    assert(nnz(idx)==NF(k), 'Frame count mismatch for condition %d.',k);
    assert(sum(T.NVortices(idx))==NV(k), ...
        'Vortex count mismatch for condition %d.',k);
end
assert(sum(T.NVortices)==81773,'Final catalog requires 81773 detections.');
result.frames=height(T);result.vortices=sum(T.NVortices);
result.framesPerCondition=NF;result.vorticesPerCondition=NV;
end
