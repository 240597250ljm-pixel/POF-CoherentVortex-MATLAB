function out = pof_preflight(cfg)
% POF_PREFLIGHT Read-only local provenance and cohort check.
% Checks existence / counts of VECs, sampled original VEC headers,
% an optionally supplied frozen final 2935-frame CSV and 81773-row catalog.
% DOES NOT execute any legacy algorithm or verify all PIV vectors.
checks=cell(0,4);
if ~isstruct(cfg),error('Configuration must be struct.');end
if isfield(cfg,'pathsConfirmedByUser') && cfg.pathsConfirmedByUser
    record('CONFIG','PATHS','PASS','Paths explicitly confirmed in local config');
else
    record('CONFIG','PATHS','NEEDS_REVIEW','Set pathsConfirmedByUser=true after checking folders');
end
if isfield(cfg,'allowProductionRecompute') && cfg.allowProductionRecompute
    record('CONFIG','PRODUCTION','FAIL','Production execution is not supported by v0.4');
else
    record('CONFIG','PRODUCTION','PASS','Legacy science execution remains disabled');
end
expectedNames=cfg.primaryFolders;
scanDay('0604',cfg.primaryRoot,expectedNames,cfg.primaryExpectedUnit, ...
    cfg.primaryPixelCalibrationUm,cfg.primaryPulseIntervalUs);
scanDay('0605',cfg.repeatRoot,cfg.repeatFolders,cfg.repeatExpectedUnit, ...
    cfg.repeatPixelCalibrationUm,cfg.repeatPulseIntervalUs);
if cfg.repeatCalibrationConfirmed && isfinite(cfg.repeatPixelCalibrationUm) ...
        && cfg.repeatPixelCalibrationUm>0 && isfinite(cfg.repeatPulseIntervalUs) ...
        && cfg.repeatPulseIntervalUs>0
    record('0605','CALIBRATION','PASS','Repeat-day calibration confirmed by user');
else
    record('0605','CALIBRATION','NEEDS_REVIEW', ...
        'Repeat-day conversion cannot be certified before original metadata check');
end
if isfield(cfg,'frozenPrimaryFrameCSV') && ~isempty(cfg.frozenPrimaryFrameCSV)
    if isfile(cfg.frozenPrimaryFrameCSV)
        try
            frameT=readtable(cfg.frozenPrimaryFrameCSV);
            validated=pof_validate_primary_cohort(frameT);
            record('FINAL','2935_81773_FRAME','PASS', ...
                sprintf('Frozen frame reference: %d frames / %d vortices', ...
                validated.frames,validated.vortices));
            if isfield(cfg,'frozenPrimaryCatalogCSV') && ...
                    ~isempty(cfg.frozenPrimaryCatalogCSV)
                checkCatalog(frameT,cfg.frozenPrimaryCatalogCSV);
            else
                record('FINAL','81773_CATALOG','NOT_CHECKED', ...
                    'Optional finalized 81773-row vortex catalog not supplied');
            end
        catch ME
            record('FINAL','2935_81773_FRAME','FAIL',ME.message);
        end
    else
        record('FINAL','2935_81773_FRAME','FAIL','Configured frozen CSV not found');
    end
else
    record('FINAL','2935_81773_FRAME','NOT_CHECKED', ...
        'Point frozenPrimaryFrameCSV to final 2935-frame table');
    record('FINAL','81773_CATALOG','NOT_CHECKED', ...
        'No catalog comparison attempted');
end
T=cell2table(checks,'VariableNames',{'Scope','Check','Status','Details'});
T=sortrows(T,{'Scope','Check'});
out=struct();out.checks=T;
out.numPass=nnz(strcmp(T.Status,'PASS'));
out.numFail=nnz(strcmp(T.Status,'FAIL'));
out.numNeedsReview=nnz(strcmp(T.Status,'NEEDS_REVIEW'));
out.numNotChecked=nnz(strcmp(T.Status,'NOT_CHECKED'));
out.scientificResultsReproduced=false;
out.productionEnabled=false;
fprintf('\nPOF v0.4 read-only preflight: PASS=%d FAIL=%d NEEDS_REVIEW=%d NOT_CHECKED=%d\n', ...
    out.numPass,out.numFail,out.numNeedsReview,out.numNotChecked);
disp(T);
if isfield(cfg,'reportDir') && ~isempty(cfg.reportDir)
    if ~isfolder(cfg.reportDir),mkdir(cfg.reportDir);end
    writetable(T,fullfile(cfg.reportDir,'POF_preflight_checks.csv'));
end

    function record(scope,key,state,detail)
        checks(end+1,:)={scope,key,state,detail}; %#ok<AGROW>
    end

    function scanDay(day,root,folders,expectVelUnit,calib,pulse)
        if ~isfolder(root)
            record(day,'ROOT','FAIL',['Root not found: ',root]);return;
        end
        record(day,'ROOT','PASS',['Found root ',root]);
        for i=1:numel(folders)
            condition=folders{i};folder=fullfile(root,condition);
            label=sprintf('%s_%s',day,condition);
            if ~isfolder(folder)
                record(label,'FOLDER','FAIL',['Missing ',folder]);continue;
            end
            d=dir(fullfile(folder,'*.vec'));
            if numel(d)==cfg.expectedFilesPerCondition
                record(label,'FILE_COUNT','PASS', ...
                    sprintf('%d VEC files (acquisition only; not QC)',numel(d)));
            else
                record(label,'FILE_COUNT','NEEDS_REVIEW', ...
                    sprintf('%d VEC files, expected %d',numel(d), ...
                    cfg.expectedFilesPerCondition));
            end
            if isempty(d),continue;end
            pick=unique([1,ceil(numel(d)/2),numel(d)]);
            for j=1:numel(pick)
                vname=d(pick(j)).name;
                try
                    m=pof_parse_vec_header(fullfile(folder,vname));
                    if ~strcmpi(m.velocityUnit,expectVelUnit)
                        record(label,'HEADER_UNIT','FAIL', ...
                            sprintf('%s: %s != expected %s',vname,m.velocityUnit,expectVelUnit));
                        continue;
                    end
                    if isfinite(calib) && isfinite(m.micrometersPerPixel)
                        if abs(m.micrometersPerPixel-calib)>1e-3*calib
                            record(label,'HEADER_CALIB','FAIL',sprintf('%s: calibration mismatch',vname));
                            continue;
                        end
                    end
                    if isfinite(pulse) && isfinite(m.pulseIntervalUs)
                        if abs(m.pulseIntervalUs-pulse)>1e-3*pulse
                            record(label,'HEADER_PULSE','FAIL',sprintf('%s: pulse mismatch',vname));
                            continue;
                        end
                    end
                    record(label,'SAMPLE_HEADER','PASS', ...
                        sprintf('%s: X/%s U/%s',vname,m.positionUnit,m.velocityUnit));
                catch ME
                    record(label,'SAMPLE_HEADER','FAIL',[vname,': ',ME.message]);
                end
            end
        end
    end

    function checkCatalog(frameT,catalogFile)
        if ~isfile(catalogFile)
            record('FINAL','81773_CATALOG','FAIL','Configured vortex catalog not found');
            return;
        end
        try
            C=readtable(catalogFile);
            assert(height(C)==81773,'Final catalog must contain exactly 81773 detections.');
            assert(all(ismember({'Ra','Frame'},C.Properties.VariableNames)), ...
                'Catalog must contain Ra and Frame columns.');
            [g,r,f]=findgroups(C.Ra,C.Frame);
            n=splitapply(@numel,C.Frame,g);
            U=table(r,f,n,'VariableNames',{'Ra','Frame','NVortices'});
            A=frameT(:,{'Ra','Frame','NVortices'});
            A.Properties.VariableNames(strcmp(A.Properties.VariableNames,'NVortices'))={'CountFromFrames'};
            U.Properties.VariableNames(strcmp(U.Properties.VariableNames,'NVortices'))={'CountFromCatalog'};
            L=innerjoin(A,U,'Keys',{'Ra','Frame'});
            assert(height(L)==2935,'Catalog missing matched frame keys.');
            assert(all(L.CountFromFrames==L.CountFromCatalog), ...
                'Frame and catalog per-frame vortex counts differ.');
            record('FINAL','81773_CATALOG','PASS', ...
                '81773 catalog rows, matched counts for all 2935 Ra+Frame keys');
        catch ME
            record('FINAL','81773_CATALOG','FAIL',ME.message);
        end
    end
end
