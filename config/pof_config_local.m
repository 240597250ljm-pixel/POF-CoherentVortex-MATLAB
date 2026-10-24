function cfg = pof_config_local()
% Portable local configuration for the POF source release.
% Put original PIV data in a separate, untracked directory and supply
% POF_DATA_ROOT (containing 0604/ and 0605/) or edit this local function.
% The code package does NOT contain experimental PIV datasets.
dataRoot = getenv('POF_DATA_ROOT');
if isempty(dataRoot)
    dataRoot = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'external_data');
end
cfg.primaryRoot = fullfile(dataRoot,'0604');
cfg.repeatRoot = fullfile(dataRoot,'0605');
cfg.primaryFolders = {'05','10','15','20','25','30'};
cfg.repeatFolders  = {'5','10','15','20','25','30'};
cfg.expectedFilesPerCondition = 500;
cfg.primaryExpectedUnit = 'm/s';
cfg.repeatExpectedUnit  = 'pixel';
cfg.primaryPixelCalibrationUm = 114.910004;
cfg.primaryPulseIntervalUs = 50;
% 0605 same calibration confirmed by the experimenter; NOT present in VEC header.
cfg.repeatPixelCalibrationUm = 114.910004;
cfg.repeatPulseIntervalUs = 50;
cfg.repeatCalibrationConfirmed = true;
cfg.pathsConfirmedByUser = false;
cfg.frozenPrimaryFrameCSV = getenv('POF_FRAME_CSV');
cfg.frozenPrimaryCatalogCSV = getenv('POF_CATALOG_CSV');
cfg.reportDir = '';
cfg.allowProductionRecompute = false;
end
