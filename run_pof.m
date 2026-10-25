function out = run_pof(command)
% Single safe entry point. No legacy scripts, raw PIV recomputations or
% manuscript figure regenerations are launched automatically.
if nargin<1 || isempty(command),command='status';end
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'io'),fullfile(root,'config'),fullfile(root,'tests'));
command=lower(char(command));
switch command
 case 'status'
    out=struct('version','1.0-rc4','implementation','reference-methods', ...
       'productionFiguresCertified',false,'scienceRegressionTested',false);
    fprintf('POF reference methods v%s | no historical computation auto-run.\n',out.version);
    fprintf('Synthetic tests: run_pof(''test''); code is NOT a certified figure reproduction.\n');
 case 'test'
    out=run_method_tests();
 case 'preflight'
    out=pof_preflight(pof_config_local());
 case 'methods'
    d=dir(fullfile(root,'+pof','*.m'));
    out={d.name}';disp(out);
 case 'example'
    out=run_synthetic_example();
 otherwise
    error('Allowed commands: status, test, preflight, methods, example.');
end
end
