function p = pof_paper_parameters()
% Fixed definitions from the 2026-10-24 POF manuscript and SI.
% Values are documentary; no unverified figure regression is implied.
p=struct();
p.nominalRa=[4.63e8;9.26e8;1.39e9;1.85e9;2.32e9;2.78e9];
p.conditionK=[5;10;15;20;25;30];
p.frozenAcceptedFrames=[498;499;499;492;470;477];
p.frozenVortexCount=81773;
p.vectorSpacing_m=32*114.910004e-6;
p.gridFull=[63 63];p.gridCropped=[57 57];
p.liutexSigmaCells=0.8;p.liutexPositivePercentile=80;
p.vortexMinCells=4;p.vortexBoundaryCells=2;
p.filterFWHM_m=[11.03;15.01;20.43;27.80;37.83;51.48]*1e-3;
p.mainFilterFWHM_m=20.43e-3;
p.nearCutoff=1.5;p.backgroundCutoff=2.0;
p.bootstrapBlock=10;p.bootstrapReplicates=3000;
p.voronoiNullReplicates=[12 24 24 12];
p.voronoiNullExclusionFactors=[0 0 0.5 0.75];
p.structureFunctionBins=14;
p.velocityPreSmoothSigmaCells=0.8; % Manuscript II E/SI 1.1 velocity pre-smoothing
p.velocityPreSmoothKernelTruncationSigma=4; % separate from 3-sigma ell filter
p.isFullThreeDimensionalFlux=false;
end
