# Public-release gates

1. Historical PIV scripts with double m/s conversion must remain outside the supported entry point. **DONE**: stored separately.
2. Main module has no hidden raw-data scaling. **DONE** via unit normalization boundary and method functions accepting m/s.
3. Revised manuscript, SI, code default and documentation must disclose the distinct 0.8-grid U/V pre-smoothing and ell-FWHM coarse graining; exact numerical output equivalence remains untested. **METHOD-LEVEL ALIGNMENT DONE**.
4. Final frozen figure CSV/producer provenance and the historical 81775 vs 81773 input for SI S6/S7 must be resolved **BEFORE claiming complete reproduction**. OPEN.
5. True raw-data-to-paper end-to-end MATLAB regression and MATLAB syntax execution of new kernels. OPEN (environment lacks MATLAB/Octave; static checks only).
6. Authors/coauthors' code redistribution approval and an explicit license. OPEN.

The repository is already public as a **reference-methods** source release. No license is granted by its public visibility alone. Do not describe it as verified end-to-end paper replication or fully scientifically certified.
