# Public-release gates

1. Historical PIV scripts with double m/s conversion must remain outside the supported entry point. **DONE**: stored separately.
2. Main module has no hidden raw-data scaling. **DONE** via unit normalization boundary and method functions accepting m/s.
3. README and methods description must disclose pre-smoothing ambiguity and record that exact numerical output equivalence is untested. **DONE**.
4. Final frozen figure CSV/producer provenance and the historical 81775 vs 81773 input for SI S6/S7 must be resolved **BEFORE claiming complete reproduction**. OPEN.
5. True raw-data-to-paper end-to-end MATLAB regression and MATLAB syntax execution of new kernels. OPEN (environment lacks MATLAB/Octave; static checks only).
6. Authors/coauthors' code redistribution approval and an explicit license. OPEN.

The private GitHub ZIP is uploadable as a software **reference-methods** repository with these truthful claims. Do not describe it as verified end-to-end paper replication or fully scientifically certified.
