# MassGap-EnergyLevels

An autonomous research program asking whether the Yang-Mills mass gap can be explained with the mathematics and concepts behind electron energy levels and energy gaps: bound-state quantization, band gaps, Mott and BCS gaps, and the spectral-gap methods used for Haldane and AKLT spin chains.

Run by a Claude Code overseer agent that dispatches teams of research subagents in a loop. See `OVERSEER_PROMPT.md` for the program and its stop conditions, `STATE.md` for current status, and `CLAIMS.md` for claims and their status.

This project does not claim to solve the Clay Millennium problem. Possible outcomes are a rigorous partial result, a quantitative model of the glueball spectrum, or a precise account of where each analogy fails. Nothing here is a conclusion unless `RESULT.md` exists.

## Layout
- `literature/`: paper notes with citations
- `data/`: lattice spectra and string tensions with sources
- `src/`: shared code
- `lean/`: Lean 4 formalizations
- `work/`: one folder per task, including dead ends
