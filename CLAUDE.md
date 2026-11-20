# Repo rules (overseer and all subagents)

Research program: can the Yang-Mills mass gap be explained with the mathematics and concepts behind electron energy levels and gaps (bound states, band gaps, Mott/BCS gaps, Haldane/AKLT spectral-gap methods)? The full program, stop conditions, and team structure are in `OVERSEER_PROMPT.md`. Read it before doing anything.

- Work only inside this repo. Do not read, list, or reference any other local folder or repository.
- `STATE.md` is the current state of the program. Read it first. Only the overseer edits it, `CLAIMS.md`, and `RESULT.md`.
- Each task lives in `work/<YYYYMMDD>-<NNN>-<slug>/` with a `README.md` (question, method, result, verdict, open issues) and the scripts needed to rerun it.
- Every proof step is marked proved, cited, or conjectured. Every claim states its setting: lattice or continuum, finite or infinite volume, dimension, gauge group.
- A numerical gap at finite volume or lattice spacing is an observation, not a proof.
- Lattice and measured numbers come from `data/` with citations in `data/SOURCES.md`. Don't hardcode them in scripts.
- Python lives in `.venv/`, pinned in `requirements.txt`. Lean 4 project lives in `lean/`.
- CPU: this program owns logical CPUs 9-11 only (3 CPUs); 0-3 are reserved for the user. Run every compute command as `python tools/pin.py <command>` and size threads/workers to 3 total. Intel iGPU, no CUDA. Long jobs run in the background and log to `work/<id>/run.log`.
- Repo is public. No personal info, credentials, or absolute local paths in committed files. Don't change the git identity, don't force-push.
- Never claim the Millennium problem is solved. Use "candidate" or "partial" unless `RESULT.md` exists.
