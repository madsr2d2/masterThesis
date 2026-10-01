# orcamon in its ecosystem

Where orcamon sits among the tools around ORCA and computational chemistry,
what it does that nothing else was found to do, what it should stop trying to
win, and what that implies for development. Written 2026-10-01 against
orcamon 0.1.0, from a survey of PyPI, GitHub and project sites made the same
day.

This document is a **development aid**, not marketing. Section 8 says exactly
how far each claim was checked, and section 9 says how to repeat the survey.
Star counts and dates are as of 2026-10-01 and will drift.

## 1. Summary

- **orcamon's defensible core is live-log triage with a stable machine
  interface**: it says whether a job that ORCA is still writing is running,
  stalled, stopped, quiet or failed, and answers in bounded `--json` with exit
  codes, plus `wait`. No tool found does this. It is also the part that cannot
  be reduced to "a script around a parser", because the input is an unfinished
  file.
- **The second asset is one report model behind two tiers**, with a skill
  generated from the same argparse parsers and version-checked against the
  installed copy. A person in the TUI and an agent on the CLI see the same
  numbers. No other tool found has this, or ships a tested skill with the
  tool.
- **The molecule and normal-mode viewer is not a differentiator any more.**
  MolTUI (147 stars, created April 2026) does terminal normal modes,
  trajectories and orbitals, and is a dedicated project. orcamon's viewer is
  worth keeping as a convenience on the job's own geometry, not as the pitch.
- **Finished-job parsing is crowded and now has an official library.** OPI,
  maintained by FACCTs (the ORCA co-developers), parses ORCA 6.1 output; cclib
  and several smaller parsers exist. orcamon should not compete there, and
  already says its numbers are monitoring values.
- **Agent tooling for ORCA exists but is shallow on monitoring.** The ORCA MCP
  servers found generate inputs or run a calculation and return a file;
  `dft-skills-market` offers prompt-level parsing and diagnosis guidance;
  DELFIN has a real correctness critic but inside a larger platform. None
  gives an agent a bounded, tested view of a running job.
- **The window is real but not permanent.** A cclib or OPI wrapper with JSON
  output is a small project for someone else. What lasts is the accumulated
  knowledge of ORCA's awkward output (section 4) and the contract's stability.

## 2. What orcamon is, in one table

Taken from `README.md` and `orcamon --help` (0.1.0).

| Aspect | orcamon |
|---|---|
| Input | The live or finished `.out` (plus the `.inp`), read from a directory tree |
| Output | A `JobReport`: `show`, `ls`, `conv`, `energies`, `geom`, `measure`, `freqs`, `input`, `errors`, `tail`, `wait`, `snapshot`; every command has `--json` (`schema: 1`) and meaningful exit codes (0/1/2/4/130) |
| Status model | `not run`, `queued`, `failed`, `finished`, `stalled?`, `running`, `quiet`, `stopped`, plus attention flags such as `opt_not_converged`, `scan_incomplete`, `irc_not_converged`, `ts_hessian`, `ts_imaginary`, `minimum_imaginary`, `qm2_errors` |
| Liveness | `process` (the process table), `slurm` (`squeue --me`, matched by working directory), `mtime` (copied trees); it never claims a job is dead from mtime alone |
| Job kinds | Opt, OptTS, frequencies, relaxed scans, IRC, NEB, QM/QM2 multilayer (energy of the QM/QM2 total, not the QM1 region) |
| Human tier | Textual TUI over ssh: every job live, convergence charts, geometry (pixel via herdr or Kitty graphics, braille text elsewhere), normal-mode animation, notifications (bell, OSC 9/777, `--on-event` hook) |
| Agent tier | Standard-library-only core; shipped skill (`orcamon skill install`, `--check` detects staleness) |
| Install | `uv tool install` / `pipx`; extras `tui`, `images`, `procs`; Python 3.10+ |
| Not in scope | Property JSON, inputs, submitting or killing jobs, spectra, orbitals, a web UI, an MCP server (see `PLAN_ORCAMON.md` section 10) |

## 3. The ecosystem, by what each tool is for

### 3.1 Live monitoring of ORCA output

| Tool | What it is | Evidence | Relation to orcamon |
|---|---|---|---|
| [MonitorOLive](https://github.com/AzuleneG/MonitorOLive) | Python module; plots SCF, Opt/OptTS and NEB-TS convergence in gnuplot windows, writes CSV, splits an `.out` by cycle or frequency block | README read; MIT, created 2026-08-04, 0 stars | Nearest direct competitor in intent. One file at a time, needs gnuplot, no status model, no JSON, no multi-job view. |
| [CCWatcher](https://ccwatcher.sourceforge.net/) | Qt GUI and CLI; live parsing and SCF energy plots for ORCA, Gaussian, Turbomole, Molpro, GAMESS, NWChem and others; multi-file comparison | Project site read; v1.3.2, site last updated 2015-04 | Dormant. Covers many programs, none deeply. |
| [SCFMonitor](https://cran.rstudio.com/web/packages/SCFMonitor/index.html) | R package; live SCF and optimization graphs | Search summary only | Gaussian logs only. |

### 3.2 Job queues and schedulers (no ORCA knowledge)

[turm](https://terminaltrove.com/turm/), [slurmtui](https://pypi.org/project/slurmtui/0.1.4),
[spilot](https://pypi.org/project/spilot/), [swatch](https://lib.rs/crates/swatch),
[slurmer](https://rust-digger.code-maven.com/crates/slurmer) and stui show the
Slurm queue and tail stdout. Search summaries only. They tell you a job is
`RUNNING`; they cannot tell you it stopped converging. orcamon's SLURM probe
is a liveness input, not a queue manager, and is complementary. Generic MCPs
for [Slurm](https://claudemarketplaces.com/mcp/io.github.iowarp/slurm-mcp)
and for [cron/timer failures](https://claudemarketplaces.com/mcp/temurkhan13/silentwatch-mcp)
exist in the same vein.

### 3.3 Finished-output parsers and the official library

| Tool | Notes | Evidence |
|---|---|---|
| [faccts/opi](https://github.com/faccts/opi) | ORCA Python Interface, by FACCTs. Creates input and parses output; JCTC paper; v2.0.0; 191 stars; requires ORCA 6.1.1+. Has an output "grepper" with patterns for markers such as the optimization-converged line. Links an MCP docs page. | README read; source layout listed, `grepper` not read |
| [cclib](https://github.com/cclib/cclib) | The standard multi-program parser; underpins aiida-orca and many pipelines | Existence confirmed by code search; not read |
| [ChemParse](https://github.com/ChemParse/ChemParse), [orca_parser](https://github.com/avanteijlingen/orca_parser), [ORCAParser](https://github.com/HenriqueCSJ/ORCAParser), [FAIRmat nomad-parser-orca](https://github.com/FAIRmat-NFDI/nomad-parser-orca) | Smaller or domain-specific parsers | Repository descriptions only |
| [radi0sus/orca_ir, orca_uv, orca_orb, ...](https://github.com/radi0sus/orca_ir) | Single-purpose plotters and extractors | Descriptions only |

All of these assume a finished output and return values. None models "is it
still going", and none is aimed at an agent.

### 3.4 Terminal molecule and mode viewers

| Tool | Notes | Evidence |
|---|---|---|
| [MolTUI](https://github.com/kszenes/moltui) | Unicode terminal viewer for geometries, **trajectories**, **orbitals** and **normal modes**, built for ssh. Reads XYZ, CIF, ZMAT, Cube, Molden, Gaussian `.fchk`, ORCA `.gbw` (needs `orca_2mkl`) and ORCA `.hess`. Geometry sidebar with distances, angles, dihedrals; PNG export. MIT, 147 stars, last push 2026-08-03. | README read |
| [briling/v](https://github.com/briling/v) | X11 viewer with vibration animation, reads outputs through cclib | README read; a GUI, not a terminal tool |
| Avogadro 2, Chemcraft, ChimeraX | The viewers [ORCA's own manual recommends](https://www.faccts.de/docs/orca/6.1/manual/contents/quickstartguide/gui.html); ORCA has no GUI | Manual read via search |

MolTUI works from **files** (`.hess`, Molden, `.gbw`), not from a job's live
log. orcamon takes modes and geometries from the running job's own output.
Whether MolTUI reads a plain ORCA `.out` for geometry or frequencies was not
established: its README lists ORCA under *inputs*.

### 3.5 Agent tooling for ORCA and for computational chemistry

| Tool | What it does | Monitors a job? | Evidence |
|---|---|---|---|
| [ORCA MCP server](https://github.com/PhelanShao/orca-mcp-server) | `generate_input_file`, `validate_input_syntax`, `suggest_keywords` | No | Description read; 19 stars |
| [AI4S-agent-tools `ORCA_tools`](https://github.com/deepmodeling/AI4S-agent-tools/tree/main/servers/ORCA_tools) | MCP server with two tools: run a calculation in a temp directory and return file paths; retrieve content from an output | No | README and tool list read |
| [dft-skills-market](https://github.com/sylcliff/dft-skills-market) | Prompt-level skills for ORCA/VASP/QE/Gaussian/CP2K: input generation, `output-parse`, `error-diagnosis`, Slurm templates | No (guidance, not a tool) | README read; 8 stars, created 2026-05 |
| [DELFIN](https://github.com/ComPlat/DELFIN) | AI-orchestrated platform: SMILES in, properties out; agents configure, submit and monitor calculations; includes a read-only **correctness critic** for finished ORCA runs (saddle points, spin contamination, unconverged SCF, wrong TS imaginary-mode count) | Inside its own workflow | README and `result_critic.py` header read |
| [catgo-LRG](https://github.com/Hello-QM/catgo-LRG) | Tauri/Svelte workbench with an AI assistant, workflow DAG and HPC submission, driving VASP, ORCA, CP2K and others; 205 stars | Inside its own platform | Description only |
| [Hello-QM ORCA skills](https://www.skillsdirectory.com/skills/hello-qm-neb-ts) | Claude Code skills routing to opt, NEB-TS and freq workflows; status comes from the CatGo engine | Via CatGo | Page read |
| [dft-autopilot-mcp](https://github.com/The66user/dft-autopilot-mcp) | MCP for QE/VASP/Gaussian over Slurm; `check_job_status`, `download_job_results` | Queue state only; no ORCA | README read |

The pattern: agents are being given ORCA as a tool to **run**, and platforms
that embed monitoring do it privately. The standalone, installable piece
that answers "what is this running job doing" is the gap.

### 3.6 Workflow engines and automation (neighbours, not competitors)

autodE, pysisyphus, AiiDA (`aiida-orca`), `chemsmart`, `AaronTools`, QMatSuite
and similar drive ORCA and parse results as part of a pipeline. They own the
job lifecycle, so they do not need an external monitor; a user running jobs by
hand, from a PBS-free homelab or a shell, does. This is also who might one day
call orcamon.

## 4. Where orcamon is distinct

Ranked by how confident the claim is.

1. **A status model for unfinished output**, with liveness from the process
   table or `squeue` and an explicit refusal to call a job dead from mtime
   alone. Not found anywhere.
2. **A stable, versioned contract**: `--json` with units in key names, a
   `schema` integer, closed flag codes, defined exit codes, and `wait` with a
   timeout under an agent shell's ceiling. This is what turns a script into an
   interface.
3. **The pitfalls of ORCA's output are handled, and tested**: four energies
   per QM/XTB geometry, the QM/QM2 total rather than the QM1 region, scan
   cycles that restart, IRC and NEB paths (NEB is rewritten every iteration),
   Hessians recomputed mid-optimization labelled rather than flagged,
   imaginary modes counted without a cutoff. Each is the sort of thing a
   generic parser or an agent's `grep` gets wrong. This accumulated knowledge
   is the hardest part for a newcomer to copy.
4. **One `JobReport` for both tiers**, so the TUI summary and `orcamon show`
   cannot disagree.
5. **A skill released with the tool**, generated from the parsers and checked
   for staleness. No equivalent found; `dft-skills-market` is hand-written
   prose with no tests against a tool.
6. **Low-friction install**: a standard-library core, no gnuplot, optional
   extras.

## 5. Where orcamon is not distinct

- **Normal-mode and trajectory viewing.** MolTUI is broader (orbitals, many
  formats, crystals) and maintained as its own project. `PLAN_ORCAMON.md`
  section 10 originally put "orbitals, densities, vibration animation" out of
  scope and handed them to Avogadro or Chemcraft; the vibration work later
  went in anyway. Treat it as a convenience that shows the job's own modes
  without leaving the monitor, and do not extend it toward orbitals or
  crystals, where MolTUI would always win.
- **Parsing a finished job into values.** OPI is official and cclib is the
  standard. orcamon's README is right that its numbers are for monitoring and
  triage, and that quoted numbers should come from ORCA's structured output.
- **Convergence plots.** MonitorOLive and CCWatcher do this; orcamon's chart
  is better integrated (scrubbing, multi-job) but not a new capability.
- **Input generation and running jobs.** Several MCP servers and platforms do
  it. orcamon has chosen not to.

### 5.1 Rendering: MolTUI against orcamon

Read from source (MolTUI `image_renderer.py` and `app.py`; orcamon `tui/raster.py`
and `tui/geometry_text.py`). **Neither was run for this comparison**, so
statements about how the output looks are inferred from the code.

**MolTUI draws one shaded 3D image and converts it to text.** A numpy
z-buffer rasteriser renders spheres with Blinn-Phong shading (ambient,
diffuse, specular) and shaded cylinders for bonds. The frame is sized
`cols*2` by `rows*4` pixels, one pixel per braille dot. Each terminal cell then
takes:

- its **shape** from the hit mask: the braille character whose dots are the
  ones that touched a surface, so the silhouette is drawn at 2x4 dots per
  cell and occlusion comes from the z-buffer;
- its **colour** from the mean of the covered dots, as 24-bit RGB.

The result is lit, depth-occluded spheres in any truecolor terminal, with no
graphics protocol, so it survives tmux, mosh and any ssh path. The limit is
that shading is per-cell colour, since a cell averages its dots into one.

**orcamon has two renderers.** The pixel path (Kitty or herdr graphics,
`raster.py`) is richer than MolTUI's: fog, depth-based silhouette halos,
separate QM and environment layers, a see-through mode, and labels and
distances drawn into the frame. The text path (`geometry_text.py`) is not a
rasteriser: bonds are braille lines, atoms are element letters, with no
shading or occlusion. It is pure Python with no numpy.

| | MolTUI | orcamon |
|---|---|---|
| Look in tmux, mosh, plain terminals | Shaded spheres via braille plus truecolor | Braille bond lines plus letters, flat |
| Needs a graphics protocol for the good picture | No | Yes (Kitty or herdr) |
| Occlusion in the text path | Yes (z-buffer hit mask) | No |
| Specular highlights, adjustable lighting | Yes | Fixed lighting, pixel path only |
| Fog, silhouette halos, see-through | No | Yes, pixel path only |
| QM and environment layers | No | Yes, both paths |
| Dependencies for the viewer | numpy | Text path none; pixel path numpy and PIL |
| Source of the geometry | A file you open (XYZ, Molden, `.hess`, and others) | The running job's own log |
| Orbitals, periodic cells | Yes | No |

**What this means for orcamon.** The text pane is the weak point, and it is
the one that matters, because ssh and tmux are where orcamon is meant to
run. `raster.render` already produces a z-buffered, lit RGB frame. Two
changes would give a shaded text mode:

1. Expose the hit mask. It is `frame.zbuf != -inf` and stays internal to
   `render` today, which returns only the PIL image.
2. Render a frame at `2*cols` by `4*rows` pixels and reduce it to braille
   cells as MolTUI does: dot mask for the shape, mean colour for the
   foreground.

It would need numpy, so it would sit behind the `images` extra with the
current pure-Python path kept as the standard-library fallback. Fog, halos
and the QM/environment layers would carry over. Open questions that only a
prototype answers:

- **Thin features.** At 2x4 dots per cell, the environment layer's thin
  sticks (0.06 radius against 0.10 for QM bonds) may vanish. MolTUI renders
  at `ssaa=1` in the TUI, so a bond there is a dot or two wide.
- **Labels.** The element letters in today's text pane say which atom is
  which. A braille image needs an overlay of label cells, as MolTUI does for
  atom numbers.
- **Colour depth.** 24-bit colour is assumed. A 256-colour fallback and a
  probe for support would be needed.

## 6. Risks

| Risk | Likelihood and effect | Mitigation |
|---|---|---|
| A competitor wraps cclib/OPI in a JSON CLI | Plausible; small effort for them. Would cover finished jobs only. | Keep the live-state and ORCA-quirk lead; keep the contract stable; publish the status model. |
| FACCTs adds monitoring to OPI or ORCA itself | Possible; they are the authority on the format. | Track OPI releases. If they add a status API, consider consuming it instead of competing. |
| Output format drift in new ORCA versions | Certain over time. orcamon parses text markers. | `python -m orcamon.validate ROOT` over real trees on each new ORCA release; fixtures from ORCA's own output in `test_monitor.py`. |
| MolTUI makes the viewer redundant | Already largely true for modes. | Do not invest further; consider interop (6.2 below). |
| The tool is mostly known by one user | Real. 0.1.0, local-only, one repository. | Decide whether to publish (section 7). |
| Name collisions | `orcamon` is **not registered on PyPI** (checked 2026-10-01: 404 for `orcamon`, `orca-mon`, `orcamonitor`). GitHub has unrelated `orcamonitor` and `orcamonsta` repositories. | Reserve the name when publishing. |

## 7. Development implications

Ordered by expected value, not effort. None is committed; each is a
suggestion for the next planning round.

### 7.1 Strengthen the core

1. **Publish the status and flag semantics as a specification**, separate from
   the code: what each status means, which evidence decides it, which flags
   exist and when each fires. This is the artifact a competitor cannot get by
   wrapping a parser, and it makes the JSON contract citable. The README table
   is a start.
2. **Cross-check the finished-job values against an independent parser.** Run
   orcamon and OPI (or cclib) over the same finished outputs in a test and
   compare final energy, convergence and frequencies. It is a check with an
   independent source, which also gives a concrete statement for the README:
   "agrees with OPI on N outputs". Keep the dependency out of the core and out
   of the shipped package; test-only.
3. **Version coverage.** OPI requires ORCA 6.1.1+. State which ORCA versions
   orcamon is tested against, and keep sample outputs from 5.x and 6.x in the
   fixtures. Being the tool that still reads 5.x logs on old clusters is a
   practical advantage worth checking, not assuming.
4. **More failure modes.** The set in `errors` and the flags is where
   expertise accumulates. Candidates: SCF non-convergence trajectories
   (oscillation versus slow creep), DIIS stalls, disk or memory errors, MPI
   aborts, jobs killed by a walltime. Each needs a real output to test
   against.

### 7.2 Agent distribution

5. **Consider a thin MCP wrapper**, optional and separate from the core. The
   plan declined an MCP server because `ssh host orcamon ls --json` is enough,
   and for Claude Code with a skill that remains true. But the surveyed ORCA
   MCP servers show where non-Claude-Code agents look. A wrapper that shells
   out to the existing commands costs little and does not touch the contract.
6. **Keep the skill as the primary agent route** and make `orcamon skill
   install` work for other agents' formats only if a user asks.

### 7.3 Learn from MolTUI, and interoperate instead of competing

7. **Prototype a shaded braille text mode** from the existing rasteriser
   (section 5.1), as a third `--graphics` option, and compare it with the
   current text pane on the same geometry before deciding. This is the one
   place MolTUI is plainly ahead of orcamon on something orcamon needs.
8. **Hand off to MolTUI for rich viewing.** orcamon already writes geometries
   and could point to a `.hess` or Molden file: a `--open-in moltui` style
   helper, or documentation of the one-line command, costs little and turns a
   competitor into a downstream tool.
9. **Emit what workflow engines can consume.** `show --json` already does.
   Documenting it for autodE/pysisyphus/AiiDA users as "a liveness probe for
   jobs you launched" is an adoption route.

### 7.4 What not to do

- Do not add input generation, job submission, orbitals, spectra or a web UI.
  Others do them, and each dilutes the claim that orcamon answers one question
  well.
- Do not import thesis-specific rules (the project memory already records
  this): the tool's value is that it is general.
- Do not market it as a parser or a viewer.

### 7.5 Publishing

If it is published, the minimum is: a repository of its own (the plan notes a
`git subtree split`), a PyPI release under a reserved name, a related-work
section in the README taken from sections 3 and 5 here, and a changelog tied
to `schema`. `pyproject.toml` declares MIT but the package directory has no
LICENSE file; add one before any release.

## 8. How far this was checked

**Read directly (this session):** orcamon's README, `--help`, `PLAN_ORCAMON.md`;
MolTUI's `image_renderer.py` and `app.py` (the rendering pipeline), and
orcamon's `tui/raster.py`, `tui/geometry_text.py` and `tui/geometry_render.py`;
READMEs of MonitorOLive, MolTUI, OPI, dft-skills-market, dft-autopilot-mcp,
DELFIN (plus the head of `result_critic.py`), and `ORCA_tools`' README and
function list; CCWatcher's project page; the Hello-QM NEB-TS skill page; the
ORCA MCP server page.

**Seen only as search results or repository descriptions:** cclib's ORCA
parser (code-search hit), ChemParse, orca_parser, ORCAParser, nomad-parser-orca,
the radi0sus plotters, catgo-LRG, QMatSuite, chemsmart, AiiDA-ORCA,
SCFMonitor and every Slurm TUI. Claims about these are limited to what their
descriptions say.

**Not examined:** source code of any competitor beyond the files named above;
the running output of MolTUI (installed or executed), so no visual comparison
has been made;
OPI's `grepper` and whether it covers live or truncated logs; any tool's
behaviour when run; GitLab, Codeberg and institutional repositories; private
lab scripts; the ORCA forum; commercial software such as Chemcraft;
Q-Chem, Psi4 and Turbomole ecosystems.

Two earlier statements were wrong and are corrected here: MolTUI **does**
support normal modes and trajectories (first-pass search summaries said
otherwise), and GitHub keyword search is a poor test of absence, since
"orca" matches many unrelated projects.

## 9. Repeating the survey

Re-run before a release, a README claim, or a thesis statement.

```bash
# repositories by name and topic
gh search repos "orca-parser" --sort stars --limit 20
gh search repos "orca-output" --sort stars --limit 20
gh search repos --topic=orca-quantum-chemistry --limit 20
gh search repos --topic=computational-chemistry --topic=mcp --limit 20
gh search repos --topic=computational-chemistry --topic=claude-code --limit 20
gh search repos "moltui"

# ORCA output markers in code (surfaces parsers and tools that read logs)
gh api -X GET search/code -f q='ORCA TERMINATED NORMALLY language:Python' -f per_page=30
gh api -X GET search/code -f q='HURRAY THE OPTIMIZATION HAS CONVERGED' -f per_page=30
gh api -X GET search/code -f q='orca_tools' -f per_page=30

# the name
curl -s -o /dev/null -w "%{http_code}\n" https://pypi.org/pypi/orcamon/json
```

Things worth checking by hand each time: whether OPI or MolTUI gained a
monitoring or status feature; whether any tool offers `wait`, exit codes or a
JSON status for a **running** ORCA job; whether `dft-skills-market` or DELFIN
split out a standalone checker.
