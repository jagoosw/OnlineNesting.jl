# OnlineNesting.jl — Agent Rules

Two-way, multi-level nesting of ocean models, built on `NumericalEarth.NestedModels` and
Oceananigans' open boundary conditions.

## Start of every session

1. Read `PLAN.md`: the design (section 3), the milestones (section 4), how sessions and agents
   are used (section 6) and the status log at the end of section 6.
2. Read `git log`.
3. Work on the block of milestones the user names. End by appending to the status log what was
   built, the check results with their numbers, and any departure from the plan.

## Rules

This package follows NumericalEarth's agent rules. The rules in `.claude/rules/` are copies of
NumericalEarth's and load automatically; `restraint-rules.md` and its self-review checklist are
mandatory for every edit. The general rules are imported here:

@~/Documents/NumericalEarth-OMIP/AGENTS.md

Where those files disagree, `AGENTS.md` wins: `$(TYPEDSIGNATURES)` in docstrings, and literal
zeros (`max(0, a)`) in kernels. Oceananigans' own rules (`~/Documents/Oceananigans.jl/AGENTS.md`)
apply to anything proposed upstream there.

## Running

- Tests: `julia --project=. -e 'using Pkg; Pkg.test()'`; one file with
  `Pkg.test(test_args=["test_nested_grid"])`. Loading Oceananigans takes minutes the first time.
- Commit locally at the end of each milestone (the user authorized this on 2026-11-10). Never push without being asked.
