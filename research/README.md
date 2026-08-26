# Multi-Agent Autonomous Discovery Lab — Shared Research Repository

This repository is the persistent memory of a coordinated research effort across five
independent sessions (S1–S5), coordinated by the lead coordinator. Goal: produce at least one
result that is **NEW + USEFUL/SIGNIFICANT + REPRODUCIBLE + PROVEN**.

## Layout

```text
/research
  /candidates      machine-readable records (candidate.json) for every serious candidate
  /experiments     experiment code + raw outputs, per session: experiments/S{n}/...
  /benchmarks      fixed benchmark instances/generators shared across sessions
  /discoveries     candidates that survived adversarial testing + novelty audit + reproduction
  /failed          archived dead candidates with post-mortems (do NOT re-investigate)
  /proofs          proof artifacts: exhaustive-search completeness args, certificates, derivations
  /literature      novelty audits: what prior work exists, with citations/links
  /reproductions   independent re-derivations of key results by agents that did not author them
  /reports         session reports and the final lab-wide report
  /leaderboard     scored ranking of all candidates (session files, merged by coordinator)
  /coordination    roadmap, mission briefs, protocol, decision log
  /sessions        working scratch space per session S1..S5 (each session owns its dir ONLY)
```

## Hard rules

1. Each session writes **only** inside `research/sessions/S{n}/`, plus exactly these global touchpoints:
   `research/candidates/<ID>.json`, `research/experiments/S{n}/...`, `research/literature/S{n}-novelty.md`,
   `research/reports/S{n}-report.md`. The coordinator owns everything else.
2. No session runs `git` commands; the coordinator commits.
3. Every number claimed in a report must be regenerable by a committed script + command line in the report.
4. Claim labels are strict: `observation < hypothesis < conjecture < experimentally_validated_result < proven_result < formally_verified_result`.
5. Nothing is called a discovery until the novelty audit is done. Failed candidates go to `/failed` with a post-mortem.
6. Experiments must be timeboxed and deterministic where possible (fixed seeds, pinned inputs).

## Scoring

Score = Novelty × Usefulness × Provability × Reproducibility × Improvement, each factor 0–5,
anchored in `coordination/PROTOCOL.md`. Coordinator computes scores; sessions provide self-assessments as input only.
