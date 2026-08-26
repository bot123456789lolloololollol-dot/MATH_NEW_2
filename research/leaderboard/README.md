# Leaderboard

Per-session self-assessments live here as `session_S{n}.json`. The coordinator merges them into
`LEADERBOARD.json` with coordinator-adjusted scores after the adversarial round.

Schema for session files:
```json
{"session": "S1", "candidates": [{"id": "C-S1-01", "title": "...", "scores": {"novelty": 0, "usefulness": 0, "provability": 0, "reproducibility": 0, "improvement": 0}, "justification": "one line per factor"}]}
```
Scores are products (see ../coordination/PROTOCOL.md). Self-scores are inputs only; the coordinator's
post-adversarial scores are authoritative.
