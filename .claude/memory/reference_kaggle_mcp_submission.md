---
name: reference-kaggle-mcp-submission
description: Kaggle MCP is configured in D:\RUNQI\cmi-flu-2026/.mcp.json for automated competition submission; how it works and the two-call flow
metadata: 
  node_type: memory
  type: reference
  originSessionId: 80382f7a-b8ba-4b02-829c-4d0d6a4f091a
  modified: 2026-07-29T18:38:32.571Z
---

**Automated Kaggle submission is set up** (2026-07-29) for the CMI-Flu project. Full
detail: `D:\RUNQI\cmi-flu-2026\docs\kaggle_mcp_submission.md`. See
[[project-cmi-flu-2026-challenge]].

**Config:** `.mcp.json` at the flu repo root (committed, no secret) runs
`npx -y mcp-remote https://www.kaggle.com/mcp --header "Authorization: Bearer ${KAGGLE_TOKEN}"`.
`KAGGLE_TOKEN` is a Windows **User**-scope env var holding the `KGAT_` token from
`~/.kaggle/access_token`; `${KAGGLE_TOKEN}` is expanded by Claude Code at launch so the
secret never enters a tracked file. Re-set after rotation with
`[Environment]::SetEnvironmentVariable('KAGGLE_TOKEN', (Get-Content ~/.kaggle/access_token -Raw).Trim(), 'User')`.

**How it works:** `mcp-remote` is a stdio↔HTTP bridge — Claude Code talks MCP over
stdin/stdout to a local npx process, which forwards JSON-RPC to Kaggle's server with the
auth header. Kaggle hosts the tools; nothing is installed but the bridge.

**CRITICAL GOTCHA: `.mcp.json` loads only at Claude Code STARTUP and only for the project
directory containing it.** Claude Code must be launched from `D:\RUNQI\cmi-flu-2026`. A
session started elsewhere (e.g. RZ_MATLAB) will NOT have the Kaggle tools, and
`ToolSearch` for them returns nothing.

**Verified:** POST to `https://www.kaggle.com/mcp` returns HTTP 200,
serverInfo `Kaggle.Web 1.0.0.0`, 30 tools.

**Submission = two calls, args always wrapped in a `request` object:**
1. `start_competition_submission_upload` — `{competitionName, fileName, contentLength,
   lastModifiedEpochSeconds}` → returns a blob token
2. `submit_to_competition` — `{competitionName, blobFileTokens, submissionDescription,
   sandbox}`. **`sandbox: true` validates format WITHOUT consuming a daily submission** —
   always use it for a format check first.
3. `search_competition_submissions` — `{competitionName, sortBy:"Date", group:"All"}` to
   read the resulting public score.

Competition name: `cmi-flu-first-prediction-challenge`.

**Operating rule: never submit without an explicit instruction naming the file.** 5/day,
2 final. Variants in `submissions/` differ materially — `part1_v4_bcf` and `part1_v4_b2`
rank Task 1.3 almost oppositely (rho = −0.366). Use the registry `experiment_id` as
`submissionDescription` so the leaderboard entry maps back to its commit and CV score.
