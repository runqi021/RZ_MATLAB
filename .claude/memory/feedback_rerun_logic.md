---
name: Keep rerun logic simple
description: When user asks for a rerun option, just delete and redo — no validation checks
type: feedback
originSessionId: 6047df83-9fd0-4783-b4b4-f85d8119b32c
---
When user asks for a "rerun" flag, implement it as unconditional: delete existing output and redo. Don't add smart checking (count files, compare params, verify completeness).

**Why:** User finds conditional skip/check logic frustrating — "if I said re run just re run". Checking defeats the purpose of a rerun flag.

**How to apply:** Rerun flags should: (1) delete existing output, (2) regenerate. No counting, no param comparison, no partial-completion detection.
