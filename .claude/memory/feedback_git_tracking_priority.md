---
name: don-t-pester-about-git-tracking-commits
description: "User doesn't care about meticulous git hygiene; what matters is the current working code on disk being good"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: be23b187-1863-428f-9dbd-81b57a9328c2
---

Don't alarm-bell about files being "untracked," "uncommitted," or "at risk." If a file is on disk and works, that IS the source of truth for this user. "Untracked" ≠ "lost" or "in danger."

**Why:** User is a neuroscientist, not a CS engineer. Doesn't get paid by commit count. Cares about whether the current code on disk is good and runs correctly, not about git bookkeeping. My "today's work is gone if you rm the folder" framing got pushback for being alarmist and condescending.

**How to apply:**
- Answer git status questions factually (yes/no pushed, what's tracked) without urgency framing
- Don't offer to commit/push unless the user asks
- Don't list commit hygiene options as if they're consequential decisions
- Don't suggest `.gitignore` cleanup, branch hygiene, or commit messages unprompted
- If the user explicitly asks to commit/push, do it cleanly without lecturing
