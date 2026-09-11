---
name: Always save scripts to disk before running
description: Never rely on unsaved MATLAB/Python editor buffers — commit every generated script to a file under DLC_tracking/ or similar so it survives session end
type: feedback
originSessionId: 8ad2ee8d-8d96-4c84-b84e-fab65cbda69e
---
When generating a script for the user (MATLAB or Python), **write it to a file on disk immediately**, not just describe or paste it. Past sessions lost scripts (including a whisker batch analysis over 30 folders) because they were held only in an unsaved editor buffer.

Conventions:
- MATLAB scripts → `C:\Users\Admin\Desktop\RZ_MATLAB\<subfolder>\<name>.m` (commit to repo via user or explicit git command).
- Python scripts → `C:\Users\Admin\Desktop\RZ_MATLAB\DLC_tracking\<name>.py`.
- When iterating, edit in place — don't spawn new files unless the user asks.
- After writing, optionally lint via `mcp__matlab__check_matlab_code` for MATLAB.

**Why:** User explicitly called this out after discovering a lost script: "why u dont save script". The pattern of generating code inline → user copies to editor → editor buffer unsaved at session end → script gone.

**How to apply:** First action after deciding on a script is `Write` to the filesystem. Only then run or demo. Never deliver code in a chat-only form if it's meant to be reused.
