# Repository Instructions

This repo stores agent-written research notes. Notes are organized by project and by the local date on which the note is added.

## Directory Layout

Research notes go under:

`projects/<project-slug>/<YYYY-MM-DD>/<note-slug>.md`

Use the local date in `America/Los_Angeles` unless the user specifies another date.

Create missing project and date folders as needed. Keep project slugs and note slugs lowercase, hyphen-separated, and descriptive.

Recommended structure:

```text
projects/
  project-slug/
    README.md
    2026-07-17/
      note-slug.md
```

## Note Format

Each note should be Markdown with YAML frontmatter:

```md
---
title:
date:
project:
agent:
status: draft
sources:
tags:
---
```

Use this body structure unless the task clearly calls for a different shape:

```md
# Summary

# Key Points

# Details

# Related Notes

# Open Questions

# Sources
```

## Linking Between Notes

Link to other notes when they influenced the reasoning, provide background, contain supporting evidence, or record a related decision.

Use GitHub-compatible relative Markdown links within the project:

```md
[related note](../2026-07-16/example-note.md)
[project overview](../README.md)
```

Prefer relative links over absolute filesystem paths. Keep links stable by avoiding unnecessary renames of note files and folders.

In a note's `Related Notes` section, include a short explanation of why each linked note matters:

```md
- [Prior experiment notes](../2026-07-16/prior-experiment.md): Baseline assumptions used for this analysis.
```

If a note depends heavily on another note, also mention that dependency in the `Summary` or `Details` section where it affects the reasoning.

## Project Indexes

Each project should have a `README.md` index. When adding a note, update the project index with a relative link to the new note.

Prefer reverse chronological order:

```md
# Project Name

## Notes

- 2026-07-17: [Note title](2026-07-17/note-slug.md)
```

## Source Handling

Prefer source-grounded notes over broad summaries. If using web sources, include links in the `Sources` section. If using local files, link to the relevant repo path when possible.

Clearly label claims that are speculative, inferred, or uncertain.

## Figures

Include small PNG figures when a plot shows a finding better than prose or a table. Do not use images for text or tables.

- Store figures in the note's companion folder and embed them with a relative link:

  ```md
  ![Validation loss by run](note-slug/loss_curves.png)

  *One sentence stating what the figure shows.*
  ```

  The companion folder is `projects/<project>/<YYYY-MM-DD>/<note-slug>/`, next to the note.
- Copy figures into the repo instead of only linking to `$SCRATCH` paths, which are purged after 90 days without modification.
- Aim for under 200 KB per image and a few figures per note. Ask the user before adding any image over 500 KB. Git keeps every binary forever, even after deletion.
- To shrink a figure, save at 100–150 dpi with a modest figure size, or quantize the PNG to 256 colors (Pillow: `Image.open(f).convert("RGB").quantize(256).save(out, optimize=True)`), which usually cuts matplotlib plots by 2–3x with no visible change. Prefer PNG over SVG; matplotlib SVGs of dense plots can be huge.
- Put a one-sentence caption under each figure stating its takeaway, so the note still reads correctly without the image.
- Record the script or command that produced each figure, per the command provenance rules.
- Do not overwrite a figure an existing note already embeds; add a new file instead.

## Agent Instruction File

`CLAUDE.md` in the repository root is a symlink to this file. Codex and Claude Code read the same contract; edit `AGENTS.md` only.

## Skill Source

The tracked source for the `write-research-note` skill lives at `skills/write-research-note/`. It is agent-neutral: Codex invokes it as `$write-research-note`, Claude Code as `/write-research-note`.

Both cross-project installs are symlinks back to the tracked source, so there is no copy to keep in sync:

```text
~/.codex/skills/write-research-note   -> skills/write-research-note/
~/.claude/skills/write-research-note  -> skills/write-research-note/
```

Editing `skills/write-research-note/SKILL.md` updates both agents immediately. If either symlink is missing, recreate it:

```bash
ln -s /home/users/diamant/repos/LM_research_notes/skills/write-research-note ~/.claude/skills/write-research-note
```

`skills/write-research-note/agents/openai.yaml` is Codex-only interface metadata; Claude Code ignores it.

## Git Workflow

On this Stanford Sherlock server, commit every completed modification to `sherlock`.

After committing, try to push `sherlock` to `origin/sherlock`. If the push fails because network access, authentication, sandboxing, or the remote is unavailable, do not block the completed work. Report the commit hash and the push failure clearly to the user.

Use small, coherent commits that describe the note or repository change. Before committing, check `git status` and avoid including unrelated local changes unless the user explicitly asks.

Do not rewrite history, amend published commits, force-push, or delete branches unless the user explicitly asks.

## Editing Rules

- Do not overwrite an existing note unless the user explicitly asks.
- Do not move or rename existing notes unless the user explicitly asks.
- Preserve useful links when editing notes.
- If a linked note is moved or renamed as part of an explicit user request, update inbound links that are easy to find with repository search.
- Keep notes concise, but include enough context that a future agent can understand the decision or finding without reconstructing the whole conversation.
