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

## Skill Source

The tracked source for the `$write-research-note` Codex skill lives at `skills/write-research-note/`.

The installed cross-project copy normally lives at `/Users/ndiamant/.codex/skills/write-research-note/`. When editing the tracked skill source, keep the installed copy in sync so Codex can use the latest workflow across projects.

## Git Workflow

Commit every completed modification to `main`.

After committing, try to push `main` to `origin/main`. If the push fails because network access, authentication, sandboxing, or the remote is unavailable, do not block the completed work. Report the commit hash and the push failure clearly to the user.

Use small, coherent commits that describe the note or repository change. Before committing, check `git status` and avoid including unrelated local changes unless the user explicitly asks.

Do not rewrite history, amend published commits, force-push, or delete branches unless the user explicitly asks.

## Editing Rules

- Do not overwrite an existing note unless the user explicitly asks.
- Do not move or rename existing notes unless the user explicitly asks.
- Preserve useful links when editing notes.
- If a linked note is moved or renamed as part of an explicit user request, update inbound links that are easy to find with repository search.
- Keep notes concise, but include enough context that a future agent can understand the decision or finding without reconstructing the whole conversation.
