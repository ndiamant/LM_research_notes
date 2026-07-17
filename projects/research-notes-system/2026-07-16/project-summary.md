---
title: Research Notes System Project Summary
date: 2026-07-16
project: research-notes-system
agent: Codex
status: draft
sources:
  - ../../../AGENTS.md
  - /Users/ndiamant/.codex/skills/write-research-note/SKILL.md
tags:
  - project-summary
  - research-notes
  - agent-workflow
---

# Summary

This project establishes a durable repository for agent-written research notes. Notes are organized by project and by the date they are added, with repo-level instructions defining the note format, linking conventions, indexing expectations, and git workflow.

# Key Points

- Research notes live under `projects/<project-slug>/<YYYY-MM-DD>/<note-slug>.md`.
- Each project has a `README.md` index linking to its notes in reverse chronological order.
- Notes should use structured Markdown, YAML frontmatter, source lists, and GitHub-compatible relative links to related notes.
- The global `$write-research-note` skill provides a reusable workflow for adding notes across projects while deferring repo-specific rules to `AGENTS.md`.
- Completed modifications should be committed to `main`; agents should attempt to push but report push failures without blocking completed local work.

# Details

The repository is intended to become a shared memory layer for research work done by agents. The root `AGENTS.md` acts as the local contract: it tells agents where notes go, how they should be formatted, how to link related notes, how to update project indexes, and how to handle commits and best-effort pushes.

The accompanying user-level Codex skill, `$write-research-note`, makes the workflow reusable outside this repository. Agents in other projects can invoke the skill to create durable research notes, provided they have access to this notes repository through a writable local checkout or another approved mechanism such as GitHub access.

# Related Notes

None yet. This is the first project note in the repository.

# Open Questions

- How should other project agents be granted routine write access to this notes repository?
- Should each external project include a short `AGENTS.md` pointer to `$write-research-note` and this repository?
- Should the repository later add templates or scripts for note creation, validation, or link checking?

# Sources

- [Repository instructions](../../../AGENTS.md)
- `$write-research-note` skill at `/Users/ndiamant/.codex/skills/write-research-note/SKILL.md`
