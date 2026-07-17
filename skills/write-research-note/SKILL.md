---
name: write-research-note
description: Create, update, and organize agent-written Markdown research notes in the LM_research_notes repository at /Users/ndiamant/repos/LM_research_notes. Use when the user asks Codex to add a research note, capture findings, summarize sources into a durable note, document a research decision, update a project note index, or preserve research context for later agents.
---

# Write Research Note

## Overview

Use this skill to add durable, source-grounded Markdown research notes to the user's dedicated notes repository. Treat that repository's `AGENTS.md` as the highest-priority local contract for directory layout, note format, linking, indexing, and git workflow.

## Target Repository

Use this notes repository by default:

```text
/Users/ndiamant/repos/LM_research_notes
```

The expected remote is:

```text
git@github.com:ndiamant/LM_research_notes.git
```

If the current task is running in a different repository, use the notes repository path above when it is available and writable. If the path is not available or not writable in the current sandbox, report that clearly. Do not silently write research notes into the current project repository unless the user explicitly asks.

## Workflow

1. Locate the target notes repository. Prefer `/Users/ndiamant/repos/LM_research_notes`; if already in that repo, use the current checkout.
2. Read the target notes repository's `AGENTS.md`.
3. Identify the target project. If the user does not specify one, infer a concise project slug from the request and existing `projects/` folders.
4. Use the repository's required date convention. If none is specified, use the local date in the user's configured timezone.
5. Create the project folder, date folder, and project `README.md` index when missing.
6. Write the note as Markdown with YAML frontmatter and the section structure required by `AGENTS.md`.
7. Link related notes using GitHub-compatible relative Markdown links when they influenced the reasoning, provide background, contain supporting evidence, or record a related decision.
8. Update the project index with a relative link to the new or changed note.
9. Review the diff for unrelated changes before committing.
10. Follow the notes repo's git workflow. Commit completed changes locally, then attempt the configured push if requested or required. If pushing fails because network access, authentication, sandboxing, or the remote is unavailable, report the commit hash and failure.

## Note Writing Rules

- Prefer concise, source-grounded notes over broad summaries.
- Preserve uncertainty. Label speculative, inferred, or unresolved claims.
- Do not overwrite, move, or rename existing notes unless the user explicitly asks.
- When editing an existing note, preserve useful links and frontmatter unless the task requires changing them.
- Use relative repo links, not absolute filesystem paths, for links intended to render on GitHub.
- Include enough context that a future agent can understand the finding or decision without reconstructing the full conversation.

## Related Notes

Before writing, search existing project notes for likely related material. Use filenames, project indexes, headings, tags, and keyword search. Link only notes that are actually relevant to the new note's reasoning or context.

In a `Related Notes` section, include a short reason for each link:

```md
- [Prior baseline](../2026-07-16/prior-baseline.md): Establishes the assumptions reused here.
```

If a related note materially affects the conclusion, also mention it in `Summary` or `Details`.

## Source Handling

When the user provides source material, cite it directly in the `Sources` section. When using web sources, include links. When using local files, link to the repo path if the file is inside the repository.

If more source gathering is needed and the user has not asked for live research, ask before browsing or making external calls unless the answer would otherwise be obviously incomplete.

## Git Handling

Before staging, run `git status` and inspect the relevant diff. Stage only files related to the note task unless the user explicitly asks otherwise.

Use small commit messages that describe the durable repository change, for example:

```text
Add note on transformer scaling assumptions
Update retrieval project notes
```

If the repository has no git remote or the push cannot complete, leave the local commit in place and tell the user exactly what happened.
