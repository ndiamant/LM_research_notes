---
title: Private dev repo kundajelab/SMFNet-dev next to the public kundajelab/SMFNet, for unreleased work with Till
date: 2026-11-09
project: smfnet
agent: Claude Code (claude-opus-5-5)
status: draft
sources:
  - Local clone ~/repos/SMFNet_public (remotes origin = kundajelab/SMFNet, dev = kundajelab/SMFNet-dev)
  - https://github.com/kundajelab/SMFNet-dev
tags: [git, github, repo-setup, collaboration, release-workflow]
---

# Summary

Unreleased SMFNet work now lives in a private GitHub repo,
[kundajelab/SMFNet-dev](https://github.com/kundajelab/SMFNet-dev), so Nathan and co-author Till
(GitHub `Till-mpimg`) can share improvements without publishing them. The dev repo keeps the full
history of the public repo `kundajelab/SMFNet`, so changes can move between the two with ordinary
`git merge` / `git cherry-pick`. In the local clone `~/repos/SMFNet_public`, pushes to the public
remote are disabled, and `main` plus all six private `nd_*` branches now track the dev remote.

# Key Points

- **History was kept, not wiped.** The public history is already public, so deleting it hides
  nothing. A shared history is what makes releasing from dev to public (and syncing public fixes
  into dev) a normal merge. A re-`git init`ed copy would have no common ancestor with public.
- **Not a GitHub fork.** GitHub does not allow private forks of public repos, so the dev repo is
  a new private repo with the existing history pushed into it.
- **Remotes in `~/repos/SMFNet_public`:** `origin` = public `kundajelab/SMFNet`, fetch only
  (push URL is `DISABLED`). `dev` = private `kundajelab/SMFNet-dev`, upstream of every local branch.
- **Branches pushed to dev:** `main`, `nd_alpha_genome`, `nd_atac`, `nd_cheri`, `nd_unet`,
  `nd_unet_fast`, `nd_unet_patch`. Five of these are checked out in worktrees (`SMFNet_atac` under
  `$SCRATCH/src`, plus `SMFNet_cheri`, `SMFNet_unet`, `SMFNet_fast`, `SMFNet_patch` under
  `~/repos`). The worktrees share the main clone's refs and config, so a plain `git push` from any
  of them goes to dev.
- **Till's access:** `Till-mpimg` has `write` on SMFNet-dev. GitHub returned no pending
  invitation, which means Till already had access, probably as a kundajelab org member
  (inferred, not checked). If so, other org members may also see the repo, depending on the org's
  default repo permission.

# Details

## Setup commands (run 2026-11-09 from `~/repos/SMFNet_public` on Sherlock)

`gh` is a Sherlock module. Nathan logged in interactively (`gh auth login`, protocol SSH). Git
traffic uses the existing `~/.ssh/id_rsa` key. The gh OAuth token (scopes `repo`, `read:org`,
`gist`, `admin:public_key`) is stored in `~/.config/gh/hosts.yml`. No other GitHub tokens were on
the machine.

```bash
ml devel gh/2.88.1
gh auth login                       # interactive, by Nathan

git remote set-url --push origin DISABLED          # no accidental pushes to public
gh repo create kundajelab/SMFNet-dev --private --description "Private dev version of SMFNet"
git remote add dev git@github.com:kundajelab/SMFNet-dev.git
git push -u dev main nd_alpha_genome nd_atac nd_cheri nd_unet nd_unet_fast nd_unet_patch
gh api -X PUT repos/kundajelab/SMFNet-dev/collaborators/Till-mpimg -f permission=push

# verification
gh api repos/kundajelab/SMFNet-dev/collaborators/Till-mpimg/permission   # -> write
gh repo view kundajelab/SMFNet-dev --json visibility                     # -> PRIVATE
```

## Till's setup

```bash
git clone git@github.com:kundajelab/SMFNet-dev.git
```

Till's clone has only `origin` = dev, so it cannot push to public by accident. The push-disable
setting lives in each clone's `.git/config` and does not carry over. If Till adds the public repo
to pull fixes from it, they should disable its push URL too:

```bash
git remote add public git@github.com:kundajelab/SMFNet.git
git remote set-url --push public DISABLED
```

Protection that doesn't depend on anyone's local setup: either don't give Till write access on
the public repo, or turn on branch protection for public `main`.

## Day-to-day and release workflow (in `~/repos/SMFNet_public`)

```bash
# pull public fixes into dev
git fetch origin && git checkout main && git merge origin/main && git push      # push -> dev

# release some or all of dev to public
git remote set-url --push origin git@github.com:kundajelab/SMFNet.git           # re-enable
git checkout -b release origin/main
git cherry-pick <commits>            # or: git merge dev/main
git push origin release:main         # or open a PR from a pushed branch
git remote set-url --push origin DISABLED                                       # disable again
```

At setup time, two untracked files sat in the public clone's working tree:
`ALPHAGENOME_IMPLEMENTATION_PLAN.md` and `ATAC_IMPLEMENTATION_PLAN.md`. They are not in either
repo. To share them with Till, commit them on a dev-tracked branch.

# Related Notes

Several earlier notes say their branch is "not pushed". That was true when they were written.
From 2026-11-09 those branches are in SMFNet-dev; the notes themselves were not edited.

- [Frozen AlphaGenome encoder](../2026-10-29/frozen-alphagenome-encoder.md): Work on `nd_alpha_genome`, now in dev.
- [Cheri blocks and Muon](../2026-11-03/cheri-muon-smfnet.md): Work on `nd_cheri`, now in dev.
- [U-Net SMF body](../2026-11-03/unet-smf-body.md): Work on `nd_unet`, now in dev.
- [U-Net training speed](../2026-11-06/unet-training-speed.md): Work on `nd_unet_fast`, now in dev.
- [Patched U-Net](../2026-11-07/patched-unet.md): Work on `nd_unet_patch`, now in dev.
- [Kernel width](../2026-11-09/kernel-width.md): Also on `nd_unet_patch`.

# Open Questions

- Does kundajelab's default repo permission let all org members see SMFNet-dev? If it should be
  only Nathan and Till, check Settings → Collaborators and teams.
- Should public `main` get branch protection as a safeguard that doesn't rely on local config?

# Sources

- https://github.com/kundajelab/SMFNet (public)
- https://github.com/kundajelab/SMFNet-dev (private dev)
- Command output from the setup session above (gh auth status, push log, collaborator permission API).
