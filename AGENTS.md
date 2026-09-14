# ALRD workspace workflow

## Canonical roles

- `/Users/zhikuan_quan/Documents/ALRD` is the development and evaluation
  workspace. Make code, documentation, manuscript, and configuration changes
  here first.
- `/Users/zhikuan_quan/Documents/GitHub/ALRD` is the GitHub-connected checkout.
  Treat it as a publication/synchronization target, not as the primary editing
  workspace.

## Required update sequence

1. Implement the requested change only in `/Users/zhikuan_quan/Documents/ALRD`.
2. Run tests and checks appropriate to the change in that workspace.
3. Report the files changed and the evaluation results to the user.
4. Ask the user whether the same evaluated update should be synchronized to
   `/Users/zhikuan_quan/Documents/GitHub/ALRD`.
5. Do not modify, copy into, commit in, or push from the GitHub-connected
   checkout until the user explicitly approves that synchronization.

Approval applies only to the specific update just reviewed. Do not assume that
approval carries forward to later updates.

## Synchronization safety

- Before an approved synchronization, run
  `scripts/compare_github_checkout.sh` and inspect both Git statuses.
- Copy only the reviewed files associated with the approved update. Do not use
  destructive mirroring or deletion commands.
- Re-run relevant tests in the GitHub-connected checkout after copying.
- Show the resulting diff before committing or pushing. A request to
  synchronize files does not by itself authorize a commit or push unless the
  user also asks for it.
- The two repositories currently have unrelated Git histories. Do not merge,
  rebase, force-push, or replace either history merely to connect them.
