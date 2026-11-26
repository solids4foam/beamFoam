# Contributing to beamFoam

Thank you for considering contributing to **beamFoam**.

## Reporting bugs and suggesting features

Open a [GitHub issue](https://github.com/solids4foam/beamFoam/issues) with a
short, descriptive title. For a bug, give the steps to reproduce it, what you
expected, the relevant log or error message, and your environment:

```plaintext
- OS version: Ubuntu 22.04.5 LTS
- OpenFOAM version: v2312
- beamFoam commit: <git rev-parse --short HEAD>
```

For a feature request, start the title with **Feature request:** and explain
what the feature is for.

## Submitting code changes

1. Create a branch from an up-to-date `main` (fork the repository first if you
   do not have write access).
2. Make the change. Follow the
   [OpenFOAM coding style](https://openfoam.org/dev/coding-style-guide/) in
   the lines you add or change, and the patterns in the surrounding code. Keep
   unrelated clean-up out of the change.
3. Build with `./Allwmake -j` and run the tutorials the change affects. If the
   results should not change, compare them with `main`.
4. Fix file headers with `bin/beamFoamUpdateHeaders <files>` and check
   Markdown with `npx markdownlint-cli2 "**/*.md"`. CI runs both checks.
5. Push the branch and open a pull request into `main`. Fill in the pull
   request template: what the problem was, what you changed, what you tested
   and what you did not.

Keep a pull request in draft until the build passes and the "Before review"
checklist in the template is done.

## Working with AI coding agents

You may use AI coding agents such as Claude Code or Codex. The repository
gives them the same instructions:

- `AGENTS.md` holds the key rules and is read at the start of every session
  by Codex and other agents. `CLAUDE.md` imports it for Claude Code.
- `.agents/skills/beamfoam/SKILL.md` is the full guide.
  `.claude/skills/beamfoam/SKILL.md` is a link to the same file, so Claude
  Code and Codex read the same guide.
- `src/wireBunchingModels/numerics/README.md` lists the shared helpers, so
  that agents and people reuse them instead of writing new ones.

Edit `AGENTS.md` and the file in `.agents/skills`; the link and the import
keep both agents in step. Keep personal agent settings, such as
`.claude/settings.local.json`, out of the repository.

When an agent writes or reviews part of a change:

- You are responsible for the change. Read every changed line, and check the
  numerical results yourself; do not rely on the agent's summary.
- Ask the agent to state what it verified and what it did not, and copy that
  into the pull request.
- Keep the change as small as the task needs. Reject unrelated refactoring,
  renaming or reformatting.
- Do not paste private data, credentials or unpublished results into a hosted
  agent.
- Say in the pull request that an agent was used, as the template asks.

## Questions

For usage questions, use the
[CFD Online forum](https://www.cfd-online.com/Forums/openfoam-cc-toolkits-fluid-structure-interaction/).
For development discussions, open an issue or discussion on GitHub.
