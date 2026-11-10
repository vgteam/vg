# Maintainability lens

Ask only: will the next person to change this code understand it and change it
safely? These findings rank lowest unless they hide a bug.

- **Fit with the code around it.** Does the change follow the patterns of the
  file and module it is in? Does it reimplement something vg or libhandlegraph
  already has (`reverse_complement`, `topological_sort`, `copy_graph`, path
  metadata parsing, the subcommand `Logger`)? Search before reporting.
- **Classes and functions.** Each class should be something a reader can
  picture; apply the `class-identity` skill to new or split classes. A function
  that needs a paragraph to explain what it does probably does too much.
- **Names.** Use the terms in the CLAUDE.md glossary for the ideas they name. A
  synonym makes a reader look for a second idea.
- **Dead code.** Unused functions, parameters, includes or options left by the
  change. Commented-out code.
- **Tests.** Is the new behavior tested at the level it lives at: a unit test
  for a class, a Bash-TAP test for a subcommand option? Is the test small and
  quick? Does it check output, not just exit status?
- **Change size.** Is unrelated reformatting or renaming mixed into the diff,
  hiding the real change?

Comments get their own lens (`comments.md`).
