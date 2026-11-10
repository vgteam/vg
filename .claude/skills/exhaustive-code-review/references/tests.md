# Tests lens

Ask only: are the tests the change adds or edits worth what they cost, and can
a reader tell by looking that each expects the right thing? Test code has no
tests of its own, so it must be correct by inspection. The `vg-code-review`
skill lists test shapes that reviewers have pushed back on; read it first.

## Worth it

- For each test, name the bug it would catch that no other test would. Suggest
  cutting tests with no answer: setup written as a test, a test of another
  subcommand, a test of which short letter an option got.
- Are the tests in proportion to the code? Tests longer or more elaborate than
  the code they test need a reason.
- Is the new behavior tested at the level it lives at: a unit test for a
  class or function, a Bash-TAP test for a subcommand option?
- Does each test check output, not just exit status, and would it fail without
  the change? Name important cases with no test.

## Readable

- Each case should read as a list of true facts about the expected output for
  one input: "B is not secondary", "B has score 60". Flag check logic shared
  across cases through flags, conditionals or ternaries; it hides what each
  input should give. Repeating the call and the checks per case is better.
  To share checks, use helper functions that each test one named condition.
- Look things up the way the expectation is stated (by name, from a map), not
  by a property the reader has to work out ("the one that is neither
  secondary nor supplementary").
- A check should not depend on an incidental fact, such as a tag being the only
  tag present, or a count of separators standing in for a parsed value. Parse
  the value and check it.
- Does the test name, `SECTION` name or comment say what the check actually
  checks? A test that says "score does not change" but checks "the best path
  wins" misleads.
- When a test changes subject partway through, it should be a separate test or
  say so.
- Is the logic of an elaborate test (several moving parts, a built graph)
  explained in a comment?
- Inputs: committed files or heredocs, not strings assembled in code (see
  `vg-code-review`).
