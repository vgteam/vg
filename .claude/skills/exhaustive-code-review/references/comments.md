# Comments lens

Ask only: is each comment the change adds or edits true, needed, and readable
by someone new to the code? Read the `intelligible-comments` and
`explanatory-comments` skills first; they are the standard. The rules in
CLAUDE.md about documentation apply too.

## True

- Does the comment describe what the code does now? Check each claim against
  the code. A stated precondition, guarantee or rule that nothing enforces is a
  correctness finding, not a comment nit.
- Is the comment narrower or broader than the code it sits on?
- Does a doc comment say what the function guarantees and what it needs from
  its caller, or only restate its name?

## Needed

- Comments that narrate the conversation or the task ("as requested", "new
  approach", "Step 2b") instead of explaining the code.
- Comments about what used to be there, or what the code does not do.
- Comments that reassure instead of explaining ("safely", "robustly",
  "to avoid silently dropping"). These often hide something nobody worked out.
- Missing comments: an invariant at a storage location with no note there, a
  non-obvious shortcut with no note of the guarantee it relies on.

## Readable

- One idea per comment; plain words; terms from the CLAUDE.md glossary.
- No AI writing patterns (inflated significance, "ensures", lists of three,
  em-dash chains). The `humanizer` skill lists them.
- `///` and `/** */` comments apply only to the next declaration.

To test a hard comment, give it with its declaration, and nothing else, to a
subagent with no context, and ask what it means. Report where its answer goes
wrong.
