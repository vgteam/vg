# Comments lens

Ask only: is each comment the change adds or edits true, needed, and readable
by someone new to the code, and does each new idea have one? Before starting,
read the `intelligible-comments` and `explanatory-comments` skills in
`.claude/skills/`; they are the standard, and this lens is not done without
them. The rules in CLAUDE.md about documentation apply too.

## Present

- A concept the change introduces, used in names or comments in several places
  (a kind of status that is passed on, a "group" of records), with no doc
  comment that says what it is and how it fits with the rest. Say where the
  overview belongs: usually the class or header doc comment.
- A term a test or comment uses that the code under test never defines.

## True

- Does the comment describe what the code does now? Check each claim against
  the code. A stated precondition, guarantee or rule that nothing enforces is a
  correctness finding, not a comment nit.
- Is the comment narrower or broader than the code it sits on?
- Does a doc comment say what the function guarantees and what it needs from
  its caller, or only restate its name?
- Does a comment explain code by what other code relies on ("callers expect the
  primary first") instead of by what this code guarantees ("put the primary
  first")? Rewrite it as the guarantee.
- Does a comment describe something that happens elsewhere instead of what the
  code under it does?

## Needed

- Comments that narrate the conversation or the task ("as requested", "new
  approach", "Step 2b") instead of explaining the code.
- Comments about what used to be there, or what the code does not do.
- Comments that answer a question no reader would ask: an edge case whose
  behavior is the obvious one ("a zero-length tail prunes nothing"), or a
  denial of something nobody would assume ("these coordinates are not relative
  to the path"). They make the reader wonder why the case was worth a note.
- Comments that reassure instead of explaining ("safely", "robustly",
  "to avoid silently dropping"). These often hide something nobody worked out.
- Missing comments: an invariant at a storage location with no note there, a
  non-obvious shortcut with no note of the guarantee it relies on.

## Too much

A comment can say too much as easily as too little. Look for both.

- A comment that argues for the code against alternatives it does not use, by
  saying how they would fail. Rewrite it to say what the code does and what it
  is for. Not "`unordered_map` would make the output order change between
  runs, so this uses `map`", but "Keep paths sorted by name for stable output".
  Do not ask for such a comparison where a comment lacks one.
- A doc comment that says how the function finds its answer: what it reads,
  how far it looks, what it compares. The caller needs what the function
  returns, what it requires, and when it fails. Cut the mechanism, unless a
  caller would act differently knowing it.
- A comment on a step of a test that argues why the expectation holds, or
  repeats the test's description string. Name the step in one line ("Make a
  graph with a bubble on each strand").
- A comment longer than the code under it, other than a doc comment. Look for
  what to cut.

A rewrite you suggest must be shorter than the comment it replaces, unless the
comment lacks a fact that a reader would get wrong without it. Name that fact.

## Readable

- One idea per comment; plain words; terms from the CLAUDE.md glossary.
- No AI writing patterns (inflated significance, "ensures", lists of three,
  em-dash chains). The `humanizer` skill lists them.
- `///` and `/** */` comments apply only to the next declaration.

The `explanatory-comments` skill's full protocol is for an author rewriting
their own comments; in a review, apply its rules without running the protocol.
To test a hard comment, give it with its declaration, and nothing else, to a
subagent with no context, and ask what it means. Report where its answer goes
wrong.
