---
name: writing-skills
description: Write a skill so that an agent with no other context can understand and follow it, and test that with new agents. Use when adding a skill, adding rules to one, or reviewing one.
---

# Skills stand alone

A skill is a `SKILL.md` file of rules, kept in `.claude/skills/<name>/` under the repository root.
It starts with a header holding its `name:` and a `description:` that says what the skill does and
when to use it. Its reader is an agent that has the skill and the other skills it names, and
nothing else: it has not seen the conversation, review or bug that led the author to write a rule.
Call that the rule's origin. The reader can follow only what the skills say, so every rule,
example and term must make sense without the origin.

## Content

- **State each rule as a general rule,** so the reader can apply it to code it has never seen.
  Bad: "Do not call the retry counter `tries2`." Good: "Name a variable for what it counts:
  `retries`, not `tries2`."
- **Use self-contained examples.** Build examples from ideas any programmer knows, such as files,
  queues, caches, or a shop's orders and stock, and explain any other idea in the same sentence.
  A name from the origin means nothing to the reader. Bad: "Don't keep a `TrackRecord*` across
  `SlotTable::grow()`." Good: "Don't keep a pointer into a growable array across a call that
  grows it, because growing can move the elements."
- **Pair each bad example with its fix,** on the same subject, so the reader sees what the rule
  changes. A rule that is easy to apply needs no example; one that is not needs a pair.
- **Give every rule its reason.** Lead with the action, then say what goes wrong without it. A
  rule with no reason cannot be stretched to a case it does not mention, so check each rule for a
  "because" or a consequence before you finish.

## Form

- **Use one name for each idea,** within a skill and across the skills it names. A reader who
  sees two names assumes two ideas. Bad: "task" in one rule and "job" in the next. Good: "task"
  in both.
- **Name another skill, by the `name:` in its header,** for an idea that skill explains, instead
  of restating the idea, because two copies drift apart. Bad: copying the rules for doc comments
  into a skill about classes. Good: "Document each class as the `doc-comments` skill says." 
- **Keep it short:** about five rules or bullets per section, and about 150 lines in all. Leave out
  rules that a careful programmer follows unprompted, such as "make sure it compiles". Readers
  skim long skills and miss rules.

## Test a skill

A new agent, in these steps, is one started with no conversation, told to read only the files you
name. It reads the skill the way a future reader will.

1. **Explain-back.** Give a new agent the skill and the skills it names, and nothing else. Ask it to restate each rule in its own words
   with a new example, to list each term, example or instruction it could not understand or could
   read two ways (with line numbers), and to list each rule the skill breaks itself.
2. **Fix the skill.** Compare each restatement with what the author meant. When you are not the
   author, ask the author, because only the author knows which reading was meant. Fix each
   misreading, unclear point and broken rule in the skill itself. Do not explain it in the test
   prompt, because the future reader will not get that prompt.
3. **Repeat** with a new agent each time, until the restatements match what the author meant and
   each remaining point is a matter of taste. A point is about meaning when two readers who
   followed the text could do different things; otherwise it is taste.
4. **Trial.** Give a new agent the skill and a small task of the kind the skill governs, such as a
   short file to write or fix. Give a second new agent the skill, the task and the result, and have
   it grade each rule pass, fail or not applicable. Make at least one task unrelated to the rules'
   origin (when you do not know the origin, any new task will do), so the skill does not fit only
   one case. After each fix to the skill, rerun the explain-back, then the trial, until every
   applicable rule passes on two different tasks.
