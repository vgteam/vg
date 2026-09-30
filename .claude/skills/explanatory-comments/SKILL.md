---
name: explanatory-comments
description: Write or review doc comments and design docs so that a reader new to the code can follow them, and test that with a context-free agent. Use when adding or rewriting doc comments across a change, writing a Markdown explanation of a subsystem, or when a reviewer says comments are hard to follow.
---

# Explanatory doc comments

A doc comment or doc section explains one concept to a reader who is new to the project and holds
only a few ideas in mind at once. This skill gives the rules such prose follows and a protocol that
tests it.

## Rules

- **One coherent abstraction per comment or section.** It says what one thing does, succinctly, in
  terms of lower-level ideas already introduced. If it cannot, the code may need renaming or
  splitting until it can.
- **Frame before detail.** A parent section, or a class comment, first names the concepts it
  contains and how they relate; only then do subsections or member comments define them. No
  "pile of nouns".
- **At most about five ideas per level.** Group more under a named sub-concept.
- **State the definition, not the history.** No statistics, results, measurements or version
  archaeology; no disavowals ("this is not X", "never used as Y") of things the reader had no
  reason to assume; no "this runs before that and after the other" unless ordering is the point.
- **Few cross-references and few circular ones.** Explain in place what the reader needs.
- **Sentence subjects can do their verbs.** A function computes, a flag selects; a "pass" does
  not "decide" unless it is the thing that decides.
- **One name per concept, one concept per name**, in code, comments and the Markdown alike. Where a
  word is overloaded (window, block, chain, pin, strand), say which sense, or rename.

## Protocol

1. **Build a concept map.** List the concepts in the change, lowest level first, with the module
   that implements each and the concepts it builds on. It should be close to acyclic. Where the
   code's structure departs from the map, rename or move code, or record why not.
2. **Extract the doc comments in concept order**: every doc comment the change adds or touches,
   grouped by concept-map row, each prefixed only by `file:line symbol`. A script that does this
   should be rerunnable.
3. **Explain back.** Give the extracted comments, and the Markdown doc if there is one, to an agent
   with no other context. Ask it to explain each concept in its own words and to list every
   question, contradiction, undefined term, sequencing puzzle and suspected error, each citing the
   `file:line` it came from.
4. **Check against the code.** For each explanation that is wrong and each real question, find the
   truth in the code and fix the comment, the doc, or the code (a rename, a split, removing dead
   code). Behaviour must not change unless that is the task: gate with a byte-identity check of the
   program's output against a binary built before the change.
5. **Repeat** from step 2 until the agent's explanation is right and its questions are about the
   domain rather than about the prose.
