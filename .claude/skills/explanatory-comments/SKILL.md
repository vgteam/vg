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

## Grice's maxims

A reader assumes the writer says as much as is needed and no more (quantity), only what is true
(quality), only what is relevant (relation), and plainly (manner). Every break makes the reader
infer a meaning the writer did not intend:

- "exact tie" implies that inexact ties exist; "strictly greater" implies that "greater or equal"
  was a live reading. Keep a qualifier only where the other reading is plausible.
- A stray "usually", "in general" or "effectively" implies exceptions the text never names.
- A negative-sense phrase ("nothing inside it is called") reads as a universal law; state the
  specific fact ("its nested sites are not called").
- Emphasis on exactness, completeness or guarantees implies that something nearby is approximate,
  partial or unguaranteed.

Write each word because the reader needs it. In step 3 of the protocol the reader applies the
maxims too, and reports every break with the meaning a reader would infer from it.

## Complexity budget

A document, or a top-level section, is too big to read in one sitting, and should be split, when it
exceeds any of these. The concept map from step 1 gives most of the counts.

- **New terms defined:** 15.
- **New symbols:** 10, in one notation table.
- **Sections, at every heading level:** 15.
- **New ideas per section:** 5.
- **Length:** about 4,000 words, or twenty minutes of reading.
- **Independent models or processes:** one. A model is independent when it has its own inputs,
  state and notation, and could be replaced without rewriting the others.

Before writing, and again before handing the document over, predict the reaction of a reader who
would rather be doing their own work. If it is dread at the length, or "which part do I need?",
the document is too big whatever the counts say.

Split along the independent models, not by length: one document per model, plus a framework
document that says how they fit together, owns the shared vocabulary, and states the reading order.
Each piece must be intelligible on its own given that vocabulary, does not redefine it, and touches
the others at a few named, linked points (what it takes in and what it hands on). A piece that
adapts a published model opens by saying how it differs from it. Run the protocol on each piece
separately.

## Protocol

1. **Build a concept map.** List the concepts in the change, lowest level first, with the module
   that implements each and the concepts it builds on. It should be close to acyclic. Where the
   code's structure departs from the map, rename or move code, or record why not.
2. **Extract the doc comments in concept order**: every doc comment the change adds or touches,
   grouped by concept-map row, each prefixed only by `file:line symbol`. A script that does this
   should be rerunnable.
3. **Explain back.** Give the extracted comments, and the Markdown doc if there is one, to an agent
   with no other context. Ask it to explain each concept in its own words and to list every
   question, contradiction, undefined term, sequencing puzzle, suspected error and break of a
   Grice maxim (with the meaning it would infer), each citing the `file:line` it came from.
4. **Check against the code.** For each explanation that is wrong and each real question, find the
   truth in the code and fix the comment, the doc, or the code (a rename, a split, removing dead
   code). Behaviour must not change unless that is the task: gate with a byte-identity check of the
   program's output against a binary built before the change.
5. **Repeat** from step 2 until the agent's explanation is right and its questions are about the
   domain rather than about the prose.
