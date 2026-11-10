---
name: plain-language
description: Rewrite text for a reader who is fluent but not native in English — short sentences, one idea each, active voice, no idioms or phrasal verbs — while keeping every technical and domain term untouched. Use for papers, protocols, vg documentation and wiki pages, help text, release notes, or email written for an international audience, and whenever asked to make writing clearer, simpler, plainer, or easier to translate. Not for making text sound less AI-generated (that is the humanizer skill), and not for explaining concepts more simply — only the language is simplified, never the content.
---

# Plain language

The reader is an expert in the field and cannot be helped by explanation. The reader
is not a native speaker of English and is slowed down by the prose. Fix the prose.
Change nothing else.

## The rules

Write for a reader who is fluent but not native in English.

- **Keep all technical and domain terminology.** Do not simplify or paraphrase
  field-specific terms.
- Simplify everything around them: **short sentences, one idea per sentence, active
  voice, subject early.**
- Avoid **idioms, metaphors, phrasal verbs, and cultural references.**
- Prefer common words over rare or literary ones ("use" not "leverage", "start" not
  "commence").
- Define an unavoidable uncommon term in a parenthetical the first time.
- **Do not explain concepts more simply than asked.** The reader is an expert. Only
  the language is the barrier.

## What you must not change

This is a language pass, not an edit. Treat everything below as frozen:

- Technical vocabulary, house style, and the program's own wording: subcommand
  names, option names and the glossary terms in CLAUDE.md. If vg says **snarl**,
  the text says **snarl**.
- Formatting that carries meaning — code spans, command lines, italic gene names.
- Numbers, coordinates, node IDs, path names, file names, identifiers, URLs,
  citations, figure anchors, cross-references, and section structure.
- Claims. A shorter sentence that says slightly less is a rewrite failure, not a
  simplification. Splitting must preserve every qualifier.

## Procedure

1. **Read the whole document first.** Sentence splitting depends on what the next
   sentence is doing.
2. **Rewrite sentence by sentence.** For each long sentence, find the ideas and give
   each one a sentence. Put the actor first and the verb next to it.
3. **Hunt the specific failure modes** — these are what actually trip a non-native
   reader, in rough order of damage:
   - A subject and its verb separated by an aside.
   - A verb separated from its complement by an aside ("copies every reference
     path — along with any haplotype paths the GBWT holds — into the new graph").
   - Phrasal verbs, especially ones with a non-literal sense: *carry across*, *turn
     on*, *land on*, *bring along*, *throw away*, *go stale*, *layer on*.
   - Metaphors that look like plain prose: *blank slate*, *stress test*, *payload*,
     *confidence layer*, *the truth*, *second look*, *natural endpoint*.
   - A pronoun ("this", "it") whose referent is a whole preceding clause.
   - Semicolon chains holding three or more independent claims.
4. **Keep an aside when the aside is doing real work.** Do not strip em dashes
   mechanically. A parenthetical that marks genuinely subordinate detail is easier to
   read than a flat sentence sequence that gives everything equal weight.
5. **Run the variety pass** (below).
6. **Verify** (below).

## The variety pass — do not skip this

The rules above, applied without judgment, produce monotony. Uniform short sentences
with the subject early means the same subject opens sentence after sentence. Measured
on a real 10,000-word rewrite: mean sentence length fell 25.1 → 16.6 words, but
distinct sentence openers per sentence fell 0.35 → 0.26, and "The" grew from 20% to
24% of all openings.

After the rewrite, reread and ask: *where has this become a list of short declaratives
that all start the same way?* Fix those places by varying the opening — a subordinate
clause of three or four words, a fronted prepositional phrase — without rebuilding the
long sentences you just took apart. Short and varied beats short and identical.

Then do one honest audit pass: *what did splitting cost here?* Answer it before
declaring the work done.

## Verify

Confirm that only the language moved:

```bash
# code spans and command-line options, with counts
diff <(grep -oE '`[^`]+`|--[a-z][a-z-]+' old.md | sort | uniq -c) \
     <(grep -oE '`[^`]+`|--[a-z][a-z-]+' new.md | sort | uniq -c)

# numbers, coordinates and PanSN path names
diff <(grep -oE '[A-Za-z0-9_.-]+#[0-9]+#[A-Za-z0-9_.-]+|[0-9][0-9,.]*' old.md | sort | uniq -c) \
     <(grep -oE '[A-Za-z0-9_.-]+#[0-9]+#[A-Za-z0-9_.-]+|[0-9][0-9,.]*' new.md | sort | uniq -c)
```

Expect the prose to grow by roughly 4–5%. Splitting a sentence repeats its subject.
If it grew much more, content was added; if it shrank, content was lost.

## Relationship to the humanizer skill

They overlap on active voice, plain vocabulary over Latinate, no filler, no hedging,
and fewer em dashes. **They contradict each other on rhythm.** The humanizer treats
uniform sentence length as a sign of soulless, machine-like writing and asks for
varied rhythm, opinions, first person, and "some mess." This skill asks for the
opposite, because uniformity is what makes a second-language reader faster.

Pick by audience, not by taste. A blog post or an announcement wants the humanizer. A
protocol, a methods paper, or a help page written for an international readership
wants this skill. Two of the humanizer's rules are actively wrong for technical
documentation and should never be applied here: boldface reduction where bold marks
command or option names the reader must find, and de-hyphenating compound modifiers
("single-base precision" is not a style tic).
