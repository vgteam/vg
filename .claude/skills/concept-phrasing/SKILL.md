---
name: concept-phrasing
description: Phrase comments, doc comments, identifiers and test names so that each idea is introduced on its own terms, by a name for what it is, rather than by what it is not or by what uses it. Use whenever writing or reviewing comments, names or test names.
---

# Introduce each idea on its own terms

A reader builds a picture from the words that first introduce an idea, and hangs every later fact
on that picture. An idea introduced by what it is not, or as a gap inside a sentence about
something else, gives the reader nothing to picture.

These rules apply to identifiers as well as to comments.

## Shapes that do not introduce an idea

Each of these is wrong where it first introduces an idea. Fix it by naming the idea with a noun for
what it is, and defining it in a sentence of its own that says what the thing is, positively.

1. **Negation.** "Nodes that bound no snarl", "a snarl that is not an ultrabubble". Name the
   positive category and say what makes a member one: "Free nodes sit directly in a snarl's net
   graph, each as a chain of its own." A new name defined by the old negation ("a free node is a
   node that bounds no snarl") has not fixed anything. When the code tests a negation
   (`type != ULTRABUBBLE`), say what the members positively are, or name the complement as a
   concept with its own description.
2. **Leftovers.** After defining one category, do not describe the rest as "the others" or "all
   other X". Name the rest, too, and say what they are.
3. **Headless noun phrases.** "What a net handle refers to", "what we reach by leaving here",
   "where a handle lies". The phrase stands where a noun for the thing should be. Say what the
   thing is: "The kind of tree element a net handle names: a snarl, a chain, a node or a sentinel."
4. **A phrase standing in for a domain concept.** "A snarl the manager owns" is a description, not
   a concept the reader knows. Either use the domain concept ("a snarl") and state the extra
   condition as a precondition, or name a real concept.
5. **A definition through its user, its later use, or a function's result.** "The data the caller
   needs", "what the record will be built from", "the stuff for the trivial chains feature", "the
   sibling_entered_by() result". Say what it holds or what it is.
6. **A definition folded into a clause inside a noun phrase.** "The endpoint of the chain that a
   handle at the given place reads into the chain through." Introduce the idea in its own sentence
   first ("A chain is entered at its start or its end."), then use it.
7. **A new term for an existing convention.** "Along or against the chain" where the domain already
   says `is_reverse` relative to the chain. Use the domain's term. A rename that changes what a
   value means (from "along" to "backward") is a code change: change every use, or do not rename.
8. **A word that already means something else here.** A chain's "bounding nodes" are its two ends,
   so the nodes between its snarls need their own name.
9. **An unbounded quantifier.** "Every comparison", "all the laws", "the known differences". Name
   the members of the set, or the property they share.
10. **A definition that is true of everything.** "A node that lies in a snarl or in the root"
    describes every node. A definition must tell its members apart from their neighbours: name
    the neighbouring categories and check that each one is excluded.
11. **A definition by a list of fields.** "Two iterators are equal if every field except the cache
    matches" restates the code. Say what the things have in common when they are equal: "Equal
    iterators stand at the same snarl of the same chain and read it the same way."
12. **The wrong kind of thing.** Name the type the reader holds: an iterator is not an element, a
    handle is not a node, a net handle is not a snarl.

Negative sentences are fine once an idea has been introduced, to state a fact a reader might not
expect ("A free node never bounds a snarl"). Special return values are facts about the return, and
read best positively: "Returns nullopt when the handle reads out of the chain."

## Procedure

1. Run the flagger from the repository root:
   - `python3 .claude/skills/concept-phrasing/flag_phrasing.py [BASE]` checks the comment lines
     changed since BASE (default `origin/master`) in C++ files under `src/`.
   - `python3 .claude/skills/concept-phrasing/flag_phrasing.py --files FILE...` checks every
     comment in the named files.

   It prints candidates, not verdicts, under these labels:
   - `negation`: shapes 1 and 2;
   - `headless noun phrase`: shape 3;
   - `definition through a user`: shape 5;
   - `unbounded quantifier`: shape 9;
   - `contract word`: words like "once" that hide a contract; say whether the thing happens on the
     first call only, or must be called only once.
2. For each comment and name you wrote, find the noun phrases that first introduce an idea. For
   each, ask: is its head a noun for a real thing in the domain? Is the thing described by what it
   is? Does its definition have a sentence of its own? Does the definition tell members apart from
   non-members?
3. Where an idea has no good name, the code may be missing a concept. If you may change the code,
   name it there (a type, a function, a variable) and use the same name in the comments.
4. Check every definition you wrote against the code, including its edge cases: empty, a single
   element, the root, a cycle that returns to its start, both orientations. A new name with a
   false definition is worse than the old phrase.
