---
name: intelligible-comments
description: Rewrite the code comments you wrote or changed so that a programmer new to the code understands each one on one reading. Decide what each comment should say (the idea of the code, the reasons for it, the guarantees its shortcuts rely on), keep only the facts a reader needs, and write them plainly, one idea per sentence. Use before finishing any task that wrote or changed comments.
---

# Comments a newcomer can read once

Your reader has not been in your session. They do not know the names you gave things while
working, and they read each sentence once, left to right. Every word that refers to something
they cannot find, and every sentence they must untangle, costs them the idea.

Go over each comment you wrote or changed, and rewrite it by this guide.

## Decide what to say

A comment gives the reader what the code cannot: the idea it carries out, the reasons for it, and
the guarantees it relies on. A reader who has those can read the code. Notes that restate
variables, indexes and cases leave the reader to rebuild them.

1. **Head each function, and each block of more than a few lines, with its idea:** what it
   computes, and why it is built the way it is, in the terms of computer science or of the field.
   "Work out which snarls and chains can reach a tip. Whether a chain can depends on its snarls,
   and whether a snarl can depends on its chains, so snarl_tips() and chain_tips() call each
   other and save each answer in tips_found."
2. **Say what a condition means, not what it tests.** The code already shows that
   `forward_loops[0]` is finite. The comment says what that means for the thing the code models:
   "A finite loop distance at the first node means a walk can enter the chain at its start, turn
   around, and leave through the start again."
3. **Give the reason for anything a reader would not expect,** including requirements from outside
   this code: "The SnarlDecomposition interface wants every parent that is not a chain traversed
   start to end." Look for reasons in commit messages, issues, and the contracts of the
   interfaces involved, not only in the code.
4. **Justify a shortcut by the guarantee that makes it correct.** Code that skips a copy, a check
   or a general path needs a comment saying why skipping is safe, checked against the source of
   the code that gives the guarantee: "int_vector::resize() keeps the old bits as a prefix of the
   new word array, so the entries need no copy." How often the shortcut runs, how big things
   get, and who calls it do not make it correct, so leave them out.
5. **Keep each comment about the code it sits on.** Facts about callers, other components or
   typical sizes belong with that code, if anywhere.
6. **Explain each representation once, where it is declared:** what an element, a bit or an
   index stands for. Later comments use it without explaining it again.
7. **Do not narrate cases that the code handles by a general rule.** A comment on what happens
   when an entry is 0, or at the last node, describes one situation and hides the rule. If the
   case matters, state the rule in the block's heading.
8. **Write down what an honest comment shows.** Justifying code often turns up a precondition the
   code needs or a case it gets wrong. Put a precondition in the doc comment, and a bug in a
   `TODO:` that says what goes wrong, and tell the person you are working for. Do not word around
   it.
9. **When rewriting comments, decide what each block needs before you reread its old comment.**
   Treat the old comment as claims to check and facts you might have missed, not as text to
   reword. A reworded bad comment is still a bad comment.
10. **Check the comments as an outline.** Read only the comments of a function, top to bottom:
    they should say what it does and why, at one level, without the code. Then read the code
    without them: each comment should tell the reader something the code does not show quickly.
    Delete a comment that only restates its line.

## Write

1. **Give the comment the form of what it sits on.**
   - A function's doc comment starts with a command: "Return the traversal chosen for the
     parent." So does a comment before a step of a procedure: "Collect the children in chain
     order."
   - A comment on a branch, a condition or a special case says what is true there: "Going right
     from the last node of a looping chain wraps around to the first node." Forced into a
     command ("Return the first node of the chain"), it hides why the branch exists.
   - A variable's, member's or type's comment starts with a noun phrase for what it holds or is:
     "The mappings in the window, sorted by node ID."
   - Doxygen tags answer other questions: `@return` describes the value, and `@throws` and `@pre`
     state a condition ("If a snarl is unary.").
2. **Make the subject the thing that acts,** and keep ordinary English order: subject, verb,
   object. Write "The correction subtracts the read term", not "Subtracted by the correction is
   the read term" or "Where the term goes is decided by the correction".
3. **Give each idea its own sentence.** Name the thing first, then say why or add a caveat, each
   in a sentence of its own: "Return the K widest traversals through the site. Order them as
   `find_traversals()` does." Do not splice clauses with semicolons or colons.
4. **Use only words the reader has:** the terms in the project's vocabulary, names in the code,
   and plain English in its ordinary sense. The vocabulary is a glossary in the project's agent
   instructions (in vg, the Terminology section of `BOTS.md`), with a list of module terms for
   words that mean something narrower, or something else, inside one module.
   - Use a vocabulary term without defining it. The reader has the vocabulary; a definition in
     each comment only repeats it, and drifts from it.
   - In a module that has its own list, use a word in that module's sense.
   - Treat a word that names a step, part or quantity, but is in neither list and is not a name
     in the code, as a possible coinage ("the direct pass", "the slot"). Replace it with a
     vocabulary term or a code name.
   - When an idea is new and needs a name, propose the term and its definition for the
     vocabulary to the person you are working for, instead of defining it in a comment.
   - Name things outright; do not point at them ("that one", "the other", "the rest").
5. **Put the head noun first:** "the strand of the read", not "the read's strand". With the
   owner first, the reader takes the owner for the subject and must reread when the real noun
   arrives.
6. **Say what things are and do.** State a negation only for a fact the reader would not expect;
   do not describe a thing by what it is not.

## Choose the facts

**Keep only the facts worth stating.** Writing costs you nothing, so you will think of more
facts than a reader can use, and a comment that gives them all buries the ones that matter.
For each fact, ask what a reader would get wrong without it: call the function wrongly, keep
a pointer too long, miss an exception. Keep only those, usually one to three. Leave out what
the declaration, the names or the code already show.

## Write like these

These sentences come from comments in vg and libbdsg. A maintainer accepted the first set and
rejected the second, for the reason given after each. Write like the accepted ones.

Accepted:
- Return the largest ID in the graph, or some larger number if the largest ID is unavailable.
- When memory mapping, use the SPP hash tables with the YOMO allocator, so they live in the
  memory map.
- We need constructors, destructors, copy, and move because we are keeping the graph we are
  proxying for in mapped memory.
- This creates a directed graph that has the same sequence space as the overlaid bidirected
  graph.
- Returns the handle in the underlying graph that corresponds to a handle in the overlay
- keep track of whether we've encountered a node in two orientations
- The sites enclosing the child: its parent, then those enclosing the parent.
- No allele has a finite likelihood, so every state has zero mass.
- Where the records go.

Rejected:
- A site left with no read to any other site stays unordered. (contorted order; "stays" is
  vague)
- The per-site call moves with the likelihoods. (vague verb: moves how?)
- Uniform unless the caller supplied the message entering the chain. (coined term "the
  message"; "entering" is vague)
- How long a batch runs on its own thread before its remaining sites become tasks, and how many
  sites a task holds. (crowded: two variables in one comment; give each its own)
- The header: the output's, with the read-likelihood genotyper's lines and this caller's.
  (contorted order, crowded)
- The ranges sorted, with overlapping and adjacent ones merged: the same node IDs, in a form
  that in_ranges and touches can search. (crowded, and says what other code does with it)
- Count crossings: an entry at one boundary followed by the other. (coined term "crossings";
  contorted)
- Reads with no strand log-odds, or with a strand from another phase set, are not counted,
  since `phase_hets` does not affect them. (about other code; contorted)
- AD and GL are taken from the site, so every block of a snarl reports the same evidence;
  INFO/SB lets a consumer avoid counting it more than once. (two ideas spliced with a
  semicolon; coined term)
- Likewise the off-panel candidates, on the window's own nodes only, since a read that reaches
  into another window is fetched with that window too. (contorted; "reaches" is vague; the
  reason is about other code)

## Voice

Write comments the way a hacker in a hurry would: someone who knows this code cold and is
leaving a note for the next person at 2am, with no time for ceremony.
