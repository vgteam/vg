---
name: intelligible-comments
description: Rewrite the code comments you wrote or changed so that a programmer new to the code understands each one on one reading. Decide what each comment should say (the idea of the code, the reasons for it, the guarantees its shortcuts rely on), keep only the facts a reader needs, and write them plainly, one idea per sentence. Use before finishing any task that wrote or changed comments.
---

# Comments a newcomer can read once

Your reader has not been in your session. They do not know the names you gave
things while working, and they read each sentence once, left to right. A word
that refers to something they cannot find, or a sentence that needs untangling,
breaks the reader's concentration.

Go over each comment you wrote or changed, and rewrite it using this guide.

When you are done, run `lint_prose.py` from this skill's directory. It checks the comments, unit
test names and Markdown you added since the last commit for phrasings this guide rules out, and
for command-line options named outside `src/subcommand/`. Give it `--base` and a branch to check
all of a branch's work, or paths to check whole files. Fix what it reports in your own prose.

## Decide what to say

A comment gives the reader the abstraction, rules, and animating idea behind
the code.

1. **Head each function, and each block of more than a few lines, with its
   idea:** what it computes, and why it is built the way it is, in the terms of
   computer science or of the field. 
2. **Say what a condition means, not what it tests.** The code already shows
   that `forward_loops[0]` is finite. The comment says what that means for the
   thing the code models: "A finite loop distance at the first node means a
   walk can enter the chain at its start, turn around, and leave through the
   start again."
3. **Give the reason for anything a reader would not expect,** including
   requirements from outside this code. Look for reasons in commit messages,
   issues, internal thinking notes, and the contracts of the interfaces
   involved, not only in the code.
4. **Justify a shortcut with the guarantee that makes it correct.** Code that
   skips a copy, a check or a general path needs a comment saying why skipping
   is safe. How often the shortcut runs, how big things get, and who calls it
   do not make it correct, so leave them out.
5. **Keep each comment about the code it is attached to.** Facts about callers
   or other components belong with that code, if anywhere.
6. **Explain each representation once, where it is declared:** what an element,
   a bit or an index stands for. Later comments use it without explaining it
   again. 
7. **Expect to find some bugs at the comment-writing stage.** Justifying code
   often turns up a precondition the code needs or a case it gets wrong. Put a
   precondition in the doc comment, and a bug in a `TODO:` that says what goes
   wrong, and tell the person you are working for. Do not word around it by
   saying that the function works in particular contrived circumstances, when
   the person you are working for probably wanted a solution for the general
   case.
8. **When rewriting comments, decide what questions each block needs answered
   before you reread its old comment.** Treat the old comment as claims to
   check and facts you might have missed, not as text to reword. A reworded bad
   comment is still a bad comment. In particular, feel free to throw out
   irrelevant information.
9. **Check the comments as an outline.** Read only the comments of a function,
   top to bottom: they should say what it does and why, at one level of
   abstraction, without the code. Then read the code without them: each
   comment should tell the reader something not immediately obvious from the
   code. Delete a comment that only restates its line.

## Write

1. **Give the comment the form of what it sits on.**
   - A function's doc comment starts with a command: "Return the traversal
     chosen for the parent." So does a comment before a step of a procedure:
     "Collect the children in chain order."
   - A comment on a branch, a condition or a special case says what is true
     there: "We are going right from the last node of a looping chain".
   - A variable's, member's or type's comment starts with a noun phrase for
     what it holds or is: "The mappings in the window, sorted by node ID."
   - Doxygen tags answer other questions: `@return` describes the value,
     `@pre` state a condition ("If a snarl is unary."), and `@throws` states a
     condition and an exception type.
2. **Make the subject the thing that acts,** and keep ordinary English order:
   subject, verb, object. If none of the things in play makes sense to
   personify enough to be the subject, use "we". Write "We subtract the
   correction from the read term", not "The correction subtracts from the read
   term" or "Subtracted by the correction is the read term".
3. **Give each idea its own sentence.** Name the thing first, then say why or
   add a caveat, each in a sentence of its own: "Return the K widest traversals
   through the site. Order them widest-first." Do not string things together
   with lots of colons and semiconons to fit more stuff into one sentence.
4. **Use only words the reader has:** the terms in the project's vocabulary,
   names in the code, and plain English in its ordinary sense. The vocabulary
   is a glossary in the project's agent instructions (in vg, the Terminology
   section of `BOTS.md`), with a list of module terms for words that mean
   something narrower, or something else, inside one module.
   - Use a vocabulary term without defining it. The reader has the vocabulary;
     a definition in each comment only repeats it, and drifts from it.
   - In a module that has its own list, use a word in that module's sense.
   - Apply the vocabulary where a sentence is about the vocabulary's idea. A
     word can mean something else in another context: its ordinary sense, or
     a meaning that the surrounding class or module defines in its comments.
     A step of an algorithm is still a step, even though the vocabulary uses
     "step" for one visit of a path to a node. Do not rewrite a word
     that is clear where it stands just because the vocabulary gives it another
     meaning.
   - Treat a word that names a step, part or quantity, but is in neither list,
     is not a name in the code, and is not defined by the surrounding class or
     module, as a possible coinage ("the direct pass", "the slot"). Replace it
     with a vocabulary term or a code name.
   - When an idea is new and needs a name, propose the term and its definition
     for the vocabulary to the person you are working for, instead of defining
     it in a comment.
   - Name things outright; do not point at them ("that one", "the other", "the
     rest") or define them by reference to other components and their behavior
     ("the order `find_traversals()` does").
5. **Put the head noun first in long possessives:** "the strand of the last
   unaligned read", not "the last unaligned read's strand". With the owner
   first, the reader takes the owner for the subject and must reread when the
   real noun arrives.
6. **Say what things are and do.** Don't (as here) say what things aren't, or
   what they don't do, unless the reader has a really good reason to assume
   something would happen and needs to be disabused of that notion. Do not
   describe a thing by what it is not.
7. **Say directly what a thing is.** Don't first describe it by its likeness
   to another type or convention and then explain it. The reader works out
   the likeness, only to find the explanation in the next sentence. Drop the
   likeness and keep the explanation. If no explanation follows, replace the
   likeness with one.
   - Not: "The two bounds of a site, oriented like the start and end of a
     Snarl. `start` reads into the site and `end` reads out of it."
   - But: "The two bounds of a site. `start` reads into the site and `end`
     reads out of it."

## Choose the facts

**Keep only the facts worth stating.** Writing costs you nothing, so you will
think of more facts than a reader can use, and a comment that gives them all
buries the ones that matter. For each fact, ask what a reader would get wrong
without it: call the function wrongly, keep a pointer too long, miss an
exception. Keep only those, usually one to three. Leave out what the
declaration, the names or the code already show.

## Write like these

These sentences come from comments in vg and libbdsg. A maintainer accepted the
first set and rejected the second, for the reason given after each. Write like
the accepted ones.

Accepted:
- Return the largest ID in the graph, or some larger number if the largest ID
  is unavailable.
- When memory mapping, use the SPP hash tables with the YOMO allocator, so they
  live in the memory map.
- We need constructors, destructors, copy, and move because we are keeping the
  graph we are proxying for in mapped memory.
- This creates a directed graph that has the same sequence space as the
  overlaid bidirected graph.
- Returns the handle in the underlying graph that corresponds to a handle in
  the overlay
- keep track of whether we've encountered a node in two orientations
- All the sites enclosing the child, from lowest to highest.
- No allele has a finite likelihood, so every state has zero mass.
- Where the records go.

Rejected:
- A site left with no read to any other site stays unordered. (contorted order;
  "stays" is vague, too many negatives)
- The per-site call moves with the likelihoods. (wrong verb, moves how?
  shouldn't sites stay still?)
- Uniform unless the caller supplied the message entering the chain. (coined
  term "the message", "entering" is vague, weirdly negated)
- How long a batch runs on its own thread before its remaining sites become
  tasks, and how many sites a task holds. (crowded: two variables in one
  comment, give each its own)
- The header: the output's, with the read-likelihood genotyper's lines and this
  caller's. (contorted order, should be several sentences)
- The ranges sorted, with overlapping and adjacent ones merged: the same node
  IDs, in a form that in_ranges and touches can search. (should be several
  sentences, defined in terms of other code)
- Count crossings: an entry at one boundary followed by the other. (should be
  several sentences, does not coherently define a new term but also doesn't
  reference a defined term, gesturing "the other")
- Reads with no strand log-odds, or with a strand from another phase set, are
  not counted, since `phase_hets` does not affect them. (about other code,
  contorted, weirdly negated)
- AD and GL are taken from the site, so every block of a snarl reports the same
  evidence, INFO/SB lets a consumer avoid counting it more than once. (should
  be multiple sentences, was wrong level of abstraction for where it was used)
- Likewise the off-panel candidates, on the window's own nodes only, since a
  read that reaches into another window is fetched with that window too.
  (contorted, should be multiple sentences, "reaches" is vague, the reason is
  about other code)
- The two bounds of a site or chain, oriented like the start and end of a
  Snarl. (indirect: describes the bounds by likeness to another type, then
  explains them anyway in the next sentence)

## Voice

Write comments the way a hacker in a hurry would: someone who knows this code
cold and is leaving a note for the next person at 2am, with no time for
ceremony.

## Format

For any substantial function, use Doxygen format for its doc comment in the
header:

```
/**
 * Brief description of the function
 *
 * A more detailed description of the function. This can be several lines and
 * use several sentences to convey several ideas, unlike the brief description
 * which should not feature things like `:` or `;`. This explains what the
 * function does and why you might need it.
 *
 * @note If there is some particular thing to note about the function, it would
 * go here.
 *
 * @pre The input must be sorted. (Or any applicable precondition.)
 *
 * @post The database will have been initialized with all the tables. (Or any
 * applicable postcondition)
 *
 * @returns The thing that the function returns, if any.
 *
 * @throws std::runtime_error when it is Tuesday.
 *
 * @param foo Description of the foo parameter, which is the nubmer of quuxes
 * to baz.
 */
```

Only use `@pre`/`@post`/`@throws`/`@note`/`@param` when justified. If the class
is a mountain and the function is `int get_height() const;` then a simple
comment with just a brief description of "Get the height of the mountain in
meters" is perfectly fine.

This format of brief and detailed descriptions also applies to class definition
doc comments.
