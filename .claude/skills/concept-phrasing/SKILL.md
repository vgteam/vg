---
name: concept-phrasing
description: Phrase comments, doc comments, identifiers and test names so that each idea is introduced on its own terms, by a name for what it is, rather than by what it is not or by what uses it. Use whenever writing or reviewing comments, names or test names.
---

# Introduce each idea on its own terms

A reader builds a picture from the words that first introduce an idea, and hangs every later fact
on that picture. An idea introduced by what it is not, or as a gap inside a sentence about
something else, gives the reader nothing to picture.

These rules apply to identifiers and test names as well as to comments. An identifier cannot
hold a definition, so give it a name for the idea (a noun, or `is_` and an adjective for a
predicate), and define the idea in the doc comment of the declaration. The examples below come
from a lending library's software, where borrowers take books home and each loan has a due
date.

## Phrasings to replace

Phrasings 1 to 12 are wrong where a file first introduces an idea, in a comment or a name.
Phrasings 13 and 14 are wrong anywhere. Fix it by
naming the idea with a noun for what it is, and defining it in a sentence of its own that says what
the thing is, positively.

1. **Negation.** Bad: "books that are not reference books". Good: "Lendable books are the books a
   borrower can take home." A new name defined by the old negation ("a lendable book is a book
   that is not a reference book") has not fixed anything. When the code tests a negation
   (`status != LOANED`), say what the selected items positively are, or name them as an idea
   with its own description.
2. **Leftovers.** After defining one category, such as lendable books, do not call the rest "the
   others" or "all other books". Name them too: "Reference books stay in the building."
3. **Headless noun phrases.** The head of a noun phrase is the noun it is about: in "the overdue
   loan's fine", the head is "fine". Bad: "what a loan refers to", "where a book is". The phrase stands
   where a noun for the thing should be. Good: "the copy of a book that the loan lends out", "the
   branch that holds the book".
4. **A phrase standing in for an idea.** Bad: "a book the catalogue knows about". Use the plain
   idea ("a book") and state the extra condition as a precondition, or name a real idea ("a
   catalogued book").
5. **A definition through its user, its later use, or a function's result.** Bad: "the data the
   overdue report needs", "what the notice will be built from", "the `find_loans()` result". Good:
   "the borrower's open loans, oldest first".
6. **A definition folded into a clause inside a noun phrase.** One short clause is fine ("the copy
   that the loan lends out"); a clause that holds another clause is not. Bad: "the branch of the library
   that the borrower who took the book out returns it through". Introduce the idea in a sentence
   of its own first ("A book is returned at its return branch, which may differ from the branch
   that lent it."), then use it ("the return branch").
7. **A new term for an existing convention.** If the code already says `is_overdue`, do not call
   the same thing "late" in comments. To replace a code term that breaks one of these rules,
   rename it in the code, at every use, and then use the new name in comments. When your task is
   limited to comments, keep the code term and add a `TODO:` suggesting the rename.
8. **A word that already means something else here.** If "copy" means one physical book, do not
   also use "copy" for a duplicated database record.
9. **An unbounded quantifier.** Bad: "every check", "all the rules", "the known problems". Name
   the members of the set, or the property they share.
10. **A definition that is true of everything.** Bad: "a book that is on a shelf or on loan"
    describes every book. A definition must tell its members apart from their neighbours: name the
    neighbouring categories, and check that each one is excluded.
11. **A definition by a list of fields.** Bad: "two loans are equal if every field except the
    timestamp matches". Good: "Equal loans lend the same copy to the same borrower for the same
    period."
12. **The wrong kind of thing.** Name the type the reader holds: a loan ID is not a loan, an
    iterator is not a book. Keep a thing apart from one pass along it. A shelf has no direction,
    but a pass along it, left to right or right to left, has one, and "the third book on the
    shelf" needs one. A name for the shelf and a name for one pass along it are different names;
    do not define one and then use it as the other.
13. **Ellipsis.** Bad: "the overdue report gives reference books none", "the overdue one's". Each
    leaves out the head noun and makes the reader recover it from an earlier sentence. Good: "the
    overdue report gives reference books no due date", "the overdue loan's fine".
14. **A vague verb or metaphor for a relation.** Bad: "the loan meets the book", "the report sees
    every branch", "the borrower touches the catalogue". Say the concrete relation: which thing
    holds, reads, lends, contains or refers to which.

## Write sentences that parse on the first reading

A reader takes in a sentence once, left to right, and should not have to go back.

- Put the subject and verb early. Bad: "Each shelf of every branch the catalogue lists, open or
  closed, with its overflow carts, is counted." Good: "The count includes each branch the catalogue
  lists, open or closed. In each branch it counts the shelves and the overflow carts."
- Give each pronoun one possible referent. Bad: "The loan names the book, so check it first."
  Good: "The loan names the book, so check the loan first."
- Put at most two words before a noun; move the rest into a clause after it. Bad: "the first unreturned renewed overdue loan". Good:
  "the first overdue loan that was renewed and is still out".
- Use a term only after it is defined.

## Give each idea room

A sentence that states two facts and a cause at once is too dense, even when every word is
right. Bad: "The overdue report gives reference books none, because it looks up due dates before
loans are recorded, and it can miss a renewed book's." Good: "The overdue report gives reference
books no due date. It can also miss the due date of a renewed book." Then add the cause in its
own sentence, if a reader of this comment would act differently knowing it.

Negative sentences are fine once an idea has been introduced, to state a fact a reader might not
expect ("A reference book is never lent"). Describe a special return value by the case that
causes it: "Returns nullopt for a borrower with no loans out."

## Procedure

1. Run the phrasing flagger, a script in this skill's folder, from the repository root:
   - `python3 .claude/skills/concept-phrasing/flag_phrasing.py [BASE]` checks the comment lines
     changed since the git commit BASE (default `origin/master`) in C++ files under `src/`.
   - `python3 .claude/skills/concept-phrasing/flag_phrasing.py --files FILE...` checks every
     comment in the named files.

   It reads only comments. It prints candidates, not verdicts. Leave a `negation`, `headless noun
   phrase` or `definition through a user` candidate alone when the phrase introduces nothing new;
   judge the other labels wherever they appear. Phrasings 4, 6, 7, 8, 10, 11 and 12 have no label, so look for
   them by reading. The labels are:

   | Label | Shape |
   |---|---|
   | `negation` | 1 and 2 |
   | `headless noun phrase` | 3 |
   | `definition through a user` | 5 |
   | `unbounded quantifier` | 9 |
   | `ellipsis` | 13 |
   | `vague relation` | 14 |
   | `contract word` | "Once", which can mean "only on the first call" or "callers must call it at most once"; write the one that is true |
   | `borrowed justification` | Behaviour or an order explained as matching another component; explain it in this code's own terms instead (see the `doc-comments` skill) |
2. For each comment and name you wrote or are reviewing, find the noun phrases that first introduce an idea. For
   each, ask: is its head a noun for a real thing? Is the thing described by what it is? Does its
   definition have a sentence of its own? Does the definition tell members apart from
   non-members?
3. Where an idea has no good name, the code may be missing an idea. Name it in the code (a type,
   a function, a variable), and use the same name in the comments, as phrasing 7 says.
4. Check every definition you wrote against the code, including its edge cases: empty, a single
   element, the first and the last element, and each direction the code can read or walk things
   in. A new name with a false definition is worse than the old phrase.
5. Reread the whole file once, top to bottom, as a newcomer would. Check that each term appears
   after its definition, that no new name is defined by a negation, and that each verb relating
   two things says what one does to the other (a loan *lends out* a copy of a book; it does not
   *touch* one).
