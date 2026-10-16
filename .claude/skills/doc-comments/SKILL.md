---
name: doc-comments
description: Write every C++ doc comment as a short summary, plus Doxygen-tagged entries for the details, parameters, return value, preconditions, postconditions and exceptions a caller needs, each checked against the code. Use when writing or reviewing doc comments on files, classes, structs, functions, methods, data members or test helpers.
---

# Doc comments

A doc comment is reference material. Its reader arrives at one declaration with one question (what
is this, may I call it here, what will I get back) and leaves when the question is answered. Give
each kind of question its own place in the comment, so the reader can find the answer
without reading the rest.

## Layout

When a declaration needs only a summary, or a summary and a few plain sentences, write `///`
lines with no tags:

```cpp
/// Get the number of items waiting in the queue.
size_t size() const;
```

A declaration needs more than a summary when a caller could misuse it, or be surprised by it,
knowing only the summary and the declaration, and the tags below would give the answer a place.
Then write a `/** */` comment that starts with
`@brief`, and add Doxygen tags in this order:

```cpp
/**
 * @brief Try to take the oldest order from the queue, waiting up to a time limit.
 *
 * @details Packers on several threads share one queue.
 *
 * @param timeout How long to wait. Zero means do not wait.
 * @return
 * - The oldest order, when the queue holds one or one arrives in time.
 * - nullopt, when the time limit passes first.
 *
 * @pre timeout is not negative.
 * @post Each order is returned by exactly one call, even when several threads take orders.
 * @throws std::logic_error If the queue has been closed and is empty.
 * @note A closed queue still hands out the orders already in it.
 */
std::optional<Order> take(std::chrono::milliseconds timeout);
```

Leave out any tag that would say nothing the reader needs.

- **@brief**: one short sentence that says what the thing is or does, and nothing else. For a
  function, an imperative verb phrase that says what it does and names the kind of thing it
  returns, if any ("Take the oldest order from the queue"). For a class, struct or data member, a noun phrase saying what it *is* ("A
  queue of orders waiting to be packed"). It covers the whole declaration, not one side fact
  about it. Why the thing exists and how it works go in `@details`. A
  summary that also says why ("A queue of orders, so that packers can work in parallel") is two
  sentences squeezed into one.
- **@details**: what a caller needs that the summary did not give. Lead with the unifying idea,
  then the specifics. Do not translate the body into English.
- **@param**: one per parameter whose meaning or allowed values are not plain from its name, its
  type and the summary. Say which values are allowed: may it be null, empty, or zero? For a
  reference or pointer that the object keeps, say that it is kept and must outlive the object.
  For a parameter that is ignored, say so.
- **@return**: what the value means, including every special value it can take (null, nullopt,
  empty, an end iterator) and what each one means. When the result depends on the kind of input,
  give each case its own bullet, and make the `@brief` name a kind of thing that covers every
  case. For a collection or a series of callbacks, say what order it comes in (see "Orders"
  below).
- **@pre**: one per condition the caller must guarantee that the types do not enforce, such as a
  non-null pointer that is dereferenced, including conditions that are asserted or whose breach
  is undefined. A condition whose breach throws is a `@throws` entry instead. Limits on one
  parameter's value can go in its `@param`.
- **@post**: one per guarantee about the result or the object's state that a caller may rely on
  and could not assume from the summary: order, uniqueness, the lifetime of a returned reference.
  Put the surprising *non*-guarantees here too, such as an arbitrary order. Special values of
  the result belong in `@return`.
- **@throws**: one per exception type, with when it is thrown.
- **@note**: anything else that would surprise a caller and is not a guarantee about the result.

## What each kind of declaration needs

- **File**: the `\file` comment says which part of the program the file holds and why a reader
  would open it ("Reading and writing orders in the warehouse's CSV format"). Do not list the
  classes and functions it defines; the reader can see those.
- **Class or struct**: the summary says what one instance is. Add `@details` for what it is for,
  who uses it, and the invariants that hold for its whole life, when those are not plain.
- **Constructor**: conditions on its arguments, including how long borrowed references must
  live, go here as `@param` and `@pre`, not in the class comment. A condition shared by many
  methods goes once in the class `@details`.
- **Data member**: its own doc comment, with no tags: a summary of what it holds, then one
  sentence for each invariant the storage keeps (sorted, nonempty, null meaning "not loaded
  yet"). State the invariant here, once, instead
  of at every place that maintains it.
- **Override**: when the base class documents a contract, write no comment unless this
  implementation narrows or departs from it. Then say only how it differs, in the terms of the
  base's contract: the inputs it refuses that the base allows, how it makes choices the base leaves
  open, and anything it does that the base forbids, with a `TODO:` in the same comment saying it
  should not. Do not paraphrase the base class's comment. When the base documents no contract, document the override fully.
- **Encoding or bit layout**: describe it once, next to the type, constants, and pack and unpack
  functions that implement it. Other comments refer to it by the name of its type.
- **Test helper that checks properties**: list each property as an English sentence, one `@post`
  each, so a reader knows the finite set of things checked.

## Orders

For any function that returns or visits several things, either promise an order and say what it
is in this function's own terms ("oldest first"), or say the order is arbitrary. "In the order
`load_all()` returns them" promises nothing a caller can use; even a promise to match another
function's order needs saying in this function's terms. When you rewrite a comment that
promised an order, keep the promise. If the code cannot keep it, fix the code, or say "arbitrary"
and add a `TODO:`.

## Keep it plain

- **Say what the declaration does not already say.** Parameters named `first` and `last` give a
  first and a last item; a comment saying so only adds words to skip.
- **Give each idea room.** In a doc comment, a bullet is one short sentence or one short phrase.
  Bad: "- a file, at its path; a socket, at its address, port and protocol; a pipe, at its other
  end". Good: a table, sub-bullets, or one sentence per item.
- **Explain behaviour in this code's own terms.** Bad: "Sort by name, which is what the old
  exporter did." Good: "Sort by name, so that two exports of one catalogue can be compared line
  by line." Mention another component only when this code promises to match it, and then state
  the promise as a contract item.
- **Lead a description of another component's defect with the rule it breaks.** Bad: "The CSV
  reader copies the field before checking for quotes, so quoted fields keep their quotes." Good:
  "RFC 4180 says a reader removes the quotes around a quoted field. The CSV reader does not always
  remove them." Explain how the defect arises only if the reader needs that, and after the
  rule.

## Words to avoid

- "Once", "only" or "lazily" without the contract. Say "loads the file on its first call; later
  calls do nothing", or "callers must call this at most once".
- "Every", "all" or "any" over a set the reader cannot list. "Checks every rule" must become the
  list of rules.
- Sentences the reader would assume without being told.
- History, and alternatives that were not chosen.

## When the honest comment is bad news

Writing the contract down often shows that the code is wrong: one function does two unrelated
things, a function can return something that does not exist, a member does nothing. Do not polish
the wording around such a problem.

- Unless your task forbids changing code, change it so an honest comment is simple. When
  reviewing, report the problem instead.
- If it does, document what the code really does, and add a `TODO:` saying what is wrong.
- A `TODO:` that excuses something, such as a test that checks less than it should, makes a claim
  about why. Check that claim against the code like any other.

## Check each comment

1. **Check it against the code.** Read the body of each function you document, and the callees
   it relies on. Each claim must be true. Check the edge cases: empty, a single element, the first
   and last element, and each direction the code can read or walk things in. Describe sibling functions in matching words, so that a difference in wording
   means a difference in behaviour. When you change code, check again every comment that
   describes it, including the class comment.
2. **Read it alone, without the code, and answer**: What is this? What must I guarantee before I
   use it? What do I get, including every special value? What could surprise me? Add any answer
   the comment lacks. Delete any part that answers a question nobody would ask.
3. **Check its phrasing** with the `concept-phrasing` skill.
4. **Count misreadings** when a change adds or edits more than five doc comments. A comment that a person must read two or three times either has more
   than one reading, so independent readers disagree about it, or has no reading at all, so a
   reader cannot say what it would mean in one concrete case. For a whole change, extract every
   comment it adds or edits, each with its `file:line` and the declaration under it. Give the
   extract to two new agents (each with no other context). Ask each, for every comment, to read it
   once, then give a plain paraphrase and one concrete instance (a specific input or situation, and
   what the comment says about it), or to answer "no reading" when it cannot do that without
   guessing or repairing the sentence. Then give a third agent both answers and the code, and ask
   it to list each comment that got "no reading", where the answers disagree, or where either one
   is false of the code. Rewrite each listed comment, and repeat until the list holds only
   comments you have checked are right and that the readers misread for lack of domain
   knowledge.
