---
name: doc-comments
description: Write every C++ doc comment as a Doxygen @brief summary, a @details description, and one tagged entry per parameter, return value, precondition, postcondition and exception that a caller needs, each checked against the code. Use when writing or reviewing doc comments on classes, structs, functions, methods, data members or test helpers.
---

# Doc comments

A doc comment is reference material. Its reader arrives at one declaration with one question (what
is this, may I call it here, what will I get back) and leaves when the question is answered. Give
each kind of question its own tagged place in the comment, so the reader can find the answer
without reading the rest.

## Layout

Use Doxygen tags, in this order:

```cpp
/**
 * @brief Get the net handle for a snarl.
 *
 * @details A longer description: what the result means and how it relates to the
 * neighbouring concepts.
 *
 * @param snarl The snarl to look up.
 * @return The snarl's start-to-end traversal.
 *
 * @pre The snarl belongs to the SnarlManager the view was made from.
 * @post The handle is start-to-end even when the snarl cannot be crossed from start to end.
 * @throws std::invalid_argument If the snarl is hidden from the decomposition.
 * @note Anything else that would surprise a caller.
 */
net_handle_t get_snarl_net(const Snarl* snarl) const;
```

Every doc comment has a `@brief`. Leave out any other tag that would say nothing the reader needs.
A one-line `/// @brief ...` is enough when nothing else is needed.

- **@brief**: one sentence. For a function, an imperative verb phrase saying what calling it does
  and naming the kind of thing it returns ("Get the sibling chain that a handle enters"). For a
  class, struct or data member, a noun phrase saying what it *is* ("A directory that lists each
  node's place in the decomposition."). It covers the whole declaration at one level of abstraction, not
  one side fact about it.
- **@details**: what a caller needs that the summary did not give. Lead with the unifying idea, then
  the specifics. Do not translate the body into English.
- **@param**: one per parameter whose meaning or allowed values are not plain from its name, its
  type and the summary. Say which values are allowed: may it be null, empty, or from another graph?
  For a reference or pointer that the object keeps, say that it is kept and that it must outlive
  the object. For a parameter that is ignored, say so.
- **@return**: what the value means, including every special value it can take (null, nullopt,
  empty, end iterator) and what each one means. When the result depends on the kind of input,
  give each case its own bullet instead of joining them with semicolons, and make the `@brief`
  name a kind of thing that covers every case.
- **@pre**: one per condition the caller must guarantee that the types do not enforce, including
  non-null pointers that are dereferenced. Say what happens when it is broken if that is defined
  (an exception, an assert).
- **@post**: one per guarantee about the result or the object's state that a caller may rely on and
  could not assume: order, uniqueness, orientation, lifetime of a returned reference. Also state the
  surprising *non*-guarantees, such as that a pointer may be null or that order is arbitrary.
- **@throws**: one per exception type, with when it is thrown.
- **@note**: anything else that would surprise a caller.

## What each kind of declaration needs

- **Class or struct**: `@brief` says what one instance is. `@details` says what it is for, who uses
  it, and the invariants that hold for its whole life. Give the reader a picture of the object; a
  list of facts about some other class does not describe this one.
- **Constructor**: the conditions on its arguments, including how long borrowed references must
  live, go here as `@param` and `@pre`, not in the class comment. A condition that applies to
  every method of a kind ("methods that take a graph need one with the manager's node IDs") goes
  once in the class `@details`.
- **Data member**: its own doc comment, saying what it holds and every invariant the storage keeps
  (sorted, nonempty, null meaning X). State the invariant here, once, instead of at every place
  that maintains it. A member that is an instance of a helper class gets a doc comment saying what
  part it plays for its owner.
- **Override**: when the base class documents a contract, write no comment unless this
  implementation adds to it; then write a `@brief` and only the added items: extra preconditions,
  choices the interface leaves open, and departures from the interface. Do not repeat what the
  interface already requires. When the base documents no contract, document the override fully.
  Cover every kind of input the method accepts; a comment that explains two cases out of five
  leaves the reader guessing about the other three.
- **Bit layout or encoding**: describe it once, in one place, next to the type, constants, and pack
  and unpack functions that implement it. Other comments refer to that place.
- **Test helper that checks properties**: list each property as an English sentence (one `@post` or
  one bullet each), so a reader knows the finite set of things checked.
- **TEST_CASE and SECTION names**: state a property that ought to be true ("A cyclic chain lists its
  shared node once", "Three nested snarls are presented the same way"), not a topic ("Cyclic
  chains", "Three nested snarls") or a category ("Known differences").

## Words to avoid

- "Once", "only", "lazily" without the contract: say "builds the index on its first call; later
  calls do nothing", or "callers must call this at most once".
- "Every", "all", "any" over a domain the reader cannot list: "every comparison" must become the
  list of comparisons.
- Sentences the reader would assume without being told. If one is a real contract item, put it in
  `@post` beside the rest of the contract, never as the whole comment.
- History, and alternatives that were not chosen. A departure from an interface, or from another
  implementation the reader will compare against, is a contract item and belongs in `@note`.

## When the honest comment is bad news

Writing the contract down often shows that the code is wrong: one function does two unrelated
things, a function can return a handle to something that does not exist, a member does nothing.
Do not polish the wording around such a problem.

- If you may change the code, change it so an honest comment is simple.
- If you may not, document what the code really does and add a `TODO:` saying what is wrong, so the
  reader sees the problem instead of a confident description of it.
- A `TODO:` that excuses something, such as a test that checks less than it should, makes a claim
  about why. Check that claim against the code like any other.

## Check each comment

1. **Check it against the code.** Read the body of every function you document, and the callees it
   relies on. Every claim must be true: values, conditions, which cases occur, numbers. Check the
   edge cases: empty, a single element, the root, a cycle, both orientations. A member whose
   feature never runs must not be described as having an effect. Describe sibling functions in
   matching words, so that a difference in wording means a difference in behaviour. When you
   change code, check again every comment that describes it, including the class comment.
2. **Read it alone, without the code, and answer**: What is this? What must I guarantee before I use
   it? What do I get, including every special value? What could surprise me? Where an answer is not
   obvious and the comment does not give it, add it. Where the comment answers a question nobody
   would ask, delete that part.
3. **Check its phrasing** with the `concept-phrasing` skill.

`///` and `/** */` comments apply only to the single declaration that follows them; give each
similar member its own comment.
