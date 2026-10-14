---
name: class-identity
description: Make every C++ class or struct something a reader can picture as an object in a machine or a character in a play, and name, build and document it as that thing. Use when adding a class or struct, splitting one up (including to satisfy a linter), or reviewing a design.
---

# Classes are things

People understand object-oriented code as a machine or a play: objects are parts or characters
that hold things, do jobs, and ask each other for help. A class that only collects the functions
and data for some feature cannot be pictured, so a reader cannot guess what it will do, what it
knows, or where new code belongs.

## Rules

1. **Say what one instance is.** Before writing a class, finish the sentence "A `Name` is a ___
   that ___." The first blank is a concrete noun for a part or a character: a directory, an index,
   a cursor, a ledger, a queue, a builder, a judge, a view, a translator, a record, a table, a set.
   The second blank says what it holds or does. "Helpers", "support", "utilities", "logic", "stuff
   for X", "context" and "data" do not fill the first blank. The sentence becomes the class's
   `@brief`.
2. **Name it as that thing.** The name is the noun from the sentence, qualified as needed
   (`NodeDirectory`, `ChainCursor`). A collection is a thing too, but name it as one (`NodeSet`,
   `CapTable`), not with a bare plural (`TrivialChains`), which names its contents and not the
   object.
3. **Describe what it knows by the domain fact, not by the feature that first needed it.** A table
   of the snarl that contains each node is a "containment table", even if the first caller only
   asks about trivial chains. A fact about the domain will have other users; a feature's data will
   not be found by them.
4. **Give it only methods the thing could do or be asked.** Phrase each as a request to that thing
   ("directory, where is node 7?"). A method that does not fit belongs to another object, or shows
   that the noun is wrong.
5. **Give it only state the thing would hold.** Build that state when the object is made, so a
   `const` object stays unchanging and needs no locks. A lazily filled `mutable` cache needs a
   reason that fits the picture and a measured cost to justify it.
6. **Make parallel data one thing.** Two containers indexed by the same key (`snarls_by_number`
   and `hidden`) describe one record per key. Make the record a struct, and the container a table
   of records with a name of its own.
7. **Split along parts, not along method lists.** When a class is too big or a linter counts too
   many methods, look for a real part of the machine that the big object hands work to. A split
   that just deals the methods into two piles makes two things nobody can picture. If no real
   part exists, simplify the design or keep the warning.
8. **A bundle is not a thing unless the bundle has a job.** Two objects that are used together
   ("pair", "bundle", "context") do not make a class. Write a function that builds them, or find
   the role the bundle plays: a test bench that builds several decompositions of one graph and
   keeps alive the inputs they borrow is a thing; a "decomposition pair" is not.
9. **Facts about a type you cannot give methods to** (a typedef of `std::vector`, a protobuf
   message) go in free functions named for that type's concept (`is_cyclic(const Chain&)`), or in
   a class that wraps the type and is itself a thing.
10. **Give each member that is another object a doc comment** saying the part it plays for its
    owner ("Answers where each node sits in the decomposition.").
11. **Group by concept.** Order members and methods by the part of the picture they belong to,
    not by which method happens to use them.

## Redesigning

- A rename is not a redesign. When a class has no identity, its members and its scope have to
  change, not just its name; a better noun on the same grab bag still fails.
- A name must not claim more than the object holds. A "containment table" must say where every
  node is contained, not only some of them.
- Keep every capability the old classes had, or delete one on purpose and say so in the change
  description and in a `TODO:` where a caller might look for it. Check that every caller can still
  do what it did.
- Where the old intent is unclear, keep the behaviour and add a `TODO:` asking the question.

## Check a design

For each main operation of the class, act out a scene of three to six lines in which objects ask
each other for things: "The view asks the directory where node 7 is. The directory says it is the
third node of chain C, read in reverse. The view asks chain C for the snarl after it." If a line
has to say "the helper does the X part", or an object answers a question it could not plausibly
know, the class has no identity yet; rename it, move methods, or redraw the parts. A plain record with no operations
needs no scene; check instead that its fields describe one thing.
