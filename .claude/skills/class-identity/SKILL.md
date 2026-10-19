---
name: class-identity
description: Make every C++ class or struct something a reader can picture as a part in a machine or a character in a play, and name, build and document it as that thing. Use when adding a class or struct, splitting one up (including to satisfy a linter), or reviewing a design.
---

# Classes are things

People understand object-oriented code as a machine or a play: objects are parts or characters
that hold things, do jobs, and ask each other for help. A class that only collects the functions
and data for some feature cannot be pictured, so a reader cannot guess what it will do, what it
knows, or where new code belongs.

## Say what it is

1. **Finish the sentence "A `Name` is a ___ that ___."** The first blank is a concrete noun for a
   part or a character: a directory, an index, a cursor, a ledger, a queue, a builder, a judge, a
   view, a translator, a record, a table, a set. The second blank says what it holds or does.
   "Helpers", "support", "utilities", "logic", "manager of X stuff", "context" and "data" do not
   fill the first blank. The sentence becomes the class's summary comment. Say what the thing is,
   not what it is not: "files that are not directories" names no thing a reader can look for,
   and "regular files" does.
2. **Name it as that thing.** The name is the noun from the sentence, qualified as needed
   (`PriceIndex`, `PageCursor`). Name a collection as one object (`OrderQueue`, `RateTable`), not
   with a bare plural (`Orders`), which names its contents. The noun must mean the thing in its
   ordinary dictionary sense, or be the field's established term for it. Check the word's senses:
   a class named `Bench` makes the reader choose between a seat, a workbench and a court; the
   field's terms for a thing that runs code under test are "fixture", "harness" and "testbench".
3. **Describe what it knows by the fact it holds, not by the feature that first needed it.** A
   table of each customer's home branch is a "branch directory", even if the first caller is the
   overdue-notice feature. Other code will look for the fact; it will not look for the feature.

## Give it what fits

4. **Give it only methods the thing could do or be asked.** Phrase each as a request to the thing
   ("directory, which branch is customer 7 at?"). A method that does not fit belongs to another
   object, or shows that the noun is wrong. Each method is one operation that one sentence can
   explain. A library's `shelf_of(handle)` that returns a book's shelf for a book handle, and the
   shelf a bookmark's book sits on for a bookmark handle, is two operations under one name, and
   a caller cannot predict it; make `shelf_of_book()` and `shelf_of_bookmark()`.
5. **Give it only state it uses.** Build state that describes a fixed input in the constructor,
   so that a `const` object never changes and needs no locks. An object whose job is to change,
   such as a builder or a cursor, changes only through its own methods. A lazily filled `mutable`
   cache needs a speed measurement that justifies it; in existing code, keep an unmeasured one and
   add a `TODO:`. Do not keep a member only to hand it back: a
   reference to the input file, kept just for a `file()` getter, is the caller's to keep. A plain
   record is the exception; holding its fields is its job.
6. **State what it promises about everything it hands out,** in the class comment: "Every handle
   a catalog returns names a book in the catalog." Enforce the promise where outside values come
   in, by checking and throwing or by a stated precondition, because otherwise every caller must
   wonder whether a handle is real.
7. **Make parallel data one thing.** Two containers indexed by the same key (`prices` and
   `discontinued`, both by product ID) describe one record per key. Make the record a struct and
   the container a table of records.
8. **Use the interface and conventions you already have.** Before wrapping an existing interface
   in new functions or classes, check whether callers can use the interface directly. A layer
   needs a job the interface does not do. Prefer the interface's callback style
   (`for_each_item(f)`) over building a vector of results to loop over afterwards. State new facts
   in the forms the code already uses: where the code gives direction as an `is_reverse` flag, a
   new "along" and "against" pair is a second idea for the reader to learn.

## Split and group

9. **Split along parts, not along method lists.** When a class is too big, or a linter counts too
   many methods, look for a real part of the machine that the big object can hand work to (a
   parser that hands tokenizing to a tokenizer). Dealing the methods into two piles makes two
   things nobody can picture. If no real part exists, simplify the design, or leave the warning
   and say why it remains.
10. **A bundle is not a thing unless the bundle has a job.** Two objects with behaviour that are
   used together ("pair", "bundle", "context") do not make a class. Write a function that builds
   them and returns a `std::pair`, or takes them as arguments, or find the role the bundle plays.
   A record of settings (`ParserOptions`) is a plain record, not a bundle.
11. **Zero, one, or many interchangeable things.** Where any number of interchangeable things
    could take part, such as implementations under test, codecs or backends, hold a collection
    and make the code work for any number. A test harness with one member for each of two
    implementations (`file_store`, `memory_store`) must change every time an implementation is
    added. Two fields with different roles, such as a range's `begin` and `end`, are fine.
12. **Facts about a type you cannot give methods to** (a typedef of `std::vector`, a generated
    message class) go in free functions named for that type's idea
    (`is_empty_basket(const Basket&)`), or in a class that wraps the type and is itself a thing.
13. **Document each member that is an object of a class with behaviour** (not a plain container
    or number) with the part it plays for its owner: "Looks up each customer's home branch for
    the scheduler." Give each nested class a summary sentence of its own, as rule 1
    says, because a nested class with none reads as a bag of its owner's helpers.
14. **Group by concept.** Within each access section, put members and methods that describe the
    same part of the object together (all the position fields, all the lookup methods), not in the
    order some method happens to use them. Across files, keep one idea in one place: a
    representation, such as an encoding of handles into bits, goes with its encode and decode
    functions and the one comment that describes the scheme. Pieces of one job strewn over
    several files, such as helpers that each write part of a text format, become one class that
    does the job and documents the format, because a reader cannot learn a scheme from fragments.

## Redesigning

- A rename fixes a class whose members already describe one thing under a bad name. When the
  members do not describe one thing, they and the class's scope have to change; a better noun on
  the same grab bag still fails.
- A name must not claim more than the object holds. A "branch directory" must know every
  customer's branch, not only some customers'.
- Keep every capability the old classes had, or delete one on purpose and say so in the commit
  message.
- Where the old intent is unclear, keep the behaviour and add a `TODO:` asking the question.

## Check a design

For each public method that does real work, act out, in your notes, a scene of three to six lines in which objects ask each other for
things: "The scheduler asks the directory which branch customer 7 uses. The directory says the
north branch. The scheduler asks the north branch's calendar for a free slot." If a line has to
say "the helper does the X part", or an object answers a question it could not plausibly know, the
class has no identity yet: rename it, move methods, or redraw the parts. A plain record with no
operations needs no scene; check instead that its fields describe one thing.
