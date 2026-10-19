---
name: class-identity
description: Make every C++ class or struct something a reader can picture as a part in a machine or a character in a play, and name, build and document it as that thing. Use when adding a class or struct, splitting one up (including to satisfy a linter), or reviewing a design.
---

# Classes are things

People understand object-oriented code as a machine or a play: objects are
parts or characters that hold things, do jobs, and ask each other for help. A
class that only collects the functions and data for some feature cannot be
pictured, so a reader cannot guess what it will do, what it knows, or where new
code belongs.

## Say what it is

1. **Finish the sentence "A `Name` is a ___ that ___."** The first blank is a
   concrete noun for a part or a character: a directory, an index, a cursor, a
   ledger, a queue, a builder, a judge, a view, a translator, a record, a
   table, a set. The second blank says what it holds or does. "Helpers",
   "support", "utilities", "logic", "manager of X stuff", "context" and "data"
   do not fill the first blank. The sentence becomes the class's core idea, and
   can be written down in its doc comment. Here it is particularly important to
   say what the thing is, not what it is not: "files that are not directories"
   names no thing a reader can look for, and "regular files" does. Never define
   the `Barber` class with something like "a man who shaves all those who do
   not shave themselves."
2. **Name it as that thing.** The name is the noun from the sentence, qualified
   as needed (`PriceIndex`, `PageCursor`). Name a collection as one object
   (`OrderQueue`, `RateTable`), not with a bare plural (`Orders`), which names
   its contents. The noun must mean the thing in its ordinary dictionary sense,
   or be the field's established term for it. Check the word's senses and
   whether you need to use a compound word: a class named `Bench` is probably
   some kind of metaphorical chair, whereas a `Testbench` is something else
   entirely.
3. **Describe the thing on its own terms, like the data it holds or the thing
   it represents, not in the terms of the feature that first needed it.** A
   table of each customer's home branch is a "branch directory", even if the
   first caller is the overdue-notice feature. When writing new code later, we
   will will look for the concept, not the other feature.

## Give it what fits

4. **Give it only methods the thing could do or be asked.** Phrase each as a
   request to the thing ("directory, which branch is customer 7 at?"). A method
   that does not fit belongs to another object, or shows that the noun is
   wrong. Each method is one operation that one sentence can explain. A
   library's `shelf_of(ref)` that needs several clauses to explain the handling
   of different `ref`s might be multiple operations under one name
   (`shelf_of_book()` and `shelf_of_bookmark()`).
5. **Give it only state it uses.** Build state that describes a fixed input in
   the constructor. An object whose job is to change, such as a builder or a
   cursor, should change only through its own methods. A lazily filled
   `mutable` cache is probably a bad idea. Do not keep a member only to hand it
   back: a reference to the input file, kept just for a `file()` getter, is a
   bad design. A plain record is the exception; holding its fields is its job.
6. **State important invariants** in the class comment: "All call numbers
   returned by the catalog correspond to books that exist. Call numbers for
   nonexistent books cannot be obtained." Make sure to enforce the invariants.
7. **Make parallel data one thing.** Two containers indexed by the same key
   (`prices` and `discontinued`, both by product ID) describe one record per
   key. Make the record a struct and the container a table of records.
   (Sometimes we need to use columnar storage instead, but if it can be avoided
   it makes the code more understandable.)
8. **Use the interface and conventions you already have.** Before wrapping an
   existing interface in new functions or classes, check whether callers can
   use the interface directly. For example, if the existing interface uses
   callback iterators (`for_each_item(f)`), use those instead of writing a
   facade that dumps everything to vectors. State new facts in the forms the
   code already uses: where the code gives direction as an `is_reverse` flag, a
   new "along" and "against" concept would be a second idea for the reader to
   learn, which will confuse a reader who already understood the first idea.

## Split and group

9. **Split along parts, not along method lists.** When a class is too big, or a
   linter counts too many methods, look for a real part of the machine that the
   big object can hand work to (a parser that hands tokenizing to a tokenizer).
   Dealing the methods into two piles makes two things nobody can picture. If
   no real part exists, simplify the design, or leave the warning and say why
   it remains.
10. **A bundle is not a thing unless the bundle has a job.** Two (or more)
    objects with behaviour that are used together ("pair", "bundle", "context")
    do not always make a good class. Find a role for the bundle, or treat the
    members as separate concepts.
11. **Zero, one, or many interchangeable things.** Where any number of
    interchangeable things could take part, such as implementations under test,
    codecs or backends, hold a collection and make the code work for any
    number. A test harness with one member for each of two implementations
    (`file_store`, `memory_store`) must change every time an implementation is
    added. Two fields with different roles, such as a range's `begin` and
    `end`, are fine.
12. **Facts about a type you cannot give methods to** (a typedef of
    `std::vector`, a generated message class) go in free functions named for
    that type's idea (`is_empty_basket(const Basket&)`), or in a class that
    wraps the type and is itself a thing, or in a long comment on the typedef
    that explains how to use it to represent a concept.
13. **Document each member that is an object of a class with behaviour** (not a
    plain container or number) with the part it plays for its owner: "Looks up
    each customer's home branch." Give each nested class a summary sentence of
    its own, as rule 1 says, because a nested class with no description is just
    a bag of its owner's helpers.
14. **Group by concept.** Within each access section, put members and methods
    that describe the same part of the object together (all the position
    fields, all the lookup methods), not in the order some method happens to
    use them. Across files, keep one idea in one place: a representation, such
    as an encoding of structured data into bits, goes with its encode and
    decode functions and the one comment that describes the scheme. Pieces of
    one concern strewn over several files, such as helpers that each write part
    of a text format, become one class that handles the concern and documents
    the format, because a reader cannot learn a scheme from widely distributed
    fragments.

## Redesigning

- A rename fixes a class whose members already describe one thing under a bad
  name. When the members do not describe one thing, they and the class's scope
  have to change; a better noun on the same grab bag is not really much better.
- A name must not imply more than the object really does. A "branch directory"
  must know every customer's branch, not only some customers'. Rather than
  narrowing the name, try making the object more powerful first.
- Keep every capability the old classes had, or delete one on purpose, if it
  isn't used, and say so in the commit message.
- Where the old intent is unclear, keep the behaviour and add a `TODO:` asking
  the question.

## Check a design

For each public method that does real work, act out, in your notes, a scene of
three to six lines in which objects ask each other for things: "The scheduler
asks the directory which branch customer 7 uses. The directory says the north
branch. The scheduler asks the north branch's calendar for a free slot." If a
line has to say "the helper does the X part", or an object answers a question
it could not plausibly know, the class has no identity yet: rename it, move
methods, or redraw the parts. A plain record with no operations needs no scene;
check instead that its fields describe one thing.
