---
name: test-design
description: Decide what each unit test checks, on what input, and against what answer. Test a unit's own contract on inputs the test sets up, test each implementation of an interface against the interface, record bugs in other code as expected-to-fail tests, check each fact once, and name tests by property. Use when writing, restructuring or reviewing unit tests.
---

# Test what the unit promises

A unit is the function or class a test is about. A test is evidence that the unit keeps its
contract, so the reader of a failing test needs to know which unit broke and which promise it
broke. Each check therefore tests one unit, on an input the test controls, against a stated
property, and a failure names the unit and the property.

## What to check

1. **Test a unit on inputs you set up.** Give the unit data, or an implementation of the
   interface it reads, that the test builds, and check its results against an answer worked out
   by hand. Bad: test a report generator by checking that it gives the same report from the SQL
   store and from the in-memory store; the stores are not the unit under test, and the check
   passes when the generator is wrong in the same way on both. Good: load three known orders into
   the in-memory store and check the report's totals.
2. **Test each implementation against its interface.** Write the interface's laws as English
   sentences ("a value that was put can be got back", "size() counts each key once"), check them
   in one helper, and run that helper on every implementation, with an `INFO` naming the
   implementation, so a failure names it.
3. **Compare implementations only on what their shared input decides,** and only on inputs too
   big or random to work out by hand. This is a test of the implementations themselves; a unit
   that reads them is still tested as rule 1 says. Two indexes built from the same random list of words must
   list the same words: no law says which words, but the input does. Compare those facts and
   nothing else, and test everything else with the laws. A failed comparison does not say which
   implementation broke, so give each one law tests too.
4. **Check each fact once in a test case.** Before adding a check, ask whether a stronger check
   in the same `TEST_CASE` implies it. Separate test cases may check the same fact. Bad: `REQUIRE(v == vector{1, 2, 3})` then `REQUIRE(v.size() == 3)`. Good: the first
   check alone.

## How to organize tests

1. **Handle any number of implementations.** Keep the implementations under test in a list
   (factory functions, or Catch2's `TEMPLATE_TEST_CASE`) and loop over it, so adding one is one
   more entry, not a new member or a copied test.
2. **Record a bug in other code as a test of that code.** When code that is not the unit under
   test (a library, another module) breaks its interface, write a test of that code on the
   smallest input that shows the bug, named for the rule it breaks, and tag it `[!mayfail]`.
   Catch2 runs such a test and reports its failure without failing the run. Put it with that
   code's own tests when the task lets you edit that code; otherwise put it beside the tests of the
   code that depends on it. Where
   another test has to avoid the bug, add a `TODO:` naming the bug test. Catch2 does not report
   a `[!mayfail]` test that starts passing, so after changing the buggy code, run its tagged tests;
   for each that passes, remove the tag and the `TODO:`s.
3. **Name tests by property.** Each `TEST_CASE` or `SECTION` name states something that ought to
   be true ("An empty cart costs nothing"), not a topic ("Empty carts") or a category ("Edge
   cases"). Use the names of the calls and constants under test where they exist:
   "`total()` of an empty cart is 0" says more than "An empty cart costs nothing". A `TEST_CASE`
   whose sections state different properties names the unit and situation they share ("A cart
   with one item"). Leave commas out of names, because Catch2 splits a name given on the command
   line at its commas.
4. **Use the interface under test as callers do.** Call its own methods and callbacks. A helper
   that turns results into strings or vectors to compare them needs a comment saying exactly what
   it compares and what it ignores.
