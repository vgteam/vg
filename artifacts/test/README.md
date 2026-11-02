# Command-line tests not run by `make test`

`vg_call_read_likelihood.t` holds the command-line tests of `vg call --read-likelihood` that were
cut from `test/t/18_vg_call.t`. The pull request that added read-likelihood genotyping wrote more
than 300 of them, which would have tripled the size of that file. To keep the tests that CI runs
manageable, `test/t/18_vg_call.t` keeps about 30 that are most likely to catch a regression, and
this file keeps the rest, runnable:

```bash
cd test
prove -v ../artifacts/test/vg_call_read_likelihood.t
```

It builds its own fixtures from the test data in `test/`, and uses `bin/vg`.

Each test here pins a case that one change to the caller once got wrong: an error message, an
option that must be refused, or a property of the mosaic or anchor file on a hand-built graph. Some
of them check library behaviour that a unit test in `artifacts/unittest/` could check more cheaply;
those that were easy to move were moved there.

In the future, we should explore ways to manage command-line tests at least an order of magnitude
more numerous and complex than the flat files of counted TAP tests in `test/t/` allow, which can
only be run all together. Then tests like these could run in CI again.
