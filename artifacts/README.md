# Artifacts

The code and write-ups in this tree were written by coding agents, and no person has yet read them
and cleared them for reading. They build and run as part of vg, and its tests cover them, but their
comments and documents may not explain them well. Read the code in `src/` to learn how vg works.

The tree mirrors `src/`:

- `artifacts/` holds library code. It goes into `libvg`, beside the code in `src/`.
- `artifacts/unittest/` holds unit tests, and `artifacts/unittest/support/` their helpers. They go
  into the `vg` binary and into a test binary of their own in `bin/unittest/`, beside the tests in
  `src/unittest/`.
- `artifacts/doc/` holds write-ups, listed in `artifacts/doc/README.md`.
- `artifacts/scripts/` holds scripts, such as `mosaic_to_path.py`, which the call tests use to check
  a mosaic against the graph.
- `artifacts/test/` holds command-line tests that `make test` does not run, described in
  `artifacts/test/README.md`.

Code in `src/` may include headers from here. No file here other than a `README.md` may share a
name with a file in `src/` or `doc/`, since objects, test binaries and includes are found by name.

A file leaves this tree for `src/` or `doc/` once a person has read it and cleared it for reading.
