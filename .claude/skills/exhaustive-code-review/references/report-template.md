# Review report template

Keep it short. Leave out any section with nothing in it.

```
## Review of <branch or PR> (<base>..<head>, N files)

**Verdict:** ready / ready after the blocking fixes / needs another round

**Ran:** <build with which warnings, which tests, which sanitizers or tools>
**Not run:** <tools or tests skipped, and why>
**Lenses:** <lenses run>; skipped <lens> because <reason>

### Blocking
1. `src/file.cpp:123` — <the problem, as the input that triggers it and the
   wrong result>. <Why it matters.> Fix: <concrete change>. (Regression /
   pre-existing, when both versions were run; "probable regression, traced"
   when settled only by reading.)

### Important
1. ...

### Minor
- `file:line` — <one line each>

### Checked and fine
<One short paragraph: what was looked at closely and found correct, so the
author knows it was covered.>
```

State confidence on any finding you did not fully trace ("didn't check all
callers"). Phrase judgment calls as questions.
