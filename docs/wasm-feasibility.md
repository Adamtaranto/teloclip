# WebAssembly / Pyodide feasibility

Investigation of the "Add wasm wheel" item in `.plans/02_code_review.md`.

## Conclusion

**No separate wasm wheel is needed, and none can be built.** teloclip is pure
Python, so the wheel it already publishes is `py3-none-any` — that *is* the
artefact a browser would install. There is no architecture-specific variant to
add.

The real question was whether teloclip's dependencies can be satisfied under
Pyodide. They can, contrary to the assumption that pysam was a blocker.

## Verified findings

Built wheel tag, from `pip wheel . --no-deps`:

```
teloclip-0.3.5.dev33+ga91feb0d9-py3-none-any.whl
```

Dependency availability, checked against the Pyodide 0.28.0 lockfile
(`pyodide-lock.json`, the authoritative list of what Pyodide ships):

| Dependency | Wheel type | In Pyodide | Notes |
|---|---|---|---|
| `click` | pure | yes (8.2.1) | |
| `rich` | pure | yes (13.9.4) | |
| `biopython` | native | yes (1.85) | Built by the Pyodide project |
| `pysam` | native | **yes (0.23.0)** | Built by the Pyodide project |
| `pyfaidx` | pure | no | Not shipped, but installs from PyPI via `micropip` because it is pure Python |

So every dependency is reachable. `pysam` being available is the load-bearing
finding: the plan assumed it was not, which would have restricted any wasm build
to a `filter`-only subset.

## Sub-command dependency footprint

Established by simulating an environment where `import pysam` fails and running
each command against `tests/data`:

| Command | Runs without pysam | Notes |
|---|---|---|
| `filter` | yes | Pure standard library, click and rich. Kept 37 records, as with pysam present |
| `extract` | yes | Uses biopython for FASTA writing; the two `import pysam` statements in `extract_io.py` are function-local and are not on the extract path |
| `extend` | no | Imports pysam and pyfaidx at module scope |

This matters beyond wasm: it means `filter` and `extract` are usable wherever
pysam cannot be built.

## Bug found and fixed

`cli.register_commands` imported all three sub-commands inside a single `try`
block. Any one failing meant `add_command` was never reached for the others, so
in a pysam-free environment **no** commands registered at all:

```
Warning: Could not import commands: No module named 'pysam'
registered subcommands: []
```

Each is now imported independently, and only the genuinely unavailable one is
skipped, with a message naming it:

```
Warning: sub-command from teloclip.commands.extend is unavailable: No module named 'pysam'
registered subcommands: ['extract', 'filter']
```

Covered by `test_one_unimportable_command_does_not_hide_the_others`.

## Not verified

I could not run teloclip inside an actual Pyodide runtime: no Node.js is
available in this environment, and Pyodide cannot be exercised from CPython. The
dependency analysis above is from the Pyodide lockfile and from simulating
pysam's absence under CPython, not from a browser run.

Anyone picking this up should confirm with a real Pyodide session before
advertising browser support:

```js
const pyodide = await loadPyodide();
await pyodide.loadPackage(["micropip", "pysam", "biopython"]);
await pyodide.runPythonAsync(`
  import micropip
  await micropip.install(["pyfaidx", "teloclip"])
  from teloclip.cli import main
`);
```

Two things to watch, neither of which the static analysis can settle:

- **File system.** Pyodide's virtual FS has no real files. `extend` requires an
  indexed BAM and FASTA on disk, so a browser workflow would need files mounted
  into the Emscripten FS first. `filter` reads a stream, so it is the natural
  entry point.
- **pysam's htslib build.** Pyodide ships pysam, but whether BAM indexing and
  random access work under wasm is a property of that build, not of teloclip.

## Recommendation

Close the "add wasm wheel" item as not-applicable: the published wheel already
serves. If browser support is wanted as a feature, the work is a Pyodide
integration test plus documentation of the `micropip` install, not a build
change. That is worth doing only if there is a concrete use case, since it adds
a runtime nobody currently tests against.
