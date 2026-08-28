# Personal Copilot CLI instructions (@atlejg)

Always-on notes loaded in every Copilot CLI session, across all repositories.

## Interaction conventions

- **Prompts starting with `q&a`** are **discussion only**: answer / explain and
  do **not** take any action or make changes (no edits, no mutating commands).

## Gotchas

### IPython `%pylab` shadows Python builtins with numpy

Interactive scripts run via `ipython -i script.py` (e.g. the AWT sandbox scripts
under `equinor/apo-r` at `sandbox/Atle/Automated_Well_Testing/`) often run in a
`%pylab` / `from numpy import *` namespace. This **replaces the builtins**
`any`, `all`, `sum`, `min`, `max` with their numpy versions, which behave
differently and cause silent logic bugs or crashes:

- `max(a, b)` → `numpy.max(a, axis=b)` → `AxisError`. Same for two-arg `min`.
- `any(<generator>)` → `numpy.any(<generator>)` returns **`True`
  unconditionally** — numpy wraps the generator as a truthy 0-d object array and
  never iterates it. `all(<generator>)` is likewise broken.

**Rules when writing/reviewing code that may run under `%pylab`:**
- Never call `any(...)` / `all(...)` on a **generator**. Use list truthiness
  instead, e.g. `failed = [x for x in items if cond(x)]; return not failed`.
- Never call `max` / `min` / `sum` with the two-arg / axis-style signature on
  scalars. Use explicit comparisons (`a if a > b else b`) or `math.fsum`.
- Symptom to watch for: diagnostics that use list-truthiness disagree with
  detection logic that uses `any()`/`all()` (e.g. "long stable run reported, but
  0 results detected").

## Technical considerations

- use 'uv run python' for running python
- i ususally run python scripts from an ipython --pylab session, as indicated above

