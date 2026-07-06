# Publishing

How releases work for `stages-thermo`: PyPI + crates.io publishes happen
automatically on `v*` tag pushes via `.github/workflows/release.yml`. End users
get the library via `pip install stages-thermo` (`import stages`) or
`cargo add stages-thermo`.

`stages-thermo` publishes **one** crate (unlike vle's two: `vle-units` +
`vle-thermo`). It reuses `vle-units` transitively through its `vle-thermo`
dependency, so there is nothing extra to publish here.

> **Milestone rule:** each milestone that adds Rust functionality also adds its
> PyO3 bindings in the same commit series (see `CLAUDE.md`) and typically cuts a
> new release tag. Tag → push → registries is the standard milestone-completion
> workflow, matching vle.

---

## How a release flows

```
git tag vX.Y.Z && git push origin vX.Y.Z
   │
   └── release.yml fires:
        ├── check-already-published → probe PyPI / crates.io / GH Releases
        │     (idempotent: a re-run skips already-published targets)
        ├── build-wheels  → cibuildwheel matrix, GitHub-hosted:
        │     • Linux x86_64      (ubuntu-latest)
        │     • Linux aarch64     (ubuntu-24.04-arm)
        │     • macOS arm64       (macos-14)
        │     • Windows AMD64     (windows-latest)
        │     All abi3-tagged (cp310-abi3-*): one wheel per (OS, arch) covers
        │     CPython 3.10+.
        ├── build-sdist   → maturin sdist
        ├── publish-pypi   → PyPI Trusted Publishing (OIDC, no token)
        ├── publish-crates → cargo publish stages-thermo (token from 1Password)
        └── gh-release     → GitHub Release with all wheels + sdist
```

---

## Cutting a release

1. **Bump the version in both places — they must match.**
   - `Cargo.toml` (workspace root) → `[workspace.package] version = "X.Y.Z"`
   - `python/pyproject.toml` → `[project] version = "X.Y.Z"`

2. **Update the docs** per `CLAUDE.md` release rules: `ROADMAP.md`, `TODO.md`,
   `README.md`, and both package READMEs (`engine/README.md`,
   `python/README.md`) — remembering the immutable-per-published-version rule
   for the two package READMEs.

3. **Run the pre-push fmt gate** (`hooks/pre-push`, wired via
   `git config core.hooksPath hooks`).

4. **Commit and push the version bump to `main` first** (YubiKey-signed).

5. **Tag and push the tag.**
   ```sh
   git tag -a vX.Y.Z -m "Release X.Y.Z"
   git push origin vX.Y.Z
   ```

6. **Watch the workflow.** In the Actions tab, approve the `crates-io`
   environment (and `pypi` if it has required reviewers) to release the publish
   steps.

7. **Verify** (below).

---

## Name-holding stubs (Milestone 0)

The `0.0.1` stubs hold the `stages-thermo` name on **both** registries (PLAN §2,
§11). **Prerequisite:** the credentials + pending publisher in "Credentials &
1Password vault" below must exist first. Then publish via **Route A** (the
tagged-release pipeline — token-free for PyPI, and the cleanest way to convert
the PyPI pending publisher) or **Route B** (manual crates.io only). Both are
spelled out at the end of this document.

crates.io + PyPI versions are **immutable** — `0.0.1` is spent once uploaded.
Double-check the version before publishing either.

---

## Post-publish verification

```sh
# crates.io
cargo search stages-thermo
cargo new --lib scratch && cd scratch && cargo add stages-thermo && cargo build

# PyPI — fresh venv
python -m venv /tmp/stages-check && source /tmp/stages-check/bin/activate
pip install stages-thermo==X.Y.Z
python -c "import stages; print(stages.__version__); print(stages.smoke_bubble_temperature())"
deactivate && rm -rf /tmp/stages-check
```

---

## Rolling back

**crates.io:** `cargo yank --version X.Y.Z stages-thermo` (hides the version
from new resolvers; existing lockfiles still resolve it). You **cannot** delete
a published version — yank + publish a patch.

**PyPI:** `twine yank stages-thermo==X.Y.Z --reason "..."` (or a PyPI removal
request). Same "cannot truly delete" rule — yank + patch.

---

## Credentials & 1Password vault (one-time, do this first)

> **Do NOT reuse the `vle-thermo-ci` vault or its `crates-io` token.** That
> crates.io token is **scoped to the vle crates** (`vle-thermo`, `vle-units`)
> and returns **HTTP 403** on a new crate name like `stages-thermo`. This
> project gets its **own** vault and its **own** credentials.

**PyPI needs no token at all** — Miguel publishes via **Trusted Publishing
(OIDC)**, the same as vle. crates.io is still token-based, so that one token is
the only registry secret we store.

### 1. New 1Password vault + crates.io token

1. Create a 1Password vault **`stages-thermo-ci`** (mirrors `vle-thermo-ci`; a
   separate vault keeps each project's publish token isolated).
2. Create a crates.io token at <https://crates.io/settings/tokens>:
   - Scopes: **`publish-new`** + **`publish-update`**.
   - Crate scope: restrict to **`stages-thermo`** (or leave unscoped). It must
     allow **publish-new** or the very first `0.0.1` upload 403s.
3. Store it in the new vault as item **`crates-io`**, field **`token`** →
   reference `op://stages-thermo-ci/crates-io/token` (this is exactly what
   `release.yml` reads).
4. For CI, create/attach a **read-only service account** on the
   `stages-thermo-ci` vault and put its token in the GitHub repo secret
   **`OP_SERVICE_ACCOUNT_TOKEN`** — the **only** GitHub-side secret. (For a
   local manual publish, `op read op://stages-thermo-ci/crates-io/token` is
   enough; no service account required.)

### 2. PyPI Trusted Publisher (no token — OIDC)

Because the `stages-thermo` project **does not exist on PyPI yet**, register a
**pending** publisher (under your *account*, not a project, since there's no
project to attach to). At
<https://pypi.org/manage/account/publishing/> → "Add a pending publisher":

- **PyPI Project Name:** `stages-thermo`
- **Owner:** `miguelju`
- **Repository name:** `stages-thermo`
- **Workflow name:** `release.yml`
- **Environment name:** `pypi`

On the first successful OIDC publish the pending publisher **auto-converts to a
normal publisher** and creates the project — no token, ever. ⚠️ **Caveat:** if
anyone else registers the name `stages-thermo` on PyPI before your first
publish, the pending publisher is invalidated — so register it and cut the
`v0.0.1` release promptly. (Ref: PyPI docs, "Creating a PyPI Project with a
Trusted Publisher".)

### 3. GitHub Environments (in the repo settings)

- **`pypi`** — for the OIDC publish (matches the pending-publisher config above;
  add required reviewers if desired).
- **`crates-io`** — add Miguel as a required reviewer so every crate publish
  pauses for a human click (crates.io versions are immutable).

---

## Publishing the 0.0.1 name-holding stub — the two routes

Once the vault + token (crates.io) and pending publisher (PyPI) exist:

**Route A — the real pipeline (recommended, token-free for PyPI):** create the
GitHub repo, push, set up the two environments, then `git tag v0.0.1 && git push
origin v0.0.1`. `release.yml` builds the wheels and publishes to **both**
registries (PyPI via OIDC, crates.io via the vault token). This is also the
end-to-end CI smoke test.

**Route B — manual (crates.io only; faster for just holding the name):**
```sh
# crates.io — token from the new vault, never printed:
CARGO_REGISTRY_TOKEN="$(op read op://stages-thermo-ci/crates-io/token)" \
  cargo publish -p stages-thermo --allow-dirty
```
PyPI has no manual token route in this setup — use Route A (a tagged release, or
a `workflow_dispatch` run with `dry_run=false`) so the pending publisher
converts. If you ever need a truly manual PyPI upload, mint a one-off
account-scoped token, but that's off the standard path.
