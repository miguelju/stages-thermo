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

Before any real release, `0.0.1` stubs are registered on **both** registries to
hold the `stages-thermo` name (PLAN §2, §11). This is a one-time manual publish,
not a tagged release:

```sh
# crates.io — needs a login token (crates.io → Account Settings → API Tokens)
cargo publish -p stages-thermo            # from the workspace root
# (add --dry-run first to sanity-check the packaged files)

# PyPI — build the sdist + one wheel and upload with twine (needs a PyPI token)
cd python && maturin sdist --out ../dist
twine upload ../dist/*
```

crates.io + PyPI versions are **immutable** — `0.0.1` is spent once uploaded.
Double-check the version before pushing either.

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

## One-time setup (do once, then document as done)

### PyPI Trusted Publishing
At <https://pypi.org/manage/account/publishing/>, add a Trusted Publisher for
project `stages-thermo` pointing at:
- Owner: `miguelju`
- Repository: `stages-thermo`
- Workflow: `release.yml`
- Environment: `pypi`

No API token needed thereafter. (For the very first `0.0.1` stub, PyPI's
"pending publisher" flow or a one-off token upload is fine, since the project
doesn't exist yet.)

### crates.io token (1Password)
1. Create a token at <https://crates.io/settings/tokens> scoped to
   `publish-new` + `publish-update` for `stages-thermo`.
2. Store it in 1Password at `stages-thermo-ci/crates-io/token`.
3. The release workflow loads it via `1password/load-secrets-action` using the
   single GitHub secret `OP_SERVICE_ACCOUNT_TOKEN` (a read-only service account
   on the `stages-thermo-ci` vault). This is the **only** GitHub-side secret.

### GitHub Environments
- `pypi` — for the OIDC publish (add required reviewers if desired).
- `crates-io` — add Miguel as a required reviewer so every crate publish pauses
  for a human click (crates.io versions are immutable).
