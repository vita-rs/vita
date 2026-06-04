# Contributing

## Prerequisites

Rust 1.85 (edition 2024) or later, with `rustfmt` and `clippy`:

```sh
rustup component add rustfmt clippy
```

## Development

All crates live under `crates/`. Before every commit:

```sh
cargo build --workspace
cargo test --workspace
cargo clippy --workspace -- -D warnings
cargo fmt --check
```

CI enforces zero clippy warnings and correct formatting.

## Commit Convention

Vita uses [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/).

```text
<type>[!](<scope>): <description>

[optional body]

[optional footer]
```

### Types

| Type       | When to use                                     |
| ---------- | ----------------------------------------------- |
| `feat`     | New public API or capability                    |
| `fix`      | Bug fix                                         |
| `perf`     | Performance improvement with no behavior change |
| `refactor` | Internal restructure with no behavior change    |
| `test`     | Add or fix tests                                |
| `docs`     | Documentation only                              |
| `chore`    | Dependency updates, config, tooling             |
| `ci`       | CI/CD pipeline changes                          |

### Scopes

Use the crate's short name (without the `vita-` prefix): `core` · `vita` · `workspace`

Additional scopes will be added as new crates are introduced.

### Breaking changes

Append `!` to the type, or add a `BREAKING CHANGE:` footer:

```text
feat!(core): Rename atom identifier type

BREAKING CHANGE: The atom ID type has been renamed. Update all call sites.
```

### Examples

```text
feat(core): Add charge query trait
fix(core): Handle malformed atom records
perf(core): Reduce allocation in atom store
chore(workspace): Bump serde to 1.0.220
ci: Add miri job
```

**Rules:**

- Description starts with a capital letter, imperative mood, no trailing period
- Body and footer are separated from the description by a blank line
- Scope is required for all types except `ci` and `workspace`-level `chore`

## Code Rules

### Safety

Every `unsafe` block must carry a `// SAFETY:` comment explaining why the invariants required by the unsafe operation are upheld:

```rust
// SAFETY: idx < u32::MAX is asserted above, so adding 1 cannot overflow,
// and the result is non-zero.
let id = unsafe { NonZeroU32::new_unchecked(idx as u32 + 1) };
```

PRs that add `unsafe` without a `SAFETY:` comment will not be merged.

<!-- TODO: Add more rules here. -->

## Pull Requests

Before opening a PR, make sure your follow the contributing guidelines.

Vita uses **merge commits**. Each merged PR produces one merge commit in `main`, marking the PR boundary. The PR title becomes the merge commit message and must follow the commit convention. Individual commits within the PR must also follow the commit convention and represent exactly one logical change. Maintainers will ask you to rebase and clean up before merge if commits are not atomic.

If you find an unrelated issue while working, fix it in a separate PR.

## Issues

<!--TODO: Add issue guidelines here. -->

## Accountability

AI-assisted development is welcome. Vibe coding is not.

You are responsible for every line you submit — its correctness, its safety properties, and its alignment with vita's design. Maintainers may ask you to explain any part of a PR in detail; inability to do so is grounds for rejection.

AI slop will not be accepted.

## Code of Conduct

Vita follows the [Contributor Covenant 3.0](CODE_OF_CONDUCT.md).
