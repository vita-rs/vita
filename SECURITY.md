# Security Policy

## Scope

Vita is a Rust framework. A security vulnerability in this context is any defect in **unsafe code** that can be triggered by **safe Rust** calling code, resulting in undefined behavior (UB) — including memory corruption, use-after-free, invalid memory access, or data races.

The following are **not** security vulnerabilities and should be reported as regular issues:

- Incorrect scientific results or numerical inaccuracies
- Unexpected panics in documented error conditions
- API design concerns or ergonomic issues
- Performance regressions

## Supported Versions

Vita has not yet reached 1.0. Only the latest release receives security fixes. No patches are backported to older versions.

## Reporting a Vulnerability

**Do not open a public GitHub issue for security vulnerabilities.**

**Preferred:** [Open a private security advisory](https://github.com/vita-rs/vita/security/advisories/new) on GitHub. This keeps the report confidential until a fix is released.

**Fallback:** Send an email to <tianchengkan@gmail.com> with a description of the issue.

Please include:

- A description of the vulnerability and its potential impact
- A minimal reproducible example or steps to trigger the issue
- Affected crate(s) and version(s), if known

## Response

This project is maintained on a best-effort basis. We aim to:

- Acknowledge receipt within **7 days**
- Provide an initial assessment within **14 days**
- Release a fix as soon as practical after confirmation

We practice coordinated disclosure: please do not publicly disclose the vulnerability until a fix has been released.

## Credit

Reporters are credited by name or pseudonym in the [GitHub Security Advisory](https://github.com/vita-rs/vita/security/advisories) published upon fixing the vulnerability. If anonymity is preferred, let us know and no identifying information will be included.
