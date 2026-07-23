# Vita IO (`vita-io`)

**I/O for the vita ecosystem: maps byte streams to vita capabilities and back.**

## Overview

`vita-io` supports the following formats:

| format | module         | read output   | write bounds                    |
| :----- | :------------- | :------------ | :------------------------------ |
| XYZ    | `vita_io::xyz` | `xyz::System` | `HasElements + HasPositions<V>` |

Reading produces a concrete type implementing exactly the capabilities the format records. Writing accepts any type bounded on those same capabilities, whatever its origin.

```rust
use vita_io::xyz;
use vita_core::units::length::Angstrom;

// any `BufRead` -> `xyz::System<f64>`
let system = xyz::read(reader).system::<f64>()?;

// any `HasElements + HasPositions<f64>` -> any `Write`
xyz::write::<_, Angstrom, _>(&mut writer, &system, &xyz::Config { comment: "" })?;
```

---

## Installation

```toml
[dependencies]
vita-core = "0.3.0"
vita-io = "0.2.0"
```

---

## Documentation

The documentation for `vita-io` is available at [docs.rs/vita-io](https://docs.rs/vita-io).
