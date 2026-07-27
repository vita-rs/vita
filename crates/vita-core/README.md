# Vita Core (`vita-core`)

**The questions a physical system can answer, not the structures that store the answers.**

## Overview

`vita-core` is a zero-dependency vocabulary for atomistic and molecular computation.

| layer              | what it provides                                                                        |
| :----------------- | :-------------------------------------------------------------------------------------- |
| **primitives**     | `Element`, `Isotope`, `SiteId`, … — what is at each site                                |
| **linear algebra** | `Scalar`, `Point3`, `Vector3`, `Matrix3` — abstract 3‑D math                            |
| **capability**     | `HasSites`, `HasPositions`, `HasElements`, … — questions a system can answer            |
| **dimensions**     | 23 physical dimension modules, each with a `Quantity<V, U>` newtype — no unit confusion |

A system is never a concrete type in this crate. Everything is a trait bound — code declares the capabilities it needs and stays blind to the storage behind them:

```rust
use vita_core::prelude::*;

fn center_of_mass<M>(m: &M) -> Point3<Length<f64, Angstrom>>
where
    M: HasSites + HasPositions<f64> + HasMasses<f64>,
{
    // ...
}
```

---

## Installation

```toml
[dependencies]
vita-core = "0.3.0"
```

---

## Documentation

The documentation for `vita-core` is available at [docs.rs/vita-core](https://docs.rs/vita-core).
