# Vita Core (`vita-core`)

**The questions a physical system can answer, not the structures that store the answers.**

## Overview

`vita-core` is a zero-dependency vocabulary for atomistic and molecular computation. Every concept it exposes falls into one of three categories — the identity of a thing, the geometry of where it sits, or a question you can ask about it — and nothing else.

| layer          | what it provides                                                                        |
| :------------- | :-------------------------------------------------------------------------------------- |
| **primitives** | `Element`, `Isotope`, `SiteId` — what is at each site                                   |
| **geometry**   | `Point3`, `Vector3`, `Matrix3`, `Lattice` — where each site sits, in what cell          |
| **capability** | `HasSites`, `HasPositions`, `HasElements`, … — questions a system can answer            |
| **dimensions** | 20 physical dimension modules, each with a `Quantity<V, U>` newtype — no unit confusion |

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
vita-core = "0.0.0"
```

---

## Documentation

The documentation for `vita-core` is available at [docs.rs/vita-core](https://docs.rs/vita-core).
