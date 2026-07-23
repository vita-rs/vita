# Vita Chem (`vita-chem`)

**What a declared molecular structure determines, never what a heuristic would guess.**

## Overview

`vita-chem` is the cheminformatics crate of the vita ecosystem: discrete chemistry derived exactly from what a molecule declares — identity is computed, values are never invented.

| layer          | what it provides                                                                                              |
| :------------- | :------------------------------------------------------------------------------------------------------------ |
| **primitives** | `BondId`, `BondOrder`, `StereoConfiguration`, … — the bond vocabulary                                         |
| **capability** | `HasBondOrders`, `HasFormalCharges`, `HasStereoConfigurations`, … — the chemical questions                    |
| **graph**      | `topology`, `canonical`, `isomorphism` — skeleton facts, portable identity, substructure search               |
| **perception** | `valence`, `aromaticity`, `conjugation`, `hybridization`, `stereo` — electron bookkeeping and declared models |
| **projection** | `composition`, `fingerprint`, `notation` — particle inventory, similarity, standard symbols                   |

A molecule is never a concrete type here either. Every algorithm bounds on exactly the capabilities it consumes and reads the molecule unchanged:

```rust
use vita_chem::{HasFormalCharges, composition::elemental, notation::formula};
use vita_core::HasElements;

fn hill<M: HasElements + HasFormalCharges>(m: &M) -> String {
    formula::write(&elemental(m))
}
```

---

## Installation

```toml
[dependencies]
vita-chem = "0.0.0"
```

---

## Documentation

The documentation for `vita-chem` is available at [docs.rs/vita-chem](https://docs.rs/vita-chem).
