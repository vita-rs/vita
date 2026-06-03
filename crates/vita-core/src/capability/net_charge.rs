use crate::Scalar;
use crate::units::charge::{Charge, ChargeUnit};

/// The net electric charge of a system.
///
/// [`net_charge`](HasNetCharge::net_charge) returns the system's total charge in any
/// requested [unit](ChargeUnit), positive for a net electron deficit.
pub trait HasNetCharge<V: Scalar> {
    /// Returns the system's net charge, in unit `U`.
    fn net_charge<U: ChargeUnit>(&self) -> Charge<V, U>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::charge::ElementaryCharge;

    struct Ion {
        charge: Charge<f64, ElementaryCharge>,
    }
    impl HasNetCharge<f64> for Ion {
        fn net_charge<U: ChargeUnit>(&self) -> Charge<f64, U> {
            self.charge.to()
        }
    }

    #[test]
    fn net_charge() {
        let ion = Ion {
            charge: Charge::new(-1.0),
        };
        assert_eq!(ion.net_charge::<ElementaryCharge>(), Charge::new(-1.0));
    }
}
