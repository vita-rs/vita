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
