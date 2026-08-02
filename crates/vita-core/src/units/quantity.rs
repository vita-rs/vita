/// Generate a physical-quantity newtype parameterized by a scalar type `V`
/// and a unit marker `U`, implementing [`Quantity`](crate::Quantity).
macro_rules! define_quantity {
    ($(#[$meta:meta])* $Qty:ident, $UnitTrait:ident) => {

        $(#[$meta])*
        #[repr(transparent)]
        pub struct $Qty<V, U>(V, ::core::marker::PhantomData<U>);

        impl<V: ::core::marker::Copy, U> ::core::marker::Copy for $Qty<V, U> {}

        impl<V: ::core::marker::Copy, U> ::core::clone::Clone for $Qty<V, U> {
            #[inline]
            fn clone(&self) -> Self { *self }
        }

        impl<V: ::core::fmt::Debug, U> ::core::fmt::Debug for $Qty<V, U> {
            fn fmt(&self, f: &mut ::core::fmt::Formatter<'_>) -> ::core::fmt::Result {
                f.debug_tuple(stringify!($Qty)).field(&self.0).finish()
            }
        }

        impl<V: ::core::cmp::PartialEq, U> ::core::cmp::PartialEq for $Qty<V, U> {
            #[inline]
            fn eq(&self, other: &Self) -> bool { self.0 == other.0 }
        }

        impl<V: ::core::cmp::PartialOrd, U> ::core::cmp::PartialOrd for $Qty<V, U> {
            #[inline]
            fn partial_cmp(
                &self,
                other: &Self,
            ) -> ::core::option::Option<::core::cmp::Ordering> {
                self.0.partial_cmp(&other.0)
            }
        }

        impl<V, U> $Qty<V, U> {
            /// Constructs a quantity holding `value` in unit `U`.
            #[inline]
            pub const fn new(value: V) -> Self {
                Self(value, ::core::marker::PhantomData)
            }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> $Qty<V, U> {
            /// Returns `self` expressed in unit `T`.
            ///
            /// Conversion is performed in `f64` precision via the `TO_CANONICAL`
            /// factors on `U` and `T`; the result is then narrowed back to `V`.
            #[inline]
            pub fn to<T: $UnitTrait>(self) -> $Qty<V, T> {
                let ratio = const { U::TO_CANONICAL / T::TO_CANONICAL };
                $Qty::new(V::from_f64(self.0.to_f64() * ratio))
            }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> $crate::Quantity for $Qty<V, U> {
            type Value = V;

            const ZERO: Self = Self::new(<V as $crate::Quantity>::ZERO);
            const INFINITY: Self = Self::new(<V as $crate::Quantity>::INFINITY);
            const NEG_INFINITY: Self = Self::new(<V as $crate::Quantity>::NEG_INFINITY);
            const NAN: Self = Self::new(<V as $crate::Quantity>::NAN);

            #[inline]
            fn value(self) -> V { self.0 }

            #[inline]
            fn from_value(value: V) -> Self { Self::new(value) }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::default::Default
            for $Qty<V, U>
        {
            /// Returns a zero quantity in unit `U`.
            #[inline]
            fn default() -> Self { Self::new(V::default()) }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::ops::Neg for $Qty<V, U> {
            type Output = Self;
            #[inline]
            fn neg(self) -> Self { Self::new(-self.0) }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::ops::Add for $Qty<V, U> {
            type Output = Self;
            #[inline]
            fn add(self, rhs: Self) -> Self { Self::new(self.0 + rhs.0) }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::ops::AddAssign
            for $Qty<V, U>
        {
            #[inline]
            fn add_assign(&mut self, rhs: Self) { self.0 += rhs.0; }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::ops::Sub for $Qty<V, U> {
            type Output = Self;
            #[inline]
            fn sub(self, rhs: Self) -> Self { Self::new(self.0 - rhs.0) }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::ops::SubAssign
            for $Qty<V, U>
        {
            #[inline]
            fn sub_assign(&mut self, rhs: Self) { self.0 -= rhs.0; }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::ops::Rem for $Qty<V, U> {
            type Output = Self;
            #[inline]
            fn rem(self, rhs: Self) -> Self { Self::new(self.0 % rhs.0) }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::ops::RemAssign
            for $Qty<V, U>
        {
            #[inline]
            fn rem_assign(&mut self, rhs: Self) { self.0 %= rhs.0; }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::ops::Mul<V>
            for $Qty<V, U>
        {
            type Output = Self;
            #[inline]
            fn mul(self, rhs: V) -> Self { Self::new(self.0 * rhs) }
        }

        impl<U: $UnitTrait> ::core::ops::Mul<$Qty<f32, U>> for f32 {
            type Output = $Qty<f32, U>;
            #[inline]
            fn mul(self, rhs: $Qty<f32, U>) -> $Qty<f32, U> { rhs * self }
        }

        impl<U: $UnitTrait> ::core::ops::Mul<$Qty<f64, U>> for f64 {
            type Output = $Qty<f64, U>;
            #[inline]
            fn mul(self, rhs: $Qty<f64, U>) -> $Qty<f64, U> { rhs * self }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::ops::MulAssign<V>
            for $Qty<V, U>
        {
            #[inline]
            fn mul_assign(&mut self, rhs: V) { self.0 *= rhs; }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::ops::Div<V>
            for $Qty<V, U>
        {
            type Output = Self;
            #[inline]
            fn div(self, rhs: V) -> Self { Self::new(self.0 / rhs) }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::ops::DivAssign<V>
            for $Qty<V, U>
        {
            #[inline]
            fn div_assign(&mut self, rhs: V) { self.0 /= rhs; }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::ops::Rem<V>
            for $Qty<V, U>
        {
            type Output = Self;
            #[inline]
            fn rem(self, rhs: V) -> Self { Self::new(self.0 % rhs) }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::ops::RemAssign<V>
            for $Qty<V, U>
        {
            #[inline]
            fn rem_assign(&mut self, rhs: V) { self.0 %= rhs; }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::ops::Div for $Qty<V, U> {
            type Output = V;
            /// Returns the dimensionless ratio `self / rhs`.
            #[inline]
            fn div(self, rhs: Self) -> V { self.0 / rhs.0 }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::iter::Sum for $Qty<V, U> {
            #[inline]
            fn sum<I: ::core::iter::Iterator<Item = Self>>(iter: I) -> Self {
                iter.fold(Self::new(V::default()), |acc, x| acc + x)
            }
        }

        impl<'__qty, V: $crate::Scalar, U: $UnitTrait>
            ::core::iter::Sum<&'__qty $Qty<V, U>> for $Qty<V, U>
        {
            #[inline]
            fn sum<I: ::core::iter::Iterator<Item = &'__qty Self>>(iter: I) -> Self {
                iter.copied().fold(Self::new(V::default()), |acc, x| acc + x)
            }
        }

        impl<V: $crate::Scalar, U: $UnitTrait> ::core::fmt::Display
            for $Qty<V, U>
        {
            fn fmt(&self, f: &mut ::core::fmt::Formatter<'_>) -> ::core::fmt::Result {
                write!(f, "{} {}", self.0, U::SYMBOL)
            }
        }
    };
}

pub(crate) use define_quantity;
