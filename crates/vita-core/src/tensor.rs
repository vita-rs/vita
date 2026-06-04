//! Fixed-size geometric primitives in three dimensions.
//!
//! [`Vector3`] is a free vector, [`Point3`] an affine point, and [`Matrix3`] a linear map;
//! each is generic over its component type. Point and vector are distinct: subtracting two
//! points yields the [`Vector3`] between them, keeping positions and displacements from
//! being confused.

mod matrix;
mod point;
mod vector;

pub use matrix::Matrix3;
pub use point::Point3;
pub use vector::Vector3;
