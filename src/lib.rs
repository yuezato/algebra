#[macro_use]
extern crate lazy_static;

pub mod erasure_code;
pub mod field;
pub mod fin_field;
pub mod matrix;
pub mod reed_solomon;
pub mod univariate_polynomial;
pub mod vandermonde;
pub mod vecteur;
pub mod bit_based_gf8;

#[cfg(test)]
extern crate rand;
