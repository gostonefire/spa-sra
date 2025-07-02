const PI: f64 = std::f64::consts::PI;   // in original code: 3.1415926535897932384626433832795028841971
const SUN_RADIUS: f64 = 0.26667;

pub const L_COUNT: usize = 6;
pub const B_COUNT: usize = 2;
pub const R_COUNT: usize = 5;
pub const Y_COUNT: usize = 63;

pub const L_MAX_SUBCOUNT: usize = 64;
pub const B_MAX_SUBCOUNT: usize = 5;
pub const R_MAX_SUBCOUNT: usize = 40;

pub const TERM_A: usize = 0;
pub const TERM_B: usize = 1;
pub const TERM_C: usize = 2;
pub const TERM_COUNT: usize = 3;

pub const TERM_X0: usize = 0;
pub const TERM_X1: usize = 1;
pub const TERM_X2: usize = 2;
pub const TERM_X3: usize = 3;
pub const TERM_X4: usize = 4;
pub const TERM_X_COUNT: usize = 5;

pub const TERM_PSI_A: usize = 0;
pub const TERM_PSI_B: usize = 1;
pub const TERM_EPS_C: usize = 2;
pub const TERM_EPS_D: usize = 3;
pub const TERM_PE_COUNT: usize = 4;

pub const JD_MINUS: usize = 0;
pub const JD_ZERO: usize = 1;
pub const JD_PLUS: usize = 2;
pub const JD_COUNT: usize = 3;

pub const SUN_TRANSIT: usize = 0;
pub const SUN_RISE: usize = 1;
pub const SUN_SET: usize = 2;
pub const SUN_COUNT: usize = 3;

pub const TERM_Y_COUNT: usize = TERM_X_COUNT;


/*

enum {TERM_A, TERM_B, TERM_C, TERM_COUNT};
enum {TERM_X0, TERM_X1, TERM_X2, TERM_X3, TERM_X4, TERM_X_COUNT};
enum {TERM_PSI_A, TERM_PSI_B, TERM_EPS_C, TERM_EPS_D, TERM_PE_COUNT};
enum {JD_MINUS, JD_ZERO, JD_PLUS, JD_COUNT};
enum {SUN_TRANSIT, SUN_RISE, SUN_SET, SUN_COUNT};

#define TERM_Y_COUNT TERM_X_COUNT
 */