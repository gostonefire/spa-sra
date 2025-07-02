use crate::constants::PI;

fn rad2deg(radians: f64) -> f64 {
    180.0/PI * radians
}

fn deg2rad(degrees: f64) -> f64 {
    PI/180.0 * degrees
}

fn integer(value: f64) -> i64 {
    value as i64
}

fn limit_degrees(degrees: f64) -> f64 {
    let deg = degrees / 360.0;
    let limited = 360.0 * (deg - deg.floor());
    
    if limited < 0.0 {
        limited + 360.0
    } else {
        limited
    }
}

fn limit_degrees180pm(degrees: f64) -> f64 {
    let deg = degrees / 360.0;
    let limited = 360.0 * (deg - deg.floor());
    
    if limited < -180.0 {
        limited + 360.0
    } else if limited > 180.0 {
        limited - 360.0
    } else {
        limited
    }
}

fn limit_degrees180(degrees: f64) -> f64 {
    let deg = degrees / 180.0;
    let limited = 180.0 * (deg - deg.floor());
    
    if limited < 0.0 {
        limited + 180.0
    } else {
        limited
    }
}

fn limit_zero2one(value: f64) -> f64 {
    let limited = value - value.floor();
    
    if limited < 0.0 {
        limited + 1.0
    } else {
        limited
    }
}

fn limit_minutes(minutes: f64) -> f64 {
    let limited = minutes;
    
    if limited < -20.0 {
        limited + 1440.0
    } else if limited > 20.0 {
        limited - 1440.0
    } else {
        limited
    }
}

fn dayfrac_to_local_hr(dayfrac: f64, timezone: f64) -> f64 {
    24.0 * limit_zero2one(dayfrac + timezone / 24.0)    
}

fn third_order_polynomial(a: f64, b: f64, c: f64, d: f64, x: f64) -> f64 {
    ((a * x + b) * x + c) * x + d
}

