use crate::functions::spa_calculate;
use crate::spa::SpaData;

mod earth_periodic_terms;
mod constants;
mod nutation_obliquity_periodic_terms;
mod functions;
pub mod spa;

pub fn spa_calc(spa: &mut SpaData) -> i64  {
    spa_calculate(spa)
}

#[cfg(test)]
mod tests {
    use crate::spa::Output;
    use super::*;

    #[test]
    fn full_spa_output() {
        let mut spa = SpaData::new();

        // enter required input values into SPA structure

        spa.year          = 2003;
        spa.month         = 10;
        spa.day           = 17;
        spa.hour          = 12;
        spa.minute        = 30;
        spa.second        = 30.0;
        spa.timezone      = -7.0;
        spa.delta_ut1     = 0.0;
        spa.delta_t       = 67.0;
        spa.longitude     = -105.1786;
        spa.latitude      = 39.742476;
        spa.elevation     = 1830.14;
        spa.pressure      = 820.0;
        spa.temperature   = 11.0;
        spa.slope         = 30.0;
        spa.azm_rotation  = -10.0;
        spa.atmos_refract = 0.5667;
        spa.function      = Output::SpaAll;

        // call the SPA calculate function and pass the SPA structure
        // test and test results according original code

        let result = spa_calc(&mut spa);
        assert_eq!(result, 0, "SPA Error Code: 0");
        
        assert_eq!(format!("{:.6}", spa.jd), "2452930.312847",        "Julian Day:    2452930.312847");
        assert_eq!(format!("{:.6e}", spa.l), "2.401826e1",            "L:             2.401826e+01 degrees");
        assert_eq!(format!("{:.6e}", spa.b), "-1.011219e-4",          "B:             -1.011219e-04 degrees");
        assert_eq!(format!("{:.6}", spa.r), "0.996542",               "R:             0.996542 AU");
        assert_eq!(format!("{:.6}", spa.h), "11.105902",              "H:             11.105902 degrees");
        assert_eq!(format!("{:.6e}", spa.del_psi), "-3.998404e-3",    "Delta Psi:     -3.998404e-03 degrees");
        assert_eq!(format!("{:.6e}", spa.del_epsilon), "1.666568e-3", "Delta Epsilon: 1.666568e-03 degrees");
        assert_eq!(format!("{:.6}", spa.epsilon), "23.440465",        "Epsilon:       23.440465 degrees");
        assert_eq!(format!("{:.6}", spa.zenith), "50.111622",         "Zenith:        50.111622 degrees");
        assert_eq!(format!("{:.6}", spa.azimuth), "194.340241",       "Azimuth:       194.340241 degrees");
        assert_eq!(format!("{:.6}", spa.incidence), "25.187000",      "Incidence:     25.187000 degrees");

        let mut min: f64 = 60.0 * (spa.sunrise - (spa.sunrise as i64) as f64);
        let mut sec: f64 = 60.0 * (min - (min as i64) as f64);
        assert_eq!(format!("{:0>2}:{:0>2}:{:0>2}", spa.sunrise as i64, min as i64, sec as i64), "06:12:43", 
                   "Sunrise: 06:12:43 Local Time");

        min = 60.0 * (spa.sunset - (spa.sunset as i64) as f64);
        sec = 60.0 * (min - (min as i64) as f64);
        assert_eq!(format!("{:0>2}:{:0>2}:{:0>2}", spa.sunset as i64, min as i64, sec as i64), "17:20:19", 
                   "Sunset: 17:20:19 Local Time");
    }
}
