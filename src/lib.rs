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
        let mut spa = SpaData {
            year: 0,
            month: 0,
            day: 0,
            hour: 0,
            minute: 0,
            second: 0.0,
            delta_ut1: 0.0,
            delta_t: 0.0,
            timezone: 0.0,
            longitude: 0.0,
            latitude: 0.0,
            elevation: 0.0,
            pressure: 0.0,
            temperature: 0.0,
            slope: 0.0,
            azm_rotation: 0.0,
            atmos_refract: 0.0,
            function: Output::SpaAll,
            jd: 0.0,
            jc: 0.0,
            jde: 0.0,
            jce: 0.0,
            jme: 0.0,
            l: 0.0,
            b: 0.0,
            r: 0.0,
            theta: 0.0,
            beta: 0.0,
            x0: 0.0,
            x1: 0.0,
            x2: 0.0,
            x3: 0.0,
            x4: 0.0,
            del_psi: 0.0,
            del_epsilon: 0.0,
            epsilon0: 0.0,
            epsilon: 0.0,
            del_tau: 0.0,
            lamda: 0.0,
            nu0: 0.0,
            nu: 0.0,
            alpha: 0.0,
            delta: 0.0,
            h: 0.0,
            xi: 0.0,
            del_alpha: 0.0,
            delta_prime: 0.0,
            alpha_prime: 0.0,
            h_prime: 0.0,
            e0: 0.0,
            del_e: 0.0,
            e: 0.0,
            eot: 0.0,
            srha: 0.0,
            ssha: 0.0,
            sta: 0.0,
            zenith: 0.0,
            azimuth_astro: 0.0,
            azimuth: 0.0,
            incidence: 0.0,
            suntransit: 0.0,
            sunrise: 0.0,
            sunset: 0.0,
        };

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
