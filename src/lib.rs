use crate::spa::{Function, SpaData};

mod earth_periodic_terms;
mod constants;
mod nutation_obliquity_periodic_terms;
mod utils;
pub mod spa;

/// Builder for creating an operational SpaData struct
/// 
pub struct SpaBuilder {
    spa: SpaData,
}

impl SpaBuilder {
    /// Creates a new SpaBuilder with default values as of original code tests, hence an immediate
    /// call to the `build` function with argument Function::SpaAll returns a fully operational
    /// SpaData struct which when invoking the XXXXX will output those test results.
    ///
    /// To change any settings, just use the builder object and call whatever change is needed and
    /// then call the `build` function to get the new SpaData struct ready for execution with the
    /// new settings.
    /// 
    pub fn new() -> Self {
        Self { spa: SpaData::new() }
    }
    
    /// Sets date and time
    /// 
    /// # Arguments
    /// 
    /// * 'year' - 4-digit year, valid range: -2000 to 6000
    /// * 'month' - 2-digit month, valid range: 1 to  12
    /// * 'day' - 2-digit day, valid range: 1 to  31
    /// * 'hour' - Observer local hour, valid range: 0 to  24
    /// * 'minute' - Observer local minute, valid range: 0 to  59
    /// * 'second' - Observer local second, valid range: 0 to <60 (accepts fraction)
    pub fn date_time(mut self, year: i64, month: i64, day: i64, hour: i64, minute: i64, second: f64) -> Self {
        self.spa.year = year;
        self.spa.month = month;
        self.spa.day = day;
        self.spa.hour = hour;
        self.spa.minute = minute;
        self.spa.second = second;
        
        self
    }
    
    /// Sets observer time zone (negative west of Greenwich)
    ///
    /// # Arguments
    /// 
    /// * 'timezone' - valid range: -18 to 18 hours
    pub fn timezone(mut self, timezone: f64) -> Self {
        self.spa.timezone = timezone;
        
        self
    }
    
    /// Sets observer latitude (negative south of equator) and longitude (negative west of Greenwich)
    /// 
    /// # Arguments
    /// 
    /// * 'latitude' - valid range: -90 to 90 degrees
    /// * 'longitude' - valid range: -180 to 180 degrees
    pub fn lat_long(mut self, latitude: f64, longitude: f64) -> Self {
        self.spa.latitude = latitude;
        self.spa.longitude = longitude;
        
        self
    }
    
    /// Sets annual average local pressure [millibars] and annual average local temperature [degrees Celsius]
    /// 
    /// # Arguments
    /// 
    /// * 'pressure' - valid range: 0 to 5000 millibars
    /// * 'temperature' - valid range: -273 to 6000 degrees Celsius
    pub fn pressure_temperature(mut self, pressure: f64, temperature: f64) -> Self {
        self.spa.pressure = pressure;
        self.spa.temperature = temperature;
        
        self
    }
    
    /// Sets atmospheric refraction at sunrise and sunset (0.5667 deg is typical)
    /// 
    /// # Arguments
    /// 
    /// * 'atmos_refract' - valid range: -5 to 5 degrees
    pub fn atmospheric_refraction(mut self, atmos_refract: f64) -> Self {
        self.spa.atmos_refract = atmos_refract;
        
        self
    }
    
    /// Sets observer elevation [meters]
    /// 
    /// # Arguments
    /// 
    /// * 'elevation' - valid range: -6500000 or higher meters
    pub fn elevation(mut self, elevation: f64) -> Self {
        self.spa.elevation = elevation;
        
        self
    }
    
    /// Sets surface slope (measured from the horizontal plane)
    /// 
    /// # Arguments
    /// 
    /// * 'slope' - valid range: -360 to 360 degrees
    pub fn slope(mut self, slope: f64) -> Self {
        self.spa.slope = slope;
        
        self
    }
    
    /// Sets surface azimuth rotation (measured from south to projection of
    /// surface normal on horizontal plane, negative east)
    /// 
    /// # Arguments
    /// 
    /// * 'azm_rotation' - -360 to 360 degrees
    pub fn azimuth_rotation(mut self, azm_rotation: f64) -> Self {
        self.spa.azm_rotation = azm_rotation;
        
        self
    }
    
    /// Fractional second difference between UTC and UT which is used
    /// to adjust UTC for earth's irregular rotation rate and is derived
    /// from observation only and is reported in this bulletin:
    /// http://maia.usno.navy.mil/ser7/ser7.dat,
    /// where delta_ut1 = DUT1
    /// 
    /// # Arguments
    /// 
    /// * 'delta_ut1' - valid range: -1 to 1 second (exclusive)
    pub fn delta_ut1(mut self, delta_ut1: f64) -> Self {
        self.spa.delta_ut1 = delta_ut1;
        
        self
    }
    
    /// Difference between earth rotation time and terrestrial time
    /// It is derived from observation only and is reported in this
    /// bulletin: http://maia.usno.navy.mil/ser7/ser7.dat,
    /// where delta_t = 32.184 + (TAI-UTC) - DUT1
    /// 
    /// # Arguments
    /// 
    /// * 'delta_t' - valid range: -8000 to 8000 seconds
    pub fn delta_t(mut self, delta_t: f64) -> Self {
        self.spa.delta_t = delta_t;
        
        self
    }

    /// Builds and return a SpaData struct ready to execute.
    /// It is important to call all builder functions before this since the default SpaData struct
    /// will operate on settings as defined in the original code test, which reflect some situation
    /// and place as of 2003.
    ///
    /// The argument to this function sets what output is expected:
    /// * Function::SpaZa    - calculate zenith and azimuth
    /// * Function::SpaZaInc - calculate zenith, azimuth, and incidence
    /// * Function::SpaZaRts - calculate zenith, azimuth, and sun rise/transit/set values
    /// * Function::SpaAll   - calculate all SPA output values
    ///
    /// # Arguments
    ///
    /// * 'function' - switch to choose functions for desired output (from enumeration `Function`)
    pub fn build(mut self, function: Function) -> SpaData {
        self.spa.function = function;

        self.spa.clone()
    }
}

pub fn spa_calc(spa: &mut SpaData) -> i64  {
    spa.spa_calculate()
    //spa_calculate(spa)
}

#[cfg(test)]
mod tests {
    use crate::spa::Function;
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
        spa.function      = Function::SpaAll;

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