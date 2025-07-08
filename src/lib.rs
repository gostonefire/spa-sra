#[cfg(feature = "chrono_0_4")]
use chrono::{DateTime, Datelike, Offset, TimeZone, Timelike};
use crate::errors::{SpaError, MESSAGES};
use crate::spa::{Function, Input, SpaData};

mod earth_periodic_terms;
mod constants;
mod nutation_obliquity_periodic_terms;
mod utils;
pub mod spa;
pub mod errors;

/// Builder for creating an operational SpaData struct
///
/// Example usage which gives the same results as in the original C program test code:
/// ```rust
/// use spa_sra::SpaBuilder;
/// use spa_sra::errors::SpaError;
/// use spa_sra::spa::{Function, SpaData};
///
/// fn main() {
///
///     match get_spa() {
///         Err(e) => { println!("{}", e) },
///         Ok(spa) => {
///             println!("Julian Day:    {:.6}", spa.spa_za.jd);
///             println!("L:             {:.6e} degrees", spa.spa_za.l);
///             println!("B:             {:.6e} degrees", spa.spa_za.b);
///             println!("R:             {:.6} AU", spa.spa_za.r);
///             println!("H:             {:.6} degrees", spa.spa_za.h);
///             println!("Delta Psi:     {:.6e} degrees", spa.spa_za.del_psi);
///             println!("Delta Epsilon: {:.6e} degrees", spa.spa_za.del_epsilon);
///             println!("Epsilon:       {:.6} degrees", spa.spa_za.epsilon);
///
///             println!("Zenith:        {:.6} degrees", spa.spa_za.zenith);
///             println!("Azimuth:       {:.6} degrees", spa.spa_za.azimuth);
///             println!("Incidence:     {:.6} degrees", spa.spa_za_inc.incidence);
///
///             let mut min: f64 = 60.0 * (spa.spa_za_rts.sunrise - (spa.spa_za_rts.sunrise as i64) as f64);
///             let mut sec: f64 = 60.0 * (min - (min as i64) as f64);
///             println!("Sunrise:       {:0>2}:{:0>2}:{:0>2} Local Time", spa.spa_za_rts.sunrise as i64, min as i64, sec as i64);
///
///             min = 60.0 * (spa.spa_za_rts.sunset - (spa.spa_za_rts.sunset as i64) as f64);
///             sec = 60.0 * (min - (min as i64) as f64);
///             println!("Sunset:        {:0>2}:{:0>2}:{:0>2} Local Time", spa.spa_za_rts.sunset as i64, min as i64, sec as i64);
///         }
///     }
/// }
///
/// fn get_spa() -> Result<SpaData, SpaError<'static>> {
///
///     let spa = SpaBuilder::new()
///         .date(2003, 10, 17)?
///         .time(12, 30, 30.0)?
///         .timezone(-7.0)?
///         .lat_long(39.742476, -105.1786)?
///         .pressure_temperature(820.0, 11.0)?
///         .atmospheric_refraction(0.5667)?
///         .elevation(1830.14)?
///         .slope(30.0)?
///         .azimuth_rotation(-10.0)?
///         .delta_ut1(0.0)?
///         .delta_t(67.0)?
///         .build_and_calculate(Function::SpaAll)?;
///
///     Ok(spa)
/// }
/// ```
///
pub struct SpaBuilder {
    input: Input,
}

impl SpaBuilder {
    /// Creates a new SpaBuilder.
    ///
    /// All values besides a few need to be set, either using builder functions or manually via the [`Input`] struct fields
    /// after calling [`SpaBuilder::build`], at least within valid ranges.
    ///
    /// These fields have quite valid defaults, at least within the neighborhood of 2025-07-10, those are:
    /// * [`SpaBuilder::atmospheric_refraction`] / [`Input::atmos_refract`] - 0.5667 (which is a typical value)
    /// * [`SpaBuilder::delta_ut1`] / [`Input::delta_ut1`] - 0.1 (doesn't change that often)
    /// * [`SpaBuilder::delta_t`] / [`Input::delta_t`] - 69.084 (doesn't change that often)
    /// 
    pub fn new() -> Self {
        Self { input: Input::new() }
    }

    /// Creates a new SpaBuilder from an existing [SpaData] struct.
    ///
    /// Useful to for instance re-run calculations after just changing one or a few settings
    ///
    /// # Arguments
    ///
    /// * 'spa' - existing [SpaData] struct
    pub fn from_spa_data(spa: SpaData) -> Self {
        Self { input: spa.input }
    }
    
    /// Sets date
    /// 
    /// # Arguments
    /// 
    /// * 'year' - 4-digit year, valid range: -2000 to 6000
    /// * 'month' - 2-digit month, valid range: 1 to  12
    /// * 'day' - 2-digit day, valid range: 1 to  31
    pub fn date(mut self, year: i64, month: i64, day: i64) -> Result<Self, SpaError<'static>> {
        if year        < -2000   || year        > 6000   { return Err(SpaError{ code: 1, message: MESSAGES[1] }) };
        if month       < 1       || month       > 12     { return Err(SpaError{ code: 2, message: MESSAGES[2] }) };
        if day         < 1       || day         > 31     { return Err(SpaError{ code: 3, message: MESSAGES[3] }) };

        self.input.year = year;
        self.input.month = month;
        self.input.day = day;

        Ok(self)
    }

    /// Sets time
    ///
    /// # Arguments
    ///
    /// * 'hour' - Observer local hour, valid range: 0 to  24
    /// * 'minute' - Observer local minute, valid range: 0 to  59
    /// * 'second' - Observer local second, valid range: 0 to <60 (accepts fraction)
    pub fn time(mut self, hour: i64, minute: i64, second: f64) -> Result<Self, SpaError<'static>> {
        if hour        < 0       || hour        > 24     { return Err(SpaError{ code: 4, message: MESSAGES[4] }) };
        if minute      < 0       || minute      > 59     { return Err(SpaError{ code: 5, message: MESSAGES[5] }) };
        if second      < 0.0     || second      >=60.0   { return Err(SpaError{ code: 6, message: MESSAGES[6] }) };
        if hour        == 24     && minute      > 0      { return Err(SpaError{ code: 5, message: MESSAGES[5] }) };
        if hour        == 24     && second      > 0.0    { return Err(SpaError{ code: 6, message: MESSAGES[6] }) };

        self.input.hour = hour;
        self.input.minute = minute;
        self.input.second = second;

        Ok(self)
    }

    /// Sets all date, time and timezone fields to the [SpaData] struct from the given date_time parameter.
    ///
    /// This method is dependent on the feature "chrono_0_4" which will include the [chrono] crate.
    ///
    /// # Arguments
    ///
    /// * 'date_time' - a [DateTime] object including the time zone
    #[cfg(feature = "chrono_0_4")]
    pub fn date_time<T: TimeZone>(mut self, date_time: DateTime<T>) -> Result<Self, SpaError<'static>> {

        self.input.year = date_time.year() as i64;
        self.input.month = date_time.month() as i64;
        self.input.day = date_time.day() as i64;
        self.input.hour = date_time.hour() as i64;
        self.input.minute = date_time.minute() as i64;
        self.input.second = date_time.second() as f64 + date_time.nanosecond() as f64 / 1_000_000_000f64;

        self.input.timezone = date_time.offset().fix().local_minus_utc() as f64 / 3600.0;

        Ok(self)
    }

    /// Sets observer time zone (negative west of Greenwich)
    ///
    /// # Arguments
    /// 
    /// * 'timezone' - valid range: -18 to 18 hours
    pub fn timezone(mut self, timezone: f64) -> Result<Self, SpaError<'static>> {
        if timezone.abs()      > 18.0       { return Err(SpaError{ code: 8, message: MESSAGES[8] }) };

        self.input.timezone = timezone;
        
        Ok(self)
    }
    
    /// Sets observer latitude (negative south of the equator) and longitude (negative west of Greenwich)
    /// 
    /// # Arguments
    /// 
    /// * 'latitude' - valid range: -90 to 90 degrees
    /// * 'longitude' - valid range: -180 to 180 degrees
    pub fn lat_long(mut self, latitude: f64, longitude: f64) -> Result<Self, SpaError<'static>> {
        if longitude.abs()     > 180.0      { return Err(SpaError{ code: 9, message: MESSAGES[9] }) };
        if latitude.abs()      > 90.0       { return Err(SpaError{ code: 10, message: MESSAGES[10] }) };

        self.input.latitude = latitude;
        self.input.longitude = longitude;
        
        Ok(self)
    }
    
    /// Sets annual average local pressure \[millibars\] and annual average local temperature \[degrees Celsius\]
    /// 
    /// # Arguments
    /// 
    /// * 'pressure' - valid range: 0 to 5000 millibars
    /// * 'temperature' - valid range: -273 to 6000 degrees Celsius
    pub fn pressure_temperature(mut self, pressure: f64, temperature: f64) -> Result<Self, SpaError<'static>> {
        if pressure    < 0.0     || pressure    > 5000.0 { return Err(SpaError{ code: 12, message: MESSAGES[12] }) };
        if temperature <= -273.0 || temperature > 6000.0 { return Err(SpaError{ code: 13, message: MESSAGES[13] }) };

        self.input.pressure = pressure;
        self.input.temperature = temperature;
        
        Ok(self)
    }
    
    /// Sets atmospheric refraction at sunrise and sunset (0.5667 deg is typical)
    /// 
    /// # Arguments
    /// 
    /// * 'atmos_refract' - valid range: -5 to 5 degrees
    pub fn atmospheric_refraction(mut self, atmos_refract: f64) -> Result<Self, SpaError<'static>> {
        if atmos_refract.abs() > 5.0        { return Err(SpaError{ code: 16, message: MESSAGES[16] }) };

        self.input.atmos_refract = atmos_refract;
        
        Ok(self)
    }
    
    /// Sets observer elevation \[meters\]
    /// 
    /// # Arguments
    /// 
    /// * 'elevation' - valid range: -6500000 or higher meters
    pub fn elevation(mut self, elevation: f64) -> Result<Self, SpaError<'static>> {
        if elevation           < -6500000.0 { return Err(SpaError{ code: 11, message: MESSAGES[11] }) };

        self.input.elevation = elevation;
        
        Ok(self)
    }
    
    /// Sets surface slope (measured from the horizontal plane).
    ///
    /// This value is used to calculate the surface incidence angle for e.g. a solar panel.
    /// The surface incidence angle is the angle between an incoming ray (like light or radar) and
    /// a line perpendicular to the surface at the point where the ray hits.
    ///
    /// No need to set this unless [Function::SpaZaInc] or [Function::SpaAll] will be used
    /// when building.
    /// 
    /// # Arguments
    /// 
    /// * 'slope' - valid range: -360 to 360 degrees
    pub fn slope(mut self, slope: f64) -> Result<Self, SpaError<'static>> {
        if slope.abs()  > 360.0 { return Err(SpaError{ code: 14, message: MESSAGES[14] }) };

        self.input.slope = slope;
        
        Ok(self)
    }
    
    /// Sets surface azimuth rotation (measured from south to projection of
    /// surface normal on horizontal plane, negative east)
    ///
    /// This value is used to calculate the surface incidence angle for e.g. a solar panel.
    /// The surface incidence angle is the angle between an incoming ray (like light or radar) and
    /// a line perpendicular to the surface at the point where the ray hits.
    ///
    /// No need to set this unless [Function::SpaZaInc] or [Function::SpaAll] will be used
    /// when building.
    ///
    /// # Arguments
    /// 
    /// * 'azm_rotation' - -360 to 360 degrees
    pub fn azimuth_rotation(mut self, azm_rotation: f64) -> Result<Self, SpaError<'static>> {
        if azm_rotation.abs() > 360.0 { return Err(SpaError{ code: 15, message: MESSAGES[15] }) };

        self.input.azm_rotation = azm_rotation;
        
        Ok(self)
    }
    
    /// Fractional second difference between UTC and UT which is used
    /// to adjust UTC for earth's irregular rotation rate and is derived
    /// from observation only and is reported in this bulletin:
    /// <http://maia.usno.navy.mil/ser7/ser7.dat>,
    /// where delta_ut1 = DUT1
    /// 
    /// # Arguments
    /// 
    /// * 'delta_ut1' - valid range: -1 to 1 second (exclusive)
    pub fn delta_ut1(mut self, delta_ut1: f64) -> Result<Self, SpaError<'static>> {
        if delta_ut1   <= -1.0   || delta_ut1   >= 1.0   { return Err(SpaError{ code: 17, message: MESSAGES[17] }) };

        self.input.delta_ut1 = delta_ut1;
        
        Ok(self)
    }
    
    /// Difference between earth rotation time and terrestrial time
    /// It is derived from observation only and is reported in this
    /// bulletin: <http://maia.usno.navy.mil/ser7/ser7.dat>,
    /// where delta_t = 32.184 + (TAI-UTC) - DUT1
    /// 
    /// # Arguments
    /// 
    /// * 'delta_t' - valid range: -8000 to 8000 seconds
    pub fn delta_t(mut self, delta_t: f64) -> Result<Self, SpaError<'static>> {
        if delta_t.abs()       > 8000.0     { return Err(SpaError{ code: 7, message: MESSAGES[7] }) };

        self.input.delta_t = delta_t;
        
        Ok(self)
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
        self.input.function = function;

        SpaData::new(self.input)
    }

    /// Builds and calculates a SpaData struct
    ///
    /// See builder method [SpaBuilder::build] for more information
    ///
    ///# Arguments
    ///
    /// * 'function' - switch to choose functions for desired output (from enumeration [Function])
    pub fn build_and_calculate(mut self, function: Function) -> Result<SpaData, SpaError<'static>> {
        self.input.function = function;
        let mut spa_data = SpaData::new(self.input);
        spa_data.spa_calculate()?;

        Ok(spa_data)
    }
}


#[cfg(test)]
mod tests {
    use crate::spa::Function;
    use super::*;

    #[test]
    fn full_spa_output() {
        // enter required input values into SPA structure

        let input = Input {
            year:          2003,
            month:         10,
            day:           17,
            hour:          12,
            minute:        30,
            second:        30.0,
            timezone:      -7.0,
            delta_ut1:     0.0,
            delta_t:       67.0,
            longitude:     -105.1786,
            latitude:      39.742476,
            elevation:     1830.14,
            pressure:      820.0,
            temperature:   11.0,
            slope:         30.0,
            azm_rotation:  -10.0,
            atmos_refract: 0.5667,
            function:      Function::SpaAll,      
        };

        let mut spa = SpaData::new(input);

        // call the SPA calculate function and pass the SPA structure
        // test and test results according original code

        if let Err(e) = spa.spa_calculate() {
            panic!("{}", e);
        }

        assert_eq!(format!("{:.6}", spa.spa_za.jd), "2452930.312847",            "Julian Day:    2452930.312847");
        assert_eq!(format!("{:.6e}", spa.spa_za.l), "2.401826e1",                "L:             2.401826e+01 degrees");
        assert_eq!(format!("{:.6e}", spa.spa_za.b), "-1.011219e-4",              "B:             -1.011219e-04 degrees");
        assert_eq!(format!("{:.6}", spa.spa_za.r), "0.996542",                   "R:             0.996542 AU");
        assert_eq!(format!("{:.6}", spa.spa_za.h), "11.105902",                  "H:             11.105902 degrees");
        assert_eq!(format!("{:.6e}", spa.spa_za.del_psi), "-3.998404e-3",        "Delta Psi:     -3.998404e-03 degrees");
        assert_eq!(format!("{:.6e}", spa.spa_za.del_epsilon), "1.666568e-3",     "Delta Epsilon: 1.666568e-03 degrees");
        assert_eq!(format!("{:.6}", spa.spa_za.epsilon), "23.440465",            "Epsilon:       23.440465 degrees");
        assert_eq!(format!("{:.6}", spa.spa_za.zenith), "50.111622",             "Zenith:        50.111622 degrees");
        assert_eq!(format!("{:.6}", spa.spa_za.azimuth), "194.340241",           "Azimuth:       194.340241 degrees");
        assert_eq!(format!("{:.6}", spa.spa_za_inc.incidence), "25.187000",      "Incidence:     25.187000 degrees");

        let mut min: f64 = 60.0 * (spa.spa_za_rts.sunrise - (spa.spa_za_rts.sunrise as i64) as f64);
        let mut sec: f64 = 60.0 * (min - (min as i64) as f64);
        assert_eq!(format!("{:0>2}:{:0>2}:{:0>2}", spa.spa_za_rts.sunrise as i64, min as i64, sec as i64), "06:12:43",
                   "Sunrise: 06:12:43 Local Time");

        min = 60.0 * (spa.spa_za_rts.sunset - (spa.spa_za_rts.sunset as i64) as f64);
        sec = 60.0 * (min - (min as i64) as f64);
        assert_eq!(format!("{:0>2}:{:0>2}:{:0>2}", spa.spa_za_rts.sunset as i64, min as i64, sec as i64), "17:20:19",
                   "Sunset: 17:20:19 Local Time");
    }
}