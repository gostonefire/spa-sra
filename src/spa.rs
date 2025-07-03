
/// Enumeration for function codes to select desired final outputs from SPA
#[derive(PartialEq, Clone)]
pub enum Output{
    /// calculate zenith and azimuth
    SpaZa,
    /// calculate zenith, azimuth, and incidence
    SpaZaInc,
    /// calculate zenith, azimuth, and sun rise/transit/set values
    SpaZaRts,
    /// calculate all SPA output values
    SpaAll,
}

#[derive(Clone)]
pub struct SpaData {
    //----------------------INPUT VALUES------------------------

    /// 4-digit year,      valid range: -2000 to 6000, error code: 1
    pub year: i64,
    /// 2-digit month,         valid range: 1 to  12,  error code: 2
    pub month: i64,
    /// 2-digit day,           valid range: 1 to  31,  error code: 3
    pub day: i64,
    /// Observer local hour,   valid range: 0 to  24,  error code: 4
    pub hour: i64,
    /// Observer local minute, valid range: 0 to  59,  error code: 5
    pub minute: i64,
    /// Observer local second, valid range: 0 to <60,  error code: 6
    pub second: f64,

    /// Fractional second difference between UTC and UT which is used
    /// to adjust UTC for earth's irregular rotation rate and is derived
    /// from observation only and is reported in this bulletin:
    /// http://maia.usno.navy.mil/ser7/ser7.dat,
    /// where delta_ut1 = DUT1
    /// 
    /// valid range: -1 to 1 second (exclusive), error code 17
    pub delta_ut1: f64,

    /// Difference between earth rotation time and terrestrial time
    /// It is derived from observation only and is reported in this
    /// bulletin: http://maia.usno.navy.mil/ser7/ser7.dat,
    /// where delta_t = 32.184 + (TAI-UTC) - DUT1
    /// 
    /// valid range: -8000 to 8000 seconds, error code: 7
    pub delta_t: f64,

    /// Observer time zone (negative west of Greenwich)
    /// 
    /// valid range: -18   to   18 hours,   error code: 8
    pub timezone: f64,

    /// Observer longitude (negative west of Greenwich)
    /// 
    /// valid range: -180  to  180 degrees, error code: 9
    pub longitude: f64,

    /// Observer latitude (negative south of equator)
    /// 
    /// valid range: -90   to   90 degrees, error code: 10
    pub latitude: f64,

    /// Observer elevation [meters]
    /// 
    /// valid range: -6500000 or higher meters,    error code: 11
    pub elevation: f64,

    /// Annual average local pressure [millibars]
    /// 
    /// valid range:    0 to 5000 millibars,       error code: 12
    pub pressure: f64,

    /// Annual average local temperature [degrees Celsius]
    /// 
    /// valid range: -273 to 6000 degrees Celsius, error code; 13
    pub temperature: f64,

    /// Surface slope (measured from the horizontal plane)
    /// 
    /// valid range: -360 to 360 degrees, error code: 14
    pub slope: f64,

    /// Surface azimuth rotation (measured from south to projection of
    /// surface normal on horizontal plane, negative east)
    /// 
    /// valid range: -360 to 360 degrees, error code: 15
    pub azm_rotation: f64,

    /// Atmospheric refraction at sunrise and sunset (0.5667 deg is typical)
    /// 
    /// valid range: -5   to   5 degrees, error code: 16
    pub atmos_refract: f64,

    /// Switch to choose functions for desired output (from enumeration)
    pub function: Output,

    //-----------------Intermediate OUTPUT VALUES--------------------

    pub(crate) jd: f64,             //Julian day
    pub(crate) jc: f64,             //Julian century

    pub(crate) jde: f64,            //Julian ephemeris day
    pub(crate) jce: f64,            //Julian ephemeris century
    pub(crate) jme: f64,            //Julian ephemeris millennium

    pub(crate) l: f64,              //earth heliocentric longitude [degrees]
    pub(crate) b: f64,              //earth heliocentric latitude [degrees]
    pub(crate) r: f64,              //earth radius vector [Astronomical Units, AU]

    pub(crate) theta: f64,          //geocentric longitude [degrees]
    pub(crate) beta: f64,           //geocentric latitude [degrees]

    pub(crate) x0: f64,             //mean elongation (moon-sun) [degrees]
    pub(crate) x1: f64,             //mean anomaly (sun) [degrees]
    pub(crate) x2: f64,             //mean anomaly (moon) [degrees]
    pub(crate) x3: f64,             //argument latitude (moon) [degrees]
    pub(crate) x4: f64,             //ascending longitude (moon) [degrees]

    pub(crate) del_psi: f64,        //nutation longitude [degrees]
    pub(crate) del_epsilon: f64,    //nutation obliquity [degrees]
    pub(crate) epsilon0: f64,       //ecliptic mean obliquity [arc seconds]
    pub(crate) epsilon: f64,        //ecliptic true obliquity  [degrees]

    pub(crate) del_tau: f64,        //aberration correction [degrees]
    pub(crate) lamda: f64,          //apparent sun longitude [degrees]
    pub(crate) nu0: f64,            //Greenwich mean sidereal time [degrees]
    pub(crate) nu: f64,             //Greenwich sidereal time [degrees]

    pub(crate) alpha: f64,          //geocentric sun right ascension [degrees]
    pub(crate) delta: f64,          //geocentric sun declination [degrees]

    pub(crate) h: f64,              //observer hour angle [degrees]
    pub(crate) xi: f64,             //sun equatorial horizontal parallax [degrees]
    pub(crate) del_alpha: f64,      //sun right ascension parallax [degrees]
    pub(crate) delta_prime: f64,    //topocentric sun declination [degrees]
    pub(crate) alpha_prime: f64,    //topocentric sun right ascension [degrees]
    pub(crate) h_prime: f64,        //topocentric local hour angle [degrees]

    pub(crate) e0: f64,             //topocentric elevation angle (uncorrected) [degrees]
    pub(crate) del_e: f64,          //atmospheric refraction correction [degrees]
    pub(crate) e: f64,              //topocentric elevation angle (corrected) [degrees]

    pub(crate) eot: f64,            //equation of time [minutes]
    pub(crate) srha: f64,           //sunrise hour angle [degrees]
    pub(crate) ssha: f64,           //sunset hour angle [degrees]
    pub(crate) sta: f64,            //sun transit altitude [degrees]

    //---------------------Final OUTPUT VALUES------------------------

    pub(crate) zenith: f64,         //topocentric zenith angle [degrees]
    pub(crate) azimuth_astro: f64,  //topocentric azimuth angle (westward from south) [for astronomers]
    pub(crate) azimuth: f64,        //topocentric azimuth angle (eastward from north) [for navigators and solar radiation]
    pub(crate) incidence: f64,      //surface incidence angle [degrees]

    pub(crate) suntransit: f64,     //local sun transit time (or solar noon) [fractional hour]
    pub(crate) sunrise: f64,        //local sunrise time (+/- 30 seconds) [fractional hour]
    pub(crate) sunset: f64,         //local sunset time (+/- 30 seconds) [fractional hour]
}

impl SpaData {
    pub fn new() -> Self {
        SpaData {
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
        }
    }
}
