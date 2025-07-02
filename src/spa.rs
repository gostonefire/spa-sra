
/// Enumeration for function codes to select desired final outputs from SPA
enum Output{
    /// calculate zenith and azimuth
    SpaZa,
    /// calculate zenith, azimuth, and incidence
    SpaZaInc,
    /// calculate zenith, azimuth, and sun rise/transit/set values
    SpaZaRts,
    /// calculate all SPA output values
    SpaAll,
}

struct SpaData {
    //----------------------INPUT VALUES------------------------

    /// 4-digit year,      valid range: -2000 to 6000, error code: 1
    year: i64,
    /// 2-digit month,         valid range: 1 to  12,  error code: 2
    month: i64,
    /// 2-digit day,           valid range: 1 to  31,  error code: 3
    day: i64,
    /// Observer local hour,   valid range: 0 to  24,  error code: 4
    hour: i64,
    /// Observer local minute, valid range: 0 to  59,  error code: 5
    minute: i64,
    /// Observer local second, valid range: 0 to <60,  error code: 6
    second: f64,

    /// Fractional second difference between UTC and UT which is used
    /// to adjust UTC for earth's irregular rotation rate and is derived
    /// from observation only and is reported in this bulletin:
    /// http://maia.usno.navy.mil/ser7/ser7.dat,
    /// where delta_ut1 = DUT1
    /// 
    /// valid range: -1 to 1 second (exclusive), error code 17
    delta_ut1: f64,

    /// Difference between earth rotation time and terrestrial time
    /// It is derived from observation only and is reported in this
    /// bulletin: http://maia.usno.navy.mil/ser7/ser7.dat,
    /// where delta_t = 32.184 + (TAI-UTC) - DUT1
    /// 
    /// valid range: -8000 to 8000 seconds, error code: 7
    delta_t: f64,

    /// Observer time zone (negative west of Greenwich)
    /// 
    /// valid range: -18   to   18 hours,   error code: 8
    timezone: f64,

    /// Observer longitude (negative west of Greenwich)
    /// 
    /// valid range: -180  to  180 degrees, error code: 9
    longitude: f64,

    /// Observer latitude (negative south of equator)
    /// 
    /// valid range: -90   to   90 degrees, error code: 10
    latitude: f64,

    /// Observer elevation [meters]
    /// 
    /// valid range: -6500000 or higher meters,    error code: 11
    elevation: f64,

    /// Annual average local pressure [millibars]
    /// 
    /// valid range:    0 to 5000 millibars,       error code: 12
    pressure: f64,

    /// Annual average local temperature [degrees Celsius]
    /// 
    /// valid range: -273 to 6000 degrees Celsius, error code; 13
    temperature: f64,

    /// Surface slope (measured from the horizontal plane)
    /// 
    /// valid range: -360 to 360 degrees, error code: 14
    slope: f64,

    /// Surface azimuth rotation (measured from south to projection of
    /// surface normal on horizontal plane, negative east)
    /// 
    /// valid range: -360 to 360 degrees, error code: 15
    azm_rotation: f64,

    /// Atmospheric refraction at sunrise and sunset (0.5667 deg is typical)
    /// 
    /// valid range: -5   to   5 degrees, error code: 16
    atmos_refract: f64,

    /// Switch to choose functions for desired output (from enumeration)
    function: Output,

    //-----------------Intermediate OUTPUT VALUES--------------------

    jd: f64,             //Julian day
    jc: f64,             //Julian century

    jde: f64,            //Julian ephemeris day
    jce: f64,            //Julian ephemeris century
    jme: f64,            //Julian ephemeris millennium

    l: f64,              //earth heliocentric longitude [degrees]
    b: f64,              //earth heliocentric latitude [degrees]
    r: f64,              //earth radius vector [Astronomical Units, AU]

    theta: f64,          //geocentric longitude [degrees]
    beta: f64,           //geocentric latitude [degrees]

    x0: f64,             //mean elongation (moon-sun) [degrees]
    x1: f64,             //mean anomaly (sun) [degrees]
    x2: f64,             //mean anomaly (moon) [degrees]
    x3: f64,             //argument latitude (moon) [degrees]
    x4: f64,             //ascending longitude (moon) [degrees]

    del_psi: f64,        //nutation longitude [degrees]
    del_epsilon: f64,    //nutation obliquity [degrees]
    epsilon0: f64,       //ecliptic mean obliquity [arc seconds]
    epsilon: f64,        //ecliptic true obliquity  [degrees]

    del_tau: f64,        //aberration correction [degrees]
    lamda: f64,          //apparent sun longitude [degrees]
    nu0: f64,            //Greenwich mean sidereal time [degrees]
    nu: f64,             //Greenwich sidereal time [degrees]

    alpha: f64,          //geocentric sun right ascension [degrees]
    delta: f64,          //geocentric sun declination [degrees]

    h: f64,              //observer hour angle [degrees]
    xi: f64,             //sun equatorial horizontal parallax [degrees]
    del_alpha: f64,      //sun right ascension parallax [degrees]
    delta_prime: f64,    //topocentric sun declination [degrees]
    alpha_prime: f64,    //topocentric sun right ascension [degrees]
    h_prime: f64,        //topocentric local hour angle [degrees]

    e0: f64,             //topocentric elevation angle (uncorrected) [degrees]
    del_e: f64,          //atmospheric refraction correction [degrees]
    e: f64,              //topocentric elevation angle (corrected) [degrees]

    eot: f64,            //equation of time [minutes]
    srha: f64,           //sunrise hour angle [degrees]
    ssha: f64,           //sunset hour angle [degrees]
    sta: f64,            //sun transit altitude [degrees]

    //---------------------Final OUTPUT VALUES------------------------

    zenith: f64,         //topocentric zenith angle [degrees]
    azimuth_astro: f64,  //topocentric azimuth angle (westward from south) [for astronomers]
    azimuth: f64,        //topocentric azimuth angle (eastward from north) [for navigators and solar radiation]
    incidence: f64,      //surface incidence angle [degrees]

    suntransit: f64,     //local sun transit time (or solar noon) [fractional hour]
    sunrise: f64,        //local sunrise time (+/- 30 seconds) [fractional hour]
    sunset: f64,         //local sunset time (+/- 30 seconds) [fractional hour]
}
