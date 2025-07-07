use crate::constants::{B_COUNT, B_SUBCOUNT, JD_COUNT, JD_MINUS, JD_PLUS, JD_ZERO, L_COUNT, L_SUBCOUNT, R_COUNT, R_SUBCOUNT, SUN_COUNT, SUN_RADIUS, SUN_RISE, SUN_SET, SUN_TRANSIT, TERM_A, TERM_B, TERM_C, TERM_COUNT, TERM_EPS_C, TERM_EPS_D, TERM_PSI_A, TERM_PSI_B, TERM_X0, TERM_X1, TERM_X2, TERM_X3, TERM_X4, TERM_X_COUNT, TERM_Y_COUNT, Y_COUNT};
use crate::earth_periodic_terms::{B_TERMS, L_TERMS, R_TERMS};
use crate::utils;
use crate::utils::{atmospheric_refraction_correction, deg2rad, geocentric_declination, geocentric_right_ascension, limit_degrees, observer_hour_angle, rad2deg, right_ascension_parallax_and_topocentric_dec, third_order_polynomial, topocentric_azimuth_angle, topocentric_azimuth_angle_astro, topocentric_elevation_angle, topocentric_elevation_angle_corrected, topocentric_local_hour_angle, topocentric_right_ascension};
use crate::nutation_obliquity_periodic_terms::{PE_TERMS, Y_TERMS};

/// Enumeration for function codes to select desired final outputs from SPA
#[derive(PartialEq, Clone)]
pub enum Function {
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
    pub function: Function,

    //-----------------Intermediate OUTPUT VALUES--------------------

    /// Julian day
    pub(crate) jd: f64,
    /// Julian century
    pub(crate) jc: f64,

    /// Julian ephemeris day
    pub(crate) jde: f64,
    /// Julian ephemeris century
    pub(crate) jce: f64,
    /// Julian ephemeris millennium
    pub(crate) jme: f64,

    /// earth heliocentric longitude [degrees]
    pub(crate) l: f64,
    /// earth heliocentric latitude [degrees]
    pub(crate) b: f64,
    /// earth radius vector [Astronomical Units, AU]
    pub(crate) r: f64,

    /// geocentric longitude [degrees]
    pub(crate) theta: f64,
    /// geocentric latitude [degrees]
    pub(crate) beta: f64,

    /// mean elongation (moon-sun) [degrees]
    pub(crate) x0: f64,
    /// mean anomaly (sun) [degrees]
    pub(crate) x1: f64,
    /// mean anomaly (moon) [degrees]
    pub(crate) x2: f64,
    /// argument latitude (moon) [degrees]
    pub(crate) x3: f64,
    /// ascending longitude (moon) [degrees]
    pub(crate) x4: f64,

    /// nutation longitude [degrees]
    pub(crate) del_psi: f64,
    /// nutation obliquity [degrees]
    pub(crate) del_epsilon: f64,
    /// ecliptic mean obliquity [arc seconds]
    pub(crate) epsilon0: f64,
    /// ecliptic true obliquity  [degrees]
    pub(crate) epsilon: f64,

    /// aberration correction [degrees]
    pub(crate) del_tau: f64,
    /// apparent sun longitude [degrees]
    pub(crate) lamda: f64,
    /// Greenwich mean sidereal time [degrees]
    pub(crate) nu0: f64,
    /// Greenwich sidereal time [degrees]
    pub(crate) nu: f64,

    /// geocentric sun right ascension [degrees]
    pub(crate) alpha: f64,
    /// geocentric sun declination [degrees]
    pub(crate) delta: f64,

    /// observer hour angle [degrees]
    pub(crate) h: f64,
    /// sun equatorial horizontal parallax [degrees]
    pub(crate) xi: f64,
    /// sun right ascension parallax [degrees]
    pub(crate) del_alpha: f64,
    /// topocentric sun declination [degrees]
    pub(crate) delta_prime: f64,
    /// topocentric sun right ascension [degrees]
    pub(crate) alpha_prime: f64,
    /// topocentric local hour angle [degrees]
    pub(crate) h_prime: f64,

    /// topocentric elevation angle (uncorrected) [degrees]
    pub(crate) e0: f64,
    /// atmospheric refraction correction [degrees]
    pub(crate) del_e: f64,
    /// topocentric elevation angle (corrected) [degrees]
    pub(crate) e: f64,

    /// equation of time [minutes]
    pub(crate) eot: f64,
    /// sunrise hour angle [degrees]
    pub(crate) srha: f64,
    /// sunset hour angle [degrees]
    pub(crate) ssha: f64,
    /// sun transit altitude [degrees]
    pub(crate) sta: f64,

    //---------------------Final OUTPUT VALUES------------------------

    /// topocentric zenith angle [degrees]
    pub zenith: f64,
    /// topocentric azimuth angle (westward from south) [for astronomers]
    pub azimuth_astro: f64,
    /// topocentric azimuth angle (eastward from north) [for navigators and solar radiation]
    pub azimuth: f64,
    /// surface incidence angle [degrees]
    pub incidence: f64,

    /// local sun transit time (or solar noon) [fractional hour]
    pub suntransit: f64,
    /// local sunrise time (+/- 30 seconds) [fractional hour]
    pub sunrise: f64,
    /// local sunset time (+/- 30 seconds) [fractional hour]
    pub sunset: f64,
}

impl SpaData {
    pub fn new() -> Self {
        SpaData {
            year: 2003,
            month: 10,
            day: 17,
            hour: 12,
            minute: 30,
            second: 30.0,
            delta_ut1: 0.0,
            delta_t: 67.0,
            timezone: -7.0,
            longitude: -105.1786,
            latitude: 39.742476,
            elevation: 1830.14,
            pressure: 820.0,
            temperature: 11.0,
            slope: 30.0,
            azm_rotation: -10.0,
            atmos_refract: 0.5667,
            function: Function::SpaAll,
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

    /// Calculate all SPA parameters and put into structure
    /// Note: All inputs values (listed in header file) must already be in structure
    ///
    pub fn spa_calculate(&mut self) -> i64 {
        let result: i64 = self.validate_inputs();

        if result == 0 {
            self.jd = julian_day(self.year, self.month, self.day, self.hour,
                                 self.minute, self.second, self.delta_ut1, self.timezone);

            self.calculate_geocentric_sun_right_ascension_and_declination();

            self.h  = observer_hour_angle(self.nu, self.longitude, self.alpha);
            self.xi = sun_equatorial_horizontal_parallax(self.r);

            right_ascension_parallax_and_topocentric_dec(self.latitude, self.elevation, self.xi,
                                                         self.h, self.delta, &mut self.del_alpha, &mut self.delta_prime);

            self.alpha_prime = topocentric_right_ascension(self.alpha, self.del_alpha);
            self.h_prime     = topocentric_local_hour_angle(self.h, self.del_alpha);

            self.e0      = topocentric_elevation_angle(self.latitude, self.delta_prime, self.h_prime);
            self.del_e   = atmospheric_refraction_correction(self.pressure, self.temperature,
                                                             self.atmos_refract, self.e0);
            self.e       = topocentric_elevation_angle_corrected(self.e0, self.del_e);

            self.zenith        = utils::topocentric_zenith_angle(self.e);
            self.azimuth_astro = topocentric_azimuth_angle_astro(self.h_prime, self.latitude,
                                                                 self.delta_prime);
            self.azimuth       = topocentric_azimuth_angle(self.azimuth_astro);

            if self.function == Function::SpaZaInc || self.function == Function::SpaAll {
                self.incidence  = surface_incidence_angle(self.zenith, self.azimuth_astro,
                                                          self.azm_rotation, self.slope);
            }

            if self.function == Function::SpaZaRts || self.function == Function::SpaAll {
                self.calculate_eot_and_sun_rise_transit_set();
            }
        }

        result
    }

    /// Validates inputs
    ///
    fn validate_inputs(&self) -> i64 {
        if self.year        < -2000   || self.year        > 6000   { return 1 };
        if self.month       < 1       || self.month       > 12     { return 2 };
        if self.day         < 1       || self.day         > 31     { return 3 };
        if self.hour        < 0       || self.hour        > 24     { return 4 };
        if self.minute      < 0       || self.minute      > 59     { return 5 };
        if self.second      < 0.0     || self.second      >=60.0   { return 6 };
        if self.pressure    < 0.0     || self.pressure    > 5000.0 { return 12 };
        if self.temperature <= -273.0 || self.temperature > 6000.0 { return 13 };
        if self.delta_ut1   <= -1.0   || self.delta_ut1   >= 1.0   { return 17 };
        if self.hour        == 24     && self.minute      > 0      { return 5 };
        if self.hour        == 24     && self.second      > 0.0    { return 6 };

        if self.delta_t.abs()       > 8000.0     { return 7 };
        if self.timezone.abs()      > 18.0       { return 8 };
        if self.longitude.abs()     > 180.0      { return 9 };
        if self.latitude.abs()      > 90.0       { return 10 };
        if self.atmos_refract.abs() > 5.0        { return 16 };
        if self.elevation           < -6500000.0 { return 11 };

        if self.function == Function::SpaZaInc || self.function == Function::SpaAll {
            if self.slope.abs()  > 360.0 { return 14 };
            if self.azm_rotation > 360.0 { return 15 };
        }

        0
    }
    /// Calculate Equation of Time (EOT) and Sun Rise, Transit, & Set (RTS)
    ///
    fn calculate_eot_and_sun_rise_transit_set(&mut self) {
        let mut sun_rts: SpaData = self.clone();

        let mut alpha: [f64;JD_COUNT] = [0.0; JD_COUNT];
        let mut delta: [f64;JD_COUNT] = [0.0; JD_COUNT];

        let mut m_rts: [f64;SUN_COUNT] = [0.0; SUN_COUNT];
        let mut nu_rts: [f64;SUN_COUNT] = [0.0; SUN_COUNT];
        let mut h_rts: [f64;SUN_COUNT] = [0.0; SUN_COUNT];

        let mut alpha_prime: [f64;SUN_COUNT] = [0.0; SUN_COUNT];
        let mut delta_prime: [f64;SUN_COUNT] = [0.0; SUN_COUNT];
        let mut h_prime: [f64;SUN_COUNT] = [0.0; SUN_COUNT];

        let h0_prime: f64 = -1.0 * (SUN_RADIUS + self.atmos_refract);

        let m: f64 = sun_mean_longitude(self.jme);
        self.eot = eot(m, self.alpha, self.del_psi, self.epsilon);

        sun_rts.hour = 0; sun_rts.minute = 0; sun_rts.second = 0.0;
        sun_rts.delta_ut1 = 0.0; sun_rts.timezone = 0.0;

        sun_rts.jd = julian_day (sun_rts.year,   sun_rts.month,  sun_rts.day,       sun_rts.hour,
                                 sun_rts.minute, sun_rts.second, sun_rts.delta_ut1, sun_rts.timezone);

        sun_rts.calculate_geocentric_sun_right_ascension_and_declination();
        let nu: f64 = sun_rts.nu;

        sun_rts.delta_t = 0.0;
        sun_rts.jd -= 1.0;

        for i in 0..JD_COUNT {
            sun_rts.calculate_geocentric_sun_right_ascension_and_declination();
            alpha[i] = sun_rts.alpha;
            delta[i] = sun_rts.delta;
            sun_rts.jd += 1.0;
        }

        m_rts[SUN_TRANSIT] = approx_sun_transit_time(alpha[JD_ZERO], self.longitude, nu);
        let h0: f64 = sun_hour_angle_at_rise_set(self.latitude, delta[JD_ZERO], h0_prime);

        if h0 >= 0.0 {
            approx_sun_rise_and_set(&mut m_rts, h0);

            for i in 0..SUN_COUNT {
                nu_rts[i]      = nu + 360.985647*m_rts[i];

                let n: f64     = m_rts[i] + self.delta_t / 86400.0;
                alpha_prime[i] = rts_alpha_delta_prime(&alpha, n);
                delta_prime[i] = rts_alpha_delta_prime(&delta, n);

                h_prime[i]     = limit_degrees180pm(nu_rts[i] + self.longitude - alpha_prime[i]);

                h_rts[i]       = rts_sun_altitude(self.latitude, delta_prime[i], h_prime[i]);
            }

            self.srha = h_prime[SUN_RISE];
            self.ssha = h_prime[SUN_SET];
            self.sta  = h_rts[SUN_TRANSIT];

            self.suntransit = dayfrac_to_local_hr(m_rts[SUN_TRANSIT] - h_prime[SUN_TRANSIT] / 360.0, self.timezone);

            self.sunrise = dayfrac_to_local_hr(sun_rise_and_set(&m_rts, &h_rts, &delta_prime,
                                                                self.latitude, &h_prime, h0_prime, SUN_RISE), self.timezone);

            self.sunset  = dayfrac_to_local_hr(sun_rise_and_set(&m_rts, &h_rts, &delta_prime,
                                                                self.latitude, &h_prime, h0_prime, SUN_SET), self.timezone);

        } else {
            self.srha       = -99999.0;
            self.ssha       = -99999.0;
            self.sta        = -99999.0;
            self.suntransit = -99999.0;
            self.sunrise    = -99999.0;
            self.sunset     = -99999.0;
        }
    }

    /// Calculate required SPA parameters to get the right ascension (alpha) and declination (delta)
    /// Note: JD must be already calculated and in structure
    ///
    /// # Arguments
    ///
    /// * 'spa' - the `SpaData` struct
    fn calculate_geocentric_sun_right_ascension_and_declination(&mut self) {
        let mut x: [f64;TERM_X_COUNT] = [0.0; TERM_X_COUNT];

        self.jc = julian_century(self.jd);

        self.jde = julian_ephemeris_day(self.jd, self.delta_t);
        self.jce = julian_ephemeris_century(self.jde);
        self.jme = julian_ephemeris_millennium(self.jce);

        self.l = earth_heliocentric_longitude(self.jme);
        self.b = earth_heliocentric_latitude(self.jme);
        self.r = earth_radius_vector(self.jme);

        self.theta = geocentric_longitude(self.l);
        self.beta  = geocentric_latitude(self.b);

        x[TERM_X0] = { self.x0 = mean_elongation_moon_sun(self.jce); self.x0 };
        x[TERM_X1] = { self.x1 = mean_anomaly_sun(self.jce);         self.x1 };
        x[TERM_X2] = { self.x2 = mean_anomaly_moon(self.jce);        self.x2 };
        x[TERM_X3] = { self.x3 = argument_latitude_moon(self.jce);   self.x3 };
        x[TERM_X4] = { self.x4 = ascending_longitude_moon(self.jce); self.x4 };

        self.nutation_longitude_and_obliquity(&x);

        self.epsilon0 = ecliptic_mean_obliquity(self.jme);
        self.epsilon  = ecliptic_true_obliquity(self.del_epsilon, self.epsilon0);

        self.del_tau   = aberration_correction(self.r);
        self.lamda     = apparent_sun_longitude(self.theta, self.del_psi, self.del_tau);
        self.nu0       = greenwich_mean_sidereal_time(self.jd, self.jc);
        self.nu        = greenwich_sidereal_time(self.nu0, self.del_psi, self.epsilon);

        self.alpha = geocentric_right_ascension(self.lamda, self.epsilon, self.beta);
        self.delta = geocentric_declination(self.beta, self.epsilon, self.lamda);
    }

    fn nutation_longitude_and_obliquity(&mut self, x: &[f64]) {
        let mut sum_psi: f64 = 0.0;
        let mut sum_epsilon: f64 = 0.0;

        for i in 0..Y_COUNT {
            let xy_term_sum  = deg2rad(xy_term_summation(i, x));
            sum_psi     += (PE_TERMS[i][TERM_PSI_A] + self.jce * PE_TERMS[i][TERM_PSI_B]) * xy_term_sum.sin();
            sum_epsilon += (PE_TERMS[i][TERM_EPS_C] + self.jce * PE_TERMS[i][TERM_EPS_D]) * xy_term_sum.cos();
        }

        self.del_psi = sum_psi / 36000000.0;
        self.del_epsilon = sum_epsilon / 36000000.0;
    }
}

fn integer(value: f64) -> f64 {
    (value as i64) as f64
}

fn limit_degrees180pm(mut degrees: f64) -> f64 {
    degrees /= 360.0;
    let limited = 360.0 * (degrees - degrees.floor());

    if limited < -180.0 {
        limited + 360.0
    } else if limited > 180.0 {
        limited - 360.0
    } else {
        limited
    }
}

fn limit_degrees180(mut degrees: f64) -> f64 {
    degrees /= 180.0;
    let limited = 180.0 * (degrees - degrees.floor());

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

fn julian_day(mut year: i64, mut month: i64, day: i64, hour: i64, minute: i64, second: f64, dut1: f64, tz: f64) -> f64 {
    let day_decimal: f64 = day as f64 + (hour as f64 - tz + (minute as f64 + (second + dut1) / 60.0) / 60.0) / 24.0;

    if month < 3 {
        month += 12;
        year -= 1;
    }

    let mut julian_day: f64 = integer(365.25 * (year as f64 + 4716.0)) + integer(30.6001 * (month as f64 + 1.0)) + day_decimal - 1524.5;

    if julian_day > 2299160.0 {
        let a = integer(year as f64 / 100.0);
        julian_day += 2.0 - a + integer(a / 4.0);
    }

    julian_day
}

fn julian_century(jd: f64) -> f64 {
    (jd - 2451545.0) / 36525.0
}

fn julian_ephemeris_day(jd: f64, delta_t: f64) -> f64 {
    jd + delta_t / 86400.0
}

fn julian_ephemeris_century(jde: f64) -> f64 {
    (jde - 2451545.0) / 36525.0
}

fn julian_ephemeris_millennium(jce: f64) -> f64 {
    jce / 10.0
}

fn earth_periodic_term_summation(terms: &[[f64;TERM_COUNT]], count: i64, jme: f64) -> f64 {
    let mut sum:f64 = 0.0;
    for i in 0..count as usize {
        sum += terms[i][TERM_A] * (terms[i][TERM_B] + terms[i][TERM_C] * jme).cos();
    }

    sum
}

fn earth_values(term_sum: &[f64], count: usize, jme: f64) -> f64 {
    let mut sum:f64 = 0.0;
    for i in 0..count {
        sum += term_sum[i] * jme.powi(i as i32);
    }

    sum /= 1.0e8;

    sum
}

fn earth_heliocentric_longitude(jme: f64) -> f64 {
    let mut sum: [f64;L_COUNT] = [0.0; L_COUNT];
    for i in 0..L_COUNT {
        sum[i] = earth_periodic_term_summation(&L_TERMS[i], L_SUBCOUNT[i], jme);
    }

    limit_degrees(rad2deg(earth_values(&sum, L_COUNT, jme)))
}

fn earth_heliocentric_latitude(jme: f64) -> f64 {
    let mut sum: [f64;B_COUNT] = [0.0; B_COUNT];
    for i in 0..B_COUNT {
        sum[i] = earth_periodic_term_summation(&B_TERMS[i], B_SUBCOUNT[i], jme);
    }

    rad2deg(earth_values(&sum, B_COUNT, jme))
}

fn earth_radius_vector(jme: f64) -> f64 {
    let mut sum: [f64;R_COUNT] = [0.0; R_COUNT];
    for i in 0..R_COUNT {
        sum[i] = earth_periodic_term_summation(&R_TERMS[i], R_SUBCOUNT[i], jme);
    }

    earth_values(&sum, R_COUNT, jme)
}

fn geocentric_longitude(l: f64) -> f64 {
    let theta: f64 = l + 180.0;

    if theta >= 360.0 {
        theta - 360.0
    } else {
        theta
    }
}

fn geocentric_latitude(b: f64) -> f64 {
    -b
}

fn mean_elongation_moon_sun(jce: f64) -> f64 {
    third_order_polynomial(1.0 / 189474.0, -0.0019142, 445267.11148, 297.85036, jce)
}

fn mean_anomaly_sun(jce: f64) -> f64 {
    third_order_polynomial(-1.0 / 300000.0, -0.0001603, 35999.05034, 357.52772, jce)
}

fn mean_anomaly_moon(jce: f64) -> f64 {
    third_order_polynomial(1.0 / 56250.0, 0.0086972, 477198.867398, 134.96298, jce)
}

fn argument_latitude_moon(jce: f64) -> f64 {
    third_order_polynomial(1.0 / 327270.0, -0.0036825, 483202.017538, 93.27191, jce)
}

fn ascending_longitude_moon(jce: f64) -> f64 {
    third_order_polynomial(1.0 / 450000.0, 0.0020708, -1934.136261, 125.04452, jce)
}

fn xy_term_summation(i: usize, x: &[f64]) -> f64 {
    let mut sum:f64 = 0.0;
    for j in 0..TERM_Y_COUNT {
        sum += x[j] * Y_TERMS[i][j] as f64;
    }

    sum
}

fn ecliptic_mean_obliquity(jme: f64) -> f64 {
    let u: f64 = jme / 10.0;

    84381.448 + u * (-4680.93 + u * (-1.55 + u * (1999.25 + u * (-51.38 + u * (-249.67 +
        u * (  -39.05 + u * ( 7.12 + u * (  27.87 + u * (  5.79 + u * 2.45)))))))))
}

fn ecliptic_true_obliquity(delta_epsilon: f64, epsilon0: f64) -> f64 {
    delta_epsilon + epsilon0 / 3600.0
}

fn aberration_correction(r: f64) -> f64 {
    -20.4898 / (3600.0 * r)
}

fn apparent_sun_longitude(theta: f64, delta_psi: f64, delta_tau: f64) -> f64 {
    theta + delta_psi + delta_tau
}

fn greenwich_mean_sidereal_time(jd: f64, jc: f64) -> f64 {
    limit_degrees(280.46061837 + 360.98564736629 * (jd - 2451545.0) +
        jc * jc * (0.000387933 - jc/38710000.0))
}

fn greenwich_sidereal_time(nu0: f64, delta_psi: f64, epsilon: f64) -> f64 {
    nu0 + delta_psi * deg2rad(epsilon).cos()
}

fn sun_equatorial_horizontal_parallax(r: f64) -> f64 {
    8.794 / (3600.0 * r)
}

fn surface_incidence_angle(zenith: f64, azimuth_astro: f64, azm_rotation: f64, slope: f64) -> f64 {
    let zenith_rad: f64 = deg2rad(zenith);
    let slope_rad: f64  = deg2rad(slope);

    rad2deg((zenith_rad.cos() * slope_rad.cos()  +
        slope_rad.sin() * zenith_rad.sin() * deg2rad(azimuth_astro - azm_rotation).cos()).acos())
}

fn sun_mean_longitude(jme: f64) -> f64 {
    limit_degrees(280.4664567 + jme * (360007.6982779 + jme * (0.03032028 +
        jme * (1.0 / 49931.0 + jme * (-1.0 / 15300.0 + jme * (-1.0 / 2000000.0))))))
}

fn eot(m: f64, alpha: f64, del_psi: f64, epsilon: f64) -> f64 {
    limit_minutes(4.0 * (m - 0.0057183 - alpha + del_psi * deg2rad(epsilon).cos()))
}

fn approx_sun_transit_time(alpha_zero: f64, longitude: f64, nu: f64) -> f64 {
    (alpha_zero - longitude - nu) / 360.0
}

fn sun_hour_angle_at_rise_set(latitude: f64, delta_zero: f64, h0_prime: f64) -> f64 {
    let mut h0: f64             = -99999.0;
    let latitude_rad: f64   = deg2rad(latitude);
    let delta_zero_rad: f64 = deg2rad(delta_zero);
    let argument: f64       = (deg2rad(h0_prime).sin() - latitude_rad.sin() * delta_zero_rad.sin()) /
        (latitude_rad.cos() * delta_zero_rad.cos());

    if argument.abs() <= 1.0 {
        h0 = limit_degrees180(rad2deg(argument.acos()));
    }

    h0
}

fn approx_sun_rise_and_set(m_rts: &mut [f64], h0: f64) {
    let h0_dfrac: f64 = h0 / 360.0;

    m_rts[SUN_RISE]    = limit_zero2one(m_rts[SUN_TRANSIT] - h0_dfrac);
    m_rts[SUN_SET]     = limit_zero2one(m_rts[SUN_TRANSIT] + h0_dfrac);
    m_rts[SUN_TRANSIT] = limit_zero2one(m_rts[SUN_TRANSIT]);
}

fn rts_alpha_delta_prime(ad: &[f64], n: f64) -> f64 {
    let mut a: f64 = ad[JD_ZERO] - ad[JD_MINUS];
    let mut b: f64 = ad[JD_PLUS] - ad[JD_ZERO];

    if a.abs() >= 2.0 {
        a = limit_zero2one(a);
    }
    if b.abs() >= 2.0 {
        b = limit_zero2one(b);
    }

    ad[JD_ZERO] + n * (a + b + (b - a) * n) / 2.0
}

fn rts_sun_altitude(latitude: f64, delta_prime: f64, h_prime: f64) -> f64 {
    let latitude_rad: f64    = deg2rad(latitude);
    let delta_prime_rad: f64 = deg2rad(delta_prime);

    rad2deg((latitude_rad.sin() * delta_prime_rad.sin() +
        latitude_rad.cos() * delta_prime_rad.cos() * deg2rad(h_prime).cos()).asin())
}

fn sun_rise_and_set(m_rts: &[f64], h_rts: &[f64], delta_prime: &[f64], latitude: f64, h_prime: &[f64], h0_prime: f64, sun: usize) -> f64 {
    m_rts[sun] + (h_rts[sun] - h0_prime) /
        (360.0 * deg2rad(delta_prime[sun]).cos() * deg2rad(latitude).cos() * deg2rad(h_prime[sun]).sin())
}
