use crate::constants::{b_subcount, l_subcount, r_subcount, B_COUNT, JD_MINUS, JD_PLUS, JD_ZERO, L_COUNT, PI, R_COUNT, SUN_RADIUS, SUN_RISE, SUN_SET, SUN_TRANSIT, TERM_A, TERM_B, TERM_C, TERM_COUNT, TERM_EPS_C, TERM_EPS_D, TERM_PSI_A, TERM_PSI_B, TERM_Y_COUNT, Y_COUNT};
use crate::earth_periodic_terms::{B_TERMS, L_TERMS, R_TERMS};
use crate::nutation_obliquity_periodic_terms::{PE_TERMS, Y_TERMS};
use crate::spa::{Output, SpaData};

fn rad2deg(radians: f64) -> f64 {
    180.0/PI * radians
}

fn deg2rad(degrees: f64) -> f64 {
    PI/180.0 * degrees
}

fn integer(value: f64) -> f64 {
    (value as i64) as f64
}

fn limit_degrees(mut degrees: f64) -> f64 {
    degrees /= 360.0;
    let limited = 360.0 * (degrees - degrees.floor());
    
    if limited < 0.0 {
        limited + 360.0
    } else {
        limited
    }
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

fn third_order_polynomial(a: f64, b: f64, c: f64, d: f64, x: f64) -> f64 {
    ((a * x + b) * x + c) * x + d
}


/// Validates inputs
/// 
/// # Arguments
/// 
/// * 'spa' - the `SpaData` struct
fn validate_inputs(spa: &SpaData) -> i64 {
    if spa.year        < -2000   || spa.year        > 6000   { return 1 };
    if spa.month       < 1       || spa.month       > 12     { return 2 };
    if spa.day         < 1       || spa.day         > 31     { return 3 };
    if spa.hour        < 0       || spa.hour        > 24     { return 4 };
    if spa.minute      < 0       || spa.minute      > 59     { return 5 };
    if spa.second      < 0.0     || spa.second      >=60.0   { return 6 };
    if spa.pressure    < 0.0     || spa.pressure    > 5000.0 { return 12 };
    if spa.temperature <= -273.0 || spa.temperature > 6000.0 { return 13 };
    if spa.delta_ut1   <= -1.0   || spa.delta_ut1   >= 1.0   { return 17 };
    if spa.hour        == 24     && spa.minute      > 0      { return 5 };
    if spa.hour        == 24     && spa.second      > 0.0    { return 6 };
    
    if spa.delta_t.abs()       > 8000.0     { return 7 };
    if spa.timezone.abs()      > 18.0       { return 8 };
    if spa.longitude.abs()     > 180.0      { return 9 };
    if spa.latitude.abs()      > 90.0       { return 10 };
    if spa.atmos_refract.abs() > 5.0        { return 16 };
    if spa.elevation           < -6500000.0 { return 11 };

    if spa.function == Output::SpaZaInc || spa.function == Output::SpaAll {
        if spa.slope.abs()  > 360.0 { return 14 };
        if spa.azm_rotation > 360.0 { return 15 };
    }
    
    0
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
        sum[i] = earth_periodic_term_summation(&L_TERMS[i], l_subcount[i], jme);
    }

    limit_degrees(rad2deg(earth_values(&sum, L_COUNT, jme)))
}

fn earth_heliocentric_latitude(jme: f64) -> f64 {
    let mut sum: [f64;B_COUNT] = [0.0; B_COUNT];
    for i in 0..B_COUNT {
        sum[i] = earth_periodic_term_summation(&B_TERMS[i], b_subcount[i], jme);
    }

    rad2deg(earth_values(&sum, B_COUNT, jme))
}

fn earth_radius_vector(jme: f64) -> f64 {
    let mut sum: [f64;R_COUNT] = [0.0; R_COUNT];
    for i in 0..R_COUNT {
        sum[i] = earth_periodic_term_summation(&R_TERMS[i], r_subcount[i], jme);
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

fn nutation_longitude_and_obliquity(jce: f64, x: &[f64], del_psi: &mut f64, del_epsilon: &mut f64) {
    let mut sum_psi: f64 = 0.0;
    let mut sum_epsilon: f64 = 0.0;
    
    for i in 0..Y_COUNT {
        let xy_term_sum  = deg2rad(xy_term_summation(i, x));
        sum_psi     += (PE_TERMS[i][TERM_PSI_A] + jce*PE_TERMS[i][TERM_PSI_B]) * xy_term_sum.sin();
        sum_epsilon += (PE_TERMS[i][TERM_EPS_C] + jce*PE_TERMS[i][TERM_EPS_D]) * xy_term_sum.cos();
    }
    
    *del_psi = sum_psi / 36000000.0;
    *del_epsilon = sum_epsilon / 36000000.0;
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

fn geocentric_right_ascension(lamda: f64, epsilon: f64, beta: f64) -> f64 {
    let lamda_rad: f64   = deg2rad(lamda);
    let epsilon_rad: f64 = deg2rad(epsilon);

    limit_degrees(rad2deg((lamda_rad.sin() * epsilon_rad.cos() - 
        deg2rad(beta).tan() * epsilon_rad.sin()).atan2(lamda_rad.cos())))
}

fn geocentric_declination(beta: f64, epsilon: f64, lamda: f64) -> f64 {
    let beta_rad: f64    = deg2rad(beta);
    let epsilon_rad: f64 = deg2rad(epsilon);

    rad2deg((beta_rad.sin() * epsilon_rad.cos() + 
        beta_rad.cos() * epsilon_rad.sin() * deg2rad(lamda).sin()).asin())
}

fn observer_hour_angle(nu: f64, longitude: f64, alpha_deg: f64) -> f64 {
    limit_degrees(nu + longitude - alpha_deg)
}

fn sun_equatorial_horizontal_parallax(r: f64) -> f64 {
    8.794 / (3600.0 * r)
}

fn right_ascension_parallax_and_topocentric_dec(latitude: f64, elevation: f64, xi: f64, h: f64, delta: f64, delta_alpha: &mut f64, delta_prime: &mut f64) {
    let lat_rad: f64   = deg2rad(latitude);
    let xi_rad: f64    = deg2rad(xi);
    let h_rad: f64     = deg2rad(h);
    let delta_rad: f64 = deg2rad(delta);
    let u: f64 = (0.99664719 * lat_rad.tan()).atan();
    let y: f64 = 0.99664719 * u.sin() + elevation * lat_rad.sin() / 6378140.0;
    let x: f64 =              u.cos() + elevation * lat_rad.cos() / 6378140.0;

    let delta_alpha_rad: f64 = (             - x * xi_rad.sin()  * h_rad.sin())
        .atan2(        delta_rad.cos() - x * xi_rad.sin()  * h_rad.cos());

    *delta_prime = rad2deg(((delta_rad.sin() - y * xi_rad.sin()) * delta_alpha_rad.cos())
        .atan2(        delta_rad.cos() - x * xi_rad.sin()  * h_rad.cos()));

    *delta_alpha = rad2deg(delta_alpha_rad);
}

fn topocentric_right_ascension(alpha_deg: f64, delta_alpha: f64) -> f64 {
    alpha_deg + delta_alpha
}

fn topocentric_local_hour_angle(h: f64, delta_alpha: f64) -> f64 {
    h - delta_alpha
}

fn topocentric_elevation_angle(latitude: f64, delta_prime: f64, h_prime: f64) -> f64 {
    let lat_rad: f64         = deg2rad(latitude);
    let delta_prime_rad: f64 = deg2rad(delta_prime);

    rad2deg((lat_rad.sin() * delta_prime_rad.sin() + 
             lat_rad.cos() * delta_prime_rad.cos() * deg2rad(h_prime).cos()).asin())
}

fn atmospheric_refraction_correction(pressure: f64, temperature: f64, atmos_refract: f64, e0: f64) -> f64 {
    let mut del_e: f64 = 0.0;

    if e0 >= -1.0 * (SUN_RADIUS + atmos_refract) {
        del_e = (pressure / 1010.0) * (283.0 / (273.0 + temperature)) * 
            1.02 / (60.0 * deg2rad(e0 + 10.3/(e0 + 5.11)).tan());
    }

    del_e
}

fn topocentric_elevation_angle_corrected(e0: f64, delta_e: f64) -> f64 {
    e0 + delta_e
}

fn topocentric_zenith_angle(e: f64) -> f64 {
    90.0 - e
}

fn topocentric_azimuth_angle_astro(h_prime: f64, latitude: f64, delta_prime: f64) -> f64 {
    let h_prime_rad: f64 = deg2rad(h_prime);
    let lat_rad: f64     = deg2rad(latitude);

    limit_degrees(rad2deg(h_prime_rad.sin().atan2(h_prime_rad.cos() * lat_rad.sin() - deg2rad(delta_prime).tan() * lat_rad.cos())))
}

fn topocentric_azimuth_angle(azimuth_astro: f64) -> f64 {
    limit_degrees(azimuth_astro + 180.0)
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

/*

////////////////////////////////////////////////////////////////////////////////////////////////
// Calculate required SPA parameters to get the right ascension (alpha) and declination (delta)
// Note: JD must be already calculated and in structure
////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_geocentric_sun_right_ascension_and_declination(spa_data *spa)
{
    double x[TERM_X_COUNT];

    spa->jc = julian_century(spa->jd);

    spa->jde = julian_ephemeris_day(spa->jd, spa->delta_t);
    spa->jce = julian_ephemeris_century(spa->jde);
    spa->jme = julian_ephemeris_millennium(spa->jce);

    spa->l = earth_heliocentric_longitude(spa->jme);
    spa->b = earth_heliocentric_latitude(spa->jme);
    spa->r = earth_radius_vector(spa->jme);

    spa->theta = geocentric_longitude(spa->l);
    spa->beta  = geocentric_latitude(spa->b);

    x[TERM_X0] = spa->x0 = mean_elongation_moon_sun(spa->jce);
    x[TERM_X1] = spa->x1 = mean_anomaly_sun(spa->jce);
    x[TERM_X2] = spa->x2 = mean_anomaly_moon(spa->jce);
    x[TERM_X3] = spa->x3 = argument_latitude_moon(spa->jce);
    x[TERM_X4] = spa->x4 = ascending_longitude_moon(spa->jce);

    nutation_longitude_and_obliquity(spa->jce, x, &(spa->del_psi), &(spa->del_epsilon));

    spa->epsilon0 = ecliptic_mean_obliquity(spa->jme);
    spa->epsilon  = ecliptic_true_obliquity(spa->del_epsilon, spa->epsilon0);

    spa->del_tau   = aberration_correction(spa->r);
    spa->lamda     = apparent_sun_longitude(spa->theta, spa->del_psi, spa->del_tau);
    spa->nu0       = greenwich_mean_sidereal_time (spa->jd, spa->jc);
    spa->nu        = greenwich_sidereal_time (spa->nu0, spa->del_psi, spa->epsilon);

    spa->alpha = geocentric_right_ascension(spa->lamda, spa->epsilon, spa->beta);
    spa->delta = geocentric_declination(spa->beta, spa->epsilon, spa->lamda);
}

////////////////////////////////////////////////////////////////////////
// Calculate Equation of Time (EOT) and Sun Rise, Transit, & Set (RTS)
////////////////////////////////////////////////////////////////////////

void calculate_eot_and_sun_rise_transit_set(spa_data *spa)
{
    spa_data sun_rts;
    double nu, m, h0, n;
    double alpha[JD_COUNT], delta[JD_COUNT];
    double m_rts[SUN_COUNT], nu_rts[SUN_COUNT], h_rts[SUN_COUNT];
    double alpha_prime[SUN_COUNT], delta_prime[SUN_COUNT], h_prime[SUN_COUNT];
    double h0_prime = -1*(SUN_RADIUS + spa->atmos_refract);
    int i;

	sun_rts  = *spa;
    m        = sun_mean_longitude(spa->jme);
    spa->eot = eot(m, spa->alpha, spa->del_psi, spa->epsilon);

    sun_rts.hour = sun_rts.minute = sun_rts.second = 0;
	sun_rts.delta_ut1 = sun_rts.timezone = 0.0;

    sun_rts.jd = julian_day (sun_rts.year,   sun_rts.month,  sun_rts.day,       sun_rts.hour,
		                     sun_rts.minute, sun_rts.second, sun_rts.delta_ut1, sun_rts.timezone);

    calculate_geocentric_sun_right_ascension_and_declination(&sun_rts);
    nu = sun_rts.nu;

    sun_rts.delta_t = 0;
    sun_rts.jd--;
    for (i = 0; i < JD_COUNT; i++) {
        calculate_geocentric_sun_right_ascension_and_declination(&sun_rts);
        alpha[i] = sun_rts.alpha;
        delta[i] = sun_rts.delta;
        sun_rts.jd++;
    }

    m_rts[SUN_TRANSIT] = approx_sun_transit_time(alpha[JD_ZERO], spa->longitude, nu);
    h0 = sun_hour_angle_at_rise_set(spa->latitude, delta[JD_ZERO], h0_prime);

    if (h0 >= 0) {

        approx_sun_rise_and_set(m_rts, h0);

        for (i = 0; i < SUN_COUNT; i++) {

            nu_rts[i]      = nu + 360.985647*m_rts[i];

            n              = m_rts[i] + spa->delta_t/86400.0;
            alpha_prime[i] = rts_alpha_delta_prime(alpha, n);
            delta_prime[i] = rts_alpha_delta_prime(delta, n);

            h_prime[i]     = limit_degrees180pm(nu_rts[i] + spa->longitude - alpha_prime[i]);

            h_rts[i]       = rts_sun_altitude(spa->latitude, delta_prime[i], h_prime[i]);
        }

        spa->srha = h_prime[SUN_RISE];
        spa->ssha = h_prime[SUN_SET];
        spa->sta  = h_rts[SUN_TRANSIT];

        spa->suntransit = dayfrac_to_local_hr(m_rts[SUN_TRANSIT] - h_prime[SUN_TRANSIT] / 360.0,
                                              spa->timezone);

        spa->sunrise = dayfrac_to_local_hr(sun_rise_and_set(m_rts, h_rts, delta_prime,
                          spa->latitude, h_prime, h0_prime, SUN_RISE), spa->timezone);

        spa->sunset  = dayfrac_to_local_hr(sun_rise_and_set(m_rts, h_rts, delta_prime,
                          spa->latitude, h_prime, h0_prime, SUN_SET),  spa->timezone);

    } else spa->srha= spa->ssha= spa->sta= spa->suntransit= spa->sunrise= spa->sunset= -99999;

}

///////////////////////////////////////////////////////////////////////////////////////////
// Calculate all SPA parameters and put into structure
// Note: All inputs values (listed in header file) must already be in structure
///////////////////////////////////////////////////////////////////////////////////////////
int spa_calculate(spa_data *spa)
{
    int result;

    result = validate_inputs(spa);

    if (result == 0)
    {
        spa->jd = julian_day (spa->year,   spa->month,  spa->day,       spa->hour,
			                  spa->minute, spa->second, spa->delta_ut1, spa->timezone);

        calculate_geocentric_sun_right_ascension_and_declination(spa);

        spa->h  = observer_hour_angle(spa->nu, spa->longitude, spa->alpha);
        spa->xi = sun_equatorial_horizontal_parallax(spa->r);

        right_ascension_parallax_and_topocentric_dec(spa->latitude, spa->elevation, spa->xi,
                                spa->h, spa->delta, &(spa->del_alpha), &(spa->delta_prime));

        spa->alpha_prime = topocentric_right_ascension(spa->alpha, spa->del_alpha);
        spa->h_prime     = topocentric_local_hour_angle(spa->h, spa->del_alpha);

        spa->e0      = topocentric_elevation_angle(spa->latitude, spa->delta_prime, spa->h_prime);
        spa->del_e   = atmospheric_refraction_correction(spa->pressure, spa->temperature,
                                                         spa->atmos_refract, spa->e0);
        spa->e       = topocentric_elevation_angle_corrected(spa->e0, spa->del_e);

        spa->zenith        = topocentric_zenith_angle(spa->e);
        spa->azimuth_astro = topocentric_azimuth_angle_astro(spa->h_prime, spa->latitude,
                                                                           spa->delta_prime);
        spa->azimuth       = topocentric_azimuth_angle(spa->azimuth_astro);

        if ((spa->function == SPA_ZA_INC) || (spa->function == SPA_ALL))
            spa->incidence  = surface_incidence_angle(spa->zenith, spa->azimuth_astro,
                                                      spa->azm_rotation, spa->slope);

        if ((spa->function == SPA_ZA_RTS) || (spa->function == SPA_ALL))
            calculate_eot_and_sun_rise_transit_set(spa);
    }

    return result;
}
 */