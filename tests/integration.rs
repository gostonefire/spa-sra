use spa_sra::spa::Function;
use spa_sra::SpaBuilder;

#[test]
fn builder_and_run() {
    let mut spa = SpaBuilder::new()
        .date(2003, 10, 17).unwrap_or_else(|e| panic!("{}", e))
        .time(12, 30, 30.0).unwrap_or_else(|e| panic!("{}", e))
        .timezone(-7.0).unwrap_or_else(|e| panic!("{}", e))
        .lat_long(39.742476, -105.1786).unwrap_or_else(|e| panic!("{}", e))
        .pressure_temperature(820.0, 11.0).unwrap_or_else(|e| panic!("{}", e))
        .atmospheric_refraction(0.5667).unwrap_or_else(|e| panic!("{}", e))
        .elevation(1830.14).unwrap_or_else(|e| panic!("{}", e))
        .slope(30.0).unwrap_or_else(|e| panic!("{}", e))
        .azimuth_rotation(-10.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_ut1(0.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_t(67.0).unwrap_or_else(|e| panic!("{}", e))
        .build(Function::SpaAll);


    // call the SPA calculate function and pass the SPA structure
    // test and test results according original code

    if let Err(e) = spa.spa_calculate() {
        panic!("{}", e);
    }

    assert_eq!(format!("{:.6}", spa.spa_za.zenith),        "50.111622",         "Zenith:        50.111622 degrees");
    assert_eq!(format!("{:.6}", spa.spa_za.azimuth),       "194.340241",       "Azimuth:       194.340241 degrees");
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

#[test]
fn build_and_calculate() {
    let spa = SpaBuilder::new()
        .date(2003, 10, 17).unwrap_or_else(|e| panic!("{}", e))
        .time(12, 30, 30.0).unwrap_or_else(|e| panic!("{}", e))
        .timezone(-7.0).unwrap_or_else(|e| panic!("{}", e))
        .lat_long(39.742476, -105.1786).unwrap_or_else(|e| panic!("{}", e))
        .pressure_temperature(820.0, 11.0).unwrap_or_else(|e| panic!("{}", e))
        .atmospheric_refraction(0.5667).unwrap_or_else(|e| panic!("{}", e))
        .elevation(1830.14).unwrap_or_else(|e| panic!("{}", e))
        .slope(30.0).unwrap_or_else(|e| panic!("{}", e))
        .azimuth_rotation(-10.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_ut1(0.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_t(67.0).unwrap_or_else(|e| panic!("{}", e))
        .build_and_calculate(Function::SpaAll).unwrap_or_else(|e| panic!("{}", e));

    // test and test results according original code

    assert_eq!(format!("{:.6}", spa.spa_za.zenith),        "50.111622",         "Zenith:        50.111622 degrees");
    assert_eq!(format!("{:.6}", spa.spa_za.azimuth),       "194.340241",       "Azimuth:       194.340241 degrees");
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

#[cfg(feature = "chrono_0_4")]
#[test]
fn feature_chrono_0_4() {
    use chrono::{FixedOffset, TimeZone, Timelike};

    let date_time = FixedOffset::east_opt(-7 * 3600)
        .unwrap()
        .with_ymd_and_hms(2003, 10, 17, 12, 30, 30)
        .unwrap();

    let spa = SpaBuilder::new()
        .date_time(date_time).unwrap_or_else(|e| panic!("{}", e))
        .lat_long(39.742476, -105.1786).unwrap_or_else(|e| panic!("{}", e))
        .pressure_temperature(820.0, 11.0).unwrap_or_else(|e| panic!("{}", e))
        .atmospheric_refraction(0.5667).unwrap_or_else(|e| panic!("{}", e))
        .elevation(1830.14).unwrap_or_else(|e| panic!("{}", e))
        .slope(30.0).unwrap_or_else(|e| panic!("{}", e))
        .azimuth_rotation(-10.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_ut1(0.0).unwrap_or_else(|e| panic!("{}", e))
        .delta_t(67.0).unwrap_or_else(|e| panic!("{}", e))
        .build_and_calculate(Function::SpaAll).unwrap_or_else(|e| panic!("{}", e));

    // test and test results according original code

    assert_eq!(format!("{:.6}", spa.spa_za.zenith),        "50.111622",         "Zenith:        50.111622 degrees");
    assert_eq!(format!("{:.6}", spa.spa_za.azimuth),       "194.340241",       "Azimuth:       194.340241 degrees");
    assert_eq!(format!("{:.6}", spa.spa_za_inc.incidence), "25.187000",      "Incidence:     25.187000 degrees");

    let mut min: f64 = 60.0 * (spa.spa_za_rts.sunrise - (spa.spa_za_rts.sunrise as i64) as f64);
    let mut sec: f64 = 60.0 * (min - (min as i64) as f64);
    assert_eq!(format!("{:0>2}:{:0>2}:{:0>2}", spa.spa_za_rts.sunrise as i64, min as i64, sec as i64), "06:12:43",
               "Sunrise: 06:12:43 Local Time");

    min = 60.0 * (spa.spa_za_rts.sunset - (spa.spa_za_rts.sunset as i64) as f64);
    sec = 60.0 * (min - (min as i64) as f64);
    assert_eq!(format!("{:0>2}:{:0>2}:{:0>2}", spa.spa_za_rts.sunset as i64, min as i64, sec as i64), "17:20:19",
               "Sunset: 17:20:19 Local Time");

    let was = spa.get_sunrise(date_time.timezone());
    let should = FixedOffset::east_opt(-7 * 3600)
        .unwrap()
        .with_ymd_and_hms(2003, 10, 17, 06, 12, 43)
        .unwrap()
        .with_nanosecond(439793424)
        .unwrap();

    assert_eq!(was, should);
}