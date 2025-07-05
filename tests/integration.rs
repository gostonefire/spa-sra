use spa_sra::spa::Function;
use spa_sra::{spa_calc, SpaBuilder};

#[test]
fn builder_and_run() {
    let mut spa = SpaBuilder::new()
        .date_time(2003, 10, 17, 12, 30, 30.0)
        .timezone(-7.0)
        .lat_long(39.742476, -105.1786)
        .pressure_temperature(820.0, 11.0)
        .atmospheric_refraction(0.5667)
        .elevation(1830.14)
        .slope(30.0)
        .azimuth_rotation(-10.0)
        .delta_ut1(0.0)
        .delta_t(67.0)
        .build(Function::SpaAll);


    // call the SPA calculate function and pass the SPA structure
    // test and test results according original code

    let result = spa_calc(&mut spa);
    assert_eq!(result, 0, "SPA Error Code: 0");

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
