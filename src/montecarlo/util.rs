use std::num;

pub struct Point2 {
    pub x: f32,
    pub y: f32,
}

pub fn calc_dis(first_point:(f64,f64), second_point:(f64,f64)) -> f64
{
    let x1 = first_point.0;
    let y1 = first_point.1;

    let x2 = second_point.0;
    let y2 = second_point.1;

    let xds = f64::powf(x2 - x1, 2.0);
    let yds = f64::powf(y2 - y1, 2.0);

    (xds + yds).sqrt()
}

pub fn avg_dis()
{

}

pub fn p()
{

}
