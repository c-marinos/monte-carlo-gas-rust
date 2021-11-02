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

pub fn avg_dis(coords: Vec<(f64,f64)>) -> f64
{
    let mut avg_d:f64 = 0.0;
    let n = coords.len();

    for i in 0..n
    {
        for j in 0..n
        {
            if i != j && j > i
            {
                avg_d += calc_dis(coords[i], coords[j]);
            }
        }
    }

    avg_d / (n as f64 * (n as f64 - 1.0) / 2.0)
}
