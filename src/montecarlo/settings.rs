
// Current state of UI buttons in window
pub struct Settings
{
    pub potential_function: PotentialFunction,
    pub temp: Temperature,
    pub r: f64,
    pub sigma: f64,
    pub epsilon: f64,
}

impl Settings
{
    pub fn new(potential_function: PotentialFunction, temp: Temperature, r:f64, sigma:f64, epsilon: f64) -> Settings
    {

        Settings
        {
            potential_function: potential_function,
            temp: temp,
            r: r,
            sigma: sigma,
            epsilon: epsilon,
        }
    }
}

pub enum Temperature
{
    OneHundred,
    ThreeHundred,
    FiveHundred
}

pub enum PotentialFunction
{
    HardSphere,
    SquareWell,
    LennardJones
}
