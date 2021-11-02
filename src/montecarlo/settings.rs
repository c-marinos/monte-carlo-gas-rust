
// Current state of UI buttons in window
pub struct Settings
{
    pub potential_function: PotentialFunction,
    pub temp: Temperature,
    pub r: i32,
    pub sigma: f64,
}

impl Settings
{
    pub fn new(potential_function: PotentialFunction, temp: Temperature, r:i32, sigma:f64) -> Settings
    {

        Settings
        {
            potential_function: potential_function,
            temp: temp,
            r: r,
            sigma: sigma,
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
