
// Current state of UI buttons in window
pub struct Settings
{
    pub potential_function: PotentialFunction,
    pub temp: Temperature,
    pub r: f64,
}

impl Settings
{
    pub fn new()
    {

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
