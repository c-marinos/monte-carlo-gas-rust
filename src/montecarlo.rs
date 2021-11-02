use rand::Rng;

pub mod settings;
pub mod util;

static KB:f64 = 0.0019872036;
pub const E: f64 = 2.71828182845904523536028747135266250f64;

pub struct State {
    pub n: i32,
    pub settings: settings::Settings,
    pub coords: Vec<(f64,f64)>,
    pub l: f32, // length of box
}

impl State
{
    pub fn new(n: i32, settings: settings::Settings, l: f32) -> State
    {
        let coord = (0.0, 0.0);
        let mut init_coords = Vec::<(f64, f64)>::new();
        for i in 0..n{
            init_coords.push(coord);
        }

        State {
            n: n,
            settings: settings,
            coords: init_coords,
            l: l,
        }
    }

    pub fn evolve()
    {

    }

    pub fn hard_sphere_energy_seed(&self) -> f64
    {
        let mut this_energy: f64;
        let mut tot_energy: f64 = 0.0;

        for i in 0..self.n
        {
            let coord_one: (f64,f64) = self.coords[i as usize];

            for j in 0..self.n
            {
                let coord_two: (f64,f64) = self.coords[j as usize];

                if coord_one != coord_two && j > i
                {
                    let d = util::calc_dis(coord_two, coord_one);

                    if d >= self.settings.sigma
                    {
                        this_energy = 0.0;
                    }
                    else
                    {
                        this_energy = std::f64::INFINITY;
                    }
                    tot_energy += this_energy
                }
            }
        }
        tot_energy
    }

    pub fn square_well_energy_seed(&self) -> f64
    {
        let mut this_energy: f64 = 0.0;
        let mut tot_energy: f64 = 0.0;

        for i in 0..self.n
        {
            let coord_one: (f64,f64) = self.coords[i as usize];

            for j in 0..self.n
            {
                let coord_two: (f64,f64) = self.coords[j as usize];

                if coord_one != coord_two && j > i
                {
                    let d = util::calc_dis(coord_two, coord_one);

                    if d < self.settings.sigma
                    {
                        this_energy = std::f64::INFINITY;
                    }
                    else if d >= self.settings.sigma && d <= self.settings.r * self.settings.sigma
                    {
                        this_energy = -1 as f64 * self.settings.epsilon;
                    }
                    else if d > self.settings.r * self.settings.sigma
                    {
                        this_energy = 0.0;
                    }
                    tot_energy += this_energy
                }
            }
        }
        tot_energy
    }

    pub fn lennard_jones_energy_seed(&self) -> f64
    {
        let mut this_energy: f64;
        let mut tot_energy: f64 = 0.0;

        for i in 0..self.n
        {
            let coord_one: (f64,f64) = self.coords[i as usize];

            for j in 0..self.n
            {
                let coord_two: (f64,f64) = self.coords[j as usize];

                if coord_one != coord_two && j > i
                {
                    let d = util::calc_dis(coord_two, coord_one);
                    this_energy = 4.0 * self.settings.epsilon * f64::powf(self.settings.sigma / d , 12.0) - f64::powf(self.settings.sigma / d, 6.0);
                    tot_energy += this_energy
                }
            }
        }
        tot_energy
    }

    pub fn hard_sphere_energy(&self, coord: (f64, f64)) -> f64
    {
        let mut this_energy: f64;
        let mut tot_energy: f64 = 0.0;

        for i in 0..self.n
        {
            let coord_one: (f64,f64) = self.coords[i as usize];

            if coord_one != coord
            {
                let d = util::calc_dis(coord, coord_one);

                if d >= self.settings.sigma
                {
                    this_energy = 0.0;
                }
                else
                {
                    this_energy = std::f64::INFINITY;
                }
                tot_energy += this_energy
            }
        }

        tot_energy
    }

    pub fn square_well_energy(&self, coord: (f64, f64)) -> f64
    {
        let mut this_energy: f64 = 0.0;
        let mut tot_energy: f64 = 0.0;

        for i in 0..self.n
        {
            let coord_one: (f64,f64) = self.coords[i as usize];

            if coord_one != coord
            {
                let d = util::calc_dis(coord, coord_one);

                if d < self.settings.sigma
                {
                    this_energy = std::f64::INFINITY;
                }
                else if d >= self.settings.sigma && d <= self.settings.r * self.settings.sigma
                {
                    this_energy = -1 as f64 * self.settings.epsilon;
                }
                else if d > self.settings.r * self.settings.sigma
                {
                    this_energy = 0.0;
                }
                tot_energy += this_energy
            }
        }

        tot_energy
    }

    pub fn lennard_jones_energy(&self, coord: (f64, f64)) -> f64
    {
        let mut this_energy: f64;
        let mut tot_energy: f64 = 0.0;

        for i in 0..self.n
        {
            let coord_one: (f64,f64) = self.coords[i as usize];

            if coord_one != coord
            {
                let d = util::calc_dis(coord, coord_one);
                this_energy = 4.0 * self.settings.epsilon * f64::powf(self.settings.sigma / d, 12.0) - f64::powf(self.settings.sigma / d, 6.0);
                tot_energy += this_energy
            }
        }

        tot_energy
    }

    pub fn probability(dE: f64, temperature: settings::Temperature) -> f64
    {
        let mut temp:i32 = 0;
        match temperature {
            settings::Temperature::OneHundred => temp = 100,
            settings::Temperature::ThreeHundred => temp = 300,
            settings::Temperature::FiveHundred => temp = 500,
            _ => println!("Error with temperature in probability function."),
        }

        let  kbt:f64 = KB * temp as f64;

        f64::powf(E, -1 as f64 * dE / kbt)
    }
}
