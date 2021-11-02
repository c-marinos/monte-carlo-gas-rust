use rand::Rng;

pub mod settings;
pub mod util;

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
        let mut this_energy: f64 = 0.0;
        let mut tot_energy: f64 = 0.0;

        for i in 0..self.n
        {
            let coord_one: (f64,f64) = self.coords[i as usize];

            for j in 0..self.n
            {
                let coord_two: (f64,f64) = self.coords[j as usize];

                if(coord_one != coord_two && j > i)
                {
                    let d = util::calc_dis(coord_two, coord_one);
                    if(d > self.settings.sigma)
                    {
                        this_energy = 0.0;
                    }
                    else
                    {
                        this_energy = std::f64::INFINITY;
                    }
                }
            }
        }
        this_energy
    }

    pub fn square_well_energy_seed()
    {

    }

    pub fn lennard_jones_energy_seed()
    {

    }

    pub fn hard_sphere_energy()
    {

    }

    pub fn square_well_energy()
    {

    }

    pub fn lennard_jones_energy()
    {

    }
}
