use rand::Rng;

pub mod settings;
pub mod util;

pub struct State {
    pub n: i32,
    pub settings: settings::Settings,
    pub coords: Vec<Vec<[f64;2]>>,
    pub l: f32, // length of box

}

impl State
{
    pub fn new(n:i32, l:f32)
    {
        let mut vec: Vec<Vec<[f64;2]>> = Vec::with_capacity(n as usize * 2 as usize);
        for i in 0..n{
            vec.push(vec!());
            for j in 0..2{
                vec[i as usize].push([0.0,0.0]);
            }
        }
        //Self{coords,..} = vec;
    }
}

pub fn evolve()
{

}

pub fn hard_sphere_energy_seed()
{

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
