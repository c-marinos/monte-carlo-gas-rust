use rand::Rng;

pub mod settings;
pub mod util;

static KB:f64 = 0.0019872036;
pub const E: f64 = 2.71828182845904523536028747135266250f64;

pub struct State {
    pub n: i32,
    pub settings: settings::Settings,
    pub coords: Vec<(f64,f64)>,
    pub l: f64, // length of box
    pub energies: Energies,
}

pub struct Energies
{
    pub hs_energy: f64,
    pub old_try_hs_energy: f64,
    pub try_hs_energy: f64,
    pub sw_energy: f64,
    pub old_try_sw_energy: f64,
    pub try_sw_energy: f64,
    pub lj_energy: f64,
    pub old_try_lj_energy: f64,
    pub try_lj_energy: f64,
    pub hs_dE: f64,
    pub sw_dE: f64,
    pub lj_dE: f64,
    pub hs_dE_list: Vec<f64>,
    pub sw_dE_list: Vec<f64>,
    pub lj_dE_list: Vec<f64>,
}

impl Energies
{
    pub fn new() -> Energies
    {
        Energies {
            hs_energy: 0.0,
            old_try_hs_energy: 0.0,
            try_hs_energy: 0.0,
            sw_energy: 0.0,
            old_try_sw_energy: 0.0,
            try_sw_energy: 0.0,
            lj_energy: 0.0,
            old_try_lj_energy: 0.0,
            try_lj_energy: 0.0,
            hs_dE: 0.0,
            sw_dE: 0.0,
            lj_dE: 0.0,
            hs_dE_list: Vec::<f64>::new(),
            sw_dE_list: Vec::<f64>::new(),
            lj_dE_list: Vec::<f64>::new(),
        }
    }
}

impl State
{
    pub fn new(n: i32, settings: settings::Settings, l: f64) -> State
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
            energies: Energies::new(),
        }
    }

    pub fn evolve(&self)
    {
        let mut this_try = Vec::<(f64, f64)>::new();
        let mut old_try = Vec::<(f64, f64)>::new();
        let randno: i32 = rand::thread_rng().gen_range(0,self.n);
        let e: (f64, f64) = self.coords[randno as usize];

        let mut changed_atom: (f64, f64) = (0.0, 0.0);
        let mut old_energy: (f64, f64) = (0.0, 0.0);

        for i in 0..self.n
        {
            let mut d = self.coords[i as usize];
            old_try.push(d);

            if i == randno
            {
                let x = e.0;
                let y = e.1;

                old_energy = (x, y);

                let xr = rand::thread_rng().gen_range(-1.0 * self.settings.step_size, self.settings.step_size);
                let yr = rand::thread_rng().gen_range(-1.0 * self.settings.step_size, self.settings.step_size);

                let mut new_x = x + xr;

                if new_x > self.l / 2.0 || new_x < -1.0 * self.l / 2.0
                {
                    new_x = -1.0*x + xr;
                }

                let mut new_y = y + yr;

                if new_y > self.l / 2.0 || new_y < -1.0 * self.l / 2.0
                {
                    new_y = -1.0*y + yr;
                }

                d = (new_x, new_y);
                changed_atom = (new_x, new_y);

            }

            this_try.push(d);
        }

        /*

		// we start with this bool set to true. if our energy goes up and our P() < alpha then this will be changed to false
		bool validmove = true;
		float alpha = UnityEngine.Random.Range(0f, 1f);

		if (hsb)
		{
			oldtryhsenergy = HardSphereEnergy(oldtry, oldenergy);
			tryhsenergy = HardSphereEnergy(thistry, changedatom);
			if (hsenergy - oldtryhsenergy + tryhsenergy <= hsenergy)
				goto Accepted;
			HSdE = -oldtryhsenergy + tryhsenergy;
			if (P(HSdE) < alpha)
				validmove = false;
		}
		else if (swb)
		{
			oldtryswenergy = SquareWellEnergy(oldtry, oldenergy);
			tryswenergy = SquareWellEnergy(thistry, changedatom);
			if (swenergy - oldtryswenergy + tryswenergy <= swenergy)
				goto Accepted;
			SWdE = -oldtryswenergy + tryswenergy;
			if (P(SWdE) < alpha)
				validmove = false;
		}
		else if (ljb)
		{
			oldtryljenergy = LennardJonesEnergy(oldtry, oldenergy);
			tryljenergy = LennardJonesEnergy(thistry, changedatom);
			if (ljenergy - oldtryljenergy + tryljenergy <= ljenergy)
				goto Accepted;
			LJdE = -oldtryljenergy + tryljenergy;
			if (P(LJdE) < alpha)
				validmove = false;
		}

		Accepted:
		if (validmove == true)
		{
			if (hsb)
			{
				hsenergy += -oldtryhsenergy + tryhsenergy;
				// Grapher doesn't show 0 values so this will not show up unless a small fudge amount is added to hsenergy in this function
				//Grapher.Log(hsenergy + 0.000000001f, "HS Energies");
				Grapher.Log(hsenergy, "HS Energies");
				HSdEL.Add(hsenergy);
				Grapher.Log(HSdEL.Average(), "Avg HS Energy");
			}
			else if (swb)
			{
				swenergy += -oldtryswenergy + tryswenergy;
				Grapher.Log(swenergy, "SW Energies");
				SWdEL.Add(swenergy);
				Grapher.Log(SWdEL.Average(), "Avg SW Energy");
			}
			else if (ljb)
			{
				ljenergy += -oldtryljenergy + tryljenergy;
				Grapher.Log(ljenergy, "LJ Energies");
				LJdEL.Add(ljenergy);
				Grapher.Log(LJdEL.Average(), "Avg LJ Energy");
			}

			accepted.Add(thistry);
			e.transform.localPosition = changedatom;

			// do distance calculation every 20 accepted steps because it slows the program and we don't really care about each individual step, just about convergence
			if (dist && accepted.Count % 20 == 0)
				Debug.Log(String.Format("Average distance: {0} Angstroms.", AvgDis()));
		}

		totsteps++;

		float value = 100f * accepted.Count / (float)totsteps;

		Grapher.Log(value, "Percent accepted");
        */
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
