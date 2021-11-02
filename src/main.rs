extern crate glutin_window;
extern crate graphics;
extern crate opengl_graphics;
extern crate piston;

mod montecarlo;

use piston::Size;
use glutin_window::GlutinWindow as Window;
use opengl_graphics::{GlGraphics, OpenGL};
use piston::event_loop::{EventSettings, Events};
use piston::input::{RenderArgs, RenderEvent, UpdateArgs, UpdateEvent};
use piston::window::WindowSettings;

const N:i32 = 150;
static KB:f32 = 0.0019872036;
static T:i32 = 300; // 100 to 500
static EPSILON:i32 = 5; // 0 to 10
static R:i32 = 5; // 0 to 10
static SIGMA:f64 = 0.5; // 0 to 1
static STEPSIZE:f32 = 0.5; // 0 to 1
static WAITTIME:f32 = 0.00000000001;
static L:f32 = 3.0; // 0 to 5
const WHITE: [f32; 4] = [0.97, 0.97, 0.97, 1.0];
const RED: [f32; 4] = [0.95, 0.05, 0.05, 1.0];
const WINDOWSIZE: Size = Size{width:800.0,height:800.0};

pub struct App {
    gl: GlGraphics, // OpenGL drawing backend.
    x: f64,  // Rotation for the square.
    y: f64,
}

impl App {
    fn render(&mut self, args: &RenderArgs) {
        use graphics::*;

        self.gl.draw(args.viewport(), |c, g| {
            // Clear the screen.
            // use this for a reset button later
            clear(WHITE, g);
});

        for i in 0..N {

            let (x, y) = (self.x, self.y);//(args.window_size[0] / 2.0, args.window_size[1] / 2.0);
            let rect = [(i as f64)*x, (i as f64)*y, (R as f64),(R as f64)]; // x, y, width, height

            self.gl.draw(args.viewport(), |c, g| {
                let transform = c
                    .transform
                    .trans(x, y);

                // Draw a box rotating around the middle of the screen.
                ellipse(RED, rect, transform, g);
                Ellipse::new_border(RED, 2.0).draw(rect, &c.draw_state, transform, g);

            });
        }
    }

    fn update(&mut self, args: &UpdateArgs) {
        self.x += 2.0 * args.dt;
        self.y += 2.0 * args.dt;
    }
}

fn main() {
    let settings = montecarlo::settings::Settings::new();
    let default_state = montecarlo::State::new(N,L);


    // Change this to OpenGL::V2_1 if not working.
    let opengl = OpenGL::V3_2;

    // Create an Glutin window.
    let mut window: Window = WindowSettings::new("Monte Carlo 2D Gas", WINDOWSIZE)
        .graphics_api(opengl)
        .exit_on_esc(true)
        .build()
        .unwrap();

    // Create a new game and run it.
    let mut app = App {
        gl: GlGraphics::new(opengl),
        x: 0.0,
        y: 0.0,
    };

    let mut events = Events::new(EventSettings::new());
    while let Some(e) = events.next(&mut window) {
        if let Some(args) = e.render_args() {
            app.render(&args);
        }

        if let Some(args) = e.update_args() {
            app.update(&args);
        }
    }
}
