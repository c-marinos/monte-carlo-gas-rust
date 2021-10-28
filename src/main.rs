extern crate glutin_window;
extern crate graphics;
extern crate opengl_graphics;
extern crate piston;

use glutin_window::GlutinWindow as Window;
use opengl_graphics::{GlGraphics, OpenGL};
use piston::event_loop::{EventSettings, Events};
use piston::input::{RenderArgs, RenderEvent, UpdateArgs, UpdateEvent};
use piston::window::WindowSettings;

//static KB:f32 = 0.0019872036;
//static T:i32 = 300;

pub struct App {
    gl: GlGraphics, // OpenGL drawing backend.
    x: f64,  // Rotation for the square.
}

impl App {
    fn render(&mut self, args: &RenderArgs) {
        use graphics::*;

        const WHITE: [f32; 4] = [0.97, 0.97, 0.97, 1.0];
        const RED: [f32; 4] = [1.0, 0.0, 0.0, 1.0];

        let (x, y) = (self.x,0.0);//(args.window_size[0] / 2.0, args.window_size[1] / 2.0);
        let rect = [x,y,50.0,50.0]; // x, y, width, height

        self.gl.draw(args.viewport(), |c, g| {
            // Clear the screen.
            clear(WHITE, g);

            let transform = c
                .transform
                .trans(x, y);

            // Draw a box rotating around the middle of the screen.
            ellipse(RED, rect, transform, g);
            Ellipse::new_border(RED, 2.0).draw(rect, &c.draw_state, transform, g);
        });
    }

    fn update(&mut self, args: &UpdateArgs) {
        self.x += 2.0 * args.dt;
    }
}

fn main() {
    // Change this to OpenGL::V2_1 if not working.
    let opengl = OpenGL::V3_2;

    // Create an Glutin window.
    let mut window: Window = WindowSettings::new("spinning-square", [200, 200])
        .graphics_api(opengl)
        .exit_on_esc(true)
        .build()
        .unwrap();

    // Create a new game and run it.
    let mut app = App {
        gl: GlGraphics::new(opengl),
        x: 0.0,
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
