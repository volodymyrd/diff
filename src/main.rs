mod myers;

use termion::{color, style};

use crate::myers::DiffOperation;
use clap::Parser;

fn main() {
    let args = Args::parse();
    println!("{} => {}", args.before, args.after);
    let before = args.before.chars().collect::<Vec<_>>();
    let after = args.after.chars().collect::<Vec<_>>();
    let diff = myers::myers_diff(&before, &after);
    println!("Found {} modifications:", diff.0);
    for d in diff.1 {
        match d {
            DiffOperation::Match(c) => {
                print!("{}{c}", style::Reset);
            }
            DiffOperation::Insertion(c) => {
                print!(
                    "{}{}{}{c}",
                    style::Bold,
                    style::Underline,
                    color::Fg(color::Green)
                );
            }
            DiffOperation::Deletion(c) => {
                print!("{}{}({c})", style::Bold, color::Fg(color::Red));
            }
        }
    }
}

#[derive(Parser, Debug)]
struct Args {
    #[arg(short, long)]
    before: String,

    #[arg(short, long)]
    after: String,
}
