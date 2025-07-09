use c_elegans_chemotaxis::analysis::*;
use c_elegans_chemotaxis::genetic_algorithm::*;
// use c_elegans_chemotaxis::simulation::*;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: cargo run --release -- <command>");
        eprintln!("Commands:");
        eprintln!("  ga-normal");
        eprintln!("  ga-constrain");
        eprintln!("  analysis");
        return;
    }

    match args[1].as_str() {
        "ga-normal" => genetic_algorithm_normal(),
        "ga-constrain" => genetic_algorithm_constrain_aiy_aiz(),
        "analysis" => analytical_control(),
        _ => {
            eprintln!("Unknown command: {}", args[1]);
            eprintln!("Available commands: ga-normal, ga-constrain, analysis");
        }
    }
}
