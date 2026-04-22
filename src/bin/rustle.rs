// Primary `rustle` binary: StringTie-compatible assembler port (Rust).
//
// Entrypoint is now rooted under `src/rustle/` as the single active codebase.
#[path = "../rustle/main.rs"]
mod rustle_cli;

fn main() -> anyhow::Result<()> {
    rustle_cli::run_cli()
}
