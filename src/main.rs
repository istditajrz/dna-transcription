use std::io::{self, Read};
use dna;

fn main() -> io::Result<()> {
    println!("Enter a template strand: ");

	let mut input = String::new();

	io::stdin().read_line(&mut input)?;
	input.make_ascii_uppercase();
	input = input.replace(" ", "");

	println!(
		"{:?}", 
		dna::Protein::new(dna::TemplateStrand::new(input.trim())?.transcribe_to_RNA()).polypeptide
			.iter()
			.map(|x| x.as_full())
			.collect::<Vec<&str>>()
	);
    Ok(())
}
