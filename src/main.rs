use std::io;

#[derive(Debug, Clone, Copy)]
enum DNABases {
    Adenine,
    Cytosine,
    Guanine,
    Thymine,
}

#[derive(Debug, Clone, Copy)]
enum RNABases {
    Adenine,
    Cytosine,
    Guanine,
    Uracil,
}

#[derive(Debug, Clone, Copy)]
struct Codon {
    base1: DNABases,
    base2: DNABases,
    base3: DNABases,
}

impl Codon {
    fn transcribeToRNA(&self) -> RNACodon {
        fn convertDNAtoRNA(base: DNABases) -> RNABases {
            match base {
                DNABases::Adenine => RNABases::Adenine,
                DNABases::Cytosine => RNABases::Cytosine,
                DNABases::Guanine => RNABases::Guanine,
                DNABases::Thymine => RNABases::Uracil
            }
        }
        RNACodon {
            base1: convertDNAtoRNA(self.base1),
            base2: convertDNAtoRNA(self.base2),
            base3: convertDNAtoRNA(self.base3),
        }
    }

    fn new(bases: &str) -> io::Result<Self> {
        Ok(Codon {
            base1: match &bases[0..1] {
                "A" => DNABases::Adenine,
                "C" => DNABases::Cytosine,
                "G" => DNABases::Guanine,
                "T" => DNABases::Thymine,
                _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid base inputted")),
            },
            base2: match &bases[1..2] { 
                "A" => DNABases::Adenine,
                "C" => DNABases::Cytosine,
                "G" => DNABases::Guanine,
                "T" => DNABases::Thymine,
                _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid base inputted")),
            },
            base3: match &bases[2..3] { 
                "A" => DNABases::Adenine,
                "C" => DNABases::Cytosine,
                "G" => DNABases::Guanine,
                "T" => DNABases::Thymine,
                _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid base inputted")),
            }
        })
    }
}

#[derive(Debug)]
struct TemplateStrand<'a> {
    // number of codons
    length: usize,
    codons: &'a mut Vec<Codon>,
}

#[derive(Debug, Clone, Copy)]
struct RNACodon {
    base1: RNABases,
    base2: RNABases,
    base3: RNABases,
}

#[derive(Debug)]
struct RNA <'a> {
    // number of codons
    length: usize,
    codons: &'a mut Vec<RNACodon>,
}

impl Clone for RNA<'_> {
    fn clone(&'_ self) -> RNA<'_> {
        RNA {
            length: self.length,
            codons: &mut self.codons.clone(),
        }
    }
}

impl TemplateStrand<'_> { 
   fn transcibeToRNA(&self) -> RNA {
        let mut rna: RNA = RNA {
            length: self.length,
            codons: &mut Vec::new()
        };
        for i in 0..self.length {
            (*rna.codons).push(self.codons[i].transcribeToRNA());
        }
        rna
    }
	
    fn from_coding_strand(string: &str) -> Result<Self, io::Error> {
		let mut template = String::new();
		for i in 0..string.len() {
			template += match &string[i..i+1] {
				"A" => "T",
				"C" => "G",
				"G" => "C",
				"T" => "A",
                _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid base inputted")),
			};	
		}
		TemplateStrand::new(&template)
	}

    fn new(string: &str) -> Result<Self, io::Error> {
        if string.len() / 3 != 0 {
            return Err(
                io::Error::new(
                    io::ErrorKind::InvalidInput, 
                    "input string not in 3 base groups (codons)"
                )
            );
        };
        let mut template: TemplateStrand = TemplateStrand {
            length: string.len() / 3,
            codons: &mut Vec::new()
        };
        for i in 0..(string.len() / 3) {
            let letter = i * 3;
            (*template.codons).push(Codon::new(&string[letter..(letter+3)])?);
        }
		Ok(template)
    }
}

#[derive(Debug, Clone, Copy)]
enum AminoAcids {
    Phe,
    Leu,
    Tyr,
    Cys,
    Pro,
    His,
    Gln,
    Arg,
    Ile,
    Met,
    Thr,
    Asn,
    Lys,
    Ser,
    Val,
    Ala,
    Asp,
    Glu,
    Gly,
    Stop,
}
#[derive(Debug)]
struct Protein<'a> {
    length: usize,
    polypeptide: &'a mut Vec<AminoAcids>,
}

impl Clone for Protein<'_> {
    fn clone(&self) -> Protein<'_> {
        Protein {
            length: self.length,
            polypeptide: &mut self.polypeptide.clone(),
        }
    }
}

impl Protein<'_> {
    fn convertRNACodonToAminoAcid(codon: RNACodon) -> AminoAcids {
        match codon.base1 {
            RNABases::Adenine => match codon.base2 {
                RNABases::Adenine => match codon.base3 {
                    RNABases::Adenine => AminoAcids::Lys,
                    RNABases::Cytosine => AminoAcids::Asn,
                    RNABases::Guanine => AminoAcids::Lys,
                    RNABases::Uracil => AminoAcids::Asn,
                },
                RNABases::Cytosine => AminoAcids::Thr,
                RNABases::Guanine => match codon.base3 {
                    RNABases::Adenine => AminoAcids::Arg,
                    RNABases::Cytosine => AminoAcids::Ser,
                    RNABases::Guanine => AminoAcids::Arg,
                    RNABases::Uracil => AminoAcids::Ser,
                },
                RNABases::Uracil => match codon.base3 {
                    RNABases::Guanine => AminoAcids::Met,
                    _ => AminoAcids::Ile,
                }
            },
            RNABases::Cytosine => match codon.base2 {
                RNABases::Adenine => match codon.base3 {
                    RNABases::Adenine => AminoAcids::Gln,
                    RNABases::Cytosine => AminoAcids::His,
                    RNABases::Guanine => AminoAcids::Gln,
                    RNABases::Uracil => AminoAcids::His,
                },
                RNABases::Cytosine => AminoAcids::Pro,
                RNABases::Guanine => AminoAcids::Arg,
                RNABases::Uracil => AminoAcids::Leu,
            },
            RNABases::Guanine => match codon.base2 {
                RNABases::Adenine => match codon.base3 {
                    RNABases::Adenine => AminoAcids::Glu,
                    RNABases::Cytosine => AminoAcids::Asp,
                    RNABases::Guanine => AminoAcids::Glu,
                    RNABases::Uracil => AminoAcids::Asp,
                },
                RNABases::Cytosine => AminoAcids::Ala,
                RNABases::Guanine => AminoAcids::Gly,
                RNABases::Uracil => AminoAcids::Val,
            },
            RNABases::Uracil => match codon.base2 {
                RNABases::Adenine => AminoAcids::Stop,
                RNABases::Cytosine => AminoAcids::Tyr,
                RNABases::Guanine => AminoAcids::Stop,
                RNABases::Uracil => AminoAcids::Tyr,
            },
        }
    }

    fn new(mRNA: RNA) -> Self {
        let mut protein: Protein = Protein {
            length: mRNA.length,
            polypeptide: &mut Vec::new()
        }; 
        for codon in 0..mRNA.length {
            (*protein.polypeptide).push(Protein::convertRNACodonToAminoAcid(mRNA.codons[codon]));
        }
        protein.clone()
    }
}

fn main() -> io::Result<()> {
    print!("Enter a template strand: ");

	let mut input = String::new();

	io::stdin().read_line(&mut input)?;

	println!("{:?}", TemplateStrand::new(&input)?);
    Ok(())
}
