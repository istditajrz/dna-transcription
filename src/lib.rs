use std::io::{self, Read};

#[derive(Debug, Clone, Copy)]
enum DNABases {
    Adenine,
    Cytosine,
    Guanine,
    Thymine,
}

impl DNABases {
	fn as_slice(&self) -> &str {
		match self {
			DNABases::Adenine => "A",
			DNABases::Cytosine => "C",
			DNABases::Guanine => "G",
			DNABases::Thymine => "T",
		}
	}
}

#[derive(Debug, Clone, Copy)]
enum RNABases {
    Adenine,
    Cytosine,
    Guanine,
    Uracil,
}

impl RNABases {
	fn as_slice(&self) -> &str {
		match self {
			RNABases::Adenine => "A",
			RNABases::Cytosine => "C",
			RNABases::Guanine => "G",
			RNABases::Uracil => "U",
		}
	}
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

	fn as_string(&self) -> String {
		format!("{}{}{}", self.base1.as_slice(), self.base2.as_slice(), self.base3.as_slice())
	}
}

#[derive(Debug, Clone)]
struct TemplateStrand {
    // number of codons
    length: usize,
    codons: Vec<Codon>,
}

#[derive(Debug, Clone, Copy)]
struct RNACodon {
    base1: RNABases,
    base2: RNABases,
    base3: RNABases,
}

impl RNACodon {
	fn as_string(&self) -> String {
		format!("{}{}{}", self.base1.as_slice(), self.base2.as_slice(), self.base3.as_slice())
	}
}

#[derive(Debug, Clone)]
struct RNA {
    // number of codons
    length: usize,
    codons: Vec<RNACodon>,
}

impl TemplateStrand { 
   fn transcibeToRNA(&self) -> RNA {	
        RNA {
            length: self.length,
            codons: self.codons
					.as_slice()
					.iter()
					.map(|x| x.transcribeToRNA())
					.collect::<Vec<RNACodon>>()
        }
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
        if string.len() % 3 != 0 {
            return Err(
                io::Error::new(
                    io::ErrorKind::InvalidInput, 
                    "input string not in 3 base groups (codons)"
                )
            );
        };
		let mut vec: Vec<Codon> = Vec::new();
        for i in 0..(string.len() / 3) {
            let letter = i * 3;
            vec.push(Codon::new(&string[letter..(letter+3)])?);
        }
		Ok(TemplateStrand {
            length: string.len() / 3,
            codons: vec
        })
    }
}

#[derive(Debug, Clone, Copy)]
enum AminoAcids {
    Phenylalanine,
    Leucine,
    Tyrosine,
    Cysteine,
    Proline,
    Histidine,
    Glutamine,
    Arginine,
    Isoleucine,
    Methionine,
    Threonine,
    Asparagine,
    Lysine,
    Serine,
    Valine,
    Alanine,
    AsparticAcid,
    GlutamicAcid,
    Glycine,
	Tryptophan,
    Stop,
	NotFound,
}

impl AminoAcids {
	fn as_slice(&self) -> &str {
		match self {
			AminoAcids::Phenylalanine => "TTT",
			AminoAcids::Serine => "TCT",
			AminoAcids::Tyrosine => "TAT",
			AminoAcids::Cysteine => "TGT",
			AminoAcids::Phenylalanine => "TTC",
			AminoAcids::Serine => "TCC",
			AminoAcids::Tyrosine => "TAC",
			AminoAcids::Cysteine => "TGC",
			AminoAcids::Leucine => "TTA",
			AminoAcids::Serine => "TCA",
			AminoAcids::Stop => "TAA",
			AminoAcids::Stop => "TGA",
			AminoAcids::Leucine => "TTG",
			AminoAcids::Serine => "TCG",
			AminoAcids::Stop => "TAG",
			AminoAcids::Tryptophan => "TGG",
			AminoAcids::Leucine => "CTT",
			AminoAcids::Proline => "CCT",
			AminoAcids::Histidine => "CAT",
			AminoAcids::Arginine => "CGT",
			AminoAcids::Leucine => "CTC",
			AminoAcids::Proline => "CCC",
			AminoAcids::Histidine => "CAC",
			AminoAcids::Arginine => "CGC",
			AminoAcids::Leucine => "CTA",
			AminoAcids::Proline => "CCA",
			AminoAcids::Glutamine => "CAA",
			AminoAcids::Arginine => "CGA",
			AminoAcids::Leucine => "CTG",
			AminoAcids::Proline => "CCG",
			AminoAcids::Glutamine => "CAG",
			AminoAcids::Arginine => "CGG",
			AminoAcids::Isoleucine => "ATT",
			AminoAcids::Threonine => "ACT",
			AminoAcids::Asparagine => "AAT",
			AminoAcids::Serine => "AGT",
			AminoAcids::Isoleucine => "ATC",
			AminoAcids::Threonine => "ACC",
			AminoAcids::Asparagine => "AAC",
			AminoAcids::Serine => "AGC",
			AminoAcids::Isoleucine => "ATA",
			AminoAcids::Threonine => "ACA",
			AminoAcids::Lysine => "AAA",
			AminoAcids::Arginine => "AGA",
			AminoAcids::Methionine => "ATG",
			AminoAcids::Threonine => "ACG",
			AminoAcids::Lysine => "AAG",
			AminoAcids::Arginine => "AGG",
			AminoAcids::Valine => "GTT",
			AminoAcids::Alanine => "GCT",
			AminoAcids::AsparticAcid => "GAT",
			AminoAcids::Glycine => "GGT",
			AminoAcids::Valine => "GTC",
			AminoAcids::Alanine => "GCC",
			AminoAcids::AsparticAcid => "GAC",
			AminoAcids::Glycine => "GGC",
			AminoAcids::Valine => "GTA",
			AminoAcids::Alanine => "GCA",
			AminoAcids::GlutamicAcid => "GAA",
			AminoAcids::Glycine => "GGA",
			AminoAcids::Valine => "GTG",
			AminoAcids::Alanine => "GCG",
			AminoAcids::GlutamicAcid => "GAG",
			AminoAcids::Glycine => "GGG",
			AminoAcids::NotFound => "XXX"
		}
	}

	fn as_char(&self) -> &str {
		match self {
			AminoAcids::Phenylalanine => "F",
			AminoAcids::Serine => "S",
			AminoAcids::Tyrosine => "Y",
			AminoAcids::Cysteine => "C",
			AminoAcids::Leucine => "L",
			AminoAcids::Stop => "*",
			AminoAcids::Tryptophan => "W",
			AminoAcids::Proline => "P",
			AminoAcids::Histidine => "H",
			AminoAcids::Arginine => "R",
			AminoAcids::Glutamine => "Q",
			AminoAcids::Isoleucine => "I",
			AminoAcids::Threonine => "T",
			AminoAcids::Asparagine => "N",
			AminoAcids::Lysine => "K",
			AminoAcids::Methionine => "M",
			AminoAcids::Valine => "V",
			AminoAcids::Alanine => "A",
			AminoAcids::AsparticAcid => "D",
			AminoAcids::Glycine => "G",
			AminoAcids::GlutamicAcid => "E",
			AminoAcids::NotFound => "X",
		}
	}

	fn as_short(&self) -> &str {
		match self {
			AminoAcids::Phenylalanine => "Phe",
			AminoAcids::Serine => "Ser",
			AminoAcids::Tyrosine => "Tyr",
			AminoAcids::Cysteine => "Cys",
			AminoAcids::Leucine => "Leu",
			AminoAcids::Stop => "Ter",
			AminoAcids::Tryptophan => "Trp",
			AminoAcids::Proline => "Pro",
			AminoAcids::Histidine => "His",
			AminoAcids::Arginine => "Arg",
			AminoAcids::Glutamine => "Gln",
			AminoAcids::Isoleucine => "Ile",
			AminoAcids::Threonine => "Thr",
			AminoAcids::Asparagine => "Asn",
			AminoAcids::Lysine => "Lys",
			AminoAcids::Methionine => "Met",
			AminoAcids::Valine => "Val",
			AminoAcids::Alanine => "Ala",
			AminoAcids::AsparticAcid => "Asp",
			AminoAcids::Glycine => "Gly",
			AminoAcids::GlutamicAcid => "Glu",
			AminoAcids::NotFound => "None",
		}
	}

	fn as_full(&self) -> &str {
		match self {
			AminoAcids::Phenylalanine => "Phenylalanine",
			AminoAcids::Serine => "Serine",
			AminoAcids::Tyrosine => "Tyrosine",
			AminoAcids::Cysteine => "Cysteine",
			AminoAcids::Leucine => "Leucine",
			AminoAcids::Stop => "Stop",
			AminoAcids::Tryptophan => "Tryptophan",
			AminoAcids::Proline => "Proline",
			AminoAcids::Histidine => "Histidine",
			AminoAcids::Arginine => "Arginine",
			AminoAcids::Glutamine => "Glutamine",
			AminoAcids::Isoleucine => "Isoleucine",
			AminoAcids::Threonine => "Threonine",
			AminoAcids::Asparagine => "Asparagine",
			AminoAcids::Lysine => "Lysine",
			AminoAcids::Methionine => "Methionine",
			AminoAcids::Valine => "Valine",
			AminoAcids::Alanine => "Alanine",
			AminoAcids::AsparticAcid => "AsparticAcid",
			AminoAcids::Glycine => "Glycine",
			AminoAcids::GlutamicAcid => "GlutamicAcid",
			AminoAcids::NotFound => "None",
		}
	}
}


#[derive(Debug)]
struct Protein {
    length: usize,
    polypeptide: Vec<AminoAcids>,
}

impl Protein {
    fn convert_RNACodon_to_AminoAcids(codon: &RNACodon) -> AminoAcids {
        match codon.as_string().as_str() {
			"UUU" => AminoAcids::Phenylalanine,
			"UCU" => AminoAcids::Serine,
			"UAU" => AminoAcids::Tyrosine,
			"UGU" => AminoAcids::Cysteine,
			"UUC" => AminoAcids::Phenylalanine,
			"UCC" => AminoAcids::Serine,
			"UAC" => AminoAcids::Tyrosine,
			"UGC" => AminoAcids::Cysteine,
			"UUA" => AminoAcids::Leucine,
			"UCA" => AminoAcids::Serine,
			"UAA" => AminoAcids::Stop,
			"UGA" => AminoAcids::Stop,
			"UUG" => AminoAcids::Leucine,
			"UCG" => AminoAcids::Serine,
			"UAG" => AminoAcids::Stop,
			"UGG" => AminoAcids::Tryptophan,
			"CUU" => AminoAcids::Leucine,
			"CCU" => AminoAcids::Proline,
			"CAU" => AminoAcids::Histidine,
			"CGU" => AminoAcids::Arginine,
			"CUC" => AminoAcids::Leucine,
			"CCC" => AminoAcids::Proline,
			"CAC" => AminoAcids::Histidine,
			"CGC" => AminoAcids::Arginine,
			"CUA" => AminoAcids::Leucine,
			"CCA" => AminoAcids::Proline,
			"CAA" => AminoAcids::Glutamine,
			"CGA" => AminoAcids::Arginine,
			"CUG" => AminoAcids::Leucine,
			"CCG" => AminoAcids::Proline,
			"CAG" => AminoAcids::Glutamine,
			"CGG" => AminoAcids::Arginine,
			"AUU" => AminoAcids::Isoleucine,
			"ACU" => AminoAcids::Threonine,
			"AAU" => AminoAcids::Asparagine,
			"AGU" => AminoAcids::Serine,
			"AUC" => AminoAcids::Isoleucine,
			"ACC" => AminoAcids::Threonine,
			"AAC" => AminoAcids::Asparagine,
			"AGC" => AminoAcids::Serine,
			"AUA" => AminoAcids::Isoleucine,
			"ACA" => AminoAcids::Threonine,
			"AAA" => AminoAcids::Lysine,
			"AGA" => AminoAcids::Arginine,
			"AUG" => AminoAcids::Methionine,
			"ACG" => AminoAcids::Threonine,
			"AAG" => AminoAcids::Lysine,
			"AGG" => AminoAcids::Arginine,
			"GUU" => AminoAcids::Valine,
			"GCU" => AminoAcids::Alanine,
			"GAU" => AminoAcids::AsparticAcid,
			"GGU" => AminoAcids::Glycine,
			"GUC" => AminoAcids::Valine,
			"GCC" => AminoAcids::Alanine,
			"GAC" => AminoAcids::AsparticAcid,
			"GGC" => AminoAcids::Glycine,
			"GUA" => AminoAcids::Valine,
			"GCA" => AminoAcids::Alanine,
			"GAA" => AminoAcids::GlutamicAcid,
			"GGA" => AminoAcids::Glycine,
			"GUG" => AminoAcids::Valine,
			"GCG" => AminoAcids::Alanine,
			"GAG" => AminoAcids::GlutamicAcid,
			"GGG" => AminoAcids::Glycine,
			_ => AminoAcids::NotFound
		}
    }

    fn new(mRNA: RNA) -> Self {
        // let mut protein: Protein = Protein {
        //     length: mRNA.length,
        //     polypeptide: &mut Vec::new()
        // }; 
        // for codon in 0..mRNA.length {
        //     (*protein.polypeptide).push(Protein::convertRNACodonToAminoAcid(mRNA.codons[codon]));
        // }
        // protein.clone()

		Protein {
			length: mRNA.length,
			polypeptide: mRNA.codons
							.iter()
							.map(|x| Protein::convert_RNACodon_to_AminoAcids(x))
							.collect::<Vec<AminoAcids>>()
		}
    }
}
