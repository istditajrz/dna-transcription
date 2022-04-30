#![allow(non_snake_case, unreachable_patterns)]
use std::{io, convert::TryFrom};

macro_rules! pub_struct {
    ($name:ident {$($field:ident: $t:ty,)*}) => {
        #[derive(Debug, Clone)] // ewww
        pub struct $name {
            $(pub $field: $t),*
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum DNABases {
    Adenine,
    Cytosine,
    Guanine,
    Thymine,
}

impl DNABases {
	pub fn as_slice(&self) -> &str {
		match self {
			DNABases::Adenine => "A",
			DNABases::Cytosine => "C",
			DNABases::Guanine => "G",
			DNABases::Thymine => "T",
		}
	}
    fn to_RNA(&self) -> RNABases {
        match self {
            DNABases::Adenine => RNABases::Adenine,
            DNABases::Cytosine => RNABases::Cytosine,
            DNABases::Guanine => RNABases::Guanine,
            DNABases::Thymine => RNABases::Uracil
        }
    }
}

impl std::convert::TryFrom<&str> for DNABases {
    type Error = io::Error;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "A" => Ok(DNABases::Adenine),
            "C" => Ok(DNABases::Cytosine),
            "G" => Ok(DNABases::Guanine),
            "T" => Ok(DNABases::Thymine),
            _ =>   Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid base inputted")),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum RNABases {
    Adenine,
    Cytosine,
    Guanine,
    Uracil,
}

impl RNABases {
	pub fn as_slice(&self) -> &str {
		match self {
			RNABases::Adenine => "A",
			RNABases::Cytosine => "C",
			RNABases::Guanine => "G",
			RNABases::Uracil => "U",
		}
	}
}


pub_struct!(Codon {
    base1: DNABases,
    base2: DNABases,
    base3: DNABases,
});

impl Codon {
    pub fn transcribeToRNA(&self) -> RNACodon {
        RNACodon {
            base1: self.base1.to_RNA(),
            base2: self.base2.to_RNA(),
            base3: self.base3.to_RNA(),
        }
    }

    pub fn new(bases: &str) -> io::Result<Self> {
        Ok(Codon {
            base1: DNABases::try_from(&bases[0..1])?,
            base2: DNABases::try_from(&bases[1..2])?,
            base3: DNABases::try_from(&bases[2..3])?,
        })
    }
}
impl std::fmt::Display for Codon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
		write!(f, "{}{}{}", self.base1.as_slice(), self.base2.as_slice(), self.base3.as_slice())
	}
}


pub_struct!(TemplateStrand {
    // number of codons
    length: usize,
    codons: Vec<Codon>,
});


pub_struct!(RNACodon {
    base1: RNABases,
    base2: RNABases,
    base3: RNABases,
});

impl std::fmt::Display for RNACodon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
		write!(f, "{}{}{}", self.base1.as_slice(), self.base2.as_slice(), self.base3.as_slice())
	}
}


pub_struct!(RNA {
    // number of codons
    length: usize,
    codons: Vec<RNACodon>,
});

impl TemplateStrand { 
    pub fn transcribe_to_RNA(&self) -> RNA {	
        RNA {
            length: self.length,
            codons: self.codons
					.as_slice()
					.iter()
					.map(|x| x.transcribeToRNA())
					.collect::<Vec<RNACodon>>()
        }
    }
	
    pub fn from_coding_strand(string: &str) -> Result<Self, io::Error> {
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

    pub fn new(string: &str) -> Result<Self, io::Error> {
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
pub enum AminoAcids {
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

impl std::convert::Into<&'static str> for AminoAcids {
    fn into(self) -> &'static str {
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
}
impl std::convert::Into<char> for AminoAcids {
    fn into(self) -> char {
		match self {
			AminoAcids::Phenylalanine   => 'F',
			AminoAcids::Serine          => 'S',
			AminoAcids::Tyrosine        => 'Y',
			AminoAcids::Cysteine        => 'C',
			AminoAcids::Leucine         => 'L',
			AminoAcids::Stop            => '*',
			AminoAcids::Tryptophan      => 'W',
			AminoAcids::Proline         => 'P',
			AminoAcids::Histidine       => 'H',
			AminoAcids::Arginine        => 'R',
			AminoAcids::Glutamine       => 'Q',
			AminoAcids::Isoleucine      => 'I',
			AminoAcids::Threonine       => 'T',
			AminoAcids::Asparagine      => 'N',
			AminoAcids::Lysine          => 'K',
			AminoAcids::Methionine      => 'M',
			AminoAcids::Valine          => 'V',
			AminoAcids::Alanine         => 'A',
			AminoAcids::AsparticAcid    => 'D',
			AminoAcids::Glycine         => 'G',
			AminoAcids::GlutamicAcid    => 'E',
			AminoAcids::NotFound        => 'X',
		}
	}
}
impl AminoAcids {
	pub fn as_short(&self) -> &str {
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

	pub fn as_full(&self) -> &str {
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

impl std::convert::From<&RNACodon> for AminoAcids {
    fn from(codon: &RNACodon) -> Self {
        match (codon.base1, codon.base2, codon.base3) {
            (RNABases::Uracil, RNABases::Uracil, RNABases::Uracil)          => AminoAcids::Phenylalanine,
            (RNABases::Uracil, RNABases::Cytosine, RNABases::Uracil)        => AminoAcids::Serine,
            (RNABases::Uracil, RNABases::Adenine, RNABases::Uracil)         => AminoAcids::Tyrosine,
            (RNABases::Uracil, RNABases::Guanine, RNABases::Uracil)         => AminoAcids::Cysteine,
            (RNABases::Uracil, RNABases::Uracil, RNABases::Cytosine)        => AminoAcids::Phenylalanine,
            (RNABases::Uracil, RNABases::Cytosine, RNABases::Cytosine)      => AminoAcids::Serine,
            (RNABases::Uracil, RNABases::Adenine, RNABases::Cytosine)       => AminoAcids::Tyrosine,
            (RNABases::Uracil, RNABases::Guanine, RNABases::Cytosine)       => AminoAcids::Cysteine,
            (RNABases::Uracil, RNABases::Uracil, RNABases::Adenine)         => AminoAcids::Leucine,
            (RNABases::Uracil, RNABases::Cytosine, RNABases::Adenine)       => AminoAcids::Serine,
            (RNABases::Uracil, RNABases::Adenine, RNABases::Adenine)        => AminoAcids::Stop,
            (RNABases::Uracil, RNABases::Guanine, RNABases::Adenine)        => AminoAcids::Stop,
            (RNABases::Uracil, RNABases::Uracil, RNABases::Guanine)         => AminoAcids::Leucine,
            (RNABases::Uracil, RNABases::Cytosine, RNABases::Guanine)       => AminoAcids::Serine,
            (RNABases::Uracil, RNABases::Adenine, RNABases::Guanine)        => AminoAcids::Stop,
            (RNABases::Uracil, RNABases::Guanine, RNABases::Guanine)        => AminoAcids::Tryptophan,
            (RNABases::Cytosine, RNABases::Uracil, RNABases::Uracil)        => AminoAcids::Leucine,
            (RNABases::Cytosine, RNABases::Cytosine, RNABases::Uracil)      => AminoAcids::Proline,
            (RNABases::Cytosine, RNABases::Adenine, RNABases::Uracil)       => AminoAcids::Histidine,
            (RNABases::Cytosine, RNABases::Guanine, RNABases::Uracil)       => AminoAcids::Arginine,
            (RNABases::Cytosine, RNABases::Uracil, RNABases::Cytosine)      => AminoAcids::Leucine,
            (RNABases::Cytosine, RNABases::Cytosine, RNABases::Cytosine)    => AminoAcids::Proline,
            (RNABases::Cytosine, RNABases::Adenine, RNABases::Cytosine)     => AminoAcids::Histidine,
            (RNABases::Cytosine, RNABases::Guanine, RNABases::Cytosine)     => AminoAcids::Arginine,
            (RNABases::Cytosine, RNABases::Uracil, RNABases::Adenine)       => AminoAcids::Leucine,
            (RNABases::Cytosine, RNABases::Cytosine, RNABases::Adenine)     => AminoAcids::Proline,
            (RNABases::Cytosine, RNABases::Adenine, RNABases::Adenine)      => AminoAcids::Glutamine,
            (RNABases::Cytosine, RNABases::Guanine, RNABases::Adenine)      => AminoAcids::Arginine,
            (RNABases::Cytosine, RNABases::Uracil, RNABases::Guanine)       => AminoAcids::Leucine,
            (RNABases::Cytosine, RNABases::Cytosine, RNABases::Guanine)     => AminoAcids::Proline,
            (RNABases::Cytosine, RNABases::Adenine, RNABases::Guanine)      => AminoAcids::Glutamine,
            (RNABases::Cytosine, RNABases::Guanine, RNABases::Guanine)      => AminoAcids::Arginine,
            (RNABases::Adenine, RNABases::Uracil, RNABases::Uracil)         => AminoAcids::Isoleucine,
            (RNABases::Adenine, RNABases::Cytosine, RNABases::Uracil)       => AminoAcids::Threonine,
            (RNABases::Adenine, RNABases::Adenine, RNABases::Uracil)        => AminoAcids::Asparagine,
            (RNABases::Adenine, RNABases::Guanine, RNABases::Uracil)        => AminoAcids::Serine,
            (RNABases::Adenine, RNABases::Uracil, RNABases::Cytosine)       => AminoAcids::Isoleucine,
            (RNABases::Adenine, RNABases::Cytosine, RNABases::Cytosine)     => AminoAcids::Threonine,
            (RNABases::Adenine, RNABases::Adenine, RNABases::Cytosine)      => AminoAcids::Asparagine,
            (RNABases::Adenine, RNABases::Guanine, RNABases::Cytosine)      => AminoAcids::Serine,
            (RNABases::Adenine, RNABases::Uracil, RNABases::Adenine)        => AminoAcids::Isoleucine,
            (RNABases::Adenine, RNABases::Cytosine, RNABases::Adenine)      => AminoAcids::Threonine,
            (RNABases::Adenine, RNABases::Adenine, RNABases::Adenine)       => AminoAcids::Lysine,
            (RNABases::Adenine, RNABases::Guanine, RNABases::Adenine)       => AminoAcids::Arginine,
            (RNABases::Adenine, RNABases::Uracil, RNABases::Guanine)        => AminoAcids::Methionine,
            (RNABases::Adenine, RNABases::Cytosine, RNABases::Guanine)      => AminoAcids::Threonine,
            (RNABases::Adenine, RNABases::Adenine, RNABases::Guanine)       => AminoAcids::Lysine,
            (RNABases::Adenine, RNABases::Guanine, RNABases::Guanine)       => AminoAcids::Arginine,
            (RNABases::Guanine, RNABases::Uracil, RNABases::Uracil)         => AminoAcids::Valine,
            (RNABases::Guanine, RNABases::Cytosine, RNABases::Uracil)       => AminoAcids::Alanine,
            (RNABases::Guanine, RNABases::Adenine, RNABases::Uracil)        => AminoAcids::AsparticAcid,
            (RNABases::Guanine, RNABases::Guanine, RNABases::Uracil)        => AminoAcids::Glycine,
            (RNABases::Guanine, RNABases::Uracil, RNABases::Cytosine)       => AminoAcids::Valine,
            (RNABases::Guanine, RNABases::Cytosine, RNABases::Cytosine)     => AminoAcids::Alanine,
            (RNABases::Guanine, RNABases::Adenine, RNABases::Cytosine)      => AminoAcids::AsparticAcid,
            (RNABases::Guanine, RNABases::Guanine, RNABases::Cytosine)      => AminoAcids::Glycine,
            (RNABases::Guanine, RNABases::Uracil, RNABases::Adenine)        => AminoAcids::Valine,
            (RNABases::Guanine, RNABases::Cytosine, RNABases::Adenine)      => AminoAcids::Alanine,
            (RNABases::Guanine, RNABases::Adenine, RNABases::Adenine)       => AminoAcids::GlutamicAcid,
            (RNABases::Guanine, RNABases::Guanine, RNABases::Adenine)       => AminoAcids::Glycine,
            (RNABases::Guanine, RNABases::Uracil, RNABases::Guanine)        => AminoAcids::Valine,
            (RNABases::Guanine, RNABases::Cytosine, RNABases::Guanine)      => AminoAcids::Alanine,
            (RNABases::Guanine, RNABases::Adenine, RNABases::Guanine)       => AminoAcids::GlutamicAcid,
            (RNABases::Guanine, RNABases::Guanine, RNABases::Guanine)       => AminoAcids::Glycine,
            _                                                               => AminoAcids::NotFound
        }
    }
}
impl std::fmt::Display for AminoAcids {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_full())
    }
}

pub_struct!(Protein {
    length: usize,
    polypeptide: Vec<AminoAcids>,
});

impl Protein {
    pub fn new(mRNA: RNA) -> Self {
		Protein {
			length: mRNA.length,
			polypeptide: mRNA.codons
							.iter()
							.map(|x: &RNACodon| x.into())
							.collect::<Vec<AminoAcids>>()
		}
    }
}
