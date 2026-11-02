
const aminoAcids = {
  ALA: {
    name: "Alanine",
    letter: "A",
    formula: "C3H7NO2",
    mass: 89.09,
    category: "Nonpolar",
    sideChain: "–CH₃",
    description: "A small amino acid that helps proteins pack together. Its simple methyl side chain is useful in many protein structures.",
    waterAffinity: "Hydrophobic",
    inProteins: "Often contributes to stable protein cores and alpha helices."
  },



  
  ARG: {
    name: "Arginine",
    letter: "R",
    formula: "C6H14N4O2",
    mass: 174.20,
    category: "Positively charged",
    sideChain: "Guanidinium-containing",
    description: "A basic amino acid that usually carries a positive charge at physiological pH. It can interact with negatively charged molecules.",
    waterAffinity: "Hydrophilic",
    inProteins: "Forms salt bridges, hydrogen bonds and interactions with DNA."
  },

  ASN: {
    name: "Asparagine",
    letter: "N",
    formula: "C4H8N2O3",
    mass: 132.12,
    category: "Polar, uncharged",
    sideChain: "Amide-containing",
    description: "A polar amino acid whose amide group can form hydrogen bonds with nearby molecules.",
    waterAffinity: "Hydrophilic",
    inProteins: "Helps stabilise structures through hydrogen bonding."
  },

  ASP: {
    name: "Aspartic acid",
    letter: "D",
    formula: "C4H7NO4",
    mass: 133.10,
    category: "Negatively charged",
    sideChain: "Carboxyl-containing",
    description: "An acidic amino acid that is usually negatively charged at physiological pH.",
    waterAffinity: "Hydrophilic",
    inProteins: "Often forms salt bridges and participates in enzyme reactions."
  },

  CYS: {
    name: "Cysteine",
    letter: "C",
    formula: "C3H7NO2S",
    mass: 121.16,
    category: "Polar, uncharged",
    sideChain: "–CH₂–SH",
    description: "Contains sulfur in a thiol group. Two cysteine residues can form a disulfide bond under suitable conditions.",
    waterAffinity: "Context-dependent",
    inProteins: "Disulfide bonds can help stabilise a protein's folded structure."
  },

  GLN: {
    name: "Glutamine",
    letter: "Q",
    formula: "C5H10N2O3",
    mass: 146.15,
    category: "Polar, uncharged",
    sideChain: "Amide-containing",
    description: "A polar amino acid with a relatively long side chain that can form hydrogen bonds.",
    waterAffinity: "Hydrophilic",
    inProteins: "Often participates in molecular recognition and hydrogen bonding."
  },

  GLU: {
    name: "Glutamic acid",
    letter: "E",
    formula: "C5H9NO4",
    mass: 147.13,
    category: "Negatively charged",
    sideChain: "Carboxyl-containing",
    description: "An acidic amino acid that usually carries a negative charge at physiological pH.",
    waterAffinity: "Hydrophilic",
    inProteins: "Can form salt bridges and participate in enzyme catalysis."
  },

  GLY: {
    name: "Glycine",
    letter: "G",
    formula: "C2H5NO2",
    mass: 75.07,
    category: "Special: flexible",
    sideChain: "–H",
    description: "The smallest amino acid. Because its side chain is just hydrogen, its backbone has greater flexibility than most other amino acids.",
    waterAffinity: "Context-dependent",
    inProteins: "Commonly found in tight turns and regions requiring flexibility."
  },

  HIS: {
    name: "Histidine",
    letter: "H",
    formula: "C6H9N3O2",
    mass: 155.16,
    category: "Polar, ionisable",
    sideChain: "Imidazole-containing",
    description: "Its imidazole group can gain or lose a proton near physiological pH, making histidine especially useful in chemical reactions.",
    waterAffinity: "Hydrophilic",
    inProteins: "Frequently involved in enzyme catalysis and metal binding."
  },

  ILE: {
    name: "Isoleucine",
    letter: "I",
    formula: "C6H13NO2",
    mass: 131.17,
    category: "Nonpolar",
    sideChain: "Branched hydrocarbon",
    description: "A hydrophobic amino acid with a branched side chain that tends to avoid water.",
    waterAffinity: "Hydrophobic",
    inProteins: "Often helps form tightly packed protein cores."
  },

  LEU: {
    name: "Leucine",
    letter: "L",
    formula: "C6H13NO2",
    mass: 131.17,
    category: "Nonpolar",
    sideChain: "Branched hydrocarbon",
    description: "A hydrophobic amino acid that is common in the interiors of folded proteins.",
    waterAffinity: "Hydrophobic",
    inProteins: "Contributes to protein packing and is frequently found in alpha helices."
  },

  LYS: {
    name: "Lysine",
    letter: "K",
    formula: "C6H14N2O2",
    mass: 146.19,
    category: "Positively charged",
    sideChain: "Terminal amino group",
    description: "A basic amino acid whose side chain usually carries a positive charge at physiological pH.",
    waterAffinity: "Hydrophilic",
    inProteins: "Often interacts with negatively charged groups and can undergo chemical modifications."
  },

  MET: {
    name: "Methionine",
    letter: "M",
    formula: "C5H11NO2S",
    mass: 149.21,
    category: "Nonpolar",
    sideChain: "Sulfur-containing thioether",
    description: "Contains sulfur and commonly serves as the first amino acid added when a protein is made.",
    waterAffinity: "Hydrophobic",
    inProteins: "Contributes to hydrophobic packing. The initial methionine may be removed after translation."
  },

  PHE: {
    name: "Phenylalanine",
    letter: "F",
    formula: "C9H11NO2",
    mass: 165.19,
    category: "Aromatic, nonpolar",
    sideChain: "Phenyl ring",
    description: "Contains a large aromatic ring that tends to avoid water and interact with other nonpolar groups.",
    waterAffinity: "Hydrophobic",
    inProteins: "Often contributes to protein core packing and aromatic interactions."
  },

  PRO: {
    name: "Proline",
    letter: "P",
    formula: "C5H9NO2",
    mass: 115.13,
    category: "Special: conformationally restricted",
    sideChain: "Ring bonded to backbone nitrogen",
    description: "Its side chain forms a ring with the backbone nitrogen, restricting how the protein chain can bend.",
    waterAffinity: "Hydrophobic",
    inProteins: "Often influences turns and can disrupt or cap alpha helices."
  },

  SER: {
    name: "Serine",
    letter: "S",
    formula: "C3H7NO3",
    mass: 105.09,
    category: "Polar, uncharged",
    sideChain: "–CH₂–OH",
    description: "Contains a hydroxyl group that can form hydrogen bonds and undergo phosphorylation.",
    waterAffinity: "Hydrophilic",
    inProteins: "Important in enzyme reactions, hydrogen bonding and cell signalling."
  },

  THR: {
    name: "Threonine",
    letter: "T",
    formula: "C4H9NO3",
    mass: 119.12,
    category: "Polar, uncharged",
    sideChain: "Hydroxyl-containing",
    description: "A polar amino acid with a hydroxyl group that can participate in hydrogen bonding and phosphorylation.",
    waterAffinity: "Hydrophilic",
    inProteins: "Helps with molecular interactions and protein regulation."
  },

  TRP: {
    name: "Tryptophan",
    letter: "W",
    formula: "C11H12N2O2",
    mass: 204.23,
    category: "Aromatic",
    sideChain: "Indole ring",
    description: "A large aromatic amino acid whose indole ring can participate in hydrophobic interactions and some hydrogen bonds.",
    waterAffinity: "Mostly hydrophobic",
    inProteins: "Often helps stabilise protein structures and molecular interactions."
  },

  TYR: {
    name: "Tyrosine",
    letter: "Y",
    formula: "C9H11NO3",
    mass: 181.19,
    category: "Aromatic, polar",
    sideChain: "Phenolic hydroxyl",
    description: "Combines an aromatic ring with a hydroxyl group, allowing both hydrophobic and polar interactions.",
    waterAffinity: "Amphipathic",
    inProteins: "Can form hydrogen bonds and be phosphorylated during cell signalling."
  },

  VAL: {
    name: "Valine",
    letter: "V",
    formula: "C5H11NO2",
    mass: 117.15,
    category: "Nonpolar",
    sideChain: "Branched hydrocarbon",
    description: "A hydrophobic amino acid with a branched side chain that often contributes to protein stability.",
    waterAffinity: "Hydrophobic",
    inProteins: "Commonly contributes to the tightly packed interiors of proteins."
  }
};

export const aminoAcidList = Object.values(aminoAcids).sort(
  (a, b) => a.name.localeCompare(b.name)
);

export const aminoAcidGroups = [
  ...new Set(aminoAcidList.map((aminoAcid) => aminoAcid.category))
];

export default aminoAcids;
