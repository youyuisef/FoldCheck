
const atomChemistry = {
  ARG: {
    CD: {
      role: "This carbon is part of arginine's flexible side chain. It connects the hydrocarbon portion to the nitrogen-containing guanidinium group.",
      bonds: ["CG", "NE"],
      significance: "Arginine's positively charged guanidinium group can form electrostatic interactions and hydrogen bonds, including with DNA. CD helps position that group but is not itself positively charged."
    },
    CZ: {
      role: "This carbon is the central carbon of arginine's guanidinium group, bonded to three nitrogen atoms.",
      bonds: ["NE", "NH1", "NH2"],
      significance: "The guanidinium group normally carries a positive charge at physiological pH and can participate in interactions with negatively charged groups."
    }
  },
  GLU: {
    CD: {
      role: "This carbon is the central carbon of glutamate's terminal carboxylate group.",
      bonds: ["CG", "OE1", "OE2"],
      significance: "The carboxylate group is normally negatively charged at physiological pH and can participate in electrostatic interactions."
    }
  },
  GLN: {
    CD: {
      role: "This carbon is the carbonyl carbon of glutamine's side-chain amide group.",
      bonds: ["CG", "OE1", "NE2"],
      significance: "The amide group is polar but normally uncharged. Its oxygen and nitrogen can participate in hydrogen bonding."
    }
  },
  LYS: {
    CD: {
      role: "This carbon is part of lysine's flexible hydrocarbon side chain, connecting CG and CE.",
      bonds: ["CG", "CE"],
      significance: "The side chain terminates in an amino group that is usually positively charged at physiological pH."
    },
    NZ: {
      role: "This is the terminal nitrogen of lysine's side chain.",
      bonds: ["CE"],
      significance: "The terminal amino group is usually protonated and positively charged at physiological pH, making it important in electrostatic interactions."
    }
  },
  PRO: {
    CD: {
      role: "This carbon is part of proline's five-membered ring, which connects the side chain back to the backbone nitrogen.",
      bonds: ["CG", "N"],
      significance: "Proline's ring restricts backbone flexibility and can influence turns and secondary structure."
    }
  },
  TYR: {
    OH: {
      role: "This oxygen belongs to the hydroxyl group attached to tyrosine's aromatic ring.",
      bonds: ["CZ"],
      significance: "The hydroxyl group can participate in hydrogen bonding and can be a site of phosphorylation."
    }
  },
  SER: {
    OG: {
      role: "This oxygen belongs to serine's side-chain hydroxyl group.",
      bonds: ["CB"],
      significance: "The hydroxyl group can participate in hydrogen bonding and may be phosphorylated in regulatory processes."
    }
  },
  THR: {
    OG1: {
      role: "This oxygen belongs to threonine's side-chain hydroxyl group.",
      bonds: ["CB"],
      significance: "The hydroxyl group can participate in hydrogen bonding and may be phosphorylated."
    }
  },
  CYS: {
    SG: {
      role: "This sulfur belongs to cysteine's thiol group.",
      bonds: ["CB"],
      significance: "Cysteine thiols can participate in redox chemistry, coordinate metals, or form disulfide bonds in suitable environments."
    }
  },
  MET: {
    SD: {
      role: "This sulfur is part of methionine's thioether group.",
      bonds: ["CG", "CE"],
      significance: "Methionine's sulfur is normally uncharged and can influence side-chain packing and chemical reactivity."
    }
  },
  ASP: {
    CG: {
      role: "This carbon is the central carbon of aspartate's carboxylate group.",
      bonds: ["CB", "OD1", "OD2"],
      significance: "The carboxylate group is normally negatively charged at physiological pH."
    }
  },
  ASN: {
    CG: {
      role: "This carbon is the carbonyl carbon of asparagine's amide group.",
      bonds: ["CB", "OD1", "ND2"],
      significance: "The polar amide group can participate in hydrogen bonding."
    }
  },
  PHE: {
    CZ: {
      role: "This carbon belongs to phenylalanine's aromatic ring.",
      bonds: ["CE1", "CE2"],
      significance: "The aromatic ring contributes to hydrophobic packing and can participate in aromatic interactions."
    }
  },
  TRP: {
    NE1: {
      role: "This nitrogen is part of tryptophan's indole ring.",
      bonds: ["CD1", "CE2"],
      significance: "The indole nitrogen can donate a hydrogen bond, while the aromatic system contributes to molecular packing."
    }
  },
  HIS: {
    ND1: {
      role: "This nitrogen belongs to histidine's imidazole ring.",
      bonds: ["CG", "CE1"],
      significance: "Histidine can participate in proton transfer, hydrogen bonding and metal coordination, depending on its environment and protonation state."
    },
    NE2: {
      role: "This nitrogen belongs to histidine's imidazole ring.",
      bonds: ["CD2", "CE1"],
      significance: "Histidine's imidazole group can change protonation state near physiological pH, making it chemically versatile."
    }
  }
};

export default atomChemistry;
