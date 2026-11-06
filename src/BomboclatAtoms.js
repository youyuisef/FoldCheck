const atomChemistry = {
  ARG: {
    CD: {
      role: "This carbon is part of arginine's flexible side chain. It connects the hydrocarbon section to the nitrogen-containing guanidinium group.",
      bonds: ["CG", "NE"],
      significance: "Arginine's guanidinium group is usually positively charged at physiological pH. It can form electrostatic interactions and hydrogen bonds, including with DNA. CD helps position this group but is not itself positively charged."
    },
    CZ: {
      role: "This is the central carbon of arginine's guanidinium group, bonded to three nitrogen atoms.",
      bonds: ["NE", "NH1", "NH2"],
      significance: "The guanidinium group normally carries a positive charge at physiological pH and can interact with negatively charged groups."
    }
  },

  GLU: {
    CD: {
      role: "This carbon is the central carbon of glutamate's terminal carboxylate group.",
      bonds: ["CG", "OE1", "OE2"],
      significance: "The carboxylate group is usually negatively charged at physiological pH and can form electrostatic interactions, including salt bridges."
    }
  },
  GLN: {
    CD: {
      role: "This carbon is the carbonyl carbon of glutamine's side-chain amide group.",
      bonds: ["CG", "OE1", "NE2"],
      significance: "The amide group is polar but normally uncharged. Its oxygen can accept hydrogen bonds, while its nitrogen can donate them."
    }
  },

  LYS: {
    CD: {
      role: "This carbon is part of lysine's flexible hydrocarbon side chain, connecting CG and CE.",
      bonds: ["CG", "CE"],
      significance: "Lysine's side chain ends in an amino group that is usually positively charged at physiological pH."
    },
    NZ: {
      role: "This is the terminal nitrogen of lysine's side chain.",
      bonds: ["CE"],
      significance: "The terminal amino group is usually protonated and positively charged at physiological pH, allowing it to form electrostatic interactions with negatively charged groups."
    }
  },

  PRO: {
    CD: {
      role: "This carbon is part of proline's five-membered ring, which connects its side chain back to the backbone nitrogen.",
      bonds: ["CG", "N"],
      significance: "Proline's ring restricts backbone movement, affecting how the protein chain bends and folds."
    }
  },
  TYR: {
    OH: {
      role: "This oxygen belongs to the hydroxyl group attached to tyrosine's aromatic ring.",
      bonds: ["CZ"],
      significance: "The hydroxyl group can form hydrogen bonds and may be phosphorylated during cell signalling."
    }
  },

  SER: {
    OG: {
      role: "This oxygen belongs to serine's side-chain hydroxyl group.",
      bonds: ["CB"],
      significance: "The hydroxyl group can form hydrogen bonds and may be phosphorylated to help regulate protein activity."
    }
  },

  THR: {
    OG1: {
      role: "This oxygen belongs to threonine's side-chain hydroxyl group.",
      bonds: ["CB"],
      significance: "The hydroxyl group can form hydrogen bonds and may be phosphorylated during protein regulation."
    }
  },
  CYS: {
    SG: {
      role: "This sulfur belongs to cysteine's thiol group.",
      bonds: ["CB"],
      significance: "Cysteine's thiol group can participate in redox reactions, bind certain metals and form disulfide bonds under suitable conditions."
    }
  },

  MET: {
    SD: {
      role: "This sulfur is part of methionine's thioether group.",
      bonds: ["CG", "CE"],
      significance: "Methionine's sulfur is normally uncharged. It can affect side-chain packing and participate in certain chemical reactions."
    }
  },

  ASP: {
    CG: {
      role: "This carbon is the central carbon of aspartate's carboxylate group.",
      bonds: ["CB", "OD1", "OD2"],
      significance: "The carboxylate group is normally negatively charged at physiological pH and can interact with positively charged groups."
    }
  },
  ASN: {
    CG: {
      role: "This carbon is the carbonyl carbon of asparagine's amide group.",
      bonds: ["CB", "OD1", "ND2"],
      significance: "The polar amide group can form hydrogen bonds, helping with protein stability and molecular recognition."
    }
  },

  PHE: {
    CZ: {
      role: "This carbon belongs to phenylalanine's aromatic ring.",
      bonds: ["CE1", "CE2"],
      significance: "The aromatic ring contributes to hydrophobic packing and can participate in interactions with other aromatic groups."
    }
  },
  TRP: {
    NE1: {
      role: "This nitrogen is part of tryptophan's indole ring.",
      bonds: ["CD1", "CE2"],
      significance: "The indole nitrogen can donate a hydrogen bond, while the aromatic ring system contributes to protein packing and molecular interactions."
    }
  },

  HIS: {
    ND1: {
      role: "This nitrogen belongs to histidine's imidazole ring.",
      bonds: ["CG", "CE1"],
      significance: "Histidine can participate in proton transfer, hydrogen bonding and metal coordination, depending on its environment and protonation state."
    },
    NE2: {
      role: "This nitrogen also belongs to histidine's imidazole ring.",
      bonds: ["CD2", "CE1"],
      significance: "Histidine's imidazole group can gain or lose a proton near physiological pH, making it useful in enzyme reactions."
    }
  }
};

export default atomChemistry;
