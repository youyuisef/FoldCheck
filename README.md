# FoldCheck

**Explore protein structures. Investigate predictions. Understand uncertainty.**

FoldCheck is an interactive protein-structure exploration project designed to make complex molecular structures easier to visualise and understand.

Its broader goal is to help users investigate AI-predicted protein structures critically, rather than assuming that every part of a prediction is equally reliable.

**Live Demo:** Coming soon — a permanent public link will be added after deployment.

---

## 1. Why I Created FoldCheck

I've always been interested in biology, chemistry and how artificial intelligence can help scientists understand complex problems.

When I learned about AlphaFold and its ability to predict the three-dimensional structures of proteins, I was fascinated by the possibilities.

However, I realised that seeing a detailed 3D protein model doesn't necessarily tell someone how reliable the prediction is.

That led me to an important question:

**How can we investigate a protein prediction instead of simply accepting what it shows?**

I created FoldCheck to explore that question.

My vision is to develop a tool that combines interactive molecular visualisation, scientific information and confidence analysis in a single accessible environment.

I also wanted to make protein science more approachable for students and people who may not have experience using specialised molecular-visualisation software.

---

## 2. What Is FoldCheck?

FoldCheck is a React-based web application focused on exploring protein structures.

The project brings together an interactive molecular viewer and a user interface designed to help people examine proteins, select individual residues and understand the scientific context behind predicted structures.

Rather than developing a new protein-folding AI model, FoldCheck focuses on **visualising and investigating existing protein structures and predictions**.

Its long-term purpose is to help users move from simply viewing a protein to asking questions about its structural reliability.

---

## 3. Current Features

### Interactive 3D Protein Viewer

FoldCheck integrates Mol*, a molecular visualisation library used for exploring three-dimensional molecular structures.

The viewer allows users to inspect proteins in 3D, rotate structures and zoom into specific regions.

During development, I used crambin (PDB ID: 1CRN) as an example experimental structure and explored an AlphaFold prediction of human p53 (UniProt: P04637).

### Residue Selection and Inspection

FoldCheck connects interactions inside the molecular viewer to the React interface.

Selecting a residue can display information about the selection in a separate details panel, providing a foundation for more advanced residue-level analysis.

### Scientific Introduction

FoldCheck includes an animated introduction designed to explain how proteins relate to the wider biological world.

It introduces concepts such as cells, amino acids, protein folding, AlphaFold predictions and scientific uncertainty.

The aim is to make the project understandable even to users who are relatively new to structural biology.

### Protein Exploration Interface

The application includes a dedicated interface for exploring proteins and navigating between the project's scientific content and interactive visualisation features.

The design aims to balance scientific information with an approachable user experience.

---

## 4. Planned Features and Future Development

FoldCheck remains under development. The following capabilities represent my goals rather than features I am claiming are fully implemented.

**Residue-Level Confidence Analysis**

I want to integrate AlphaFold confidence information, particularly pLDDT, so users can investigate how confidence varies across a predicted structure.

**Predicted Aligned Error (PAE)**

I plan to explore ways of displaying PAE information to help users understand uncertainty in the relative arrangement of protein regions.

**Structural Investigation**

I want to develop tools that help users identify regions worth investigating further, including potentially uncertain or unusual structural features.

**Comparison With Experimental Structures**

A future goal is to allow users to compare predicted protein structures with experimentally determined structures when suitable reference data is available.

**Accessible Scientific Explanations**

I want to provide explanations of confidence scores, structural features and scientific terminology so that FoldCheck can be useful for learning as well as exploration.

**Integrated Investigation Workspace**

My longer-term vision is to connect the molecular viewer, residue selections, confidence information and other scientific evidence into one interactive environment.

---

## 5. The Science Behind FoldCheck

### What Are Proteins?

Proteins are biological molecules made from chains of amino acids.

Their three-dimensional structures are important because molecular shape and interactions contribute to biological function.

Understanding protein structures can help researchers investigate biological processes, molecular interactions and disease-related mechanisms.

### What Is AlphaFold?

AlphaFold is an artificial intelligence system developed to predict protein structures from amino-acid sequences.

Its predictions have contributed significantly to structural biology by making predicted protein structures available for a large number of sequences.

However, a predicted structure should not automatically be treated as experimentally confirmed.

### What Is pLDDT?

The predicted Local Distance Difference Test, or pLDDT, is a measure of confidence in the local structure predicted by AlphaFold.

It is typically represented on a scale from 0 to 100.

Higher scores generally indicate greater confidence in the predicted local structure, while lower scores suggest greater uncertainty.

Importantly, pLDDT is a prediction-confidence measure, not direct experimental proof that a structure is correct.

### What Is PAE?

Predicted Aligned Error, or PAE, estimates uncertainty in the relative positions of residues within a predicted structure.

It can help researchers investigate whether the relative arrangement of different protein regions is confidently predicted.

### Why Does Structural Confidence Matter?

Protein predictions can contain regions with different levels of confidence.

A visually convincing model does not necessarily mean that every region has the same reliability.

FoldCheck's central motivation is to make this distinction easier to investigate and understand.

---

## 6. Technologies Used

FoldCheck has been developed using:

- **React** — for building the interactive user interface.
- **JavaScript** — for application behaviour and interactions.
- **Vite** — for frontend development and build tooling.
- **Mol*** — for three-dimensional molecular visualisation.
- **HTML and CSS** — for layout, styling and visual presentation.

The project also explores scientific concepts and publicly available protein-structure resources.

---

## 7. Development Journey

FoldCheck has progressed through several stages.

**Stage 1 — Initial Concept**

I began by exploring how to make protein-structure predictions more understandable and interactive.

**Stage 2 — Interface Design**

I developed the initial React interface, visual styling and scientific presentation.

**Stage 3 — Scientific Introduction**

I worked on an animated introduction explaining proteins, folding and AI predictions.

**Stage 4 — Molecular Visualisation**

I integrated Mol* and began displaying real protein structures inside the application.

**Stage 5 — Interactive Residue Selection**

I connected residue selection in the 3D viewer with information displayed elsewhere in the interface.

**Next Stage — Confidence Analysis**

My next major objective is to integrate confidence data and develop more advanced tools for investigating predicted structures.

---

## 8. Current Limitations

FoldCheck is an ongoing student-developed project, not a completed scientific validation platform.

Although the project provides interactive protein visualisation and an initial foundation for residue inspection, several advanced analysis features remain in development.

FoldCheck does not currently guarantee the correctness of a protein structure, independently verify AlphaFold predictions or replace experimental structural analysis.

It should not be used as the sole basis for scientific, clinical or medical decisions.

---

## 9. Live Demo

A permanent public demonstration link will be added once FoldCheck has been deployed.

The live demo will allow visitors to explore the application directly through their browsers without installing development software.

---

## 10. Scientific Resources and Acknowledgements

The following resources have been valuable for understanding the scientific concepts and technologies related to FoldCheck:

- [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/)
- [RCSB Protein Data Bank](https://www.rcsb.org/)
- [UniProt](https://www.uniprot.org/)
- [Mol* Molecular Visualisation](https://molstar.org/)

FoldCheck is an independent student project and is not affiliated with or officially endorsed by these organisations.

---

## 11. Copyright and Usage

**© 2026 Yousef Khanfar. All rights reserved.**

FoldCheck is an original student-developed project.

The source code is publicly accessible for project demonstration, inspection and review.

Unless otherwise permitted by an applicable licence or law, no permission is granted to copy, modify, redistribute or republish my original source code without prior authorisation.

Third-party libraries, scientific resources and other external materials remain subject to their respective licences and terms.

---

## 12. About the Creator

FoldCheck was developed as part of my participation in **Hack Club Stardance**.

The project combines my interests in biology, chemistry, artificial intelligence, programming and interactive design.

I hope to continue improving FoldCheck as I learn more about structural biology and software development.

**My ultimate goal is simple: make it easier to question, explore and understand protein predictions — not just look at them.**
