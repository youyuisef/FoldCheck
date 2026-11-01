FOLDCHECK

Proteins are some of the most important molecules in our bodies, responsible for multiple processes that keep us alive. Their shapes play a huge role in how they function, and thanks to AI systems like AlphaFold, we can now predict what these complicated 3D structures might look like. But what if those predictions aren't completely accurate? How would we even know which parts to trust?



Proteins are made from amino-acid chains that fold into three-dimensional structures. Their shapes and interactions are important for biological functions such as chemical reactions, transport, signalling and structural support. Understanding protein structures can help scientists investigate how proteins function and how structural changes might affect them.



Hi! My name is Yousef, and I created FoldCheck for Hack Club Stardance! FoldCheck is a web application that allows you to explore protein structures in 3D, investigate individual residues and atoms, and examine how reliable AI-predicted protein structures actually are.



I've always been interested in biology, chemistry and AI, so I wanted to combine all three into one project. When I first learned about AlphaFold, I thought it was incredible that AI could predict something as complicated as a protein's structure. However, I also realised that just because a prediction looks convincing, it doesn't necessarily mean every part of it is correct.



That's where FoldCheck comes in!



LIVE WEBSITE

https://fold-check-seven.vercel.app

You can open FoldCheck directly in your browser without installing anything!



WHAT CAN FOLDCHECK DO



Interactive 3D protein viewer



One of the biggest parts of FoldCheck is its interactive molecular viewer, which I built using Mol*, a molecular visualisation library.



You can rotate proteins, zoom into different areas and inspect their structures in 3D.

During development, I used crambin (PDB ID: 1CRN) as an experimental protein structure and human p53 (UniProt ID: P04637) as an example AlphaFold prediction.



Residue and atom inspection



Proteins are made from chains of amino acids, and each amino acid within a protein chain is generally called a residue.

FoldCheck allows users to interact with protein structures and inspect residues and atoms through the viewer.

I connected the Mol* viewer with React so information about selections can appear elsewhere in the interface.

Getting the interactions between the viewer and the rest of the application working was definitely one of the more challenging parts of development 😭



pLDDT confidence analysis



FoldCheck includes pLDDT confidence analysis, which helps users investigate how confident AlphaFold is about different regions of a predicted protein structure.

pLDDT stands for predicted Local Distance Difference Test and uses scores from 0 to 100.

Generally, scores above 90 indicate very high confidence, 70–90 indicate confidence, 50–70 indicate low confidence, and below 50 indicate very low confidence.

However, even a very high score doesn't prove that a structure is experimentally correct!



PAE analysis



FoldCheck also includes PAE heatmaps.

PAE stands for Predicted Aligned Error, and it helps estimate uncertainty in the relative positions of different residues or regions.

For example, AlphaFold might be confident about the individual shapes of two regions, but less confident about their positions relative to each other.

PAE helps investigate these differences.



Structural geometry



I also worked on structural geometry investigation, allowing users to examine different aspects of protein structures beyond just their overall appearance.

The goal is to make it easier to investigate the structure itself instead of simply rotating a 3D model.



Functional regions



FoldCheck includes functional region exploration, allowing users to interact with highlighted regions of a protein and investigate them.

This helps connect the actual structure with biological information about different parts of the protein.



Experimental structure comparison



Another feature is experimental structure comparison.

This allows users to investigate predicted structures alongside experimentally determined structures, including highlighting linked residues.

I wanted FoldCheck to provide more than just confidence scores, since comparing predictions with experimental evidence is another important way of investigating protein structures.



Scientific introduction



I also created an animated scientific introduction explaining proteins, amino acids, folding, AlphaFold and uncertainty.

I wanted people who aren't familiar with structural biology to understand what FoldCheck is actually investigating.



DEVELOPMENT



I built FoldCheck using React, JavaScript, Vite, HTML and CSS, with Mol* handling the molecular visualisation.

I also used scientific resources and data from UniProt, the AlphaFold Protein Structure Database and the Protein Data Bank.

I started by developing the interface and scientific introduction, before moving on to integrating the actual 3D protein viewer.

One of the biggest challenges was connecting Mol* interactions with React, especially getting residue selection, highlighting and camera behaviour working together without breaking other features.

After getting the viewer working, I continued developing confidence analysis, PAE visualisation, structural investigation tools and experimental comparisons.

There were definitely moments when fixing one thing seemed to break something else (which was REALLY fun 😭), but eventually I managed to get the project built and deployed!

The science behind FoldCheck



These resources were useful throughout development:





AlphaFold Protein Structure Database



RCSB Protein Data Bank



UniProt



Mol* Molecular Visualisation





WHAT IS ALPHAFOLD

AlphaFold is an AI system that predicts protein structures from amino-acid sequences.

It has made predicted structures available for enormous numbers of proteins, but predictions aren't automatically equivalent to experimentally determined structures.

That's why investigating confidence and comparing different forms of evidence is important.

What's the difference between pLDDT and PAE?

pLDDT measures confidence in the predicted local structure, while PAE estimates uncertainty in the relative positioning of residues.

A protein can have regions with high local confidence while still having uncertainty about how those regions are arranged relative to each other.

Both measures are useful, but neither independently proves that a prediction is correct.



IMPORTANT NOTE



FoldCheck is a completed Stardance project, but that doesn't mean it can scientifically prove whether a protein prediction is correct.

Its purpose is to help users investigate protein structures, confidence scores and scientific evidence.

It should not be used as the sole basis for scientific, clinical or medical decisions.

Scientific resources





FoldCheck is an independent student project and isn't officially affiliated with organisations.

Copyright and usage

© 2026 Yousef Khanfar. All rights reserved.

FoldCheck is my student-developed project. The source code is publicly available for viewing and review, while third-party libraries and resources remain subject to their respective licences.



FINAL THOUGHTS

I wanted FoldCheck to be more than just a website that displays protein structures. I wanted to create something that helps people actually investigate what they're looking at and understand why confidence matters. I'm also planning on introducing FoldCheck to my school as a possible learning resource, especially for students interested in biology and chemistry!

It's been a pretty challenging project, especially working with molecular visualisation and scientific data, but I've learned alot from building it!





Thanks for checking out FoldCheck!!!
