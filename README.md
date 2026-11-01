# FOLDCHECK



Hi! I'm Yousef, and I made FoldCheck for Hack Club Stardance!



I've always been interested in biology, chemistry and AI, and I wanted to build something that combined all three. When I first learned about AlphaFold, I was amazed that AI could predict something as complicated as a protein's 3D structure-but I also wondered, how could scientists tell whether those predictions were actually reliable?



That gave me the idea for FoldCheck!

I wanted to create a website where people could explore protein structures, investigate how confident AI is about its predictions and learn about the science behind them.





LIVE WEBSITE



https://fold-check-seven.vercel.app



You can use FoldCheck in your browser!





WHAT DOES FOLDCHECK DO?



FoldCheck is a web application for exploring proteins and investigating AI-predicted structures.



I built several features, including:



- An interactive 3D protein viewer where you can rotate, zoom and explore structures.

- Residue and atom inspection, allowing you to select individual parts of a protein.

- pLDDT confidence analysis to investigate how confident AlphaFold is about different regions.

- PAE heatmaps to explore uncertainty in how different parts of a protein are positioned.

- Structural geometry tools for examining details such as backbone geometry.

- Functional region exploration and experimental structure comparisons.

- An animated scientific introduction for people who are unfamiliar with proteins and AlphaFold.



I used crambin (PDB ID: 1CRN) and human p53 (UniProt ID: P04637) while developing and testing FoldCheck.






HOW I STARTED



I started by working on the interface and the animated scientific introduction.



I wanted the website to be understandable even if someone didn't know much about structural biology. I didn't want people to open FoldCheck, see a complicated protein model and have absolutely no idea what they were looking at 😭



After working on the introduction, I moved on to the actual protein viewer.








THE HARDEST PART



The 3D viewer was probably my biggest challenge.



I used Mol*, which is a molecular visualisation library, but getting it to work properly with React took quite a bit of effort.



It wasn't just about displaying a protein. I wanted people to click on residues and atoms, see information about what they'd selected, highlight different parts of the structure and move the camera around.



Getting all those things to work together was difficult.



Sometimes I'd finally fix one interaction, only to discover that another feature had stopped working. That was definitely frustrating, especially when something that had been working perfectly suddenly wasn't anymore.



I spent a lot of time working on the connection between the viewer and the rest of the application, particularly selections, highlighting and camera behaviour.



But when those features finally started working together, it was VERY satisfying. It felt like FoldCheck was becoming an actual application!





BUILDING THE REST OF FOLDCHECK



Once I had the viewer working, I continued adding the tools I wanted people to use.



I worked on pLDDT confidence scores, PAE heatmaps, structural geometry, functional regions and experimental comparisons.



This part was interesting because I was learning about the science while also trying to figure out how to represent it in the application.



Before starting FoldCheck, I knew about AlphaFold, but I didn't understand nearly as much about how scientists assess the reliability of protein predictions.



I found it really interesting learning about things like pLDDT and PAE, especially the fact that a prediction can be highly confident in one way but still uncertain in another.






WHAT I ENJOYED



Even though FoldCheck was challenging, I genuinely enjoyed building it!



My favourite part was being able to combine subjects I'm interested in instead of working on something completely unrelated to what I enjoy learning.



I liked seeing the website gradually improve as I added more features, and I especially enjoyed being able to interact with actual protein structures. There was something really satisfying about going from an idea in my head to a working website that I could open and use.



I also enjoyed learning things that I wouldn't normally cover in school yet, particularly the more complicated parts of protein structures and molecular visualisation. Finally getting FoldCheck deployed after all the problems I ran into was a pretty good feeling!






WHAT I LEARNED



I built FoldCheck using React, JavaScript, Vite, HTML and CSS, with Mol* handling the 3D molecular viewer.



I learned a lot about connecting different parts of an application, managing interactions between React and another library, and dealing with problems that appeared when I changed existing features. I also learned that making something work once doesn't necessarily mean it's finished. Sometimes I had to revisit parts of the application several times before everything behaved the way I wanted.



On the science side, I learned much more about proteins, AlphaFold, confidence scores and how predicted structures can be investigated.






RESOURCES I USED



These resources helped me with the scientific information and molecular visualisation:



- [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/)

- [RCSB Protein Data Bank](https://www.rcsb.org/)

- [UniProt](https://www.uniprot.org/)

- [Mol* Molecular Visualisation](https://molstar.org/)





IMPORTANT NOTE



FoldCheck is a completed Stardance project, but it cannot prove whether a protein prediction is scientifically correct.



It's intended for exploring protein structures and investigating confidence and scientific evidence, not for making medical or clinical decisions.



FoldCheck is an independent student project and isn't officially affiliated with the organisations listed above.





FINAL THOUGHTS



FoldCheck has been one of the most challenging projects I've worked on, but also one of the most enjoyable.



I started with an interest in biology, chemistry and AI, and ended up learning so much more about both programming and protein structures than I expected.



I'm really happy that I managed to turn my idea into a working website, even if getting there involved quite a few frustrating moments!



I'm also planning on introducing FoldCheck to my school as a possible learning resource, especially for students interested in biology and chemistry.



I'd love to keep improving it in the future, but for now I'm happy with how far it's come.



Thanks for checking out FoldCheck!!!
