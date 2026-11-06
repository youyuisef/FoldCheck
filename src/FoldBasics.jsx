import {useEffect,useRef,useState} from "react";
import "./FoldStyle.css";
import ProteinViewer from "./ProteinViewer";
import aminoZoo from "./AmazingAminos";
import atomLore from "./BomboclatAtoms";
import RegionDetails from "./RegionDetails";
import PAEHeatmap from "./PAEHeatmap";

const aminoPick = Object.entries(aminoZoo).map(([code,amino])=>[
  code,
  amino.name,
]);

const topDoors = [
  ["finder","Find a protein"],
  ["basics","Basics"],
  ["limitations","Limitations"],
  ["about","About Me"],
  ["sources","Sources"],
];

const viewTabs = [
  ["workspace","3D Structure"],
  ["confidence","Confidence"],
  ["pae","PAE"],
  ["geometry","Geometry"],
  ["experimental","Experimental Comparison"],
];

const critterChips = [
  ["Human","organism_id:9606","human"],
  ["Mouse","organism_id:10090","mouse"],
  ["Yeast","organism_id:559292","yeast"],
  ["E. coli","organism_id:83333","E. coli"],
  ["Plant","organism_id:3702","plant"],
  ["Any","",""],
];

const jobChips = [
  ["Carries oxygen","keyword:KW-0561","helps carry oxygen around"],
  ["Speeds up reactions","ec:*","speeds up chemical reactions"],
  ["Fights germs","keyword:KW-0391","helps defend against germs"],
  ["Sends messages","keyword:KW-0372","sends messages between cells"],
  ["Stops cancer","keyword:KW-0043","helps stop cells growing out of control"],
  ["Moves muscles","keyword:KW-0514","helps muscles move"],
  ["Gives cells shape","keyword:KW-0206","helps cells keep their shape"],
  ["Switches genes on and off","keyword:KW-0805","helps switch genes on and off"],
  ["Not sure","",""],
];

const termBits = [
  ["Amino acids","Proteins are built from 20 commonly used amino acids. Each has an amino group,carboxyl group and a side chain that affects its chemical behaviour."],
  ["Alpha carbon (Cα)","The central carbon of each amino acid,bonded to the backbone nitrogen,the carbonyl carbon,a hydrogen and the side chain. FoldCheck uses one Cα per residue as that residue's position when it compares structures."],
  ["Alpha helices","A common secondary structure where the backbone coils into a helix. Hydrogen bonds between backbone groups help stabilise it."],
  ["AlphaFold model","AlphaFold predicts a possible 3D protein structure using AI. A prediction is not direct experimental proof,and proteins can move or change shape."],
  ["Backbone","The repeating N–Cα–C chain running through every residue,with the carbonyl oxygen attached. Side chains hang off it,and phi and psi angles describe how it twists."],
  ["Beta sheets","A secondary structure formed when extended stretches of the protein backbone align and form hydrogen bonds."],
  ["Chains","In a structure file,each polypeptide gets a letter (chain A,B…). Experimental entries often hold several copies or partners,so FoldCheck picks the chain that best matches the AlphaFold sequence before comparing."],
  ["Confidence scores","A confidence score indicates how reliable a model is expected to be in a particular respect. It is not the probability that the whole protein is correct."],
  ["Disordered regions","Some regions do not settle into one stable 3D shape. Low model confidence can reflect genuine flexibility or disorder,not necessarily a modelling failure."],
  ["Disulfide bonds","Covalent links between cysteine side chains that can help stabilise a protein,especially outside cells."],
  ["Domains","A domain is a section of a protein that can fold into a distinct structural unit and often carries a particular function."],
  ["Enzymes","Proteins that speed up chemical reactions by lowering activation energy. Their active sites bind substrates and help reactions happen."],
  ["Evidence codes","UniProt tags each annotation with how it's supported: experimental evidence in a publication,a curator's inference,or similarity to another protein. The Region Explorer shows these,and experimental evidence is usually the strongest."],
  ["Experimental structures","Scientists determine structures using methods such as X-ray crystallography,cryo-electron microscopy and NMR. Experimental models also have limitations and may only show part of a protein."],
  ["Folding","The process by which an amino-acid chain adopts its three-dimensional shape. The sequence,environment and sometimes helper proteins influence folding."],
  ["Hydrogen bonds","Attractions involving hydrogen and electronegative atoms. Many hydrogen bonds help stabilise alpha helices,beta sheets and other protein interactions."],
  ["Hydrophobic effect","Nonpolar side chains often pack away from water,helping many soluble proteins form a stable interior."],
  ["Ionic interactions","Attractions between oppositely charged groups,such as certain amino-acid side chains. Their strength depends on the surrounding environment."],
  ["Ligands","Molecules or ions that bind to a protein. Some help proteins work,while others inhibit or change their activity."],
  ["Membrane proteins","Proteins embedded in or associated with cell membranes. They can transport substances,detect signals or connect cells."],
  ["Mutations","Changes in DNA can alter a protein sequence or how much protein is produced. Effects range from harmless to important,depending on the change and context."],
  ["Organelles and localisation","Where a protein works matters. Proteins may function in the nucleus,cytoplasm,mitochondria,membranes or other compartments,and some move between locations."],
  ["PAE","Predicted aligned error estimates uncertainty in the relative positions of residues or domains. It helps assess whether the arrangement of two regions is reliable."],
  ["PDB","The Protein Data Bank,an archive of experimentally determined structures. Each entry has a four-character ID,like 1CRN,the demo structure on the homepage. UniProt links proteins to their PDB entries."],
  ["Peptide bonds","Covalent bonds linking amino acids into a polypeptide chain. They form between the carboxyl group of one amino acid and the amino group of another."],
  ["Phi and psi angles","Two rotations per residue along the backbone: phi (φ) around the N–Cα bond and psi (ψ) around the Cα–C bond. Together they set the local shape. The first residue has no phi and the last has no psi,which is why the geometry table shows N/A at chain ends."],
  ["pLDDT","AlphaFold’s local confidence score,from 0 to 100. Above 90 is very high confidence:70–90 confident:50–70 low:and below 50 very low."],
  ["Primary structure","The order of amino acids in a protein,usually written from the N-terminus to the C-terminus. This sequence influences how the protein folds."],
  ["Protein function","A protein’s job may include catalysing reactions,carrying molecules,providing structure,sending signals or regulating genes. Structure helps,but does not alone prove function."],
  ["Quaternary structure","The arrangement of multiple protein chains into one functional complex. Many proteins work as single chains and do not have quaternary structure."],
  ["Ramachandran plot","A scatter plot of each residue's phi (across) against psi (up). Helix residues cluster near φ ≈ −60°,ψ ≈ −45°,and strand residues near φ ≈ −120°,ψ ≈+130°,because most other combinations make atoms clash. Glycine is more flexible and proline more restricted. FoldCheck draws the points but doesn't yet mark favoured or outlier zones."],
  ["Residues","An amino acid is called a residue once incorporated into a protein chain. Residue numbers help locate positions in sequences and structures."],
  ["Resolution","How finely an experimental structure is measured,in Å. Smaller is sharper. Around 2 Å or better,side chains are usually clear,and around 3 Å or worse their positions are less certain. NMR entries usually have no resolution value."],
  ["RMSD","Root-mean-square deviation,the typical distance between matching atoms in two structures after they've been optimally rotated and shifted onto each other. FoldCheck uses Cα atoms and shows it in Å. Lower means closer,but it's an average,so one floppy loop can raise it a lot."],
  ["Secondary structure","Local patterns of the protein backbone,especially alpha helices and beta sheets,stabilised largely by hydrogen bonds."],
  ["Sequence alignment","A method of matching amino-acid positions between sequences or structures. It helps compare related proteins and identify substitutions or missing regions."],
  ["Side chains","The variable parts of amino acids. Their size,charge and polarity affect interactions and protein behaviour."],
  ["Signal peptides","A short stretch at the start of some proteins that directs them to the secretory pathway or certain compartments,and is usually cut off afterwards. UniProt marks these as \"Signal\",so the mature protein may lack those first residues."],
  ["Superposition","Rotating and shifting one structure to sit as closely as possible over another,so the remaining differences aren't just orientation. FoldCheck does this mathematically before calculating RMSD. The two 3D viewers on screen aren't superimposed."],
  ["Termini","The two ends of a chain. The N-terminus has a free amino group and the C-terminus a free carboxyl group,which is where the extra OXT oxygen in the Atom Explorer comes from. Residues are numbered from N to C."],
  ["Tertiary structure","The overall three-dimensional arrangement of one protein chain,held together by many interactions including hydrophobic packing,hydrogen bonds and sometimes disulfide bonds."],
  ["UniProt","A database of protein sequences and biological annotations. Its descriptions and region information are distinct from AlphaFold’s predicted coordinates."],
];

const caveatBits = [
  ["AlphaFold confidence (pLDDT)","pLDDT is AlphaFold's own estimate of its confidence in the local structure around a residue. It isn't a measured accuracy. A high score doesn't guarantee correct side chains or interactions,and a low score doesn't prove disorder."],
  ["Predicted aligned error (PAE)","PAE is AlphaFold's predicted uncertainty in the relative position of two residues,in Å. It isn't a measured distance or experimental evidence. The protein-wide averages include lots of neighbouring pairs,which pull them down,so being below the average doesn't mean a pair is well predicted."],
  ["Predicted structures","An AlphaFold model is one predicted shape,not an observation. Real proteins move and work with partners or ligands that the model may not show."],
  ["Residue-level interpretation","A residue's pLDDT says how much AlphaFold trusts the local prediction,not whether the residue is important. What a residue does depends on its neighbours in 3D,so don't read exact contacts from low-confidence regions or claim a binding site from pLDDT."],
  ["Atom-level interpretation","Atom explanations describe the element and where the atom sits in the amino acid. They don't show that a particular hydrogen bond or charge interaction happens here,and atom positions are shakier where pLDDT is low."],
  ["Functional region annotations","Regions come from UniProt annotations. They describe biology but don't verify binding or the accuracy of the displayed coordinates. No annotation doesn't mean no function,and the evidence can be experimental,a curator's inference or similarity to another protein."],
  ["Experimental structure comparisons","A PDB entry may cover only part of the protein,include other molecules,and reflects its own experimental conditions. The two viewers aren't superimposed on screen. FoldCheck matches residues by sequence and calculates a Cα RMSD,but check the chain and matching before treating it as proof of agreement. A low RMSD doesn't validate the prediction,and a high one doesn't disprove it."],
];

const pctBit = (part,total)=>(total?(part/total)*100 : 0);

const pickTitle = (protein)=>
  protein.proteinDescription?.recommendedName?.fullName?.value||
  protein.proteinDescription?.submissionNames?.[0]?.fullName?.value||
  protein.primaryAccession;

const funcBlurb = (comments = [])=>
  comments
    .filter((comment)=>comment.commentType==="FUNCTION")
    .flatMap((comment)=>comment.texts||[])
    .map((text)=>text.value)
    .join(" ")||"";

function lineUpSeq(predicted,experimental) {
  const n = predicted.length;
  const m = experimental.length;
  if (n===0||m===0) return [];

  const MATCH = 2;
  const MISMATCH  = -1;
  const GAP =-2;

  const scores = Array.from({length: n+1},()=>new Int32Array(m+1));
  for (let i = 1; i<=n; i++) scores[i][0] = i*GAP;
  for (let j = 1; j<=m; j++) scores[0][j] = j*GAP;

  for (let i = 1; i<=n; i++) {
    for (let j = 1; j<=m; j++) {
      const diagonal =
        scores[i-1][j-1]+
        (predicted[i-1].code===experimental[j-1].code?MATCH : MISMATCH);
      const up = scores[i-1][j]+GAP;
      const left = scores[i][j-1]+GAP;
      scores[i][j] = Math.max(diagonal,up,left);
}
}

  const matches = [];
  let i  =  n;
  let j  =  m;

  while (i>0&&j>0) {
    const sameCode = predicted[i-1].code===experimental[j-1].code;
    const diagonal = scores[i-1][j-1]+(sameCode?MATCH : MISMATCH);

    if (sameCode&&scores[i][j]===diagonal) {
      matches.push({
        predicted: predicted[i-1],
        experimental: experimental[j-1],
});
}

    if (scores[i][j]===diagonal) {
      i--;
      j--;
} else if (scores[i][j]===scores[i-1][j]+GAP) {
      i--;
} else {
      j--;
}
}

  return matches.reverse();
}

function superRmsd(pairs) {
  if (pairs.length<3) return null;

  const n = pairs.length;

  const mean = (key)=>
    [0,1,2].map(
      (axis)=>pairs.reduce((sum,pair)=>sum+pair[key][axis],0)/n
);

  const predictedCenter  =  mean("predicted");
  const experimentalCenter = mean("experimental");

  const covariance = Array.from({length: 3},()=>[0,0,0]);

  for (const pair of pairs) {
    const p = pair.predicted.map((v,i)=>v-predictedCenter[i]);
    const q = pair.experimental.map((v,i)=>v-experimentalCenter[i]);

    for (let i = 0; i<3; i++) {
      for (let j = 0; j<3; j++) {
        covariance[i][j]+= p[i]*q[j];
}
}
}

  const [
    [sxx,sxy,sxz],
    [syx,syy,syz],
    [szx,szy,szz],
] = covariance;

  const K = [
    [sxx+syy+szz,syz-szy,szx-sxz,sxy-syx],
    [syz-szy,sxx-syy-szz,sxy+syx,szx+sxz],
    [szx-sxz,sxy+syx,-sxx+syy-szz,syz+szy],
    [sxy-syx,szx+sxz,syz+szy,-sxx-syy+szz],
];

  
  const A = K.map((row)=>[...row]);
  const V = Array.from({length: 4},(_,i)=>
    Array.from({length: 4},(_,j)=>(i===j?1 : 0))
);

  for (let iteration = 0; iteration<100; iteration++) {
    let a  =  0;
    let b = 1;
    let largest = 0;

    for (let i  =  0; i<4; i++) {
      for (let j = i+1; j<4; j++) {
        if (Math.abs(A[i][j])>largest) {
          largest = Math.abs(A[i][j]);
          a = i;
          b = j;
}
}
}

    if (largest<1e-12) break;

    const angle = 0.5*Math.atan2(2*A[a][b],A[b][b]-A[a][a]);
    const c = Math.cos(angle);
    const s = Math.sin(angle);

    for (let i = 0; i<4; i++) {
      if (i===a||i===b) continue;
      const x = A[i][a];
      const y  =  A[i][b];
      A[i][a] = A[a][i] = c*x-s*y;
      A[i][b] = A[b][i] = s*x+c*y;
}

    const aa = A[a][a];
    const bb = A[b][b];
    const ab = A[a][b];

    A[a][a] = c*c*aa-2*c*s*ab+s*s*bb;
    A[b][b] = s*s*aa+2*c*s*ab+c*c*bb;
    A[a][b] = A[b][a] = 0;

    for (let i = 0; i<4; i++) {
      const x = V[i][a];
      const y = V[i][b];
      V[i][a]  =  c*x-s*y;
      V[i][b] = s*x+c*y;
}
}

  const eigenvalues = [0,1,2,3].map((i)=>A[i][i]);
  const best = eigenvalues.indexOf(Math.max(...eigenvalues));
  const [w,x,y,z] = V.map((row)=>row[best]);

  const rotation = [
    [1-2*(y*y+z*z),2*(x*y-w*z),2*(x*z+w*y)],
    [2*(x*y+w*z),1-2*(x*x+z*z),2*(y*z-w*x)],
    [2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y)],
];

  let squaredError = 0;

  for (const pair of pairs) {
    const p = pair.predicted.map((v,i)=>v-predictedCenter[i]);
    const q = pair.experimental.map((v,i)=>v-experimentalCenter[i]);

    const rotated = rotation.map((row)=>
      row.reduce((sum,value,i)=>sum+value*p[i],0)
);

    squaredError+= rotated.reduce(
      (sum,value,i)=>sum+(value-q[i])**2,
      0
);
}

  return Math.sqrt(squaredError/n);
}

function phiPsiTwist(p1,p2,p3,p4) {
  const subtract = (a,b)=>a.map((v,i)=>v-b[i]);
  const dot = (a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
  const cross = (a,b)=>[
    a[1]*b[2]-a[2]*b[1],
    a[2]*b[0]-a[0]*b[2],
    a[0]*b[1]-a[1]*b[0],
];

  const b0 = subtract(p1,p2);
  const b1 = subtract(p3,p2);
  const b2 = subtract(p4,p3);

  const length = Math.sqrt(dot(b1,b1));
  if (length===0) return null;

  const axis = b1.map((v)=>v/length);
  const v = subtract(b0,axis.map((n)=>dot(b0,axis)*n));
  const w = subtract(b2,axis.map((n)=>dot(b2,axis)*n));

  const angle = Math.atan2(dot(cross(axis,v),w),dot(v,w));
  return (angle*180)/Math.PI;
}

function TopStrip({setActivePage}) {
  return (
<header className="nav">
<button className="brand brandButton" onClick={()=>setActivePage("home")}>
<span>FoldCheck</span>
</button>

<div className="navLinks">
        {topDoors.map(([page,label])=>(
<button
            key={page}
            className="navTextButton"
            onClick={()=>setActivePage(page)}
>
            {label}
</button>
))}
</div>

<button className="navButton" onClick={()=>setActivePage("home")}>
        Explore a protein
</button>
</header>
);
}

function LensStrip({activePage,setActivePage}) {
  return (
<nav className = "analysisTabs">
      {viewTabs.map(([page,label])=>(
<button
          key={page}
          className={`analysisTab ${activePage === page ? "active" : ""}`}
          onClick={()=>setActivePage(page)}
>
          {label}
</button>
))}
</nav>
);
}

function FoldBasics() {

  const [geomRows,setGeomRows] = useState([]);
  const [labRows,setLabRows] = useState([]);
  const [linkedSelection,setLinkedSelection] = useState(null);

  const experimentalChains = [...new Set(labRows.map((r)=>r.chain))];

  const predictedBackbone = geomRows.filter((r)=>r.backbone?.CA);

  const selectedExperimentalChain =
    experimentalChains
      .map((chain)=>{
        const chainResidues = labRows.filter(
          (r)=>r.chain===chain&&r.backbone?.CA
);
        const matches = lineUpSeq(predictedBackbone,chainResidues);

        return {
          chain,
          matchedCount: matches.length,
          coverage: chainResidues.length?matches.length/chainResidues.length : 0,
};
})
      .sort((a,b)=>b.matchedCount-a.matchedCount||b.coverage-a.coverage)[0]
?.chain??null;

  const asSeqEntry = (residue)=>({
    number: residue.number,
    code: residue.code,
    chain: residue.chain,
    CA: residue.backbone.CA,
});

  const predictedSequence = predictedBackbone.map(asSeqEntry);

  const experimentalSequence = labRows
    .filter((r)=>r.chain===selectedExperimentalChain&&r.backbone?.CA)
    .map(asSeqEntry);

  const pairedResidues  =  lineUpSeq(
    predictedSequence,
    experimentalSequence
).map(({predicted,experimental})=>({
    number: predicted.number,
    code: predicted.code,
    predictedCA: predicted.CA,
    experimentalCA: experimental.CA,
    predictedNumber: predicted.number,
    predictedChain: predicted.chain,
    experimentalNumber: experimental.number,
    experimentalChain: experimental.chain,
}));

  const mirrorPick = (side,residue)=>{
    const number = Number(residue?.residueNumber);
    const numberKey = side==="alphafold"?"predictedNumber" : "experimentalNumber";
    const chainKey  =  side==="alphafold"?"predictedChain" : "experimentalChain";

    const match = pairedResidues.find(
      (item)=>
        item[numberKey]===number&&
        (!residue.chain||item[chainKey]===residue.chain)
);

    if (!match) {
      setLinkedSelection(null);
      return;
}

    setLinkedSelection({
      alphafold: {number: match.predictedNumber,chain: match.predictedChain},
      pdb: {number: match.experimentalNumber,chain: match.experimentalChain},
});
};

  const rmsdPairs = pairedResidues
    .filter(
      (residue)=>
        residue.predictedCA?.length===3&&
        residue.experimentalCA?.length===3&&
        [...residue.predictedCA,...residue.experimentalCA].every(Number.isFinite)
)
    .map((residue)=>({
      predicted: residue.predictedCA,
      experimental: residue.experimentalCA,
}));

  const superposedRmsd =
    rmsdPairs.length>=3?superRmsd(rmsdPairs) : null;

  const [showIntro,setShowIntro] = useState(true);
  const [introStage,setIntroStage] = useState(0);
  const [pageNow,setPageNow] = useState("home");
  const [glossaryQuery,setGlossaryQuery] = useState("");
  const [glossaryLetter,setGlossaryLetter] = useState("All");

  useEffect(()=>{
    window.scrollTo(0,0);
},[pageNow]);

  useEffect(()=>{
    const title = document.querySelector(".newHomeHero h1");
    if (!title||showIntro||window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
      return;
}

    const blobs = [
      {x: 0,y: 0,vx: 0,vy: 0,fx: 0.00095,fy: 0.00132,px: 0.0,py: 1.3},
      {x: 0,y: 0,vx: 0,vy: 0,fx: 0.00078,fy: 0.00104,px: 2.1,py: 0.4},
      {x: 0,y: 0,vx: 0,vy: 0,fx: 0.00112,fy: 0.00086,px: 4.2,py: 2.6},
];
    let pointer  =  null;
    let frame = 0;
    let ready = false;

    const onMove = (event)=>{
      const box = title.getBoundingClientRect();
      pointer  =  {x: event.clientX-box.left,y: event.clientY-box.top};
};
    const onLeave = ()=>{
      pointer = null;
};

    title.addEventListener("mousemove",onMove);
    title.addEventListener("mouseleave",onLeave);

    const tick = (time)=>{
      const w = title.offsetWidth;
      const h  =  title.offsetHeight;

      blobs.forEach((blob,index)=>{
        const targetX = w*(0.5+0.42*Math.sin(time*blob.fx+blob.px));
        const targetY  =  h*(0.5+0.55*Math.sin(time*blob.fy+blob.py));

        if (!ready) {
          blob.x = targetX;
          blob.y = targetY;
}

        blob.vx+= (targetX-blob.x)*0.0045;
        blob.vy+= (targetY-blob.y)*0.0045;

        if (pointer) {
          const dx = blob.x-pointer.x;
          const dy = blob.y-pointer.y;
          const dist = Math.hypot(dx,dy)+1;

          if (dist<260) {
            const push = (1-dist/260)*2.2;
            blob.vx+= (dx/dist)*push;
            blob.vy+= (dy/dist)*push*0.6;
}
}

        blob.vx *= 0.93;
        blob.vy *= 0.93;
        blob.x+= blob.vx;
        blob.y+= blob.vy;

        title.style.setProperty(`--b${index}x`, `${(blob.x / w) * 100}%`);
        title.style.setProperty(`--b${index}y`, `${(blob.y / h) * 100}%`);
});

      ready = true;
      frame  =  requestAnimationFrame(tick);
};

    frame = requestAnimationFrame(tick);

    return ()=>{
      cancelAnimationFrame(frame);
      title.removeEventListener("mousemove",onMove);
      title.removeEventListener("mouseleave",onLeave);
};
},[pageNow,showIntro]);

  const [searchBoxText,setHomeQuery] = useState("");
  const [proteinSearchLoading,setProteinSearchLoading] = useState(false);
  const [proteinSearchError,setProteinSearchError] = useState("");
  const [searchHits,setHomeHits]  =  useState([]);
  const [showProteinSuggestions,setShowProteinSuggestions] = useState(false);
  const [activeProteinSuggestionIndex,setActiveProteinSuggestionIndex] = useState(-1);

  useEffect(()=>{
    const query = searchBoxText.trim();

    if (query.length<2) {
      setHomeHits([]);
      return;
}

    const controller  =  new AbortController();

    const timer = setTimeout(async ()=>{
      try {
        const response = await fetch(
          `https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(
            `(${query}) AND reviewed:true`
          )}&fields=accession,protein_name,gene_names,organism_name&size=8`,
          {signal: controller.signal}
);

        if (!response.ok) return;

        const data = await response.json();

        if (!controller.signal.aborted) {
          setHomeHits(data.results||[]);
}
} catch (error) {
        if (error.name!=="AbortError") {
          setHomeHits([]);
}
}
},350);

    return ()=>{
      clearTimeout(timer);
      controller.abort();
};
},[searchBoxText]);

  const visitProtein = async (accession)=>{
    setProteinSearchLoading(true);
    setProteinSearchError("");

    try {
      const response = await fetch(
        `https://alphafold.ebi.ac.uk/api/prediction/${encodeURIComponent(accession)}`
);

      if (response.status===404) {
        setProteinSearchError(
          "This protein exists in UniProt,but no AlphaFold prediction is available."
);
        return;
}

      if (!response.ok) {
        throw new Error("AlphaFold availability check failed");
}

      setShowProteinSuggestions(false);
      setActiveProteinSuggestionIndex(-1);
      setSelRes(null);
      setSideTab("residue");
      setActiveAcc(accession);
      setPageNow("workspace");
} catch {
      setProteinSearchError(
        "Could not check AlphaFold availability. Please try again."
);
} finally {
      setProteinSearchLoading(false);
}
};

  const lookUpProtein = async ()=>{
    const query = searchBoxText.trim();
    if (!query||proteinSearchLoading) return;

    setProteinSearchLoading(true);
    setProteinSearchError("");

    try {
      const accession = query.toUpperCase();

      const response = await fetch(
        `https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(
          `(${query}) AND reviewed:true`
        )}&fields=accession,protein_name,gene_names,organism_name&size=10`
);

      if (!response.ok) throw new Error("UniProt search failed");

      const data = await response.json();

      const exactMatch = data.results?.find(
        (protein)=>protein.primaryAccession===accession
);
      const protein = exactMatch||data.results?.[0];

      if (!protein) {
        setProteinSearchError("No matching protein found in UniProt.");
        return;
}

      await visitProtein(protein.primaryAccession);
} catch {
      setProteinSearchError("Could not search UniProt. Please try again.");
} finally {
      setProteinSearchLoading(false);
}
};

  const [finderSpecies,setFinderSpecies] = useState(0);
  const [finderRole,setFinderRole] = useState(null);
  const [finderMore,setFinderMore]  =  useState(false);
  const [finderWord,setFinderWord] = useState("");
  const [finderSize,setFinderSize] = useState("any");
  const [finderStructure,setFinderStructure] = useState(false);
  const [finderPicks,setFinderPicks] = useState(null);
  const [finderLoading,setFinderLoading] = useState(false);
  const [finderError,setFinderError] = useState("");

  const findProteins = async ()=>{
    if (finderRole===null||finderLoading) return;

    setFinderLoading(true);
    setFinderError("");
    setFinderPicks(null);
    setProteinSearchError("");

    const parts  =  ["reviewed:true"];
    if (critterChips[finderSpecies][1]) parts.push(critterChips[finderSpecies][1]);
    if (jobChips[finderRole][1]) parts.push(jobChips[finderRole][1]);
    if (finderWord.trim()) parts.push(`(${finderWord.trim()})`);
    if (finderStructure) parts.push("structure_3d:true");
    if (finderSize==="small") parts.push("length:[1 TO 200]");
    if (finderSize==="large") parts.push("length:[500 TO*]");

    try {
      const response  =  await fetch(
        `https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(
          parts.join(" AND ")
)}&fields=accession,protein_name,gene_names,organism_name&sort=${encodeURIComponent(
          "annotation_score desc"
        )}&size=4`
);

      if (!response.ok) throw new Error("finder search failed");

      const data = await response.json();

      const found = (data.results||[]).map((protein)=>({
        id: protein.primaryAccession,
        name: pickTitle(protein),
        species: protein.organism?.commonName||protein.organism?.scientificName||"",
}));

      if (!found.length) {
        setFinderError("Nothing matched those answers. Try fewer options or pick Any.");
}
      setFinderPicks(found);
} catch {
      setFinderError("Could not search UniProt. Please try again.");
} finally {
      setFinderLoading(false);
}
};

  const [activeAcc,setActiveAcc] = useState("P04637");
  const [protCard,setProtCard] = useState(null);
  const [proteinInfoLoading,setProteinInfoLoading] = useState(false);
  const [proteinInfoError,setProteinInfoError]  =  useState("");
  const [pdbEntries,setPdbEntries] = useState([]);
  const [experimentalLoading,setExperimentalLoading] = useState(false);
  const [experimentalError,setExperimentalError]  =  useState("");
  const [selectedExperimentalId,setSelectedExperimentalId] = useState("");
  const [partnerInfo,setPartnerInfo] = useState({});
  const [selRegion,setSelRegion]  =  useState(null);
  const [selRes,setSelRes] = useState(null);
  const [modelScore,setModelScore] = useState(null);
  const [errorGrid,setErrorGrid] = useState(null);
  const [paeReferenceResidue,setPaeReferenceResidue] = useState(null);
  const [paeSelectedPair,setPaeSelectedPair] = useState(null);
  const [paeLoading,setPaeLoading] = useState(false);
  const [paeError,setPaeError]  =  useState("");
  const [panelW,setPanelW] = useState(320);
  const [showConfidenceInfo,setShowConfidenceInfo] = useState(false);
  const [sideTab,setSideTab] = useState("residue");
  const [jumpQuery,setJumpQuery] = useState("");
  const [proteinResidues,setProteinResidues] = useState([]);
  const [showStructureSuggestions,setShowStructureSuggestions] = useState(false);
  const [activeSuggestionIndex,setActiveSuggestionIndex]  =  useState(-1);
  const [goToRes,setGoToRes] = useState(null);
  const [resetCam,setResetCam]  =  useState(0);
  const [jumpFeedback,setJumpFeedback] = useState(null);

  useEffect(()=>{
    setGeomRows([]);
},[activeAcc]);

  useEffect(()=>{
    const partners = protCard?.interactionPartners||[];
    const regions = protCard?.regions||[];

    const partnerIds  =  [
      ...new Set(partners.filter((partner)=>partner.id).map((partner)=>partner.id)),
];

    const regionPartnerNames = regions
      .map((region)=>region.name.match(/^Interaction with\s+(.+)$/i)?.[1]?.trim())
      .filter(Boolean);

    const missingPartnerNames = regionPartnerNames.filter(
      (name)=>
        !partners.some((partner)=>partner.name.toLowerCase()===name.toLowerCase())
);

    const controller  =  new AbortController();

    const loadPartners = async ()=>{
      const matchingPartners = partnerIds.filter((id)=>
        regions.some((region)=>
          partners.some(
            (partner)=>
              partner.id===id&&
              partner.name.length>2&&
              region.name.toLowerCase().includes(partner.name.toLowerCase())
)
)
);

      const missingPartnerResults = await Promise.all(
        missingPartnerNames.slice(0,10).map(async (name)=>{
          try {
            const query = `gene_exact:${name} AND organism_id:9606 AND reviewed:true`;

            const response = await fetch(
              `https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(
                query
              )}&size=1`,
              {signal: controller.signal}
);

            if (!response.ok) return null;

            const data = await response.json();
            const protein = data.results?.[0];
            if (!protein) return null;

            return {
              id: protein.primaryAccession,
              gene: name,
              name: protein.proteinDescription?.recommendedName?.fullName?.value||name,
              function: funcBlurb(protein.comments),
};
} catch {
            return null;
}
})
);

      const results = await Promise.all(
        matchingPartners.slice(0,10).map(async (id)=>{
          try {
            const response = await fetch(
              `https://rest.uniprot.org/uniprotkb/${encodeURIComponent(id)}.json`,
              {signal: controller.signal}
);

            if (!response.ok) return null;

            const data = await response.json();

            return {
              id,
              name: data.proteinDescription?.recommendedName?.fullName?.value||id,
              function: funcBlurb(data.comments),
};
} catch {
            return null;
}
})
);

      if (controller.signal.aborted) return;

      setPartnerInfo(
        Object.fromEntries(
          [...results,...missingPartnerResults]
            .filter(Boolean)
            .flatMap((partner)=>[
              [partner.id,partner],
              ...(partner.gene?[[partner.gene.toUpperCase(),partner]] : []),
])
)
);
};

    setPartnerInfo({});
    loadPartners();

    return ()=>controller.abort();
},[protCard]);

  useEffect(()=>{
    const controller = new AbortController();

    const loadPAE = async ()=>{
      setErrorGrid(null);
      setPaeSelectedPair(null);
      setPaeReferenceResidue(null);
      setPaeError("");
      setPaeLoading(true);

      try {
        const response = await fetch(
          `https://alphafold.ebi.ac.uk/files/AF-${activeAcc}-F1-predicted_aligned_error_v6.json`,
          {signal: controller.signal}
);

        if (!response.ok) {
          throw new Error("PAE data unavailable");
}

        const data = await response.json();
        const entry = Array.isArray(data)?data[0] : data;

        let matrix  =  entry?.predicted_aligned_error;

        
        if (!matrix&&Array.isArray(entry?.residue1)&&Array.isArray(entry?.distance)) {
          const size = Math.max(...entry.residue1,...entry.residue2);
          matrix = Array.from({length: size},()=>Array(size).fill(30));

          for (let i  =  0; i<entry.distance.length; i++) {
            const row = entry.residue1[i]-1;
            const col  =  entry.residue2[i]-1;

            if (row>=0&&col>=0&&row<size&&col<size) {
              matrix[row][col]  =  entry.distance[i];
}
}
}

        if (!Array.isArray(matrix)||!Array.isArray(matrix[0])) {
          throw new Error("Invalid PAE data");
}

        if (!controller.signal.aborted) {
          setErrorGrid(matrix);
}
} catch {
        if (!controller.signal.aborted) {
          setPaeError("PAE data could not be loaded for this protein.");
}
} finally {
        if (!controller.signal.aborted) {
          setPaeLoading(false);
}
}
};

    loadPAE();

    return ()=>controller.abort();
},[activeAcc]);

  useEffect(()=>{
    const controller = new AbortController();

    const loadProteinInfo = async ()=>{
      setProtCard(null);
      setProteinInfoError("");
      setProteinInfoLoading(true);
      setPdbEntries([]);
      setSelectedExperimentalId("");
      setExperimentalError("");
      setExperimentalLoading(true);

      try {
        const response = await fetch(
          `https://rest.uniprot.org/uniprotkb/${encodeURIComponent(activeAcc)}.json`,
          {signal: controller.signal}
);

        if (!response.ok) {
          throw new Error("Protein information unavailable");
}

        const data = await response.json();
        if (controller.signal.aborted) return;

        const propertyValue = (reference,key)=>
          reference.properties?.find((property)=>property.key===key)?.value||
          "Not specified";

        const pdbReferences = (data.uniProtKBCrossReferences||[])
          .filter((reference)=>reference.database==="PDB")
          .map((reference)=>({
            id: reference.id,
            method: propertyValue(reference,"Method"),
            resolution: propertyValue(reference,"Resolution"),
            chains: propertyValue(reference,"Chains"),
}));

        setPdbEntries(pdbReferences);
        setSelectedExperimentalId(pdbReferences[0]?.id||"");

        const comments = data.comments||[];
        const interactionComments = comments.filter(
          (comment)=>comment.commentType==="INTERACTION"
);

        setProtCard({
          name:
            data.proteinDescription?.recommendedName?.fullName?.value||
            data.proteinDescription?.submissionNames?.[0]?.fullName?.value||
            activeAcc,
          gene: data.genes?.[0]?.geneName?.value||"",
          aliases: [
            ...(data.genes?.[0]?.synonyms||[]).map((item)=>item.value),
            ...(data.proteinDescription?.recommendedName?.shortNames||[]).map(
              (item)=>item.value
),
].filter(Boolean),
          organism: data.organism?.scientificName||"",
          organismCommonName: data.organism?.commonName||"",
          length: data.sequence?.length||null,
          sequence: data.sequence?.value||"",
          subcellularLocations: comments
            .filter((comment)=>comment.commentType==="SUBCELLULAR LOCATION")
            .flatMap((comment)=>comment.subcellularLocations||[])
            .map((item)=>({
              name: item.location?.value||"Unknown",
              topology: item.topology?.value||"",
              orientation: item.orientation?.value||"",
}))
            .filter((item)=>item.name),
          function: comments
            .find((item)=>item.commentType==="FUNCTION")
?.texts?.map((item)=>item.value)
            .join(" ")
            .replace(/\s*\(PubMed:[^)]+\)/g,"")
            .replace(/\s+/g," "),
          interactionPartners: interactionComments
            .flatMap((comment)=>comment.interactions||[])
            .map((interaction)=>({
              name: interaction.interactantTwo?.geneName||"",
              id: interaction.interactantTwo?.uniProtKBAccession||"",
              experiments: interaction.numberOfExperiments??null,
              intActId: interaction.interactantTwo?.intActId||"",
              organismDiffer: interaction.organismDiffer??false,
})),
          interactionComments: interactionComments.flatMap(
            (comment)=>comment.interactions||[]
),
          regions: (data.features||[])
            .filter((feature)=>
              ["Domain","Region","Chain","Peptide","Propeptide","Signal"].includes(
                feature.type
)
)
            .map((feature)=>({
              name: feature.description||feature.type,
              type: feature.type,
              description: [
                feature.description,
                ...(feature.ligand
                  ? [`Ligand: ${feature.ligand.name || feature.ligand.id || "specified"}`]
                  : []),
]
                .filter(Boolean)
                .join(". "),
              start: feature.location?.start?.value,
              end: feature.location?.end?.value,
              evidence: feature.evidences||[],
})),
});
} catch (error) {
        if (error.name==="AbortError") return;

        setProteinInfoError("Could not retrieve protein information from UniProt.");
        setExperimentalError("Could not retrieve experimental structure references.");
} finally {
        if (!controller.signal.aborted) {
          setProteinInfoLoading(false);
          setExperimentalLoading(false);
}
}
};

    loadProteinInfo();

    return ()=>controller.abort();
},[activeAcc]);

  const workspaceRef = useRef(null);
  const resizingRef = useRef(false);

  const grabHandle = (event)=>{
    event.preventDefault();
    resizingRef.current = true;
    event.currentTarget.setPointerCapture(event.pointerId);
};

  const dragHandle = (event)=>{
    if (!resizingRef.current||!workspaceRef.current) return;

    const bounds = workspaceRef.current.getBoundingClientRect();
    const availableWidth  =  bounds.width-10;
    const requestedWidth = bounds.right-event.clientX;
    const maximumWidth  =  Math.max(280,availableWidth-350);

    setPanelW(Math.max(280,Math.min(requestedWidth,maximumWidth)));
};

  const dropHandle = ()=>{
    resizingRef.current = false;
};

  const structureSuggestions = (()=>{
    const query = String(jumpQuery??"").toUpperCase().trimStart();
    if (!query) return [];

    const atomMatch = query.match(/^([A-Z]{3})\s+(\d+)\s+([A-Z0-9']*)$/);
    if (atomMatch) {
      const [,code,number,atomInput] = atomMatch;
      const residue  =  proteinResidues.find(
        (item)=>item.code===code&&item.number===Number(number)
);

      if (!residue) return [];

      return (residue.atoms||[])
        .filter((atom)=>atom.startsWith(atomInput))
        .slice(0,12)
        .map((atom) => [`${code} ${number} ${atom}`, `Atom · ${residue.chain}`]);
}

    const match = query.match(/^([A-Z]+)\s+(\d*)$/);
    if (match) {
      const [,inputCode,inputNumber] = match;
      const amino  =  aminoPick.find(
        ([code,name])=>code===inputCode||name.toUpperCase()===inputCode
);

      if (!amino) return [];

      return proteinResidues
        .filter(
          (residue)=>
            residue.code===amino[0]&&String(residue.number).startsWith(inputNumber)
)
        .slice(0,8)
        .map((residue) => [`${residue.code} ${residue.number}`, `Chain ${residue.chain}`]);
}

    return aminoPick.filter(
      ([code,name])=>code.startsWith(query)||name.toUpperCase().startsWith(query)
).slice(0,6);
})();

  const findStructureTarget = ()=>{
    const input = jumpQuery.trim().toUpperCase();
    const match = input.match(/^(?:([A-Z]{3})\s+)?(\d+)(?:\s+([A-Z][A-Z0-9']*))?$/);

    if (!match) {
      setJumpFeedback({
        success: false,
        message: "Try 158,ARG 158,or ARG 158 CD.",
});
      return;
}

    const residueNumber = Number(match[2]);

    if (!proteinResidues.some((residue)=>residue.number===residueNumber)) {
      setJumpFeedback({
        success: false,
        message: "That residue number is not present in the loaded protein structure.",
});
      return;
}

    setJumpFeedback(null);
    setGoToRes({
      residueCode: match[1]||null,
      residueNumber,
      atomName: match[3]||null,
      requestId: Date.now(),
});
};

  const handleStructureSearchResult = (result)=>{
    setJumpFeedback(result);

    if (result.success) {
      setSideTab(result.atomSearch?"atom" : "residue");
}
};

  const replayIntro = ()=>{
    setIntroStage(0);
    setShowIntro(true);
};

  const selectedAtomName = selRes?.atom?.trim().toUpperCase()||"";
  const spineAtoms = ["N","CA","C","O","OXT"];

  const atomElement =
    selectedAtomName==="N"||selectedAtomName.startsWith("N")
?"Nitrogen"
      : selectedAtomName==="O"||selectedAtomName.startsWith("O")
?"Oxygen"
        : selectedAtomName==="S"||selectedAtomName.startsWith("S")
?"Sulfur"
          : selectedAtomName.startsWith("C")
?"Carbon"
            : selectedAtomName.startsWith("H")
?"Hydrogen"
              : null;

  const atomLocation = !selectedAtomName
?null
    : spineAtoms.includes(selectedAtomName)
?"Protein backbone"
      : "Amino acid side chain";

  const atomDossier =
    atomLore[selRes?.residueCode]?.[selectedAtomName]||null;

  const sideChainPositions = {
    B: "beta",
    G: "gamma",
    D: "delta",
    E: "epsilon",
    Z: "zeta",
    H: "eta",
};
  const atomPosition = sideChainPositions[selectedAtomName.charAt(1)]||null;

  const atomResidueName =
    aminoZoo[selRes?.residueCode]?.name||
    selRes?.residueCode||
    "this residue";

  const atomElementNote =
    {
      Carbon:
        "Side-chain carbons make up the carbon skeleton of the residue. Carbons bonded only to other carbons and hydrogens are usually nonpolar,while carbons next to nitrogen or oxygen are partly polarised.",
      Nitrogen:
        "Side-chain nitrogens sit in groups such as amines,amides and rings. Depending on the group,they can carry a charge or take part in hydrogen bonds.",
      Oxygen:
        "Side-chain oxygens sit in groups such as hydroxyls,carboxylates and amides,and can often accept hydrogen bonds.",
      Sulfur:
        "Sulfur only appears in the side chains of cysteine and methionine. In cysteine it can sometimes form a disulfide bond with another cysteine.",
      Hydrogen:
        "A hydrogen is bonded to a heavier atom. Hydrogens on nitrogen or oxygen can act as hydrogen bond donors.",
}[atomElement]||"";

  const atomPlddt =
    typeof selRes?.plddt==="number"?selRes.plddt : null;

  const bandWord = (score)=>
    score>=90
?"so this atom sits in a confidently predicted stretch"
      : score>=70
?"so this stretch is predicted fairly confidently"
        : score>=50
?"so this atom's position is uncertain"
          : "so this atom's position is very uncertain";

  const atomRole = !selectedAtomName
?null
    : selectedAtomName==="N"
?Number(selRes?.residueNumber)===1
?"As the first residue,this nitrogen is the free amino end (N-terminus) of the chain,so it is not bonded to a previous residue."
        : selRes?.residueCode==="PRO"
?"This backbone nitrogen forms the peptide bond to the carbonyl carbon of the previous residue. Proline's nitrogen has no hydrogen because it is part of the side-chain ring,so unlike most residues it can't donate a hydrogen bond."
          : "This backbone nitrogen forms the peptide bond to the carbonyl carbon of the previous residue. Its N–H group can donate a hydrogen bond in secondary structures such as helices and sheets."
      : selectedAtomName==="CA"
        ? `The alpha carbon is the central carbon of ${atomResidueName}. It is bonded to the backbone nitrogen, the carbonyl carbon and the side chain (glycine only has a hydrogen here). Alpha carbon positions are the reference points FoldCheck uses for the RMSD comparison.`
        : selectedAtomName==="C"
?"This is the carbonyl carbon of the backbone. It is double bonded to the backbone oxygen and forms the peptide bond to the nitrogen of the next residue."
          : selectedAtomName==="O"
?"This is the carbonyl oxygen of the backbone. It can accept hydrogen bonds,which helps hold helices and sheets together."
            : selectedAtomName==="OXT"
?"OXT is the extra oxygen at the carboxyl end (C-terminus) of the chain,so it only appears on the last residue."
              : spineAtoms.includes(selectedAtomName)
?"This atom belongs to the protein backbone."
                : `${selectedAtomName} belongs to the side chain of ${atomResidueName}${
                    atomPosition
                      ? `, at the ${atomPosition} position counted outward from the alpha carbon`
                      : ""
}. ${atomElementNote} ${
                    atomPlddt!==null
                      ? `Residue ${selRes.residueNumber} has a pLDDT of ${atomPlddt.toFixed(1)}, ${bandWord(atomPlddt)}. pLDDT is scored on the backbone, so side-chain atoms are less certain than that number suggests. `
                      : ""
                  }This page doesn't measure whether the atom makes a hydrogen bond or forms a salt bridge.`;

  const pickedPlddt = typeof selRes?.plddt==="number"?selRes.plddt : null;

  const trustBucket =
    pickedPlddt===null
?null
      : pickedPlddt>=90
?"Very high"
        : pickedPlddt>=70
?"Confident"
          : pickedPlddt>=50
?"Low"
            : "Very low";

  const seqSpot  =  Number(selRes?.residueNumber);

  const regionsHere  =  Number.isInteger(seqSpot)
?(protCard?.regions||[]).filter(
        (region)=>
          region.type!=="Chain"&&
          region.start!=null&&
          region.end!=null&&
          Number(region.start)<=seqSpot&&
          Number(region.end)>=seqSpot
)
    : [];

  const confidenceDifference =
    pickedPlddt!==null&&modelScore
?pickedPlddt-modelScore.mean
      : null;

  const chosenAmino = selRes?aminoZoo[selRes.residueCode] : null;

  const spotStory =
    selRes&&chosenAmino&&protCard?.length&&Number.isInteger(seqSpot)
      ? `${chosenAmino.name} is classed as ${chosenAmino.category.toLowerCase()} (${chosenAmino.waterAffinity.toLowerCase()}${
          chosenAmino.sideChain
            ? `; side chain: ${
/^[A-Z][a-z]/.test(chosenAmino.sideChain)
?chosenAmino.sideChain.charAt(0).toLowerCase()+
                    chosenAmino.sideChain.slice(1)
                  : chosenAmino.sideChain
              }`
            : ""
}). Residue ${seqSpot} is about ${Math.ceil(
          (seqSpot/protCard.length)*100
)}% of the way along the ${protCard.length}-residue chain${
          seqSpot<=20||seqSpot>protCard.length-20
?",close to one end,where flexible stretches can be more common"
            : ""
        }.`
      : "";

  const proofTag = (region)=>{
    const code = region?.evidence?.[0]?.evidenceCode;
    return code==="ECO:0000269"
?"experimental evidence in a publication"
      : code==="ECO:0000305"
?"curator inference"
        : code==="ECO:0000250"
?"similarity to another protein"
          : "";
};

 
  const pdbHitCount = (start,end)=>
    pdbEntries.filter((structure)=>
      String(structure.chains||"")
        .split(",")
        .some((part)=>{
          const range = part.match(/=\s*(\d+)\s*-\s*(\d+)/);
          return range&&Number(range[1])<=start&&Number(range[2])>=end;
})
).length;

  const bandCount = (value)=>
    value>=90
?modelScore.veryHigh
      : value>=70
?modelScore.confident
        : value>=50
?modelScore.low
          : modelScore.veryLow;

  const bandSlice =
    pickedPlddt!==null&&modelScore
?pctBit(bandCount(pickedPlddt),modelScore.residueCount)
      : null;

  const topTierShare = modelScore
?pctBit(modelScore.veryHigh,modelScore.residueCount)
    : 0;
  const bottomTierShare = modelScore
?pctBit(modelScore.veryLow,modelScore.residueCount)
    : 0;

  const residueDisordered = regionsHere.some((region)=>
    region.name.toLowerCase().includes("disordered")
);

  const nearbyResidues = Number.isInteger(seqSpot)
?geomRows.filter(
        (item)=>typeof item.plddt==="number"&&Math.abs(item.number-seqSpot)<=5
)
    : [];

  const residueNeighbourText =
    nearbyResidues.length>1
      ? ` (their pLDDT runs from ${Math.min(...nearbyResidues.map((i) => i.plddt)).toFixed(
          1
)} to ${Math.max(...nearbyResidues.map((i)=>i.plddt)).toFixed(
          1
)},averaging ${(
          nearbyResidues.reduce((sum,item)=>sum+item.plddt,0)/nearbyResidues.length
        ).toFixed(1)})`
      : "";

  const residuePdbText = (()=>{
    if (!Number.isInteger(seqSpot)) return "";
    if (pdbEntries.length===0) return "UniProt lists no PDB entries for this protein.";

    const covering = pdbHitCount(seqSpot,seqSpot);
    const entryWord = pdbEntries.length===1?"entry" : "entries";

    return covering>0
      ? `UniProt lists ${pdbEntries.length} PDB ${entryWord}, and ${covering} of them ${
          covering===1?"reports" : "report"
        } a chain covering residue ${seqSpot}, so the Experimental Comparison page can show how an experimental structure compares here.`
      : `UniProt lists ${pdbEntries.length} PDB ${entryWord}, but none report a chain covering residue ${seqSpot}, so the Experimental Comparison page may have nothing to compare here.`;
})();

  const looseRegions = (protCard?.regions||[]).filter(
    (region)=>
      region.name?.toLowerCase().includes("disordered")&&
      region.start!=null&&
      region.end!=null
);

  const disorderedResidueCount = (()=>{
    const covered = new Set();
    looseRegions.forEach((region)=>{
      for (let i = Number(region.start); i<=Number(region.end); i++) covered.add(i);
});
    return covered.size;
})();

  const labVerdicts = (structure)=>{
    if (!pairedResidues.length) {
      return [
        "No residues could be matched between the two structures,so there is nothing to compare yet.",
];
}

    const average = (list)=>list.reduce((sum,value)=>sum + value,0)/list.length;
    const numbers = pairedResidues.map((item)=>item.predictedNumber);
    const scoreLookup = new Map(
      geomRows.map((residue) => [`${residue.chain}:${residue.number}`, residue.plddt])
);
    const matchedScores = pairedResidues
      .map((item) => scoreLookup.get(`${item.predictedChain}:${item.predictedNumber}`))
      .filter((value)=>typeof value==="number");
    const allScores = geomRows
      .map((residue)=>residue.plddt)
      .filter((value)=>typeof value==="number");

    const lines  =  [];

    let first = `Chain ${selectedExperimentalChain} of ${structure.id} lines up with residues ${Math.min(
      ...numbers
)}–${Math.max(...numbers)} of ${protCard?.name||activeAcc}: ${
      rmsdPairs.length
} Cα pairs${
      protCard?.length
        ? `, ${((pairedResidues.length / protCard.length) * 100).toFixed(1)}% of the ${
            protCard.length
          }-residue protein`
        : ""
    }.`;

    if (superposedRmsd!==null&&Number.isFinite(superposedRmsd)) {
      first += ` After superposition the RMSD is ${superposedRmsd.toFixed(2)} Å${
        superposedRmsd<2
?",so the matched backbone overlaps closely"
          : superposedRmsd<4
?",a broadly similar fold with some parts sitting differently"
            : ",so the backbones differ a lot over this stretch"
      }.`;
} else {
      first+=  " Fewer than three pairs had usable Cα coordinates,so no RMSD could be calculated.";
}
    lines.push(first);

    if (matchedScores.length&&allScores.length) {
      const matchedMean = average(matchedScores);
      const wholeMean = average(allScores);
      const difference = matchedMean-wholeMean;

      lines.push(
        `The matched residues have a mean pLDDT of ${matchedMean.toFixed(
          1
)},against ${wholeMean.toFixed(1)} for the whole model${
          difference>=3
?",so this comparison covers a more confidently predicted part than average."
            : difference<=-3
?",so this comparison covers a less confident part than average,and a larger RMSD there is less surprising."
              : ",so the matched part is about as confident as the rest."
        }`
);
}

    if (looseRegions.length) {
      const inDisorder = numbers.filter((number)=>
        looseRegions.some(
          (region)=>number>=Number(region.start)&&number<=Number(region.end)
)
).length;

      lines.push(
        inDisorder
          ? `${inDisorder} matched ${
              inDisorder===1?"residue sits" : "residues sit"
            } inside UniProt-annotated disordered regions, which are often poorly defined in experiments too.`
          : "None of the matched residues fall inside a UniProt-annotated disordered region."
);
}

    const method = String(structure.method||"");

    if (/nmr/i.test(method)) {
      lines.push(
        `${method} entries are usually an ensemble of models and report no resolution, so one RMSD may not capture the spread between models.`
);
} else if (Number.isFinite(parseFloat(structure.resolution))) {
      lines.push(
        `Reported resolution is ${structure.resolution} (${method}). ${
          parseFloat(structure.resolution)>=3
?"Fine detail such as side chains is less certain at this resolution,so small differences shouldn't be over-read."
            : "Backbone positions are usually well defined at this resolution,though flexible loops can still be poorly resolved."
        }`
);
} else {
      lines.push(`No resolution is reported for this ${method || "entry"}.`);
}

    lines.push(
      "RMSD only compares backbone positions in this one entry. A low value doesn't validate the prediction and a high one doesn't disprove it."
);

    return lines;
};

  const regionFacts = (()=>{
    const start = Number(selRegion?.start);
    const end = Number(selRegion?.end);
    if (!selRegion||!Number.isFinite(start)||!Number.isFinite(end)) return [];

    const facts = [];

    const others = (protCard?.regions||[]).filter(
      (region)=>
        region.type!=="Chain"&&
        !(
          region.name===selRegion.name&&
          region.start===selRegion.start&&
          region.end===selRegion.end
)&&
        Number(region.start)<=end&&
        Number(region.end)>=start
);

    facts.push(
      others.length
        ? `Overlaps with ${others.length} other annotated ${
            others.length===1?"feature" : "features"
},including ${others
            .slice(0,3)
            .map((region)=>region.name)
            .join(", ")}.`
        : "No other annotated feature overlaps this one."
);

    if (protCard?.length&&!selRegion.name.toLowerCase().includes("disordered")) {
      facts[0] += ` This region makes up ${Math.max(
        1,
        Math.round(((end-start+1)/protCard.length)*100)
      )}% of the ${protCard.length}-residue sequence.`;
}

    const scores = geomRows
      .filter((item)=>typeof item.plddt==="number"&&item.number>=start&&item.number<=end)
      .map((item)=>item.plddt);

    if (scores.length) {
      facts.push(
        `AlphaFold's pLDDT across this region averages ${(
          scores.reduce((a,b)=>a + b,0)/scores.length
).toFixed(1)} (lowest ${Math.min(...scores).toFixed(1)},highest ${Math.max(
          ...scores
).toFixed(1)})${
          modelScore ? `, compared with ${modelScore.mean.toFixed(1)} for the whole protein` : ""
        }.`
);
}

    if (pdbEntries.length) {
      const covering  =  pdbHitCount(start,end);
      facts.push(
        `${covering === 0 ? "None" : covering} of the ${pdbEntries.length} PDB ${
          pdbEntries.length===1?"entry" : "entries"
        } listed in UniProt ${covering === 1 ? "reports" : "report"} a chain covering this whole region.`
);
}

    return facts;
})();

  useEffect(()=>{
    if (!showIntro) return;

    const stageTimes = [3000,9000,15000,21000,27000,33000];
    const timers = stageTimes.map((time,index)=>
      setTimeout(()=>setIntroStage(index + 1),time)
);
    timers.push(setTimeout(()=>setShowIntro(false),40000));

    return ()=>timers.forEach(clearTimeout);
},[showIntro]);

  const errorGridStats = (()=>{
    if (!errorGrid?.length) return null;

    let total = 0;
    let count = 0;
    let low = 0;
    let moderate  =  0;
    let high = 0;
    let maximum = 0;

    for (let i = 0; i<errorGrid.length; i++) {
      const row = errorGrid[i];
      if (!Array.isArray(row)) continue;

      for (let j = 0; j<row.length; j++) {
        if (i===j) continue;

        const value  =  Number(row[j]);
        if (!Number.isFinite(value)) continue;

        total+= value;
        count++;
        maximum = Math.max(maximum,value);

        if (value<5) low++;
        else if (value<15) moderate++;
        else high++;
}
}

    if (!count) return null;

    return {
      proteinLength: errorGrid.length,
      mean: total/count,
      maximum,
      lowPercent: pctBit(low,count),
      moderatePercent: pctBit(moderate,count),
      highPercent: pctBit(high,count),
      pairCount: count,
};
})();

  

  

  if (pageNow==="confidence") {
    return (
<div className = "app">
<TopStrip activePage={pageNow} setActivePage={setPageNow}/>

<main className = "workspacePage">
<section className="workspaceTop">
<div>
<h1>Local Prediction Confidence</h1>
</div>
</section>

<LensStrip activePage={pageNow} setActivePage={setPageNow}/>

<section className="paeSection">
<h2>AlphaFold Confidence (pLDDT)</h2>
<p>
              pLDDT estimates confidence in the predicted local structure of each amino
              acid,on a scale from 0 to 100.
</p>

<div className="confidenceLegend">
<div className = "plddtGradient"/>
<div className="plddtTicks">
<span>0</span>
<span>50</span>
<span>70</span>
<span>90</span>
<span>100</span>
</div>
</div>

<div
              className="viewerShell"
              style={{
                height: "480px",
                marginTop: "24px",
                marginBottom: "28px",
                borderRadius: "12px",
                overflow: "hidden",
}}
>
<ProteinViewer
                source="alphafold"
                proteinId={activeAcc}
                onProteinConfidence={setModelScore}
                onResidueSelect={setSelRes}
/>
</div>

            {modelScore&&(
<div style={{marginBottom: "28px"}}>
<h2>Confidence Distribution</h2>
<p>Percentage of residues in each AlphaFold confidence category.</p>

                {[
                  ["Very high (90–100)",modelScore.veryHigh,"#0053d6"],
                  ["Confident (70–90)",modelScore.confident,"#65cbf3"],
                  ["Low (50–70)",modelScore.low,"#ffdb13"],
                  ["Very low (0–50)",modelScore.veryLow,"#ff7d45"],
].map(([label,count,color])=>{
                  const share = pctBit(count,modelScore.residueCount);

                  return (
<div key={label} style={{marginBottom: "18px"}}>
<div
                        style={{
                          display: "flex",
                          justifyContent: "space-between",
                          marginBottom: "8px",
}}
>
<strong>{label}</strong>
<span>{share.toFixed(1)}%</span>
</div>

<div
                        style={{
                          height: "14px",
                          background: "#1c303b",
                          borderRadius: "8px",
                          overflow: "hidden",
}}
>
<div
                          style={{
                            width: `${share}%`,
                            height: "100%",
                            background: color,
                            borderRadius: "8px",
}}
/>
</div>
</div>
);
})}
</div>
)}

<h3>Whole-protein confidence</h3>
            {modelScore?(
<>
<h2>{modelScore.mean.toFixed(1)}/100</h2>
<p>Mean pLDDT across {modelScore.residueCount} residues.</p>
</>
) : (
<p>Confidence data will load with the protein structure.</p>
)}
</section>
</main>
</div>
);
}

  if (pageNow==="pae") {
    const statSlab = (stat)=>(
<div
        key = {stat.label}
        style={{
          padding: "20px",
          background: "#0a0a0a",
          border: "1px solid #222",
          borderRadius: "12px",
}}
>
<p style={{color: "#e8eef0",fontWeight: 600,fontSize: "13px",marginTop: 0}}>
          {stat.label}
</p>
<h3 style={{fontSize: "27px",margin: "8px 0"}}>{stat.value}</h3>
<p style={{color: "#cfd8dc",fontSize: "12px",marginBottom: 0}}>
          {stat.detail}
</p>
</div>
);

    const pairVsWholeBlurb = ()=>(
<div style={{marginTop: "20px"}}>
<h3>How does this pair compare with the whole protein?</h3>

<p>
          Selected pair PAE:<strong>{Number(paeSelectedPair.error).toFixed(2)} Å</strong>
</p>

<p>
          Whole-protein mean PAE:<strong>{errorGridStats.mean.toFixed(2)} Å</strong>
</p>

<p>
          Difference:{" "}
<strong>
            {Math.abs(Number(paeSelectedPair.error)-errorGridStats.mean).toFixed(2)} Å{" "}
            {Number(paeSelectedPair.error)<errorGridStats.mean
?"below"
              : Number(paeSelectedPair.error)>errorGridStats.mean
?"above"
                : "equal to"}{" "}
            the protein average
</strong>
</p>

<p style={{color: "#b7cbd5",lineHeight: 1.7}}>
          {Number(paeSelectedPair.error)<errorGridStats.mean
?"This residue pair has lower predicted relative-position uncertainty than the average pair in this protein."
            : Number(paeSelectedPair.error)>errorGridStats.mean
?"This residue pair has higher predicted relative-position uncertainty than the average pair in this protein."
              : "This residue pair has predicted relative-position uncertainty equal to the protein average."}
</p>
</div>
);

    const paeSideBySide = ()=>(
<div
        style={{
          marginTop: "20px",
          padding: "22px",
          background: "#0a0a0a",
          border: "1px solid #222",
          borderRadius: "12px",
}}
>
<h3>Visual PAE Comparison</h3>

<div
          style={{
            display: "flex",
            justifyContent: "space-between",
            marginBottom: "10px",
            color: "#b7cbd5",
}}
>
<span>0 Å</span>
<span>32 Å</span>
</div>

<div
          style = {{
            position: "relative",
            height: "24px",
            background: "linear-gradient(to right,#1464d2,#8a997d,#ffcd28)",
            borderRadius: "8px",
}}
>
<div
            style={{
              position: "absolute",
              left: `${Math.min(100, Math.max(0, (Number(paeSelectedPair.error) / 32) * 100))}%`,
              top: "-6px",
              width: "4px",
              height: "36px",
              background: "#ffffff",
              borderRadius: "3px",
              transform: "translateX(-50%)",
}}
/>
<div
            style={{
              position: "absolute",
              left: `${Math.min(100, (errorGridStats.mean / 32) * 100)}%`,
              top: "-6px",
              width: "4px",
              height: "36px",
              background: "#ff647c",
              borderRadius: "3px",
              transform: "translateX(-50%)",
}}
/>
</div>

<div
          style = {{
            display: "flex",
            gap: "24px",
            flexWrap: "wrap",
            marginTop: "20px",
}}
>
<span>
<strong style={{color: "#ffffff"}}>●</strong>Selected pair:{" "}
            {Number(paeSelectedPair.error).toFixed(2)} Å
</span>
<span>
<strong style={{color: "#ff647c"}}>●</strong>Protein average:{" "}
            {errorGridStats.mean.toFixed(2)} Å
</span>
</div>

<p style={{color: "#9fb8c7",fontSize: "13px"}}>
          The white marker represents the selected residue pair. The pink marker
          represents the whole-protein average. Lower PAE indicates lower
          predicted relative-position uncertainty.
</p>
</div>
);

    const modelWindow = ()=>(
<div
        className="viewerShell"
        style = {{
          height: "650px",
          minHeight: "650px",
          overflow: "hidden",
          borderRadius: "12px",
}}
>
<ProteinViewer
          source="alphafold"
          proteinId={activeAcc}
          onResidueSelect={setSelRes}
          paeMatrix={errorGrid}
          paeReferenceResidue={paeReferenceResidue}
/>
</div>
);

    const predictionTakeaways = ()=>(
<div>
<h3>What does this prediction suggest?</h3>

<p>
          {Number(paeSelectedPair.error)<5
?"Low predicted aligned error: AlphaFold predicts the relative positioning of these residues with comparatively high confidence."
            : Number(paeSelectedPair.error)<15
?"Moderate predicted aligned error: AlphaFold indicates some uncertainty in their relative positioning."
              : "High predicted aligned error: AlphaFold indicates substantial uncertainty in their relative positioning."}
</p>

<p>
          {Math.abs(paeSelectedPair.x-paeSelectedPair.y)<=5
?"These residues are close together in the amino-acid sequence. Low PAE between nearby residues is common and does not establish confidence in the protein's overall fold."
            : "These residues are separated in the amino-acid sequence. Their PAE provides information about the predicted uncertainty in their relative structural arrangement."}
</p>
</div>
);

    return (
<div className="app">
<TopStrip activePage = {pageNow} setActivePage={setPageNow}/>

<main className="workspacePage">
<section className="workspaceTop">
<div>
<h1>Predicted Aligned Error</h1>
</div>
</section>

<LensStrip activePage={pageNow} setActivePage={setPageNow}/>

<section
            className="paeSection"
            style={{
              padding: "28px",
              background: "#0a0a0a",
              border: "1px solid #222",
              borderRadius: "16px",
}}
>
<div
              style = {{
                display: "grid",
                gridTemplateColumns: "repeat(auto-fit,minmax(min(100%,450px),1fr))",
                gap: "32px",
                alignItems: "start",
}}
>
<div style={{minWidth: 0}}>
                {paeLoading?(
<p>Loading AlphaFold PAE data...</p>
) : errorGrid?(
<PAEHeatmap
                    key={activeAcc}
                    matrix = {errorGrid}
                    onReferenceChange={setPaeReferenceResidue}
                    onPairSelect={setPaeSelectedPair}
/>
) : paeError?(
<p>{paeError}</p>
) : null}
</div>

<div style={{minWidth: 0}}>
<h2>3D Protein Structure</h2>

<p style={{color: "#9fb8c7"}}>
                  Explore the predicted structure alongside its relative-position uncertainty.
</p>

                {modelWindow()}
</div>
</div>
</section>

          {errorGridStats&&(
<section
              className="paeSection"
              style = {{
                marginTop: "24px",
                padding: "28px",
                background: "#0a0a0a",
                border: "1px solid #222",
                borderRadius: "16px",
}}
>
<h2>PAE Confidence Statistics</h2>

<p style={{color: "#9fb8c7"}}>
                Calculated for {protCard?.name||activeAcc} ({activeAcc}).
</p>

<div
                style = {{
                  display: "grid",
                  gridTemplateColumns: "repeat(auto-fit,minmax(190px,1fr))",
                  gap: "16px",
                  marginTop: "24px",
}}
>
                {[
                  {
                    label: "Mean PAE",
                    value: `${errorGridStats.mean.toFixed(2)} Å`,
                    detail: "Average across non-self residue pairs",
},
                  {
                    label: "Low Error",
                    value: `${errorGridStats.lowPercent.toFixed(1)}%`,
                    detail: "Predicted error below 5 Å",
},
                  {
                    label: "Moderate Error",
                    value: `${errorGridStats.moderatePercent.toFixed(1)}%`,
                    detail: "Predicted error from 5 to below 15 Å",
},
                  {
                    label: "High Error",
                    value: `${errorGridStats.highPercent.toFixed(1)}%`,
                    detail: "Predicted error of 15 Å or greater",
},
                  {
                    label: "Maximum PAE",
                    value: `${errorGridStats.maximum.toFixed(2)} Å`,
                    detail: "Highest value in the analysed matrix",
},
                  {
                    label: "Residue Pairs",
                    value: errorGridStats.pairCount.toLocaleString(),
                    detail: "Ordered pairs,excluding self-comparisons",
},
].map((stat)=>statSlab(stat))}
</div>
</section>
)}

          {paeSelectedPair&&errorGrid&&(
<section
              className="paeSection"
              style={{
                marginTop: "24px",
                padding: "28px",
                background: "#0a0a0a",
                border: "1px solid #222",
                borderRadius: "16px",
}}
>
<div>
<h2>Selected Residue Pair</h2>

<h3>
                  Residue {paeSelectedPair.y} relative to residue {paeSelectedPair.x}
</h3>

<p>
                  Predicted aligned error:<strong>{Number(paeSelectedPair.error).toFixed(2)} Å</strong>
</p>

                {predictionTakeaways()}

                {errorGridStats&&Number.isFinite(Number(paeSelectedPair.error))&&pairVsWholeBlurb()}

                {paeSelectedPair&&errorGridStats&&paeSideBySide()}

<h3 style={{marginTop: "24px"}}>Annotated biological regions</h3>

                {[paeSelectedPair.y,paeSelectedPair.x].map((residue,index)=>{
                  const matchingRegions  =  (protCard?.regions||[]).filter(
                    (region)=>
                      Number(region.start)<=residue&&
                      Number(region.end)>=residue&&
                      region.type!=="Chain"
);

                  return (
<div key={index} style={{marginTop: "12px"}}>
<strong>Residue {residue}</strong>

                      {matchingRegions.length>0?(
                        matchingRegions.map((region,regionIndex)=>(
<p key={regionIndex}>
                            {region.name} ({region.type},residues {region.start}–{region.end})
</p>
))
) : (
<p style = {{color: "#9fb8c7"}}>
                          No specific domain or region annotation available from UniProt for
                          this residue.
</p>
)}
</div>
);
})}
</div>
</section>
)}

<section
            style = {{
              marginTop: "24px",
              padding: "28px",
              background: "#0a0a0a",
              border: "1px solid #222",
              borderRadius: "16px",
}}
>
<h2>What does Predicted Aligned Error tell us?</h2>

<p>
              Predicted Aligned Error (PAE) estimates how uncertain AlphaFold is about the
              position of one residue relative to another. It is measured in ångströms (Å).
</p>

<p style={{marginTop: "20px"}}>
              Low PAE (blue): AlphaFold is fairly sure how these two residues sit relative to
              each other. Blue blocks often mark a part of the protein that is well
              constrained internally.
</p>

<p>
              High PAE (yellow): AlphaFold is less sure of their relative position,which is
              common between flexible regions or domains with an uncertain orientation.
</p>

<p>
              The diagonal: Each residue compared with itself. Low values here are expected
              and say nothing about the overall fold.
</p>

<p>
              Off-diagonal blocks: Blue blocks away from the diagonal suggest two distant
              parts of the sequence are confidently placed relative to each other. Yellow
              between blocks suggests uncertain domain positioning.
</p>

<h3 style={{marginTop: "30px"}}>PAE versus pLDDT: what's the difference?</h3>

<p>
              pLDDT estimates confidence in the local structure around individual residues.
              PAE estimates uncertainty in the relative positions of pairs of residues.
</p>

<p>
              A protein can have highly confident individual domains but uncertain
              orientations between those domains. That's why examining both pLDDT and PAE
              provides a more complete understanding of an AlphaFold prediction.
</p>
</section>
</main>
</div>
);
}

  if (pageNow==="geometry") {
    const backboneAngles = geomRows.map((residue,index)=>{
      const previous  =  geomRows[index-1];
      const next = geomRows[index + 1];
      const current = residue.backbone;

      const phi =
        previous?.chain===residue.chain&&
        previous.number===residue.number-1&&
        previous.backbone?.C&&
        current?.N&&
        current?.CA&&
        current?.C
?phiPsiTwist(previous.backbone.C,current.N,current.CA,current.C)
          : null;

      const psi =
        next?.chain===residue.chain&&
        next.number===residue.number+1&&
        current?.N&&
        current?.CA&&
        current?.C&&
        next.backbone?.N
?phiPsiTwist(current.N,current.CA,current.C,next.backbone.N)
          : null;

      return {number: residue.number,chain: residue.chain,phi,psi};
});

    const helixLike = backboneAngles.filter(
      (r)=>
        Number.isFinite(r.phi)&&
        Number.isFinite(r.psi)&&
        r.phi>-160&&
        r.phi<-20&&
        r.psi>-80&&
        r.psi<0
).length;

    const strandLike = backboneAngles.filter(
      (r)=>
        Number.isFinite(r.phi)&&
        Number.isFinite(r.psi)&&
        r.phi<-45&&
        (r.psi>90||r.psi<-170)
).length;

    const completeBackbone  =  geomRows.filter(
      (r)=>r.backbone?.N&&r.backbone?.CA&&r.backbone?.C
).length;

    const phiPsiRows = ()=>(
<tbody>
        {backboneAngles.map((residue,index)=>(
          <tr key={`${residue.chain}-${residue.number}-${index}`}>
<td>
              {residue.chain} · {residue.number}
</td>
<td>
              {residue.phi === null ? "N/A" : `${residue.phi.toFixed(1)}°`}
</td>
<td>
              {residue.psi === null ? "N/A" : `${residue.psi.toFixed(1)}°`}
</td>
</tr>
))}
</tbody>
);

    const angleLedger = ()=>(
<table>
<thead>
<tr>
<th>Residue</th>
<th>Phi (φ)</th>
<th>Psi (ψ)</th>
</tr>
</thead>
        {phiPsiRows()}
</table>
);

    return (
<div className="app">
<TopStrip activePage={pageNow} setActivePage={setPageNow}/>

<main className="workspacePage">
<div
            style={{
              position: "absolute",
              left: "-10000px",
              top: 0,
              width: "600px",
              height: "500px",
              pointerEvents: "none",
}}
>
<ProteinViewer
              source="alphafold"
              proteinId = {activeAcc}
              onGeometryLoaded={setGeomRows}
              onProteinConfidence={setModelScore}
/>
</div>

<LensStrip activePage={pageNow} setActivePage={setPageNow}/>

<section
            className="workspaceTop"
            style={{
              display: "flex",
              flexDirection: "column",
              alignItems: "stretch",
              gap: "24px",
              width: "100%",
              maxWidth: "1100px",
              margin: "0 auto",
}}
>
<h1>Protein Geometry Analysis</h1>

<p style={{fontSize: "1.15rem"}}>
              {protCard?.name||activeAcc} has{" "}
              {protCard?.length??"an unknown number of"} residues
              {modelScore?.mean!=null
                ? ` and a mean pLDDT of ${modelScore.mean.toFixed(
                    1
                  )} / 100 (AlphaFold's confidence, not a measure of geometric quality)`
                : ""}
              .
</p>

<p style = {{fontSize: "1.15rem"}}>
              {geomRows.length
                ? `${completeBackbone} of the ${geomRows.length} residues here have complete backbone coordinates. Roughly ${helixLike} fall in a helix-like φ/ψ (phi/psi) range and ${strandLike} in a strand-like range (rough boxes only, not a proper classification).`
                : "Waiting for structural coordinates..."}
</p>

<div className="geometryMetric">
<h3>Residue Backbone Angles</h3>
<div style = {{maxHeight: "350px",overflowY: "auto"}}>
                {angleLedger()}
</div>
</div>

<div className="geometryMetric">
<h3>Ramachandran Plot</h3>
<p>Backbone torsion angles for residues with both φ and ψ available.</p>

<svg
                viewBox = "0 0 400 400"
                style={{
                  width: "100%",
                  maxWidth: "440px",
                  background: "#0b1624",
                  borderRadius: "12px",
}}
>
<rect x="40" y="20" width="340" height="340" fill="#101f30" stroke="#64748b"/>
<line x1="210" y1="20" x2="210" y2="360" stroke="#475569"/>
<line x1 = "40" y1="190" x2="380" y2="190" stroke="#475569"/>
<text x="190" y="390" fill="white">
                  Phi (φ)
</text>
<text x="8" y="195" fill="white" fontSize="12">
                  ψ
</text>

                {backboneAngles
                  .filter((r)=>Number.isFinite(r.phi)&&Number.isFinite(r.psi))
                  .map((r,index)=>(
<circle
                      key={`${r.chain}-${r.number}-${index}`}
                      cx={40+((r.phi+180)/360)*340}
                      cy={360-((r.psi+180)/360)*340}
                      r="3"
                      fill="#38bdf8"
                      opacity = "0.75"
/>
))}
</svg>

<small>One point represents one residue.</small>
</div>
</section>
</main>
</div>
);
}

  if (pageNow==="experimental") {
    const pdbEntryCard = (structure)=>(
<div key={structure.id}>
<div
          style={{
            display: "grid",
            gridTemplateColumns: "repeat(auto-fit,minmax(220px,1fr))",
            gap: "16px",
            marginBottom: "28px",
}}
>
          {[
            {
              label: "Selected Protein",
              value: protCard?.gene||activeAcc,
              detail: `${protCard?.length ?? "Unknown"} amino acids in UniProt`,
},
            {
              label: "Predicted Structure",
              value: "AlphaFold",
              detail:
                modelScore?.mean!=null
                  ? `Mean pLDDT: ${modelScore.mean.toFixed(1)} / 100`
                  : "Local confidence score unavailable",
},
            {
              label: "Experimental Structure",
              value: structure.id,
              detail: structure.method,
},
            {
              label: "Experimental Resolution",
              value: structure.resolution,
              detail: "Reported by UniProt,where available",
},
].map((item)=>(
<div
              key = {item.label}
              style={{
                padding: "22px",
                background: "#0a0a0a",
                border: "1px solid #222",
                borderRadius: "12px",
}}
>
<p style={{color: "#e8eef0",fontWeight: 600,marginTop: 0}}>
                {item.label}
</p>
<h3 style = {{fontSize: "24px",margin: "10px 0"}}>{item.value}</h3>
<p style={{color: "#cfd8dc",fontSize: "13px"}}>{item.detail}</p>
</div>
))}
</div>

<p className="foldcheckLinkedHint">
          Select a residue on either 3D model to highlight the corresponding residue
          in both structures,when a matching residue is present. The experimental
          structure may cover only part of the protein.
</p>

<div
          style={{
            display: "grid",
            gridTemplateColumns: "repeat(auto-fit,minmax(min(100%,420px),1fr))",
            gap: "24px",
            alignItems: "start",
            marginTop: "24px",
}}
>
<section style={{minWidth: 0}}>
<h3>AlphaFold Prediction</h3>
<div
              className="viewerShell"
              style={{
                height: "500px",
                overflow: "hidden",
                borderRadius: "12px",
                marginTop: "16px",
}}
>
<ProteinViewer
                source="alphafold"
                proteinId={activeAcc}
                onGeometryLoaded={setGeomRows}
                linkedResidue={linkedSelection?.alphafold}
                onResidueSelect={(residue)=>
                  mirrorPick("alphafold",residue)
}
/>
</div>
<p>
<strong>Protein:</strong>{protCard?.name||activeAcc}
</p>
<p>
<strong>Sequence length:</strong>{" "}
              {protCard?.length??"Unknown"} amino acids
</p>
<p>
<strong>Mean pLDDT:</strong>{" "}
              {modelScore?.mean!=null
?modelScore.mean.toFixed(1)
                : "Unavailable"}
</p>
</section>

<section style={{minWidth: 0}}>
<h3>Experimental PDB {structure.id}</h3>
<div
              className="viewerShell"
              style={{
                height: "500px",
                overflow: "hidden",
                borderRadius: "12px",
                marginTop: "16px",
}}
>
<ProteinViewer
                source="pdb"
                proteinId = {structure.id}
                onGeometryLoaded={setLabRows}
                linkedResidue={linkedSelection?.pdb}
                onResidueSelect={(residue)=>mirrorPick("pdb",residue)}
/>
</div>
<p>
<strong>Method:</strong>{structure.method}
</p>
<p>
<strong>Resolution:</strong>{structure.resolution}
</p>
<p>
<strong>Chains/coverage:</strong>{structure.chains}
</p>
</section>
</div>

<a
          href={`https://www.rcsb.org/structure/${encodeURIComponent(structure.id)}`}
          target = "_blank"
          rel="noopener noreferrer"
>
          View PDB {structure.id} at RCSB ↗
</a>

<div
          style={{
            marginTop: "28px",
            padding: "28px",
            background: "#0a0a0a",
            border: "1px solid #222",
            borderRadius: "14px",
}}
>
<div
            style = {{
              padding: "22px",
              marginBottom: "24px",
              background: "#0a0a0a",
              border: "1px solid #222",
              borderRadius: "12px",
}}
>
<h3>Cα Backbone RMSD</h3>

<div
              style={{
                display: "flex",
                gap: "12px",
                flexWrap: "wrap",
                margin: "14px 0",
}}
>
<div>
<strong>Selected PDB chain:</strong>
                {selectedExperimentalChain??"Unavailable"}
</div>

<div>
<strong>Matched residues:</strong>
                {pairedResidues.length}
</div>

<div>
<strong>Experimental coverage:</strong>
                {experimentalSequence.length>0
                  ? `${pctBit(
                      pairedResidues.length,
                      experimentalSequence.length
                    ).toFixed(1)}%`
                  : "Unavailable"}
</div>
</div>

            {superposedRmsd!==null&&Number.isFinite(superposedRmsd)?(
<>
<h2 style={{fontSize: "34px",margin: "12px 0"}}>
                  {superposedRmsd.toFixed(2)} Å
</h2>
<p>
                  Calculated after optimal rigid-body superposition using{" "}
<strong>{rmsdPairs.length}</strong>matched Cα atom pairs.
</p>
</>
) : (
<p>
                RMSD unavailable. At least three matched residues with valid Cα
                coordinates are required.
</p>
)}
</div>

<h2>Interpreting the Evidence</h2>

<div
            style={{
              display: "grid",
              gridTemplateColumns: "repeat(auto-fit,minmax(240px,1fr))",
              gap: "20px",
              marginTop: "20px",
}}
>
<div>
<h3>AlphaFold Prediction</h3>
<p>
                AlphaFold predicts the protein's three-dimensional structure from
                its amino acid sequence.
</p>
<p>
                {modelScore?.mean!=null
                  ? `The loaded model has a mean pLDDT of ${modelScore.mean.toFixed(
                      1
                    )} / 100. This describes predicted local confidence, not experimentally measured accuracy.`
                  : "Local confidence information is not currently available."}
</p>
</div>

<div>
<h3>Experimental Structure</h3>
<p>
                PDB {structure.id} is linked to this protein through UniProt. Its
                reported determination method is {structure.method}.
</p>
<p>
                Reported resolution: {structure.resolution}. The experimental entry
                may represent only a portion of the protein or contain additional
                molecular components.
</p>
</div>
</div>

<h3 style={{marginTop: "24px"}}>What can we conclude?</h3>

          {labVerdicts(structure).map((line,index)=>(
<p key={index}>{line}</p>
))}
</div>
</div>
);

    const pdbPicker = ()=>(
<select
        id="experimentalPdbSelect"
        value={selectedExperimentalId}
        onChange={(event)=>setSelectedExperimentalId(event.target.value)}
        style={{
          display: "block",
          width: "100%",
          maxWidth: "500px",
          margin: "14px 0 24px",
          padding: "14px",
          background: "#0a0a0a",
          color: "#ffffff",
          border: "1px solid #222",
          borderRadius: "10px",
}}
>
        {pdbEntries.map((structure)=>(
<option key={structure.id} value={structure.id}>
            {structure.id} · {structure.method} · {structure.resolution}
</option>
))}
</select>
);

    return (
<div className="app">
<TopStrip activePage={pageNow} setActivePage={setPageNow}/>

<main className="workspacePage">
<section className = "workspaceTop">
<div>
<h1>Experimental Structure Comparison</h1>
</div>
</section>

<LensStrip activePage = {pageNow} setActivePage={setPageNow}/>

<section className = "paeSection" style={{marginTop: "24px"}}>
<h2>Available PDB Structures</h2>

<p>
              These structures are linked to {protCard?.name||activeAcc} through
              UniProt (may cover only part of the protein or include other molecules).
</p>

            {experimentalLoading?(
<p>Searching UniProt for experimental structures...</p>
) : experimentalError?(
<p role="alert">{experimentalError}</p>
) : pdbEntries.length===0?(
<p>No PDB cross-references were found for this protein in UniProt.</p>
) : (
<>
<label htmlFor="experimentalPdbSelect">Select an experimental structure</label>

                {pdbPicker()}

                {pdbEntries
                  .filter((structure)=>structure.id===selectedExperimentalId)
                  .map((structure)=>pdbEntryCard(structure))}
</>
)}
</section>
</main>
</div>
);
}

  if (pageNow==="basics") {
    const letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".split("");

    const visibleTopics = termBits.filter(([title,description])=>{
      const matchesSearch = `${title} ${description}`
        .toLowerCase()
        .includes(glossaryQuery.trim().toLowerCase());
      const matchesLetter =
        glossaryLetter==="All"||title.toUpperCase().startsWith(glossaryLetter);
      return matchesSearch&&matchesLetter;
});

    return (
<div className="app">
<TopStrip activePage={pageNow} setActivePage={setPageNow}/>

<main
          className="workspacePage"
          style={{maxWidth: 1050,margin: "0 auto",padding: "40px 24px"}}
>
<h1>Protein Basics</h1>
<p>
            Explore protein architecture,cellular localisation and the evidence behind
            structure predictions.
</p>

<div style={{margin: "30px 0 20px"}}>
<label htmlFor="basicsSearch" style={{display: "block",marginBottom: 10}}>
              Search the glossary
</label>
<input
              id="basicsSearch"
              type = "search"
              value={glossaryQuery}
              onChange={(event)=>setGlossaryQuery(event.target.value)}
              placeholder="Search a term or explanation..."
              style={{
                width: "100%",
                padding: "13px 16px",
                background: "#111d26",
                color: "#fff",
                border: "1px solid #29414d",
                borderRadius: 8,
}}
/>
</div>

<div
            aria-label="Filter glossary by letter"
            style = {{display: "flex",flexWrap: "wrap",gap: 7,marginBottom: 30}}
>
            {["All",...letters].map((letter)=>(
<button
                key={letter}
                type="button"
                onClick={()=>setGlossaryLetter(letter)}
                aria-pressed={glossaryLetter===letter}
                disabled = {
                  letter!=="All"&&
                  !termBits.some(([title])=>title.toUpperCase().startsWith(letter))
}
                style={{
                  padding: "8px 10px",
                  background: glossaryLetter===letter?"#28535a" : "transparent",
                  color: "#d6e7ec",
                  border: "1px solid #29414d",
                  borderRadius: 5,
                  cursor: "pointer",
}}
>
                {letter}
</button>
))}
</div>

<section aria-label="Protein science glossary">
            {visibleTopics.length?(
              visibleTopics.map(([title,description])=>(
<article key={title} style={{padding: "18px 0"}}>
<h2 style = {{margin: "0 0 8px",fontSize: "1.2rem"}}>{title}</h2>
<p style={{margin: 0,lineHeight: 1.75}}>{description}</p>
</article>
))
) : (
<p>No matching terms. Try another search or letter.</p>
)}
</section>
</main>
</div>
);
}

  if (pageNow==="limitations") {
    return (
<div className="app">
<TopStrip activePage={pageNow} setActivePage={setPageNow}/>

<main
          className="workspacePage"
          style={{maxWidth: 1050,margin: "0 auto",padding: "40px 24px"}}
>
<h1>Limitations</h1>
<p>What FoldCheck's predictions,confidence scores and comparisons can and cannot tell you.</p>

<section aria-label="Scientific limitations" style={{marginTop: 30}}>
            {caveatBits.map(([title,description])=>(
<article key={title} style={{padding: "18px 0"}}>
<h2 style={{margin: "0 0 8px",fontSize: "1.2rem"}}>{title}</h2>
<p style={{margin: 0,lineHeight: 1.75}}>{description}</p>
</article>
))}
</section>
</main>
</div>
);
}

  if (pageNow==="finder") {
    const pick = finderPicks?.[0];
    const speciesWord = critterChips[finderSpecies][2];
    const roleWords = finderRole!==null?jobChips[finderRole][2] : "";

    const finderSiblings = ()=>(
<div className="finderOthers">
<span>Also try</span>
        {finderPicks.slice(1).map((other)=>(
<button
            key={other.id}
            type="button"
            className="finderChip"
            onClick={()=>visitProtein(other.id)}
>
            {other.name}
</button>
))}
</div>
);

    return (
<div className="app">
<TopStrip activePage={pageNow} setActivePage={setPageNow}/>

<main className="finderPage">
<h1>Find a protein</h1>

<h2>Which living thing?</h2>
<div className="finderChips">
            {critterChips.map(([label],index)=>(
<button
                key={label}
                type="button"
                className={finderSpecies===index?"finderChip active" : "finderChip"}
                onClick={()=>setFinderSpecies(index)}
>
                {label}
</button>
))}
</div>

<h2>What should it do?</h2>
<div className="finderChips">
            {jobChips.map(([label],index)=>(
<button
                key={label}
                type="button"
                className = {finderRole===index?"finderChip active" : "finderChip"}
                onClick={()=>setFinderRole(index)}
>
                {label}
</button>
))}
</div>

<button
            type = "button"
            className="finderMoreButton"
            onClick={()=>setFinderMore(!finderMore)}
>
            {finderMore?"Fewer options ↑" : "More options ↓"}
</button>

          {finderMore&&(
<div className="finderMore">
<h2>Any word it should include?</h2>
<input
                type="text"
                className = "finderInput"
                value = {finderWord}
                onChange={(event)=>setFinderWord(event.target.value)}
                placeholder = "e.g. heart,blood,liver"
/>

<h2>How big?</h2>
<div className="finderChips">
                {[
                  ["small","Small"],
                  ["any","Any size"],
                  ["large","Large"],
].map(([value,label])=>(
<button
                    key = {value}
                    type="button"
                    className={finderSize===value?"finderChip active" : "finderChip"}
                    onClick={()=>setFinderSize(value)}
>
                    {label}
</button>
))}
</div>

<label className="finderCheck">
<input
                  type = "checkbox"
                  checked={finderStructure}
                  onChange={(event)=>setFinderStructure(event.target.checked)}
/>
                Only proteins scientists have already imaged in a lab
</label>
</div>
)}

<button
            type = "button"
            className="finderGo"
            disabled={finderRole===null||finderLoading}
            onClick={findProteins}
>
            {finderLoading?"Searching..." : "Find my protein →"}
</button>

          {finderError&&(
<p className="finderError" role="alert">
              {finderError}
</p>
)}
          {proteinSearchError&&(
<p className="finderError" role="alert">
              {proteinSearchError}
</p>
)}

          {pick&&(
<div className="finderResult">
<h2>{pick.name}</h2>
<p>
                {roleWords
                  ? `${pick.name} is a ${speciesWord ? speciesWord + " " : ""}protein that ${roleWords}.`
                  : `${pick.name} is a well-studied protein${
                      pick.species?" from "+pick.species : ""
                    }.`}
</p>
<button
                type="button"
                className="finderGo"
                disabled = {proteinSearchLoading}
                onClick={()=>visitProtein(pick.id)}
>
                {proteinSearchLoading ? "Opening..." : `Explore ${pick.name} →`}
</button>

              {finderPicks.length>1&&finderSiblings()}
</div>
)}
</main>
</div>
);
}

  if (pageNow==="sources") {
    return (
<div className="app">
<TopStrip activePage={pageNow} setActivePage={setPageNow}/>

<main className="sourcesPage">
<h1>Sources</h1>
<p>Where FoldCheck gets its protein information,and the tools it uses.</p>

<div className="sourcesLinks">
            <a href="https://www.uniprot.org/" target="_blank" rel="noopener noreferrer">
              UniProt ↗
</a>
            <a href="https://alphafold.ebi.ac.uk/" target="_blank" rel="noopener noreferrer">
              AlphaFold Protein Structure Database ↗
</a>
            <a href="https://www.rcsb.org/" target="_blank" rel="noopener noreferrer">
              RCSB Protein Data Bank ↗
</a>
            <a href="https://europepmc.org/" target="_blank" rel="noopener noreferrer">
              Europe PMC ↗
</a>
            <a href="https://molstar.org/" target="_blank" rel="noopener noreferrer">
              Mol*3D Viewer ↗
</a>
</div>
</main>
</div>
);
}

  if (pageNow==="about") {
    return (
<div className="app">
<TopStrip activePage = {pageNow} setActivePage={setPageNow}/>

<main className="aboutPage">
<button className="aboutBack" onClick={()=>setPageNow("home")}>
            ← Back to home
</button>
<h1>About Me</h1>

<div className="aboutStory">
<p>
              Hey! My name is Yousef,and I'm a young student from Australia. I love science,
              especially biology and chemistry,and recently I've also started getting into
              coding and circuitry. I created FoldCheck to get an ESP32,a pair of Nothing
              headphones,and maybe even a 3D printer from Hack Club. Okay,I'm kidding. Those
              are just some of the perks.
</p>
<p>
              It all started on a "moony" evening (like a sunny afternoon,but with the moon). I
              was studying proteins when I suddenly wondered: How do we actually know what
              proteins look like?So,I did some research and came across AlphaFold,an AI that
              can predict the 3D structures of proteins. I thought that was crazy! AI could
              predict how a protein folds?But then I had another question. How does it know if
              it's correct?
</p>
<p>
              I found out that AlphaFold gives confidence scores for its predictions,but I
              wanted to understand what those scores meant and how much we could actually trust
              them. So,I decided to (attempt to) help the world out. That's how FoldCheck
              started. I wanted to make a website where people could look at proteins,explore
              their structures,and understand how confident AlphaFold is in its predictions.
</p>
<p>
              I'm still learning,and there's definitely a lot more I want to add to FoldCheck.
              One of my biggest goals is to eventually bring FoldCheck into schools,so students
              can explore proteins themselves instead of just looking at pictures in a textbook.
              And yeah,hopefully I'll get that 3D printer too.
</p>
</div>

<button className="primaryButton" onClick={()=>setPageNow("home")}>
            Explore FoldCheck →
</button>
</main>
</div>
);
}

  if (pageNow==="workspace") {
    const functionTeaser = ()=>(
<summary>
<span className="functionPreviewText">
          {protCard.function.length>280
            ? `${protCard.function.slice(0, 280).trimEnd()}...`
            : protCard.function}
</span>

        {protCard.function.length>280&&(
<span className="functionExpandHint">
<span className="showMoreText">More details ↓</span>
<span className = "showLessText">Show less ↑</span>
</span>
)}
</summary>
);

    const regionRowBit = (index,region)=>(
<div className="functionalRegionRow" key={index}>
<span className = "functionalRegionName">{region.name}</span>

<button
          type="button"
          className="functionalRegionRange"
          disabled = {region.start==null||region.end==null}
          onClick={()=>{
            setSelRegion(region);
            setSelRes(null);
            setJumpQuery("");
            setJumpFeedback(null);
            setSideTab("region");

            document.querySelector(".structureWorkspace")?.scrollIntoView({
              behavior: "smooth",
              block: "start",
});
}}
          title="Explore this region in 3D"
>
          Residues {region.start??"?"}–{region.end??"?"}
<span aria-hidden="true">↗</span>
</button>
</div>
);

    const searchFieldWrap = ()=>(
<div className="structureSearchInputWrap">
<input
          type="text"
          value={jumpQuery}
          onChange={(event)=>{
            setJumpQuery(event.target.value);
            setShowStructureSuggestions(true);
}}
          onFocus={()=>setShowStructureSuggestions(true)}
          onKeyDown={(event)=>{
            if (event.key==="Escape") {
              setShowStructureSuggestions(false);
              setActiveSuggestionIndex(-1);
}

            if (showStructureSuggestions&&structureSuggestions.length>0) {
              if (event.key==="ArrowDown") {
                event.preventDefault();
                setActiveSuggestionIndex((index)=>
                  Math.min(index+1,structureSuggestions.length-1)
);
}

              if (event.key==="ArrowUp") {
                event.preventDefault();
                setActiveSuggestionIndex((index)=>Math.max(index-1,0));
}

              if (event.key==="Enter"&&activeSuggestionIndex>=0) {
                event.preventDefault();
                const [code] = structureSuggestions[activeSuggestionIndex];
                setJumpQuery(code.includes(" ") ? code : `${code} `);
                setShowStructureSuggestions(false);
                setActiveSuggestionIndex(-1);
}
}
}}
          placeholder="e.g. ARG 158 CD"
          aria-label = "Search residue or atom"
          autoComplete = "off"
/>

        {showStructureSuggestions&&structureSuggestions.length>0&&(
<div className="structureSuggestions">
            {structureSuggestions.map(([code,name],index)=>(
<button
                key={code}
                type = "button"
                className={`structureSuggestion ${
                  activeSuggestionIndex===index?"isActive" : ""
                }`}
                onClick={()=>{
                  setJumpQuery(code.includes(" ") ? code : `${code} `);
                  setShowStructureSuggestions(false);
}}
>
<strong>{code}</strong>
<span>{name}</span>
</button>
))}
</div>
)}
</div>
);

    const resetViewBtn = ()=>(
<button
        type="button"
        className="structureResetButton"
        onClick={()=>{
          setResetCam((value)=>value+1);
          setSelRegion(null);
          setSelRes(null);
          setJumpQuery("");
          setShowStructureSuggestions(false);
          setSideTab("residue");
          setJumpFeedback(null);
}}
>
        Reset View
</button>
);

    const confidenceBuckets = ()=>(
<>
<p>
          Calculated from {modelScore.residueCount} residues in the AlphaFold
          model of {protCard?.name||activeAcc} ({activeAcc}).
</p>

<div className="proteinConfidenceDistribution">
          {[
            ["Very high (≥90)",modelScore.veryHigh],
            ["Confident (70–90)",modelScore.confident],
            ["Low (50–70)",modelScore.low],
            ["Very low (&lt;50)",modelScore.veryLow],
].map(([label,count])=>(
<div key = {label}>
<span>{label}</span>
<strong>
                {pctBit(count,modelScore.residueCount).toFixed(1)}%
</strong>
</div>
))}
</div>

<div className="proteinConfidenceInterpretation">
<h4>What does this tell us about this protein?</h4>
<p>
            This model scores {modelScore.residueCount} residues with a mean
            pLDDT of {modelScore.mean.toFixed(1)}. {modelScore.veryHigh}{" "}
            of them ({topTierShare.toFixed(1)}%) score 90 or above and{" "}
            {modelScore.veryLow} ({bottomTierShare.toFixed(1)}%) score below 50,
            {topTierShare>=15&&bottomTierShare>=15
?" so the mean hides a big gap between its best-and worst-predicted parts."
              : " so the mean is only a rough summary of how the scores are spread."}
</p>

<h4>Why might its confidence vary?</h4>
<p>
            {looseRegions.length
              ? `UniProt annotates ${looseRegions.length} "Disordered" ${
                  looseRegions.length===1?"region" : "regions"
}${
                  protCard?.length
                    ? ` covering ${disorderedResidueCount} residues (${Math.round(
                        (disorderedResidueCount/protCard.length)*100
                      )}% of the sequence)`
                    : ""
}. Low pLDDT often goes with disorder${
                  bottomTierShare>=10?",and that fits here" : ""
                }, but the scores alone can't prove why a stretch is uncertain.`
              : "UniProt has no disordered-region annotation for this protein,so the data loaded here can't say why some stretches score low."}
</p>

<h4>What should you investigate?</h4>
<p>
            Pick a stretch scoring above 90 and one below 50 and compare their shapes
            in the viewer. Then check the PAE page to see whether the confident parts
            are also confidently placed relative to each other.{" "}
            {pdbEntries.length
              ? `UniProt lists ${pdbEntries.length} PDB ${
                  pdbEntries.length===1?"structure" : "structures"
                } for this protein, so you can test the prediction on the Experimental Comparison page.`
              : "UniProt lists no PDB structures for this protein,so there is no experimental structure to compare against here."}
</p>
</div>
</>
);

    const residueVsAverage = ()=>(
<div className="residueConfidenceComparison">
<strong>
          {confidenceDifference>=0?"+" : ""}
          {confidenceDifference.toFixed(1)} pLDDT points
</strong>

<p>
          {chosenAmino?.name||selRes.residueCode}{" "}
          {seqSpot} has{" "}
          {confidenceDifference>=0?"higher" : "lower"} local
          prediction confidence than the average residue (
          {modelScore.mean.toFixed(1)}) in the AlphaFold model of{" "}
          {protCard?.name||activeAcc} ({activeAcc}).
</p>
</div>
);

    const residueSentence = ()=>(
<p>
        Residue {seqSpot} is{" "}
        {/^[aeiou]/i.test(chosenAmino?.name||selRes.residueCode)
?"an"
          : "a"}{" "}
        {(chosenAmino?.name||selRes.residueCode).toLowerCase()}{" "}
        with a pLDDT of {pickedPlddt.toFixed(1)},in the{" "}
        {trustBucket} band (
        {pickedPlddt>=90
?"90–100"
          : pickedPlddt>=70
?"70–90"
            : pickedPlddt>=50
?"50–70"
              : "0–50"}
).
        {confidenceDifference!==null&&
          ` That is ${Math.abs(confidenceDifference).toFixed(1)} points ${
            confidenceDifference>=0?"above" : "below"
} this protein’s mean of ${modelScore.mean.toFixed(1)}${
            bandSlice!==null
              ? `, and ${bandSlice.toFixed(1)}% of this model’s residues fall in the same band`
              : ""
          }.`}
        {pickedPlddt>=90
?" It means AlphaFold is confident about the local backbone here,not about the side chain or any interaction."
          : pickedPlddt>=70
?" AlphaFold is fairly confident about the local backbone,but details like side-chain positions are less certain."
            : pickedPlddt>=50
?" The local shape is uncertain,so don't read exact distances or contacts here."
              : " AlphaFold has very low confidence here. The stretch may be flexible or disordered,but the score alone can't say which."}
</p>
);

    const contextSentence = ()=>(
<p>
        {spotStory}
        {regionsHere.length
          ? ` UniProt places position ${seqSpot} inside ${
              regionsHere.length
} annotated ${
              regionsHere.length===1?"feature" : "features"
},including ${regionsHere
              .slice(0,2)
              .map((region) => `${region.name} (residues ${region.start}–${region.end})`)
              .join(" and ")}. Those are sequence annotations${
              proofTag(regionsHere[0])
                ? ` (the first is backed by ${proofTag(regionsHere[0])})`
                : ""
}. ${
              regionsHere
                .slice(0,2)
                .some((region)=>region.name.toLowerCase().includes("interaction"))
?"They link this stretch to those partners (pLDDT can't show whether a partner actually binds here)."
                : "They describe the stretch biologically (pLDDT can't confirm any function or binding here)."
            }`
          : ` UniProt has no annotated feature at position ${seqSpot} (that doesn't mean the residue has no function).`}
        {residueDisordered&&
          pickedPlddt<50&&
          " UniProt annotates this stretch as disordered,which fits the low score."}
        {residueDisordered&&
          pickedPlddt>=90&&
          " UniProt annotates this stretch as disordered,yet AlphaFold models it confidently. That mismatch is worth treating cautiously."}
</p>
);

    const neighbourHint = ()=>(
<p>
        Look at residues {Math.max(1,seqSpot-5)}–
        {Math.min(
          protCard?.length||seqSpot+5,
          seqSpot+5
)}{" "}
        in the 3D viewer and see whether the colour changes
        {residueNeighbourText}. Then open the PAE page and click the
        heatmap at residue {seqSpot} on the horizontal axis to see
        which other regions AlphaFold places confidently relative to it.{" "}
        {residuePdbText}
</p>
);

    const residuePanel = ()=>(
<div className="residueInspector">
<h3>{chosenAmino?.name||selRes.residueCode}</h3>

<ul className="residueFacts">
<li>Code: {selRes.residueCode}</li>
<li>Residue {selRes.residueNumber}</li>
<li>Chain {selRes.chain}</li>
</ul>

<div className="inspectorSection">
<h4>AlphaFold prediction confidence</h4>

          {pickedPlddt!==null?(
<>
<div className = "localConfidenceScale">
<div className="localConfidenceReadout">
<span>SELECTED RESIDUE · pLDDT</span>
<strong>{pickedPlddt.toFixed(1)}/100</strong>
</div>

<div className="localConfidenceTrack">
<div
                    className="localConfidenceMarker"
                    style={{ left: `${pickedPlddt}%` }}
/>
</div>

<div className="localConfidenceTicks">
<span>0</span>
<span>50</span>
<span>70</span>
<span>90</span>
<span>100</span>
</div>

<p className="localConfidenceCategory">
                  Classification:<strong>{trustBucket}</strong>
</p>
</div>

<p className="inspectorNote">
                pLDDT estimates AlphaFold's confidence in the local structure around{" "}
                {chosenAmino?.name||selRes.residueCode}{" "}
                {seqSpot}.
</p>

<div className="confidenceInterpretation">
                {confidenceDifference!==null&&residueVsAverage()}

<h4>What does this score tell us about this residue?</h4>
                {residueSentence()}

<h4>What might explain its confidence?</h4>
                {contextSentence()}

<h4>What should you check next?</h4>
                {neighbourHint()}
</div>
</>
) : (
<p>Confidence information is unavailable for this selection.</p>
)}
</div>

</div>
);

    const disorderNote = ()=>(
<p>
        UniProt annotates residues {selRegion.start}–
        {selRegion.end} as disordered (
        {selRegion.end-selRegion.start+1} residues
        {protCard?.length
          ? `, ${Math.round(
              ((selRegion.end-selRegion.start + 1)/
                protCard.length)*
                100
            )}% of the sequence`
          : ""}
). Disordered stretches don't settle into one fixed shape,and
        AlphaFold usually scores them low and draws them as loose coils,so
        don't read these coordinates as a true structure. Disordered doesn't
        mean useless: short binding motifs often sit in these stretches.
</p>
);

    const regionBlurbPick = ()=>(
<p>
        {protCard?.function
?.split(/(?<=[.!?])\s+/)
          .filter((sentence)=>{
            const terms  =  selRegion.name
              .toLowerCase()
              .replace(/\b(region|interaction|with|binding)\b/g,"")
              .trim()
              .split(/\s+/)
              .filter((term)=>term.length>=4);

            return terms.some((term)=>
              sentence.toLowerCase().includes(term)
);
})
          .slice(0,2)
          .join(" ")||
          `No specific biological explanation is available for ${selRegion.name} in the retrieved protein function annotation.`}
</p>
);

    const regionPanel = ()=>(
<div className="residueInspector">
<h3>{selRegion.name}</h3>

<div className = "inspectorSection">
          {!selRegion.name.toLowerCase().includes("interaction")&&(
<h4>
              {selRegion.type==="Domain"
                ? `What does the ${selRegion.name} domain do?`
                : `What does ${selRegion.name} do?`}
</h4>
)}

          {selRegion.type==="Domain"?(
<p>
              {protCard?.function
?.split(/(?<=[.!?])\s+/)
                .filter((sentence)=>{
                  const terms = selRegion.name
                    .toLowerCase()
                    .replace(/\b(domain|region|motif)\b/g,"")
                    .trim()
                    .split(/\s+/)
                    .filter((term)=>term.length>=4);

                  return terms.some((term)=>
                    sentence.toLowerCase().includes(term)
);
})
                .slice(0,2)
                .join(" ")||
                `No domain-specific functional explanation was found in the retrieved UniProt annotation for ${selRegion.name}.`}
</p>
) : selRegion.type==="Region"?(
<>
              {selRegion.name.toLowerCase().includes("interaction")?(
<></>
) : selRegion.name.toLowerCase().includes("disordered")?disorderNote() : regionBlurbPick()}
</>
) : selRegion.type==="Signal"?(
<p>
              This signal peptide is a sequence that helps direct the protein into the
              cellular pathway responsible for secretion or delivery to certain
              compartments. It spans residues {selRegion.start}–
              {selRegion.end} ({selRegion.end-selRegion.start+1}{" "}
              amino acids) and is usually cleaved off once the protein has been
              targeted,so it may be missing from the mature protein.
</p>
) : selRegion.type==="Chain"?(
<p>
              This annotation identifies a protein chain,which may represent the
              mature protein after processing.
              {protCard?.length&&
                ` Here it spans residues ${selRegion.start}–${selRegion.end} of the ${protCard.length}-residue sequence${
                  Number(selRegion.start)<=1&&
                  Number(selRegion.end)>=protCard.length
?",so it covers the whole sequence."
                    : ",so part of the original sequence is not included in this chain."
                }`}
</p>
) : (
<p>
              UniProt classifies this feature as{" "}
<strong>{selRegion.type||"a functional region"}</strong>. Its
              annotation describes a sequence feature of the protein that spans
              residues {selRegion.start}–{selRegion.end} (
              {selRegion.end-selRegion.start+1} amino acids).
</p>
)}

          {regionFacts.map((fact,index)=>(
<p key={index}>{fact}</p>
))}

<RegionDetails
            region={selRegion}
            proteinInfo = {protCard}
            partnerInfo={partnerInfo}
/>

          {selRegion.evidence
?.filter((evidence)=>evidence.source==="PubMed"&&evidence.id)
            .map((evidence,index,list)=>(
<a
                key={index}
                href={`https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(
                  evidence.id
                )}/`}
                target="_blank"
                rel = "noopener noreferrer"
>
                {list.length>1
                  ? `Supporting publication ${index + 1}`
                  : "View supporting publication"}{" "}
                ↗
</a>
))}

          {selRegion.description&&
            selRegion.description!==selRegion.name&&(
<div className="inspectorNote">
<strong>Region-specific function</strong>
<p>{selRegion.description}</p>
</div>
)}
</div>

<div className="inspectorSection">
<h4>Region information</h4>

<div className="inspectorData">
<div>
<span>Start residue</span>
<strong>{selRegion.start}</strong>
</div>

<div>
<span>End residue</span>
<strong>{selRegion.end}</strong>
</div>

<div>
<span>Region length</span>
<strong>{selRegion.end-selRegion.start+1} aa</strong>
</div>
</div>
</div>
</div>
);

    const atomPanel = ()=>(
<div className="residueInspector">
<h3>{selectedAtomName||"No atom selected"}</h3>

<ul className="residueFacts">
<li>{chosenAmino?.name||selRes.residueCode}</li>
<li>Residue {selRes.residueNumber}</li>
<li>Chain {selRes.chain}</li>
</ul>

<div className="inspectorSection">
<h4>Atom Explorer</h4>

          {selectedAtomName?(
<>
<div className="inspectorData">
<div>
<span>Atom identifier</span>
<strong>{selectedAtomName}</strong>
</div>
<div>
<span>Chemical element</span>
<strong>{atomElement||"Unknown"}</strong>
</div>
<div>
<span>Structural location</span>
<strong>{atomLocation}</strong>
</div>
</div>

<h4>What does this atom do?</h4>
<p>{atomDossier?.role||atomRole}</p>

              {atomDossier&&(
<>
<h4>Normally bonded to</h4>
<p>{atomDossier.bonds.join(" · ")}</p>

<h4>Why is this atom important?</h4>
<p>{atomDossier.significance}</p>
</>
)}
</>
) : (
<p>Select an atom in the 3D viewer.</p>
)}
</div>
</div>
);

    const functionLongList = ()=>(
<details className="biologicalFunctionDetails">
        {functionTeaser()}

        {protCard.function.length>280&&(
<ul className="biologicalFunctionFull">
            {protCard.function
              .split(/(?<=[.!?;])\s+|;\s*/)
              .filter((sentence)=>sentence.trim().length>0)
              .map((sentence,index)=>(
<li key={index}>{sentence.trim()}</li>
))}
</ul>
)}
</details>
);

    const findForm = ()=>(
<form
        className="structureSearchForm"
        onSubmit={(event)=>{
          event.preventDefault();
          setSelRegion(null);
          findStructureTarget();
}}
>
        {searchFieldWrap()}

<button type="submit">Find →</button>
        {resetViewBtn()}
</form>
);

    const foldedModelBox = ()=>(
<ProteinViewer
        source="alphafold"
        proteinId={activeAcc}
        onResidueSelect={setSelRes}
        onProteinConfidence={setModelScore}
        searchRequest = {goToRes}
        onSearchResult={handleStructureSearchResult}
        resetViewRequest={resetCam}
        onSequenceLoaded={setProteinResidues}
        onGeometryLoaded={setGeomRows}
        selectedRegion={selRegion}
/>
);

    const legendHeadline = ()=>(
<div className="confidenceLegendTitle">
<span>AlphaFold confidence (pLDDT)</span>

<button
          className="confidenceInfoButton"
          onClick={()=>setShowConfidenceInfo(!showConfidenceInfo)}
          aria-label="About AlphaFold confidence"
>
          i
</button>
</div>
);

    const legendExplainer = ()=>(
<div className="confidenceInfo">
<strong>What does this mean?</strong>

<p>
          These colours show AlphaFold's confidence in the predicted local structure
          around each residue. Higher confidence does not guarantee biological
          correctness,and lower confidence does not automatically mean the prediction
          is wrong.
</p>

<button className="confidenceLearnMore">
          Learn more about confidence and how it is determined in The Basics →
</button>
</div>
);

    const meanConfidenceHead = ()=>(
<div className="proteinConfidenceHeader">
<div>
<h3>Mean predicted local confidence</h3>
</div>

<strong className = "proteinConfidenceValue">
          {modelScore?modelScore.mean.toFixed(1) : "—"}
<small>/100</small>
</strong>
</div>
);

    const nothingPickedPanel = ()=>(
<>
        {sideTab==="residue"&&selRes&&residuePanel()}

        {sideTab==="region"&&!selRegion&&(
<div className="residueInspector">
<h3>No region selected</h3>
<p className="residueSubtitle">
              Select a functional region from the protein information section to
              explore its biological role and highlight it in the 3D structure.
</p>
</div>
)}

        {sideTab==="region"&&selRegion&&regionPanel()}

        {sideTab==="atom"&&selRes&&atomPanel()}
</>
);

    return (
<div className="app">
<TopStrip activePage={pageNow} setActivePage={setPageNow}/>

<main className = "workspacePage">
<section className = "workspaceTop">
<div>
<h1>{protCard?.name||activeAcc}</h1>
<p>
                {protCard
?protCard.organismCommonName||protCard.organism
                  : proteinInfoLoading
?"Loading protein information..."
                    : proteinInfoError||"Protein information unavailable"}
</p>
</div>
</section>

<LensStrip activePage={pageNow} setActivePage={setPageNow}/>

<section className="proteinOverview">
<h2>What does this protein do?</h2>

<p>
<strong>UniProt ID:</strong>
              {activeAcc}
</p>

<p>
<strong>Sequence length:</strong>
              {protCard?.length ? `${protCard.length} amino acids` : "Unavailable"}
</p>

<h3>Biological function</h3>
<div className="functionDescription">
              {proteinInfoLoading?(
<p>Loading protein information...</p>
) : protCard?.function?functionLongList() : (
<p>No functional description available from UniProt.</p>
)}
</div>

<h3 className="functionalRegionsHeading">Annotated functional regions</h3>

<details className="functionalRegionsDetails">
<summary className="functionalRegionsToggle">
<span className="regionsShowText">
                  Show all {protCard?.regions?.length??0} regions ↓
</span>
<span className="regionsHideText">Show fewer regions ↑</span>
</summary>

              {protCard?.regions?.length?(
                protCard.regions.map((region,index)=>regionRowBit(index,region))
) : (
<p>No domain or region annotations available.</p>
)}
</details>
</section>

<section
            className="structureWorkspace"
            id="structure"
            ref={workspaceRef}
            style={{ "--panel-width": `${panelW}px` }}
>
<div className="structureViewer">
<div className="viewerHeading">
<div>
<h2>Explore the protein</h2>
</div>
</div>

<div className="structureSearch">
<div className="structureSearchHeading">
<strong>Find a residue or atom</strong>
</div>

                {findForm()}

                {jumpFeedback&&(
<p
                    className={
                      jumpFeedback.success?"structureSearchSuccess" : "structureSearchError"
}
>
                    {jumpFeedback.message}
</p>
)}
</div>

<div className="viewerShell">
                {foldedModelBox()}
</div>

<div className="confidenceLegend" id="confidence">
                {legendHeadline()}

<div className="plddtGradient"/>

<div className="plddtTicks">
<span>0</span>
<span>50</span>
<span>70</span>
<span>90</span>
<span>100</span>
</div>

<div className="plddtCategories">
<span>Very low</span>
<span>Low</span>
<span>Confident</span>
<span>Very high</span>
</div>
</div>

              {showConfidenceInfo && legendExplainer()}

<div className="proteinConfidenceSummary">
                {meanConfidenceHead()}

                {modelScore?confidenceBuckets() : (
<p>Waiting for confidence data from the AlphaFold model.</p>
)}
</div>
</div>

<div
              className="workspaceResizeHandle"
              role="separator"
              aria-label = "Resize Selection Details panel"
              aria-orientation="vertical"
              tabIndex={0}
              onPointerDown={grabHandle}
              onPointerMove={dragHandle}
              onPointerUp={dropHandle}
              onPointerCancel={dropHandle}
              onLostPointerCapture={dropHandle}
              onKeyDown={(event)=>{
                if (event.key==="ArrowLeft") {
                  event.preventDefault();
                  setPanelW((width)=>Math.min(width+20,800));
}

                if (event.key==="ArrowRight") {
                  event.preventDefault();
                  setPanelW((width)=>Math.max(width-20,280));
}
}}
/>

<aside className="evidencePanel">
<div className="inspectorTabs">
                {[
                  ["residue","Residue Explorer"],
                  ["atom","Atom Explorer"],
                  ["region","Region Explorer"],
].map(([tab,label])=>(
<button
                    key={tab}
                    type="button"
                    className={`inspectorTab ${sideTab === tab ? "active" : ""}`}
                    onClick={()=>setSideTab(tab)}
>
                    {label}
</button>
))}
</div>

              {selRes||(sideTab==="region"&&selRegion)?nothingPickedPanel() : (
<div className="emptySelection">
<h3>Select part of the protein</h3>
</div>
)}
</aside>
</section>
</main>
</div>
);
}

  
  

  const storyCards = [
    ["humanScene","/human-anatomy.png","Human anatomy","Your body depends on proteins."],
    ["cellScene","/animal-cell.png","Animal cell","Cells follow instructions stored in DNA."],
    ["sequenceScene","/amino-acids.png","Amino acids linked into a chain","A new protein starts as a chain of amino acids."],
    ["predictionScene","/alphafold-prediction.png","AlphaFold protein structure prediction","AlphaFold predicts a possible 3D structure."],
    ["uncertaintyScene","/alphafold-uncertainty.png","Different confidence levels in a predicted structure","Not every part of a prediction is equally reliable."],
];

  const dragZoomHints = ()=>(
<div className="viewerBasicsBar">
<span>
<strong>Drag</strong>to rotate
</span>
<span>
<strong>Scroll</strong>to zoom
</span>
<span>
<strong>Click</strong>to select
</span>
</div>
);

  const homeSearchField = ()=>(
<input
      type = "text"
      placeholder="Try insulin,TP53 or P04637..."
      value={searchBoxText}
      onChange={(event)=>{
        setHomeQuery(event.target.value);
        setShowProteinSuggestions(true);
        setActiveProteinSuggestionIndex(-1);
}}
      onKeyDown={(event)=>{
        if (event.key==="Escape") {
          setShowProteinSuggestions(false);
          setActiveProteinSuggestionIndex(-1);
          return;
}

        if (event.key==="ArrowDown"&&searchHits.length>0) {
          event.preventDefault();
          setShowProteinSuggestions(true);
          setActiveProteinSuggestionIndex((index)=>
            Math.min(index + 1,searchHits.length-1)
);
          return;
}

        if (event.key==="ArrowUp"&&searchHits.length>0) {
          event.preventDefault();
          setActiveProteinSuggestionIndex((index)=>Math.max(index-1,0));
          return;
}

        if (event.key==="Enter") {
          event.preventDefault();

          if (
            showProteinSuggestions&&
            activeProteinSuggestionIndex>=0&&
            searchHits[activeProteinSuggestionIndex]
) {
            visitProtein(
              searchHits[activeProteinSuggestionIndex].primaryAccession
);
} else {
            lookUpProtein();
}
}
}}
/>
);

  const homeSuggestionsDrop = ()=>(
<div className="proteinSuggestions">
      {searchHits.map((protein,index)=>{
        const accession = protein.primaryAccession;
        const name = pickTitle(protein);
        const scientificName =
          protein.organism?.scientificName||"Unknown organism";
        const commonName = protein.organism?.commonName||"";
        const organism = commonName
          ? `${scientificName} — ${commonName}`
          : scientificName;

        return (
<button
            key = {accession}
            type="button"
            className={`proteinSuggestion ${
              activeProteinSuggestionIndex===index?"isActive" : ""
            }`}
            onClick={()=>{
              setHomeQuery(name);
              visitProtein(accession);
}}
>
<strong>{name}</strong>
<span>
              {organism} · {accession}
</span>
</button>
);
})}
</div>
);

  const quickPicks = ()=>(
<div className="examples">
<span>Try an example</span>

<button className="exampleProtein" onClick={()=>visitProtein("P04637")}>
<strong>p53</strong>
</button>

<button className="exampleProtein" onClick={()=>visitProtein("P02144")}>
<strong>Myoglobin</strong>
</button>
</div>
);

  return (
<div className="app">
      {showIntro&&(
        <div className={`intro introStage${introStage}`}>
<button className="skipIntro" onClick={()=>setShowIntro(false)}>
            SKIP INTRO
</button>

<div className="introProgress" aria-label="Intro progress">
            {Array.from({length: 7},(_,index)=>(
<span key={index} className={introStage>=index?"active" : ""}/>
))}
</div>

<div className="introOpening">
<h1>Let’s show you how proteins are made.</h1>
</div>

          {storyCards.map(([sceneClass,src,alt,text])=>(
            <div key={sceneClass} className={`introScene ${sceneClass}`}>
<div className="introStoryVisual">
<img src={src} alt={alt} className="introStoryImage"/>
</div>
<div className="introCopy">
<h1>{text}</h1>
</div>
</div>
))}

<div className="introScene foldcheckScene">
<div className="introCopy">
<h1>
                Now you can explore the evidence...
<br/>
                with foldcheck!
</h1>
</div>
</div>
</div>
)}

      {!showIntro&&<TopStrip activePage={pageNow} setActivePage={setPageNow}/>}

<main>
        {pageNow==="home"&&(
<section className="hero newHomeHero" id="top">
<div className="heroContent">
<h1>Welcome to... FoldCheck.</h1>
<div className = "homeAction">
<span className="handNote">click here to start or scroll down!</span>
<svg className="handArrow" viewBox="0 0 110 65" fill="none" aria-hidden="true">
<path
                    d="M5 7 C30 1 54 12 68 30 S81 47 101 46 M85 34 L102 47 L84 57"
                    stroke="currentColor"
                    strokeWidth="2.4"
                    strokeLinecap="round"
                    strokeLinejoin="round"
/>
</svg>
<button className="primaryButton" onClick={()=>setPageNow("search")}>
                  Explore a protein →
</button>
<button className="secondaryButton" onClick={replayIntro}>
                  Replay intro ↻
</button>
</div>
</div>
</section>
)}

        {(pageNow==="home"||pageNow==="search")&&(
<section className="exploreSection" id="explore">
<div className="exploreLayout">
<div className="demoColumn">
<h2>Try the 3D viewer</h2>

<p className="sectionText">
                  Explore a real experimental protein structure before starting your own
                  investigation.
</p>

<div className="demoViewer" id="demo">
<div className="viewerShell">
<ProteinViewer/>
</div>

                  {dragZoomHints()}
</div>
</div>

<div className="searchColumn">
<h2>Which protein do you want to explore?</h2>

<p className="sectionText">
                  Search by protein name,gene,or UniProt accession. FoldCheck will retrieve
                  available structural data and help you investigate the prediction.
</p>

<div className="searchBox">
                  {homeSearchField()}

<button onClick={lookUpProtein}>Search</button>

                  {showProteinSuggestions&&searchHits.length>0&&homeSuggestionsDrop()}
</div>

                {proteinSearchError&&(
<p className="proteinSearchError" role="alert">
                    {proteinSearchError}
</p>
)}

                {quickPicks()}
</div>
</div>
</section>
)}

        {pageNow==="home"&&(
<section className="homeAbout">
<h2>Why I made FoldCheck</h2>
<p>
              I'm Yousef,a student who likes biology,chemistry and coding. After learning about
              AlphaFold,I started wondering how we can tell whether a predicted protein structure
              is reliable. That question turned into FoldCheck.
</p>
<button
              className="aboutLink"
              onClick={()=>{
                setPageNow("about");
                window.scrollTo(0,0);
}}
>
              Read my story →
</button>
</section>
)}
</main>
</div>
);
}

export default FoldBasics;

/* if you see this can i have freeesp32 :) */
