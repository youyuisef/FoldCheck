

import { useEffect, useRef, useState } from "react";
import "./App.css";
import ProteinViewer from "./ProteinViewer";
import aminoAcids from "./aminoAcids";
import atomChemistry from "./atomChemistry";
import RegionDetails from "./RegionDetails";
import PAEHeatmap from "./PAEHeatmap";

const AMINO_ACID_SUGGESTIONS = [
  ["ALA", "Alanine"],
  ["ARG", "Arginine"],
  ["ASN", "Asparagine"],
  ["ASP", "Aspartic acid"],
  ["CYS", "Cysteine"],
  ["GLN", "Glutamine"],
  ["GLU", "Glutamic acid"],
  ["GLY", "Glycine"],
  ["HIS", "Histidine"],
  ["ILE", "Isoleucine"],
  ["LEU", "Leucine"],
  ["LYS", "Lysine"],
  ["MET", "Methionine"],
  ["PHE", "Phenylalanine"],
  ["PRO", "Proline"],
  ["SER", "Serine"],
  ["THR", "Threonine"],
  ["TRP", "Tryptophan"],
  ["TYR", "Tyrosine"],
  ["VAL", "Valine"]
];

function alignResidueSequences(predicted, experimental) {
  const n = predicted.length;
  const m = experimental.length;

  if (n === 0 || m === 0) return [];

  const matchScore = 2;
  const mismatchScore = -1;
  const gapPenalty = -2;

  const scores = Array.from(
    { length: n + 1 },
    () => new Int32Array(m + 1)
  );

  for (let i = 1; i <= n; i++) {
    scores[i][0] = i * gapPenalty;
  }

  for (let j = 1; j <= m; j++) {
    scores[0][j] = j * gapPenalty;
  }

  for (let i = 1; i <= n; i++) {
    for (let j = 1; j <= m; j++) {
      const diagonal =
        scores[i - 1][j - 1] +
        (predicted[i - 1].code === experimental[j - 1].code
          ? matchScore
          : mismatchScore);

      const up = scores[i - 1][j] + gapPenalty;
      const left = scores[i][j - 1] + gapPenalty;

      scores[i][j] = Math.max(diagonal, up, left);
    }
  }

  const matches = [];
  let i = n;
  let j = m;

  while (i > 0 && j > 0) {
    const sameCode =
      predicted[i - 1].code === experimental[j - 1].code;

    const diagonal =
      scores[i - 1][j - 1] +
      (sameCode ? matchScore : mismatchScore);

    if (scores[i][j] === diagonal) {
      if (sameCode) {
        matches.push({
          predicted: predicted[i - 1],
          experimental: experimental[j - 1],
        });
      }

      i--;
      j--;
    } else if (
      scores[i][j] === scores[i - 1][j] + gapPenalty
    ) {
      i--;
    } else {
      j--;
    }
  }

  return matches.reverse();
}

function testRMSDCalculation() {
  const original = [
    [0, 0, 0],
    [2, 0, 0],
    [0, 3, 0],
    [0, 0, 4],
    [1, 2, 3],
  ];

  const transformed = original.map(([x, y, z]) => [
    -y + 10,
    x - 7,
    z + 5,
  ]);

  const pairs = original.map((predicted, index) => ({
    predicted,
    experimental: transformed[index],
  }));

  const rmsd = calculateAlignedRMSD(pairs);

  console.assert(
    rmsd !== null && Math.abs(rmsd) < 1e-6,
    "RMSD rotation/translation test FAILED:",
    rmsd
  );

  console.log(
    "FoldCheck RMSD validation:",
    rmsd !== null && Math.abs(rmsd) < 1e-6
      ? "PASS"
      : "FAIL",
    rmsd
  );
}

function calculateAlignedRMSD(pairs) {
  if (pairs.length < 3) return null;

  const n = pairs.length;

  const mean = (key) =>
    [0, 1, 2].map(
      (axis) =>
        pairs.reduce((sum, pair) => sum + pair[key][axis], 0) / n
    );

  const predictedCenter = mean("predicted");
  const experimentalCenter = mean("experimental");

  const covariance = Array.from({ length: 3 }, () => [0, 0, 0]);

  for (const pair of pairs) {
    const p = pair.predicted.map((v, i) => v - predictedCenter[i]);
    const q = pair.experimental.map((v, i) => v - experimentalCenter[i]);

    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        covariance[i][j] += p[i] * q[j];
      }
    }
  }

  const [
    [sxx, sxy, sxz],
    [syx, syy, syz],
    [szx, szy, szz],
  ] = covariance;

  const K = [
    [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
    [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
    [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
    [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
  ];

  
  const A = K.map((row) => [...row]);
  const V = Array.from({ length: 4 }, (_, i) =>
    Array.from({ length: 4 }, (_, j) => (i === j ? 1 : 0))
  );

  for (let iteration = 0; iteration < 100; iteration++) {
    let a = 0;
    let b = 1;
    let largest = 0;

    for (let i = 0; i < 4; i++) {
      for (let j = i + 1; j < 4; j++) {
        if (Math.abs(A[i][j]) > largest) {
          largest = Math.abs(A[i][j]);
          a = i;
          b = j;
        }
      }
    }

    if (largest < 1e-12) break;

    const angle =
      0.5 * Math.atan2(2 * A[a][b], A[b][b] - A[a][a]);

    const c = Math.cos(angle);
    const s = Math.sin(angle);

    for (let i = 0; i < 4; i++) {
      if (i === a || i === b) continue;

      const x = A[i][a];
      const y = A[i][b];

      A[i][a] = A[a][i] = c * x - s * y;
      A[i][b] = A[b][i] = s * x + c * y;
    }

    const aa = A[a][a];
    const bb = A[b][b];
    const ab = A[a][b];

    A[a][a] = c * c * aa - 2 * c * s * ab + s * s * bb;
    A[b][b] = s * s * aa + 2 * c * s * ab + c * c * bb;
    A[a][b] = A[b][a] = 0;

    for (let i = 0; i < 4; i++) {
      const x = V[i][a];
      const y = V[i][b];

      V[i][a] = c * x - s * y;
      V[i][b] = s * x + c * y;
    }
  }

  const eigenvalues = [0, 1, 2, 3].map((i) => A[i][i]);
  const best = eigenvalues.indexOf(Math.max(...eigenvalues));

  const [w, x, y, z] = V.map((row) => row[best]);

  const rotation = [
    [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
    [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
    [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
  ];

  let squaredError = 0;

  for (const pair of pairs) {
    const p = pair.predicted.map((v, i) => v - predictedCenter[i]);
    const q = pair.experimental.map((v, i) => v - experimentalCenter[i]);

    const rotated = rotation.map((row) =>
      row.reduce((sum, value, i) => sum + value * p[i], 0)
    );

    squaredError += rotated.reduce(
      (sum, value, i) => sum + (value - q[i]) ** 2,
      0
    );
  }

  return Math.sqrt(squaredError / n);
}
function calculateDihedral(p1, p2, p3, p4) {
  const subtract = (a, b) => a.map((v, i) => v - b[i]);
  const dot = (a, b) => a.reduce((s, v, i) => s + v * b[i], 0);
  const cross = (a, b) => [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];

  const b0 = subtract(p1, p2);
  const b1 = subtract(p3, p2);
  const b2 = subtract(p4, p3);

  const length = Math.sqrt(dot(b1, b1));
  if (length === 0) return null;

  const axis = b1.map((v) => v / length);
  const v = subtract(b0, axis.map((n) => dot(b0, axis) * n));
  const w = subtract(b2, axis.map((n) => dot(b2, axis) * n));

  const angle = Math.atan2(
    dot(cross(axis, v), w),
    dot(v, w)
  );

  return (angle * 180) / Math.PI;
}

function FoldCheckNav({ activePage, setActivePage }) {
  return (
    <header className="nav">
      <button
        className="brand brandButton"
        onClick={() => setActivePage("home")}
      >
        <div className="brandMark">F</div>
        <span>FoldCheck</span>
      </button>

      <div className="navLinks">
  <button className="navTextButton" onClick={() => setActivePage("home")}>
    Search
  </button>
  <button className="navTextButton" onClick={() => setActivePage("basics")}>
    Basics
  </button>
</div>

      <button className="navButton" onClick={() => setActivePage("home")}>
        Explore a protein
      </button>
    </header>
  );
}

const AMINO_CHEMISTRY = {
  ALA: ["Nonpolar", "Neutral", "Small"], ARG: ["Basic", "Positive", "Large"],
  ASN: ["Polar", "Neutral", "Medium"], ASP: ["Acidic", "Negative", "Medium"],
  CYS: ["Polar", "Neutral", "Small; disulfide-capable"], GLN: ["Polar", "Neutral", "Large"],
  GLU: ["Acidic", "Negative", "Large"], GLY: ["Nonpolar", "Neutral", "Very small; flexible"],
  HIS: ["Basic", "Variable near physiological pH", "Medium"],
  ILE: ["Nonpolar", "Neutral", "Large"], LEU: ["Nonpolar", "Neutral", "Large"],
  LYS: ["Basic", "Positive", "Large"], MET: ["Nonpolar", "Neutral", "Large"],
  PHE: ["Aromatic nonpolar", "Neutral", "Large"], PRO: ["Nonpolar", "Neutral", "Rigid"],
  SER: ["Polar", "Neutral", "Small"], THR: ["Polar", "Neutral", "Medium"],
  TRP: ["Aromatic nonpolar", "Neutral", "Very large"],
  TYR: ["Aromatic polar", "Neutral", "Large"], VAL: ["Nonpolar", "Neutral", "Medium"]
};
const ONE_TO_THREE = {
  A: "ALA", R: "ARG", N: "ASN", D: "ASP", C: "CYS", Q: "GLN", E: "GLU",
  G: "GLY", H: "HIS", I: "ILE", L: "LEU", K: "LYS", M: "MET", F: "PHE",
  P: "PRO", S: "SER", T: "THR", W: "TRP", Y: "TYR", V: "VAL"
};

function App() {
  const [geometryResidues, setGeometryResidues] = useState([]);
  const [experimentalResidues, setExperimentalResidues] = useState([]);
  const [linkedSelection, setLinkedSelection] = useState(null);
  const [mutationPosition, setMutationPosition] = useState("175");
  const [mutationTarget, setMutationTarget] = useState("HIS");

  const experimentalChains = [
  ...new Set(experimentalResidues.map((residue) => residue.chain)),
];

const selectedExperimentalChain =
  experimentalChains
    .map((chain) => {
      const chainResidues = experimentalResidues
        .filter(
          (residue) =>
            residue.chain === chain &&
            residue.backbone?.CA
        );

      const predictedResidues = geometryResidues
        .filter((residue) => residue.backbone?.CA);

      const matches = alignResidueSequences(
        predictedResidues,
        chainResidues
      );

      return {
        chain,
        matchedCount: matches.length,
        coverage:
          chainResidues.length > 0
            ? matches.length / chainResidues.length
            : 0,
      };
    })
    .sort(
      (a, b) =>
        b.matchedCount - a.matchedCount ||
        b.coverage - a.coverage
    )[0]?.chain ?? null;
const predictedSequence = geometryResidues
  .filter((residue) => residue.backbone?.CA)
  .map((residue) => ({
    number: residue.number,
    code: residue.code,
    chain: residue.chain,
    CA: residue.backbone.CA,
  }));

const experimentalSequence = experimentalResidues
  .filter(
    (residue) =>
      residue.chain === selectedExperimentalChain &&
      residue.backbone?.CA
  )
  .map((residue) => ({
    number: residue.number,
    code: residue.code,
    chain: residue.chain,
    CA: residue.backbone.CA,
  }));
  const matchedResidues = alignResidueSequences(
  predictedSequence,
  experimentalSequence
).map(({ predicted, experimental }) => ({
  number: predicted.number,
  code: predicted.code,
  predictedCA: predicted.CA,
  experimentalCA: experimental.CA,
  predictedNumber: predicted.number,
  predictedChain: predicted.chain,
  experimentalNumber: experimental.number,
  experimentalChain: experimental.chain,
}));
const handleLinkedSelection = (side, residue) => {
  const number = Number(residue?.residueNumber);
  const match = matchedResidues.find((item) =>
    side === "alphafold"
      ? item.predictedNumber === number && (!residue.chain || item.predictedChain === residue.chain)
      : item.experimentalNumber === number && (!residue.chain || item.experimentalChain === residue.chain)
  );
  if (!match) { setLinkedSelection(null); return; }
  setLinkedSelection({
    alphafold: { number: match.predictedNumber, chain: match.predictedChain },
    pdb: { number: match.experimentalNumber, chain: match.experimentalChain }
  });
};
const rmsdCoordinates = matchedResidues
  .filter(
    (residue) =>
      Array.isArray(residue.predictedCA) &&
      Array.isArray(residue.experimentalCA) &&
      residue.predictedCA.length === 3 &&
      residue.experimentalCA.length === 3 &&
      [...residue.predictedCA, ...residue.experimentalCA].every(
        Number.isFinite
      )
  )
  .map((residue) => ({
    predicted: residue.predictedCA,
    experimental: residue.experimentalCA,
  }));

const alignedRMSD =
  rmsdCoordinates.length >= 3
    ? calculateAlignedRMSD(rmsdCoordinates)
    : null;

  const [showIntro, setShowIntro] = useState(true);
  const [introStage, setIntroStage] = useState(0);
  const [activePage, setActivePage] = useState("home");
  useEffect(() => {
  window.scrollTo(0, 0);
}, [activePage]);
  const [proteinQuery, setProteinQuery] = useState("");
  const [proteinSearchLoading, setProteinSearchLoading] = useState(false);
const [proteinSearchError, setProteinSearchError] = useState("");
const [proteinSuggestions, setProteinSuggestions] = useState([]);
const [showProteinSuggestions, setShowProteinSuggestions] = useState(false);
const [activeProteinSuggestionIndex, setActiveProteinSuggestionIndex] = useState(-1);

useEffect(() => {
  const query = proteinQuery.trim();

  if (query.length < 2) {
    setProteinSuggestions([]);
    return;
  }

  const controller = new AbortController();

  const timer = setTimeout(async () => {
    try {
      const response = await fetch(
        `https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(
          `(${query}) AND reviewed:true`
        )}&fields=accession,protein_name,gene_names,organism_name&size=8`,
        { signal: controller.signal }
      );

      if (!response.ok) return;

      const data = await response.json();

      if (!controller.signal.aborted) {
        setProteinSuggestions(data.results || []);
      }
    } catch (error) {
      if (error.name !== "AbortError") {
        setProteinSuggestions([]);
      }
    }
  }, 350);

  return () => {
    clearTimeout(timer);
    controller.abort();
  };
}, [proteinQuery]);

const openProtein = async (accession) => {
  setProteinSearchLoading(true);
  setProteinSearchError("");

  try {
    const response = await fetch(
      `https://alphafold.ebi.ac.uk/api/prediction/${encodeURIComponent(accession)}`
    );

    if (response.status === 404) {
      setProteinSearchError(
        "This protein exists in UniProt, but no AlphaFold prediction is available."
      );
      return;
    }

    if (!response.ok) {
      throw new Error("AlphaFold availability check failed");
    }

    setShowProteinSuggestions(false);
setActiveProteinSuggestionIndex(-1);
setSelectedResidue(null);
setInspectorTab("residue");
setActiveProteinId(accession);
setActivePage("workspace");
  } catch {
    setProteinSearchError(
      "Could not check AlphaFold availability. Please try again."
    );
  } finally {
    setProteinSearchLoading(false);
  }
};

const searchProtein = async () => {
  const query = proteinQuery.trim();
  if (!query || proteinSearchLoading) return;

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
      (protein) => protein.primaryAccession === accession
    );

    const protein = exactMatch || data.results?.[0];

    if (!protein) {
      setProteinSearchError("No matching protein found in UniProt.");
      return;
    }

await openProtein(protein.primaryAccession);
  } catch {
    setProteinSearchError("Could not search UniProt. Please try again.");
  } finally {
    setProteinSearchLoading(false);
  }
};
  const [activeProteinId, setActiveProteinId] = useState("P04637");
  useEffect(() => {
  setGeometryResidues([]);
}, [activeProteinId]);
  const [proteinInfo, setProteinInfo] = useState(null);
  const [experimentalStructures, setExperimentalStructures] = useState([]);
const [experimentalLoading, setExperimentalLoading] = useState(false);
const [experimentalError, setExperimentalError] = useState("");
const [selectedExperimentalId, setSelectedExperimentalId] = useState("");
  const [partnerInfo, setPartnerInfo] = useState({});
  const [selectedRegion, setSelectedRegion] = useState(null);
const [interactionResearch, setInteractionResearch] = useState([]);
const [interactionResearchLoading, setInteractionResearchLoading] = useState(false);

useEffect(() => {
  const partner = proteinInfo?.interactionPartners?.find(
    (item) =>
      item.name.length > 2 &&
      selectedRegion?.name?.toLowerCase().includes(item.name.toLowerCase())
  );

  if (!partner || !selectedRegion) {
    setInteractionResearch([]);
    setInteractionResearchLoading(false);
    return;
  }

  const controller = new AbortController();

  const loadResearch = async () => {
    setInteractionResearch([]);
    setInteractionResearchLoading(true);

    try {
      const proteinName = proteinInfo?.gene || activeProteinId;
      const regionName = selectedRegion?.name || "";

      const query = `TITLE_ABS:("${proteinName}" AND "${partner.name}" AND (interaction OR binding OR complex OR mechanism))`;

      const response = await fetch(
        `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=${encodeURIComponent(query)}&format=json&resultType=core&pageSize=5`,
        { signal: controller.signal }
      );

      if (!response.ok) throw new Error("Research search failed");

      const data = await response.json();

      if (!controller.signal.aborted) {
        const papers = data.resultList?.result || [];
        const papersWithAbstracts = await Promise.all(
  papers.map(async (paper) => {
    if (!paper.id || !paper.source) return paper;

    try {

const abstractResponse = await fetch(
  `https://www.ebi.ac.uk/europepmc/webservices/rest/${encodeURIComponent(paper.source)}/${encodeURIComponent(paper.id)}?format=json`,
  { signal: controller.signal }
);

if (!abstractResponse.ok) {
  return {
    ...paper,
    abstract: (paper.abstractText || "")
      .replace(/<[^>]*>/g, " ")
      .replace(/\s+/g, " ")
      .trim()
  };
}

const article = await abstractResponse.json();

return {
  ...paper,
  abstract: (article.result?.abstractText || paper.abstractText || "")
  .replace(/<[^>]*>/g, " ")
  .replace(/\s+/g, " ")
  .trim()
};
    } catch {
  return {
    ...paper,
    abstract: (paper.abstractText || "")
      .replace(/<[^>]*>/g, " ")
      .replace(/\s+/g, " ")
      .trim()
  };
}
  })
);

if (!controller.signal.aborted) {
  setInteractionResearch(
  papersWithAbstracts.filter((paper) => {
    const abstract = (paper.abstract || "").toLowerCase();

    const protein = (proteinInfo?.gene || "").toLowerCase();
    const partnerName = partner.name.toLowerCase();

    return (
      abstract.length > 100 &&
      protein.length > 0 &&
      [
  proteinInfo?.gene,
  proteinInfo?.name,
  activeProteinId
].filter((name) => name && name.length > 2)
  .some((name) => abstract.includes(name.toLowerCase())) &&
      abstract.includes(partnerName)
    );
  })
);
}
      }
    } catch {
      if (!controller.signal.aborted) {
        setInteractionResearch([]);
      }
    } finally {
      if (!controller.signal.aborted) {
        setInteractionResearchLoading(false);
      }
    }
  };

  loadResearch();

  return () => controller.abort();
}, [selectedRegion, proteinInfo, activeProteinId]);

useEffect(() => {
 const partnerIds = [
  ...new Set(
    (proteinInfo?.interactionPartners || [])
      .filter((partner) => partner.id)
      .map((partner) => partner.id)
  )
];

const regionPartnerNames = (proteinInfo?.regions || [])
  .map((region) =>
    region.name.match(/^Interaction with\s+(.+)$/i)?.[1]?.trim()
  )
  .filter(Boolean);

const missingPartnerNames = regionPartnerNames.filter(
  (name) =>
    !(proteinInfo?.interactionPartners || []).some(
      (partner) => partner.name.toLowerCase() === name.toLowerCase()
    )
);
  const controller = new AbortController();

  const loadPartners = async () => {
    const matchingPartners = partnerIds.filter((id) =>
      (proteinInfo?.regions || []).some((region) =>
        (proteinInfo?.interactionPartners || []).some(
          (partner) =>
            partner.id === id &&
            partner.name.length > 2 &&
            region.name.toLowerCase().includes(partner.name.toLowerCase())
        )
      )
    );
const missingPartnerResults = await Promise.all(
  missingPartnerNames.slice(0, 10).map(async (name) => {
    try {
      const query = `gene_exact:${name} AND organism_id:9606 AND reviewed:true`;

      const response = await fetch(
        `https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(query)}&size=1`,
        { signal: controller.signal }
      );

      if (!response.ok) return null;

      const data = await response.json();
      const protein = data.results?.[0];

      if (!protein) return null;

      return {
        id: protein.primaryAccession,
        gene: name,
        name:
          protein.proteinDescription?.recommendedName?.fullName?.value ||
          name,
        function:
          protein.comments
            ?.filter((comment) => comment.commentType === "FUNCTION")
            .flatMap((comment) => comment.texts || [])
            .map((text) => text.value)
            .join(" ") || ""
      };
    } catch {
      return null;
    }
  })
);

    const results = await Promise.all(
      matchingPartners.slice(0, 10).map(async (id) => {
        try {
          const response = await fetch(
            `https://rest.uniprot.org/uniprotkb/${encodeURIComponent(id)}.json`,
            { signal: controller.signal }
          );

          if (!response.ok) return null;

          const data = await response.json();

          return {
            id,
            name:
              data.proteinDescription?.recommendedName?.fullName?.value ||
              id,
            function:
              data.comments
                ?.filter((comment) => comment.commentType === "FUNCTION")
                .flatMap((comment) => comment.texts || [])
                .map((text) => text.value)
                .join(" ") || ""
          };
        } catch {
          return null;
        }
      })
    );

    if (!controller.signal.aborted) {
      setPartnerInfo(
  Object.fromEntries(
    [...results, ...missingPartnerResults]
      .filter(Boolean)
      .flatMap((partner) => [
        [partner.id, partner],
        ...(partner.gene
          ? [[partner.gene.toUpperCase(), partner]]
          : [])
      ])
  )
);
    }
  };

  setPartnerInfo({});
  loadPartners();

  return () => controller.abort();
}, [proteinInfo]);
  const [proteinInfoLoading, setProteinInfoLoading] = useState(false);
  const [proteinInfoError, setProteinInfoError] = useState("");
  const [selectedResidue, setSelectedResidue] = useState(null);
  const [proteinConfidence, setProteinConfidence] = useState(null);
  const [paeData, setPaeData] = useState(null);
  const [paeReferenceResidue, setPaeReferenceResidue] = useState(null);
  const [paeSelectedPair, setPaeSelectedPair] = useState(null);
const [paeLoading, setPaeLoading] = useState(false);
const [paeError, setPaeError] = useState("");
useEffect(() => {
  const controller = new AbortController();

  const loadPAE = async () => {
  setPaeData(null);
  setPaeSelectedPair(null);
  setPaeReferenceResidue(null);
  setPaeError("");
  setPaeLoading(true);

    try {
      const response = await fetch(
        `https://alphafold.ebi.ac.uk/files/AF-${activeProteinId}-F1-predicted_aligned_error_v6.json`,
        { signal: controller.signal }
      );

      if (!response.ok) {
        throw new Error("PAE data unavailable");
      }

      const data = await response.json();
      const entry = Array.isArray(data) ? data[0] : data;

let matrix = entry?.predicted_aligned_error;

if (!matrix && Array.isArray(entry?.residue1) && Array.isArray(entry?.distance)) {
  const size = Math.max(...entry.residue1, ...entry.residue2);

  matrix = Array.from({ length: size }, () =>
    Array(size).fill(30)
  );

  for (let i = 0; i < entry.distance.length; i++) {
    const row = entry.residue1[i] - 1;
    const col = entry.residue2[i] - 1;

    if (row >= 0 && col >= 0 && row < size && col < size) {
      matrix[row][col] = entry.distance[i];
    }
  }
}

      if (!Array.isArray(matrix) || !Array.isArray(matrix[0])) {
        throw new Error("Invalid PAE data");
      }

      if (!controller.signal.aborted) {
        setPaeData(matrix);
      }
    } catch (error) {
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

  return () => controller.abort();
}, [activeProteinId]);
  const [panelWidth, setPanelWidth] = useState(320);
  const [showConfidenceInfo, setShowConfidenceInfo] = useState(false);
  const [inspectorTab, setInspectorTab] = useState("residue");
  const [structureQuery, setStructureQuery] = useState("");
  const [proteinResidues, setProteinResidues] = useState([]);
  const [showStructureSuggestions, setShowStructureSuggestions] = useState(false);
  const [activeSuggestionIndex, setActiveSuggestionIndex] = useState(-1);

const structureSuggestions = (() => {
  const query =  String(structureQuery ?? "").toUpperCase().trimStart();

  if (!query) return [];

const atomMatch = query.match(/^([A-Z]{3})\s+(\d+)\s+([A-Z0-9']*)$/);

if (atomMatch) {
  const [, code, number, atomInput] =atomMatch;

  const residue = proteinResidues.find(
    (item) =>
      item.code === code &&
    item.number === Number(number)
  );

  if (!residue) return [];

  return (residue.atoms || [])
  .filter((atom) => atom.startsWith(atomInput))
  .slice(0, 12)
  .map((atom) => [
    `${code} ${number} ${atom}`,
    `Atom · ${residue.chain}`
  ]);
}

  const match = query.match(/^([A-Z]+)\s+(\d*)$/);

  if (match) {
    const [, inputCode, inputNumber] = match;
    const amino = AMINO_ACID_SUGGESTIONS.find(
      ([code, name]) =>
        code === inputCode ||
      name.toUpperCase() === inputCode
    );

    if (!amino) return [];

    return proteinResidues
    .filter(
      (residue) =>
        residue.code === amino[0] &&
      String(residue.number).startsWith(inputNumber)
    )

    .slice(0, 8)
    .map((residue) => [
      `${residue.code} ${residue.number}`,
      `Chain ${residue.chain}`
    ]);
  }

return AMINO_ACID_SUGGESTIONS.filter(
  ([code, name]) =>
    code.startsWith(query) ||
  name.toUpperCase().startsWith(query)
).slice(0, 6);

})();

const [searchRequest, setSearchRequest] = useState(null);
const [resetViewRequest, setResetViewRequest] = useState(0);
const [searchFeedback, setSearchFeedback] = useState(null);

  const selectedAtomName = selectedResidue?.atom?.trim().toUpperCase() || "";

const backboneAtoms = ["N", "CA", "C", "O", "OXT"];

const atomElement =
  selectedAtomName === "N" || selectedAtomName.startsWith("N")
    ? "Nitrogen"
    : selectedAtomName === "O" || selectedAtomName.startsWith("O")
      ? "Oxygen"
      : selectedAtomName === "S" || selectedAtomName.startsWith("S")
        ? "Sulfur"
        : selectedAtomName.startsWith("C")
          ? "Carbon"
          : selectedAtomName.startsWith("H")
            ? "Hydrogen"
            : null;

const atomLocation = !selectedAtomName
  ? null
  : backboneAtoms.includes(selectedAtomName)
    ? "Protein backbone"
    : "Amino acid side chain";

const specificAtomInfo =
  atomChemistry[selectedResidue?.residueCode]?.[selectedAtomName] || null;

const atomRole = !selectedAtomName
  ? null
  : selectedAtomName === "N"
    ? "This backbone nitrogen is part of the peptide chain and can participate in hydrogen bonding."
    : selectedAtomName === "CA"
      ? "The alpha carbon connects the backbone to the amino acid's side chain. It is a key reference point for protein geometry."
      : selectedAtomName === "C"
        ? "This backbone carbon is part of the peptide group and is bonded to the carbonyl oxygen."
        : selectedAtomName === "O"
          ? "This backbone oxygen is part of a carbonyl group and can accept hydrogen bonds."
          : backboneAtoms.includes(selectedAtomName)
            ? "This atom belongs to the protein backbone."
            : "This atom belongs to the amino acid's side chain. Its chemical role depends on the residue and its molecular environment.";

useEffect(() => {
  const controller = new AbortController();

  const loadProteinInfo = async () => {
    setProteinInfo(null);
    setProteinInfoError("");
    setProteinInfoLoading(true);
    setExperimentalStructures([]);
setSelectedExperimentalId("");
setExperimentalError("");
setExperimentalLoading(true);

    try {
      const response = await fetch(
  `https://rest.uniprot.org/uniprotkb/${encodeURIComponent(activeProteinId)}.json`,
  { signal: controller.signal }
);
      if (!response.ok) {
        throw new Error("Protein information unavailable");
      }

      const data = await response.json();

      if (controller.signal.aborted) return;

const pdbReferences = (data.uniProtKBCrossReferences || [])
  .filter((reference) => reference.database === "PDB")
  .map((reference) => ({
    id: reference.id,
    method:
      reference.properties?.find(
        (property) => property.key === "Method"
      )?.value || "Not specified",
    resolution:
      reference.properties?.find(
        (property) => property.key === "Resolution"
      )?.value || "Not specified",
    chains:
      reference.properties?.find(
        (property) => property.key === "Chains"
      )?.value || "Not specified"
  }));

setExperimentalStructures(pdbReferences);
setSelectedExperimentalId(pdbReferences[0]?.id || "");

      setProteinInfo({
        name:
          data.proteinDescription?.recommendedName?.fullName?.value ||
          data.proteinDescription?.submissionNames?.[0]?.fullName?.value ||
          activeProteinId,
        gene: data.genes?.[0]?.geneName?.value || "",
        aliases: [
  ...(data.genes?.[0]?.synonyms || []).map((item) => item.value),
  ...(data.proteinDescription?.recommendedName?.shortNames || []).map(
    (item) => item.value
  )
].filter(Boolean),
        organism: data.organism?.scientificName || "",
organismCommonName: data.organism?.commonName || "",
        length: data.sequence?.length || null,
        sequence: data.sequence?.value || "",
        subcellularLocations: (data.comments || [])
          .filter((comment) => comment.commentType === "SUBCELLULAR LOCATION")
          .flatMap((comment) => comment.subcellularLocations || [])
          .map((item) => ({
            name: item.location?.value || "Unknown",
            topology: item.topology?.value || "",
            orientation: item.orientation?.value || ""
          }))
          .filter((item) => item.name),
        function:
          data.comments?.find((item) => item.commentType === "FUNCTION")
            ?.texts?.map((item) => item.value)
  .join(" ")
  .replace(/\s*\(PubMed:[^)]+\)/g, "")
  .replace(/\s+/g, " "),

interactionPartners: (data.comments || [])
  .filter((comment) => comment.commentType === "INTERACTION")
  .flatMap((comment) => comment.interactions || [])
  .map((interaction) => ({
    name: interaction.interactantTwo?.geneName || "",
id: interaction.interactantTwo?.uniProtKBAccession || "",
experiments: interaction.numberOfExperiments ?? null,
intActId: interaction.interactantTwo?.intActId || "",
organismDiffer: interaction.organismDiffer ?? false
  })),

  interactionComments: (data.comments || [])
  .filter((comment) => comment.commentType === "INTERACTION")
  .flatMap((comment) => comment.interactions || []),
        regions: (data.features || [])
          .filter((feature) =>
            ["Domain", "Region", "Chain", "Peptide", "Propeptide", "Signal"].includes(feature.type)
          )
          .map((feature) => ({
 name: feature.description || feature.type,
type: feature.type,
description: [
  feature.description,
  ...(feature.ligand ? [`Ligand: ${feature.ligand.name || feature.ligand.id || "specified"}`] : [])
].filter(Boolean).join(". "),
  start: feature.location?.start?.value,
  end: feature.location?.end?.value,
  evidence: feature.evidences || []
}))
      });
        } catch (error) {
      if (error.name === "AbortError") return;

      setProteinInfoError(
        "Could not retrieve protein information from UniProt."
      );
      setExperimentalError(
        "Could not retrieve experimental structure references."
      );
    } finally {
      if (!controller.signal.aborted) {
        setProteinInfoLoading(false);
        setExperimentalLoading(false);
      }
    }
  };

  loadProteinInfo();

  return () => controller.abort();
}, [activeProteinId]);

  const workspaceRef = useRef(null);
  const resizingRef = useRef(false);

  const startPanelResize = (event) => {
    event.preventDefault();
    resizingRef.current = true;
    event.currentTarget.setPointerCapture(event.pointerId);
  };

  const movePanelResize = (event) => {
    if (!resizingRef.current || !workspaceRef.current) return;

    const bounds = workspaceRef.current.getBoundingClientRect();
    const availableWidth = bounds.width - 10;
    const requestedWidth = bounds.right - event.clientX;
    const maximumWidth = Math.max(280, availableWidth - 350);

    setPanelWidth(Math.max(280, Math.min(requestedWidth, maximumWidth)));
  };

  const stopPanelResize = () => {
    resizingRef.current = false;
  };

  const selectedPlddt =
    typeof selectedResidue?.plddt === "number"
      ? selectedResidue.plddt
      : null;

  const confidenceCategory =
    selectedPlddt === null
      ? null
      : selectedPlddt >= 90
        ? "Very high"
        : selectedPlddt >= 70
          ? "Confident"
          : selectedPlddt >= 50
            ? "Low"
            : "Very low";

  const confidenceExplanation =
    selectedPlddt === null
      ? null
      : selectedPlddt >= 90
        ? "Very high confidence suggests AlphaFold has strong support for this residue's predicted local structure. A well-constrained folded region may contribute, but the score alone cannot establish why confidence is high."
        : selectedPlddt >= 70
          ? "AlphaFold predicts this local region with generally good confidence, although some structural details may remain uncertain. Flexibility or fewer structural constraints are possible explanations, not confirmed causes."
          : selectedPlddt >= 50
            ? "Lower confidence may reflect flexibility, an exposed loop, intrinsic disorder, or limited structural information. Additional evidence is needed to determine which explanation applies."
            : "Very low confidence may indicate a flexible or intrinsically disordered region, or a region AlphaFold could not model confidently. It does not prove the residue is disordered or that the prediction is incorrect.";

  const confidenceSignificance =
    selectedPlddt === null
      ? null
      : selectedPlddt >= 90
        ? "The local backbone arrangement is predicted with high confidence, but this does not guarantee accurate side-chain positions, biological interactions, or experimental agreement."
        : selectedPlddt >= 70
          ? "The predicted local fold may be useful for interpretation, but detailed molecular interactions should still be evaluated cautiously."
          : selectedPlddt >= 50
            ? "Avoid relying on precise atomic positions or specific contacts in this region without additional supporting evidence."
            : "Treat the displayed coordinates as highly uncertain. Avoid drawing strong conclusions about this residue's exact position, contacts, or interactions.";

  const residuePosition = Number(selectedResidue?.residueNumber);

  const p53Region =
    !Number.isInteger(residuePosition) ||
    residuePosition < 1 ||
    residuePosition > 393
      ? null
      : residuePosition <= 61
        ? {
            name: "N-terminal transactivation region",
            description:
              "This region helps p53 interact with proteins involved in regulating gene expression. Much of the N-terminal region is intrinsically disordered."
          }
        : residuePosition <= 93
          ? {
              name: "Proline-rich region",
              description:
                "This region contains many proline residues and contributes to interactions and regulatory functions of p53."
            }
          : residuePosition <= 292
            ? {
                name: "DNA-binding core domain",
                description:
                  "This folded domain is central to p53's ability to recognise DNA. Many cancer-associated TP53 mutations occur in this region."
              }
            : residuePosition <= 325
              ? {
                  name: "Linker region",
                  description:
                    "This region connects the DNA-binding core with the oligomerisation domain and may allow flexibility between them."
                }
              : residuePosition <= 356
                ? {
                    name: "Oligomerisation domain",
                    description:
                      "This region contributes to the formation of p53 tetramers, which are important for its biological activity."
                  }
                : {
                    name: "C-terminal regulatory region",
                    description:
                      "This region participates in regulating p53 activity and contains sites of post-translational modification. It is often flexible."
                  };

  const confidenceDifference =
    selectedPlddt !== null && proteinConfidence
      ? selectedPlddt - proteinConfidence.mean
      : null;

  const selectedAminoAcid = selectedResidue
    ? aminoAcids[selectedResidue.residueCode]
    : null;

  useEffect(() => {
    if (!showIntro) return;

    const timers = [
      setTimeout(() => setIntroStage(1), 6500),
      setTimeout(() => setIntroStage(2), 13500),
      setTimeout(() => setIntroStage(3), 20500),
      setTimeout(() => setIntroStage(4), 28000),
      setTimeout(() => setIntroStage(5), 36000),
      setTimeout(() => setShowIntro(false), 43000)
    ];

    return () => timers.forEach(clearTimeout);
  }, [showIntro]);

const findStructureTarget = () => {
  const input = structureQuery.trim().toUpperCase();

  const match = input.match(
    /^(?:([A-Z]{3})\s+)?(\d+)(?:\s+([A-Z][A-Z0-9']*))?$/
  );

  if (!match) {
    setSearchFeedback({
      success: false,
      message: "Try 158, ARG 158, or ARG 158 CD."
    });
    return;
  }

  const residueNumber = Number(match[2]);

  if (!proteinResidues.some((residue) => residue.number === residueNumber)) {
    setSearchFeedback({
      success: false,
      message: "That residue number is not present in the loaded protein structure."
    });
    return;
  }

  setSearchFeedback(null);

  setSearchRequest({
    residueCode: match[1] || null,
    residueNumber,
    atomName: match[3] || null,
    requestId: Date.now()
  });
};

const handleStructureSearchResult = (result) => {
  setSearchFeedback(result);

  if (result.success) {
    setInspectorTab(result.atomSearch ? "atom" : "residue");
  }
};

  const replayIntro = () => {
    setIntroStage(0);
    setShowIntro(true);
  };

  const paeStats = (() => {
  if (!paeData?.length) return null;

  let total = 0;
  let count = 0;
  let low = 0;
  let moderate = 0;
  let high = 0;
  let maximum = 0;

  for (let i = 0; i < paeData.length; i++) {
    const row = paeData[i];

    if (!Array.isArray(row)) continue;

    for (let j = 0; j < row.length; j++) {
      if (i === j) continue;

      const value = Number(row[j]);
      if (!Number.isFinite(value)) continue;

      total += value;
      count++;
      maximum = Math.max(maximum, value);

      if (value < 5) low++;
      else if (value < 15) moderate++;
      else high++;
    }
  }

  if (!count) return null;

  return {
    proteinLength: paeData.length,
    mean: total / count,
    maximum,
    lowPercent: (low / count) * 100,
    moderatePercent: (moderate / count) * 100,
    highPercent: (high / count) * 100,
    pairCount: count
  };
})();

if (activePage === "confidence") {
  return (
    <div className="app">
      <FoldCheckNav activePage={activePage} setActivePage={setActivePage} />

      <main className="workspacePage">
        <section className="workspaceTop">
          <div>
            <p className="sectionLabel">STRUCTURAL CONFIDENCE ANALYSIS</p>
            <h1>Local Prediction Confidence</h1>
            <p>{proteinInfo?.name || activeProteinId} · {activeProteinId}</p>
          </div>
        </section>

        <nav className="analysisTabs">
  <button
    className="analysisTab"
    onClick={() => setActivePage("workspace")}
  >
    3D Structure
  </button>

  <button className="analysisTab active">
    Confidence
  </button>

  <button
    className="analysisTab"
    onClick={() => setActivePage("pae")}
  >
    PAE
  </button>

  <button className="analysisTab" onClick={() => setActivePage("geometry")}>
    Geometry
  </button>

  <button className="analysisTab" onClick={() => setActivePage("experimental")}>
    Experimental Comparison
  </button>
</nav>

        <section className="paeSection">
          <h2>AlphaFold Confidence (pLDDT)</h2>
          <p>
            pLDDT estimates confidence in the predicted local structure
            of each amino acid, on a scale from 0 to 100.
          </p>

          <div className="confidenceLegend">
            <div className="plddtGradient" />
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
    overflow: "hidden"
  }}
>
  <ProteinViewer
    source="alphafold"
    proteinId={activeProteinId}
    onProteinConfidence={setProteinConfidence}
    onResidueSelect={setSelectedResidue}
  />
</div>

{proteinConfidence && (
  <div style={{ marginBottom: "28px" }}>
    <h2>Confidence Distribution</h2>
    <p>
      Percentage of residues in each AlphaFold confidence category.
    </p>

    {[
      ["Very high (90–100)", proteinConfidence.veryHigh, "#0053d6"],
      ["Confident (70–90)", proteinConfidence.confident, "#65cbf3"],
      ["Low (50–70)", proteinConfidence.low, "#ffdb13"],
      ["Very low (0–50)", proteinConfidence.veryLow, "#ff7d45"]
    ].map(([label, count, color]) => {
      const percentage =
        (count / proteinConfidence.residueCount) * 100;

      return (
        <div key={label} style={{ marginBottom: "18px" }}>
          <div style={{
            display: "flex",
            justifyContent: "space-between",
            marginBottom: "8px"
          }}>
            <strong>{label}</strong>
            <span>{percentage.toFixed(1)}%</span>
          </div>

          <div style={{
            height: "14px",
            background: "#1c303b",
            borderRadius: "8px",
            overflow: "hidden"
          }}>
            <div style={{
              width: `${percentage}%`,
              height: "100%",
              background: color,
              borderRadius: "8px"
            }} />
          </div>
        </div>
      );
    })}
  </div>
)}

          <h3>Whole-protein confidence</h3>
          {proteinConfidence ? (
            <>
              <h2>{proteinConfidence.mean.toFixed(1)} / 100</h2>
              <p>
                Mean pLDDT across {proteinConfidence.residueCount} residues.
              </p>
            </>
          ) : (
            <p>Confidence data will load with the protein structure.</p>
          )}
        </section>
      </main>
    </div>
  );
}

  if (activePage === "pae") {
  return (
    <div className="app">
      <FoldCheckNav activePage={activePage} setActivePage={setActivePage} />



      <main className="workspacePage">
        <section className="workspaceTop">
          <div>
            <p className="sectionLabel">STRUCTURAL CONFIDENCE ANALYSIS</p>
            <h1>Predicted Aligned Error</h1>
            <p>
              {proteinInfo?.name || activeProteinId} · {activeProteinId}
            </p>
          </div>
        </section>

<nav className="analysisTabs">
  <button
    className="analysisTab"
    onClick={() => setActivePage("workspace")}
  >
    3D Structure
  </button>

  <button
  className="analysisTab"
  onClick={() => setActivePage("confidence")}
>
  Confidence
</button>

  <button className="analysisTab active">
  PAE
</button>

  <button className="analysisTab" onClick={() => setActivePage("geometry")}>
    Geometry
  </button>

  <button className="analysisTab" onClick={() => setActivePage("experimental")}>
    Experimental Comparison
  </button>
</nav>

        <div
  style={{
    marginBottom: "18px",
    padding: "12px 18px",
    background: "#10242b",
    border: "1px solid #28535a",
    borderRadius: "10px",
    textAlign: "center"
  }}
>
  <a
    href="#pae-explanation"
    style={{
      color: "#a8e8e0",
      fontSize: "14px",
      fontWeight: 600,
      textDecoration: "none"
    }}
  >
    ↓ Scroll down to understand PAE and explore the scientific explanations
  </a>
</div>

        <section
          className="paeSection"
          style={{
            padding: "28px",
            background: "#0d141b",
            border: "1px solid #263843",
            borderRadius: "16px"
          }}
        >
          <div
  style={{
    display: "grid",
    gridTemplateColumns: "repeat(auto-fit, minmax(min(100%, 450px), 1fr))",
    gap: "32px",
    alignItems: "start"
  }}
>
  <div style={{ minWidth: 0 }}>
    {paeLoading ? (
      <p>Loading AlphaFold PAE data...</p>
    ) : paeData ? (
      <PAEHeatmap
  key={activeProteinId}
  matrix={paeData}
  onReferenceChange={setPaeReferenceResidue}
  onPairSelect={setPaeSelectedPair}
/>
    ) : paeError ? (
      <p>{paeError}</p>
    ) : null}
  </div>

  <div style={{ minWidth: 0 }}>
    <h2>3D Protein Structure</h2>

    <p style={{ color: "#9fb8c7" }}>
      Explore the predicted structure alongside its relative-position
      uncertainty.
    </p>

    <div
      className="viewerShell"
      style={{
        height: "650px",
        minHeight: "650px",
        overflow: "hidden",
        borderRadius: "12px"
      }}
    >
      <ProteinViewer
  source="alphafold"
  proteinId={activeProteinId}
  onResidueSelect={setSelectedResidue}
  paeMatrix={paeData}
  paeReferenceResidue={paeReferenceResidue}
/>
    </div>

    <p style={{ color: "#9fb8c7", fontSize: "13px" }}>
      The 3D structure is coloured by predicted aligned error relative to your selected reference residue. Click different positions on the heatmap to update the colouring. Blue indicates lower predicted error, while yellow indicates higher predicted error.
    </p>
  </div>
</div>
                </section>

{paeStats && (
  <section
    className="paeSection"
    style={{
      marginTop: "24px",
      padding: "28px",
      background: "#0d141b",
      border: "1px solid #263843",
      borderRadius: "16px"
    }}
  >

    <p className="sectionLabel">PROTEIN-SPECIFIC ANALYSIS</p>
    <h2>PAE Confidence Statistics</h2>

    <p style={{ color: "#9fb8c7" }}>
      Calculated from the AlphaFold PAE matrix for{" "}
      <strong>{proteinInfo?.name || activeProteinId}</strong>{" "}
      ({activeProteinId}).
    </p>

    <div
      style={{
        display: "grid",
        gridTemplateColumns:
          "repeat(auto-fit, minmax(190px, 1fr))",
        gap: "16px",
        marginTop: "24px"
      }}
    >
      {[
        {
          label: "Mean PAE",
          value: `${paeStats.mean.toFixed(2)} Å`,
          detail: "Average across non-self residue pairs"
        },
        {
          label: "Low Error",
          value: `${paeStats.lowPercent.toFixed(1)}%`,
          detail: "Predicted error below 5 Å"
        },
        {
          label: "Moderate Error",
          value: `${paeStats.moderatePercent.toFixed(1)}%`,
          detail: "Predicted error from 5 to below 15 Å"
        },
        {
          label: "High Error",
          value: `${paeStats.highPercent.toFixed(1)}%`,
          detail: "Predicted error of 15 Å or greater"
        },
        {
          label: "Maximum PAE",
          value: `${paeStats.maximum.toFixed(2)} Å`,
          detail: "Highest value in the analysed matrix"
        },
        {
          label: "Residue Pairs",
          value: paeStats.pairCount.toLocaleString(),
          detail: "Ordered pairs, excluding self-comparisons"
        }
      ].map((stat) => (
        <div
          key={stat.label}
          style={{
            padding: "20px",
            background: "#111d26",
            border: "1px solid #29414d",
            borderRadius: "12px"
          }}
        >
          <p
            style={{
              color: "#9fb8c7",
              fontSize: "13px",
              marginTop: 0
            }}
          >
            {stat.label}
          </p>

          <h3 style={{ fontSize: "27px", margin: "8px 0" }}>
            {stat.value}
          </h3>

          <p
            style={{
              color: "#9fb8c7",
              fontSize: "12px",
              marginBottom: 0
            }}
          >
            {stat.detail}
          </p>
        </div>
      ))}
    </div>

    <p
      style={{
        marginTop: "22px",
        color: "#9fb8c7",
        fontSize: "13px"
      }}
    >
      These statistics describe AlphaFold's predicted relative-position
      uncertainty. They are not experimental accuracy measurements.
      Nearby residues are included and can contribute lower PAE values.
    </p>
  </section>
)}

        <section
          id="pae-explanation"
className="paeSection"
          style={{
            marginTop: "24px",
            padding: "28px",
            background: "#0d141b",
            border: "1px solid #263843",
            borderRadius: "16px"
          }}
        >

{paeSelectedPair && paeData && (
  <div>
  
    <p className="sectionLabel">INTERACTIVE RESIDUE ANALYSIS</p>
    <h2>Selected Residue Pair</h2>

    <p style={{ color: "#9fb8c7" }}>
      Analysis for {proteinInfo?.name || activeProteinId}
    </p>

    <h3>
      Residue {paeSelectedPair.y} relative to residue {paeSelectedPair.x}
    </h3>

    <p>
      Predicted aligned error:{" "}
      <strong>
        {Number(paeSelectedPair.error).toFixed(2)} Å
      </strong>
    </p>

    <div
      style={{
        padding: "20px",
        background: "#111d26",
        border: "1px solid #29414d",
        borderRadius: "12px"
      }}
    >
      <h3>What does this prediction suggest?</h3>

      <p>
        {Number(paeSelectedPair.error) < 5
          ? "Low predicted aligned error: AlphaFold predicts the relative positioning of these residues with comparatively high confidence."
          : Number(paeSelectedPair.error) < 15
            ? "Moderate predicted aligned error: AlphaFold indicates some uncertainty in their relative positioning."
            : "High predicted aligned error: AlphaFold indicates substantial uncertainty in their relative positioning."}
      </p>

      <p>
        {Math.abs(paeSelectedPair.x - paeSelectedPair.y) <= 5
          ? "These residues are close together in the amino-acid sequence. Low PAE between nearby residues is common and does not establish confidence in the protein's overall fold."
          : "These residues are separated in the amino-acid sequence. Their PAE provides information about the predicted uncertainty in their relative structural arrangement."}
      </p>
    </div>

{paeStats && Number.isFinite(Number(paeSelectedPair.error)) && (
  <div
    style={{
      marginTop: "20px",
      padding: "20px",
      background: "#111d26",
      border: "1px solid #29414d",
      borderRadius: "12px"
    }}
  >
    <h3>How does this pair compare with the whole protein?</h3>

    <p>
      Selected pair PAE:{" "}
      <strong>{Number(paeSelectedPair.error).toFixed(2)} Å</strong>
    </p>

    <p>
      Whole-protein mean PAE:{" "}
      <strong>{paeStats.mean.toFixed(2)} Å</strong>
    </p>

    <p>
      Difference:{" "}
      <strong>
        {Math.abs(Number(paeSelectedPair.error) - paeStats.mean).toFixed(2)} Å{" "}
        {Number(paeSelectedPair.error) < paeStats.mean
          ? "below"
          : Number(paeSelectedPair.error) > paeStats.mean
            ? "above"
            : "equal to"}{" "}
        the protein average
      </strong>
    </p>

    <p style={{ color: "#b7cbd5", lineHeight: 1.7 }}>
      {Number(paeSelectedPair.error) < paeStats.mean
        ? "This residue pair has lower predicted relative-position uncertainty than the average pair in this protein."
        : Number(paeSelectedPair.error) > paeStats.mean
          ? "This residue pair has higher predicted relative-position uncertainty than the average pair in this protein."
          : "This residue pair has predicted relative-position uncertainty equal to the protein average."}
    </p>

    <p style={{ color: "#9fb8c7", fontSize: "13px" }}>
      This comparison is descriptive. Being below the protein average
      does not necessarily mean the pair has low absolute PAE or that
      its structure has been experimentally validated.
    </p>
  </div>
)}

{paeSelectedPair && paeStats && (
  <div
    style={{
      marginTop: "20px",
      padding: "22px",
      background: "#111d26",
      border: "1px solid #29414d",
      borderRadius: "12px"
    }}
  >
    <h3>Visual PAE Comparison</h3>

    <div
      style={{
        display: "flex",
        justifyContent: "space-between",
        marginBottom: "10px",
        color: "#b7cbd5"
      }}
    >
      <span>0 Å</span>
      <span>32 Å</span>
    </div>

    <div
      style={{
        position: "relative",
        height: "24px",
        background:
          "linear-gradient(to right, #1464d2, #8a997d, #ffcd28)",
        borderRadius: "8px"
      }}
    >
      <div
        style={{
          position: "absolute",
          left: `${Math.min(100, Math.max(0, Number(paeSelectedPair.error) / 32 * 100))}%`,
          top: "-6px",
          width: "4px",
          height: "36px",
          background: "#ffffff",
          borderRadius: "3px",
          transform: "translateX(-50%)"
        }}
      />

      <div
        style={{
          position: "absolute",
          left: `${Math.min(100, paeStats.mean / 32 * 100)}%`,
          top: "-6px",
          width: "4px",
          height: "36px",
          background: "#ff647c",
          borderRadius: "3px",
          transform: "translateX(-50%)"
        }}
      />
    </div>

    <div
      style={{
        display: "flex",
        gap: "24px",
        flexWrap: "wrap",
        marginTop: "20px"
      }}
    >
      <span>
        <strong style={{ color: "#ffffff" }}>●</strong>{" "}
        Selected pair: {Number(paeSelectedPair.error).toFixed(2)} Å
      </span>

      <span>
        <strong style={{ color: "#ff647c" }}>●</strong>{" "}
        Protein average: {paeStats.mean.toFixed(2)} Å
      </span>
    </div>

    <p style={{ color: "#9fb8c7", fontSize: "13px" }}>
      The white marker represents the selected residue pair.
      The pink marker represents the whole-protein average.
      Lower PAE indicates lower predicted relative-position uncertainty.
    </p>
  </div>
)}

    <h3 style={{ marginTop: "24px" }}>Annotated biological regions</h3>

    {[paeSelectedPair.y, paeSelectedPair.x].map((residue, index) => {
      const matchingRegions = (proteinInfo?.regions || []).filter(
        (region) =>
          Number(region.start) <= residue &&
          Number(region.end) >= residue &&
          region.type !== "Chain"
      );

      return (
        <div
          key={index}
          style={{
            padding: "16px 20px",
            marginTop: "12px",
            background: "#111d26",
            border: "1px solid #29414d",
            borderRadius: "12px"
          }}
        >
          <strong>Residue {residue}</strong>

          {matchingRegions.length > 0 ? (
            matchingRegions.map((region, regionIndex) => (
              <p key={regionIndex}>
                {region.name} ({region.type}, residues{" "}
                {region.start}–{region.end})
              </p>
            ))
          ) : (
            <p style={{ color: "#9fb8c7" }}>
              No specific domain or region annotation available
              from UniProt for this residue.
            </p>
          )}
        </div>
      );
    })}

    <p
      style={{
        color: "#9fb8c7",
        fontSize: "13px",
        marginTop: "20px"
      }}
    >
      Biological region information comes from UniProt sequence
      annotations. PAE describes predicted relative-position
      uncertainty; it does not directly measure physical distance,
      molecular interactions, or biological function.
    </p>
  </div>
)}
</section>
<section
    style={{
      marginTop: "24px",
      padding: "28px",
      background: "#0d141b",
      border: "1px solid #263843",
      borderRadius: "16px"
    }}
  >

          <p className="sectionLabel">UNDERSTANDING THE PREDICTION</p>
          <h2>What does Predicted Aligned Error tell us?</h2>

          <p>
            Predicted Aligned Error (PAE) estimates how uncertain AlphaFold
            is about the position of one residue relative to another.
            It is measured in ångströms (Å).
          </p>

          <div
            style={{
              display: "grid",
              gridTemplateColumns: "repeat(auto-fit, minmax(240px, 1fr))",
              gap: "18px",
              marginTop: "24px"
            }}
          >
            {[
              {
                title: "Low PAE — Blue",
                description:
                  "AlphaFold predicts the relative positions of these residues with greater confidence. Blue blocks can indicate regions whose internal arrangement is well constrained."
              },
              {
                title: "High PAE — Yellow",
                description:
                  "AlphaFold is less certain about the relative positioning of these residues. This can occur between flexible regions or domains with uncertain orientations."
              },
              {
                title: "The Diagonal",
                description:
                  "The diagonal compares each residue with itself. Low values near the diagonal are expected and do not establish that the entire protein structure is reliable."
              },
              {
                title: "Off-Diagonal Blocks",
                description:
                  "Blue blocks away from the diagonal can suggest confidently predicted relationships between distant sequence regions. Yellow areas between blocks can suggest uncertain domain positioning."
              }
            ].map((item) => (
              <div
                key={item.title}
                style={{
                  padding: "20px",
                  background: "#111d26",
                  border: "1px solid #29414d",
                  borderRadius: "12px"
                }}
              >
                <h3 style={{ marginTop: 0 }}>{item.title}</h3>
                <p style={{ color: "#b7cbd5", lineHeight: 1.65 }}>
                  {item.description}
                </p>
              </div>
            ))}
          </div>

          <h3 style={{ marginTop: "30px" }}>
            PAE versus pLDDT: what's the difference?
          </h3>

          <p>
            <strong>pLDDT</strong> estimates confidence in the local
            structure around individual residues.
            <strong> PAE</strong> estimates uncertainty in the relative
            positions of pairs of residues.
          </p>

          <p>
            A protein can have highly confident individual domains but
            uncertain orientations between those domains. That's why
            examining both pLDDT and PAE provides a more complete
            understanding of an AlphaFold prediction.
          </p>

          <p style={{ color: "#9fb8c7", fontSize: "13px" }}>
            Scientific limitation: PAE is a prediction of structural
            uncertainty, not a measured physical distance between
            residues or direct experimental evidence of accuracy.
          </p>
        </section>
      </main>
    </div>
  );
}

if (activePage === "geometry") {
  const backboneAngles = geometryResidues.map((residue, index) => {
  const previous = geometryResidues[index - 1];
  const next = geometryResidues[index + 1];

  const current = residue.backbone;

  const phi =
    previous?.chain === residue.chain &&
    previous.number === residue.number - 1 &&
    previous.backbone?.C &&
    current?.N && current?.CA && current?.C
      ? calculateDihedral(
          previous.backbone.C,
          current.N,
          current.CA,
          current.C
        )
      : null;

  const psi =
    next?.chain === residue.chain &&
    next.number === residue.number + 1 &&
    current?.N && current?.CA && current?.C &&
    next.backbone?.N
      ? calculateDihedral(
          current.N,
          current.CA,
          current.C,
          next.backbone.N
        )
      : null;

  return { number: residue.number, chain: residue.chain, phi, psi };
});
    return (
    <div className="app">
      <FoldCheckNav activePage={activePage} setActivePage={setActivePage} />

      <main className="workspacePage">

<div
  style={{
    position: "absolute",
    left: "-10000px",
    top: 0,
    width: "600px",
    height: "500px",
    pointerEvents: "none"
  }}
>
  <ProteinViewer
    source="alphafold"
    proteinId={activeProteinId}
    onGeometryLoaded={setGeometryResidues}
    onProteinConfidence={setProteinConfidence}
  />
</div>

<div className="analysisTabs">
  <button
    className="analysisTab"
    onClick={() => setActivePage("workspace")}
  >
    3D Structure
  </button>

  <button
    className="analysisTab"
    onClick={() => setActivePage("confidence")}
  >
    Confidence
  </button>

  <button
    className="analysisTab"
    onClick={() => setActivePage("pae")}
  >
    PAE
  </button>

  <button
    className="analysisTab active"
    onClick={() => setActivePage("geometry")}
  >
    Geometry
  </button>

  <button
    className="analysisTab"
    onClick={() => setActivePage("experimental")}
  >
    Experimental Comparison
  </button>
</div>

        <section
  className="workspaceTop"
  style={{
  display: "flex",
  flexDirection: "column",
  alignItems: "stretch",
  gap: "24px",
  width: "100%",
  maxWidth: "1100px",
  margin: "0 auto"
}}
>
          <h1>Protein Geometry Analysis</h1>
          <p>
            Explore the structural geometry and backbone
            characteristics of your selected protein.
          </p>
<div className="geometryMetric">
  <h3>Protein Chain Length</h3>
  <p>
    {proteinInfo?.length ?? "Not available"} amino acid residues
  </p>
  <small>
    Chain length provides context for interpreting the protein's
    structural geometry.
  </small>
</div>

<div className="geometryMetric">
  <h3>Average Structure Confidence</h3>
  <p>
  {proteinConfidence?.mean != null
    ? `${proteinConfidence.mean.toFixed(1)} / 100`
    : "Not available"}
</p>
  <small>
    Average pLDDT summarises AlphaFold's confidence
    in the predicted local structure. It is not
    a direct measurement of geometric quality.
  </small>
</div>

<div className="geometryMetric">
  <h3>Backbone Torsion Angles</h3>
  <p>Phi (φ) and Psi (ψ)</p>
  <small>
    Phi and psi describe rotations along the protein backbone.
    Their combinations help determine whether a residue adopts
    an alpha-helix, beta-sheet, or another conformation.
  </small>
  <p>
  {geometryResidues.length
    ? `${geometryResidues.filter((r) =>
        r.backbone?.N && r.backbone?.CA && r.backbone?.C
      ).length} residues have complete backbone coordinates available.`
    : "Waiting for structural coordinates..."}
</p>
</div>

<div className="geometryMetric">
  <h3>Residue Backbone Angles</h3>
  <div style={{ maxHeight: "350px", overflowY: "auto" }}>
    <table>
      <thead>
        <tr>
          <th>Residue</th>
          <th>Phi (φ)</th>
          <th>Psi (ψ)</th>
        </tr>
      </thead>
      <tbody>
        {backboneAngles.map((residue, index) => (
          <tr key={`${residue.chain}-${residue.number}-${index}`}>
            <td>{residue.chain} · {residue.number}</td>
            <td>{residue.phi === null ? "N/A" : `${residue.phi.toFixed(1)}°`}</td>
            <td>{residue.psi === null ? "N/A" : `${residue.psi.toFixed(1)}°`}</td>
          </tr>
        ))}
      </tbody>
    </table>
  </div>
</div>

<div className="geometryMetric">
  <h3>Ramachandran Plot</h3>
  <p>Backbone torsion angles for residues with both φ and ψ available.</p>

  <svg
    viewBox="0 0 400 400"
    style={{
      width: "100%",
      maxWidth: "440px",
      background: "#0b1624",
      borderRadius: "12px"
    }}
  >
    <rect x="40" y="20" width="340" height="340" fill="#101f30" stroke="#64748b" />

    <line x1="210" y1="20" x2="210" y2="360" stroke="#475569" />
    <line x1="40" y1="190" x2="380" y2="190" stroke="#475569" />

    <text x="190" y="390" fill="white">Phi (φ)</text>
    <text x="8" y="195" fill="white" fontSize="12">ψ</text>

    {backboneAngles
      .filter((r) => Number.isFinite(r.phi) && Number.isFinite(r.psi))
      .map((r, index) => (
        <circle
          key={`${r.chain}-${r.number}-${index}`}
          cx={40 + ((r.phi + 180) / 360) * 340}
          cy={360 - ((r.psi + 180) / 360) * 340}
          r="3"
          fill="#38bdf8"
          opacity="0.75"
        />
      ))}
  </svg>

  <small>
    Each point represents one residue. This plot displays calculated
    backbone angles; it does not yet classify favoured or outlier regions.
  </small>
</div>

          
        </section>
        
      </main>
    </div>
  );
}

if (activePage === "experimental") {
  return (
    <div className="app">
      <FoldCheckNav
        activePage={activePage}
        setActivePage={setActivePage}
      />

      <main className="workspacePage">
        <section className="workspaceTop">
          <div>
            <p className="sectionLabel">EXPERIMENTAL EVIDENCE</p>
            <h1>Experimental Structure Comparison</h1>
            <p>
              {proteinInfo?.name || activeProteinId} · {activeProteinId}
            </p>
          </div>
        </section>
<nav className="analysisTabs">
  <button
    className="analysisTab"
    onClick={() => setActivePage("workspace")}
  >
    3D Structure
  </button>

  <button
    className="analysisTab"
    onClick={() => setActivePage("confidence")}
  >
    Confidence
  </button>

  <button
    className="analysisTab"
    onClick={() => setActivePage("pae")}
  >
    PAE
  </button>

  <button
    className="analysisTab"
    onClick={() => setActivePage("geometry")}
  >
    Geometry
  </button>

  <button className="analysisTab active">
    Experimental Comparison
  </button>
</nav>
        <section className="paeSection" style={{ marginTop: "24px" }}>
  <p className="sectionLabel">EXPERIMENTAL STRUCTURE DATABASE</p>
  <h2>Available PDB Structures</h2>

  <p>
    These structures are linked to {proteinInfo?.name || activeProteinId}
    {" "}through UniProt. Experimental entries may cover only part
    of the protein or include other molecules.
  </p>

  {experimentalLoading ? (
    <p>Searching UniProt for experimental structures...</p>
  ) : experimentalError ? (
    <p role="alert">{experimentalError}</p>
  ) : experimentalStructures.length === 0 ? (
    <p>
      No PDB cross-references were found for this protein in UniProt.
    </p>
  ) : (
    <>
      <label htmlFor="experimentalPdbSelect">
        Select an experimental structure
      </label>

      <select
        id="experimentalPdbSelect"
        value={selectedExperimentalId}
        onChange={(event) =>
          setSelectedExperimentalId(event.target.value)
        }
        style={{
          display: "block",
          width: "100%",
          maxWidth: "500px",
          margin: "14px 0 24px",
          padding: "14px",
          background: "#111d26",
          color: "#ffffff",
          border: "1px solid #29414d",
          borderRadius: "10px"
        }}
      >
        {experimentalStructures.map((structure) => (
          <option key={structure.id} value={structure.id}>
            {structure.id} · {structure.method} · {structure.resolution}
          </option>
        ))}
      </select>

      {experimentalStructures
        .filter((structure) => structure.id === selectedExperimentalId)
        .map((structure) => (
          <div key={structure.id}>
<div
  style={{
    display: "grid",
    gridTemplateColumns: "repeat(auto-fit, minmax(220px, 1fr))",
    gap: "16px",
    marginBottom: "28px"
  }}
>
  {[
    {
      label: "Selected Protein",
      value: proteinInfo?.gene || activeProteinId,
      detail: `${proteinInfo?.length ?? "Unknown"} amino acids in UniProt`
    },
    {
      label: "Predicted Structure",
      value: "AlphaFold",
      detail: proteinConfidence?.mean != null
        ? `Mean pLDDT: ${proteinConfidence.mean.toFixed(1)} / 100`
        : "Local confidence score unavailable"
    },
    {
      label: "Experimental Structure",
      value: structure.id,
      detail: structure.method
    },
    {
      label: "Experimental Resolution",
      value: structure.resolution,
      detail: "Reported by UniProt, where available"
    }
  ].map((item) => (
    <div
      key={item.label}
      style={{
        padding: "22px",
        background: "#111d26",
        border: "1px solid #29414d",
        borderRadius: "12px"
      }}
    >
      <p style={{ color: "#9fb8c7", marginTop: 0 }}>
        {item.label}
      </p>
      <h3 style={{ fontSize: "24px", margin: "10px 0" }}>
        {item.value}
      </h3>
      <p style={{ color: "#b7cbd5", fontSize: "13px" }}>
        {item.detail}
      </p>
    </div>
  ))}
</div>

            <p className="foldcheckLinkedHint">Select a residue on either 3D model to highlight the corresponding residue in both structures, when a matching residue is present. The experimental structure may cover only part of the protein.</p>

            <div
              style={{
                display: "grid",
                gridTemplateColumns: "repeat(auto-fit, minmax(min(100%, 420px), 1fr))",
                gap: "24px",
                alignItems: "start",
                marginTop: "24px"
              }}
            >
              <section style={{ minWidth: 0 }}>
                <h3>AlphaFold Prediction</h3>
                <div className="viewerShell" style={{ height: "500px", overflow: "hidden", borderRadius: "12px", marginTop: "16px" }}>
                  <ProteinViewer
                    source="alphafold"
                    proteinId={activeProteinId}
                    onGeometryLoaded={setGeometryResidues}
                    linkedResidue={linkedSelection?.alphafold}
                    onResidueSelect={(residue) => handleLinkedSelection("alphafold", residue)}
                  />
                </div>
                <p><strong>Protein:</strong> {proteinInfo?.name || activeProteinId}</p>
                <p><strong>Sequence length:</strong> {proteinInfo?.length ?? "Unknown"} amino acids</p>
                <p><strong>Mean pLDDT:</strong> {proteinConfidence?.mean != null ? proteinConfidence.mean.toFixed(1) : "Unavailable"}</p>
              </section>
              <section style={{ minWidth: 0 }}>
                <h3>Experimental PDB {structure.id}</h3>
                <div className="viewerShell" style={{ height: "500px", overflow: "hidden", borderRadius: "12px", marginTop: "16px" }}>
                  <ProteinViewer
                    source="pdb"
                    proteinId={structure.id}
                    onGeometryLoaded={setExperimentalResidues}
                    linkedResidue={linkedSelection?.pdb}
                    onResidueSelect={(residue) => handleLinkedSelection("pdb", residue)}
                  />
                </div>
                <p><strong>Method:</strong> {structure.method}</p>
                <p><strong>Resolution:</strong> {structure.resolution}</p>
                <p><strong>Chains / coverage:</strong> {structure.chains}</p>
              </section>
            </div>

            <p style={{ color: "#9fb8c7", fontSize: "13px" }}>
              This viewer displays the experimentally determined
              PDB entry. The two structures are not yet
              superimposed, so their visual orientations cannot
              be used to measure structural agreement.
            </p>

            <a
              href={`https://www.rcsb.org/structure/${encodeURIComponent(structure.id)}`}
              target="_blank"
              rel="noopener noreferrer"
            >
              View PDB {structure.id} at RCSB ↗
            </a>
            <div
  style={{
    marginTop: "28px",
    padding: "28px",
    background: "#111d26",
    border: "1px solid #29414d",
    borderRadius: "14px"
  }}
>
  <p className="sectionLabel">STRUCTURAL COMPARISON ANALYSIS</p>
<div
  style={{
    padding: "22px",
    marginBottom: "24px",
    background: "#10242b",
    border: "1px solid #28535a",
    borderRadius: "12px",
  }}
>
  <p className="sectionLabel">QUANTITATIVE STRUCTURAL COMPARISON</p>

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
    <strong>Selected PDB chain: </strong>
    {selectedExperimentalChain ?? "Unavailable"}
  </div>

  <div>
    <strong>Matched residues: </strong>
    {matchedResidues.length}
  </div>

  <div>
    <strong>Experimental coverage: </strong>
    {experimentalSequence.length > 0
      ? `${(
          (matchedResidues.length / experimentalSequence.length) *
          100
        ).toFixed(1)}%`
      : "Unavailable"}
  </div>
</div>

  {alignedRMSD !== null && Number.isFinite(alignedRMSD) ? (
    <>
      <h2 style={{ fontSize: "34px", margin: "12px 0" }}>
        {alignedRMSD.toFixed(2)} Å
      </h2>

      <p>
        Calculated after optimal rigid-body superposition using{" "}
        <strong>{rmsdCoordinates.length}</strong> matched Cα atom pairs.
      </p>

      <p style={{ color: "#9fb8c7", fontSize: "13px" }}>
        Preliminary result. Sequence correspondence, experimental chain
        selection, and numerical validation must be checked before
        interpreting this value as structural agreement.
      </p>
    </>
  ) : (
    <p>
      RMSD unavailable. At least three matched residues with valid
      Cα coordinates are required.
    </p>
  )}
</div>
  <h2>Interpreting the Evidence</h2>

  <div
    style={{
      display: "grid",
      gridTemplateColumns: "repeat(auto-fit, minmax(240px, 1fr))",
      gap: "20px",
      marginTop: "20px"
    }}
  >
    <div>
      <h3>AlphaFold Prediction</h3>
      <p>
        AlphaFold predicts the protein's three-dimensional structure
        from its amino acid sequence.
      </p>
      <p>
        {proteinConfidence?.mean != null
          ? `The loaded model has a mean pLDDT of ${proteinConfidence.mean.toFixed(1)} / 100. This describes predicted local confidence, not experimentally measured accuracy.`
          : "Local confidence information is not currently available."}
      </p>
    </div>

    <div>
      <h3>Experimental Structure</h3>
      <p>
        PDB {structure.id} is linked to this protein through UniProt.
        Its reported determination method is {structure.method}.
      </p>
      <p>
        Reported resolution: {structure.resolution}.
        The experimental entry may represent only a portion of
        the protein or contain additional molecular components.
      </p>
    </div>
  </div>

  <h3 style={{ marginTop: "24px" }}>
    What can we conclude?
  </h3>

  <p>
    Comparing predicted and experimentally determined structures
    can help investigate whether a prediction agrees with observed
    molecular conformations. However, this page currently displays
    the structures separately.
  </p>

  <p>
    A quantitative structural comparison requires matching
    corresponding residues, aligning their atomic coordinates,
    and calculating a measure such as backbone RMSD.
    Different experimental conditions, missing residues,
    and molecular interactions can also affect the comparison.
  </p>

  <p style={{ color: "#9fb8c7", fontSize: "13px" }}>
    Scientific limitation: A high AlphaFold confidence score
    does not guarantee agreement with a particular experimental
    structure. An aligned Cα RMSD is calculated when enough matching residues are available; the two 3D viewers are not visually superimposed.
  </p>
</div>
          </div>
        ))}
    </>
  )}
</section>
      </main>
    </div>
  );
}

if (activePage === "basics") {
  const topics = [
    ["Amino acids", "Proteins are built from amino acids linked by peptide bonds. Side chains determine many of their chemical properties."],
    ["Primary structure", "The amino-acid sequence, read from the N-terminus to the C-terminus."],
    ["Secondary structure", "Alpha helices and beta sheets are local backbone patterns stabilised by hydrogen bonds."],
    ["Tertiary structure", "The 3D arrangement of a polypeptide, influenced by hydrophobic packing, ionic interactions, hydrogen bonds and sometimes disulfide bonds."],
    ["Quaternary structure", "The assembly of multiple protein chains into a complex. Not every protein forms one."],
    ["Domains", "A domain is a region of a protein that can form a distinct structural or functional unit."],
    ["Organelles and localisation", "A protein's cellular location helps explain its role. Some proteins work in the nucleus, others in mitochondria, the endoplasmic reticulum, membranes or cytoplasm. Localisation can change with conditions."],
    ["AlphaFold model", "An AlphaFold structure is a computational prediction, not an experimental measurement or a complete simulation of a protein inside a cell."],
    ["pLDDT", "Predicted local distance difference test estimates local model confidence. Above 90: very high; 70–90: confident; 50–70: low; below 50: very low."],
    ["PAE", "Predicted aligned error describes uncertainty in relative positions of residues or domains, complementing pLDDT."],
    ["Experimental structures", "PDB entries can be determined using X-ray crystallography, cryo-electron microscopy or NMR. They may cover only part of a protein."],
    ["Mutations", "A missense variant substitutes one amino acid. Effects depend on the position, molecular environment and biological context; sequence identity alone cannot establish pathogenicity."]
  ];
  return (
    <>
      <FoldCheckNav activePage={activePage} setActivePage={setActivePage} />
      <main className="workspacePage" style={{ maxWidth: 1200, margin: "0 auto", padding: "40px 24px" }}>
        <p className="sectionLabel">LEARN THE SCIENCE</p>
        <h1>Protein Basics</h1>
        <p>Explore protein architecture, cellular localisation and the evidence behind structure predictions.</p>
        <section className="foldcheckLearningGrid">
          {topics.map(([title, description]) => (
            <article className="foldcheckLearningCard" key={title}>
              <h2>{title}</h2>
              <p>{description}</p>
            </article>
          ))}
        </section>
        <section className="foldcheckDetailPanel">
          <p className="sectionLabel">CURRENTLY LOADED PROTEIN</p>
          <h2>{proteinInfo?.name || activeProteinId}</h2>
          <p><strong>Organism:</strong> {proteinInfo?.organism || "Loading UniProt information"}</p>
          <p><strong>Gene:</strong> {proteinInfo?.gene || "Not specified"}</p>
          <p><strong>Length:</strong> {proteinInfo?.length ? `${proteinInfo.length} amino acids` : "Not available"}</p>
          <h3>Where is it found in the cell?</h3>
          {proteinInfo?.subcellularLocations?.length ? (
            <div className="foldcheckTagList">
              {[...new Set(proteinInfo.subcellularLocations.map((item) => item.name))].map((name) => (
                <span className="foldcheckTag" key={name}>{name}</span>
              ))}
            </div>
          ) : (
            <p>UniProt has no subcellular-location annotation available for this entry, or the annotation is still loading.</p>
          )}
          <h3>Annotated domains and regions</h3>
          {proteinInfo?.regions?.length ? (
            <div className="foldcheckRegionList">
              {proteinInfo.regions.slice(0, 18).map((region, index) => (
                <button key={`${region.name}-${index}`} type="button" onClick={() => {
                  setSelectedRegion(region);
                  setActivePage("workspace");
                }}>
                  <strong>{region.name}</strong>
                  <span>{region.start}–{region.end} · {region.type}</span>
                </button>
              ))}
            </div>
          ) : <p>No region annotations available for this entry.</p>}
          <h3>Function</h3>
          <p>{proteinInfo?.function || "Function annotation unavailable for this entry."}</p>
          <button className="navButton" onClick={() => setActivePage("workspace")}>Inspect the actual 3D model</button>
          <p className="foldcheckSourceNote">Protein descriptions, localisation and regions are retrieved from UniProt. The 3D model is obtained separately from AlphaFold.</p>
        </section>
      </main>
    </>
  );
}

if (activePage === "workspace") {
  return(
      <div className="app">
        <FoldCheckNav activePage={activePage} setActivePage={setActivePage} />

        <main className="workspacePage">
          <section className="workspaceTop">
            <div>
              <p className="sectionLabel">PROTEIN WORKSPACE</p>
              <h1>{proteinInfo?.name || activeProteinId}</h1>
<p>
  {proteinInfo
    ? [
  proteinInfo.gene,
  proteinInfo.organismCommonName
    ? `${proteinInfo.organism} — ${proteinInfo.organismCommonName}`
    : proteinInfo.organism
]
        .filter(Boolean)
        .join(" · ")
    : proteinInfoLoading
      ? "Loading protein information..."
      : proteinInfoError || "Protein information unavailable"}
</p>
            </div>

            <div className="workspaceAccession">
              <span>UNIPROT</span>
              <strong>{activeProteinId}</strong>
            </div>
          </section>

          <nav className="analysisTabs">
            <button
  className="analysisTab active"
  onClick={() => setActivePage("workspace")}
>
  3D Structure
</button>
            <button
  className="analysisTab"
  onClick={() => setActivePage("confidence")}
>
  Confidence
</button>
            <button
  className="analysisTab"
  onClick={() => setActivePage("pae")}
>
  PAE
</button>
            <button className="analysisTab" onClick={() => setActivePage("geometry")}>Geometry</button>
            <button className="analysisTab" onClick={() => setActivePage("experimental")}>Experimental Comparison</button>
          </nav>

          <section className="proteinOverview">
  <p className="sectionLabel">PROTEIN OVERVIEW</p>
  <h2>What does this protein do?</h2>

  <p>
    <strong>Sequence length: </strong>
    {proteinInfo?.length
      ? `${proteinInfo.length} amino acids`
      : "Unavailable"}
  </p>

  <h3>Biological function</h3>
  <p>
    {proteinInfoLoading ? (
  <p>Loading protein information...</p>
) : proteinInfo?.function ? (
<details className="biologicalFunctionDetails">
  <summary>
    <span className="functionPreviewText">
      {proteinInfo.function.length > 280
        ? `${proteinInfo.function.slice(0, 280).trimEnd()}...`
        : proteinInfo.function}
    </span>

    {proteinInfo.function.length > 280 && (
      <span className="functionExpandHint">
        <span className="showMoreText">More details ↓</span>
        <span className="showLessText">Show less ↑</span>
      </span>
    )}
  </summary>

  {proteinInfo.function.length > 280 && (
    <ul className="biologicalFunctionFull">
  {proteinInfo.function
    .split(/(?<=[.!?])\s+/)
    .filter((sentence) => sentence.trim().length > 0)
    .map((sentence, index) => (
      <li key={index}>{sentence.trim()}</li>
    ))}
</ul>
  )}
</details>
) : (
  <p>No functional description available from UniProt.</p>
)}
  </p>

  <h3 className="functionalRegionsHeading">
  Annotated functional regions
</h3>
<p className="foldcheckRegionHint">Click any residue-range button below to highlight that region in the 3D protein viewer.</p>

<details className="functionalRegionsDetails">
  <summary className="functionalRegionsToggle">
  <span className="regionsShowText">
    Show all {proteinInfo?.regions?.length ?? 0} regions ↓
  </span>
  <span className="regionsHideText">
    Show fewer regions ↑
  </span>
</summary>


  {proteinInfo?.regions?.length ? (
    proteinInfo.regions.map((region, index) => (
      <div className="functionalRegionRow" key={index}>
  <span className="functionalRegionName">
    {region.name}
  </span>

  <button
  type="button"
  className="functionalRegionRange"
  disabled={region.start == null || region.end == null}
 onClick={() => {
  setSelectedRegion(region);
  setSelectedResidue(null);
  setStructureQuery("");
  setSearchFeedback(null);
  setInspectorTab("region");

  document.querySelector(".structureWorkspace")?.scrollIntoView({
    behavior: "smooth",
    block: "start",
  });
}}
  title="Explore this region in 3D"
>
  <span className="functionalRegionAction">CLICK TO HIGHLIGHT</span>
  Residues {region.start ?? "?"}–{region.end ?? "?"}
  <span aria-hidden="true"> ↗</span>
</button>
</div>
    ))
  ) : (
    <p>No domain or region annotations available.</p>
  )}
  </details>
</section>

          

          <section
            className="structureWorkspace"
            id="structure"
            ref={workspaceRef}
            style={{ "--panel-width": `${panelWidth}px` }}
          >
            <div className="structureViewer">
              <div className="viewerHeading">
                <div>
                  <span className="viewerEyebrow">
                    INTERACTIVE STRUCTURE
                  </span>
                  <h2>Explore the protein</h2>
                </div>

                <button className="viewerHelp">
                  ? How do I use this?
                </button>
              </div>

              <div className="structureSearch">
  <div className="structureSearchHeading">
    <strong>Find a residue or atom</strong>
    <span>Search within the loaded protein structure</span>
  </div>

  <form
    className="structureSearchForm"
    onSubmit={(event) => {
  event.preventDefault();
  setSelectedRegion(null);
  findStructureTarget();
}}
  >
    <div className="structureSearchInputWrap">
  <input
    type="text"
    value={structureQuery}
    onChange={(event) => {
      setStructureQuery(event.target.value);
      setShowStructureSuggestions(true);
    }}
    onFocus={() => setShowStructureSuggestions(true)}
    onKeyDown={(event) => {
  if (event.key === "Escape") {
    setShowStructureSuggestions(false);
    setActiveSuggestionIndex(-1);
  }

  if (showStructureSuggestions && structureSuggestions.length > 0) {
    if (event.key === "ArrowDown") {
      event.preventDefault();
      setActiveSuggestionIndex((index) =>
        Math.min(index + 1, structureSuggestions.length - 1)
      );
    }

    if (event.key === "ArrowUp") {
      event.preventDefault();
      setActiveSuggestionIndex((index) => Math.max(index - 1, 0));
    }

    if (event.key === "Enter" && activeSuggestionIndex >= 0) {
      event.preventDefault();
      const [code] = structureSuggestions[activeSuggestionIndex];
      setStructureQuery(
        code.includes(" ") ? code : `${code} `
      );
      setShowStructureSuggestions(false);
      setActiveSuggestionIndex(-1);
    }
  }
}}
    placeholder="e.g. ARG 158 CD"
    aria-label="Search residue or atom"
    autoComplete="off"
  />

  {showStructureSuggestions && structureSuggestions.length > 0 && (
    <div className="structureSuggestions">
      {structureSuggestions.map(([code, name], index) => (
        <button
          key={code}
          type="button"
          className={`structureSuggestion ${
  activeSuggestionIndex === index ? "isActive" : ""
}`}
          onClick={() => {
            setStructureQuery(
              code.includes(" ") ? code : `${code} `
            );
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

    <button type="submit">Find →</button>
    <button
  type="button"
  className="structureResetButton"
  onClick={() => {
    setResetViewRequest((value) => value + 1);
    setSelectedRegion(null);
    setSelectedResidue(null);
    setStructureQuery("");
setShowStructureSuggestions(false);
    setInspectorTab("residue");
    setSearchFeedback(null);
  }}
>
  Reset View
</button>
  </form>

  {searchFeedback && (
    <p
      className={
        searchFeedback.success
          ? "structureSearchSuccess"
          : "structureSearchError"
      }
    >
      {searchFeedback.message}
    </p>
  )}
</div>

              <div className="viewerShell">
                <ProteinViewer
                  source="alphafold"
                  proteinId={activeProteinId}
                  onResidueSelect={setSelectedResidue}
                  onProteinConfidence={setProteinConfidence}
                  searchRequest={searchRequest}
onSearchResult={handleStructureSearchResult}
resetViewRequest={resetViewRequest}
onSequenceLoaded={setProteinResidues}
onGeometryLoaded={setGeometryResidues}
selectedRegion={selectedRegion}
                />
              </div>

              <div className="confidenceLegend" id="confidence">
                <div className="confidenceLegendTitle">
                  <span>AlphaFold confidence (pLDDT)</span>

                  <button
                    className="confidenceInfoButton"
                    onClick={() =>
                      setShowConfidenceInfo(!showConfidenceInfo)
                    }
                    aria-label="About AlphaFold confidence"
                  >
                    i
                  </button>
                </div>

                <div className="plddtGradient" />

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

              {showConfidenceInfo && (
                <div className="confidenceInfo">
                  <strong>What does this mean?</strong>

                  <p>
                    These colours show AlphaFold's confidence in the
                    predicted local structure around each residue.
                    Higher confidence does not guarantee biological
                    correctness, and lower confidence does not
                    automatically mean the prediction is wrong.
                  </p>

                  <button className="confidenceLearnMore">
                    Learn more about confidence and how it is
                    determined in The Basics →
                  </button>
                </div>
              )}

              <div className="proteinConfidenceSummary">
                <div className="proteinConfidenceHeader">
                  <div>
                    <span className="guideLabel">
                      WHOLE-PROTEIN CONFIDENCE
                    </span>
                    <h3>Mean predicted local confidence</h3>
                  </div>

                  <strong className="proteinConfidenceValue">
                    {proteinConfidence
                      ? proteinConfidence.mean.toFixed(1)
                      : "—"}
                    <small> / 100</small>
                  </strong>
                </div>

                {proteinConfidence ? (
                  <>
                    <p>
                      Calculated from{" "}
                      {proteinConfidence.residueCount} residues in
                      the loaded AlphaFold model.
                    </p>

                    <div className="proteinConfidenceDistribution">
                      <div>
                        <span>Very high (≥90)</span>
                        <strong>
                          {(
                            (proteinConfidence.veryHigh /
                              proteinConfidence.residueCount) *
                            100
                          ).toFixed(1)}
                          %
                        </strong>
                      </div>

                      <div>
                        <span>Confident (70–90)</span>
                        <strong>
                          {(
                            (proteinConfidence.confident /
                              proteinConfidence.residueCount) *
                            100
                          ).toFixed(1)}
                          %
                        </strong>
                      </div>

                      <div>
                        <span>Low (50–70)</span>
                        <strong>
                          {(
                            (proteinConfidence.low /
                              proteinConfidence.residueCount) *
                            100
                          ).toFixed(1)}
                          %
                        </strong>
                      </div>

                      <div>
                        <span>Very low (&lt;50)</span>
                        <strong>
                          {(
                            (proteinConfidence.veryLow /
                              proteinConfidence.residueCount) *
                            100
                          ).toFixed(1)}
                          %
                        </strong>
                      </div>
                    </div>

                    <div className="proteinConfidenceInterpretation">
                      <h4>What does this tell us about this protein?</h4>

                      <p>
                        The loaded AlphaFold model contains a
                        mixture of confidently predicted and
                        uncertain regions. Its mean pLDDT is{" "}
                        <strong>
                          {proteinConfidence.mean.toFixed(1)}
                        </strong>
                        , but{" "}
                        <strong>
                          {(
                            (proteinConfidence.veryLow /
                              proteinConfidence.residueCount) *
                            100
                          ).toFixed(1)}
                          %
                        </strong>{" "}
                        of its residues have very low confidence,
                        while{" "}
                        <strong>
                          {(
                            (proteinConfidence.veryHigh /
                              proteinConfidence.residueCount) *
                            100
                          ).toFixed(1)}
                          %
                        </strong>{" "}
                        have very high confidence.
                      </p>

                      <h4>Why might confidence vary?</h4>

                     <p>
  Confidence can vary between regions because of differences
  in structural constraints, flexibility, or limitations in
  the prediction. Low pLDDT may be consistent with intrinsic
  disorder, but it does not establish the cause of uncertainty.
</p>

                      <h4>Why does it matter?</h4>

                      <p>
                        The mean score does not represent every
                        part of the protein equally. Its highly confident
                        regions can be investigated with greater
                        structural confidence, while low-confidence
                        regions require more cautious interpretation
                        and additional evidence.
                      </p>
                    </div>
                  </>
                ) : (
                  <p>
                    Waiting for confidence data from the AlphaFold
                    model.
                  </p>
                )}
              </div>
            </div>

            <div
              className="workspaceResizeHandle"
              role="separator"
              aria-label="Resize Selection Details panel"
              aria-orientation="vertical"
              tabIndex={0}
              onPointerDown={startPanelResize}
              onPointerMove={movePanelResize}
              onPointerUp={stopPanelResize}
              onPointerCancel={stopPanelResize}
              onLostPointerCapture={stopPanelResize}
              onKeyDown={(event) => {
                if (event.key === "ArrowLeft") {
                  event.preventDefault();
                  setPanelWidth((width) =>
                    Math.min(width + 20, 800)
                  );
                }

                if (event.key === "ArrowRight") {
                  event.preventDefault();
                  setPanelWidth((width) =>
                    Math.max(width - 20, 280)
                  );
                }
              }}
            />

            <aside className="evidencePanel">
              <div className="evidencePanelHeader">
                <span>SELECTION DETAILS</span>
                <div className="liveDot" />
              </div>

              <div className="inspectorTabs">
  <button
    type="button"
    className={`inspectorTab ${inspectorTab === "residue" ? "active" : ""}`}
    onClick={() => setInspectorTab("residue")}
  >
    Residue Explorer
  </button>

  <button
    type="button"
    className={`inspectorTab ${inspectorTab === "atom" ? "active" : ""}`}
    onClick={() => setInspectorTab("atom")}
  >
    Atom Explorer
  </button>
  <button
  type="button"
  className={`inspectorTab ${inspectorTab === "region" ? "active" : ""}`}
  onClick={() => setInspectorTab("region")}
>
  Region Explorer
</button>
</div>

              {selectedResidue || (inspectorTab === "region" && selectedRegion) ? (
  <>
    {inspectorTab === "residue" && (
      <div className="residueInspector">
                  <span className="guideLabel">
                    SELECTED RESIDUE
                  </span>

                  <h3>
                    {selectedAminoAcid?.name ||
                      selectedResidue.residueCode}
                  </h3>

                  <p className="residueSubtitle">
                    {selectedResidue.residueCode} · Residue{" "}
                    {selectedResidue.residueNumber} · Chain{" "}
                    {selectedResidue.chain}
                  </p>

                  {activeProteinId === "P04637" && p53Region && (
                    <div className="inspectorSection">
                      <h4>Where is this residue in p53?</h4>

                      <strong className="residueRegionName">
                        {p53Region.name}
                      </strong>

                      <p>{p53Region.description}</p>

                      <p className="inspectorNote">
                        Functional region assignments are
                        approximate and based on p53 sequence
                        position, not inferred from the 3D
                        coordinates.
                      </p>
                    </div>
                  )}

                  <div className="inspectorSection">
                    <h4>AlphaFold prediction confidence</h4>

                    {selectedPlddt !== null ? (
                      <>
                        <div className="localConfidenceScale">
                          <div className="localConfidenceReadout">
                            <span>
                              SELECTED RESIDUE · pLDDT
                            </span>

                            <strong>
                              {selectedPlddt.toFixed(1)} / 100
                            </strong>
                          </div>

                          <div className="localConfidenceTrack">
                            <div
                              className="localConfidenceMarker"
                              style={{
                                left: `${selectedPlddt}%`
                              }}
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
                            Classification:{" "}
                            <strong>{confidenceCategory}</strong>
                          </p>
                        </div>

                        <p className="inspectorNote">
                          pLDDT estimates AlphaFold's confidence
                          in the predicted local structure around
                          this residue. It is not a direct
                          measurement of experimental accuracy.
                        </p>

                        <div className="confidenceInterpretation">
                          {confidenceDifference !== null && (
                            <div className="residueConfidenceComparison">
                              <span>
                                Compared with protein's mean confidence
                              </span>

                              <strong>
                                {confidenceDifference >= 0
                                  ? "+"
                                  : ""}
                                {confidenceDifference.toFixed(1)}{" "}
                                pLDDT points
                              </strong>

                              <p>
                                This residue has{" "}
                                {confidenceDifference >= 0
                                  ? "higher"
                                  : "lower"}{" "}
                                local prediction confidence than
                                the average residue in the loaded
                                protein model.
                              </p>
                            </div>
                          )}

                          <h4>
                            Why might confidence be{" "}
                            {confidenceCategory.toLowerCase()}?
                          </h4>

                          <p>{confidenceExplanation}</p>

                          <h4>Why does it matter?</h4>
                          <p>{confidenceSignificance}</p>

                          <h4>
                            What should you investigate next?
                          </h4>

                          <p>
                            Compare neighbouring residues,
                            examine local structural features,
                            and check available PAE or
                            experimental evidence before drawing
                            conclusions about this region.
                          </p>
                        </div>
                      </>
                    ) : (
                      <p>
                        Confidence information is unavailable
                        for this selection.
                      </p>
                    )}
                  </div>

                  {selectedAminoAcid ? (
                    <>
                      <div className="inspectorSection">
                        <h4>About this amino acid</h4>
                        <p>
                          {selectedAminoAcid.description}
                        </p>

                        <div className="inspectorData">
                          <div>
                            <span>One-letter code</span>
                            <strong>
                              {selectedAminoAcid.letter}
                            </strong>
                          </div>

                          <div>
                            <span>Chemical class</span>
                            <strong>
                              {selectedAminoAcid.category}
                            </strong>
                          </div>

                          <div>
                            <span>Side chain</span>
                            <strong>
                              {selectedAminoAcid.sideChain}
                            </strong>
                          </div>
                        </div>
                      </div>

                      <div className="inspectorSection">
                        <h4>Chemical composition</h4>

                        <div className="inspectorData">
                          <div>
                            <span>Free amino acid formula</span>
                            <strong>
                              {selectedAminoAcid.formula}
                            </strong>
                          </div>

                          <div>
                            <span>Free amino acid molar mass</span>
                            <strong>
                              {selectedAminoAcid.mass.toFixed(2)}{" "}
                              g/mol
                            </strong>
                          </div>
                        </div>

                        <p className="inspectorNote">
                          These values describe the free amino
                          acid. Its composition differs when
                          incorporated into a protein chain.
                        </p>
                      </div>
                    </>
                  ) : (
                    <p>
                      Reference information is not available
                      for this residue type.
                    </p>
                  )}

                 
                </div>
)}

{inspectorTab === "region" && !selectedRegion && (
  <div className="residueInspector">
    <span className="guideLabel">REGION EXPLORER</span>
    <h3>No region selected</h3>
    <p className="residueSubtitle">
      Select a functional region from the protein information
      section to explore its biological role and highlight it
      in the 3D structure.
    </p>
  </div>
)}
{inspectorTab === "region" && selectedRegion && (
  <div className="residueInspector">
    <span className="guideLabel">SELECTED FUNCTIONAL REGION</span>

    <h3>{selectedRegion.name}</h3>

<p className="inspectorNote">
  UniProt feature type: <strong>{selectedRegion.type}</strong>
</p>

    <p className="residueSubtitle">
      Residues {selectedRegion.start}–{selectedRegion.end}
      {" · "}
      {selectedRegion.end - selectedRegion.start + 1} amino acids
    </p>

   <div className="inspectorSection">
  <h4>
  {!selectedRegion.name.toLowerCase().includes("interaction") && (
  <h4>
    {selectedRegion.type === "Domain"
      ? `What does the ${selectedRegion.name} domain do?`
      : `What does ${selectedRegion.name} do?`}
  </h4>
)}
</h4>

  

  {selectedRegion.type === "Domain" ? (
    <p>
  {proteinInfo?.function
    ?.split(/(?<=[.!?])\s+/)
    .filter((sentence) => {
      const terms = selectedRegion.name
        .toLowerCase()
        .replace(/\b(domain|region|motif)\b/g, "")
        .trim()
        .split(/\s+/)
        .filter((term) => term.length >= 4);

      return terms.some((term) =>
        sentence.toLowerCase().includes(term)
      );
    })
    .slice(0, 2)
    .join(" ") ||
    `No domain-specific functional explanation was found in the retrieved UniProt annotation for ${selectedRegion.name}.`}
</p>
  ) : selectedRegion.type === "Region" ? (
    <>
      

      {selectedRegion.name.toLowerCase().includes("interaction") ? (
  <></>
) : selectedRegion.name.toLowerCase().includes("disordered") ? (
        <p>
          This region is annotated as disordered, meaning
          it may not maintain one stable three-dimensional
          structure. Such regions can remain flexible and
          participate in important biological interactions.
        </p>
      ) : (
        <p>
  {proteinInfo?.function
    ?.split(/(?<=[.!?])\s+/)
    .filter((sentence) => {
      const terms = selectedRegion.name
        .toLowerCase()
        .replace(/\b(region|interaction|with|binding)\b/g, "")
        .trim()
        .split(/\s+/)
        .filter((term) => term.length >= 4);

      return terms.some((term) =>
        sentence.toLowerCase().includes(term)
      );
    })
    .slice(0, 2)
    .join(" ") ||
    `No specific biological explanation is available for ${selectedRegion.name} in the retrieved protein function annotation.`}
</p>
      )}
    </>
  ) : selectedRegion.type === "Signal" ? (
    <p>
      This signal peptide is a sequence that helps direct
      the protein into the cellular pathway responsible
      for secretion or delivery to certain compartments.
    </p>
  ) : selectedRegion.type === "Chain" ? (
    <p>
      This annotation identifies a protein chain, which
      may represent the mature protein after processing.
    </p>
  ) : (
    <p>
      UniProt classifies this feature as{" "}
      <strong>{selectedRegion.type || "a functional region"}</strong>.
      Its annotation describes a sequence feature of the protein.
    </p>
  )}

  <RegionDetails
  region={selectedRegion}
  proteinInfo={proteinInfo}
  partnerInfo={partnerInfo}
/>



{selectedRegion.description &&
  selectedRegion.description !== selectedRegion.name && (
    <div className="inspectorNote">
      <strong>Region-specific function</strong>
      <p>{selectedRegion.description}</p>
    </div>
  )}

  {selectedRegion.evidence?.length > 0 && (
  <h4>What evidence supports this annotation?</h4>
)}

  {selectedRegion.evidence?.length ? (
  <div className="inspectorData">
    {selectedRegion.evidence.map((evidence, index) => (
      <div key={index}>
        <span>Evidence {index + 1}</span>

        <strong>
          {evidence.source || "UniProt evidence"}
          {evidence.id ? ` · ${evidence.id}` : ""}
        </strong>

        <p className="inspectorNote">
          Evidence: {evidence.evidenceCode === "ECO:0000269"
  ? "Experimental evidence reported in a publication"
  : evidence.evidenceCode === "ECO:0000305"
    ? "Curator inference based on available information"
    : evidence.evidenceCode === "ECO:0000250"
      ? "Inferred from similarity to another protein"
      : evidence.evidenceCode || "Not specified"}
        </p>

        {evidence.source === "PubMed" && evidence.id && (
          <a
            href={`https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(evidence.id)}/`}
            target="_blank"
            rel="noopener noreferrer"
          >
            View supporting publication ↗
          </a>
        )}
      </div>
    ))}
  </div>
) : null}

  <p className="inspectorNote">
    UniProt sequence annotations describe biological features.
    They do not, by themselves, verify molecular binding or
    the accuracy of the displayed AlphaFold coordinates.
  </p>
</div>

    <div className="inspectorSection">
      <h4>Region information</h4>

      <div className="inspectorData">
        <div>
          <span>Start residue</span>
          <strong>{selectedRegion.start}</strong>
        </div>

        <div>
          <span>End residue</span>
          <strong>{selectedRegion.end}</strong>
        </div>

        <div>
          <span>Region length</span>
          <strong>
            {selectedRegion.end - selectedRegion.start + 1} aa
          </strong>
        </div>
      </div>
    </div>

    
  </div>
)}

{inspectorTab === "atom" && (
  <div className="residueInspector">
    <span className="guideLabel">SELECTED ATOM</span>

    <h3>{selectedAtomName || "No atom selected"}</h3>

    <p className="residueSubtitle">
      {selectedAminoAcid?.name || selectedResidue.residueCode}
      {" · "}Residue {selectedResidue.residueNumber}
      {" · "}Chain {selectedResidue.chain}
    </p>

    <div className="inspectorSection">
      <h4>Atom Explorer</h4>

      {selectedAtomName ? (
        <>
          <div className="inspectorData">
            <div>
              <span>Atom identifier</span>
              <strong>{selectedAtomName}</strong>
            </div>
            <div>
              <span>Chemical element</span>
              <strong>{atomElement || "Unknown"}</strong>
            </div>
            <div>
              <span>Structural location</span>
              <strong>{atomLocation}</strong>
            </div>
          </div>

          <h4>What does this atom do?</h4>
          <p>{specificAtomInfo?.role || atomRole}</p>

{specificAtomInfo && (
  <>
    <h4>Normally bonded to</h4>
    <p>{specificAtomInfo.bonds.join(" · ")}</p>

    <h4>Why is this atom important?</h4>
    <p>{specificAtomInfo.significance}</p>
  </>
)}

          
        </>
      ) : (
        <p>Select an atom in the 3D viewer.</p>
      )}
    </div>
  </div>
)}
</>
) : (
                <div className="emptySelection">
                  

                  <h3>Select part of the protein</h3>

                  <p>
                    Click a residue or region in the 3D
                    structure to investigate it in FoldCheck.
                  </p>
                </div>
              )}
            </aside>
</section>



</main>
      </div>
    );
  }

  return (
    <div className="app">
      {showIntro && (
        <div className={`intro introStage${introStage}`}>
          <button
            className="skipIntro"
            onClick={() => setShowIntro(false)}
          >
            SKIP INTRO
          </button>

          <div className="introProgress">
            <span className={introStage >= 0 ? "active" : ""} />
            <span className={introStage >= 1 ? "active" : ""} />
            <span className={introStage >= 2 ? "active" : ""} />
            <span className={introStage >= 3 ? "active" : ""} />
            <span className={introStage >= 4 ? "active" : ""} />
            <span className={introStage >= 5 ? "active" : ""} />
          </div>

          <div className="introScene humanScene">
            <p className="introIndex">01 / HUMAN BODY</p>

            <div className="scienceVisual humanScienceVisual">
              <img
                src="/human-anatomy.png"
                alt="Human body anatomy"
                className="scienceImage humanImage"
              />
            </div>

            <div className="introCopy">
              <p className="introOverline">
                FROM BODY TO MOLECULE
              </p>

              <h1>
                Your body is made of trillions of cells.
              </h1>

              <p>
                Cells form tissues and organs throughout the
                body. Inside every cell, proteins carry out
                many of the processes needed for life.
              </p>
            </div>
          </div>

          <div className="introScene cellScene">
            <p className="introIndex">02 / CELL</p>

            <div className="scienceVisual cellScienceVisual">
              <img
                src="/animal-cell.png"
                alt="Animal cell showing major organelles"
                className="scienceImage cellImage"
              />
            </div>

            <div className="introCopy">
              <p className="introOverline">
                INSIDE A CELL
              </p>

              <h1>Cells depend on proteins.</h1>

              <p>
                Proteins perform many different jobs inside
                and around cells. Their function is closely
                connected to their three-dimensional structure.
              </p>
            </div>
          </div>

          <div className="introScene sequenceScene">
            <p className="introIndex">03 / AMINO ACIDS</p>

            <div className="scienceVisual aminoScienceVisual">
              <img
                src="/amino-acids.png"
                alt="Amino acids forming a protein chain"
                className="scienceImage aminoImage"
              />
            </div>

            <div className="introCopy">
              <p className="introOverline">
                BUILDING A PROTEIN
              </p>

              <h1>
                Proteins begin as chains of amino acids.
              </h1>

              <p>
                Amino acids are joined in a specific sequence.
                The chain can then fold into a complex
                three-dimensional structure.
              </p>
            </div>
          </div>

          <div className="introScene predictionScene">
            <p className="introIndex">
              04 / STRUCTURE PREDICTION
            </p>

            <div className="predictionVisual">
              <img
                src="/alphafold-prediction.png"
                alt="Illustration of a predicted folded protein structure"
                className="predictionImage"
              />

              <div className="predictionReadout">
                <span className="readoutDot" />
                PREDICTED STRUCTURE
              </div>
            </div>

            <div className="introCopy">
              <p className="introOverline">ALPHAFOLD</p>

              <h1>
                AlphaFold predicts a protein's 3D structure.
              </h1>

              <p>
                Starting from an amino-acid sequence, AlphaFold
                can generate a structural model and confidence
                information for the prediction.
              </p>
            </div>
          </div>

          <div className="introScene uncertaintyScene">
            <p className="introIndex">05 / UNCERTAINTY</p>

            <div className="uncertaintyVisual">
              <img
                src="/alphafold-uncertainty.png"
                alt="Illustration showing varying confidence across a predicted protein structure"
                className="uncertaintyImage"
              />
            </div>

            <div className="introCopy">
              <p className="introOverline">
                CONFIDENCE IS NOT UNIFORM
              </p>

              <h1>
                Not every part of the prediction is equally
                certain.
              </h1>

              <p>
                Some regions can have lower local confidence,
                while the relative positions of different
                parts of the protein can also be uncertain.
                A convincing 3D model still needs to be
                interpreted carefully.
              </p>
            </div>
          </div>

          <div className="introScene foldcheckScene">
            <div className="foldcheckFinalVisual">
              <div className="introFinalMark">F</div>

              <div className="finalEvidence">
                <div>
                  <strong>pLDDT</strong>
                  <span>Local confidence</span>
                </div>

                <div>
                  <strong>PAE</strong>
                  <span>Position uncertainty</span>
                </div>

                <div>
                  <strong>PDB</strong>
                  <span>Experimental evidence</span>
                </div>

                <div>
                  <strong>3D</strong>
                  <span>Structural checks</span>
                </div>
              </div>
            </div>

            <p className="introOverline">FOLDCHECK</p>

            <h1>Investigate the prediction.</h1>

            <p>
              FoldCheck brings different evidence together so
              you can explore a predicted structure, find
              regions worth investigating, and understand why
              they may need a closer look.
            </p>
          </div>
        </div>
      )}

      {!showIntro && (
  <FoldCheckNav activePage={activePage} setActivePage={setActivePage} />
)}

      <main>
        <section className="hero" id="top">
          <div className="heroContent">
            <p className="eyebrow">
              ALPHAFOLD STRUCTURE EXPLORER
            </p>

            <h1>FoldCheck</h1>

            <h2>Explore protein predictions in 3D.</h2>

            <p className="heroText">
              Search for a protein, explore its AlphaFold
              structure, and investigate where the prediction
              is confident — and where it might be uncertain.
              A protein model doesn't tell you everything,
              so FoldCheck brings multiple forms of evidence
              together to help you investigate it.
            </p>

            <p className="heroFeatures">
              Explore interactive 3D structures, examine
              local confidence (pLDDT), investigate predicted
              aligned error (PAE), inspect structural
              geometry, and compare predictions with
              available experimental evidence.
            </p>

            <div className="heroButtons">
              <a className="primaryButton" href="#explore">
                Explore a protein →
              </a>

              <button
                className="secondaryButton"
                onClick={replayIntro}
              >
                ▶ Replay Intro
              </button>
            </div>
          </div>
        </section>

        <section className="exploreSection" id="explore">
          <div className="exploreLayout">
            <div className="demoColumn">
              <p className="sectionLabel">
                INTERACTIVE DEMONSTRATION
              </p>

              <h2>Try the 3D viewer</h2>

              <p className="sectionText">
                Explore a real experimental protein structure
                before starting your own investigation.
              </p>

              <div className="demoViewer" id="demo">
                <div className="demoViewerHeader">
                  <span>TEST STRUCTURE</span>
                  <strong>PDB · 1CRN</strong>
                </div>

                <div className="viewerShell">
                  <ProteinViewer />
                </div>

                <div className="viewerBasicsBar">
                  <span>
                    <strong>Drag</strong> to rotate
                  </span>

                  <span>
                    <strong>Scroll</strong> to zoom
                  </span>

                  <span>
                    <strong>Click</strong> to select
                  </span>
                </div>
              </div>
            </div>

            <div className="searchColumn">
              <p className="sectionLabel">
                BEGIN AN INVESTIGATION
              </p>

              <h2>
                Which protein do you want to explore?
              </h2>

              <p className="sectionText">
                Search by protein name, gene, or UniProt
                accession. FoldCheck will retrieve available
                structural data and help you investigate
                the prediction.
              </p>

              <div className="searchBox">
                <input
                  type="text"
                  placeholder="Try insulin, TP53 or P04637..."
                  onKeyDown={(event) => {
  if (event.key === "Escape") {
    setShowProteinSuggestions(false);
    setActiveProteinSuggestionIndex(-1);
    return;
  }

  if (event.key === "ArrowDown" && proteinSuggestions.length > 0) {
    event.preventDefault();
    setShowProteinSuggestions(true);
    setActiveProteinSuggestionIndex((index) =>
      Math.min(index + 1, proteinSuggestions.length - 1)
    );
    return;
  }

  if (event.key === "ArrowUp" && proteinSuggestions.length > 0) {
    event.preventDefault();
    setActiveProteinSuggestionIndex((index) => Math.max(index - 1, 0));
    return;
  }

  if (event.key === "Enter") {
    event.preventDefault();

    if (
      showProteinSuggestions &&
      activeProteinSuggestionIndex >= 0 &&
      proteinSuggestions[activeProteinSuggestionIndex]
    ) {
      const protein = proteinSuggestions[activeProteinSuggestionIndex];
      
      openProtein(protein.primaryAccession);
    } else {
      searchProtein();
    }
  }
}}
                  value={proteinQuery}
                  onChange={(event) => {
  setProteinQuery(event.target.value);
  setShowProteinSuggestions(true);
  setActiveProteinSuggestionIndex(-1);
  
}}
                />

                <button
                  onClick={searchProtein}
                >
                  Search
                </button>
                {showProteinSuggestions && proteinSuggestions.length > 0 && (
  <div className="proteinSuggestions">
    {proteinSuggestions.map((protein, index) => {
      const accession = protein.primaryAccession;
      const name =
        protein.proteinDescription?.recommendedName?.fullName?.value ||
        protein.proteinDescription?.submissionNames?.[0]?.fullName?.value ||
        accession;

      const scientificName = protein.organism?.scientificName || "Unknown organism";
const commonName = protein.organism?.commonName || "";

const organism = commonName
  ? `${scientificName} — ${commonName}`
  : scientificName;

      return (
        <button
          key={accession}
          type="button"
          className={`proteinSuggestion ${
  activeProteinSuggestionIndex === index ? "isActive" : ""
}`}
          onClick={() => {
  setProteinQuery(name);
  openProtein(accession);
}}
        >
          <strong>{name}</strong>
          <span>{organism} · {accession}</span>
        </button>
      );
    })}
  </div>
)}
              </div>

{proteinSearchError && (
  <p className="proteinSearchError" role="alert">
    {proteinSearchError}
  </p>
)}

              <div className="examples">
                <span>Try an example</span>

                <button
                  className="exampleProtein"
                  onClick={() => openProtein("P04637")}
                >
                  <strong>p53</strong>
                  <small>TP53 · Tumour suppressor</small>
                </button>

                <button
  className="exampleProtein"
  onClick={() => openProtein("P02144")}
>
                  <strong>Myoglobin</strong>
                  <small>MB · Oxygen binding</small>
                </button>
              </div>
            </div>
          </div>
        </section>
      </main>
    </div>
  );
}

export default App;
