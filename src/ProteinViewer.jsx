import { useEffect, useRef, useState } from "react";
import { Viewer } from "molstar/lib/apps/viewer/app";
import { lociLabel } from "molstar/lib/mol-theme/label";
import { StructureElement, StructureProperties } from "molstar/lib/mol-model/structure";
import { OrderedSet } from "molstar/lib/mol-data/int";
import { MolScriptBuilder as MS } from "molstar/lib/mol-script/language/builder";
import { StructureSelection } from "molstar/lib/mol-model/structure";
import { StructureQueryHelper } from "molstar/lib/mol-plugin-state/helpers/structure-query";
import "molstar/build/viewer/molstar.css";
import { Color } from "molstar/lib/mol-util/color";
import { ParamDefinition as PD } from "molstar/lib/mol-util/param-definition";

// PAE colours for the selected residue
const paeThemeParams = {
  referenceResidue: PD.Numeric(1, { min: 1, max: 10000, step: 1 })
};

const FoldCheckPAETheme = {
  name: "foldcheck-pae",
  label: "FoldCheck PAE",
  colorThemeCategory: "Custom",
  factory: (_ctx, props) => ({
    factory: FoldCheckPAETheme,
    granularity: "group",
    color: (location) => {
      if (!StructureElement.Location.is(location)) return Color(0x808080);
      const residue = StructureProperties.residue.auth_seq_id(location);
      const reference = Number(props.referenceResidue);
      const matrix = FoldCheckPAETheme.currentMatrix;
      const error = matrix?.[reference - 1]?.[residue - 1];
      if (!Number.isFinite(error)) return Color(0x808080);
      const t = Math.max(0, Math.min(1, error / 30));
      const r = Math.round(20 + 235 * t);
      const g = Math.round(100 + 105 * t);
      const b = Math.round(210 - 170 * t);
      return Color((r << 16) | (g << 8) | b);
    },
    props,
    description: "Predicted aligned error relative to the selected reference residue."
  }),
  getParams: () => paeThemeParams,
  defaultValues: PD.getDefaultValues(paeThemeParams),
  isApplicable: (ctx) => !!ctx.structure,
  currentMatrix: null
};

function ProteinViewer({
  source = "pdb",
  proteinId = "1CRN",
  onResidueSelect,
  onProteinConfidence,
  onSequenceLoaded,
  onGeometryLoaded,
  searchRequest,
  onSearchResult,
  resetViewRequest,
  selectedRegion,
  linkedResidue,
  paeMatrix,
  paeReferenceResidue
}) {
    const onResidueSelectRef = useRef(onResidueSelect);
    onResidueSelectRef.current = onResidueSelect;
    const viewerRef = useRef(null);
    const molstarRef = useRef(null);
    const [viewerReady, setViewerReady] = useState(0);
    const overlayRef = useRef(null);
    const atomOverlayRef = useRef(null);
    const regionOverlayRef = useRef(null);

    // Load the 3D model
    useEffect(() => {
      let cancelled = false;
      let viewer = null;
      let clickSubscription = null;
      const container = viewerRef.current;
      if (!container) return;

      const host = document.createElement("div");
      host.style.width = "100%";
      host.style.height = "100%";
      container.appendChild(host);

      const startViewer = async () => {
        try {
          viewer = await Viewer.create(host, {
            layoutIsExpanded: false,
            layoutShowControls: false,
            layoutShowSequence: false,
            layoutShowLog: false,
            layoutShowLeftPanel: false,
            layoutShowRemoteState: false,
            viewportShowExpand: false,
            viewportShowControls: true,
            viewportShowSelectionMode: false,
            viewportBackgroundColor: "#05090d"
          });
          if (cancelled) {
            viewer.dispose();
            host.remove();
            return;
          }
          molstarRef.current = viewer;
        viewer.plugin.canvas3d?.setProps({
          renderer: { backgroundColor: Color(0x07131f) }
        });
          const registry = viewer.plugin.representation.structure.themes.colorThemeRegistry;
          if (!registry.types.some(([name]) => name === "foldcheck-pae")) {
            registry.add(FoldCheckPAETheme);
          }
          if (source === "alphafold") {
            await viewer.loadAlphaFoldDb(proteinId);
          } else {
            await viewer.loadPdb(proteinId);
          }
          if (cancelled) return;
          setViewerReady((n) => n + 1);

// Get the sequence and backbone positions
const loadedStructure =
  viewer.plugin.managers.structure.hierarchy.current.structures[0]
    ?.cell?.obj?.data;

if (loadedStructure) {
  const residues = [];
  const residueMap = new Map();
  

  for (const unit of loadedStructure.units) {
    if (unit.kind !== 0) continue;

    for (const element of unit.elements) {
      const location = StructureElement.Location.create(
        loadedStructure,
        unit,
        element
      );

      const number = StructureProperties.residue.auth_seq_id(location);
      const code = StructureProperties.atom.label_comp_id(location);
      const chain = StructureProperties.chain.auth_asym_id(location);
      const key = `${chain}:${number}`;

      const atomName = StructureProperties.atom.label_atom_id(location);
      const coordinates = [
  unit.conformation.x(element),
  unit.conformation.y(element),
  unit.conformation.z(element)
];

      if (!residueMap.has(key)) {
        const residue = {
          number,
          code,
          chain,
          atoms: [],
          plddt: null,
backbone: {}
        };

        residueMap.set(key, residue);
        residues.push(residue);
      }

const residue = residueMap.get(key);
if (source === "alphafold" && residue.plddt === null) {
  const score = StructureProperties.atom.B_iso_or_equiv(location);
  if (Number.isFinite(score) && score >= 0 && score <= 100) {
    residue.plddt = score;
  }
}

if (!residue.atoms.includes(atomName)) {
  residue.atoms.push(atomName);
}

if (["N", "CA", "C"].includes(atomName)) {
  residue.backbone[atomName] = coordinates;
}

    }
  }

  residues.sort((a, b) => a.number - b.number);
  onGeometryLoaded?.(residues);
  onSequenceLoaded?.(residues);
}

// Get AlphaFold confidence scores
if (source === "alphafold") {
  const structures = viewer.plugin.managers.structure.hierarchy.current.structures;
  const structure = structures[0]?.cell?.obj?.data;

  if (structure) {
    const scores = [];

    for (const unit of structure.units) {
      if (unit.kind !== 0) continue;

      const seenResidues = new Set();

      for (let i = 0; i < unit.elements.length; i++) {
        const location = StructureElement.Location.create(
          structure,
          unit,
          unit.elements[i]
        );

        const residueIndex = StructureProperties.residue.key(location);

        if (seenResidues.has(residueIndex)) continue;
        seenResidues.add(residueIndex);

        const score = StructureProperties.atom.B_iso_or_equiv(location);

        if (Number.isFinite(score) && score >= 0 && score <= 100) {
          scores.push(score);
        }
      }
    }

    if (scores.length > 0) {
      const mean = scores.reduce((sum, score) => sum + score, 0) / scores.length;

      onProteinConfidence?.({
        mean,
        residueCount: scores.length,
        veryHigh: scores.filter(score => score >= 90).length,
        confident: scores.filter(score => score >= 70 && score < 90).length,
        low: scores.filter(score => score >= 50 && score < 70).length,
        veryLow: scores.filter(score => score < 50).length
      });
    }
  }
}

// Handle clicks on the model
clickSubscription = viewer.subscribe(
  viewer.plugin.behaviors.interaction.click,
  ({ current }) => {
    if (current?.loci?.kind !== "element-loci") return;

    const label = lociLabel(current.loci, {
      reverse: true,
      hidePrefix: true,
      htmlStyling: false
    });

    const parts = label.split("|").map(part => part.trim());

let atom = null;
let residueCode = "";
let residueNumber = "";
let chain = "";

if (parts.length >= 3) {
  const atomMatch = parts[0].match(/^([A-Z0-9']+)/);
  const residueMatch = parts[1].match(/^([A-Z]{3})\s+(-?\d+)/);

  atom = atomMatch?.[1] || null;
  residueCode = residueMatch?.[1] || "";
  residueNumber = residueMatch?.[2] || "";
  chain = parts[2] || "";
} else {
  const residueMatch = parts[0].match(/^([A-Z]{3})\s+(-?\d+)/);

  residueCode = residueMatch?.[1] || "";
  residueNumber = residueMatch?.[2] || "";
  chain = parts[1] || "";
}

let plddt = null;

const loci = current.loci;

if (StructureElement.Loci.is(loci) && loci.elements.length > 0) {
  const first = loci.elements[0];
  const index =
  OrderedSet.size(first.indices) > 0
    ? OrderedSet.getAt(first.indices, 0)
    : undefined;

  if (index !== undefined) {
    const location = StructureElement.Location.create(
      loci.structure,
      first.unit,
      first.unit.elements[index]
    );

    const value = StructureProperties.atom.B_iso_or_equiv(location);

    if (Number.isFinite(value) && value >= 0 && value <= 100 && source === "alphafold") {
      plddt = value;
    }
  }
}

console.log("FoldCheck residue confidence:", {
  residueCode,
  residueNumber,
  plddt
});

onResidueSelectRef.current?.({
  residueCode,
  residueNumber,
  chain,
  atom,
  plddt
});
  }
);

        } catch (error) {
          if (!cancelled) console.error("FoldCheck Mol* initialization failed:", error);
        }
      };
      startViewer();
      return () => {
        cancelled = true;
        clickSubscription?.unsubscribe();
        if (molstarRef.current === viewer) molstarRef.current = null;
        overlayRef.current = null;
        atomOverlayRef.current = null;
        regionOverlayRef.current = null;
        if (viewer) {
          viewer.dispose();
          host.remove();
        } else {
          host.style.display = "none";
        }
      };
    }, [source, proteinId]);

// Search and zoom to an atom or residue
useEffect(() => {
  if (!searchRequest) return;

  const viewer = molstarRef.current;

  if (!viewer) {
    onSearchResult?.({
      success: false,
      message: "The 3D structure is still loading."
    });
    return;
  }

  const structure =
    viewer.plugin.managers.structure.hierarchy.current.structures[0]
      ?.cell?.obj?.data;

  if (!structure) {
    onSearchResult?.({
      success: false,
      message: "No loaded structure was found."
    });
    return;
  }

  const { residueNumber, residueCode, atomName } = searchRequest;

  let match = null;
  let foundResidue = false;

  for (const unit of structure.units) {
    if (unit.kind !== 0) continue;

    for (let i = 0; i < unit.elements.length; i++) {
      const location = StructureElement.Location.create(
        structure,
        unit,
        unit.elements[i]
      );

      const number = StructureProperties.residue.auth_seq_id(location);
      const code = StructureProperties.atom.label_comp_id(location);
      const name = StructureProperties.atom.label_atom_id(location);

      if (number !== residueNumber) continue;
      if (residueCode && code !== residueCode) continue;

      foundResidue = true;

      if (atomName && name !== atomName) continue;

      match = { unit, index: i, location, code, name, number };
      break;
    }

    if (match) break;
  }

  if (!match) {
    onSearchResult?.({
      success: false,
      message: foundResidue
        ? `Residue ${residueNumber} was found, but atom ${atomName} was not.`
        : `No matching residue ${residueNumber} was found.`
    });
    return;
  }

  const plugin = viewer.plugin;

const addResidueOverlay = async () => {
  if (overlayRef.current) {
    await plugin.state.data
      .build()
      .delete(overlayRef.current)
      .commit();

    overlayRef.current = null;
  }

  const structureRef =
    plugin.managers.structure.hierarchy.current.structures[0];

  if (!structureRef) return;

  const expression = MS.struct.generator.atomGroups({
    "residue-test": MS.core.rel.eq([
      MS.struct.atomProperty.macromolecular.auth_seq_id(),
      residueNumber
    ]),
    "group-by": MS.struct.atomProperty.macromolecular.residueKey()
  });

  const component =
    await plugin.builders.structure.tryCreateComponentFromExpression(
      structureRef.cell,
      expression,
      "foldcheck-residue-overlay",
      { label: `${match.code} ${match.number}` }
    );

  if (!component) return;

    if (atomOverlayRef.current) {
    await plugin.state.data
      .build()
      .delete(atomOverlayRef.current)
      .commit();

    atomOverlayRef.current = null;
  }

  overlayRef.current = component.ref;

  await plugin.builders.structure.representation.addRepresentation(
    component,
    {
      type: "ball-and-stick",
      color: "element-symbol"
    }
  );

    if (atomName) {
    const atomExpression = MS.struct.generator.atomGroups({
      "atom-test": MS.core.logic.and([
        MS.core.rel.eq([
          MS.struct.atomProperty.macromolecular.auth_seq_id(),
          residueNumber
        ]),
        MS.core.rel.eq([
          MS.struct.atomProperty.macromolecular.label_atom_id(),
          atomName
        ])
      ])
    });

    const atomComponent =
      await plugin.builders.structure.tryCreateComponentFromExpression(
        structureRef.cell,
        atomExpression,
        "foldcheck-target-atom",
        { label: `${match.code} ${match.number} ${atomName}` }
      );

    if (atomComponent) {
      atomOverlayRef.current = atomComponent.ref;

      await plugin.builders.structure.representation.addRepresentation(
        atomComponent,
        {
          type: "ball-and-stick",
          color: "uniform",
          colorParams: {
            value: 0xff8c00
          },
          size: "uniform",
          sizeParams: {
            value: 0.45
          }
        }
      );
    }
  }

};

addResidueOverlay().catch(error => {
  console.error("FoldCheck residue overlay failed:", error);
});

  const loci = StructureElement.Loci(structure, [
    {
      unit: match.unit,
      indices: OrderedSet.ofSingleton(match.index)
    }
  ]);

viewer.plugin.managers.interactivity.lociSelects.deselectAll();

viewer.plugin.managers.interactivity.lociSelects.select({
  loci
});

viewer.plugin.managers.camera.focusLoci(loci, {
  minRadius: atomName ? 8 : 12,
  extraRadius: atomName ? 6 : 8,
  durationMs: 350
});

  const score = StructureProperties.atom.B_iso_or_equiv(
    match.location
  );

  const plddt =
    source === "alphafold" &&
    Number.isFinite(score) &&
    score >= 0 &&
    score <= 100
      ? score
      : null;

  onResidueSelectRef.current?.({
    residueCode: match.code,
    residueNumber: String(match.number),
    chain: StructureProperties.chain.auth_asym_id(match.location),
    atom: atomName ? match.name : null,
    plddt
  });

  onSearchResult?.({
    success: true,
    atomSearch: Boolean(atomName),
    message: atomName
      ? `Found ${match.code} ${match.number} · ${match.name}`
      : `Found ${match.code} ${match.number}`
  });
}, [searchRequest, source]);
// Highlight a selected region
useEffect(() => {
  if (!selectedRegion) return;

  const viewer = molstarRef.current;
  if (!viewer) return;

  const plugin = viewer.plugin;

  const structureRef =
    plugin.managers.structure.hierarchy.current.structures[0];

  const structure = structureRef?.cell?.obj?.data;
  if (!structure) return;

  const start = Number(selectedRegion.start);
  const end = Number(selectedRegion.end);

  if (!Number.isInteger(start) || !Number.isInteger(end) || start > end) {
    return;
  }

  let cancelled = false;

  const highlightRegion = async () => {
    
    if (regionOverlayRef.current) {
      await plugin.state.data
        .build()
        .delete(regionOverlayRef.current)
        .commit();

      regionOverlayRef.current = null;
    }

    if (cancelled) return;

    const expression = MS.struct.generator.atomGroups({
  "residue-test": MS.core.rel.inRange([
    MS.struct.atomProperty.macromolecular.auth_seq_id(),
    start,
    end
  ]),
  "group-by": MS.struct.atomProperty.macromolecular.residueKey()
});
    
    const selection = StructureQueryHelper.createAndRun(
      structure,
      expression
    );

    const loci = StructureSelection.toLociWithSourceUnits(
      selection.selection
    );

    plugin.managers.interactivity.lociSelects.deselectAll();

plugin.managers.interactivity.lociSelects.select({ loci });

    plugin.managers.camera.focusLoci(loci, {
      minRadius: 12,
      extraRadius: 8,
      durationMs: 450
    });
  };

  highlightRegion().catch((error) => {
    console.error("FoldCheck region highlight failed:", error);
  });

  return () => {
    cancelled = true;
  };
}, [selectedRegion, proteinId]);

// Update PAE colours
useEffect(() => {
  const viewer = molstarRef.current;
  if (!viewer || !paeMatrix?.length || !paeReferenceResidue) return;
  const plugin = viewer.plugin;
  const structures = plugin.managers.structure.hierarchy.current.structures;
  const components = structures.flatMap((structure) => structure.components);
  if (!components.length) return;

  FoldCheckPAETheme.currentMatrix = paeMatrix;
  plugin.managers.structure.component.updateRepresentationsTheme(components, {
    color: "foldcheck-pae",
    colorParams: { referenceResidue: Number(paeReferenceResidue) }
  }).catch((error) => console.error("FoldCheck PAE colouring failed:", error));
}, [paeMatrix, paeReferenceResidue, viewerReady]);


  useEffect(() => {
    if (!linkedResidue || !viewerReady) return;
    const plugin = molstarRef.current?.plugin;
    const structure = plugin?.managers.structure.hierarchy.current.structures[0]?.cell?.obj?.data;
    const number = Number(linkedResidue.number ?? linkedResidue.residueNumber);
    if (!structure || !Number.isFinite(number)) return;
    try {
      const tests = {
        "residue-test": MS.core.rel.eq([
          MS.struct.atomProperty.macromolecular.auth_seq_id(),
          number
        ]),
        "group-by": MS.struct.atomProperty.macromolecular.residueKey()
      };
      if (linkedResidue.chain) {
        tests["chain-test"] = MS.core.rel.eq([
          MS.struct.atomProperty.macromolecular.auth_asym_id(),
          linkedResidue.chain
        ]);
      }
      const expression = MS.struct.generator.atomGroups(tests);
      const result = StructureQueryHelper.createAndRun(structure, expression);
      const loci = StructureSelection.toLociWithSourceUnits(result.selection);
      if (!loci.elements.length) {
        console.warn("FoldCheck: linked residue not found", linkedResidue);
        return;
      }
      plugin.managers.interactivity.lociSelects.deselectAll();
      plugin.managers.interactivity.lociSelects.select({ loci });
      plugin.managers.camera.focusLoci(loci, {
        minRadius: 12,
        extraRadius: 8,
        durationMs: 350
      });
    } catch (error) {
      console.error("FoldCheck linked selection failed:", error);
    }
  }, [linkedResidue, viewerReady]);

useEffect(() => {
  if (!resetViewRequest) return;

  const viewer = molstarRef.current;
  if (!viewer) return;

 viewer.plugin.managers.interactivity.lociHighlights.clearHighlights();
viewer.plugin.managers.interactivity.lociSelects.deselectAll();
if (overlayRef.current) {
  const overlay = overlayRef.current;
  overlayRef.current = null;

  viewer.plugin.state.data
    .build()
    .delete(overlay)
    .commit()
    .catch(error => {
      console.error("FoldCheck overlay reset failed:", error);
    });
}

if (atomOverlayRef.current) {
  const atomOverlay = atomOverlayRef.current;
  atomOverlayRef.current = null;

  viewer.plugin.state.data
    .build()
    .delete(atomOverlay)
    .commit()
    .catch(error => {
      console.error("FoldCheck atom overlay reset failed:", error);
    });
}

if (regionOverlayRef.current) {
  const regionOverlay = regionOverlayRef.current;
  regionOverlayRef.current = null;

  viewer.plugin.state.data
    .build()
    .delete(regionOverlay)
    .commit()
    .catch(error => {
      console.error("FoldCheck region overlay reset failed:", error);
    });
}

viewer.plugin.canvas3d?.requestCameraReset();
}, [resetViewRequest]);

    return (
        <div
        ref={viewerRef}
        style={{
            position: "relative",
            width: "100%",
            height: "600px",
            overflow: "hidden",
            borderRadius: "24px"
        }}
        />
    );
}

export default ProteinViewer;
