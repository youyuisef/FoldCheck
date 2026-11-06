import { useEffect, useRef, useState } from "react";
import { Viewer } from "molstar/lib/apps/viewer/app";
import { lociLabel } from "molstar/lib/mol-theme/label";
import {
  StructureElement,
  StructureProperties,
  StructureSelection
} from "molstar/lib/mol-model/structure";
import { OrderedSet } from "molstar/lib/mol-data/int";
import { MolScriptBuilder as MS } from "molstar/lib/mol-script/language/builder";
import { StructureQueryHelper } from "molstar/lib/mol-plugin-state/helpers/structure-query";
import "molstar/build/viewer/molstar.css";
import { Color } from "molstar/lib/mol-util/color";
import { ParamDefinition as PD } from "molstar/lib/mol-util/param-definition";

const spursPaeKnobs = {
  referenceResidue:PD.Numeric(1,{min:1,max:10000,step:1})
};

const LilywhiteTint = {
  name:"foldcheck-pae",
  label:"FoldCheck PAE",
  colorThemeCategory:"Custom",
  factory:(_ctx,props)=>({
    factory:LilywhiteTint,
    granularity:"group",
    color:(location)=>{
      if (!StructureElement.Location.is(location)) return Color(0x808080);

      const residue = StructureProperties.residue.auth_seq_id(location);
      const reference = Number(props.referenceResidue);
      const matrix = LilywhiteTint.currentMatrix;
      const error = matrix?.[reference -1]?.[residue -1];

      if (!Number.isFinite(error)) return Color(0x808080);

      
      const t = Math.max(0,Math.min(1,error/30));
      const r = Math.round(20 +235*t);
      const g = Math.round(100 +105*t);
      const b = Math.round(210 -170*t);

      return Color((r << 16) | (g << 8) | b);
    },
    props,
    description:"Predicted aligned error relative to the selected reference residue."
  }),
  getParams:()=>spursPaeKnobs,
  defaultValues:PD.getDefaultValues(spursPaeKnobs),
  isApplicable:(ctx)=>!!ctx.structure,
  currentMatrix:null
};

const binTheOverlay = (plugin,slot)=>
  plugin.state.data.build().delete(slot).commit();

const grabFold = (engine)=>
  engine.plugin.managers.structure.hierarchy.current.structures[0]?.cell?.obj?.data;

const walkTheChain = (foldedStruct,source)=>{
  const residueRows = [];
  const byKey = new Map();

  for (const unit of foldedStruct.units) {
    if (unit.kind!==0) continue; 

    for (const element of unit.elements) {
      const location = StructureElement.Location.create(
        foldedStruct,
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

      if (!byKey.has(key)) {
        const residue = {
          number,
          code,
          chain,
          atoms:[],
          plddt:null,
          backbone:{}
        };

        byKey.set(key,residue);
        residueRows.push(residue);
      }

      const residue = byKey.get(key);

      
      if (source==="alphafold"&&residue.plddt===null) {
        const score = StructureProperties.atom.B_iso_or_equiv(location);
        if (Number.isFinite(score)&&score >= 0&&score <= 100) {
          residue.plddt = score;
        }
      }

      if (!residue.atoms.includes(atomName)) {
        residue.atoms.push(atomName);
      }

      if (["N","CA","C"].includes(atomName)) {
        residue.backbone[atomName] = coordinates;
      }
    }
  }

  residueRows.sort((a,b)=>a.number -b.number);
  return residueRows;
};

const plddtTally = (structure)=>{
  const plddtList = [];

  for (const unit of structure.units) {
    if (unit.kind!==0) continue;

    const doneResidues = new Set();

    for (let i = 0; i < unit.elements.length; i++) {
      const location = StructureElement.Location.create(
        structure,
        unit,
        unit.elements[i]
      );

      const residueIndex = StructureProperties.residue.key(location);
      if (doneResidues.has(residueIndex)) continue;
      doneResidues.add(residueIndex);

      const score = StructureProperties.atom.B_iso_or_equiv(location);
      if (Number.isFinite(score)&&score >= 0&&score <= 100) {
        plddtList.push(score);
      }
    }
  }

  if (plddtList.length===0) return null;

  const mean = plddtList.reduce((sum,score)=>sum +score,0)/plddtList.length;

  return {
    mean,
    residueCount:plddtList.length,
    veryHigh:plddtList.filter((s)=>s >= 90).length,
    confident:plddtList.filter((s)=>s >= 70&&s < 90).length,
    low:plddtList.filter((s)=>s >= 50&&s < 70).length,
    veryLow:plddtList.filter((s)=>s < 50).length
  };
};

const readClickTag = (label)=>{
  const parts = label.split("|").map((part)=>part.trim());

  let atom = null;
  let residueCode = "";
  let residueNumber = "";
  let chain = "";

  if (parts.length >= 3) {
    
    const atomMatch = parts[0].match(/^([A-Z0-9']+)/);
    const residueMatch = parts[1].match(/^([A-Z]{3})\s+(-?\d+)/);

    atom = atomMatch?.[1]||null;
    residueCode = residueMatch?.[1]||"";
    residueNumber = residueMatch?.[2]||"";
    chain = parts[2]||"";
  } else {
    
    const residueMatch = parts[0].match(/^([A-Z]{3})\s+(-?\d+)/);

    residueCode = residueMatch?.[1]||"";
    residueNumber = residueMatch?.[2]||"";
    chain = parts[1]||"";
  }

  return {atom,residueCode,residueNumber,chain};
};

const plddtAtClick = (loci,source)=>{
  let plddt = null;

  if (StructureElement.Loci.is(loci)&&loci.elements.length > 0) {
    const first = loci.elements[0];
    const index =
      OrderedSet.size(first.indices) > 0
        ? OrderedSet.getAt(first.indices,0)
        :undefined;

    if (index!==undefined) {
      const location = StructureElement.Location.create(
        loci.structure,
        first.unit,
        first.unit.elements[index]
      );

      const value = StructureProperties.atom.B_iso_or_equiv(location);

      if (Number.isFinite(value)&&value >= 0&&value <= 100&&source==="alphafold") {
        plddt = value;
      }
    }
  }

  return plddt;
};

const huntRes = (structure,residueNumber,residueCode,atomName)=>{
  let match = null;
  let foundResidue = false;

  for (const unit of structure.units) {
    if (unit.kind!==0) continue;

    for (let i = 0; i < unit.elements.length; i++) {
      const location = StructureElement.Location.create(
        structure,
        unit,
        unit.elements[i]
      );

      const number = StructureProperties.residue.auth_seq_id(location);
      const code = StructureProperties.atom.label_comp_id(location);
      const name = StructureProperties.atom.label_atom_id(location);

      if (number!==residueNumber) continue;
      if (residueCode&&code!==residueCode) continue;

      foundResidue = true;

      if (atomName&&name!==atomName) continue;

      match = {unit,index:i,location,code,name,number};
      break;
    }

    if (match) break;
  }

  return {match,foundResidue};
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
  
  const lastClickRef = useRef(onResidueSelect);
  lastClickRef.current = onResidueSelect;

  const viewportBox = useRef(null);
  const molEngine = useRef(null);
  const [molReady,setMolReady] = useState(0);

  const resHighlight = useRef(null);
  const atomHighlight = useRef(null);
  const regionHighlight = useRef(null);

  
  
  
  useEffect(()=>{
    let quitEarly = false;
    let engine = null;
    let clickSub = null;

    const holder = viewportBox.current;
    if (!holder) return;

    const mount = document.createElement("div");
    mount.style.width = "100%";
    mount.style.height = "100%";
    holder.appendChild(mount);

    const startMol = async ()=>{
      try {
        engine = await Viewer.create(mount,{
          layoutIsExpanded:false,
          layoutShowControls:false,
          layoutShowSequence:false,
          layoutShowLog:false,
          layoutShowLeftPanel:false,
          layoutShowRemoteState:false,
          viewportShowExpand:false,
          viewportShowControls:true,
          viewportShowSelectionMode:false,
          viewportBackgroundColor:"#05090d"
        });

        
        if (quitEarly) {
          engine.dispose();
          mount.remove();
          return;
        }

        molEngine.current = engine;

        engine.plugin.canvas3d?.setProps({
          renderer:{backgroundColor:Color(0x07131f)}
        });

        
        const themeShelf = engine.plugin.representation.structure.themes.colorThemeRegistry;
        if (!themeShelf.types.some(([name])=>name==="foldcheck-pae")) {
          themeShelf.add(LilywhiteTint);
        }

        if (source==="alphafold") {
          await engine.loadAlphaFoldDb(proteinId);
        } else {
          await engine.loadPdb(proteinId);
        }

        if (quitEarly) return;
        setMolReady((n)=>n +1);

        
        const foldedStruct = grabFold(engine);

        if (foldedStruct) {
          const residueRows = walkTheChain(foldedStruct,source);
          onGeometryLoaded?.(residueRows);
          onSequenceLoaded?.(residueRows);
        }

        
        if (source==="alphafold") {
          const structure = grabFold(engine);

          const summary = structure&&plddtTally(structure);

          if (summary) onProteinConfidence?.(summary);
        }

        
        clickSub = engine.subscribe(
          engine.plugin.behaviors.interaction.click,
          ({current})=>{
            if (current?.loci?.kind!=="element-loci") return;

            const label = lociLabel(current.loci,{
              reverse:true,
              hidePrefix:true,
              htmlStyling:false
            });

            const {atom,residueCode,residueNumber,chain} = readClickTag(label);
            const plddt = plddtAtClick(current.loci,source);

            lastClickRef.current?.({
              residueCode,
              residueNumber,
              chain,
              atom,
              plddt
            });
          }
        );
      } catch (error) {
        if (!quitEarly) console.error("FoldCheck Mol*initialization failed:",error);
      }
    };

    startMol();

    
    return ()=>{
      quitEarly = true;
      clickSub?.unsubscribe();

      if (molEngine.current===engine) molEngine.current = null;

      resHighlight.current = null;
      atomHighlight.current = null;
      regionHighlight.current = null;

      if (engine) {
        engine.dispose();
        mount.remove();
      } else {
        mount.style.display = "none";
      }
    };
  },[source,proteinId]);

  
  
  
  useEffect(()=>{
    if (!searchRequest) return;

    const engine = molEngine.current;

    if (!engine) {
      onSearchResult?.({
        success:false,
        message:"The 3D structure is still loading."
      });
      return;
    }

    const structure = grabFold(engine);

    if (!structure) {
      onSearchResult?.({
        success:false,
        message:"No loaded structure was found."
      });
      return;
    }

    const {residueNumber,residueCode,atomName} = searchRequest;

    const {match,foundResidue} = huntRes(structure,residueNumber,residueCode,atomName);

    if (!match) {
      onSearchResult?.({
        success:false,
        message:foundResidue
          ? `Residue ${residueNumber} was found, but atom ${atomName} was not.`
          : `No matching residue ${residueNumber} was found.`
      });
      return;
    }

    const plugin = engine.plugin;

    
    
    const addResidueOverlay = async ()=>{
      if (resHighlight.current) {
        await binTheOverlay(plugin,resHighlight.current);
        resHighlight.current = null;
      }

      const structureRef = plugin.managers.structure.hierarchy.current.structures[0];
      if (!structureRef) return;

      const expression = MS.struct.generator.atomGroups({
        "residue-test":MS.core.rel.eq([
          MS.struct.atomProperty.macromolecular.auth_seq_id(),
          residueNumber
        ]),
        "group-by":MS.struct.atomProperty.macromolecular.residueKey()
      });

      const component = await plugin.builders.structure.tryCreateComponentFromExpression(
        structureRef.cell,
        expression,
        "foldcheck-residue-overlay",
        { label: `${match.code} ${match.number}` }
      );

      if (!component) return;

      if (atomHighlight.current) {
        await binTheOverlay(plugin,atomHighlight.current);
        atomHighlight.current = null;
      }

      resHighlight.current = component.ref;

      await plugin.builders.structure.representation.addRepresentation(component,{
        type:"ball-and-stick",
        color:"element-symbol"
      });

      if (atomName) {
        const atomExpression = MS.struct.generator.atomGroups({
          "atom-test":MS.core.logic.and([
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

        const atomComponent = await plugin.builders.structure.tryCreateComponentFromExpression(
          structureRef.cell,
          atomExpression,
          "foldcheck-target-atom",
          { label: `${match.code} ${match.number} ${atomName}` }
        );

        if (atomComponent) {
          atomHighlight.current = atomComponent.ref;

          await plugin.builders.structure.representation.addRepresentation(atomComponent,{
            type:"ball-and-stick",
            color:"uniform",
            colorParams:{value:0xff8c00},
            size:"uniform",
            sizeParams:{value:0.45}
          });
        }
      }
    };

    addResidueOverlay().catch((error)=>{
      console.error("FoldCheck residue overlay failed:",error);
    });

    
    const loci = StructureElement.Loci(structure,[
      {
        unit:match.unit,
        indices:OrderedSet.ofSingleton(match.index)
      }
    ]);

    engine.plugin.managers.interactivity.lociSelects.deselectAll();
    engine.plugin.managers.interactivity.lociSelects.select({loci});

    engine.plugin.managers.camera.focusLoci(loci,{
      minRadius:atomName ? 8 :12,
      extraRadius:atomName ? 6 :8,
      durationMs:350
    });

    const score = StructureProperties.atom.B_iso_or_equiv(match.location);

    const plddt =
      source==="alphafold"&&Number.isFinite(score)&&score >= 0&&score <= 100
        ? score
        :null;

    lastClickRef.current?.({
      residueCode:match.code,
      residueNumber:String(match.number),
      chain:StructureProperties.chain.auth_asym_id(match.location),
      atom:atomName ? match.name :null,
      plddt
    });

    onSearchResult?.({
      success:true,
      atomSearch:Boolean(atomName),
      message:atomName
        ? `Found ${match.code} ${match.number} · ${match.name}`
        : `Found ${match.code} ${match.number}`
    });
  },[searchRequest,source]);

  
  
  
  useEffect(()=>{
    if (!selectedRegion) return;

    const engine = molEngine.current;
    if (!engine) return;

    const plugin = engine.plugin;
    const structure = grabFold(engine);
    if (!structure) return;

    const start = Number(selectedRegion.start);
    const end = Number(selectedRegion.end);

    if (!Number.isInteger(start)||!Number.isInteger(end)||start > end) return;

    let quitEarly = false;

    const highlightRegion = async ()=>{
      if (regionHighlight.current) {
        await binTheOverlay(plugin,regionHighlight.current);
        regionHighlight.current = null;
      }

      if (quitEarly) return;

      const expression = MS.struct.generator.atomGroups({
        "residue-test":MS.core.rel.inRange([
          MS.struct.atomProperty.macromolecular.auth_seq_id(),
          start,
          end
        ]),
        "group-by":MS.struct.atomProperty.macromolecular.residueKey()
      });

      const selection = StructureQueryHelper.createAndRun(structure,expression);
      const loci = StructureSelection.toLociWithSourceUnits(selection.selection);

      plugin.managers.interactivity.lociSelects.deselectAll();
      plugin.managers.interactivity.lociSelects.select({loci});

      plugin.managers.camera.focusLoci(loci,{
        minRadius:12,
        extraRadius:8,
        durationMs:450
      });
    };

    highlightRegion().catch((error)=>{
      console.error("FoldCheck region highlight failed:",error);
    });

    return ()=>{
      quitEarly = true;
    };
  },[selectedRegion,proteinId]);

  
  
  
  useEffect(()=>{
    const engine = molEngine.current;
    if (!engine||!paeMatrix?.length||!paeReferenceResidue) return;

    const plugin = engine.plugin;
    const structures = plugin.managers.structure.hierarchy.current.structures;
    const components = structures.flatMap((structure)=>structure.components);

    if (!components.length) return;

    LilywhiteTint.currentMatrix = paeMatrix;

    plugin.managers.structure.component
      .updateRepresentationsTheme(components,{
        color:"foldcheck-pae",
        colorParams:{referenceResidue:Number(paeReferenceResidue)}
      })
      .catch((error)=>console.error("FoldCheck PAE colouring failed:",error));
  },[paeMatrix,paeReferenceResidue,molReady]);

  
  
  
  useEffect(()=>{
    if (!linkedResidue||!molReady) return;

    const plugin = molEngine.current?.plugin;
    const structure = plugin?.managers.structure.hierarchy.current.structures[0]?.cell?.obj?.data;
    const number = Number(linkedResidue.number ?? linkedResidue.residueNumber);

    if (!structure||!Number.isFinite(number)) return;

    try {
      const tests = {
        "residue-test":MS.core.rel.eq([
          MS.struct.atomProperty.macromolecular.auth_seq_id(),
          number
        ]),
        "group-by":MS.struct.atomProperty.macromolecular.residueKey()
      };

      if (linkedResidue.chain) {
        tests["chain-test"] = MS.core.rel.eq([
          MS.struct.atomProperty.macromolecular.auth_asym_id(),
          linkedResidue.chain
        ]);
      }

      const expression = MS.struct.generator.atomGroups(tests);
      const result = StructureQueryHelper.createAndRun(structure,expression);
      const loci = StructureSelection.toLociWithSourceUnits(result.selection);

      if (!loci.elements.length) {
        console.warn("FoldCheck:linked residue not found",linkedResidue);
        return;
      }

      plugin.managers.interactivity.lociSelects.deselectAll();
      plugin.managers.interactivity.lociSelects.select({loci});

      plugin.managers.camera.focusLoci(loci,{
        minRadius:12,
        extraRadius:8,
        durationMs:350
      });
    } catch (error) {
      console.error("FoldCheck linked selection failed:",error);
    }
  },[linkedResidue,molReady]);

  
  
  
  useEffect(()=>{
    if (!resetViewRequest) return;

    const engine = molEngine.current;
    if (!engine) return;

    const plugin = engine.plugin;

    plugin.managers.interactivity.lociHighlights.clearHighlights();
    plugin.managers.interactivity.lociSelects.deselectAll();

    
    const overlayRefs = [
      ["residue overlay",resHighlight],
      ["atom overlay",atomHighlight],
      ["region overlay",regionHighlight]
    ];

    for (const [label,slot] of overlayRefs) {
      if (!slot.current) continue;

      const overlay = slot.current;
      slot.current = null;

      binTheOverlay(plugin,overlay).catch((error)=>{
        console.error(`FoldCheck ${label} reset failed:`, error);
      });
    }

    plugin.canvas3d?.requestCameraReset();
  },[resetViewRequest]);

  return (
    <div
      ref={viewportBox}
      style={{
        position:"relative",
        width:"100%",
        height:"600px",
        overflow:"hidden",
        borderRadius:"24px"
      }}
    />
  );
}

export default ProteinViewer;
