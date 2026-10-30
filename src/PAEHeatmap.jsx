
import { useEffect, useRef, useState } from "react";

export default function PAEHeatmap({ matrix, onReferenceChange, onPairSelect }) {
  const canvasRef = useRef(null);
  const [hover, setHover] = useState(null);
const [selectedPair, setSelectedPair] = useState(null);

  useEffect(() => {
    if (!matrix?.length || !canvasRef.current) return;

    const size = matrix.length;
    const canvas = canvasRef.current;
    const context = canvas.getContext("2d");

    canvas.width = size;
    canvas.height = size;

    const image = context.createImageData(size, size);

    for (let y = 0; y < size; y++) {
      for (let x = 0; x < size; x++) {
        const error = Math.max(
          0,
          Math.min(30, Number(matrix[y]?.[x] ?? 30))
        );

        const t = error / 30;
        const index = (y * size + x) * 4;

        image.data[index] = Math.round(20 + 235 * t);
        image.data[index + 1] = Math.round(100 + 105 * t);
        image.data[index + 2] = Math.round(210 - 170 * t);
        image.data[index + 3] = 255;
      }
    }

    context.putImageData(image, 0, 0);
  }, [matrix]);

  if (!matrix?.length) return null;

  const size = matrix.length;

  return (
    <div className="paeHeatmap">
      <h2>Predicted Aligned Error (PAE)</h2>

      <p>
        Explore AlphaFold's predicted error in the relative
        positions of pairs of residues.
      </p>

      <p>
        <strong>Protein length:</strong> {size} residues
      </p>

      <div
  style={{
    width: "100%",
    maxWidth: "550px",
    display: "flex",
    justifyContent: "space-between",
    marginBottom: "8px",
    color: "#9fb8c7",
    fontSize: "12px"
  }}
>
  <span>Residue 1</span>
  <span>Residue {Math.round(size / 2)}</span>
  <span>Residue {size}</span>
</div>

<div
  style={{
    position: "relative",
    width: "100%",
    maxWidth: "550px"
  }}
>
  <div
    style={{
      position: "absolute",
      left: "-65px",
      top: "0",
      bottom: "0",
      display: "flex",
      flexDirection: "column",
      justifyContent: "space-between",
      color: "#9fb8c7",
      fontSize: "12px",
      textAlign: "right",
      width: "55px"
    }}
  >
    <span>1</span>
    <span>{Math.round(size / 2)}</span>
    <span>{size}</span>
  </div>

      <canvas
        ref={canvasRef}
        style={{
          display: "block",
          width: "100%",
          maxWidth: "550px",
          height: "auto",
          aspectRatio: "1",
          cursor: "crosshair"
        }}
        onMouseMove={(event) => {
          const bounds =
            event.currentTarget.getBoundingClientRect();

          const x = Math.min(
            size - 1,
            Math.max(
              0,
              Math.floor(
                ((event.clientX - bounds.left) /
                  bounds.width) * size
              )
            )
          );

          const y = Math.min(
            size - 1,
            Math.max(
              0,
              Math.floor(
                ((event.clientY - bounds.top) /
                  bounds.height) * size
              )
            )
          );

          setHover({
            x: x + 1,
            y: y + 1,
            error: matrix[y]?.[x]
          });
        }}
        onClick={(event) => {
  const bounds = event.currentTarget.getBoundingClientRect();

  const x = Math.min(
    size - 1,
    Math.max(0, Math.floor(((event.clientX - bounds.left) / bounds.width) * size))
  );

  const y = Math.min(
    size - 1,
    Math.max(0, Math.floor(((event.clientY - bounds.top) / bounds.height) * size))
  );

  setSelectedPair({
  x: x + 1,
  y: y + 1,
  error: matrix[y]?.[x]
});

onPairSelect?.({
  x: x + 1,
  y: y + 1,
  error: matrix[y]?.[x]
});

onReferenceChange?.(x + 1);
}}
        onMouseLeave={() => setHover(null)}
      />
</div>

      <div
  style={{
    maxWidth: "550px",
    display: "grid",
    gridTemplateColumns: "1fr",
    gap: "12px",
    marginTop: "16px"
  }}
>
  {[
    { title: "Selected Pair Information", pair: selectedPair },
    { title: "Hovering Pair Information", pair: hover }
  ].map(({ title, pair }) => (
    <div
      key={title}
      style={{
        padding: "18px",
        background: "#111d26",
        border: "1px solid #29414d",
        borderRadius: "12px",
        minHeight: "150px"
      }}
    >
      <h3 style={{ marginTop: 0 }}>{title}</h3>

      {pair ? (
        <>
          <p>
            <strong>Residue {pair.y}</strong> relative to{" "}
            <strong>residue {pair.x}</strong>
          </p>

          <p>
            Predicted aligned error:{" "}
            <strong>{Number(pair.error).toFixed(2)} Å</strong>
          </p>

          <p>
            {Number(pair.error) < 5
              ? "Low predicted error: relatively confident positioning."
              : Number(pair.error) < 15
              ? "Moderate predicted error: some positional uncertainty."
              : "High predicted error: substantial positional uncertainty."}
          </p>
        </>
      ) : (
        <p style={{ color: "#9fb8c7" }}>
          {title === "Selected Pair Information"
            ? "Click the heatmap to select a residue pair."
            : "Hover over the heatmap to inspect a residue pair."}
        </p>
      )}
    </div>
  ))}
</div>

<div style={{ maxWidth: "550px", marginTop: "20px" }}>
  <div
    style={{
      height: "16px",
      borderRadius: "5px",
      background:
        "linear-gradient(to right, rgb(20,100,210), rgb(138,153,125), rgb(255,205,40))",
      border: "1px solid #34444f"
    }}
  />

  <div
    style={{
      display: "flex",
      justifyContent: "space-between",
      marginTop: "8px",
      color: "#a9bdc9",
      fontSize: "12px"
    }}
  >
    <span>0 Å</span>
    <span>15 Å</span>
    <span>30+ Å</span>
  </div>

  <p style={{ color: "#a9bdc9", fontSize: "13px" }}>
    Predicted aligned error (Å). Lower values indicate greater confidence
    in the relative positioning of residue pairs.
  </p>
</div>

      <p>
        <strong>Blue:</strong> lower predicted error ·{" "}
        <strong>Yellow:</strong> higher predicted error
      </p>
    </div>
  );
}
