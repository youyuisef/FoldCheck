
export default function RegionDetails({
  region,
  proteinInfo,
  partnerInfo = {}
}) {
  if (!region) return null;

  const proteinName = proteinInfo?.gene || "this protein";

  const partnerName = region.name
    ?.match(/^Interaction with\s+(.+)$/i)?.[1]
    ?.trim();

  const matchingInteraction = proteinInfo?.interactionPartners?.find(
  (item) =>
    partnerName &&
    item.name?.toLowerCase() === partnerName.toLowerCase()
);

const partner = partnerName
  ? (
      (matchingInteraction?.id
        ? partnerInfo[matchingInteraction.id]
        : null) ||
      partnerInfo[partnerName.toUpperCase()] ||
      Object.values(partnerInfo).find(
        (item) =>
          item?.gene?.toLowerCase() === partnerName.toLowerCase() ||
          item?.name?.toLowerCase() === partnerName.toLowerCase()
      )
    )
  : null;

const experimentCount = matchingInteraction?.experiments;

const intactId = matchingInteraction?.intActId;

  const names = [
    proteinName,
    ...(proteinInfo?.aliases || [])
  ].filter(Boolean);

  const relevantSentences = (partner?.function || "")
    .split(/(?<=[.!?])\s+/)
    .filter((sentence) =>
      names.some(
        (name) =>
          name.length > 2 &&
          sentence.toLowerCase().includes(name.toLowerCase())
      )
    );

  return (
    <div className="regionDetails">
        {(experimentCount > 0 || intactId) && (
  <section>
    <h4>Interaction evidence</h4>

    {experimentCount > 0 && (
      <p>
        {experimentCount} experiment{experimentCount === 1 ? "" : "s"} reported
        for the interaction between {proteinName} and {partnerName}.
      </p>
    )}

    {intactId && (
      <p>
        IntAct identifier: <strong>{intactId}</strong>
      </p>
    )}
  </section>
)}

{partner && (
  <section>
    <h4>Interaction partner</h4>

    <p>
      <strong>{partner.name || partnerName}</strong>
      {partner.gene &&
        partner.gene.toLowerCase() !==
          (partner.name || "").toLowerCase() && (
          <> · {partner.gene}</>
        )}
    </p>

    {partner.id && (
      <p>
        UniProt accession: <strong>{partner.id}</strong>
      </p>
    )}
  </section>
)}

      {relevantSentences.length > 0 && (
        <section>
          <h4>Biological findings</h4>
          <ul>
            {relevantSentences.map((sentence, index) => (
              <li key={index}>{sentence.trim()}</li>
            ))}
          </ul>
        </section>
      )}

      {partner?.id && relevantSentences.length > 0 && (
        <a
          href={`https://www.uniprot.org/uniprotkb/${encodeURIComponent(partner.id)}/entry`}
          target="_blank"
          rel="noopener noreferrer"
        >
          Explore {partnerName} ↗
        </a>
      )}
    </div>
  );
}
