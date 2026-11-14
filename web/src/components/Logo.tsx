/** mzDesk logo, "kneehole cube": an isometric cube whose front face is carved with two kneeholes (a desk, and the three stems of an m),
 *  a z on the side face, and the top face as the desk surface. Colours come from CSS (.logo .lg-*), so one shape serves the light page,
 *  the dark page and the blue sidebar. Geometry generated once (isometric 30°); public/favicon.svg is the same drawing. */
const PARTS: [string, string][] = [
  ["back", "23.88,36.80 30.81,40.80 30.81,64.00 23.88,60.00"],
  ["wall", "9.33,45.20 23.88,36.80 23.88,60.00 9.33,68.40"],
  ["ceil", "9.33,45.20 16.26,49.20 30.81,40.80 23.88,36.80"],
  ["back", "37.73,44.80 44.66,48.80 44.66,72.00 37.73,68.00"],
  ["wall", "23.18,53.20 37.73,44.80 37.73,68.00 23.18,76.40"],
  ["ceil", "23.18,53.20 30.11,57.20 44.66,48.80 37.73,44.80"],
  ["front", "2.40,22.40 37.04,42.40 37.04,84.40 30.11,80.40 30.11,57.20 23.18,53.20 23.18,76.40 16.26,72.40 16.26,49.20 9.33,45.20 9.33,68.40 2.40,64.40"],
  ["side", "37.04,42.40 71.68,22.40 71.68,64.40 37.04,84.40"],
  ["top", "37.04,2.40 71.68,22.40 37.04,42.40 2.40,22.40"],
  ["letter", "45.35,46.80 63.37,36.40 63.37,41.20 51.24,61.80 63.37,54.80 63.37,59.60 45.35,70.00 45.35,65.20 57.48,44.60 45.35,51.60"]
];

export function Logo({ className = "" }: { className?: string }) {
  return (
    <svg className={"logo " + className} viewBox="0 0 74.08 86.80" aria-hidden="true" focusable="false">
      {PARTS.map(([cls, points], i) => <polygon key={i} className={"lg-" + cls} points={points} />)}
    </svg>
  );
}
