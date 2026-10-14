// Sky orientation: reading celestial WCS from FITS headers and building the
// 2x2 matrices that map BAYES-LOSVD (x, y) offsets to the displayed frame.
//
// BAYES-LOSVD instrument readers define x = (column - c0) * pixsize and
// y = (row - r0) * pixsize, i.e. offsets along the cube's pixel axes. The
// spatial CD matrix of the cube therefore maps (x, y) to sky offsets:
//
//   [ΔE]          [CD1_1 CD1_2] [x]        ΔE = offset towards East (ΔRA cos δ)
//   [ΔN] = 1/s ·  [CD2_1 CD2_2] [y]        ΔN = offset towards North
//
// with s = sqrt(|det CD|) so the matrix W is orthonormal (a rotation, possibly
// with a reflection).

export const ORIENTATIONS = ["data", "sky", "major"];

const BLOCK = 2880;

function parseCard(card) {
  const key = card.slice(0, 8).trim();
  if (card.slice(8, 10) !== "= ") return [key, undefined];
  let rest = card.slice(10).trim();
  if (rest.startsWith("'")) {
    let value = "";
    let i = 1;
    while (i < rest.length) {
      if (rest[i] === "'") {
        if (rest[i + 1] === "'") {
          value += "'";
          i += 2;
          continue;
        }
        break;
      }
      value += rest[i];
      i += 1;
    }
    return [key, value.trimEnd()];
  }
  const slash = rest.indexOf("/");
  if (slash >= 0) rest = rest.slice(0, slash).trim();
  if (rest === "T") return [key, true];
  if (rest === "F") return [key, false];
  const num = Number(rest.replace(/D/i, "E"));
  return [key, Number.isFinite(num) ? num : rest];
}

async function readHeader(blob, offset) {
  const cards = {};
  let pos = offset;
  for (let guard = 0; guard < 400; guard += 1) {
    if (pos + BLOCK > blob.size) return null;
    const text = await blob.slice(pos, pos + BLOCK).text();
    pos += BLOCK;
    for (let i = 0; i < BLOCK; i += 80) {
      const card = text.slice(i, i + 80);
      const key = card.slice(0, 8).trim();
      if (key === "END") return {cards, end: pos};
      if (!key) continue;
      const [k, value] = parseCard(card);
      if (value !== undefined && !(k in cards)) cards[k] = value;
    }
  }
  return null;
}

function dataSize(cards) {
  const naxis = Number(cards.NAXIS) || 0;
  if (!naxis) return 0;
  let n = 1;
  for (let i = 1; i <= naxis; i += 1) n *= Number(cards[`NAXIS${i}`]) || 0;
  const pcount = Number(cards.PCOUNT) || 0;
  const gcount = Number(cards.GCOUNT) || 1;
  const bytes = (Math.abs(Number(cards.BITPIX) || 8) / 8) * gcount * (pcount + n);
  return Math.ceil(bytes / BLOCK) * BLOCK;
}

// Spatial CD matrix (deg/pixel) from header cards, or null.
export function wcsFromCards(cards) {
  const axes = [];
  for (let i = 1; i <= 4; i += 1) {
    const ctype = String(cards[`CTYPE${i}`] || "").toUpperCase();
    if (ctype.startsWith("RA")) axes.push({axis: i, kind: "ra"});
    else if (ctype.startsWith("DEC")) axes.push({axis: i, kind: "dec"});
  }
  const ra = axes.find((a) => a.kind === "ra");
  const dec = axes.find((a) => a.kind === "dec");
  if (!ra || !dec) return null;
  const pix = [ra.axis, dec.axis].sort((a, b) => a - b); // pixel axes of the image plane
  const num = (key, fallback) => (Number.isFinite(Number(cards[key])) ? Number(cards[key]) : fallback);
  const cdAt = (world, pixel) => {
    if (`CD${world}_${pixel}` in cards) return num(`CD${world}_${pixel}`, 0);
    const hasPc = Object.keys(cards).some((k) => /^PC\d_\d$/.test(k) || /^PC0\d{2}0\d{2}$/.test(k));
    const cdelt = num(`CDELT${world}`, 1);
    if (hasPc) {
      const pc = num(`PC${world}_${pixel}`, num(`PC00${world}00${pixel}`, world === pixel ? 1 : 0));
      return cdelt * pc;
    }
    // Legacy CROTA2 convention
    const rho = num("CROTA2", num("CROTA1", 0)) * Math.PI / 180;
    if (world === pixel) return cdelt * Math.cos(rho);
    const sign = world === ra.axis ? -1 : 1;
    return sign * num(`CDELT${pixel}`, 1) * Math.sin(rho);
  };
  // Rows: East (RA axis), North (Dec axis). Columns: image pixel axes 1, 2.
  const cd = [
    cdAt(ra.axis, pix[0]), cdAt(ra.axis, pix[1]),
    cdAt(dec.axis, pix[0]), cdAt(dec.axis, pix[1])
  ];
  if (!cd.every(Number.isFinite) || Math.abs(cd[0] * cd[3] - cd[1] * cd[2]) === 0) return null;
  return {cd, pixelScaleArcsec: Math.sqrt(Math.abs(cd[0] * cd[3] - cd[1] * cd[2])) * 3600};
}

// Reads headers HDU by HDU (only the header blocks are read, so very large
// cubes are fine) and returns the first celestial WCS found.
export async function readFitsWcs(blob, name = "FITS file") {
  let offset = 0;
  for (let hdu = 0; hdu < 12 && offset < blob.size; hdu += 1) {
    const header = await readHeader(blob, offset);
    if (!header) break;
    if (hdu === 0 && header.cards.SIMPLE !== true) throw new Error(`${name} is not a FITS file`);
    const wcs = wcsFromCards(header.cards);
    if (wcs) return {...wcs, hdu, extname: header.cards.EXTNAME || null, source: name};
    offset = header.end + dataSize(header.cards);
  }
  throw new Error(`No celestial WCS (RA/DEC axes) found in the headers of ${name}`);
}

// Orthonormal raw (x, y) -> (ΔE, ΔN) matrix from a CD matrix.
export function skyMatrix(cd) {
  const s = Math.sqrt(Math.abs(cd[0] * cd[3] - cd[1] * cd[2]));
  return cd.map((v) => v / s);
}

// Position angle (deg, East of North) of the data +x and +y axes, and parity.
export function describeSkyMatrix(W) {
  const pa = (e, n) => ((Math.atan2(e, n) * 180 / Math.PI) + 360) % 360;
  return {
    paX: pa(W[0], W[2]),
    paY: pa(W[1], W[3]),
    standard: W[0] * W[3] - W[1] * W[2] < 0 // det < 0: x is West of y, like an image with East left
  };
}

export const mul = (A, B) => [
  A[0] * B[0] + A[1] * B[2], A[0] * B[1] + A[1] * B[3],
  A[2] * B[0] + A[3] * B[2], A[2] * B[1] + A[3] * B[3]
];

export const invert = (M) => {
  const det = M[0] * M[3] - M[1] * M[2];
  return [M[3] / det, -M[1] / det, -M[2] / det, M[0] / det];
};

export const apply = (M, x, y) => [M[0] * x + M[1] * y, M[2] * x + M[3] * y];

// (ΔE, ΔN) -> display frame for the sky-based orientations. In the major-axis
// orientation the end of the major axis that keeps North in the upper half of
// the map goes to +x (flip = true picks the other end).
export function skyToDisplay(orientation, galaxyPaDeg, flip = false) {
  if (orientation === "sky") return [-1, 0, 0, 1]; // East to the left, North up
  let pa = galaxyPaDeg;
  if (Math.sin(pa * Math.PI / 180) > 1e-9) pa += 180;
  if (flip) pa += 180;
  const t = pa * Math.PI / 180;
  const s = Math.sin(t);
  const c = Math.cos(t);
  // +x' along the major axis (PA t), a proper rotation of the North-up/East-left image
  return [s, c, c, -s];
}

// Full raw -> display matrix for one result. Returns {M, mode} where mode is
// the orientation actually used (falls back to "data" without WCS / PA).
export function displayMatrix({orientation, dataRotationDeg, galaxyPaDeg, flip}, W) {
  const t = (Number(dataRotationDeg) || 0) * Math.PI / 180;
  const data = [Math.cos(t), -Math.sin(t), Math.sin(t), Math.cos(t)];
  if (orientation === "data" || !W) return {M: data, mode: "data"};
  if (orientation === "major" && !Number.isFinite(galaxyPaDeg)) return {M: mul(skyToDisplay("sky"), W), mode: "sky"};
  return {M: mul(skyToDisplay(orientation, galaxyPaDeg, flip), W), mode: orientation};
}
