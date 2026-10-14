// Colormaps as evenly spaced hex stops, expanded into 256-entry lookup tables.
// Stops are sampled from the matplotlib definitions (and plotbin's "sauron").

const STOPS = {
  RdBu_r: ["#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"],
  coolwarm: ["#3b4cc0", "#5977e3", "#7b9ff9", "#9ebeff", "#c0d4f5", "#dddcdc", "#f2cbb7", "#f7ac8e", "#ee8468", "#d65244", "#b40426"],
  RdYlBu_r: ["#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#ffffbf", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"],
  viridis: ["#440154", "#482878", "#3e4989", "#31688e", "#26828e", "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725"],
  magma: ["#000004", "#180f3d", "#440f76", "#721f81", "#9e2f7f", "#cd4071", "#f1605d", "#fd9668", "#feca8d", "#fcfdbf"],
  inferno: ["#000004", "#1b0c41", "#4a0c6b", "#781c6d", "#a52c60", "#cf4446", "#ed6925", "#fb9b06", "#f7d13d", "#fcffa4"],
  cividis: ["#00224e", "#123570", "#3b496c", "#575d6d", "#707173", "#8a8678", "#a59c74", "#c3b369", "#e1cc55", "#fee838"],
  gray: ["#000000", "#ffffff"]
};

// plotbin's sauron colormap (Cappellari), defined on an uneven grid.
const SAURON = {
  x: [1.0, 43.5, 86.0, 106.0, 118.5, 128.5, 138.5, 151.0, 171.0, 213.5, 256.0],
  r: [0.01, 0.0, 0.4, 0.5, 0.3, 0.0, 0.7, 1.0, 1.0, 1.0, 0.9],
  g: [0.01, 0.0, 0.85, 1.0, 1.0, 0.9, 1.0, 1.0, 0.85, 0.0, 0.9],
  b: [0.01, 1.0, 1.0, 1.0, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9]
};

export const COLORMAP_NAMES = ["RdBu_r", "coolwarm", "RdYlBu_r", "sauron", "viridis", "magma", "inferno", "cividis", "gray"];

const cache = new Map();

function hexToRgb(hex) {
  const s = hex.replace("#", "");
  return [parseInt(s.slice(0, 2), 16), parseInt(s.slice(2, 4), 16), parseInt(s.slice(4, 6), 16)];
}

function buildFromStops(stops) {
  const rgb = stops.map(hexToRgb);
  const lut = new Uint8ClampedArray(256 * 3);
  const segments = rgb.length - 1;
  for (let i = 0; i < 256; i += 1) {
    const t = (i / 255) * segments;
    const k = Math.min(segments - 1, Math.floor(t));
    const f = t - k;
    for (let c = 0; c < 3; c += 1) {
      lut[i * 3 + c] = Math.round(rgb[k][c] + f * (rgb[k + 1][c] - rgb[k][c]));
    }
  }
  return lut;
}

function buildSauron() {
  const lut = new Uint8ClampedArray(256 * 3);
  const channels = [SAURON.r, SAURON.g, SAURON.b];
  for (let i = 0; i < 256; i += 1) {
    const xv = 1 + (i / 255) * 255;
    let k = 0;
    while (k < SAURON.x.length - 2 && xv > SAURON.x[k + 1]) k += 1;
    const f = (xv - SAURON.x[k]) / (SAURON.x[k + 1] - SAURON.x[k]);
    channels.forEach((ch, c) => {
      lut[i * 3 + c] = Math.round(255 * (ch[k] + f * (ch[k + 1] - ch[k])));
    });
  }
  return lut;
}

export function getLut(name, reversed = false) {
  const key = `${name}|${reversed}`;
  if (cache.has(key)) return cache.get(key);
  let lut = name === "sauron" ? buildSauron() : buildFromStops(STOPS[name] || STOPS.viridis);
  if (reversed) {
    const rev = new Uint8ClampedArray(lut.length);
    for (let i = 0; i < 256; i += 1) {
      rev.set(lut.subarray((255 - i) * 3, (255 - i) * 3 + 3), i * 3);
    }
    lut = rev;
  }
  cache.set(key, lut);
  return lut;
}

export function lutIndex(value, vmin, vmax) {
  const t = (value - vmin) / (vmax - vmin || 1);
  return Math.max(0, Math.min(255, Math.round(t * 255)));
}

export function cssColor(lut, value, vmin, vmax) {
  const i = lutIndex(value, vmin, vmax) * 3;
  return `rgb(${lut[i]},${lut[i + 1]},${lut[i + 2]})`;
}
