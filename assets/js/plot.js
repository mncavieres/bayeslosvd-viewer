// Small canvas plotting toolkit: HiDPI setup, nice ticks, axes, colorbars.

export function setupCanvas(canvas) {
  const rect = canvas.getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;
  const w = Math.max(1, Math.round(rect.width * dpr));
  const h = Math.max(1, Math.round(rect.height * dpr));
  if (canvas.width !== w) canvas.width = w;
  if (canvas.height !== h) canvas.height = h;
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return {ctx, width: rect.width, height: rect.height, dpr};
}

export function theme() {
  const css = getComputedStyle(document.documentElement);
  const v = (name) => css.getPropertyValue(name).trim();
  return {
    plotBg: v("--plot-bg"),
    ink: v("--ink"),
    muted: v("--muted"),
    axis: v("--axis"),
    grid: v("--grid-line"),
    accent: v("--accent"),
    empty: v("--empty-cell"),
    font: getComputedStyle(document.body).fontFamily
  };
}

export function niceStep(range, target = 5) {
  if (!(range > 0)) return 1;
  const raw = range / target;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const f = raw / mag;
  const nice = f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10;
  return nice * mag;
}

export function niceTicks(lo, hi, target = 5) {
  const step = niceStep(hi - lo, target);
  const start = Math.ceil(lo / step - 1e-9) * step;
  const ticks = [];
  for (let v = start; v <= hi + step * 1e-9; v += step) ticks.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  return {ticks, step};
}

export function formatTick(v, step) {
  if (!Number.isFinite(v)) return "";
  const decimals = Math.max(0, Math.min(6, -Math.floor(Math.log10(step) + 1e-9)));
  const abs = Math.abs(v);
  if ((abs >= 1e5 || (abs > 0 && abs < 1e-3)) && step < abs) return v.toExponential(1);
  return v.toFixed(decimals);
}

export function formatValue(v, unit = "") {
  if (!Number.isFinite(v)) return "–";
  // Treat round-off noise (e.g. 1e-14 km/s for a symmetric LOSVD) as zero.
  const abs = Math.abs(v) < 1e-9 ? 0 : Math.abs(v);
  if (abs === 0) v = 0;
  let text;
  if (abs !== 0 && (abs >= 1e5 || abs < 1e-3)) text = v.toExponential(2);
  else text = v.toFixed(abs >= 100 ? 1 : abs >= 10 ? 2 : 3);
  return unit ? `${text} ${unit}` : text;
}

// Draw a framed axes box with ticks. plot = {left, top, width, height, xmin, xmax, ymin, ymax}
export function drawAxes(ctx, plot, t, {xlabel = "", ylabel = "", xTicks = 6, yTicks = 5, gridLines = false, xFormat = null, yFormat = null} = {}) {
  ctx.save();
  ctx.lineWidth = 1;
  ctx.font = `11px ${t.font}`;
  const xs = (v) => plot.left + (v - plot.xmin) / (plot.xmax - plot.xmin) * plot.width;
  const ys = (v) => plot.top + (plot.ymax - v) / (plot.ymax - plot.ymin) * plot.height;
  const xt = niceTicks(plot.xmin, plot.xmax, xTicks);
  const yt = niceTicks(plot.ymin, plot.ymax, yTicks);
  if (gridLines) {
    ctx.strokeStyle = t.grid;
    ctx.beginPath();
    xt.ticks.forEach((v) => {
      const x = Math.round(xs(v)) + 0.5;
      ctx.moveTo(x, plot.top);
      ctx.lineTo(x, plot.top + plot.height);
    });
    yt.ticks.forEach((v) => {
      const y = Math.round(ys(v)) + 0.5;
      ctx.moveTo(plot.left, y);
      ctx.lineTo(plot.left + plot.width, y);
    });
    ctx.stroke();
  }
  ctx.strokeStyle = t.axis;
  ctx.strokeRect(Math.round(plot.left) + 0.5, Math.round(plot.top) + 0.5, Math.round(plot.width), Math.round(plot.height));
  ctx.fillStyle = t.muted;
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  ctx.beginPath();
  xt.ticks.forEach((v) => {
    const x = Math.round(xs(v)) + 0.5;
    ctx.moveTo(x, plot.top + plot.height);
    ctx.lineTo(x, plot.top + plot.height + 4);
    ctx.fillText(xFormat ? xFormat(v, xt.step) : formatTick(v, xt.step), x, plot.top + plot.height + 6);
  });
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  yt.ticks.forEach((v) => {
    const y = Math.round(ys(v)) + 0.5;
    ctx.moveTo(plot.left - 4, y);
    ctx.lineTo(plot.left, y);
    ctx.fillText(yFormat ? yFormat(v, yt.step) : formatTick(v, yt.step), plot.left - 6, y);
  });
  ctx.stroke();
  ctx.fillStyle = t.ink;
  ctx.font = `12px ${t.font}`;
  if (xlabel) {
    ctx.textAlign = "center";
    ctx.textBaseline = "bottom";
    ctx.fillText(xlabel, plot.left + plot.width / 2, plot.top + plot.height + 38);
  }
  if (ylabel) {
    ctx.save();
    ctx.translate(plot.left - 44, plot.top + plot.height / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.fillText(ylabel, 0, 0);
    ctx.restore();
  }
  ctx.restore();
}

export function drawColorbar(canvas, lut, vmin, vmax, label, t, {top = 14, bottom = 44} = {}) {
  const {ctx, width, height} = setupCanvas(canvas);
  ctx.clearRect(0, 0, width, height);
  const left = 6;
  const barW = 14;
  const barH = Math.max(10, height - top - bottom);
  const img = ctx.createImageData(1, 256);
  for (let i = 0; i < 256; i += 1) {
    const src = (255 - i) * 3;
    img.data[i * 4] = lut[src];
    img.data[i * 4 + 1] = lut[src + 1];
    img.data[i * 4 + 2] = lut[src + 2];
    img.data[i * 4 + 3] = 255;
  }
  const off = document.createElement("canvas");
  off.width = 1;
  off.height = 256;
  off.getContext("2d").putImageData(img, 0, 0);
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(off, left, top, barW, barH);
  ctx.strokeStyle = t.axis;
  ctx.lineWidth = 1;
  ctx.strokeRect(left + 0.5, top + 0.5, barW, barH);
  const {ticks, step} = niceTicks(vmin, vmax, 6);
  ctx.fillStyle = t.muted;
  ctx.font = `11px ${t.font}`;
  ctx.textAlign = "left";
  ctx.textBaseline = "middle";
  ctx.beginPath();
  ticks.forEach((v) => {
    const y = Math.round(top + barH * (1 - (v - vmin) / (vmax - vmin))) + 0.5;
    ctx.moveTo(left + barW, y);
    ctx.lineTo(left + barW + 4, y);
    ctx.fillText(formatTick(v, step), left + barW + 6, y);
  });
  ctx.stroke();
  if (label) {
    ctx.save();
    ctx.translate(width - 4, top + barH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.textBaseline = "bottom";
    ctx.fillStyle = t.ink;
    ctx.font = `12px ${t.font}`;
    ctx.fillText(label, 0, 0);
    ctx.restore();
  }
}
