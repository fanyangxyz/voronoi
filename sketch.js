// Canvas
const WIDTH = 1000;
const HEIGHT = 1000;

// Interaction
const SCALE_STEP = 0.02;
const MIN_SCALE = 0.8;
const MAX_SCALE = 2.0;

// Geometry styling
const MIN_CELL_AREA = 9;
const PATCH_ALPHA = 0.02;

// Color clustering
const LIGHT_TONE_THRESHOLD = 0.72;
const DARK_TONE_THRESHOLD = 0.26;
const CLUSTER_GROW_PROBABILITY = 0.76;
const HUE_JITTER = 4;
const SAT_JITTER = 5;
const LIGHTNESS_JITTER = 5;

const LAYERS = [
  { count: 50000, clusterMin: 100, clusterMax: 240, accentRatio: 0.06 },
];

const PALETTE = {
  light: [
    [16, 65, 86],
    [22, 70, 82],
    [334, 58, 86],
    [45, 62, 84],
    [350, 42, 90],
  ],
  mid: [
    [12, 76, 66],
    [18, 78, 62],
    [322, 64, 66],
    [278, 44, 60],
    [300, 56, 56],
  ],
  dark: [
    [224, 44, 26],
    [238, 36, 22],
    [208, 52, 30],
    [318, 66, 34],
  ],
  accent: [
    [343, 90, 42],
    [210, 68, 38],
    [286, 58, 44],
  ],
};

// Runtime state
let cellScale = 30;
let composition = [];

// Initialize canvas and generate the first composition.
function setup() {
  const cnv = createCanvas(WIDTH, HEIGHT);
  cnv.parent(document.querySelector("main"));
  noLoop();
  generateComposition();
}

// Handle keyboard controls for regenerate and scale tuning.
function keyPressed() {
  if (key === "r" || key === "R") {
    generateComposition();
    redraw();
  } else if (key === "[" || key === "{") {
    cellScale = max(MIN_SCALE, cellScale - SCALE_STEP);
    redraw();
  } else if (key === "]" || key === "}") {
    cellScale = min(MAX_SCALE, cellScale + SCALE_STEP);
    redraw();
  }
}

// Render the full scene: patches, grain texture, and HUD frame.
function draw() {
  background(17, 22, 31);

  colorMode(HSL, 360, 100, 100, 1);
  for (const patch of composition) {
    drawPatch(patch);
  }
  colorMode(RGB, 255, 255, 255, 255);

  drawGrain();
  drawFrame();
}

// Build one full painterly composition from sampled Voronoi cells.
function generateComposition() {
  composition = [];
  for (let layerIndex = 0; layerIndex < LAYERS.length; layerIndex += 1) {
    const layer = LAYERS[layerIndex];
    const points = sampleRandomPoints(layer.count, width, height);

    const tess = buildVoronoi(points);
    const cells = extractCells(tess.voronoi, points);
    const neighbors = buildAdjacency(tess.delaunay, cells);
    // Neighbor-aware clustering assigns shared colors to adjacent cells.
    assignMergedColors(cells, neighbors, layer);

    for (const cell of cells) {
      composition.push({
        layerIndex,
        color: cell.color,
        center: cell.center,
        vertices: cell.vertices,
        area: cell.area,
      });
    }
  }

  composition.sort((a, b) => {
    if (a.layerIndex !== b.layerIndex) return a.layerIndex - b.layerIndex;
    return b.area - a.area;
  });
}

// Draw a single paint patch with scale inflation.
function drawPatch(patch) {
  const scaled = patch.vertices.map((v) => scaleAround(v, patch.center, cellScale));
  fill(patch.color.h, patch.color.s, patch.color.l, PATCH_ALPHA);
  noStroke();
  drawPolygon(scaled);
}

// Draw a closed polygon.
function drawPolygon(verts) {
  if (!verts || verts.length < 3) return;
  beginShape();
  for (const v of verts) {
    vertex(v.x, v.y);
  }
  endShape(CLOSE);
}

// Uniformly sample random points across the canvas.
function sampleRandomPoints(count, w, h) {
  const points = [];
  for (let i = 0; i < count; i += 1) {
    points.push({ x: random(w), y: random(h) });
  }
  return points;
}

// Build Delaunay triangulation and Voronoi diagram from point sites.
function buildVoronoi(points) {
  const coords = points.map((p) => [p.x, p.y]);
  const delaunay = d3.Delaunay.from(coords);
  return { delaunay, voronoi: delaunay.voronoi([0, 0, width, height]) };
}

// Convert Voronoi cells into renderable polygon records.
function extractCells(voronoi, points) {
  const cells = [];
  for (let i = 0; i < points.length; i += 1) {
    const poly = voronoi.cellPolygon(i);
    if (!poly || poly.length < 4) continue;

    const vertices = poly.slice(0, -1).map(([x, y]) => ({ x, y }));
    const center = polygonCentroid(poly);
    if (!center) continue;

    const area = abs(polygonArea(vertices));
    if (area < MIN_CELL_AREA) continue;

    cells.push({
      index: i,
      vertices,
      center,
      area,
      tone: "mid",
      color: { h: 0, s: 0, l: 0 },
    });
  }
  return cells;
}

// Build cell-neighbor graph using Delaunay connectivity.
function buildAdjacency(delaunay, cells) {
  const byIndex = new Map();
  for (let i = 0; i < cells.length; i += 1) {
    byIndex.set(cells[i].index, i);
  }

  const neighbors = Array.from({ length: cells.length }, () => []);
  for (let i = 0; i < cells.length; i += 1) {
    const src = cells[i].index;
    for (const n of delaunay.neighbors(src)) {
      const mapped = byIndex.get(n);
      if (mapped !== undefined) neighbors[i].push(mapped);
    }
  }
  return neighbors;
}

// Assign tones, then grow local color clusters across neighboring cells.
function assignMergedColors(cells, neighbors, layer) {
  for (const cell of cells) {
    // Noise-driven tone map establishes broad value structure.
    const n = noise(cell.center.x * 0.005, cell.center.y * 0.005);
    if (n > LIGHT_TONE_THRESHOLD) cell.tone = "light";
    else if (n < DARK_TONE_THRESHOLD) cell.tone = "dark";
    else cell.tone = "mid";

    if (random() < layer.accentRatio) cell.tone = "accent";
  }

  const unassigned = new Set();
  for (let i = 0; i < cells.length; i += 1) unassigned.add(i);

  while (unassigned.size > 0) {
    const start = pickFromSet(unassigned);
    const targetSize = floor(random(layer.clusterMin, layer.clusterMax + 1));
    const baseTone = cells[start].tone;
    const base = pickColor(baseTone);

    // BFS-like growth makes contiguous regions share the same color.
    const queue = [start];
    let painted = 0;

    while (queue.length > 0 && painted < targetSize) {
      const current = queue.shift();
      if (!unassigned.has(current)) continue;

      unassigned.delete(current);
      cells[current].color = base;
      painted += 1;

      const shuffled = shuffleArray(neighbors[current].slice());
      for (const nb of shuffled) {
        if (!unassigned.has(nb)) continue;
        if (!toneCompatible(baseTone, cells[nb].tone)) continue;
        if (random() < CLUSTER_GROW_PROBABILITY) queue.push(nb);
      }
    }
  }
}

// Restrict cluster growth to compatible tone families.
function toneCompatible(a, b) {
  if (a === b) return true;
  if ((a === "light" && b === "mid") || (a === "mid" && b === "light")) return true;
  if ((a === "dark" && b === "mid") || (a === "mid" && b === "dark")) return true;
  return false;
}

// Pick a palette color for a tone, then add slight jitter for variation.
function pickColor(tone) {
  const list = PALETTE[tone] || PALETTE.mid;
  const [h, s, l] = random(list);
  return {
    h: (h + random(-HUE_JITTER, HUE_JITTER) + 360) % 360,
    s: constrain(s + random(-SAT_JITTER, SAT_JITTER), 0, 100),
    l: constrain(l + random(-LIGHTNESS_JITTER, LIGHTNESS_JITTER), 0, 100),
  };
}

// Compute centroid of a closed polygon (shoelace-based).
function polygonCentroid(poly) {
  let twiceArea = 0;
  let cx = 0;
  let cy = 0;

  for (let i = 0; i < poly.length - 1; i += 1) {
    const [x0, y0] = poly[i];
    const [x1, y1] = poly[i + 1];
    const cross = x0 * y1 - x1 * y0;
    twiceArea += cross;
    cx += (x0 + x1) * cross;
    cy += (y0 + y1) * cross;
  }

  if (abs(twiceArea) < 1e-9) return null;
  return { x: cx / (3 * twiceArea), y: cy / (3 * twiceArea) };
}

// Compute signed polygon area (positive or negative by winding).
function polygonArea(verts) {
  let area = 0;
  for (let i = 0; i < verts.length; i += 1) {
    const a = verts[i];
    const b = verts[(i + 1) % verts.length];
    area += a.x * b.y - b.x * a.y;
  }
  return area * 0.5;
}

// Uniformly scale a point around an origin.
function scaleAround(p, origin, scale) {
  return {
    x: origin.x + (p.x - origin.x) * scale,
    y: origin.y + (p.y - origin.y) * scale,
  };
}

// Pick one random element from a Set.
function pickFromSet(setObj) {
  const idx = floor(random(setObj.size));
  let i = 0;
  for (const item of setObj) {
    if (i === idx) return item;
    i += 1;
  }
  return 0;
}

// In-place Fisher-Yates array shuffle.
function shuffleArray(arr) {
  for (let i = arr.length - 1; i > 0; i -= 1) {
    const j = floor(random(i + 1));
    const tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
  }
  return arr;
}

// Overlay subtle bright/dark grain for canvas-like texture.
function drawGrain() {
  noStroke();
  for (let i = 0; i < 9000; i += 1) {
    const x = random(width);
    const y = random(height);
    const a = random(10, 24);
    fill(255, 255, 255, a);
    rect(x, y, 1, 1);
  }

  for (let i = 0; i < 6500; i += 1) {
    const x = random(width);
    const y = random(height);
    const a = random(8, 18);
    fill(0, 0, 0, a);
    rect(x, y, 1, 1);
  }
}

// Draw border and controls/status text.
function drawFrame() {
  noFill();
  stroke(255, 255, 255, 20);
  strokeWeight(2);
  rect(9, 9, width - 18, height - 18);

  noStroke();
  fill(255, 220);
  textSize(14);
  textAlign(LEFT, TOP);
  text(`Painterly Voronoi | Scale: ${cellScale.toFixed(2)} ([ / ]) | R regenerate`, 18, 16);
}
