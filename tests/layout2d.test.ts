import { describe, it, expect } from 'vitest';
import { Mol } from '../src/chem/mol';
import { parseSmiles, writeSmiles } from '../src/chem/smiles';
import { perceiveRings } from '../src/chem/rings';
import { perceiveStereo2D } from '../src/chem/stereo2d';
import { layoutMol } from '../src/chem/layout2d';
import { clean2D } from '../src/chem/clean2d';
import { TEMPLATE_GROUPS, templateMol, allTemplates } from '../src/chem/templates';

// Geometry thresholds (bond length unit = 1):
//  • ordinary molecules (chains, fused/spiro rings, macrocycles): every bond within 0.9–1.1;
//  • bridged/cage ring systems (two SSSR rings sharing ≥ 3 atoms – norbornane, adamantane, cubane,
//    morphine, strychnine, taxanes, porphyrins…): 0.6–1.45, since perspective/strained drawings
//    cannot keep all bonds at unit length;
//  • no two non-bonded atoms closer than 0.4, for every molecule.
const NORMAL = [0.9, 1.1];
const BRIDGED = [0.6, 1.45];
const MIN_NONBONDED = 0.4;

const TRICKY: [string, string][] = [
  ['cholesterol', 'C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C'],
  ['morphine', 'CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O'],
  ['strychnine', 'C1CN2CC3=CCO[C@H]4CC(=O)N5[C@H]6[C@H]4[C@H]3C[C@H]2[C@@]61C7=CC=CC=C75'],
  ['paclitaxel', 'CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C'],
  ['taxadiene', 'CC1=C2CC[C@@]3(CCC=C([C@H]3C[C@@H](C2(C)C)CC1)C)C'],
  ['cubane', 'C12C3C4C1C5C2C3C45'],
  ['adamantane', 'C1C2CC3CC1CC(C2)C3'],
  ['[2.2]paracyclophane', 'C1CC2=CC=C(CCC3=CC=C1C=C3)C=C2'],
  ['porphine', 'C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2'],
  ['cyclododecane', 'C1CCCCCCCCCCC1'],
  ['beta-carotene', 'CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C2=C(CCCC2(C)C)C)\\C)\\C)/C)/C'],
  ['triacontane', 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'],
  ['quinine', 'COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O'],
  ['cocaine', 'CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)c1ccccc1)C(=O)OC'],
  ['alpha-pinene', 'CC1=CCC2CC1C2(C)C'],
  ['artemisinin', 'C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C'],
  ['erythromycin', 'CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O'],
  ['civetone (Z, 17-ring)', 'C1CCC/C=C\\CCCCCCCC(=O)CCC1'],
  ['(E)-cyclododecene', 'C1CCCCC/C=C/CCCC1'],
  ['alpha-cyclodextrin', 'C([C@@H]1[C@@H]2[C@@H]([C@H]([C@H](O1)O[C@@H]3[C@H](O[C@@H]([C@@H]([C@H]3O)O)O[C@@H]4[C@H](O[C@@H]([C@@H]([C@H]4O)O)O[C@@H]5[C@H](O[C@@H]([C@@H]([C@H]5O)O)O[C@@H]6[C@H](O[C@@H]([C@@H]([C@H]6O)O)O[C@@H]7[C@H](O[C@H](O2)[C@@H]([C@H]7O)O)CO)CO)CO)CO)CO)O)O)O'],
  ['reserpine', 'CO[C@H]1[C@@H](C[C@@H]2CN3CCC4=C([C@H]3C[C@@H]2[C@@H]1C(=O)OC)NC5=C4C=CC(=C5)OC)OC(=O)C6=CC(=C(C(=C6)OC)OC)OC'],
  ['pyrene', 'c1cc2ccc3cccc4ccc(c1)c2c34'],
  ['spiro[4.5]decane', 'C1CCC2(CC1)CCCC2'],
  ['allene + alkyne', 'CC=C=CC#CC'],
  ['phosphate/sulfonyl', 'CCOP(=O)(OCC)OCC.CS(=O)(=O)N(C)C'],
];

interface Geo { minBond: number; maxBond: number; minNB: number }

function geometry(m: Mol): Geo {
  let minBond = Infinity, maxBond = 0, minNB = Infinity;
  for (const b of m.bonds) {
    const d = Math.hypot(m.atoms[b.a].x - m.atoms[b.b].x, m.atoms[b.a].y - m.atoms[b.b].y);
    minBond = Math.min(minBond, d);
    maxBond = Math.max(maxBond, d);
  }
  for (let i = 0; i < m.atoms.length; i++) {
    for (let j = i + 1; j < m.atoms.length; j++) {
      if (m.bondBetween(i, j) >= 0) continue;
      minNB = Math.min(minNB, Math.hypot(m.atoms[i].x - m.atoms[j].x, m.atoms[i].y - m.atoms[j].y));
    }
  }
  return { minBond, maxBond, minNB };
}

/**
 * True for bridged/cage ring systems: two SSSR rings sharing ≥ 3 atoms, or an atom completely
 * surrounded by rings whose regular-polygon angles cannot tile the plane (cubane: 3 × 90°).
 */
function isBridged(m: Mol): boolean {
  const info = perceiveRings(m);
  for (const sys of info.systems) {
    for (let p = 0; p < sys.length; p++) {
      for (let q = p + 1; q < sys.length; q++) {
        const a = new Set(info.rings[sys[p]]);
        if (info.rings[sys[q]].filter((x) => a.has(x)).length >= 3) return true;
      }
    }
  }
  for (let i = 0; i < m.atoms.length; i++) {
    const rings = info.atomRings[i];
    if (rings.length < 3) continue;
    const ringNbrs = m.neighbors(i).filter((j) => rings.some((r) => info.rings[r].includes(j)));
    if (ringNbrs.length !== rings.length) continue;
    const sum = rings.reduce((s, r) => s + 180 - 360 / info.rings[r].length, 0);
    if (Math.abs(sum - 360) > 15) return true;
  }
  return false;
}

/** SMILES re-derived from the drawing (wedges + coordinates), restricted to the specified stereo. */
function drawnSmiles(m: Mol, centres: Set<number>, dbBonds: Set<number>): string {
  const c = m.clone();
  perceiveStereo2D(c);
  c.tetra = c.tetra.filter((t) => centres.has(t.center));
  c.dbStereo = c.dbStereo.filter((d) => dbBonds.has(d.bond));
  return writeSmiles(c);
}

function checkLayout(name: string, smiles: string, m: Mol): void {
  const g = geometry(m);
  const [lo, hi] = isBridged(m) ? BRIDGED : NORMAL;
  if (m.bonds.length) {
    expect(g.minBond, `${name}: shortest bond`).toBeGreaterThanOrEqual(lo);
    expect(g.maxBond, `${name}: longest bond`).toBeLessThanOrEqual(hi);
  }
  if (m.atoms.length > 1) expect(g.minNB, `${name}: closest non-bonded pair`).toBeGreaterThanOrEqual(MIN_NONBONDED);
  for (const a of m.atoms) {
    expect(Number.isFinite(a.x) && Number.isFinite(a.y), `${name}: finite coordinates`).toBe(true);
  }
  void smiles;
}

function stereoRoundTrip(name: string, smiles: string): void {
  const m = parseSmiles(smiles);
  const ref = writeSmiles(m);
  const centres = new Set(m.tetra.map((t) => t.center));
  const dbBonds = new Set(m.dbStereo.map((d) => d.bond));
  layoutMol(m);
  expect(drawnSmiles(m, centres, dbBonds), `${name}: stereo after layout`).toBe(ref);
}

/** Seeded random generator for reproducible perturbations. */
function rng(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

describe('layoutMol: templates', () => {
  for (const group of TEMPLATE_GROUPS) {
    it(`lays out every "${group.name}" template cleanly with correct stereo`, () => {
      for (const item of group.items) {
        const m = templateMol(item);
        checkLayout(item.name, item.smiles, m);
        stereoRoundTrip(item.name, item.smiles);
      }
    });
  }

  it('template SMILES all parse and the library is reasonably complete', () => {
    const items = allTemplates();
    expect(items.length).toBeGreaterThan(150);
    const cosmetic = TEMPLATE_GROUPS.find((g) => g.name.startsWith('Cosmetic'))!;
    expect(cosmetic.items.length).toBeGreaterThanOrEqual(40);
    for (const it2 of items) expect(() => parseSmiles(it2.smiles), it2.name).not.toThrow();
    expect(TEMPLATE_GROUPS.find((g) => g.name === 'Amino acids')!.items.length).toBe(20);
  });

  it('templateMol returns independent copies', () => {
    const item = TEMPLATE_GROUPS[0].items[6]; // benzene
    const a = templateMol(item);
    a.atoms[0].x += 10;
    const b = templateMol(item);
    expect(b.atoms[0].x).not.toBe(a.atoms[0].x);
  });
});

describe('layoutMol: tricky molecules', () => {
  for (const [name, smiles] of TRICKY) {
    it(name, () => {
      const m = parseSmiles(smiles);
      layoutMol(m);
      checkLayout(name, smiles, m);
      stereoRoundTrip(name, smiles);
    });
  }
});

describe('layoutMol: geometry conventions', () => {
  const angleAt = (m: Mol, a: number, c: number, b: number) => {
    const A = m.atoms[a], C = m.atoms[c], B = m.atoms[b];
    const v1 = [A.x - C.x, A.y - C.y], v2 = [B.x - C.x, B.y - C.y];
    const cos = (v1[0] * v2[0] + v1[1] * v2[1]) / (Math.hypot(v1[0], v1[1]) * Math.hypot(v2[0], v2[1]));
    return (Math.acos(Math.max(-1, Math.min(1, cos))) * 180) / Math.PI;
  };

  it('draws chains as 120° trans zig-zags with the long axis horizontal', () => {
    const m = parseSmiles('CCCCCCCCCC');
    layoutMol(m);
    for (let i = 1; i < 9; i++) expect(angleAt(m, i - 1, i, i + 1)).toBeCloseTo(120, 1);
    // trans: atoms i and i+3 far apart (2.65 for a zig-zag)
    for (let i = 0; i + 3 < 10; i++) {
      const d = Math.hypot(m.atoms[i].x - m.atoms[i + 3].x, m.atoms[i].y - m.atoms[i + 3].y);
      expect(d).toBeCloseTo(Math.sqrt(7), 2);
    }
    const bb = m.bbox();
    expect(bb.maxX - bb.minX).toBeGreaterThan(3 * (bb.maxY - bb.minY));
  });

  it('keeps triple bonds and cumulenes linear', () => {
    const m = parseSmiles('CC#CC');
    layoutMol(m);
    expect(angleAt(m, 0, 1, 2)).toBeCloseTo(180, 3);
    expect(angleAt(m, 1, 2, 3)).toBeCloseTo(180, 3);
    const a = parseSmiles('CC=C=C=CC');
    layoutMol(a);
    for (const i of [2, 3]) expect(angleAt(a, i - 1, i, i + 1)).toBeCloseTo(180, 3);
  });

  it('draws rings as regular polygons with radial substituents', () => {
    const m = parseSmiles('Cc1ccccc1');
    layoutMol(m);
    for (let k = 1; k <= 6; k++) {
      const a = k, b = k === 6 ? 1 : k + 1;
      const d = Math.hypot(m.atoms[a].x - m.atoms[b].x, m.atoms[a].y - m.atoms[b].y);
      expect(d).toBeCloseTo(1, 6);
    }
    expect(angleAt(m, 0, 1, 2)).toBeCloseTo(120, 3);
    expect(angleAt(m, 0, 1, 6)).toBeCloseTo(120, 3);
  });

  it('isolated benzene has vertical side bonds', () => {
    const m = parseSmiles('c1ccccc1');
    layoutMol(m);
    const vertical = m.bonds.filter((b) => Math.abs(m.atoms[b.a].x - m.atoms[b.b].x) < 1e-6).length;
    expect(vertical).toBe(2);
  });

  it('uses the cross layout for 4-coordinate chain atoms', () => {
    const m = parseSmiles('CC(C)(C)C');
    layoutMol(m);
    for (const j of [2, 3, 4]) {
      const ang = angleAt(m, 0, 1, j);
      expect([90, 180].some((t) => Math.abs(ang - t) < 1)).toBe(true);
    }
  });

  it('honours E/Z specifications', () => {
    for (const s of ['C/C=C/C', 'C/C=C\\C', 'CC/C=C(/C)CC', 'OC/C=C\\C=C\\C=C/CO', 'C1CCCCC/C=C/CCCC1', 'C1CCCCC/C=C\\CCCC1']) {
      stereoRoundTrip(s, s);
    }
  });

  it('is deterministic', () => {
    const s = 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O';
    const a = parseSmiles(s), b = parseSmiles(s);
    layoutMol(a);
    layoutMol(b);
    expect(a.atoms.map((p) => [p.x, p.y])).toEqual(b.atoms.map((p) => [p.x, p.y]));
  });

  it('places explicit hydrogens like other substituents', () => {
    const m = parseSmiles('[H]C([H])([H])C([H])([H])O[H]');
    layoutMol(m);
    checkLayout('ethanol with H', '', m);
  });

  it('arranges components left to right, counter-ion next to its partner', () => {
    const m = parseSmiles('CCO.c1ccccc1.O');
    layoutMol(m);
    const comps = m.components();
    const boxes = comps.map((c) => ({
      minX: Math.min(...c.map((a) => m.atoms[a].x)),
      maxX: Math.max(...c.map((a) => m.atoms[a].x)),
      cy: (Math.min(...c.map((a) => m.atoms[a].y)) + Math.max(...c.map((a) => m.atoms[a].y))) / 2,
    }));
    boxes.sort((p, q) => p.minX - q.minX);
    for (let k = 1; k < boxes.length; k++) {
      expect(boxes[k].minX - boxes[k - 1].maxX).toBeCloseTo(1.5, 6);
      expect(boxes[k].cy).toBeCloseTo(boxes[0].cy, 6);
    }
    const salt = parseSmiles('CCCCCCCCCCCC(=O)[O-].[Na+]');
    layoutMol(salt);
    const o = salt.atoms.findIndex((a) => a.charge === -1);
    const na = salt.atoms.findIndex((a) => a.el === 'Na');
    const d = Math.hypot(salt.atoms[o].x - salt.atoms[na].x, salt.atoms[o].y - salt.atoms[na].y);
    expect(d).toBeLessThan(2);
    expect(geometry(salt).minNB).toBeGreaterThan(1.2);
  });

  it('draws wedges that reproduce every stereocentre', () => {
    const m = parseSmiles('C[C@@H](O)[C@H](N)C(=O)O');
    layoutMol(m);
    expect(m.bonds.filter((b) => b.style === 'wedge' || b.style === 'hash').length).toBe(2);
    stereoRoundTrip('threonine-like', 'C[C@@H](O)[C@H](N)C(=O)O');
  });
});

describe('layoutMol: fixed atoms', () => {
  it('keeps fixed coordinates and places new atoms around them', () => {
    const m = parseSmiles('c1ccccc1C(=O)O');
    layoutMol(m);
    const before = m.atoms.map((a) => ({ x: a.x, y: a.y }));
    const fixed = new Set(m.atoms.map((_, i) => i));
    // attach an ethyl ester and a methoxy group (new, unplaced atoms)
    const o = m.atoms.findIndex((a, i) => a.el === 'O' && m.degree(i) === 1 && m.bonds[m.adj[i][0]].order === 1);
    const c1 = m.addAtom({ el: 'C' });
    const c2 = m.addAtom({ el: 'C' });
    m.addBond(o, c1);
    m.addBond(c1, c2);
    const om = m.addAtom({ el: 'O' });
    const cm = m.addAtom({ el: 'C' });
    m.addBond(2, om);
    m.addBond(om, cm);
    layoutMol(m, { fixed });
    before.forEach((p, i) => {
      expect(m.atoms[i].x).toBe(p.x);
      expect(m.atoms[i].y).toBe(p.y);
    });
    const g = geometry(m);
    expect(g.minBond).toBeGreaterThan(0.85);
    expect(g.maxBond).toBeLessThan(1.15);
    expect(g.minNB).toBeGreaterThan(0.6);
  });

  it('expands an abbreviation (phenyl) from a fixed attachment atom', () => {
    const m = parseSmiles('CC(=O)N');
    layoutMol(m);
    const fixed = new Set([0, 1, 2, 3]);
    const before = m.atoms.map((a) => ({ x: a.x, y: a.y }));
    const ring: number[] = [];
    for (let k = 0; k < 6; k++) ring.push(m.addAtom({ el: 'C' }));
    for (let k = 0; k < 6; k++) m.addBond(ring[k], ring[(k + 1) % 6], k % 2 ? 2 : 1);
    m.addBond(3, ring[0]);
    layoutMol(m, { fixed });
    for (const i of fixed) expect([m.atoms[i].x, m.atoms[i].y]).toEqual([before[i].x, before[i].y]);
    const g = geometry(m);
    expect(g.minBond).toBeGreaterThan(0.9);
    expect(g.maxBond).toBeLessThan(1.1);
    expect(g.minNB).toBeGreaterThan(0.8);
  });

  it('completes a partially drawn ring system around fixed atoms', () => {
    const m = parseSmiles('c1ccc2ccccc2c1');
    layoutMol(m);
    const fixed = new Set([0, 1, 2, 3, 4]);
    const before = m.atoms.map((a) => ({ x: a.x, y: a.y }));
    for (let i = 5; i < m.atoms.length; i++) { m.atoms[i].x = 0; m.atoms[i].y = 0; }
    layoutMol(m, { fixed });
    for (const i of fixed) expect([m.atoms[i].x, m.atoms[i].y]).toEqual([before[i].x, before[i].y]);
    const g = geometry(m);
    expect(g.minBond).toBeGreaterThan(0.95);
    expect(g.maxBond).toBeLessThan(1.05);
  });

  it('keeps wedges on fixed stereocentres', () => {
    const m = parseSmiles('C[C@@H](O)C(=O)O');
    layoutMol(m);
    const wedge = m.bonds.findIndex((b) => b.style === 'wedge' || b.style === 'hash');
    const snapshot = { ...m.bonds[wedge] };
    const fixed = new Set(m.atoms.map((_, i) => i));
    const c = m.addAtom({ el: 'C' });
    m.addBond(0, c);
    layoutMol(m, { fixed });
    expect(m.bonds[wedge]).toEqual(snapshot);
  });

  it('places a new fragment beside the fixed drawing', () => {
    const m = parseSmiles('CCO');
    layoutMol(m);
    const fixed = new Set([0, 1, 2]);
    const off = m.append(parseSmiles('c1ccccc1'));
    layoutMol(m, { fixed });
    const maxFixed = Math.max(...[0, 1, 2].map((i) => m.atoms[i].x));
    for (let i = off; i < m.atoms.length; i++) expect(m.atoms[i].x).toBeGreaterThan(maxFixed + 1);
  });
});

describe('clean2D', () => {
  const CASES: [string, string][] = [
    ['cholesterol', TRICKY[0][1]],
    ['sucrose', 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O'],
    ['retinol', 'CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/CO)/C)/C'],
    ['L-isoleucine', 'CC[C@H](C)[C@@H](C(=O)O)N'],
    ['octocrylene', 'CCCCC(CC)COC(=O)C(=C(C1=CC=CC=C1)C2=CC=CC=C2)C#N'],
    ['bakuchiol', 'CC(=CCC[C@@](C)(C=C)/C=C/C1=CC=C(C=C1)O)C'],
    ['ceramide NP', 'CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@@H]([C@@H](CCCCCCCCCCCCCC)O)O'],
    ['ibuprofen', 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O'],
    ['(E)-cyclododecene', 'C1CCCCC/C=C/CCCC1'],
  ];

  for (const [name, smiles] of CASES) {
    it(`restores a perturbed drawing of ${name}`, () => {
      const m = parseSmiles(smiles);
      const ref = writeSmiles(m);
      const centres = new Set(m.tetra.map((t) => t.center));
      const dbBonds = new Set(m.dbStereo.map((d) => d.bond));
      layoutMol(m);
      const r = rng(name.length * 7919);
      for (const a of m.atoms) {
        a.x += (r() - 0.5) * 0.5;
        a.y += (r() - 0.5) * 0.5;
      }
      const c0 = centroid(m);
      clean2D(m);
      const g = geometry(m);
      expect(g.minBond, name).toBeGreaterThan(0.9);
      expect(g.maxBond, name).toBeLessThan(1.1);
      expect(g.minNB, name).toBeGreaterThan(MIN_NONBONDED);
      const c1 = centroid(m);
      expect(Math.hypot(c1.x - c0.x, c1.y - c0.y), name).toBeLessThan(0.5);
      expect(drawnSmiles(m, centres, dbBonds), `${name} stereo`).toBe(ref);
    });
  }

  it('keeps overall orientation and position of a rotated, scaled drawing', () => {
    const m = parseSmiles('CCCCCCO');
    layoutMol(m);
    // rotate by 40°, scale ×1.6, translate
    const t = (40 * Math.PI) / 180;
    for (const a of m.atoms) {
      const x = a.x * 1.6, y = a.y * 1.6;
      a.x = Math.cos(t) * x - Math.sin(t) * y + 10;
      a.y = Math.sin(t) * x + Math.cos(t) * y + 5;
    }
    const c0 = centroid(m);
    const dir0 = Math.atan2(m.atoms[6].y - m.atoms[0].y, m.atoms[6].x - m.atoms[0].x);
    clean2D(m);
    const c1 = centroid(m);
    expect(Math.hypot(c1.x - c0.x, c1.y - c0.y)).toBeLessThan(0.05);
    const dir1 = Math.atan2(m.atoms[6].y - m.atoms[0].y, m.atoms[6].x - m.atoms[0].x);
    expect(Math.abs(dir1 - dir0)).toBeLessThan(0.1);
    expect(geometry(m).maxBond).toBeLessThan(1.02);
  });

  it('regularises a crude hand drawing of a substituted ring', () => {
    // toluene-like hexagon drawn with sloppy coordinates
    const m = new Mol();
    const pts = [[0, 0], [1.3, 0.2], [1.9, 1.1], [1.2, 2.2], [0.1, 1.8], [-0.4, 0.9], [-0.6, -0.9]];
    for (const [x, y] of pts) m.addAtom({ el: 'C', x, y });
    for (let k = 0; k < 6; k++) m.addBond(k, (k + 1) % 6, k % 2 ? 2 : 1);
    m.addBond(0, 6);
    clean2D(m);
    for (let k = 0; k < 6; k++) {
      const a = m.atoms[k], b = m.atoms[(k + 1) % 6];
      expect(Math.hypot(a.x - b.x, a.y - b.y)).toBeCloseTo(1, 1);
    }
    // methyl radial: 120° to both ring neighbours
    const ang = (i: number, c: number, j: number) => {
      const A = m.atoms[i], C = m.atoms[c], B = m.atoms[j];
      return (Math.acos(((A.x - C.x) * (B.x - C.x) + (A.y - C.y) * (B.y - C.y)) / (Math.hypot(A.x - C.x, A.y - C.y) * Math.hypot(B.x - C.x, B.y - C.y))) * 180) / Math.PI;
    };
    expect(ang(6, 0, 1)).toBeCloseTo(120, -1);
    expect(ang(6, 0, 5)).toBeCloseTo(120, -1);
  });

  it('moves only the requested atoms', () => {
    const m = parseSmiles('CCCCCC');
    layoutMol(m);
    m.atoms[5].x += 0.4;
    m.atoms[5].y -= 0.3;
    const before = m.atoms.map((a) => ({ x: a.x, y: a.y }));
    clean2D(m, { atoms: new Set([4, 5]) });
    for (let i = 0; i < 4; i++) expect([m.atoms[i].x, m.atoms[i].y]).toEqual([before[i].x, before[i].y]);
    const d = Math.hypot(m.atoms[4].x - m.atoms[5].x, m.atoms[4].y - m.atoms[5].y);
    expect(d).toBeCloseTo(1, 1);
  });

  it('restores wedge-defined configuration if cleaning disturbs it', () => {
    const smiles = 'N[C@@H](Cc1ccccc1)C(=O)O';
    const m = parseSmiles(smiles);
    const ref = writeSmiles(m);
    layoutMol(m);
    const r = rng(42);
    for (const a of m.atoms) { a.x += (r() - 0.5) * 0.8; a.y += (r() - 0.5) * 0.8; }
    // the perturbed drawing still encodes the configuration
    expect(drawnSmiles(m, new Set([1]), new Set())).toBe(ref);
    clean2D(m);
    expect(drawnSmiles(m, new Set([1]), new Set())).toBe(ref);
  });
});

describe('performance', () => {
  it('lays out ~100 heavy atoms quickly and handles 500+ atoms', () => {
    // a 100-heavy-atom drug-like construct: peptide of mixed residues
    const pep = 'N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCSC)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H]([C@@H](C)CC)C(=O)O';
    const probe = parseSmiles(pep);
    expect(probe.atoms.length).toBeGreaterThanOrEqual(90);
    for (let k = 0; k < 3; k++) layoutMol(parseSmiles(pep)); // warm-up
    const times: number[] = [];
    for (let k = 0; k < 5; k++) {
      const m = parseSmiles(pep);
      const t0 = performance.now();
      layoutMol(m);
      times.push(performance.now() - t0);
      checkLayout('peptide', pep, m);
    }
    times.sort((a, b) => a - b);
    // target is < 30 ms on a desktop machine; leave headroom for slow CI runners
    expect(times[2]).toBeLessThan(90);

    const big = 'C' + 'COCCOC'.repeat(60) + 'c1ccc(cc1)' + 'CCCCCCCCCC'.repeat(14) + 'O';
    const m = parseSmiles(big);
    expect(m.atoms.length).toBeGreaterThan(500);
    const t0 = performance.now();
    layoutMol(m);
    expect(performance.now() - t0).toBeLessThan(3000);
    checkLayout('500+ atoms', big, m);
  });
});

function centroid(m: Mol): { x: number; y: number } {
  let x = 0, y = 0;
  for (const a of m.atoms) { x += a.x; y += a.y; }
  return { x: x / m.atoms.length, y: y / m.atoms.length };
}
