import { Mol } from './mol';

export interface RingInfo {
  /** Smallest set of smallest rings; each ring is an ordered cycle of atom indices. */
  rings: number[][];
  /** Bond indices of each ring (ringBonds[k][j] joins rings[k][j] and rings[k][j+1]). */
  ringBonds: number[][];
  /** Ring indices each atom belongs to. */
  atomRings: number[][];
  /** Ring indices each bond belongs to. */
  bondRings: number[][];
  /** Ring systems: groups of ring indices connected through shared atoms (fused, spiro, bridged). */
  systems: number[][];
  /** True if the atom lies on any cycle (including large ones not in the SSSR). */
  inRing: boolean[];
  /** True if the bond lies on any cycle. */
  bondInRing: boolean[];
}

/** Bonds that lie on at least one cycle (non-bridges), via Tarjan's bridge-finding algorithm. */
export function cyclicBonds(mol: Mol): boolean[] {
  const n = mol.atoms.length;
  const disc = new Array(n).fill(-1);
  const low = new Array(n).fill(0);
  const isBridge = new Array(mol.bonds.length).fill(false);
  let t = 0;
  for (let s = 0; s < n; s++) {
    if (disc[s] >= 0) continue;
    // iterative DFS: stack of [vertex, parentBond, adjacency cursor]
    const stack: [number, number, number][] = [[s, -1, 0]];
    disc[s] = low[s] = t++;
    while (stack.length) {
      const top = stack[stack.length - 1];
      const [v, pb] = top;
      const adj = mol.adj[v];
      if (top[2] < adj.length) {
        const bi = adj[top[2]++];
        if (bi === pb) continue;
        const w = mol.other(bi, v);
        if (disc[w] < 0) {
          disc[w] = low[w] = t++;
          stack.push([w, bi, 0]);
        } else {
          low[v] = Math.min(low[v], disc[w]);
        }
      } else {
        stack.pop();
        if (stack.length) {
          const u = stack[stack.length - 1][0];
          low[u] = Math.min(low[u], low[v]);
          if (low[v] > disc[u]) isBridge[pb] = true;
        }
      }
    }
  }
  return isBridge.map((b) => !b);
}

/** Perceives the SSSR using Horton candidate cycles and GF(2) Gaussian elimination. */
export function perceiveRings(mol: Mol): RingInfo {
  const nA = mol.atoms.length;
  const nB = mol.bonds.length;
  const bondInRing = cyclicBonds(mol);
  const inRing = new Array(nA).fill(false);
  mol.bonds.forEach((b, i) => {
    if (bondInRing[i]) {
      inRing[b.a] = true;
      inRing[b.b] = true;
    }
  });

  // Cyclic subgraph adjacency
  const cadj: number[][] = mol.atoms.map(() => []);
  let cyclicEdgeCount = 0;
  mol.bonds.forEach((b, i) => {
    if (!bondInRing[i]) return;
    cadj[b.a].push(i);
    cadj[b.b].push(i);
    cyclicEdgeCount++;
  });
  const ringAtoms = [];
  for (let i = 0; i < nA; i++) if (inRing[i]) ringAtoms.push(i);

  // Components of cyclic subgraph
  const comp = new Array(nA).fill(-1);
  let nComp = 0;
  for (const s of ringAtoms) {
    if (comp[s] >= 0) continue;
    const st = [s];
    comp[s] = nComp;
    while (st.length) {
      const v = st.pop()!;
      for (const bi of cadj[v]) {
        const w = mol.other(bi, v);
        if (comp[w] < 0) {
          comp[w] = nComp;
          st.push(w);
        }
      }
    }
    nComp++;
  }
  const nRingsWanted = cyclicEdgeCount - ringAtoms.length + nComp;

  const rings: number[][] = [];
  const ringBonds: number[][] = [];
  if (nRingsWanted > 0) {
    // Candidate cycles (Horton): for every vertex v and edge (x,y): P(v,x) + (x,y) + P(y,v)
    const words = Math.ceil(nB / 32);
    const candidates = new Map<string, { atoms: number[]; bonds: number[]; vec: Uint32Array }>();
    for (const v of ringAtoms) {
      // BFS tree from v
      const dist = new Map<number, number>();
      const parentBond = new Map<number, number>();
      dist.set(v, 0);
      const q = [v];
      for (let qi = 0; qi < q.length; qi++) {
        const u = q[qi];
        for (const bi of cadj[u]) {
          const w = mol.other(bi, u);
          if (!dist.has(w)) {
            dist.set(w, dist.get(u)! + 1);
            parentBond.set(w, bi);
            q.push(w);
          }
        }
      }
      const pathTo = (x: number): { atoms: number[]; bonds: number[] } => {
        const atoms = [x];
        const bonds: number[] = [];
        let c = x;
        while (c !== v) {
          const pb = parentBond.get(c)!;
          bonds.push(pb);
          c = mol.other(pb, c);
          atoms.push(c);
        }
        return { atoms: atoms.reverse(), bonds: bonds.reverse() }; // v ... x
      };
      for (let bi = 0; bi < nB; bi++) {
        if (!bondInRing[bi]) continue;
        const { a: x, b: y } = mol.bonds[bi];
        if (!dist.has(x) || !dist.has(y)) continue;
        if (parentBond.get(x) === bi || parentBond.get(y) === bi) continue;
        const dx = dist.get(x)!;
        const dy = dist.get(y)!;
        if (Math.abs(dx - dy) > 1) continue;
        const px = pathTo(x);
        const py = pathTo(y);
        // paths must only share v
        const setX = new Set(px.atoms);
        let ok = true;
        for (let k = 1; k < py.atoms.length; k++) if (setX.has(py.atoms[k])) { ok = false; break; }
        if (!ok) continue;
        const atoms = [...px.atoms, ...py.atoms.slice(1).reverse()];
        const bonds = [...px.bonds, bi, ...py.bonds.slice().reverse()];
        if (atoms.length < 3) continue;
        const key = [...bonds].sort((p, q2) => p - q2).join(',');
        if (candidates.has(key)) continue;
        const vec = new Uint32Array(words);
        for (const b of bonds) vec[b >>> 5] |= 1 << (b & 31);
        candidates.set(key, { atoms, bonds, vec });
      }
    }
    const sorted = [...candidates.values()].sort((p, q2) => p.atoms.length - q2.atoms.length);
    // Gaussian elimination over GF(2)
    const basis: { vec: Uint32Array; pivot: number }[] = [];
    const lowestBit = (v: Uint32Array): number => {
      for (let w = 0; w < v.length; w++) if (v[w]) return w * 32 + (31 - Math.clz32(v[w] & -v[w]));
      return -1;
    };
    for (const c of sorted) {
      if (rings.length >= nRingsWanted) break;
      const v = c.vec.slice();
      // reduce
      let changed = true;
      while (changed) {
        changed = false;
        const p = lowestBit(v);
        if (p < 0) break;
        for (const bv of basis) {
          if (bv.pivot === p) {
            for (let w = 0; w < v.length; w++) v[w] ^= bv.vec[w];
            changed = true;
            break;
          }
        }
      }
      const p = lowestBit(v);
      if (p < 0) continue; // dependent
      basis.push({ vec: v, pivot: p });
      rings.push(c.atoms);
      ringBonds.push(c.bonds);
    }
  }

  const atomRings: number[][] = mol.atoms.map(() => []);
  const bondRings: number[][] = mol.bonds.map(() => []);
  rings.forEach((r, k) => {
    for (const a of r) atomRings[a].push(k);
    for (const b of ringBonds[k]) bondRings[b].push(k);
  });

  // ring systems via union-find on shared atoms
  const parent = rings.map((_, i) => i);
  const find = (x: number): number => (parent[x] === x ? x : (parent[x] = find(parent[x])));
  for (const list of atomRings) {
    for (let k = 1; k < list.length; k++) {
      const a = find(list[0]);
      const b = find(list[k]);
      if (a !== b) parent[a] = b;
    }
  }
  const sysMap = new Map<number, number[]>();
  rings.forEach((_, i) => {
    const r = find(i);
    if (!sysMap.has(r)) sysMap.set(r, []);
    sysMap.get(r)!.push(i);
  });

  return { rings, ringBonds, atomRings, bondRings, systems: [...sysMap.values()], inRing, bondInRing };
}

/** Atoms of a ring system (union of its rings' atoms). */
export function systemAtoms(info: RingInfo, system: number[]): number[] {
  const s = new Set<number>();
  for (const r of system) for (const a of info.rings[r]) s.add(a);
  return [...s];
}

/** Size of the smallest SSSR ring containing the bond (0 if acyclic). */
export function smallestRingSizeOfBond(info: RingInfo, bi: number): number {
  let m = 0;
  for (const r of info.bondRings[bi]) {
    const n = info.rings[r].length;
    if (!m || n < m) m = n;
  }
  return m;
}
