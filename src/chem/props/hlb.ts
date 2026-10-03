// Hydrophile–lipophile balance for amphiphiles (surfactants, emulsifiers, fatty alcohols/acids…).
//
// Davies (1957) group contributions:  HLB = 7 + Σ hydrophilic − Σ lipophilic, with
//   −SO4⁻ (Na⁺) 38.7 · −COO⁻K⁺ 21.1 · −COO⁻Na⁺ (any other counter-ion) 19.1 · sulfonate −SO3⁻ 11.0 ·
//   amine N (tertiary amine value, also used for 1°/2° amines and quaternary ammonium) 9.4 ·
//   ester on a sorbitan ring moiety 6.8 · free ester 2.4 · −COOH 2.1 · free −OH 1.9 · ether −O− 1.3 ·
//   −OH on a sorbitan ring moiety 0.5 · −(CH2CH2O)− 0.33 per unit ·
//   −(CH(CH3)CH2O)− −0.15 per unit · every remaining CH/CH2/CH3/=CH− (and C) −0.475 · CF2/CF3 −0.870.
//   Amides are not in Davies' table; they are approximated as a free ester (2.4).
//   "Sorbitan ring moiety" = a saturated oxolane (5-ring with one O) whose carbons bear O substituents,
//   plus O-bearing sp3 carbons attached to it (C5/C6 side chain of sorbitan). Sugar pyranose rings are
//   treated as free OH/ether groups.
//
// Griffin (1949/1954), non-ionic surfactants only:  HLB = 20 · Mh / M.
//   Partition (Mh = hydrophilic mass):
//   1. ethylene-oxide units −CH2CH2O− (C2H4O each; a chain of n CH2CH2 links between n+1 oxygens
//      takes n of them, leaving the terminal −OH or the ether O on the far side out — Griffin's E/5);
//   2. polyol moieties: ≥ 2 connected sp3 carbons that each carry an O (glycerol, sorbitan — also when
//      ethoxylated as in polysorbates — sugar residues), counted with their H atoms, their free −OH
//      groups and the ring ether O of the moiety; ester O atoms and acyl groups stay on the lipophilic
//      side (equivalent to Griffin's 20·(1 − S/A) for polyol fatty esters; (E + P)/5 for polysorbates);
//   3. only when neither 1 nor 2 is present: the free polar head groups −OH (O+H) and −COOH.
//   Propylene-oxide units count as lipophilic (as in Griffin's EO-weight-% formula).
//
// HLB is only reported for amphiphilic molecules: a contiguous tail of ≥ 6 carbons bearing no
// heteroatom plus at least one hydrophilic group.
import { MatchGraph } from './smarts';

export interface HLBResult {
  griffin: number | null;
  davies: number | null;
  note: string;
}

const MASS = (z: number) => (z === 1 ? 1.00794 : z === 6 ? 12.0107 : z === 7 ? 14.0067 : z === 8 ? 15.9994 : 0);
const ALKALI = new Set([3, 11, 19, 37, 55]);

export function computeHLB(g: MatchGraph, mw: number): HLBResult | null {
  const n = g.n;
  const isC = (i: number) => g.z[i] === 6;
  const heavyNbrs = (i: number) => g.nbr[i].filter((j) => g.z[j] > 1);
  const isHetero = (i: number) => g.z[i] > 1 && g.z[i] !== 6;
  const bondOrder = (a: number, b: number) => {
    const t = g.nbr[a].indexOf(b);
    return t < 0 ? 0 : g.bondArom[g.nbrBond[a][t]] ? 1.5 : g.bondOrder[g.nbrBond[a][t]];
  };
  const fCount = (i: number) => g.nbr[i].filter((j) => g.z[j] === 9).length;

  // ── Tail detection ────────────────────────────────────────────────
  const tailC = new Array(n).fill(false);
  for (let i = 0; i < n; i++) {
    if (!isC(i) || g.arom[i]) continue;
    const het = heavyNbrs(i).filter((j) => isHetero(j));
    // hydrocarbon carbon, or perfluoro carbon (only F as heteroatoms, ≥ 2 F)
    if (het.length === 0 || (het.every((j) => g.z[j] === 9) && het.length >= 2)) tailC[i] = true;
  }
  let maxTail = 0;
  const seen = new Array(n).fill(false);
  for (let s = 0; s < n; s++) {
    if (!tailC[s] || seen[s]) continue;
    let size = 0;
    const st = [s];
    seen[s] = true;
    while (st.length) {
      const v = st.pop()!;
      size++;
      for (const w of g.nbr[v]) if (tailC[w] && !seen[w]) { seen[w] = true; st.push(w); }
    }
    maxTail = Math.max(maxTail, size);
  }
  if (maxTail < 6) return null;

  // ── Group classification ──────────────────────────────────────────
  const used = new Array(n).fill(false); // atoms consumed by a hydrophilic group
  let davies = 7;
  let ionic = false;
  let hasHead = false;
  let approxAmide = false;
  const unknown = new Set<string>();
  const hasK = g.z.some((z, i) => z === 19 && g.charge[i] >= 0);
  const oxygens = (i: number) => g.nbr[i].filter((j) => g.z[j] === 8);
  const isOH = (o: number) => g.z[o] === 8 && g.hTotal[o] >= 1 && heavyNbrs(o).length === 1;
  const isAnionicO = (o: number) => g.charge[o] < 0 || heavyNbrs(o).some((j) => ALKALI.has(g.z[j]));
  for (let i = 0; i < n; i++) if (g.charge[i] !== 0 && g.z[i] > 1) ionic = true;

  // Sorbitan ring moiety (see header)
  const sorbitan = new Array(n).fill(false);
  for (const ring of g.rings.rings) {
    if (ring.length !== 5) continue;
    const gi = ring.map((a) => g.fromMol[a]);
    if (gi.some((k) => k < 0 || g.arom[k])) continue;
    const ringO = gi.filter((k) => g.z[k] === 8);
    const ringC = gi.filter((k) => isC(k));
    if (ringO.length !== 1 || ringC.length !== 4) continue;
    const oBearing = ringC.filter((k) => oxygens(k).some((o) => !gi.includes(o)));
    if (oBearing.length < 2) continue;
    for (const k of gi) sorbitan[k] = true;
    // O-bearing sp3 side-chain carbons (up to two bonds away)
    let frontier = ringC;
    for (let depth = 0; depth < 2; depth++) {
      const next: number[] = [];
      for (const c of frontier)
        for (const w of g.nbr[c])
          if (isC(w) && !sorbitan[w] && !g.arom[w] && g.degree[w] + g.hImplicit[w] === 4 && oxygens(w).length) {
            sorbitan[w] = true;
            next.push(w);
          }
      frontier = next;
    }
  }
  const onSorbitan = (o: number) => heavyNbrs(o).some((c) => sorbitan[c] && isC(c));

  // Sulfates, sulfonates, carboxylates/acids, esters, amides
  for (let i = 0; i < n; i++) {
    if (g.z[i] === 16) {
      const os = oxygens(i);
      const cs = heavyNbrs(i).filter(isC);
      if (os.length === 4) {
        davies += 38.7; // sulfate (ester) salt
        ionic = true;
        hasHead = true;
        used[i] = true;
        os.forEach((o) => (used[o] = true));
      } else if (os.length === 3 && cs.length === 1) {
        davies += 11.0; // sulfonate
        ionic = true;
        hasHead = true;
        used[i] = true;
        os.forEach((o) => (used[o] = true));
      }
    }
    if (!isC(i)) continue;
    const os = oxygens(i);
    const dblO = os.find((o) => bondOrder(i, o) === 2);
    if (dblO === undefined) continue;
    const single = os.filter((o) => o !== dblO);
    const ns = heavyNbrs(i).filter((j) => g.z[j] === 7);
    if (single.length === 1 && heavyNbrs(i).length <= 3) {
      const o = single[0];
      if (isAnionicO(o)) {
        davies += hasK ? 21.1 : 19.1;
        ionic = true;
      } else if (isOH(o)) davies += 2.1;
      else davies += onSorbitan(o) ? 6.8 : 2.4; // ester (O also bonded to a carbon)
      hasHead = true;
      used[i] = used[o] = used[dblO] = true;
    } else if (single.length === 0 && ns.length === 1) {
      davies += 2.4; // amide ≈ free ester
      approxAmide = true;
      hasHead = true;
      used[i] = used[dblO] = used[ns[0]] = true;
    }
  }

  // EO / PO chains: C–C links between two oxygens
  const links: { c1: number; c2: number; o1: number; o2: number; po: boolean; methyl: number }[] = [];
  for (let a = 0; a < n; a++) {
    if (!isC(a) || used[a] || g.arom[a]) continue;
    for (const b of g.nbr[a]) {
      if (b <= a || !isC(b) || used[b] || bondOrder(a, b) !== 1) continue;
      // link oxygens may belong to an ester or sulfate (laureth sulfate, polysorbate laurate);
      // such oxygens are never taken by the EO unit itself
      const oa = oxygens(a).filter((o) => bondOrder(a, o) === 1);
      const ob = oxygens(b).filter((o) => bondOrder(b, o) === 1);
      if (oa.length !== 1 || ob.length !== 1 || oa[0] === ob[0]) continue;
      const ha = g.hTotal[a], hb = g.hTotal[b];
      const extraA = heavyNbrs(a).filter((j) => j !== b && j !== oa[0]);
      const extraB = heavyNbrs(b).filter((j) => j !== a && j !== ob[0]);
      if (ha === 2 && hb === 2 && !extraA.length && !extraB.length) {
        links.push({ c1: a, c2: b, o1: oa[0], o2: ob[0], po: false, methyl: -1 });
      } else {
        // propylene oxide: one CH2 and one CH carrying a methyl
        const [ch, ch2, ex, ex2] = ha === 1 ? [a, b, extraA, extraB] : [b, a, extraB, extraA];
        if (g.hTotal[ch] === 1 && g.hTotal[ch2] === 2 && ex.length === 1 && !ex2.length && g.hTotal[ex[0]] === 3 && isC(ex[0]))
          links.push({ c1: a, c2: b, o1: oa[0], o2: ob[0], po: true, methyl: ex[0] });
      }
    }
  }
  // assign one oxygen per link, walking chains from an end so that a terminal OH stays free
  const linkOf = new Map<number, number[]>(); // oxygen → link indices
  links.forEach((l, k) => {
    for (const o of [l.o1, l.o2]) {
      if (!linkOf.has(o)) linkOf.set(o, []);
      linkOf.get(o)!.push(k);
    }
  });
  const linkDone = new Array(links.length).fill(false);
  const oTaken = new Set<number>();
  let eo = 0, po = 0;
  let eoMass = 0;
  const assignChain = (startO: number) => {
    let o = startO;
    for (;;) {
      const k = (linkOf.get(o) ?? []).find((x) => !linkDone[x]);
      if (k === undefined) break;
      linkDone[k] = true;
      const l = links[k];
      const far = l.o1 === o ? l.o2 : l.o1;
      // the unit takes the oxygen on the start side unless it is a free OH or already assigned;
      // a CH2CH2 between two ester/sulfate oxygens (glycol diesters) is not an EO unit
      const free = (x: number) => !oTaken.has(x) && !used[x];
      const take = free(o) && !isOH(o) ? o : free(far) && !isOH(far) ? far : free(far) ? far : free(o) ? o : -1;
      if (take < 0) {
        o = far;
        continue;
      }
      oTaken.add(take);
      used[l.c1] = used[l.c2] = used[take] = true;
      if (l.po) {
        po++;
        used[l.methyl] = true;
      } else {
        eo++;
        eoMass += 2 * MASS(6) + 4 * MASS(1) + MASS(8);
      }
      o = far;
    }
  };
  // chain ends: oxygens with a single link; start from an end whose oxygen the unit may take
  // (ether O on the hydrophobe side) so that a terminal OH stays free
  const ends = [...linkOf.entries()].filter(([, ks]) => ks.length === 1).map(([o]) => o);
  const endRank = (o: number) => (isOH(o) ? 2 : used[o] ? 1 : 0);
  ends.sort((p, q) => endRank(p) - endRank(q));
  for (const o of ends) assignChain(o);
  for (const [o] of linkOf) assignChain(o); // rings of EO units (crown ethers)
  davies += 0.33 * eo - 0.15 * po;
  if (eo) hasHead = true;

  // Polyol moieties (Griffin): ≥ 2 connected sp3 carbons that each carry an oxygen (free OH, ring
  // ether, ester or EO linkage) — glycerol, sorbitan (also when ethoxylated), sugar residues
  const polyolC = new Array(n).fill(false);
  const cand = (i: number) => isC(i) && !used[i] && !g.arom[i] && g.degree[i] + g.hImplicit[i] === 4 && oxygens(i).length > 0;
  const visited = new Array(n).fill(false);
  for (let s = 0; s < n; s++) {
    if (!cand(s) || visited[s]) continue;
    const comp: number[] = [];
    const st = [s];
    visited[s] = true;
    while (st.length) {
      const v = st.pop()!;
      comp.push(v);
      for (const w of g.nbr[v]) if (!visited[w] && cand(w)) { visited[w] = true; st.push(w); }
    }
    if (comp.length >= 2) for (const c of comp) polyolC[c] = true;
  }

  // Hydroxyls, ethers, amines
  let headMass = 0; // Griffin fallback (free OH / COOH heads)
  for (let i = 0; i < n; i++) {
    if (used[i]) continue;
    if (g.z[i] === 8) {
      const hn = heavyNbrs(i);
      if (isOH(i) && isC(hn[0])) {
        davies += onSorbitan(i) ? 0.5 : 1.9;
        headMass += MASS(8) + MASS(1);
        hasHead = true;
      } else if (hn.length === 2 && hn.every(isC)) {
        davies += 1.3;
        hasHead = true;
      } else if (!(hn.length === 1 && bondOrder(i, hn[0]) === 2 && isC(hn[0]))) {
        // ketone/aldehyde C=O has no Davies value and is ignored silently
        if (hn.length) unknown.add('O');
      }
      used[i] = true;
    } else if (g.z[i] === 7) {
      davies += 9.4;
      hasHead = true;
      used[i] = true;
    } else if (g.z[i] === 15 || g.z[i] === 16) {
      unknown.add(g.z[i] === 15 ? 'P' : 'S');
    }
  }
  // COOH heads for the Griffin fallback
  for (let i = 0; i < n; i++) {
    if (!isC(i)) continue;
    const os = oxygens(i);
    if (os.length === 2 && os.some((o) => bondOrder(i, o) === 2) && os.some((o) => isOH(o))) headMass += MASS(6) + 2 * MASS(8) + MASS(1);
  }
  // Lipophilic carbons
  for (let i = 0; i < n; i++) {
    if (!isC(i) || used[i]) continue;
    davies += fCount(i) >= 2 ? -0.87 : -0.475;
  }
  if (!hasHead) return null;

  // ── Griffin ──────────────────────────────────────────────────────
  let griffin: number | null = null;
  if (!ionic) {
    let mh = eoMass;
    let polyMass = 0;
    for (let i = 0; i < n; i++) {
      if (!polyolC[i]) continue;
      polyMass += MASS(6) + g.hTotal[i] * MASS(1);
      for (const o of oxygens(i)) if (isOH(o) && !oTaken.has(o)) polyMass += MASS(8) + MASS(1);
    }
    // ring ether oxygens inside a polyol moiety (sorbitan/sugar ring O)
    for (let i = 0; i < n; i++) {
      if (g.z[i] !== 8 || oTaken.has(i) || g.ringCount[i] === 0) continue;
      const hn = heavyNbrs(i);
      if (hn.length === 2 && hn.every((c) => polyolC[c])) polyMass += MASS(8);
    }
    mh += polyMass;
    if (mh === 0) mh = headMass;
    if (mh > 0 && mw > 0) griffin = (20 * mh) / mw;
  }

  const main = griffin ?? davies;
  const use =
    main < 3 ? 'oil-soluble; antifoam / co-emulsifier' :
    main < 7 ? 'W/O emulsifier range (3–6)' :
    main < 9 ? 'wetting / spreading agent (7–9)' :
    main < 13 ? 'O/W emulsifier range (8–16)' :
    main < 16 ? 'O/W emulsifier / detergent (13–16)' :
    'solubiliser / hydrotrope (16–18+)';
  const parts = [
    `Davies (7 + Σ group numbers) = ${davies.toFixed(1)}${ionic ? ' (ionic head)' : ''}`,
    griffin !== null
      ? `Griffin (20·Mh/M) = ${griffin.toFixed(1)}`
      : ionic
        ? 'Griffin not applicable to ionic surfactants'
        : 'Griffin not defined (no EO, polyol or free OH/COOH head)',
    `HLB ≈ ${main.toFixed(1)}: ${use}`,
  ];
  if (eo) parts.push(`${eo} EO unit${eo > 1 ? 's' : ''}`);
  if (po) parts.push(`${po} PO unit${po > 1 ? 's' : ''}`);
  if (approxAmide) parts.push('amide approximated as ester (not in Davies table)');
  if (unknown.size) parts.push(`groups with ${[...unknown].join(', ')} not covered by Davies' table`);
  return { griffin: griffin === null ? null : round(griffin, 2), davies: round(davies, 3), note: parts.join('; ') };
}

function round(x: number, d: number): number {
  const f = 10 ** d;
  return Math.round(x * f) / f;
}
