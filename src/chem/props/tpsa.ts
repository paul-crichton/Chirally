// Topological polar surface area: P. Ertl, B. Rohde, P. Selzer, J. Med. Chem. 2000, 43, 3714–3717.
// Fragment contributions of N and O (RDKit default: S and P excluded), classified by heavy-atom
// neighbour count, hydrogen count, formal charge, bond types (aromatic bonds counted separately)
// and three-membered-ring membership — the same scheme as RDKit's CalcTPSA, so values agree.
import { MatchGraph } from './smarts';

/** Per graph atom TPSA contribution (graph built with hMode 'none'). */
export function tpsaContribs(g: MatchGraph): number[] {
  const out = new Array(g.n).fill(0);
  for (let i = 0; i < g.n; i++) {
    const z = g.z[i];
    if (z !== 7 && z !== 8) continue;
    let nbrs = 0, sing = 0, doub = 0, trip = 0, arom = 0;
    for (let t = 0; t < g.nbr[i].length; t++) {
      const j = g.nbr[i][t];
      if (g.z[j] === 1) continue; // explicit H atoms are already in hTotal
      nbrs++;
      const b = g.nbrBond[i][t];
      if (g.bondArom[b]) arom++;
      else if (g.bondOrder[b] === 1) sing++;
      else if (g.bondOrder[b] === 2) doub++;
      else if (g.bondOrder[b] === 3) trip++;
    }
    const h = g.hTotal[i];
    const chg = g.charge[i];
    const in3 = g.ringSizes[i].includes(3);
    let v = -1;
    if (z === 7) {
      switch (nbrs) {
        case 1:
          if (h === 0 && chg === 0 && trip === 1) v = 23.79; // N#
          else if (h === 1 && chg === 0 && doub === 1) v = 23.85; // [NH]=
          else if (h === 2 && chg === 0 && sing === 1) v = 26.02; // [NH2]-
          else if (h === 2 && chg === 1 && doub === 1) v = 25.59; // [NH2+]=
          else if (h === 3 && chg === 1 && sing === 1) v = 27.64; // [NH3+]-
          break;
        case 2:
          if (h === 0 && chg === 0 && sing === 1 && doub === 1) v = 12.36; // =N-
          else if (h === 0 && chg === 0 && trip === 1 && doub === 1) v = 13.6; // =N# (azide middle)
          else if (h === 1 && chg === 0 && sing === 2 && in3) v = 21.94; // aziridine NH
          else if (h === 1 && chg === 0 && sing === 2) v = 12.03; // -NH-
          else if (h === 0 && chg === 1 && trip === 1 && sing === 1) v = 4.36; // -N+#
          else if (h === 1 && chg === 1 && doub === 1 && sing === 1) v = 13.97; // -[NH+]=
          else if (h === 2 && chg === 1 && sing === 2) v = 16.61; // -[NH2+]-
          else if (h === 0 && chg === 0 && arom === 2) v = 12.89; // :n:
          else if (h === 1 && chg === 0 && arom === 2) v = 15.79; // :[nH]:
          else if (h === 1 && chg === 1 && arom === 2) v = 14.14; // :[nH+]:
          break;
        case 3:
          if (h === 0 && chg === 0 && sing === 3 && in3) v = 3.01; // aziridine N
          else if (h === 0 && chg === 0 && sing === 3) v = 3.24; // tertiary amine
          else if (h === 0 && chg === 0 && sing === 1 && doub === 2) v = 11.68; // nitro (pentavalent form)
          else if (h === 0 && chg === 1 && sing === 2 && doub === 1) v = 3.01; // -N+(=)- (nitro, charge-separated)
          else if (h === 1 && chg === 1 && sing === 3) v = 4.44; // -[NH+](-)-
          else if (h === 0 && chg === 0 && arom === 3) v = 4.41; // :n(:):
          else if (h === 0 && chg === 0 && sing === 1 && arom === 2) v = 4.93; // -n(:):
          else if (h === 0 && chg === 0 && doub === 1 && arom === 2) v = 8.39; // =n(:):
          else if (h === 0 && chg === 1 && arom === 3) v = 4.1; // :[n+](:):
          else if (h === 0 && chg === 1 && sing === 1 && arom === 2) v = 3.88; // -[n+](:):
          break;
        case 4:
          if (h === 0 && sing === 4 && chg === 1) v = 0; // quaternary ammonium
          break;
      }
      // fallback for unusual N environments (as in RDKit)
      if (v < 0) v = Math.max(0, 30.5 - nbrs * 8.2 + h * 1.5);
    } else {
      switch (nbrs) {
        case 1:
          if (h === 0 && chg === 0 && doub === 1) v = 17.07; // O=
          else if (h === 1 && chg === 0 && sing === 1) v = 20.23; // -OH
          else if (h === 0 && chg === -1 && sing === 1) v = 23.06; // -O(-)
          break;
        case 2:
          if (h === 0 && chg === 0 && sing === 2 && in3) v = 12.53; // epoxide
          else if (h === 0 && chg === 0 && sing === 2) v = 9.23; // ether
          else if (h === 0 && chg === 0 && arom === 2) v = 13.14; // :o:
          break;
      }
      if (v < 0) v = Math.max(0, 28.5 - nbrs * 8.6 + h * 1.5);
    }
    out[i] = v;
  }
  return out;
}
