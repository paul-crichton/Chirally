// Numerical terms and multiplying prefixes (IUPAC 2013, P-14.2).

const UNITS = ['', 'hen', 'do', 'tri', 'tetra', 'penta', 'hexa', 'hepta', 'octa', 'nona'];
const TENS = ['', 'deca', 'icosa', 'triaconta', 'tetraconta', 'pentaconta', 'hexaconta', 'heptaconta', 'octaconta', 'nonaconta'];
const HUNDREDS = ['', 'hecta', 'dicta', 'tricta', 'tetracta', 'pentacta', 'hexacta', 'heptacta', 'octacta', 'nonacta'];
const THOUSANDS = ['', 'kilia', 'dilia', 'trilia', 'tetralia', 'pentalia', 'hexalia', 'heptalia', 'octalia', 'nonalia'];

/**
 * Basic numerical term with its final "a" (e.g. 5 → "penta", 11 → "undeca", 21 → "henicosa",
 * 23 → "tricosa", 100 → "hecta"). Used to build alkane stems and multiplying prefixes.
 */
export function numeralTerm(n: number): string {
  if (n <= 0 || n >= 10000 || !Number.isInteger(n)) throw new Error('numeralTerm: unsupported ' + n);
  if (n === 1) return 'mono';
  if (n === 2) return 'di';
  if (n === 3) return 'tri';
  if (n === 11) return 'undeca';
  const u = n % 10;
  const t = Math.floor(n / 10) % 10;
  const h = Math.floor(n / 100) % 10;
  const k = Math.floor(n / 1000);
  let unit = UNITS[u];
  if (n < 10) unit = u === 1 ? 'mono' : u === 2 ? 'di' : unit;
  // "undeca" for 11 inside larger numbers too (111 = undecahecta)
  let tens = TENS[t];
  if (t === 1 && u === 1) {
    unit = 'un';
  }
  // "icosa": the initial i is elided after a vowel (tricosa, docosa) but kept after "hen"
  if (t === 2 && unit && /[aeiou]$/.test(unit)) tens = 'cosa';
  let s = unit + tens + HUNDREDS[h] + THOUSANDS[k];
  if (n === 20) s = 'icosa';
  return s;
}

/** Alkane stem without the final "a": 1 meth, 2 eth, 3 prop, 4 but, 5 pent, 11 undec, 20 icos, 21 henicos. */
export function alkaneStem(n: number): string {
  if (n === 1) return 'meth';
  if (n === 2) return 'eth';
  if (n === 3) return 'prop';
  if (n === 4) return 'but';
  const t = numeralTerm(n);
  return t.endsWith('a') ? t.slice(0, -1) : t;
}

/** Simple multiplying prefix: 1 '', 2 di, 3 tri, 4 tetra, … */
export function multiplier(n: number): string {
  if (n <= 1) return '';
  return numeralTerm(n);
}

/** Multiplying prefix for complex entities: 2 bis, 3 tris, 4 tetrakis, … */
export function multiplierComplex(n: number): string {
  if (n <= 1) return '';
  if (n === 2) return 'bis';
  if (n === 3) return 'tris';
  return numeralTerm(n) + 'kis';
}
