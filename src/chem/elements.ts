// Periodic table data used throughout the chemistry engine.
// Columns: Z | symbol | name | standard atomic weight | monoisotopic mass (most abundant isotope)
//          | CPK colour (Jmol) | covalent radius (Å, Cordero 2008) | Pauling electronegativity
//          | IUPAC group (0 = f-block) | period

export interface ElementData {
  z: number;
  symbol: string;
  name: string;
  mass: number;
  mono: number;
  color: string;
  covRadius: number;
  en: number;
  group: number;
  period: number;
}

const RAW = `
1 H Hydrogen 1.00794 1.0078250319 FFFFFF 0.31 2.20 1 1
2 He Helium 4.002602 4.0026032 D9FFFF 0.28 0 18 1
3 Li Lithium 6.941 7.0160034 CC80FF 1.28 0.98 1 2
4 Be Beryllium 9.012182 9.0121822 C2FF00 0.96 1.57 2 2
5 B Boron 10.811 11.0093054 FFB5B5 0.84 2.04 13 2
6 C Carbon 12.0107 12.0 909090 0.76 2.55 14 2
7 N Nitrogen 14.0067 14.0030740052 3050F8 0.71 3.04 15 2
8 O Oxygen 15.9994 15.9949146221 FF0D0D 0.66 3.44 16 2
9 F Fluorine 18.9984032 18.9984032 90E050 0.57 3.98 17 2
10 Ne Neon 20.1797 19.9924402 B3E3F5 0.58 0 18 2
11 Na Sodium 22.98976928 22.98976966 AB5CF2 1.66 0.93 1 3
12 Mg Magnesium 24.305 23.98504187 8AFF00 1.41 1.31 2 3
13 Al Aluminium 26.9815386 26.98153841 BFA6A6 1.21 1.61 13 3
14 Si Silicon 28.0855 27.97692649 F0C8A0 1.11 1.90 14 3
15 P Phosphorus 30.973762 30.97376151 FF8000 1.07 2.19 15 3
16 S Sulfur 32.065 31.97207069 FFFF30 1.05 2.58 16 3
17 Cl Chlorine 35.453 34.96885271 1FF01F 1.02 3.16 17 3
18 Ar Argon 39.948 39.9623831 80D1E3 1.06 0 18 3
19 K Potassium 39.0983 38.9637069 8F40D4 2.03 0.82 1 4
20 Ca Calcium 40.078 39.9625912 3DFF00 1.76 1.00 2 4
21 Sc Scandium 44.955912 44.9559102 E6E6E6 1.70 1.36 3 4
22 Ti Titanium 47.867 47.9479471 BFC2C7 1.60 1.54 4 4
23 V Vanadium 50.9415 50.9439637 A6A6AB 1.53 1.63 5 4
24 Cr Chromium 51.9961 51.9405119 8A99C7 1.39 1.66 6 4
25 Mn Manganese 54.938045 54.9380496 9C7AC7 1.39 1.55 7 4
26 Fe Iron 55.845 55.9349421 E06633 1.32 1.83 8 4
27 Co Cobalt 58.933195 58.9332002 F090A0 1.26 1.88 9 4
28 Ni Nickel 58.6934 57.9353479 50D050 1.24 1.91 10 4
29 Cu Copper 63.546 62.9296011 C88033 1.32 1.90 11 4
30 Zn Zinc 65.38 63.9291466 7D80B0 1.22 1.65 12 4
31 Ga Gallium 69.723 68.925581 C28F8F 1.22 1.81 13 4
32 Ge Germanium 72.64 73.9211782 668F8F 1.20 2.01 14 4
33 As Arsenic 74.9216 74.9215964 BD80E3 1.19 2.18 15 4
34 Se Selenium 78.96 79.9165218 FFA100 1.20 2.55 16 4
35 Br Bromine 79.904 78.9183376 A62929 1.20 2.96 17 4
36 Kr Krypton 83.798 83.911507 5CB8D1 1.16 3.00 18 4
37 Rb Rubidium 85.4678 84.9117893 702EB0 2.20 0.82 1 5
38 Sr Strontium 87.62 87.9056143 00FF00 1.95 0.95 2 5
39 Y Yttrium 88.90585 88.9058479 94FFFF 1.90 1.22 3 5
40 Zr Zirconium 91.224 89.9047037 94E0E0 1.75 1.33 4 5
41 Nb Niobium 92.90638 92.9063775 73C2C9 1.64 1.6 5 5
42 Mo Molybdenum 95.96 97.9054078 54B5B5 1.54 2.16 6 5
43 Tc Technetium 98 97.907216 3B9E9E 1.47 1.9 7 5
44 Ru Ruthenium 101.07 101.9043495 248F8F 1.46 2.2 8 5
45 Rh Rhodium 102.9055 102.905504 0A7D8C 1.42 2.28 9 5
46 Pd Palladium 106.42 105.903483 006985 1.39 2.20 10 5
47 Ag Silver 107.8682 106.905093 C0C0C0 1.45 1.93 11 5
48 Cd Cadmium 112.411 113.9033581 FFD98F 1.44 1.69 12 5
49 In Indium 114.818 114.903878 A67573 1.42 1.78 13 5
50 Sn Tin 118.71 119.9021966 668080 1.39 1.96 14 5
51 Sb Antimony 121.76 120.903818 9E63B5 1.39 2.05 15 5
52 Te Tellurium 127.6 129.9062228 D47A00 1.38 2.1 16 5
53 I Iodine 126.90447 126.904468 940094 1.39 2.66 17 5
54 Xe Xenon 131.293 131.9041545 429EB0 1.40 2.6 18 5
55 Cs Caesium 132.9054519 132.905447 57178F 2.44 0.79 1 6
56 Ba Barium 137.327 137.905241 00C900 2.15 0.89 2 6
57 La Lanthanum 138.90547 138.906348 70D4FF 2.07 1.10 0 6
58 Ce Cerium 140.116 139.905434 FFFFC7 2.04 1.12 0 6
59 Pr Praseodymium 140.90765 140.907648 D9FFC7 2.03 1.13 0 6
60 Nd Neodymium 144.242 141.907719 C7FFC7 2.01 1.14 0 6
61 Pm Promethium 145 144.912744 A3FFC7 1.99 1.13 0 6
62 Sm Samarium 150.36 151.919728 8FFFC7 1.98 1.17 0 6
63 Eu Europium 151.964 152.921226 61FFC7 1.98 1.2 0 6
64 Gd Gadolinium 157.25 157.924101 45FFC7 1.96 1.20 0 6
65 Tb Terbium 158.92535 158.925343 30FFC7 1.94 1.2 0 6
66 Dy Dysprosium 162.5 163.929171 1FFFC7 1.92 1.22 0 6
67 Ho Holmium 164.93032 164.930319 00FF9C 1.92 1.23 0 6
68 Er Erbium 167.259 165.93029 00E675 1.89 1.24 0 6
69 Tm Thulium 168.93421 168.934211 00D452 1.90 1.25 0 6
70 Yb Ytterbium 173.054 173.938858 00BF38 1.87 1.1 0 6
71 Lu Lutetium 174.9668 174.940768 00AB24 1.87 1.27 3 6
72 Hf Hafnium 178.49 179.946549 4DC2FF 1.75 1.3 4 6
73 Ta Tantalum 180.94788 180.947996 4DA6FF 1.70 1.5 5 6
74 W Tungsten 183.84 183.950933 2194D6 1.62 2.36 6 6
75 Re Rhenium 186.207 186.955751 267DAB 1.51 1.9 7 6
76 Os Osmium 190.23 191.961479 266696 1.44 2.2 8 6
77 Ir Iridium 192.217 192.962924 175487 1.41 2.20 9 6
78 Pt Platinum 195.084 194.964774 D0D0E0 1.36 2.28 10 6
79 Au Gold 196.966569 196.966552 FFD123 1.36 2.54 11 6
80 Hg Mercury 200.59 201.970626 B8B8D0 1.32 2.00 12 6
81 Tl Thallium 204.3833 204.974412 A6544D 1.45 1.62 13 6
82 Pb Lead 207.2 207.976636 575961 1.46 2.33 14 6
83 Bi Bismuth 208.9804 208.980383 9E4FB5 1.48 2.02 15 6
84 Po Polonium 209 208.982416 AB5C00 1.40 2.0 16 6
85 At Astatine 210 209.987131 754F45 1.50 2.2 17 6
86 Rn Radon 222 222.0175705 428296 1.50 0 18 6
87 Fr Francium 223 223.0197307 420066 2.60 0.7 1 7
88 Ra Radium 226 226.0254026 007D00 2.21 0.9 2 7
89 Ac Actinium 227 227.027747 70ABFA 2.15 1.1 0 7
90 Th Thorium 232.03806 232.0380504 00BAFF 2.06 1.3 0 7
91 Pa Protactinium 231.03588 231.0358789 00A1FF 2.00 1.5 0 7
92 U Uranium 238.02891 238.0507826 008FFF 1.96 1.38 0 7
93 Np Neptunium 237 237.0481673 0080FF 1.90 1.36 0 7
94 Pu Plutonium 244 244.064198 006BFF 1.87 1.28 0 7
95 Am Americium 243 243.0613727 545CF2 1.80 1.3 0 7
96 Cm Curium 247 247.070347 785CE3 1.69 1.3 0 7
97 Bk Berkelium 247 247.070299 8A4FE3 1.68 1.3 0 7
98 Cf Californium 251 251.07958 A136D4 1.68 1.3 0 7
99 Es Einsteinium 252 252.08297 B31FD4 1.65 1.3 0 7
100 Fm Fermium 257 257.095099 B31FBA 1.67 1.3 0 7
101 Md Mendelevium 258 258.098425 B30DA6 1.73 1.3 0 7
102 No Nobelium 259 259.10102 BD0D87 1.76 1.3 0 7
103 Lr Lawrencium 266 262.10969 C70066 1.61 1.3 3 7
104 Rf Rutherfordium 267 267.12179 CC0059 1.57 0 4 7
105 Db Dubnium 268 268.12567 D1004F 1.49 0 5 7
106 Sg Seaborgium 269 269.12863 D90045 1.43 0 6 7
107 Bh Bohrium 270 270.13336 E00038 1.41 0 7 7
108 Hs Hassium 269 269.13375 E6002E 1.34 0 8 7
109 Mt Meitnerium 278 278.15631 EB0026 1.29 0 9 7
110 Ds Darmstadtium 281 281.16451 EB0026 1.28 0 10 7
111 Rg Roentgenium 282 282.16912 EB0026 1.21 0 11 7
112 Cn Copernicium 285 285.17712 EB0026 1.22 0 12 7
113 Nh Nihonium 286 286.18221 EB0026 1.36 0 13 7
114 Fl Flerovium 289 289.19042 EB0026 1.43 0 14 7
115 Mc Moscovium 290 290.19598 EB0026 1.62 0 15 7
116 Lv Livermorium 293 293.20449 EB0026 1.75 0 16 7
117 Ts Tennessine 294 294.21046 EB0026 1.65 0 17 7
118 Og Oganesson 294 294.21392 EB0026 1.57 0 18 7
`;

export const ELEMENTS: ElementData[] = [];
export const BY_SYMBOL = new Map<string, ElementData>();

for (const line of RAW.trim().split('\n')) {
  const p = line.trim().split(/\s+/);
  const e: ElementData = {
    z: +p[0],
    symbol: p[1],
    name: p[2],
    mass: +p[3],
    mono: +p[4],
    color: '#' + p[5],
    covRadius: +p[6],
    en: +p[7],
    group: +p[8],
    period: +p[9],
  };
  ELEMENTS[e.z] = e;
  BY_SYMBOL.set(e.symbol, e);
}

/** Deuterium / tritium are accepted as element aliases on input. */
export const ISOTOPE_ALIASES: Record<string, { el: string; isotope: number }> = {
  D: { el: 'H', isotope: 2 },
  T: { el: 'H', isotope: 3 },
};

/** Exact masses of selected isotopes, used when an atom carries an explicit mass number. */
export const ISOTOPE_MASS: Record<string, number> = {
  'H1': 1.0078250319, 'H2': 2.0141017778, 'H3': 3.0160492777,
  'C12': 12.0, 'C13': 13.0033548378, 'C14': 14.003241989,
  'N14': 14.0030740048, 'N15': 15.0001088982,
  'O16': 15.99491461956, 'O17': 16.9991317, 'O18': 17.999161,
  'F18': 18.000938, 'F19': 18.99840322,
  'P31': 30.97376163, 'P32': 31.97390727,
  'S32': 31.972071, 'S33': 32.97145876, 'S34': 33.9678669, 'S35': 34.96903216,
  'Cl35': 34.96885268, 'Cl37': 36.96590259,
  'Br79': 78.9183371, 'Br81': 80.9162906,
  'I125': 124.9046302, 'I127': 126.904473, 'I131': 130.9061246,
};

export function element(symbol: string): ElementData | undefined {
  return BY_SYMBOL.get(symbol);
}

export function isElement(symbol: string): boolean {
  return BY_SYMBOL.has(symbol);
}

export function atomicNumber(symbol: string): number {
  return BY_SYMBOL.get(symbol)?.z ?? 0;
}

/** Number of valence electrons for main-group elements (s/p block). Transition metals return group number capped at 12 (d electrons + s). */
export function valenceElectrons(symbol: string): number {
  const e = BY_SYMBOL.get(symbol);
  if (!e) return 0;
  if (e.group >= 13) return e.group - 10;
  if (e.group === 0) return 3;
  return e.group;
}

/** Elements that receive implicit hydrogens when drawn without an explicit H count. */
export const IMPLICIT_H_ELEMENTS = new Set([
  'H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Ge', 'As', 'Se', 'Br', 'Te', 'I', 'At',
]);

/** Allowed (neutral) valences of main-group elements by group and period. */
export function allowedValences(group: number, period: number): number[] {
  switch (group) {
    case 1: return [1];
    case 2: return [2];
    case 13: return [3];
    case 14: return [4];
    case 15: return period >= 3 ? [3, 5] : [3];
    case 16: return period >= 3 ? [2, 4, 6] : [2];
    case 17: return period >= 3 ? [1, 3, 5, 7] : [1];
    case 18: return [0];
    default: return [];
  }
}

/** Darker per-element colours suitable for 2D drawings on a white background. */
export const LABEL_COLORS: Record<string, string> = {
  H: '#404040', C: '#000000', N: '#2140d9', O: '#e00000', F: '#2f9e2f', Cl: '#1b9e1b',
  Br: '#a52a2a', I: '#8a008a', S: '#b38600', P: '#e06c00', B: '#c05a5a', Si: '#8a6d3b',
  Se: '#c47800', Na: '#7b3fbf', K: '#6a2fb0', Li: '#9b59d0', Mg: '#3c8f00', Ca: '#2f8f00',
  Fe: '#c0501f', Cu: '#a8642a', Zn: '#5a5e8f', Al: '#8f7676', Sn: '#4f6b6b', Pt: '#707088',
  Pd: '#005266', Au: '#b38f00', Ag: '#7a7a7a', Hg: '#6a6a88', As: '#8a4fbf', Ti: '#707378',
};

export function labelColor(symbol: string): string {
  return LABEL_COLORS[symbol] ?? '#333333';
}

/** Rows of the periodic table for the element picker: [symbol | '' (gap)] laid out 18 columns wide. */
export const PERIODIC_LAYOUT: (string | '')[][] = [
  ['H', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'He'],
  ['Li', 'Be', '', '', '', '', '', '', '', '', '', '', 'B', 'C', 'N', 'O', 'F', 'Ne'],
  ['Na', 'Mg', '', '', '', '', '', '', '', '', '', '', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
  ['K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'],
  ['Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe'],
  ['Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'],
  ['Fr', 'Ra', 'Ac', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'],
  ['', '', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', '', ''],
  ['', '', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', '', ''],
];
