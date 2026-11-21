// Liquid properties at the melting point (values of Table 1 of the accompanying
// article; Slag6: CaO-Al2O3-SiO2 slag), and water at 20 C as a familiar reference.
// density [kg/m^3], viscosity [mPa s], surface tension [N/m], T [K]
export const MATERIALS = [
  { key: 'Ni', T: 1728, rho: 7890, mu: 4.7, sigma: 1.739 },
  { key: 'Pt', T: 2041, rho: 19200, mu: 6.59, sigma: 1.798 },
  { key: 'Ti', T: 1941, rho: 4170, mu: 4.4, sigma: 1.557 },
  { key: 'V', T: 2183, rho: 5460, mu: 4.3, sigma: 1.935 },
  { key: 'Zr', T: 2128, rho: 6210, mu: 4.7, sigma: 1.5 },
  { key: 'Nb', T: 2750, rho: 7730, mu: 4.5, sigma: 1.937 },
  { key: 'Mo', T: 2896, rho: 9110, mu: 5.6, sigma: 2.29 },
  { key: 'Ru', T: 2607, rho: 10750, mu: 6.1, sigma: 2.256 },
  { key: 'Hf', T: 2504, rho: 11820, mu: 5.2, sigma: 1.614 },
  { key: 'Ta', T: 3290, rho: 14750, mu: 8.6, sigma: 2.154 },
  { key: 'W', T: 3695, rho: 16430, mu: 6.9, sigma: 2.477 },
  { key: 'Re', T: 3441, rho: 18650, mu: 7.9, sigma: 2.71 },
  { key: 'Ir', T: 2719, rho: 19870, mu: 7.0, sigma: 2.241 },
  { key: 'Si', T: 1687, rho: 2520, mu: 0.574, sigma: 0.732 },
  { key: 'Ti-6Al-4V', T: 1933, rho: 4087, mu: 4.0, sigma: 1.493 },
  { key: 'Al2O3', T: 2327, rho: 3000, mu: 30, sigma: 0.65 },
  { key: 'Tb2O3', T: 2683, rho: 7450, mu: 27, sigma: 0.73 },
  { key: 'Slag6', T: 1792, rho: 2824, mu: 250, sigma: 0.548 },
  { key: 'Water', T: 293, rho: 998, mu: 1.0, sigma: 0.0728 },
];
