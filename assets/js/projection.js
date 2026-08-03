/**
 * projection.js — 1D projection machinery for the TrigPoly Fourier maps
 * produced by AngularExpression.compute().
 *
 * Extracted verbatim from tools/angular-projection.html so the same θ/φ
 * marginal formulas are shared by the Angular Projection tool and the
 * Dalitz Plot tool (which uses the θ₁-integrated angular distribution).
 *
 * Dependencies: surd.js (Surd, SurdSum)
 */

// ══════════════════════════════════════════════════════════
// INTEGRAL TABLES
// ══════════════════════════════════════════════════════════

// ── Integral tables for θ variables ──
// I_c(k) = ∫₀^π cos(k·θ/2) sin θ dθ  (as rational, no π)
// I_s(k) = ∫₀^π sin(k·θ/2) sin θ dθ  (as {num, den, piExp})
var _THETA_COS_INT = {
  0: { n: 2, d: 1 },
  1: { n: 4, d: 3 },
  2: { n: 0, d: 1 },
  3: { n: -4, d: 5 },
  4: { n: -2, d: 3 },
  5: { n: -4, d: 21 },
  6: { n: 0, d: 1 },
  7: { n: -4, d: 45 },
  8: { n: -2, d: 15 },
  9: { n: -4, d: 77 },
  10: { n: 0, d: 1 }
};

var _THETA_SIN_INT = {
  0: { n: 0, d: 1, pi: 0 },
  1: { n: 4, d: 3, pi: 0 },
  2: { n: 1, d: 2, pi: 1 },    // π/2
  3: { n: 4, d: 5, pi: 0 },
  4: { n: 0, d: 1, pi: 0 },
  5: { n: -4, d: 21, pi: 0 },
  6: { n: 0, d: 1, pi: 0 },
  7: { n: 4, d: 45, pi: 0 },
  8: { n: 0, d: 1, pi: 0 },
  9: { n: -4, d: 77, pi: 0 },
  10: { n: 0, d: 1, pi: 0 }
};

// Extended for larger k via formula
function _lookupInt(table, k) {
  if (table[k]) return table[k];
  // Compute on-the-fly for k > 10 (rare in practice)
  if (table === _THETA_COS_INT) {
    return _computeThetaCosInt(k);
  } else {
    return _computeThetaSinInt(k);
  }
}

function _computeThetaCosInt(k) {
  // I_c(k) = (J(k+2) - J(k-2)) / 2
  // J(n) = 2/n * (1 - cos(nπ/2)) for n>0, J(-n) = -J(n), J(0) = 0
  var jp = _J(k + 2);
  var jm = _J(k - 2);
  if (jp === 0 && jm === 0) return { n: 0, d: 1 };
  // Convert to rational
  var num = (jp - jm);
  var den = 2;
  // Reduce by gcd
  var g = _gcd(Math.abs(num), den);
  return { n: num / g, d: den / g };
}

function _computeThetaSinInt(k) {
  // I_s(k) = (K(k-2) - K(k+2)) / 2
  // K(n) = 2/n * sin(nπ/2) for n≠0, K(0) = π
  var km = k - 2;
  var kp = k + 2;
  if (km === 0) {
    // K(0) = π — need to handle π separately
    var Kp = _K(kp);
    return { n: -Kp, d: 2, pi: 1 };  // (π - Kp) / 2
  }
  if (kp === 0) {
    var Km = _K(km);
    return { n: Km, d: 2, pi: 1 };  // (Km - π) / 2
  }
  var Kv = _K(km) - _K(kp);
  if (Kv === 0) return { n: 0, d: 1, pi: 0 };
  // Kv is integer
  var g = _gcd(Math.abs(Kv), 2);
  return { n: Kv / g, d: 2 / g, pi: 0 };
}

function _J(n) {
  if (n === 0) return 0;
  if (n < 0) return -_J(-n);
  // J(n) = 2/n * (1 - cos(nπ/2))
  var cosVal;
  if (n % 2 === 0) {
    // Even n: cos(nπ/2) = (-1)^(n/2)
    cosVal = (n / 2) % 2 === 0 ? 1 : -1;
  } else {
    cosVal = 0;
  }
  if (cosVal === 0) return 2 / n;
  if (cosVal === 1) return 0;
  // cosVal === -1, 1 - (-1) = 2
  return 4 / n;
}

function _K(n) {
  if (n === 0) return 0; // π case handled separately
  if (n < 0) return _K(-n); // cos is even
  // K(n) = 2/n * sin(nπ/2)
  if (n % 2 === 0) return 0;
  var rem = n % 4;
  return rem === 1 ? (2 / n) : (-2 / n);
}

// Precompute tables up to k=30
(function _buildIntTables() {
  for (var k = 11; k <= 30; k++) {
    if (!_THETA_COS_INT[k]) _THETA_COS_INT[k] = _computeThetaCosInt(k);
    if (!_THETA_SIN_INT[k]) _THETA_SIN_INT[k] = _computeThetaSinInt(k);
  }
})();

function _gcd(a, b) {
  while (b !== 0) { var t = b; b = a % b; a = t; }
  return a;
}

/**
 * Check whether a variable name is a θ type (half-angle Fourier).
 */
function _isThetaVar(name) {
  return name.indexOf('theta_') === 0;
}

/**
 * Sort variable names in reverse order: highest index first, chi last.
 * This puts the physically deepest cascade angles at the top.
 */
function _reverseVarSort(a, b) {
  if (a === 'chi') return -1;
  if (b === 'chi') return 1;
  // Extract numeric index from "theta_2" or "phi_1"
  var aIdx = parseInt(a.replace(/\D/g, ''), 10);
  var bIdx = parseInt(b.replace(/\D/g, ''), 10);
  if (aIdx !== bIdx) return bIdx - aIdx;  // descending index
  // Same index: theta before phi (theta deeper in physical sense)
  if (a.indexOf('theta') !== -1 && b.indexOf('phi') !== -1) return -1;
  if (a.indexOf('phi') !== -1 && b.indexOf('theta') !== -1) return 1;
  return a.localeCompare(b);
}

/**
 * Get the integral coefficient and π exponent for a Fourier factor {func, k}
 * over its full range.
 * Returns { n: int, d: int, pi: int }
 *   meaning the integral = (n/d) × π^pi
 *   n=0 means the integral is zero.
 */
function _integralFactor(varName, func, k) {
  if (k === 0) {
    // cos(0) = 1, sin(0) = 0
    if (func === 'sin') return { n: 0, d: 1, pi: 0 };
    // cos(0) = 1 → ∫ 1 d(measure)
    return _emptyIntegral(varName);
  }
  if (_isThetaVar(varName)) {
    if (func === 'cos') {
      var ci = _lookupInt(_THETA_COS_INT, k);
      return { n: ci.n, d: ci.d, pi: 0 };
    } else {
      var si = _lookupInt(_THETA_SIN_INT, k);
      return { n: si.n, d: si.d, pi: si.pi || 0 };
    }
  } else {
    // φ/χ variable with half-angle: func(k·var/2) over [-π, π]
    // cos(k·var/2): ∫_{-π}^{π} cos(k·var/2) d(var) = 0 for even k ≠ 0, ±4/k for odd k
    if (func === 'cos') {
      if (k % 2 === 0) return { n: 0, d: 1, pi: 0 };
      // odd k: sin(kπ/2) = 1 for k≡1 mod4, -1 for k≡3 mod4
      var sign = (k % 4 === 1) ? 1 : -1;
      var g = _gcd(4, Math.abs(k));
      return { n: sign * 4 / g, d: Math.abs(k) / g, pi: 0 };
    }
    // sin(k·var/2): ∫_{-π}^{π} sin(k·var/2) d(var) = 0 for all k ≠ 0 (symmetric range)
    return { n: 0, d: 1, pi: 0 };
  }
}

/**
 * Integral over a variable that does NOT appear in the Fourier factor.
 * Returns { n, d, pi } where the integral = (n/d) × π^pi.
 */
function _emptyIntegral(varName) {
  if (_isThetaVar(varName)) {
    // ∫ sin θ dθ from 0 to π = 2
    return { n: 2, d: 1, pi: 0 };
  } else {
    // ∫ dv from -π to π = 2π
    return { n: 2, d: 1, pi: 1 };
  }
}

/**
 * Combine two integral factors: multiply both n/d and add π exponents.
 */
function _mulIntFactors(a, b) {
  if (a.n === 0 || b.n === 0) return { n: 0, d: 1, pi: 0 };
  var num = a.n * b.n;
  var den = a.d * b.d;
  var pi = (a.pi || 0) + (b.pi || 0);
  var g = _gcd(Math.abs(num), den);
  return { n: num / g, d: den / g, pi: pi };
}

/**
 * Parse a Fourier key into an array of {name, func, k}.
 */
function _parseFourierKey(key) {
  if (key === '1') return [];
  var parts = key.split(',');
  var result = [];
  for (var i = 0; i < parts.length; i++) {
    var m = parts[i].match(/([\w\d_]+):(\w+)\|(-?[\d.]+)/);
    if (m) result.push({ name: m[1], func: m[2], k: parseFloat(m[3]) });
  }
  return result;
}

/**
 * Build a canonical Fourier key from sorted factor descriptors.
 */
function _makeFourierKey(factors) {
  if (!factors || factors.length === 0) return '1';
  // Sort by name, then k, then func
  var sorted = factors.slice().sort(function(a, b) {
    if (a.name !== b.name) return a.name < b.name ? -1 : 1;
    if (a.k !== b.k) return a.k - b.k;
    return a.func < b.func ? -1 : (a.func > b.func ? 1 : 0);
  });
  return sorted.map(function(f) { return f.name + ':' + f.func + '|' + f.k; }).join(',');
}

/**
 * Extract the unique variable names from a Fourier map (the keys of lsMap).
 * Returns: array of variable name strings.
 */
function _extractVarNames(fourierMap) {
  var names = {};
  for (var key in fourierMap) {
    if (key === '1') continue;
    var factors = _parseFourierKey(key);
    for (var i = 0; i < factors.length; i++) {
      names[factors[i].name] = true;
    }
  }
  return Object.keys(names).sort();
}

// ══════════════════════════════════════════════════════════
// SurdSum helpers
// ══════════════════════════════════════════════════════════

/**
 * Scale all terms in a SurdSum by rational n/d.
 */
function _surSumScale(sum, n, d) {
  if (n === 0 || sum.isEmpty()) return new SurdSum();
  var result = new SurdSum();
  var terms = sum.terms();
  for (var i = 0; i < terms.length; i++) {
    var s = Surd.scale(terms[i], n, d);
    if (!s.isZero()) result.add(s);
  }
  return result;
}

/**
 * Add all terms from 'other' SurdSum into 'target' SurdSum.
 * Returns target for chaining.
 */
function _surSumAdd(target, other) {
  if (other.isEmpty()) return target;
  var terms = other.terms();
  for (var i = 0; i < terms.length; i++) {
    target.add(terms[i]);
  }
  return target;
}

/**
 * Extract the rational value of a SurdSum (all terms must have r=1, no sqrt).
 * Returns { n: num, d: den } or null if terms have sqrt roots.
 */
function _surSumToRational(sum) {
  if (sum.isEmpty()) return null;
  var num = 0, den = 1;
  var terms = sum.terms();
  for (var i = 0; i < terms.length; i++) {
    var t = terms[i];
    if (t.r !== 1) return null;  // has sqrt, can't use rational
    var tn = t.s * t.p;
    var td = t.q;
    num = num * td + tn * den;
    den = den * td;
    var g = _gcd(Math.abs(num), den);
    num /= g; den /= g;
  }
  return { n: num, d: den };
}

/**
 * Convert a Surd to a floating-point number.
 */
function _surdFloat(s) {
  var v = s.s * s.p / s.q;
  if (s.r > 1) v *= Math.sqrt(s.r);
  return v;
}

/**
 * Convert a SurdSum to a floating-point number.
 */
function _surSumFloat(sum) {
  var t = 0;
  var terms = sum.terms();
  for (var i = 0; i < terms.length; i++) t += _surdFloat(terms[i]);
  return t;
}

// ══════════════════════════════════════════════════════════
// PROJECTION COMPUTATION
// ══════════════════════════════════════════════════════════

/**
 * Build the full list of variable names from nDecays and phiCombine.
 */
function _allVarNames(nDecays, phiCombine) {
  var names = [];
  // Reverse order: highest index first (phi_0 / theta_0 are least interesting)
  for (var i = nDecays - 1; i >= 0; i--) {
    names.push('theta_' + i);
  }
  for (var i = nDecays - 1; i >= 0; i--) {
    if (phiCombine && (i === phiCombine.fixIdx || i === phiCombine.chiIdx)) continue;
    names.push('phi_' + i);
  }
  if (phiCombine) {
    names.push('chi');
  }
  return names;
}

/**
 * Compute the total integral of |T|² over all angles: ∫ I(θ,φ) dΩ.
 * This gives the overall scale factor that changes with helicity filters.
 * Returns { piExp: SurdSum } — e.g. { 1: SurdSum(4) } = 4π.
 */
function computeTotalIntegral(fourierMap, allVars) {
  var groups = {};
  for (var key in fourierMap) {
    var sum = fourierMap[key];
    if (sum.isEmpty()) continue;
    var factors = _parseFourierKey(key);
    var totalInt = {n: 1, d: 1, pi: 0};
    var kill = false;
    for (var avi = 0; avi < allVars.length; avi++) {
      var av = allVars[avi];
      var found = false;
      for (var fi = 0; fi < factors.length; fi++) {
        if (factors[fi].name === av) {
          var intF = _integralFactor(av, factors[fi].func, factors[fi].k);
          if (intF.n === 0) { kill = true; break; }
          totalInt = _mulIntFactors(totalInt, intF);
          found = true; break;
        }
      }
      if (kill) break;
      if (!found) totalInt = _mulIntFactors(totalInt, _emptyIntegral(av));
    }
    if (kill || totalInt.n === 0) continue;
    var scaled = _surSumScale(sum, totalInt.n, totalInt.d);
    if (scaled.isEmpty()) continue;
    var piExp = totalInt.pi;
    if (!groups[piExp]) groups[piExp] = new SurdSum();
    _surSumAdd(groups[piExp], scaled);
  }
  return groups;
}

/**
 * Compute 1D projections of a Fourier map onto each variable.
 *
 * @param {Object} fourierMap  { key: SurdSum } from AngularExpression
 * @param {number} nDecays     Number of decay vertices
 * @param {Object|null} phiCombine  J=0 combine info
 * @return {Object}  { varName: { 0: {key: SurdSum}, 1: {key: SurdSum}, ... } }
 *   where the outer group key is the π exponent.
 */
function computeProjections(fourierMap, nDecays, phiCombine) {
  var allVars = _allVarNames(nDecays, phiCombine);
  var results = {};

  for (var vi = 0; vi < allVars.length; vi++) {
    var projVar = allVars[vi];
    // Group by π exponent: { piExp: { fourierKey: SurdSum } }
    var groups = {};

    for (var key in fourierMap) {
      var sum = fourierMap[key];
      if (sum.isEmpty()) continue;

      var factors = _parseFourierKey(key);
      var totalInt = { n: 1, d: 1, pi: 0 };
      var projFactor = null;
      var otherFactors = [];
      var kill = false;

      for (var fi = 0; fi < factors.length; fi++) {
        var f = factors[fi];
        if (f.name === projVar) {
          projFactor = f;
        } else {
          var intF = _integralFactor(f.name, f.func, f.k);
          if (intF.n === 0) { kill = true; break; }
          totalInt = _mulIntFactors(totalInt, intF);
        }
      }
      if (kill) continue;

      // Integrate over non-appearing variables (that are NOT the projection variable)
      for (var avi = 0; avi < allVars.length; avi++) {
        var av = allVars[avi];
        if (av === projVar) continue;
        var found = false;
        for (var fi = 0; fi < factors.length; fi++) {
          if (factors[fi].name === av) { found = true; break; }
        }
        if (!found) {
          var emptyI = _emptyIntegral(av);
          totalInt = _mulIntFactors(totalInt, emptyI);
        }
      }

      if (totalInt.n === 0) continue;

      // Scale SurdSum by totalInt.n/totalInt.d (integral over other variables)
      // No further normalization — the projection is the actual marginal I(x) = ∫ I(x,others) d(others),
      // so ∫ I_proj(x) dx over the full variable range equals the total integral.
      var piExp = totalInt.pi;
      if (!groups[piExp]) groups[piExp] = {};

      // Build the resulting Fourier key (only the projection variable's factor, or constant)
      var outKey;
      if (projFactor) {
        outKey = _makeFourierKey([projFactor]);
      } else {
        outKey = '1';
      }

      // Scale SurdSum by n/d
      var scaled = _surSumScale(sum, totalInt.n, totalInt.d);
      if (scaled.isEmpty()) continue;

      if (!groups[piExp][outKey]) groups[piExp][outKey] = new SurdSum();
      _surSumAdd(groups[piExp][outKey], scaled);
      if (groups[piExp][outKey].isEmpty()) delete groups[piExp][outKey];
    }

    results[projVar] = groups;
  }

  return results;
}

/**
 * Evaluate a 1D projection result (one variable's groups from
 * computeProjections) at a numeric value of that variable.
 * Each group is a Fourier map in the half-angle basis func(k·x/2) with an
 * overall π^piExp factor.
 *
 * @param {Object} projGroups  projections[varName] = { piExp: { key: SurdSum } }
 * @param {number} x           value of the projection variable (e.g. θ)
 * @returns {number}
 */
function evaluateProjection(projGroups, x) {
  var total = 0;
  for (var piExp in projGroups) {
    var pi = parseFloat(piExp);
    var piScale = (pi === 0) ? 1 : Math.pow(Math.PI, pi);
    var groups = projGroups[piExp];
    for (var key in groups) {
      var sum = groups[key];
      if (sum.isEmpty()) continue;
      var coeff = _surSumFloat(sum);
      if (coeff === 0) continue;
      var fac = 1;
      if (key !== '1') {
        var factors = _parseFourierKey(key);
        for (var i = 0; i < factors.length; i++) {
          var f = factors[i];
          fac *= (f.func === 'cos' ? Math.cos(f.k * x / 2) : Math.sin(f.k * x / 2));
        }
      }
      total += coeff * piScale * fac;
    }
  }
  return total;
}

// Node/browser export
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    computeProjections: computeProjections,
    computeTotalIntegral: computeTotalIntegral,
    evaluateProjection: evaluateProjection,
    _allVarNames: _allVarNames,
    _parseFourierKey: _parseFourierKey,
    _surSumFloat: _surSumFloat
  };
}
