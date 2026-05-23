/**
 * Wigner d-function Calculator
 * Uses Algebrite for exact symbolic computation
 *
 * Formula:
 *   d^{j}_{m1,m2}(beta) = sum_{l=0}^{2j} w_l^{(j,m1,m2)} * sin^l(beta/2) * cos^{2j-l}(beta/2)
 *
 * where the weight w_l^{(j,m1,m2)} is:
 *   w_l = (-1)^{m1-m2+k} * sqrt((j+m1)!(j-m1)!(j+m2)!(j-m2)!) /
 *         ((j-m1-k)!(j+m2-k)!(m1-m2+k)!k!)
 *
 * with k = (l + m2 - m1)/2, and k must satisfy:
 *   max(0, m2-m1) <= k <= min(j-m1, j+m2)
 *
 * Otherwise w_l = 0.
 */

// ============================================================================
// INPUT PARSING
// ============================================================================

/**
 * Parse a string input to a number
 * Accepts: "1/2", "3/2", "-1/2", "0.5", "1", "2", etc.
 * Returns: {numerator, denominator, value, algebriteForm}
 */
function parseQuantumNumber(input) {
  const s = String(input).trim();

  // Fraction format "a/b" — parse numerator/denominator as integers
  const slashIdx = s.indexOf('/');
  if (slashIdx !== -1) {
    const num = parseInt(s.substring(0, slashIdx), 10);
    const den = parseInt(s.substring(slashIdx + 1), 10);
    if (isNaN(num) || isNaN(den) || den === 0) {
      throw new Error(`Invalid quantum number: ${s}`);
    }
    return {
      numerator: num,
      denominator: den,
      value: num / den,
      algebriteForm: `(${num}/${den})`
    };
  }

  // Integer-only (no decimals — use fraction syntax for half-integers)
  const n = parseInt(s, 10);
  if (isNaN(n)) {
    throw new Error(`Invalid quantum number: ${s}. Use fraction syntax e.g. "1/2" for half-integers.`);
  }
  return {
    numerator: n,
    denominator: 1,
    value: n,
    algebriteForm: n < 0 ? `(${n})` : String(n)
  };
}

/**
 * Validate that a quantum number is valid (integer or half-integer).
 */
function isValidQuantumNumber(parsed, allowNegative) {
  if (allowNegative === undefined) allowNegative = false;
  const twiceVal = parsed.value * 2;
  if (Math.abs(twiceVal - Math.round(twiceVal)) > 1e-10) return false;
  if (!allowNegative && parsed.value < 0) return false;
  return true;
}

// ============================================================================
// WIGNER d-FUNCTION CALCULATION
// ============================================================================

/**
 * Compute a single Wigner d-function element d^{j}_{m1,m2}(beta).
 *
 * @param {string} jStr - angular momentum j (e.g. "1/2", "1", "3/2")
 * @param {string} m1Str - projection m1
 * @param {string} m2Str - projection m2
 * @param {string} betaStr - angle beta (e.g. "pi/2", "0.5", "pi/4"). If empty/null, treat as symbolic.
 * @returns {{decimal: number, symbolic: string, latex: string, error: string}}
 */
function computeWignerD(jStr, m1Str, m2Str, betaStr) {
  if (typeof Algebrite === 'undefined') {
    return { error: 'Algebrite library is not loaded. Cannot perform calculation.' };
  }

  // Parse inputs
  var pj, pm1, pm2;
  try {
    pj = parseQuantumNumber(jStr);
    pm1 = parseQuantumNumber(m1Str);
    pm2 = parseQuantumNumber(m2Str);
  } catch (e) {
    return { error: e.message };
  }

  // Validate
  if (!isValidQuantumNumber(pj)) {
    return { error: 'j = ' + jStr + ' is not a valid quantum number (must be integer or half-integer >= 0)' };
  }
  if (!isValidQuantumNumber(pm1, true)) {
    return { error: 'm1 = ' + m1Str + ' is not a valid quantum number (must be integer or half-integer)' };
  }
  if (!isValidQuantumNumber(pm2, true)) {
    return { error: 'm2 = ' + m2Str + ' is not a valid quantum number (must be integer or half-integer)' };
  }

  var jVal = pj.value;
  var m1Val = pm1.value;
  var m2Val = pm2.value;

  // Check |m| <= j
  if (Math.abs(m1Val) - jVal > 1e-10) {
    return { error: '|m1| <= j violated: |' + m1Val + '| > ' + jVal };
  }
  if (Math.abs(m2Val) - jVal > 1e-10) {
    return { error: '|m2| <= j violated: |' + m2Val + '| > ' + jVal };
  }

  // Check j +/- m are integers
  var jPlusM1 = Math.round(jVal + m1Val);
  var jMinusM1 = Math.round(jVal - m1Val);
  var jPlusM2 = Math.round(jVal + m2Val);
  var jMinusM2 = Math.round(jVal - m2Val);

  if (Math.abs((jVal + m1Val) - jPlusM1) > 1e-10) {
    return { error: 'j + m1 = ' + (jVal + m1Val) + ' is not an integer' };
  }
  if (Math.abs((jVal + m2Val) - jPlusM2) > 1e-10) {
    return { error: 'j + m2 = ' + (jVal + m2Val) + ' is not an integer' };
  }

  // Precompute the factorial prefactor
  // sqrt((j+m1)!(j-m1)!(j+m2)!(j-m2)!)
  var prefactorExpr = 'sqrt((' + jPlusM1 + ')! * (' + jMinusM1 + ')! * (' + jPlusM2 + ')! * (' + jMinusM2 + ')! )';

  // Build the sum over l = 0..2j
  var twoJ = Math.round(2 * jVal);
  var terms = [];

  for (var l = 0; l <= twoJ; l++) {
    // k = (l + m2 - m1) / 2
    var k = (l + m2Val - m1Val) / 2;
    var kInt = Math.round(k);

    // Check if k is integer and within bounds
    if (Math.abs(k - kInt) > 1e-10) continue;
    if (kInt < Math.max(0, m2Val - m1Val) - 1e-10) continue;
    if (kInt > Math.min(jMinusM1, jPlusM2) + 1e-10) continue;

    // Sign: (-1)^{m1-m2+k}
    var sign = ((Math.round(m1Val - m2Val) + kInt) % 2 === 0) ? '' : '-';

    // Denominator: (j-m1-k)! (j+m2-k)! (m1-m2+k)! k!
    var denom = '(' + (jMinusM1 - kInt) + ')! * (' + (jPlusM2 - kInt) + ')! * (' + (Math.round(m1Val - m2Val) + kInt) + ')! * (' + kInt + ')!';

    // Weight w_l
    var weightExpr = sign + prefactorExpr + ' / (' + denom + ')';

    // sin^l(beta/2) * cos^{2j-l}(beta/2)
    var sinPart, cosPart;
    if (betaStr && betaStr.trim() !== '') {
      // Numeric beta: evaluate sin(beta/2) and cos(beta/2)
      sinPart = 'sin(' + betaStr + '/2)^' + l;
      cosPart = 'cos(' + betaStr + '/2)^' + (twoJ - l);
    } else {
      // Symbolic beta
      sinPart = 'sin(beta/2)^' + l;
      cosPart = 'cos(beta/2)^' + (twoJ - l);
    }

    var termExpr = weightExpr + ' * ' + sinPart + ' * ' + cosPart;
    terms.push(termExpr);
  }

  if (terms.length === 0) {
    return { decimal: 0, symbolic: '0', latex: '0' };
  }

  var sumExpr = terms.join(' + ').replace(/\+ -/g, '- ');

  try {
    var simplified = Algebrite.run('simplify(' + sumExpr + ')');

    // Get decimal value
    var decimalStr = Algebrite.run('float(' + simplified + ')');
    var decimalVal = parseFloat(decimalStr);

    // Get LaTeX
    var latex = '';
    try {
      latex = Algebrite.run('printlatex(' + simplified + ')').trim();
    } catch (_) {
      latex = '';
    }

    // If printlatex failed or returned empty, build LaTeX from symbolic
    if (!latex) {
      latex = simplified.trim()
        .replace(/\*/g, ' \\cdot ')
        .replace(/\^/g, '^{') + '}';
    }

    // Algebrite's printlatex outputs plain "cos", "sin", "beta" — fix to LaTeX commands
    latex = latex.replace(/(?<!\\)cos/g, '\\cos');
    latex = latex.replace(/(?<!\\)sin/g, '\\sin');
    latex = latex.replace(/(?<!\\)beta/g, '\\beta');

    return {
      decimal: decimalVal,
      symbolic: simplified.trim(),
      latex: latex
    };
  } catch (e) {
    return { error: 'Algebrite computation error: ' + (e.message || e) };
  }
}

/**
 * Compute the full Wigner d-matrix for a given j and beta.
 * Returns a 2D array of results indexed by m1, m2.
 *
 * @param {string} jStr
 * @param {string} betaStr
 * @returns {{matrix: Array<Array>, mValues: Array<number>, error: string}}
 */
/**
 * Convert a numeric quantum-number value to fraction string.
 * e.g. 0.5 → "1/2", -0.5 → "-1/2", 1 → "1", -1 → "-1"
 */
function numberToFractionStr(val) {
  var twice = Math.round(val * 2);
  if (twice % 2 === 0) return String(twice / 2);
  var sign = twice < 0 ? '-' : '';
  return sign + Math.abs(twice) + '/2';
}

function computeWignerDMatrix(jStr, betaStr) {
  if (typeof Algebrite === 'undefined') {
    return { error: 'Algebrite library is not loaded.' };
  }

  var pj;
  try {
    pj = parseQuantumNumber(jStr);
  } catch (e) {
    return { error: e.message };
  }

  if (!isValidQuantumNumber(pj)) {
    return { error: 'j = ' + jStr + ' is not a valid quantum number' };
  }

  var jVal = pj.value;
  var twoJ = Math.round(2 * jVal);

  // Generate m values as fraction strings (j, j-1, ..., -j)
  var mValues = [];
  for (var i = 0; i <= twoJ; i++) {
    mValues.push(numberToFractionStr(jVal - i));
  }

  var matrix = [];
  for (var row = 0; row < mValues.length; row++) {
    var rowData = [];
    for (var col = 0; col < mValues.length; col++) {
      var result = computeWignerD(jStr, mValues[row], mValues[col], betaStr);
      rowData.push(result);
    }
    matrix.push(rowData);
  }

  return { matrix: matrix, mValues: mValues };
}

// ============================================================================
// SANITY CHECKS
// ============================================================================

function runWignerDSanityChecks() {
  console.log('%c=== Wigner d-Function Sanity Checks ===', 'color: #06b6d4; font-weight: bold; font-size: 14px;');
  console.log('');

  if (typeof Algebrite === 'undefined') {
    console.log('%c✗ Algebrite library not loaded — cannot run tests', 'color: #ef4444; font-weight: bold;');
    console.log('');
    return { passed: 0, failed: 1 };
  }

  var tests = [
    // j=1/2: d^{1/2}_{1/2,1/2}(beta) = cos(beta/2)
    { j: '1/2', m1: '1/2', m2: '1/2', beta: '', expectedSymbolic: 'cos(1/2*beta)', description: 'd^{1/2}_{1/2,1/2}(beta) = cos(beta/2)' },
    // j=1/2: d^{1/2}_{1/2,-1/2}(beta) = -sin(beta/2)
    { j: '1/2', m1: '1/2', m2: '-1/2', beta: '', expectedSymbolic: '-sin(1/2*beta)', description: 'd^{1/2}_{1/2,-1/2}(beta) = -sin(beta/2)' },
    // j=1/2: d^{1/2}_{-1/2,1/2}(beta) = sin(beta/2)
    { j: '1/2', m1: '-1/2', m2: '1/2', beta: '', expectedSymbolic: 'sin(1/2*beta)', description: 'd^{1/2}_{-1/2,1/2}(beta) = sin(beta/2)' },
    // j=1/2: d^{1/2}_{-1/2,-1/2}(beta) = cos(beta/2)
    { j: '1/2', m1: '-1/2', m2: '-1/2', beta: '', expectedSymbolic: 'cos(1/2*beta)', description: 'd^{1/2}_{-1/2,-1/2}(beta) = cos(beta/2)' },
    // j=1: d^{1}_{1,1}(beta) = cos^2(beta/2)
    { j: '1', m1: '1', m2: '1', beta: '', expectedSymbolic: 'cos(1/2*beta)^2', description: 'd^{1}_{1,1}(beta) = cos^2(beta/2)' },
    // j=1: d^{1}_{1,0}(beta) = -1/sqrt(2) * sin(beta)
    { j: '1', m1: '1', m2: '0', beta: '', expectedSymbolic: '-2^(1/2)*sin(1/2*beta)*cos(1/2*beta)', description: 'd^{1}_{1,0}(beta) = -sin(beta)/sqrt(2)' },
    // j=1: d^{1}_{1,-1}(beta) = sin^2(beta/2)
    { j: '1', m1: '1', m2: '-1', beta: '', expectedSymbolic: 'sin(1/2*beta)^2', description: 'd^{1}_{1,-1}(beta) = sin^2(beta/2)' },
    // j=1: d^{1}_{0,0}(beta) = cos(beta)
    { j: '1', m1: '0', m2: '0', beta: '', expectedSymbolic: 'cos(beta)', description: 'd^{1}_{0,0}(beta) = cos(beta)' },
    // Numeric test: j=1/2, m1=1/2, m2=1/2, beta=pi -> cos(pi/2) = 0
    { j: '1/2', m1: '1/2', m2: '1/2', beta: 'pi', expectedDecimal: 0, description: 'd^{1/2}_{1/2,1/2}(pi) = 0' },
    // Numeric test: j=1/2, m1=1/2, m2=-1/2, beta=pi -> -sin(pi/2) = -1
    { j: '1/2', m1: '1/2', m2: '-1/2', beta: 'pi', expectedDecimal: -1, description: 'd^{1/2}_{1/2,-1/2}(pi) = -1' },
    // Numeric test: j=1, m1=0, m2=0, beta=pi/2 -> cos(pi/2) = 0
    { j: '1', m1: '0', m2: '0', beta: 'pi/2', expectedDecimal: 0, description: 'd^{1}_{0,0}(pi/2) = 0' },
    // Numeric test: j=1, m1=1, m2=1, beta=pi/2 -> cos^2(pi/4) = 1/2
    { j: '1', m1: '1', m2: '1', beta: 'pi/2', expectedDecimal: 0.5, description: 'd^{1}_{1,1}(pi/2) = 1/2' },
  ];

  var passed = 0;
  var failed = 0;

  tests.forEach(function(test, index) {
    var result = computeWignerD(test.j, test.m1, test.m2, test.beta);

    if (result.error) {
      console.log('%c✗ Test ' + (index + 1) + ': ' + test.description, 'color: #ef4444;');
      console.log('  Error: ' + result.error);
      failed++;
      return;
    }

    var ok = false;
    if (test.expectedDecimal !== undefined) {
      ok = Math.abs(result.decimal - test.expectedDecimal) < 1e-10;
    } else if (test.expectedSymbolic !== undefined) {
      // Compare symbolic forms loosely
      var sym = result.symbolic.replace(/\s+/g, '');
      var exp = test.expectedSymbolic.replace(/\s+/g, '');
      ok = sym === exp || sym.indexOf(exp) !== -1 || exp.indexOf(sym) !== -1;
    }

    if (ok) {
      console.log('%c✓ Test ' + (index + 1) + ': ' + test.description, 'color: #10b981;');
      console.log('  Symbolic: ' + result.symbolic);
      console.log('  Decimal: ' + result.decimal.toFixed(12));
      passed++;
    } else {
      console.log('%c✗ Test ' + (index + 1) + ': ' + test.description, 'color: #ef4444;');
      console.log('  Symbolic: ' + result.symbolic);
      console.log('  Decimal: ' + result.decimal.toFixed(12));
      if (test.expectedSymbolic) console.log('  Expected symbolic: ' + test.expectedSymbolic);
      if (test.expectedDecimal !== undefined) console.log('  Expected decimal: ' + test.expectedDecimal);
      failed++;
    }
    console.log('');
  });

  console.log('%c========================================', 'color: #06b6d4;');
  console.log('%cResults: ' + passed + ' passed, ' + failed + ' failed',
    failed === 0 ? 'color: #10b981; font-weight: bold;' : 'color: #ef4444; font-weight: bold;');
  console.log('');

  return { passed: passed, failed: failed };
}

// ============================================================================
// STRUCTURED WIGNER-d TERMS (pre-simplification)
// ============================================================================

/**
 * Get raw Wigner-d terms BEFORE Algebrite simplification.
 * Each term has the form: weight * sin^l(beta/2) * cos^{2j-l}(beta/2)
 *
 * @param {number} jVal - angular momentum (float, can be half-integer)
 * @param {number} m1Val - projection m1
 * @param {number} m2Val - projection m2
 * @returns {{terms: Array<{weight: number, sinPow: number, cosPow: number}>}} 
 *          weight is the exact numerical coefficient.
 */
function getWignerDTerms(jVal, m1Val, m2Val) {
  var twoJ = Math.round(2 * jVal);
  if (Math.abs(2 * jVal - twoJ) > 1e-10) {
    return { error: 'j = ' + jVal + ' is not integer or half-integer' };
  }

  var jPlusM1 = Math.round(jVal + m1Val);
  var jMinusM1 = Math.round(jVal - m1Val);
  var jPlusM2 = Math.round(jVal + m2Val);
  var jMinusM2 = Math.round(jVal - m2Val);

  if (Math.abs((jVal + m1Val) - jPlusM1) > 1e-10 ||
      Math.abs((jVal + m2Val) - jPlusM2) > 1e-10) {
    return { error: 'j ± m must be integer' };
  }
  if (Math.abs(m1Val) - jVal > 1e-10 || Math.abs(m2Val) - jVal > 1e-10) {
    return { error: '|m| <= j violated' };
  }

  // Precompute factorials
  var maxN = Math.max(jPlusM1, jMinusM1, jPlusM2, jMinusM2, twoJ + 2);
  var facts = _factorials(maxN);

  var prefactor = Math.sqrt(facts[jPlusM1] * facts[jMinusM1] *
                            facts[jPlusM2] * facts[jMinusM2]);

  var terms = [];
  for (var l = 0; l <= twoJ; l++) {
    var k = (l + m2Val - m1Val) / 2;
    if (Math.abs(k - Math.round(k)) > 1e-10) continue;
    k = Math.round(k);
    if (k < Math.max(0, m2Val - m1Val) - 1e-10) continue;
    if (k > Math.min(jMinusM1, jPlusM2) + 1e-10) continue;

    var sign = ((Math.round(m1Val - m2Val) + k) % 2 === 0) ? 1 : -1;

    var d1 = jMinusM1 - k;
    var d2 = jPlusM2 - k;
    var d3 = Math.round(m1Val - m2Val) + k;
    var d4 = k;
    if (d1 < 0 || d2 < 0 || d3 < 0 || d4 < 0) continue;

    var denom = facts[d1] * facts[d2] * facts[d3] * facts[d4];
    if (denom === 0) continue;
    var weight = sign * prefactor / denom;

    terms.push({
      weight: weight,
      sinPow: l,
      cosPow: twoJ - l
    });
  }

  return { terms: terms };
}

// Factorial cache
var _factCache = [1];
function _factorials(n) {
  for (var i = _factCache.length; i <= n; i++) {
    _factCache[i] = _factCache[i - 1] * i;
  }
  return _factCache;
}

// ============================================================================
// HALF-ANGLE EXPANSION — exact symbolic formula
// ============================================================================

/**
 * Integer GCD for reducing fractions.
 */
function _intGcd(a, b) {
  a = Math.abs(a); b = Math.abs(b);
  while (b) { var t = b; b = a % b; a = t; }
  return a;
}

/**
 * Format integer fraction as Algebrite-compatible string.
 */
function _fracStr(num, den) {
  if (den < 0) { num = -num; den = -den; }
  if (num === 0) return '0';
  var g = _intGcd(Math.abs(num), den);
  num /= g; den /= g;
  if (den === 1) return num === 1 ? '1' : String(num);
  return num + '/' + den;
}

/**
 * Expand sin^a(θ/2) * cos^b(θ/2) into exact {cos(kθ/2), sin(kθ/2), 1} basis.
 *
 * Uses the closed-form Fourier expansion derived from Euler's formula:
 *
 *   sin^a(θ/2) cos^b(θ/2) = δ_{a,even}·s₀/2^N
 *     + Σ_{m} s_m/2^{N-1} · {cos|sin}(mθ/2)
 *
 * where N = a+b, s_m = Σ_{p+q=(N+m)/2} binom(a,p) binom(b,q) (-1)^{a-p},
 * and the cos/sin choice and sign follow from the factor 1/i^a:
 *
 *   a mod 4 = 0 → +cos,   a mod 4 = 1 → +sin,
 *   a mod 4 = 2 → −cos,   a mod 4 = 3 → −sin
 *
 * All coefficients are exact rational strings — no floats.
 *
 * @param {number} sinPow - exponent of sin(θ/2)
 * @param {number} cosPow - exponent of cos(θ/2)
 * @returns {Array<{s: string, func: string, k: number}>}
 */
function expandHalfAngleBasis(sinPow, cosPow) {
  var a = sinPow, b = cosPow, N = a + b;
  var facts = _factorials(N);
  var result = [];
  var aMod4 = a % 4;
  
  // Iterate over m with same parity as N
  for (var m = (N % 2); m <= N; m += 2) {
    var r = (N + m) / 2;
    // Compute s_m = sum over p+q=r
    var sm = 0;
    var pMin = Math.max(0, r - b);
    var pMax = Math.min(a, r);
    for (var p = pMin; p <= pMax; p++) {
      var q = r - p;
      sm += ((a - p) % 2 === 0 ? 1 : -1)
          * _binom(a, p, facts)
          * _binom(b, q, facts);
    }
    if (sm === 0) continue;
    
    if (m === 0) {
      // Constant term: exists only for even a (odd a gives s₀=0)
      if (a % 2 !== 0) continue;
      var constSign = (aMod4 === 0) ? 1 : -1;
      result.push({ s: _fracStr(constSign * sm, Math.pow(2, N)),
                    func: '1', k: 0 });
    } else {
      // Fourier mode: cos if a even, sin if a odd
      var func = (a % 2 === 0) ? 'cos' : 'sin';
      var fourierSign = (aMod4 === 0 || aMod4 === 1) ? 1 : -1;
      result.push({ s: _fracStr(fourierSign * sm, Math.pow(2, N - 1)),
                    func: func, k: m });
    }
  }
  
  return result;
}

function _binom(n, k, facts) {
  if (k < 0 || k > n) return 0;
  if (!facts) facts = _factorials(n);
  return Math.round(facts[n] / (facts[k] * facts[n - k]));
}

// ============================================================================
// EXACT WIGNER-d WEIGHTS (symbolic strings, no floats)
// ============================================================================

/**
 * Compute exact Wigner-d weight strings for a given (j,m1,m2).
 * Each entry: {weightStr, sinPow, cosPow}
 * Adapted for standalone use from get-angle.js.
 */
function _getExactWignerDWeights(J, m1, m2) {
  var twoJ = Math.round(2 * J);
  var jpm1 = Math.round(J + m1), jmm1 = Math.round(J - m1);
  var jpm2 = Math.round(J + m2), jmm2 = Math.round(J - m2);

  // Precompute factorials
  var maxN = Math.max(jpm1, jmm1, jpm2, jmm2, twoJ + 2, 10);
  var facts = [1];
  for (var i = 1; i <= maxN; i++) facts[i] = facts[i-1] * i;

  // Numerator product under sqrt: (j+m1)! (j-m1)! (j+m2)! (j-m2)!
  var numNum = 1;
  for (var i = 1; i <= jpm1; i++) numNum *= i;
  for (var i = 1; i <= jmm1; i++) numNum *= i;
  for (var i = 1; i <= jpm2; i++) numNum *= i;
  for (var i = 1; i <= jmm2; i++) numNum *= i;

  var weights = [];
  for (var l = 0; l <= twoJ; l++) {
    var k = (l + m2 - m1) / 2;
    if (Math.abs(k - Math.round(k)) > 1e-10) continue;
    k = Math.round(k);
    if (k < Math.max(0, m2 - m1)) continue;
    if (k > Math.min(jmm1, jpm2)) continue;

    var sign = ((Math.round(m1 - m2) + k) % 2 === 0) ? '' : '-';

    // Denominator product
    var denom = 1;
    if (jmm1 - k >= 0) { for (var i = 1; i <= jmm1 - k; i++) denom *= i; }
    if (jpm2 - k >= 0) { for (var i = 1; i <= jpm2 - k; i++) denom *= i; }
    if (Math.round(m1 - m2) + k >= 0) { for (var i = 1; i <= Math.round(m1 - m2) + k; i++) denom *= i; }
    if (k >= 0) { for (var i = 1; i <= k; i++) denom *= i; }

    // Extract perfect squares from numNum
    var p = 1, r = numNum;
    for (var i = 2; i * i <= r; i++) {
      while (r % (i * i) === 0) { r /= (i * i); p *= i; }
    }

    // Reduce p/denom
    var g = _intGcd(p, denom); p /= g; denom /= g;

    // Build weight string
    var wStr;
    if (r === 1 && p === 1 && denom === 1) wStr = '1';
    else if (r === 1 && denom === 1) wStr = String(p);
    else if (p === 1 && r === 1) wStr = '1/' + denom;
    else if (p === 1 && denom === 1) wStr = 'sqrt(' + r + ')';
    else if (denom === 1 && r === 1) wStr = String(p);
    else if (p === 1) wStr = 'sqrt(' + r + ')/' + denom;
    else if (denom === 1) wStr = p + '*sqrt(' + r + ')';
    else if (r === 1) wStr = p + '/' + denom;
    else wStr = p + '*sqrt(' + r + ')/' + denom;

    if (wStr !== '1' && wStr !== '0') wStr = sign + wStr;
    else if (wStr === '1') wStr = sign + '1';

    weights.push({ weightStr: wStr, sinPow: l, cosPow: twoJ - l });
  }
  return weights;
}

// ============================================================================
// SIMPLIFIED WIGNER-d — half-angle Fourier expansion
// ============================================================================

/**
 * Compute Wigner d-function using the half-angle Fourier expansion,
 * producing a simplified expression in the {cos(kθ/2), sin(kθ/2)} basis.
 *
 * Uses the exact formula:
 *   d^j_{m1,m2}(β) = Σ_l w_l · sin^l(β/2) cos^{2j-l}(β/2)
 * where each sin^a cos^b term is expanded via expandHalfAngleBasis().
 *
 * @param {string} jStr
 * @param {string} m1Str
 * @param {string} m2Str
 * @returns {{symbolic: string, latex: string, groups: Array, error: string}}
 */
function computeWignerDSimplified(jStr, m1Str, m2Str) {
  if (typeof Algebrite === 'undefined') {
    return { error: 'Algebrite library is not loaded.' };
  }

  var pj, pm1, pm2;
  try {
    pj = parseQuantumNumber(jStr);
    pm1 = parseQuantumNumber(m1Str);
    pm2 = parseQuantumNumber(m2Str);
  } catch (e) {
    return { error: e.message };
  }

  var jVal = pj.value, m1Val = pm1.value, m2Val = pm2.value;
  var twoJ = Math.round(2 * jVal);

  // Get exact symbolic weights
  var weights = _getExactWignerDWeights(jVal, m1Val, m2Val);
  if (weights.length === 0) {
    return { symbolic: '0', latex: '0', groups: [] };
  }

  // Build coefficient map: key "func_k" → string array of coefficient terms
  var coeffMap = {};
  var zeroStr = '0';

  for (var w = 0; w < weights.length; w++) {
    var wt = weights[w];
    var expansion = expandHalfAngleBasis(wt.sinPow, wt.cosPow);
    for (var e = 0; e < expansion.length; e++) {
      var ex = expansion[e];
      var key = ex.func + '_' + ex.k;
      if (!coeffMap[key]) {
        coeffMap[key] = { coeffStrs: [], func: ex.func, k: ex.k };
      }
      // Combine weight string with expansion coefficient string
      var combined;
      if (ex.s === '0') continue;
      if (ex.s === '1') combined = wt.weightStr;
      else if (ex.s === '-1') combined = '-(' + wt.weightStr + ')';
      else combined = '(' + wt.weightStr + ')*(' + ex.s + ')';
      coeffMap[key].coeffStrs.push(combined);
    }
  }

  // Simplify each group coefficient and build expression
  var terms = [];
  var groups = [];

  for (var key in coeffMap) {
    var grp = coeffMap[key];
    var func = grp.func;
    var k = grp.k;
    var coeffStrs = grp.coeffStrs;

    // Sum all coefficient contributions
    var sumStr;
    if (coeffStrs.length === 0) continue;
    if (coeffStrs.length === 1) sumStr = coeffStrs[0];
    else sumStr = '(' + coeffStrs.join(')+(') + ')';

    // Simplify via Algebrite
    var simplified;
    try {
      simplified = Algebrite.run('simplify(' + sumStr + ')').trim();
    } catch(e) {
      simplified = sumStr;
    }
    if (simplified === '0') continue;

    // For k=0 (constant), k=1 (half-angle), else full simplification
    var trigExpr;
    if (func === '1') {
      trigExpr = simplified;
    } else if (k === 1) {
      trigExpr = simplified + '*' + func + '(beta/2)';
    } else if (k % 2 === 0) {
      // Even k: cos(k*beta/2) = cos((k/2)*beta) — full-angle form
      var fullK = k / 2;
      trigExpr = simplified + '*' + func + '(' + (fullK === 1 ? '' : fullK + '*') + 'beta)';
    } else {
      // Odd k: stays as half-angle
      trigExpr = simplified + '*' + func + '((' + k + '/2)*beta)';
    }

    terms.push(trigExpr);
    groups.push({ coeff: simplified, func: func, k: k });
  }

  if (terms.length === 0) {
    return { symbolic: '0', latex: '0', groups: [] };
  }

  // Build symbolic expression
  var symExpr = terms.join(' + ').replace(/\+ -/g, '- ');

  // Build LaTeX via Algebrite printlatex
  var latexExpr = symExpr;
  try {
    latexExpr = Algebrite.run('simplify(' + symExpr + ')').trim();
    // Expand (only helps when there are nested products)
    latexExpr = Algebrite.run('expand(' + latexExpr + ')').trim();
    latexExpr = Algebrite.run('printlatex(' + latexExpr + ')').trim();
    // Fix LaTeX commands
    latexExpr = latexExpr.replace(/(?<!\\)cos/g, '\\cos');
    latexExpr = latexExpr.replace(/(?<!\\)sin/g, '\\sin');
    latexExpr = latexExpr.replace(/(?<!\\)beta/g, '\\beta');
  } catch(e) {
    latexExpr = symExpr;
  }

  return {
    symbolic: symExpr,
    latex: latexExpr,
    groups: groups
  };
}


// Export for use in other modules
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    computeWignerD: computeWignerD,
    computeWignerDMatrix: computeWignerDMatrix,
    runWignerDSanityChecks: runWignerDSanityChecks,
    getWignerDTerms: getWignerDTerms,
    expandHalfAngleBasis: expandHalfAngleBasis,
    computeWignerDSimplified: computeWignerDSimplified
  };
}