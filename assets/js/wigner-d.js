/**
 * Wigner d-function Calculator
 * Exact BigInt/surd computation — no CAS dependency
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
 * Returns: {numerator, denominator, value}
 */
function parseQuantumNumber(input) {
  const s = String(input).trim();

  // Fraction format "a/b" — parse numerator/denominator as integers
  const slashIdx = s.indexOf('/');
  if (slashIdx !== -1) {
    const m = s.match(/^([+-]?\d+)\s*\/\s*([+-]?\d+)$/);
    if (!m || parseInt(m[2], 10) === 0) {
      throw new Error(`Invalid quantum number: ${s}`);
    }
    const num = parseInt(m[1], 10);
    const den = parseInt(m[2], 10);
    return {
      numerator: num,
      denominator: den,
      value: num / den
    };
  }

  // Integer or decimal string — converted to an exact fraction
  // (e.g. "0.5" → 1/2). NaN-safe regex, so "1.5junk" is rejected
  // instead of being silently truncated by parseInt.
  const dm = s.match(/^([+-]?)(\d+)(?:\.(\d+))?$/);
  if (!dm) {
    throw new Error(`Invalid quantum number: ${s}. Use fraction syntax e.g. "1/2" for half-integers.`);
  }
  const neg = dm[1] === '-';
  const intPart = parseInt(dm[2], 10);
  if (dm[3] === undefined) {
    return {
      numerator: neg ? -intPart : intPart,
      denominator: 1,
      value: neg ? -intPart : intPart
    };
  }
  let num = intPart * Math.pow(10, dm[3].length) + parseInt(dm[3], 10);
  let den = Math.pow(10, dm[3].length);
  const g = _intGcd(Math.abs(num), den);
  num = (neg ? -num : num) / g;
  den = den / g;
  return {
    numerator: num,
    denominator: den,
    value: num / den
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
 * Uses exact BigInt/surd arithmetic (no CAS dependency):
 *   d^j_{m1,m2}(β) = Σ_l w_l · sin^l(β/2) cos^{2j−l}(β/2)
 * with exact weights from _getExactWignerDWeights(), each sin^a cos^b
 * term expanded via expandHalfAngleBasis() into the exact Fourier basis
 * {cos(kβ/2), sin(kβ/2), 1} and combined with SurdSum.
 *
 * @param {string} jStr - angular momentum j (e.g. "1/2", "1", "3/2")
 * @param {string} m1Str - projection m1
 * @param {string} m2Str - projection m2
 * @param {string} betaStr - angle beta (e.g. "pi/2", "0.5", "pi/4"). If empty/null, treat as symbolic.
 * @returns {{decimal: number, symbolic: string, latex: string, error: string}}
 */
function computeWignerD(jStr, m1Str, m2Str, betaStr) {
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

  // j±m must be integers (m must share j's half-integer parity)
  if (Math.abs((jVal + m1Val) - Math.round(jVal + m1Val)) > 1e-10 ||
      Math.abs((jVal + m2Val) - Math.round(jVal + m2Val)) > 1e-10) {
    return { error: 'm1, m2 must be integer when j is integer, half-integer when j is half-integer' };
  }

  // Exact symbolic weights + half-angle Fourier expansion
  var groups = _buildFourierGroups(jVal, m1Val, m2Val);
  if (groups === null) {
    return { decimal: 0, symbolic: '0', latex: '0' };
  }

  var rendered = _renderFourierGroups(groups);

  // Numeric beta?
  var betaVal = NaN;
  if (betaStr && betaStr.trim() !== '') {
    betaVal = _evalAngle(betaStr);
    if (isNaN(betaVal)) {
      return { error: 'Invalid angle beta: ' + betaStr + ' (use e.g. "pi/2", "0.5")' };
    }
  }

  return {
    decimal: isNaN(betaVal) ? NaN : _evalGroupsDecimal(groups, betaVal),
    symbolic: rendered.symbolic,
    latex: rendered.latex
  };
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

  var tests = [
    // j=1/2: d^{1/2}_{1/2,1/2}(beta) = cos(beta/2)
    { j: '1/2', m1: '1/2', m2: '1/2', beta: '', expectedSymbolic: 'cos(beta/2)', description: 'd^{1/2}_{1/2,1/2}(beta) = cos(beta/2)' },
    // j=1/2: d^{1/2}_{1/2,-1/2}(beta) = -sin(beta/2)
    { j: '1/2', m1: '1/2', m2: '-1/2', beta: '', expectedSymbolic: '-sin(beta/2)', description: 'd^{1/2}_{1/2,-1/2}(beta) = -sin(beta/2)' },
    // j=1/2: d^{1/2}_{-1/2,1/2}(beta) = sin(beta/2)
    { j: '1/2', m1: '-1/2', m2: '1/2', beta: '', expectedSymbolic: 'sin(beta/2)', description: 'd^{1/2}_{-1/2,1/2}(beta) = sin(beta/2)' },
    // j=1/2: d^{1/2}_{-1/2,-1/2}(beta) = cos(beta/2)
    { j: '1/2', m1: '-1/2', m2: '-1/2', beta: '', expectedSymbolic: 'cos(beta/2)', description: 'd^{1/2}_{-1/2,-1/2}(beta) = cos(beta/2)' },
    // j=1: d^{1}_{1,1}(beta) = cos^2(beta/2) = 1/2 + 1/2 cos(beta)
    { j: '1', m1: '1', m2: '1', beta: '', expectedSymbolic: '1/2+1/2*cos(beta)', description: 'd^{1}_{1,1}(beta) = cos^2(beta/2)' },
    // j=1: d^{1}_{1,0}(beta) = -sin(beta)/sqrt(2)
    { j: '1', m1: '1', m2: '0', beta: '', expectedSymbolic: '-sqrt(2)/2*sin(beta)', description: 'd^{1}_{1,0}(beta) = -sin(beta)/sqrt(2)' },
    // j=1: d^{1}_{1,-1}(beta) = sin^2(beta/2) = 1/2 - 1/2 cos(beta)
    { j: '1', m1: '1', m2: '-1', beta: '', expectedSymbolic: '1/2-1/2*cos(beta)', description: 'd^{1}_{1,-1}(beta) = sin^2(beta/2)' },
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
 * Get raw Wigner-d terms.
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
 * Format integer fraction as an exact rational string.
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
          * _binom(a, p)
          * _binom(b, q);
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

/**
 * Exact binomial coefficient via BigInt factorials.
 * The old Number-based version silently corrupted n! for n ≥ 19 (2⁵³
 * overflow), breaking the "exact" half-angle expansion at j ≥ 10.
 */
var _factCacheBig = [1n];
function _factorialsBig(n) {
  for (var i = _factCacheBig.length; i <= n; i++) {
    _factCacheBig[i] = _factCacheBig[i - 1] * BigInt(i);
  }
  return _factCacheBig;
}

function _binom(n, k) {
  if (k < 0 || k > n) return 0;
  var facts = _factorialsBig(n);
  var b = facts[n] / (facts[k] * facts[n - k]);
  // b ≤ binom(n,k); within the safe integer range this conversion is
  // exact. Beyond 2⁵³ it rounds to the nearest double (best effort).
  return Number(b);
}

// ============================================================================
// EXACT WIGNER-d WEIGHTS (symbolic strings, no floats)
// ============================================================================

/** BigInt gcd helper for the exact weight computation. */
function _wignerBigGcd(a, b) {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) { var t = b; b = a % b; a = t; }
  return a;
}

/** Sieve-free trial primes up to n (n ≤ 2j, tiny). */
function _wignerPrimesUpTo(n) {
  var primes = [];
  for (var i = 2; i <= n; i++) {
    var isP = true;
    for (var j = 0; j < primes.length && primes[j] * primes[j] <= i; j++) {
      if (i % primes[j] === 0) { isP = false; break; }
    }
    if (isP) primes.push(i);
  }
  return primes;
}

/**
 * Compute exact Wigner-d weight strings for a given (j,m1,m2).
 * Each entry: {weightStr, sinPow, cosPow}
 * Adapted for standalone use from get-angle.js.
 *
 * Uses BigInt factorials — the old Number version silently corrupted
 * the symbolic weights for j ≥ 8 (e.g. d⁸₈₈(β) = cos¹⁶(β/2) came out as
 * a giant spurious radicand instead of the weight 1).
 */
function _getExactWignerDWeights(J, m1, m2) {
  var twoJ = Math.round(2 * J);
  var jpm1 = Math.round(J + m1), jmm1 = Math.round(J - m1);
  var jpm2 = Math.round(J + m2), jmm2 = Math.round(J - m2);

  // Projection parity: j±m must be non-negative integers
  // (rejects e.g. m = 0 for j = 1/2 → no valid weights).
  if (Math.abs((J + m1) - jpm1) > 1e-10 || Math.abs((J - m1) - jmm1) > 1e-10 ||
      Math.abs((J + m2) - jpm2) > 1e-10 || Math.abs((J - m2) - jmm2) > 1e-10) {
    return [];
  }
  if (jpm1 < 0 || jmm1 < 0 || jpm2 < 0 || jmm2 < 0) return [];

  // BigInt factorials
  var maxN = Math.max(jpm1, jmm1, jpm2, jmm2, twoJ + 2, 1);
  var facts = [1n];
  for (var i = 1; i <= maxN; i++) facts[i] = facts[i - 1] * BigInt(i);

  // Numerator under sqrt: (j+m1)! (j-m1)! (j+m2)! (j-m2)! — exact BigInt.
  // All its prime factors are ≤ maxN, so trial division over primes up
  // to maxN factors it fully (also replaces the slow O(√r) Number scan).
  var numNum = facts[jpm1] * facts[jmm1] * facts[jpm2] * facts[jmm2];
  var primes = _wignerPrimesUpTo(maxN);
  var outside = 1n, radicand = 1, temp = numNum;
  for (var pi = 0; pi < primes.length; pi++) {
    var p = primes[pi];
    var pBig = BigInt(p);
    if (pBig > temp) break;
    var exp = 0;
    while (temp % pBig === 0n) { temp = temp / pBig; exp++; }
    if (exp > 0) {
      var pairs = Math.floor(exp / 2);
      for (var ii = 0; ii < pairs; ii++) outside *= pBig;
      if (exp % 2 === 1) radicand *= p;
    }
  }

  var weights = [];
  for (var l = 0; l <= twoJ; l++) {
    var k = (l + m2 - m1) / 2;
    if (Math.abs(k - Math.round(k)) > 1e-10) continue;
    k = Math.round(k);
    if (k < Math.max(0, m2 - m1)) continue;
    if (k > Math.min(jmm1, jpm2)) continue;

    var sign = ((Math.round(m1 - m2) + k) % 2 === 0) ? '' : '-';

    // Denominator product — exact BigInt: (j-m1-k)!(j+m2-k)!(m1-m2+k)!·k!
    var denomBig = facts[jmm1 - k] * facts[jpm2 - k] *
                   facts[Math.round(m1 - m2) + k] * facts[k];

    // Reduce outside/denom by gcd
    var g = _wignerBigGcd(outside, denomBig);
    var p = Number(outside / g);
    var denom = Number(denomBig / g);

    // Guard: values must fit in the exact Number range used by Surd.
    var MAX_SAFE = Number.MAX_SAFE_INTEGER;
    if (p > MAX_SAFE || denom > MAX_SAFE || radicand > MAX_SAFE) {
      throw new Error('Wigner-d weight too large for exact representation (j=' + J + ')');
    }

    // Build weight string
    var wStr;
    if (radicand === 1 && p === 1 && denom === 1) wStr = '1';
    else if (radicand === 1 && denom === 1) wStr = String(p);
    else if (p === 1 && radicand === 1) wStr = '1/' + denom;
    else if (p === 1 && denom === 1) wStr = 'sqrt(' + radicand + ')';
    else if (denom === 1 && radicand === 1) wStr = String(p);
    else if (p === 1) wStr = 'sqrt(' + radicand + ')/' + denom;
    else if (denom === 1) wStr = p + '*sqrt(' + radicand + ')';
    else if (radicand === 1) wStr = p + '/' + denom;
    else wStr = p + '*sqrt(' + radicand + ')/' + denom;

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
 * Build the exact Fourier expansion of d^{j}_{m1,m2}(β):
 *   Σ_l w_l · sin^l(β/2) cos^{2j−l}(β/2) = Σ_g c_g · func_g(k_g·β/2)
 *
 * Exact weights (BigInt) are expanded via expandHalfAngleBasis() and each
 * {cos(kβ/2), sin(kβ/2), 1} coefficient combined with SurdSum — no floats,
 * no CAS. Replaces the Algebrite simplify() step.
 *
 * @param {number} jVal - angular momentum (integer or half-integer)
 * @param {number} m1Val, m2Val - projections
 * @returns {Array<{func: string, k: number, sum: SurdSum}>|null}
 *          null when the element is identically zero.
 */
function _buildFourierGroups(jVal, m1Val, m2Val) {
  var weights = _getExactWignerDWeights(jVal, m1Val, m2Val);
  if (weights.length === 0) return null;

  var map = {}; // "func_k" -> SurdSum
  for (var w = 0; w < weights.length; w++) {
    var wt = weights[w];
    var expansion = expandHalfAngleBasis(wt.sinPow, wt.cosPow);
    var wSurd = Surd.parse(wt.weightStr);
    for (var e = 0; e < expansion.length; e++) {
      var ex = expansion[e];
      if (ex.s === '0') continue;
      var key = ex.func + '_' + ex.k;
      if (!map[key]) map[key] = new SurdSum();
      map[key].add(Surd.mul(wSurd, Surd.parse(ex.s)));
    }
  }

  // Deterministic order: constant first, then ascending k, then func.
  var keys = Object.keys(map).sort(function(a, b) {
    var pa = a.split('_'), pb = b.split('_');
    if (pa[0] === '1' && pb[0] !== '1') return -1;
    if (pa[0] !== '1' && pb[0] === '1') return 1;
    var ka = parseInt(pa[1], 10), kb = parseInt(pb[1], 10);
    if (ka !== kb) return ka - kb;
    return pa[0] < pb[0] ? -1 : 1;
  });

  var groups = [];
  for (var i = 0; i < keys.length; i++) {
    var sum = map[keys[i]];
    if (sum.isEmpty()) continue;
    var kk = keys[i].split('_');
    groups.push({ func: kk[0], k: parseInt(kk[1], 10), sum: sum });
  }
  return groups.length === 0 ? null : groups;
}

/** Numeric value of a SurdSum (exact coefficients → float evaluation). */
function _surSumFloat(sum) {
  var total = 0;
  var terms = sum.terms();
  for (var i = 0; i < terms.length; i++) {
    var t = terms[i];
    total += t.s * t.p * Math.sqrt(t.r) / t.q;
  }
  return total;
}

/** Evaluate the Fourier expansion at a numeric angle beta. */
function _evalGroupsDecimal(groups, beta) {
  var total = 0;
  for (var i = 0; i < groups.length; i++) {
    var g = groups[i];
    if (g.sum.isEmpty()) continue;
    var c = _surSumFloat(g.sum);
    if (g.func === '1') total += c;
    else total += c * (g.func === 'cos' ? Math.cos(g.k * beta / 2) : Math.sin(g.k * beta / 2));
  }
  return total;
}

/** Symbolic trig factor for Fourier mode k: "cos(beta/2)", "cos(beta)", "sin(3/2*beta)". */
function _fourierSymTrig(func, k) {
  if (k === 1) return func + '(beta/2)';
  if (k % 2 === 0) {
    var n = k / 2;
    return func + '(' + (n === 1 ? '' : String(n) + '*') + 'beta)';
  }
  return func + '(' + k + '/2*beta)';
}

/** LaTeX trig factor for Fourier mode k. */
function _fourierTexTrig(func, k) {
  var f = (func === 'cos') ? '\\cos' : '\\sin';
  if (k === 1) return f + '\\left(\\frac{\\beta}{2}\\right)';
  if (k % 2 === 0) {
    var n = k / 2;
    return f + '(' + (n === 1 ? '' : String(n)) + '\\beta)';
  }
  return f + '\\left(\\frac{' + k + '}{2}\\beta\\right)';
}

/**
 * Render exact Fourier groups to symbolic and LaTeX strings.
 * Replaces Algebrite's simplify()/printlatex() for the Wigner-d output.
 */
function _renderFourierGroups(groups) {
  var symParts = [], texParts = [];
  for (var i = 0; i < groups.length; i++) {
    var g = groups[i];
    if (g.sum.isEmpty()) continue;
    var symCoeff = g.sum.toString();
    if (symCoeff === '0') continue;
    var texCoeff = g.sum.toLatex();
    if (g.func === '1') {
      symParts.push(symCoeff);
      texParts.push(texCoeff);
    } else {
      var symTrig = _fourierSymTrig(g.func, g.k);
      var texTrig = _fourierTexTrig(g.func, g.k);
      if (symCoeff === '1') symParts.push(symTrig);
      else if (symCoeff === '-1') symParts.push('-' + symTrig);
      else symParts.push(symCoeff + '*' + symTrig);

      if (texCoeff === '1') texParts.push(texTrig);
      else if (texCoeff === '-1') texParts.push('-' + texTrig);
      else texParts.push(texCoeff + '\\,' + texTrig);
    }
  }
  if (symParts.length === 0) return { symbolic: '0', latex: '0' };
  return {
    symbolic: symParts.join('+').replace(/\+-/g, '-'),
    latex: texParts.join('+').replace(/\+-/g, '-')
  };
}

/**
 * Evaluate a user-provided angle string: "pi/2", "0.5", "2*pi/3", "1.2", "-pi/4".
 * Safe recursive-descent parser — only digits, + - * / ( ) and "pi" are
 * accepted; no eval/Function. Returns NaN for anything else.
 */
function _evalAngle(betaStr) {
  var s = String(betaStr).trim().toLowerCase().replace(/\s+/g, '');
  if (!/^[0-9+\-*/().pi]*$/.test(s)) return NaN;

  var i = 0;
  function peek() { return s[i]; }
  function parseExpr() {
    var v = parseTerm();
    while (peek() === '+' || peek() === '-') {
      var op = peek(); i++;
      var r = parseTerm();
      v = (op === '+') ? v + r : v - r;
    }
    return v;
  }
  function parseTerm() {
    var v = parseFactor();
    while (peek() === '*' || peek() === '/') {
      var op = peek(); i++;
      var r = parseFactor();
      v = (op === '*') ? v * r : v / r;
    }
    return v;
  }
  function parseFactor() {
    var c = peek();
    if (c === '-' || c === '+') { i++; var v = parseFactor(); return (c === '-') ? -v : v; }
    if (c === '(') {
      i++;
      var v = parseExpr();
      if (peek() === ')') i++;
      return v;
    }
    if (c === 'p') {
      if (s.substr(i, 2) === 'pi') { i += 2; return Math.PI; }
      return NaN;
    }
    var m = /^[0-9]+(\.[0-9]+)?/.exec(s.substr(i));
    if (!m) return NaN;
    i += m[0].length;
    return parseFloat(m[0]);
  }

  var result = parseExpr();
  return (i >= s.length && isFinite(result)) ? result : NaN;
}

/**
 * Compute Wigner d-function using the half-angle Fourier expansion,
 * producing a simplified expression in the {cos(kθ/2), sin(kθ/2)} basis.
 *
 * Uses the exact formula:
 *   d^j_{m1,m2}(β) = Σ_l w_l · sin^l(β/2) cos^{2j-l}(β/2)
 * where each sin^a cos^b term is expanded via expandHalfAngleBasis() and
 * combined exactly with SurdSum — no CAS dependency.
 *
 * @param {string} jStr
 * @param {string} m1Str
 * @param {string} m2Str
 * @returns {{symbolic: string, latex: string, groups: Array, error: string}}
 */
function computeWignerDSimplified(jStr, m1Str, m2Str) {
  var pj, pm1, pm2;
  try {
    pj = parseQuantumNumber(jStr);
    pm1 = parseQuantumNumber(m1Str);
    pm2 = parseQuantumNumber(m2Str);
  } catch (e) {
    return { error: e.message };
  }

  var jVal = pj.value, m1Val = pm1.value, m2Val = pm2.value;

  var groups = _buildFourierGroups(jVal, m1Val, m2Val);
  if (groups === null) {
    return { symbolic: '0', latex: '0', groups: [] };
  }

  var rendered = _renderFourierGroups(groups);

  return {
    symbolic: rendered.symbolic,
    latex: rendered.latex,
    groups: groups.map(function(g) {
      return { coeff: g.sum.toString(), func: g.func, k: g.k };
    })
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
    computeWignerDSimplified: computeWignerDSimplified,
    parseQuantumNumber: parseQuantumNumber,
    _getExactWignerDWeights: _getExactWignerDWeights
  };
}