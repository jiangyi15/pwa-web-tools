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
  const trimmed = String(input).trim();

  // Handle fraction format "a/b"
  if (trimmed.includes('/')) {
    const parts = trimmed.split('/');
    if (parts.length !== 2) {
      throw new Error(`Invalid fraction format: ${trimmed}`);
    }
    const num = parseFloat(parts[0]);
    const den = parseFloat(parts[1]);
    if (isNaN(num) || isNaN(den) || den === 0) {
      throw new Error(`Invalid fraction: ${trimmed}`);
    }
    return {
      numerator: num,
      denominator: den,
      value: num / den,
      algebriteForm: `(${num}/${den})`
    };
  }

  // Handle decimal or integer
  const val = parseFloat(trimmed);
  if (isNaN(val)) {
    throw new Error(`Invalid number: ${trimmed}`);
  }

  const twiceVal = val * 2;
  const isHalfInteger = Math.abs(twiceVal - Math.round(twiceVal)) < 1e-10;

  if (Number.isInteger(val)) {
    return {
      numerator: val,
      denominator: 1,
      value: val,
      algebriteForm: val < 0 ? `(${val})` : String(val)
    };
  }

  if (isHalfInteger) {
    const num = Math.round(twiceVal);
    return {
      numerator: num,
      denominator: 2,
      value: val,
      algebriteForm: `(${num}/2)`
    };
  }

  if (typeof Algebrite !== 'undefined') {
    const rational = Algebrite.run(`rationalize(${val})`);
    return {
      numerator: null,
      denominator: null,
      value: val,
      algebriteForm: rational.trim()
    };
  }

  throw new Error(`Cannot convert ${val} to a rational number`);
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

  // Generate m values: j, j-1, ..., -j
  var mValues = [];
  for (var i = 0; i <= twoJ; i++) {
    mValues.push(jVal - i);
  }

  var matrix = [];
  for (var row = 0; row < mValues.length; row++) {
    var rowData = [];
    for (var col = 0; col < mValues.length; col++) {
      var m1 = mValues[row];
      var m2 = mValues[col];
      var result = computeWignerD(jStr, String(m1), String(m2), betaStr);
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
// HALF-ANGLE EXPANSION TABLES
// ============================================================================

/**
 * Pre-computed expansion of sin^b(θ/2) * cos^c(θ/2) into
 * basis of {1, cos(kθ/2), sin(kθ/2)} for k=1..(b+c).
 * Each entry: { coeff, func: "cos"|"sin"|"1", k }
 * k=0 means constant (func="1"), k>0 is the multiplier of θ/2.
 */
var _halfAngleTable = {
  '0,0': [ {s:'1', f:'1', k:0} ],
  '0,1': [ {s:'1', f:'cos', k:1} ],
  '1,0': [ {s:'1', f:'sin', k:1} ],
  '0,2': [ {s:'1/2', f:'1', k:0},  {s:'1/2', f:'cos', k:2} ],
  '1,1': [ {s:'1/2', f:'sin', k:2} ],
  '2,0': [ {s:'1/2', f:'1', k:0},  {s:'-1/2', f:'cos', k:2} ],
  '0,3': [ {s:'3/4', f:'cos', k:1}, {s:'1/4', f:'cos', k:3} ],
  '1,2': [ {s:'1/4', f:'sin', k:1}, {s:'1/4', f:'sin', k:3} ],
  '2,1': [ {s:'1/4', f:'sin', k:1}, {s:'-1/4', f:'sin', k:3} ],
  '3,0': [ {s:'3/4', f:'sin', k:1}, {s:'-1/4', f:'sin', k:3} ],
  '0,4': [ {s:'3/8', f:'1', k:0}, {s:'1/2',  f:'cos', k:2}, {s:'1/8', f:'cos', k:4} ],
  '1,3': [ {s:'1/4', f:'sin', k:2}, {s:'1/8', f:'sin', k:4} ],
  '2,2': [ {s:'1/8', f:'1', k:0}, {s:'-1/8', f:'cos', k:4} ],
  '3,1': [ {s:'1/4', f:'sin', k:2}, {s:'-1/8', f:'sin', k:4} ],
  '4,0': [ {s:'3/8', f:'1', k:0}, {s:'-1/2', f:'cos', k:2}, {s:'1/8', f:'cos', k:4} ],
};

/**
 * Expand sin^b(θ/2) * cos^c(θ/2) into basis functions.
 * Returns array of {coeff: number, func: "cos"|"sin"|"1", k: number}
 * Uses precomputed table for small powers (up to b+c=6), 
 * falls back to binomial computation for higher powers.
 */
function expandHalfAngleBasis(sinPow, cosPow) {
  var key = sinPow + ',' + cosPow;
  if (_halfAngleTable[key]) {
    return _halfAngleTable[key].map(function(e) {
      return { s: e.s, func: e.f, k: e.k };
    });
  }
  return _computeHalfAngleBasis(sinPow, cosPow);
}

function _computeHalfAngleBasis(sinPow, cosPow) {
  // Binomial expansion fallback for high powers
  var b = sinPow, c = cosPow, N = b + c;
  var denom = Math.pow(2, N);
  
  // Compute e^{imx} coefficients
  var eCoeffs = {};
  var binomFacts = _factorials(N);
  for (var j = 0; j <= b; j++) {
    for (var vk = 0; vk <= c; vk++) {
      var m = N - 2*j - 2*vk;
      var sign = (j % 2 === 0) ? 1 : -1;
      var coeff = sign * _binom(b, j, binomFacts) * _binom(c, vk, binomFacts);
      eCoeffs['m' + m] = (eCoeffs['m' + m] || 0) + coeff;
    }
  }
  
  // Determine basis type and extra sign from i^b
  var extraSign = 1;
  if (b % 4 === 2 || b % 4 === 3) extraSign = -1;
  
  var isSine = (b % 2 === 1);
  
  var result = [];
  for (var m = 0; m <= N; m += 2) {
    var coeff;
    if (m === 0) {
      coeff = (eCoeffs['m0'] || 0) / denom * extraSign;
    } else {
      var cp = eCoeffs['m' + m] || 0;
      var cm = eCoeffs['m' + (-m)] || 0;
      if (isSine) {
        coeff = (cp - cm) / denom * extraSign * (b % 4 === 1 ? 1 : -1);
      } else {
        coeff = (cp + cm) / denom * extraSign;
      }
    }
    
    if (Math.abs(coeff) < 1e-14) continue;
    
    result.push({
      coeff: coeff,
      func: m === 0 ? '1' : (isSine ? 'sin' : 'cos'),
      k: m,
      halfAngle: true
    });
  }
  
  return result;
}

function _binom(n, k, facts) {
  if (k < 0 || k > n) return 0;
  if (!facts) facts = _factorials(n);
  return Math.round(facts[n] / (facts[k] * facts[n - k]));
}


// Export for use in other modules
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    computeWignerD: computeWignerD,
    computeWignerDMatrix: computeWignerDMatrix,
    runWignerDSanityChecks: runWignerDSanityChecks,
    getWignerDTerms: getWignerDTerms,
    expandHalfAngleBasis: expandHalfAngleBasis
  };
}