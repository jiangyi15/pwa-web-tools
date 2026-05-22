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

// Export for use in other modules
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    computeWignerD: computeWignerD,
    computeWignerDMatrix: computeWignerDMatrix,
    runWignerDSanityChecks: runWignerDSanityChecks
  };
}