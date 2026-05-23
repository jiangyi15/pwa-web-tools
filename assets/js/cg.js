/**
 * Clebsch-Gordan Coefficient Calculator
 * Uses Algebrite for exact symbolic computation
 * 
 * Formula:
 * ⟨j1 m1 j2 m2 | J M⟩ = δ(m1+m2, M) · √( (2J+1) · Δ(j1,j2,J) ) ·
 *   √( (j1+m1)!(j1-m1)!(j2+m2)!(j2-m2)!(J+M)!(J-M)! ) ·
 *   Σ_k (-1)^k / [ k! (j1+j2-J-k)! (j1-m1-k)! (j2+m2-k)! (J-j2+m1+k)! (J-j1-m2+k)! ]
 * 
 * where Δ(a,b,c) = (a+b-c)!(a-b+c)!(-a+b+c)! / (a+b+c+1)!
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
 * For j (angular momentum magnitude): must be non-negative; pass allowNegative=false.
 * For m (projection): may be negative; pass allowNegative=true.
 */
function isValidQuantumNumber(parsed, allowNegative = false) {
  const twiceVal = parsed.value * 2;
  // 2*value must be an integer (so value is integer or half-integer)
  if (Math.abs(twiceVal - Math.round(twiceVal)) > 1e-10) return false;
  if (!allowNegative && parsed.value < 0) return false;
  return true;
}

// ============================================================================
// SELECTION RULES
// ============================================================================

/**
 * Check selection rules for CG coefficients
 * Returns {valid: boolean, message: string, parsed: object}
 */
function checkSelectionRules(j1, m1, j2, m2, J, M) {
  // Check if Algebrite is loaded
  if (typeof Algebrite === 'undefined') {
    return { valid: false, message: 'Algebrite library is not loaded. Cannot perform calculation.' };
  }
  
  // Parse all inputs
  let pj1, pm1, pj2, pm2, pJ, pM;
  
  try {
    pj1 = parseQuantumNumber(j1);
    pm1 = parseQuantumNumber(m1);
    pj2 = parseQuantumNumber(j2);
    pm2 = parseQuantumNumber(m2);
    pJ = parseQuantumNumber(J);
    pM = parseQuantumNumber(M);
  } catch (e) {
    return { valid: false, message: e.message };
  }
  
  // Check that all j and J values are valid quantum numbers
  if (!isValidQuantumNumber(pj1)) {
    return { valid: false, message: `j₁ = ${j1} is not a valid quantum number (must be integer or half-integer ≥ 0)` };
  }
  if (!isValidQuantumNumber(pj2)) {
    return { valid: false, message: `j₂ = ${j2} is not a valid quantum number (must be integer or half-integer ≥ 0)` };
  }
  if (!isValidQuantumNumber(pJ)) {
    return { valid: false, message: `J = ${J} is not a valid quantum number (must be integer or half-integer ≥ 0)` };
  }
  
  // Check m values are valid for their respective j values (m may be negative)
  if (!isValidQuantumNumber(pm1, true)) {
    return { valid: false, message: `m₁ = ${m1} is not a valid quantum number (must be integer or half-integer)` };
  }
  if (!isValidQuantumNumber(pm2, true)) {
    return { valid: false, message: `m₂ = ${m2} is not a valid quantum number (must be integer or half-integer)` };
  }
  if (!isValidQuantumNumber(pM, true)) {
    return { valid: false, message: `M = ${M} is not a valid quantum number (must be integer or half-integer)` };
  }
  
  const j1Val = pj1.value;
  const m1Val = pm1.value;
  const j2Val = pj2.value;
  const m2Val = pm2.value;
  const JVal = pJ.value;
  const MVal = pM.value;
  
  // Check that j1+m1, j1-m1, etc. are non-negative integers
  const checkConsistency = (j, m, jName, mName) => {
    const sum = j + m;
    const diff = j - m;
    const twiceSum = sum * 2;
    const twiceDiff = diff * 2;
    if (twiceSum < 0 || Math.abs(twiceSum - Math.round(twiceSum)) > 1e-10) {
      return { valid: false, message: `${jName} + ${mName} = ${sum.toFixed(3)} is not a non-negative integer or half-integer` };
    }
    if (twiceDiff < 0 || Math.abs(twiceDiff - Math.round(twiceDiff)) > 1e-10) {
      return { valid: false, message: `${jName} − ${mName} = ${diff.toFixed(3)} is not a non-negative integer or half-integer` };
    }
    return { valid: true };
  };
  
  let check = checkConsistency(j1Val, m1Val, 'j₁', 'm₁');
  if (!check.valid) return check;
  
  check = checkConsistency(j2Val, m2Val, 'j₂', 'm₂');
  if (!check.valid) return check;
  
  check = checkConsistency(JVal, MVal, 'J', 'M');
  if (!check.valid) return check;
  
  // Check |m| ≤ j
  if (Math.abs(m1Val) - j1Val > 1e-10) {
    return { valid: false, message: `Selection rule violated: |m₁| ≤ j₁ (|${m1Val}| > ${j1Val})` };
  }
  if (Math.abs(m2Val) - j2Val > 1e-10) {
    return { valid: false, message: `Selection rule violated: |m₂| ≤ j₂ (|${m2Val}| > ${j2Val})` };
  }
  if (Math.abs(MVal) - JVal > 1e-10) {
    return { valid: false, message: `Selection rule violated: |M| ≤ J (|${MVal}| > ${JVal})` };
  }
  
  // Check m1 + m2 = M
  if (Math.abs(m1Val + m2Val - MVal) > 1e-10) {
    return { valid: false, message: `Selection rule violated: m₁ + m₂ ≠ M (${m1Val} + ${m2Val} = ${m1Val + m2Val} ≠ ${MVal})` };
  }
  
  // Check triangle inequality: |j1 - j2| ≤ J ≤ j1 + j2
  if (JVal < Math.abs(j1Val - j2Val) - 1e-10 || JVal > j1Val + j2Val + 1e-10) {
    return { valid: false, message: `Selection rule violated: |j₁ − j₂| ≤ J ≤ j₁ + j₂ (${Math.abs(j1Val - j2Val)} ≤ ${JVal} ≤ ${j1Val + j2Val})` };
  }
  
  return { valid: true, parsed: { pj1, pm1, pj2, pm2, pJ, pM } };
}

// ============================================================================
// CLEBSCH-GORDAN CALCULATION WITH ALGEBRITE
// ============================================================================

/**
 * Compute Clebsch-Gordan coefficient using Algebrite
 * Returns { decimal, symbolic, error }
 */
function computeCG(j1, m1, j2, m2, J, M) {
  // Check selection rules
  const rules = checkSelectionRules(j1, m1, j2, m2, J, M);
  if (!rules.valid) {
    return { error: rules.message };
  }
  
  const { pj1, pm1, pj2, pm2, pJ, pM } = rules.parsed;
  
  const j1Val = pj1.value;
  const m1Val = pm1.value;
  const j2Val = pj2.value;
  const m2Val = pm2.value;
  const JVal = pJ.value;
  const MVal = pM.value;
  
  // Convert to integers for the sum range calculation
  // These quantities are always integers (j±m is always integer)
  const j1PlusM1 = Math.round(j1Val + m1Val);
  const j1MinusM1 = Math.round(j1Val - m1Val);
  const j2PlusM2 = Math.round(j2Val + m2Val);
  const j2MinusM2 = Math.round(j2Val - m2Val);
  const JPlusM = Math.round(JVal + MVal);
  const JMinusM = Math.round(JVal - MVal);
  
  const j1PlusJ2MinusJ = Math.round(j1Val + j2Val - JVal);
  const JminusJ2PlusM1 = Math.round(JVal - j2Val + m1Val);
  const JminusJ1MinusM2 = Math.round(JVal - j1Val - m2Val);
  
  const kMin = Math.max(0, -JminusJ2PlusM1, -JminusJ1MinusM2);
  const kMax = Math.min(j1PlusJ2MinusJ, j1MinusM1, j2PlusM2);
  
  if (kMin > kMax) {
    // Sum is empty, coefficient is zero
    return {
      decimal: 0,
      symbolic: '0'
    };
  }
  
  try {
    // Build the Algebrite expression
    // We'll construct it piece by piece
    
    // Get Algebrite forms for the quantum numbers
    const aj1 = pj1.algebriteForm;
    const am1 = pm1.algebriteForm;
    const aj2 = pj2.algebriteForm;
    const am2 = pm2.algebriteForm;
    const aJ = pJ.algebriteForm;
    const aM = pM.algebriteForm;
    
    // Triangle coefficient Δ(j1,j2,J)
    // Δ = (j1+j2-J)! * (j1-j2+J)! * (-j1+j2+J)! / (j1+j2+J+1)!
    const triangleExpr = `((${aj1}+${aj2}-${aJ})! * (${aj1}-${aj2}+${aJ})! * (-${aj1}+${aj2}+${aJ})!) / (${aj1}+${aj2}+${aJ}+1)!`;
    
    // Prefactor: sqrt((2J+1) * Δ(j1,j2,J) * factorial_product)
    const factorialProduct = `(${aj1}+${am1})! * (${aj1}-${am1})! * (${aj2}+${am2})! * (${aj2}-${am2})! * (${aJ}+${aM})! * (${aJ}-${aM})!`;
    
    // Build the sum over k
    let sumTerms = [];
    for (let k = kMin; k <= kMax; k++) {
      const sign = (k % 2 === 0) ? '' : '-';
      const term = `${sign}1 / (${k}! * (${j1PlusJ2MinusJ - k})! * (${j1MinusM1 - k})! * (${j2PlusM2 - k})! * (${JminusJ2PlusM1 + k})! * (${JminusJ1MinusM2 + k})!)`;
      sumTerms.push(term);
    }
    
    const sumExpr = sumTerms.join(' + ').replace(/\+ -/g, '- ');
    
    // Full expression
    // CG = sqrt((2*J+1) * Δ) * sqrt(factorial_product) * sum
    const fullExpr = `sqrt((2*${aJ}+1) * ${triangleExpr}) * sqrt(${factorialProduct}) * (${sumExpr})`;
    
    // Simplify with Algebrite
    const simplified = Algebrite.run(`simplify(${fullExpr})`);

    // Get decimal value
    const decimalStr = Algebrite.run(`float(${simplified})`);
    const decimalVal = parseFloat(decimalStr);

    // Get LaTeX form (best-effort; Algebrite's printlatex sometimes returns
    // identical infix text — that's still usable as basic LaTeX).
    let latex = '';
    try {
      latex = Algebrite.run(`printlatex(${simplified})`).trim();
    } catch (_) {
      latex = '';
    }

    return {
      decimal: decimalVal,
      symbolic: simplified.trim(),
      latex: latex
    };
    
  } catch (e) {
    return { error: `Algebrite computation error: ${e.message || e}` };
  }
}

/**
 * Format symbolic result for display
 */
function formatSymbolic(sym) {
  if (!sym) return '0';
  
  // Clean up the symbolic output
  let result = sym.trim();
  
  // Replace common patterns for better display
  // Algebrite uses ^ for exponentiation, sqrt for square root
  // We'll display it as-is in monospace
  
  return result;
}

// ============================================================================
// SANITY CHECKS / UNIT TESTS
// ============================================================================

function runSanityChecks() {
  console.log('%c=== Clebsch-Gordan Coefficient Sanity Checks ===', 'color: #06b6d4; font-weight: bold; font-size: 14px;');
  console.log('');
  
  if (typeof Algebrite === 'undefined') {
    console.log('%c✗ Algebrite library not loaded — cannot run tests', 'color: #ef4444; font-weight: bold;');
    console.log('');
    return { passed: 0, failed: 1 };
  }
  
  const tests = [
    // ⟨j1 m1 j2 m2 | J M⟩
    { j1: '1/2', m1: '1/2', j2: '1/2', m2: '1/2', J: '1', M: '1', expected: 1, description: '⟨½ ½ ½ ½ | 1 1⟩' },
    { j1: '1/2', m1: '1/2', j2: '1/2', m2: '-1/2', J: '1', M: '0', expected: 1/Math.sqrt(2), description: '⟨½ ½ ½ −½ | 1 0⟩' },
    { j1: '1/2', m1: '-1/2', j2: '1/2', m2: '1/2', J: '1', M: '0', expected: 1/Math.sqrt(2), description: '⟨½ −½ ½ ½ | 1 0⟩' },
    { j1: '1/2', m1: '1/2', j2: '1/2', m2: '-1/2', J: '0', M: '0', expected: 1/Math.sqrt(2), description: '⟨½ ½ ½ −½ | 0 0⟩' },
    { j1: '1/2', m1: '-1/2', j2: '1/2', m2: '1/2', J: '0', M: '0', expected: -1/Math.sqrt(2), description: '⟨½ −½ ½ ½ | 0 0⟩' },
    { j1: '1', m1: '0', j2: '1', m2: '0', J: '2', M: '0', expected: Math.sqrt(2/3), description: '⟨1 0 1 0 | 2 0⟩' },
    { j1: '1', m1: '1', j2: '1', m2: '-1', J: '0', M: '0', expected: 1/Math.sqrt(3), description: '⟨1 1 1 −1 | 0 0⟩' },
    { j1: '1', m1: '1', j2: '1', m2: '0', J: '2', M: '1', expected: 1/Math.sqrt(2), description: '⟨1 1 1 0 | 2 1⟩' },
    { j1: '1', m1: '0', j2: '1', m2: '0', J: '0', M: '0', expected: -1/Math.sqrt(3), description: '⟨1 0 1 0 | 0 0⟩' },
    { j1: '3/2', m1: '1/2', j2: '1/2', m2: '1/2', J: '2', M: '1', expected: 1, description: '⟨3/2 ½ ½ ½ | 2 1⟩' },
  ];
  
  let passed = 0;
  let failed = 0;
  
  tests.forEach((test, index) => {
    const result = computeCG(test.j1, test.m1, test.j2, test.m2, test.J, test.M);
    
    if (result.error) {
      console.log(`%c✗ Test ${index + 1}: ${test.description}`, 'color: #ef4444;');
      console.log(`  Error: ${result.error}`);
      failed++;
      return;
    }
    
    const diff = Math.abs(result.decimal - test.expected);
    const tolerance = 1e-10;
    
    if (diff < tolerance) {
      console.log(`%c✓ Test ${index + 1}: ${test.description}`, 'color: #10b981;');
      console.log(`  Symbolic: ${result.symbolic}`);
      console.log(`  Decimal: ${result.decimal.toFixed(12)}`);
      console.log(`  Expected: ${test.expected.toFixed(12)}`);
      passed++;
    } else {
      console.log(`%c✗ Test ${index + 1}: ${test.description}`, 'color: #ef4444;');
      console.log(`  Symbolic: ${result.symbolic}`);
      console.log(`  Decimal: ${result.decimal.toFixed(12)}`);
      console.log(`  Expected: ${test.expected.toFixed(12)}`);
      console.log(`  Difference: ${diff.toExponential(3)}`);
      failed++;
    }
    console.log('');
  });
  
  console.log('%c========================================', 'color: #06b6d4;');
  console.log(`%cResults: ${passed} passed, ${failed} failed`, 
    failed === 0 ? 'color: #10b981; font-weight: bold;' : 'color: #ef4444; font-weight: bold;');
  console.log('');
  
  return { passed, failed };
}

// Export for use in other modules
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    computeCG,
    checkSelectionRules,
    parseQuantumNumber,
    runSanityChecks
  };
}
