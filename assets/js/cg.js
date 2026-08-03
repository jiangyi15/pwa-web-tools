/**
 * Clebsch-Gordan Coefficient Calculator
 * Exact BigInt arithmetic via the Racah formula
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
  const g = Surd._gcd(Math.abs(num), den);
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

/**
 * Check that m is a valid projection of j: both j+m and j−m must be
 * non-negative integers. This is equivalent to requiring m to be an
 * integer when j is integer, and a half-integer when j is half-integer
 * (m = −j, −j+1, …, j). Rejects e.g. m = 0 for j = 1/2.
 */
function _validProjectionParity(j, m) {
  const sum = j + m;
  const diff = j - m;
  return sum >= -1e-10 && diff >= -1e-10 &&
         Math.abs(sum - Math.round(sum)) < 1e-10 &&
         Math.abs(diff - Math.round(diff)) < 1e-10;
}

// ============================================================================
// SELECTION RULES
// ============================================================================

/**
 * Check selection rules for CG coefficients
 * Returns {valid: boolean, message: string, parsed: object}
 */
function checkSelectionRules(j1, m1, j2, m2, J, M) {
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
  // (m must have the same half-integer parity as j — i.e. m = −j,…,+j)
  const checkConsistency = (j, m, jName, mName) => {
    const sum = j + m;
    const diff = j - m;
    if (sum < -1e-10 || Math.abs(sum - Math.round(sum)) > 1e-10) {
      return { valid: false, message: `${jName} + ${mName} = ${sum.toFixed(3)} is not a non-negative integer (m must be integer when j is integer, half-integer when j is half-integer)` };
    }
    if (diff < -1e-10 || Math.abs(diff - Math.round(diff)) > 1e-10) {
      return { valid: false, message: `${jName} − ${mName} = ${diff.toFixed(3)} is not a non-negative integer (m must be integer when j is integer, half-integer when j is half-integer)` };
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
// EXACT CLEBSCH-GORDAN COMPUTATION (BIGINT)
// ============================================================================

/**
 * Helper: Compute gcd of two BigInt values.
 */
function _bigIntGcd(a, b) {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) {
    var t = b;
    b = a % b;
    a = t;
  }
  return a;
}

/**
 * Helper: Check if a number is prime (simple trial division for small n).
 */
function _isPrime(n) {
  if (n < 2) return false;
  if (n === 2) return true;
  if (n % 2 === 0) return false;
  for (var i = 3; i * i <= n; i += 2) {
    if (n % i === 0) return false;
  }
  return true;
}

/**
 * Extract perfect squares from num/den fraction.
 * Returns { outside: BigInt, radicand: Number, denom: BigInt }
 * where sqrt(num/den) = outside * sqrt(radicand) / denom
 */
function _extractSqrt(num, den, maxN) {
  // sqrt(num/den) = sqrt(num*den) / den
  var product = num * den;
  var outside = 1n;
  var radicand = 1;
  var temp = product;
  
  // Trial factor using primes up to maxN
  for (var p = 2; p <= maxN; p++) {
    if (!_isPrime(p)) continue;
    var exp = 0;
    var pBig = BigInt(p);
    while (temp % pBig === 0n) {
      temp = temp / pBig;
      exp++;
    }
    if (exp > 0) {
      // outside *= p^(exp/2)
      var pairs = Math.floor(exp / 2);
      for (var i = 0; i < pairs; i++) {
        outside = outside * pBig;
      }
      // radicand *= p if exp is odd
      if (exp % 2 === 1) {
        radicand *= p;
      }
    }
  }
  
  // Any remaining prime factor in temp (should be rare for CG coefficients)
  if (temp > 1n) {
    // temp is a single prime factor (could be large)
    // We can't extract pairs, so it contributes to radicand
    // But temp might be larger than Number.MAX_SAFE_INTEGER
    // For safety, check if it fits in Number
    if (temp <= BigInt(Number.MAX_SAFE_INTEGER)) {
      radicand *= Number(temp);
    } else {
      // Large prime remaining - this shouldn't happen for typical CG coefficients
      // but we handle it by leaving it as-is (outside stays the same, radicand gets large)
      throw new Error('Large prime factor in CG coefficient sqrt extraction');
    }
  }
  
  return { outside: outside, radicand: radicand, denom: den };
}

/**
 * Compute Clebsch-Gordan coefficient exactly using BigInt arithmetic.
 * Returns a Surd object (sign × p × sqrt(r) / q).
 * 
 * Uses the Racah formula:
 * CG = sqrt( (2J+1) × Δ² × ∏ aᵢ! )  ×  Σₖ (-1)ᵏ / C(k)
 * 
 * No external dependencies - pure BigInt arithmetic.
 */
function computeCGExact(j1, m1, j2, m2, J, M) {
  // Parse all inputs
  var pj1, pm1, pj2, pm2, pJ, pM;
  
  try {
    pj1 = parseQuantumNumber(j1);
    pm1 = parseQuantumNumber(m1);
    pj2 = parseQuantumNumber(j2);
    pm2 = parseQuantumNumber(m2);
    pJ = parseQuantumNumber(J);
    pM = parseQuantumNumber(M);
  } catch (e) {
    return Surd.ZERO;
  }
  
  // Validate quantum numbers
  if (!isValidQuantumNumber(pj1)) return Surd.ZERO;
  if (!isValidQuantumNumber(pj2)) return Surd.ZERO;
  if (!isValidQuantumNumber(pJ)) return Surd.ZERO;
  if (!isValidQuantumNumber(pm1, true)) return Surd.ZERO;
  if (!isValidQuantumNumber(pm2, true)) return Surd.ZERO;
  if (!isValidQuantumNumber(pM, true)) return Surd.ZERO;
  
  var j1Val = pj1.value;
  var m1Val = pm1.value;
  var j2Val = pj2.value;
  var m2Val = pm2.value;
  var JVal = pJ.value;
  var MVal = pM.value;
  
  // Check triangle inequality
  if (JVal < Math.abs(j1Val - j2Val) - 1e-10 || JVal > j1Val + j2Val + 1e-10) {
    return Surd.ZERO;
  }
  
  // Check |m| ≤ j
  if (Math.abs(m1Val) - j1Val > 1e-10) return Surd.ZERO;
  if (Math.abs(m2Val) - j2Val > 1e-10) return Surd.ZERO;
  if (Math.abs(MVal) - JVal > 1e-10) return Surd.ZERO;
  
  // Check m1 + m2 = M
  if (Math.abs(m1Val + m2Val - MVal) > 1e-10) return Surd.ZERO;

  // Check that each m is a valid projection of its j: j±m must be
  // non-negative integers (rejects e.g. m = 0 for j = 1/2, which would
  // otherwise be silently rounded to a nonzero but wrong coefficient).
  if (!_validProjectionParity(j1Val, m1Val)) return Surd.ZERO;
  if (!_validProjectionParity(j2Val, m2Val)) return Surd.ZERO;
  if (!_validProjectionParity(JVal, MVal)) return Surd.ZERO;
  
  // Compute integer values (always integers for these expressions)
  var j1pm1 = Math.round(j1Val + m1Val);
  var j1mm1 = Math.round(j1Val - m1Val);
  var j2pm2 = Math.round(j2Val + m2Val);
  var j2mm2 = Math.round(j2Val - m2Val);
  var JpM = Math.round(JVal + MVal);
  var JmM = Math.round(JVal - MVal);
  
  var j1pj2mJ = Math.round(j1Val + j2Val - JVal);
  var j1mj2pJ = Math.round(j1Val - j2Val + JVal);
  var mj1pj2pJ = Math.round(-j1Val + j2Val + JVal);
  var j1pj2pJp1 = Math.round(j1Val + j2Val + JVal + 1);
  
  var Jmj2pm1 = Math.round(JVal - j2Val + m1Val);
  var Jmj1mm2 = Math.round(JVal - j1Val - m2Val);
  
  // Check all values are non-negative
  if (j1pm1 < 0 || j1mm1 < 0 || j2pm2 < 0 || j2mm2 < 0 ||
      JpM < 0 || JmM < 0 || j1pj2mJ < 0 || j1mj2pJ < 0 ||
      mj1pj2pJ < 0 || j1pj2pJp1 < 0) {
    return Surd.ZERO;
  }
  
  // Compute k bounds
  var kMin = Math.max(0, -Jmj2pm1, -Jmj1mm2);
  var kMax = Math.min(j1pj2mJ, j1mm1, j2pm2);
  
  if (kMin > kMax) {
    return Surd.ZERO;
  }
  
  // Precompute factorials up to max needed
  var maxN = Math.max(j1pm1, j1mm1, j2pm2, j2mm2, JpM, JmM,
                      j1pj2mJ, j1mj2pJ, mj1pj2pJ, j1pj2pJp1,
                      kMax);
  var facts = [1n];
  for (var i = 1; i <= maxN; i++) {
    facts[i] = facts[i - 1] * BigInt(i);
  }
  
  // Compute sqrt part numerator and denominator as BigInt
  // N = (2J+1) × (j1pj2mJ)! × (j1mj2pJ)! × (mj1pj2pJ)! × (j1pm1)! × (j1mm1)! × (j2pm2)! × (j2mm2)! × (JpM)! × (JmM)!
  // D = (j1pj2pJp1)!
  var twoJplus1 = Math.round(2 * JVal + 1);
  var N_big = BigInt(twoJplus1);
  N_big = N_big * facts[j1pj2mJ];
  N_big = N_big * facts[j1mj2pJ];
  N_big = N_big * facts[mj1pj2pJ];
  N_big = N_big * facts[j1pm1];
  N_big = N_big * facts[j1mm1];
  N_big = N_big * facts[j2pm2];
  N_big = N_big * facts[j2mm2];
  N_big = N_big * facts[JpM];
  N_big = N_big * facts[JmM];
  
  var D_big = facts[j1pj2pJp1];
  
  // Extract sqrt part
  var sqrtExtract;
  try {
    sqrtExtract = _extractSqrt(N_big, D_big, maxN);
  } catch (e) {
    return Surd.ZERO;
  }
  
  var sqrt_outside = sqrtExtract.outside;
  var sqrt_radicand = sqrtExtract.radicand;
  var sqrt_denom = sqrtExtract.denom;
  
  // Compute sum over k
  var sum_num = 0n;
  var sum_den = 1n;
  
  for (var k = kMin; k <= kMax; k++) {
    var denom_k = facts[k];
    denom_k = denom_k * facts[j1pj2mJ - k];
    denom_k = denom_k * facts[j1mm1 - k];
    denom_k = denom_k * facts[j2pm2 - k];
    denom_k = denom_k * facts[Jmj2pm1 + k];
    denom_k = denom_k * facts[Jmj1mm2 + k];
    
    var sign_k = (k % 2 === 0) ? 1n : -1n;
    
    // sum = sum + sign_k / denom_k
    // = sum_num/sum_den + sign_k/denom_k
    // = (sum_num * denom_k + sign_k * sum_den) / (sum_den * denom_k)
    sum_num = sum_num * denom_k + sign_k * sum_den;
    sum_den = sum_den * denom_k;
  }
  
  // Reduce sum by gcd
  if (sum_num === 0n) return Surd.ZERO;
  
  var sumSign = 1;
  if (sum_num < 0n) {
    sum_num = -sum_num;
    sumSign = -1;
  }
  
  var gcd = _bigIntGcd(sum_num, sum_den);
  sum_num = sum_num / gcd;
  sum_den = sum_den / gcd;
  
  // Combine sqrt part and sum
  // sqrt_outside × sqrt(radicand) / sqrt_denom × sum_num / sum_den
  // = (sqrt_outside × sum_num) × sqrt(radicand) / (sqrt_denom × sum_den)
  var final_numer = sqrt_outside * sum_num;
  var final_denom = sqrt_denom * sum_den;
  
  // Reduce final fraction
  gcd = _bigIntGcd(final_numer, final_denom);
  final_numer = final_numer / gcd;
  final_denom = final_denom / gcd;
  
  // Check if values fit in Number (<= 2^53)
  var MAX_SAFE = BigInt(Number.MAX_SAFE_INTEGER);
  if (final_numer > MAX_SAFE || final_denom > MAX_SAFE || sqrt_radicand > Number.MAX_SAFE_INTEGER) {
    throw new Error('CG result too large for Surd');
  }
  
  return new Surd(sumSign, Number(final_numer), Number(final_denom), sqrt_radicand);
}

// ============================================================================
// CLEBSCH-GORDAN CALCULATION (exact BigInt, no external deps)
// ============================================================================

/**
 * Compute Clebsch-Gordan coefficient using exact BigInt arithmetic.
 * Returns { decimal, symbolic, latex, error }
 */
function computeCG(j1, m1, j2, m2, J, M) {
  // Check selection rules
  const rules = checkSelectionRules(j1, m1, j2, m2, J, M);
  if (!rules.valid) {
    return { error: rules.message };
  }
  
  try {
    const surd = computeCGExact(j1, m1, j2, m2, J, M);
    
    if (surd.isZero()) {
      return {
        decimal: 0,
        symbolic: '0',
        latex: '0'
      };
    }
    
    const decimal = surd.s * surd.p * Math.sqrt(surd.r) / surd.q;
    const symbolic = surd.toString();
    const latex = surd.toLatex();
    
    return {
      decimal: decimal,
      symbolic: symbolic,
      latex: latex
    };
    
  } catch (e) {
    return { error: `Computation error: ${e.message || e}` };
  }
}

/**
 * Format symbolic result for display
 */
function formatSymbolic(sym) {
  if (!sym) return '0';
  
  // Clean up the symbolic output
  let result = sym.trim();
  return result;
}

// ============================================================================
// SANITY CHECKS / UNIT TESTS
// ============================================================================

function runSanityChecks() {
  console.log('%c=== Clebsch-Gordan Coefficient Sanity Checks ===', 'color: #06b6d4; font-weight: bold; font-size: 14px;');
  console.log('');
  
  if (typeof computeCGExact === 'undefined') {
    console.log('%c✗ computeCGExact not available — cannot run tests', 'color: #ef4444; font-weight: bold;');
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
    { j1: '3/2', m1: '1/2', j2: '1/2', m2: '1/2', J: '2', M: '1', expected: Math.sqrt(3)/2, description: '⟨3/2 ½ ½ ½ | 2 1⟩' },
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

/**
 * Run sanity checks for exact CG computation
 */
function runSanityChecksExact() {
  console.log('%c=== Exact CG Computation Sanity Checks ===', 'color: #06b6d4; font-weight: bold; font-size: 14px;');
  console.log('');
  
  var tests = [
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
    { j1: '3/2', m1: '1/2', j2: '1/2', m2: '1/2', J: '2', M: '1', expected: Math.sqrt(3)/2, description: '⟨3/2 ½ ½ ½ | 2 1⟩' },
  ];
  
  var passed = 0;
  var failed = 0;
  
  tests.forEach(function(test, index) {
    var result = computeCGExact(test.j1, test.m1, test.j2, test.m2, test.J, test.M);
    
    if (result.isZero() && test.expected !== 0) {
      console.log('%c✗ Test ' + (index + 1) + ': ' + test.description, 'color: #ef4444;');
      console.log('  Result is zero, expected non-zero');
      failed++;
      return;
    }
    
    // Compute decimal value from Surd
    var decimal = result.s * result.p * Math.sqrt(result.r) / result.q;
    var diff = Math.abs(decimal - test.expected);
    var tolerance = 1e-10;
    
    if (diff < tolerance) {
      console.log('%c✓ Test ' + (index + 1) + ': ' + test.description, 'color: #10b981;');
      console.log('  Surd: ' + result.toString());
      console.log('  Decimal: ' + decimal.toFixed(12));
      console.log('  Expected: ' + test.expected.toFixed(12));
      passed++;
    } else {
      console.log('%c✗ Test ' + (index + 1) + ': ' + test.description, 'color: #ef4444;');
      console.log('  Surd: ' + result.toString());
      console.log('  Decimal: ' + decimal.toFixed(12));
      console.log('  Expected: ' + test.expected.toFixed(12));
      console.log('  Difference: ' + diff.toExponential(3));
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
    computeCG,
    computeCGExact,
    checkSelectionRules,
    parseQuantumNumber,
    runSanityChecks,
    runSanityChecksExact
  };
}
