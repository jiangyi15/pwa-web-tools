/**
 * Get Angle Formula - Helicity Amplitude Angular Formula Calculator
 * 
 * Uses (real, imag) structure for complex numbers to avoid sqrt(-1) issues.
 * Complex operations:
 *   (r1,i1) + (r2,i2) = (r1+r2, i1+i2)
 *   (r1,i1) * (r2,i2) = (r1*r2 - i1*i2, r1*i2 + i1*r2)
 */

// ============================================================================
// COMPLEX NUMBER OPERATIONS
// ============================================================================

/**
 * Create complex number: real + i*imag
 */
function complex(real, imag) {
  return { real: real, imag: imag };
}

/**
 * Complex from real number: r + i*0
 */
function complexReal(r) {
  return { real: r, imag: '0' };
}

/**
 * Complex add: (r1,i1) + (r2,i2)
 */
function complexAdd(c1, c2) {
  return {
    real: '(' + c1.real + ') + (' + c2.real + ')',
    imag: '(' + c1.imag + ') + (' + c2.imag + ')'
  };
}

/**
 * Complex multiply: (r1,i1) * (r2,i2) = (r1*r2 - i1*i2, r1*i2 + i1*r2)
 */
function complexMul(c1, c2) {
  return {
    real: '(' + c1.real + ')*(' + c2.real + ') - (' + c1.imag + ')*(' + c2.imag + ')',
    imag: '(' + c1.real + ')*(' + c2.imag + ') + (' + c1.imag + ')*(' + c2.real + ')'
  };
}

/**
 * Complex scale: s * (r, i)
 */
function complexScale(s, c) {
  return {
    real: '(' + s + ')*(' + c.real + ')',
    imag: '(' + s + ')*(' + c.imag + ')'
  };
}

/**
 * Complex conjugate: (r, i)* = (r, -i)
 */
function complexConj(c) {
  return { real: c.real, imag: '-(' + c.imag + ')' };
}

/**
 * Simplify complex number (both parts)
 * Uses expand then simplify to help with trigonometric cancellations
 */
function complexSimplify(c) {
  try {
    var realExpanded = Algebrite.run('expand(' + c.real + ')').trim();
    var imagExpanded = Algebrite.run('expand(' + c.imag + ')').trim();
    return {
      real: Algebrite.run('simplify(' + realExpanded + ')').trim(),
      imag: Algebrite.run('simplify(' + imagExpanded + ')').trim()
    };
  } catch (e) {
    return c;
  }
}

/**
 * Convert complex to display string
 */
function complexToString(c) {
  var r = c.real;
  var i = c.imag;
  if (i === '0') return r;
  if (r === '0') return i + '*i';
  return '(' + r + ' + ' + i + '*i)';
}

/**
 * Convert complex to LaTeX
 */
function complexToLatex(c) {
  var realPart = c.real;
  var imagPart = c.imag;
  
  // Expand and simplify real part to LaTeX
  try {
    realPart = Algebrite.run('simplify(expand(' + realPart + '))').trim();
    realPart = Algebrite.run('printlatex(' + realPart + ')').trim();
    realPart = realPart.replace(/(?<!\\)cos/g, '\\cos').replace(/(?<!\\)sin/g, '\\sin');
    realPart = realPart.replace(/(?<!\\)phi/g, '\\phi').replace(/(?<!\\)theta/g, '\\theta');
  } catch(e) {}
  
  // Check if imag is zero
  if (imagPart === '0' || imagPart === '') {
    return realPart;
  }
  
  // Expand and simplify imag part to LaTeX
  try {
    imagPart = Algebrite.run('simplify(expand(' + imagPart + '))').trim();
    imagPart = Algebrite.run('printlatex(' + imagPart + ')').trim();
    imagPart = imagPart.replace(/(?<!\\)cos/g, '\\cos').replace(/(?<!\\)sin/g, '\\sin');
    imagPart = imagPart.replace(/(?<!\\)phi/g, '\\phi').replace(/(?<!\\)theta/g, '\\theta');
  } catch(e) {}
  
  // Check again if imag became 0 after simplification
  if (imagPart === '0' || imagPart === '') {
    return realPart;
  }
  
  // Purely imaginary: (imagPart)i
  if (realPart === '0' || realPart === '') {
    return '(' + imagPart + ')i';
  }
  
  // Complex: realPart + (imagPart)i
  return realPart + ' + (' + imagPart + ')i';
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

function toSpin(val) {
  if (typeof val === 'object' && val.num !== undefined) {
    return { num: val.num, den: val.den || 1, value: val.num / (val.den || 1) };
  }
  if (typeof val === 'number') {
    var twice = val * 2;
    if (Math.abs(twice - Math.round(twice)) < 1e-10) {
      return { num: Math.round(twice), den: 2, value: val };
    }
    return { num: Math.round(val), den: 1, value: val };
  }
  if (typeof val === 'string') {
    if (val.includes('/')) {
      var parts = val.split('/');
      var num = parseInt(parts[0], 10);
      var den = parseInt(parts[1], 10);
      return { num: num, den: den, value: num / den };
    }
    var n = parseInt(val, 10);
    if (isNaN(n)) throw new Error('Invalid spin value: ' + val);
    return { num: n, den: 1, value: n };
  }
  throw new Error('Invalid spin value: ' + val);
}

function formatSpin(s) {
  if (s.den === 1) return String(s.num);
  if (s.num === 1 && s.den === 2) return '1/2';
  if (s.num === -1 && s.den === 2) return '-1/2';
  return s.num + '/' + s.den;
}

function getHelicities(spin) {
  var s = toSpin(spin);
  // For spin = num/den, the helicity values range from -spin to +spin in steps of 1
  // Using 2*spin as integer for clean arithmetic
  var spin2 = 2 * s.value; // This is always an integer
  spin2 = Math.round(spin2); // Ensure it's exactly integer
  
  var values = [];
  // h2 goes from -2*spin to +2*spin in steps of 2
  // So helicity = h2/2 ranges from -spin to +spin in steps of 1
  for (var h2 = -spin2; h2 <= spin2; h2 += 2) {
    values.push({ num: h2, den: 2, value: h2 / 2 });
  }
  return values;
}

function getLSCombinations(Ja, Jb, Jc) {
  var a = toSpin(Ja);
  var b = toSpin(Jb);
  var c = toSpin(Jc);
  
  // Use 2*J as integers for clean arithmetic
  var Ja2 = Math.round(2 * a.value);
  var Jb2 = Math.round(2 * b.value);
  var Jc2 = Math.round(2 * c.value);
  
  var s2_min = Math.abs(Jb2 - Jc2);
  var s2_max = Jb2 + Jc2;
  
  var pairs = [];
  
  for (var s2 = s2_min; s2 <= s2_max; s2 += 2) {
    var l_min_2 = Math.abs(Ja2 - s2);
    var l_max_2 = Ja2 + s2;
    var l_min = Math.ceil(l_min_2 / 2);
    var l_max = Math.floor(l_max_2 / 2);
    
    for (var l = l_min; l <= l_max; l++) {
      if (Math.abs(2 * l - s2) <= Ja2 && Ja2 <= 2 * l + s2) {
        pairs.push({ l: l, s: s2 / 2 });
      }
    }
  }
  
  return pairs;
}

// ============================================================================
// WIGNER D-FUNCTION (returns complex {real, imag})
// ============================================================================

/**
 * Compute Wigner D-conjugate: D^{J*}_{m,mp}(phi, theta)
 * D^{J*}_{m,mp} = exp(i*m*phi) * d^J_{m,mp}(theta)
 * 
 * Returns {real, imag} complex structure.
 * exp(i*m*phi) = cos(m*phi) + i*sin(m*phi)
 * So D* = cos(m*phi)*d + i*sin(m*phi)*d
 */
function computeWignerDConj(J, m, mp, phiVar, thetaVar) {
  J = toSpin(J).value;
  m = toSpin(m).value;
  mp = toSpin(mp).value;
  
  // Compute d-matrix element (real for real theta)
  var dResult = computeWignerD(String(J), String(m), String(mp), thetaVar);
  if (dResult.error) return null;
  if (dResult.decimal === 0 && dResult.symbolic === '0') return complex('0', '0');
  
  var dExpr = dResult.symbolic;
  
  // If m = 0, exp(i*0*phi) = 1, so D* = d (purely real)
  if (Math.abs(m) < 1e-10) {
    return complex(dExpr, '0');
  }
  
  // exp(i*m*phi) = cos(m*phi) + i*sin(m*phi)
  var mStr = String(m);
  var cosPart = 'cos((' + mStr + ')*' + phiVar + ')';
  var sinPart = 'sin((' + mStr + ')*' + phiVar + ')';
  
  // D* = (cos + i*sin) * d = (cos*d, sin*d)
  return {
    real: '(' + cosPart + ')*(' + dExpr + ')',
    imag: '(' + sinPart + ')*(' + dExpr + ')'
  };
}

// ============================================================================
// HELICITY AMPLITUDE (returns complex {real, imag})
// ============================================================================

/**
 * Compute helicity amplitude for a specific (l, s) combination.
 * Returns {real, imag} complex structure.
 */
function computeHelicityAmp(Ja, Jb, Jc, la, lb, lc, phiVar, thetaVar, l, s) {
  Ja = toSpin(Ja).value;
  Jb = toSpin(Jb).value;
  Jc = toSpin(Jc).value;
  
  var delta = lb - lc;
  if (Math.abs(delta) > Ja + 1e-10) return null;
  
  // CG1: <Jb, lb; Jc, -lc | s, delta> (real)
  var cg1 = computeCG(
    formatSpin(toSpin(Jb)), formatSpin(toSpin(lb)),
    formatSpin(toSpin(Jc)), formatSpin(toSpin(-lc)),
    formatSpin(toSpin(s)), formatSpin(toSpin(delta))
  );
  if (cg1.error || (cg1.decimal === 0 && cg1.symbolic === '0')) return null;
  
  // CG2: <l, 0; s, delta | Ja, delta> (real)
  var cg2 = computeCG(
    String(l), '0',
    formatSpin(toSpin(s)), formatSpin(toSpin(delta)),
    formatSpin(toSpin(Ja)), formatSpin(toSpin(delta))
  );
  if (cg2.error || (cg2.decimal === 0 && cg2.symbolic === '0')) return null;
  
  // Coefficient sqrt((2l+1)/(2Ja+1)) (real)
  var coeff = 'sqrt(' + (2*l + 1) + '/' + (Math.round(2*Ja) + 1) + ')';
  
  // Full real coefficient: coeff * cg1 * cg2
  var realCoeff = '(' + coeff + ')*(' + cg1.symbolic + ')*(' + cg2.symbolic + ')';
  
  // Wigner D conjugate (complex)
  var Dconj = computeWignerDConj(Ja, la, delta, phiVar, thetaVar);
  if (!Dconj) return null;
  
  // Amplitude = realCoeff * Dconj
  var result = complexScale(realCoeff, Dconj);
  
  return complexSimplify(result);
}

// ============================================================================
// DECAY STRUCTURE
// ============================================================================

function countDecayVertices(node) {
  if (!node.children) return 0;
  return 1 + countDecayVertices(node.children[0]) + countDecayVertices(node.children[1]);
}

function getVertices(node, idx) {
  if (idx === undefined) idx = { v: 0 };
  if (!node.children) return [];
  
  var vertices = [{
    index: idx.v,
    Ja: node.j,
    Jb: node.children[0].j,
    Jc: node.children[1].j,
    child0: node.children[0],
    child1: node.children[1]
  }];
  idx.v++;
  
  vertices = vertices.concat(getVertices(node.children[0], idx));
  vertices = vertices.concat(getVertices(node.children[1], idx));
  
  return vertices;
}

// ============================================================================
// MAIN COMPUTATION
// ============================================================================

/**
 * Build amplitude for a single decay vertex.
 * Returns {ls_key: {helicity_key: {real, imag}}}
 */
function buildVertexAmp(vertex, phiVar, thetaVar) {
  var Ja = toSpin(vertex.Ja);
  var Jb = toSpin(vertex.Jb);
  var Jc = toSpin(vertex.Jc);
  
  var las = getHelicities(Ja);
  var lbs = getHelicities(Jb);
  var lcs = getHelicities(Jc);
  var lsPairs = getLSCombinations(Ja, Jb, Jc);
  
  var amp = {};
  
  for (var i = 0; i < lsPairs.length; i++) {
    var ls = lsPairs[i];
    var lsKey = ls.l + ',' + ls.s;
    amp[lsKey] = {};
    
    for (var ia = 0; ia < las.length; ia++) {
      for (var ib = 0; ib < lbs.length; ib++) {
        for (var ic = 0; ic < lcs.length; ic++) {
          var la = las[ia].value;
          var lb = lbs[ib].value;
          var lc = lcs[ic].value;
          
          var result = computeHelicityAmp(Ja, Jb, Jc, la, lb, lc, phiVar, thetaVar, ls.l, ls.s);
          
          if (result && (result.real !== '0' || result.imag !== '0')) {
            var hKey = la + ',' + lb + ',' + lc;
            amp[lsKey][hKey] = result;
          }
        }
      }
    }
    
    if (Object.keys(amp[lsKey]).length === 0) {
      delete amp[lsKey];
    }
  }
  
  return amp;
}

/**
 * Main function: compute angle formula for a decay tree.
 */
function getAngleFormula(decayTree) {
  if (typeof Algebrite === 'undefined') {
    return { error: 'Algebrite library is not loaded.' };
  }
  
  if (!decayTree.children) {
    return { error: 'Decay tree must have at least one decay.' };
  }
  
  var vertices = getVertices(decayTree);
  var nDecays = vertices.length;
  
  // Build amplitude for each vertex
  var amps = [];
  for (var i = 0; i < vertices.length; i++) {
    var phiVar = 'phi_' + i;
    var thetaVar = 'theta_' + i;
    amps.push({
      vertex: vertices[i],
      phi: phiVar,
      theta: thetaVar,
      amp: buildVertexAmp(vertices[i], phiVar, thetaVar)
    });
  }
  
  // Combine amplitudes
  var result = combineAmps(amps, decayTree);
  
  // Build angle list
  var angles = [];
  for (var i = 0; i < nDecays; i++) {
    angles.push('phi_' + i + ', theta_' + i);
  }
  
  return {
    helicities: result,
    nDecays: nDecays,
    angles: angles
  };
}

// ============================================================================
// AMPLITUDE COMBINATION (complex multiplication)
// ============================================================================

function combineAmps(amps, tree) {
  if (amps.length === 1) {
    return formatSingleAmp(amps[0].amp);
  }
  return combineRecursive(amps, tree, 0);
}

function formatSingleAmp(amp) {
  var result = {};
  for (var ls in amp) {
    for (var h in amp[ls]) {
      if (!(h in result)) result[h] = {};
      result[h][ls] = amp[ls][h];
    }
  }
  return result;
}

function combineRecursive(amps, node, idx) {
  var amp = amps[idx];
  var result = {};
  
  var child0Decays = node.children && node.children[0].children;
  var child1Decays = node.children && node.children[1].children;
  
  for (var ls in amp.amp) {
    for (var h in amp.amp[ls]) {
      var hparts = h.split(',').map(parseFloat);
      var la = hparts[0], lb = hparts[1], lc = hparts[2];
      var expr = amp.amp[ls][h]; // {real, imag}
      
      if (child0Decays && child1Decays) {
        var amp0 = combineRecursive(amps, node.children[0], idx + 1);
        var child0Size = countDecayVertices(node.children[0]);
        var amp1 = combineRecursive(amps, node.children[1], idx + 1 + child0Size);
        
        for (var h0 in amp0) {
          var h0parts = h0.split(',').map(parseFloat);
          if (Math.abs(h0parts[0] - lb) > 1e-10) continue;
          
          for (var h1 in amp1) {
            var h1parts = h1.split(',').map(parseFloat);
            if (Math.abs(h1parts[0] - lc) > 1e-10) continue;
            
            for (var ls0 in amp0[h0]) {
              for (var ls1 in amp1[h1]) {
                var fullLs = ls + ';' + ls0 + ';' + ls1;
                var fullH = la + ',' + h0parts.slice(1).join(',') + ',' + h1parts.slice(1).join(',');
                
                // Complex multiplication: expr * amp0 * amp1
                var prod1 = complexMul(expr, amp0[h0][ls0]);
                var combined = complexMul(prod1, amp1[h1][ls1]);
                combined = complexSimplify(combined);
                
                // SUM over intermediate helicities
                if (!(fullH in result)) result[fullH] = {};
                if (!(fullLs in result[fullH])) {
                  result[fullH][fullLs] = combined;
                } else {
                  result[fullH][fullLs] = complexAdd(result[fullH][fullLs], combined);
                  result[fullH][fullLs] = complexSimplify(result[fullH][fullLs]);
                }
              }
            }
          }
        }
      } else if (child0Decays) {
        var amp0 = combineRecursive(amps, node.children[0], idx + 1);
        
        for (var h0 in amp0) {
          var h0parts = h0.split(',').map(parseFloat);
          if (Math.abs(h0parts[0] - lb) > 1e-10) continue;
          
          for (var ls0 in amp0[h0]) {
            var fullLs = ls + ';' + ls0;
            var fullH = la + ',' + h0parts.slice(1).join(',') + ',' + lc;
            
            var combined = complexMul(expr, amp0[h0][ls0]);
            combined = complexSimplify(combined);
            
            // SUM over intermediate helicities
            if (!(fullH in result)) result[fullH] = {};
            if (!(fullLs in result[fullH])) {
              result[fullH][fullLs] = combined;
            } else {
              result[fullH][fullLs] = complexAdd(result[fullH][fullLs], combined);
              result[fullH][fullLs] = complexSimplify(result[fullH][fullLs]);
            }
          }
        }
      } else if (child1Decays) {
        var child0Size = countDecayVertices(node.children[0]);
        var amp1 = combineRecursive(amps, node.children[1], idx + 1 + child0Size);
        
        for (var h1 in amp1) {
          var h1parts = h1.split(',').map(parseFloat);
          if (Math.abs(h1parts[0] - lc) > 1e-10) continue;
          
          for (var ls1 in amp1[h1]) {
            var fullLs = ls + ';' + ls1;
            var fullH = la + ',' + lb + ',' + h1parts.slice(1).join(',');
            
            var combined = complexMul(expr, amp1[h1][ls1]);
            combined = complexSimplify(combined);
            
            // SUM over intermediate helicities
            if (!(fullH in result)) result[fullH] = {};
            if (!(fullLs in result[fullH])) {
              result[fullH][fullLs] = combined;
            } else {
              result[fullH][fullLs] = complexAdd(result[fullH][fullLs], combined);
              result[fullH][fullLs] = complexSimplify(result[fullH][fullLs]);
            }
          }
        }
      } else {
        var fullH = la + ',' + lb + ',' + lc;
        if (!(fullH in result)) result[fullH] = {};
        result[fullH][ls] = expr;
      }
    }
  }
  
  return result;
}

// ============================================================================
// LATEX OUTPUT
// ============================================================================

function resultToLatex(result) {
  var lines = [];
  lines.push('\\begin{array}{lll}');
  lines.push('\\hline');
  lines.push('\\lambda & LS & T(\\phi_i, \\theta_i) \\\\');
  lines.push('\\hline');
  
  var hs = Object.keys(result.helicities).sort();
  for (var i = 0; i < hs.length; i++) {
    var h = hs[i];
    var lsDict = result.helicities[h];
    var lss = Object.keys(lsDict).sort();
    
    for (var j = 0; j < lss.length; j++) {
      var ls = lss[j];
      var c = lsDict[ls]; // {real, imag}
      
      var latexExpr = complexToLatex(c);
      
      var hLatex = '(' + h.split(',').map(function(x) {
        var v = parseFloat(x);
        var s = toSpin(v);
        if (s.den === 2) return '\\tfrac{' + s.num + '}{2}';
        return String(v);
      }).join(', ') + ')';
      
      var lsLatex = '(' + ls.split(';').map(function(x) {
        return '(' + x + ')';
      }).join(', ') + ')';
      
      lines.push('$' + hLatex + '$ & $' + lsLatex + '$ & $' + latexExpr + '$ \\\\');
    }
  }
  
  lines.push('\\hline');
  lines.push('\\end{array}');
  
  return lines.join('\n');
}

// ============================================================================
// STRUCTURED ANGULAR SIMPLIFICATION
// ============================================================================
// Uses structured {coeff, sinPow, cosPow, phiFunc} terms built from CG and 
// Wigner-d weights (no string parsing). Applies half-angle expansion for theta
// and keeps phi as independent sin/cos products.
// ============================================================================

// --- Structured term constructors ---
// Coefficients are ALGEBRITE STRINGS, never floats.

function _st(coeffStr, thetaIdx, sinPow, cosPow, phiIdx, phiFunc, mVal) {
  return {
    s: coeffStr,         // Algebrite string like "sqrt(3)/3" or "-1/2" or "1"
    ti: thetaIdx, sp: sinPow, cp: cosPow,
    pi: phiIdx, pf: phiFunc || '1', pm: mVal || 0,
    im: false
  };
}

function _sti(coeffStr, thetaIdx, sinPow, cosPow, phiIdx, phiFunc, mVal) {
  var t = _st(coeffStr, thetaIdx, sinPow, cosPow, phiIdx, phiFunc, mVal);
  t.im = true;
  return t;
}

// --- Exact Wigner-d weight as Algebrite string ---
// Returns array of {weightStr, sinPow, cosPow} using exact symbolic expressions.

function _getExactWignerDWeights(J, m1, m2) {
  var twoJ = Math.round(2 * J);
  var jpm1 = Math.round(J + m1), jmm1 = Math.round(J - m1);
  var jpm2 = Math.round(J + m2), jmm2 = Math.round(J - m2);
  
  // Precompute factorials
  var maxN = Math.max(jpm1, jmm1, jpm2, jmm2, twoJ + 2, 10);
  var facts = [1];
  for (var i = 1; i <= maxN; i++) facts[i] = facts[i-1] * i;
  
  // Numerator under sqrt: (j+m1)!(j-m1)!(j+m2)!(j-m2)!
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
    if (jmm1 - k >= 0) { for (var i = 1; i <= jmm1-k; i++) denom *= i; }
    if (jpm2 - k >= 0) { for (var i = 1; i <= jpm2-k; i++) denom *= i; }
    if (Math.round(m1-m2)+k >= 0) { for (var i = 1; i <= Math.round(m1-m2)+k; i++) denom *= i; }
    if (k >= 0) { for (var i = 1; i <= k; i++) denom *= i; }
    
    // Extract perfect squares from numNum
    var p = 1, r = numNum;
    for (var i = 2; i * i <= r; i++) {
      while (r % (i * i) === 0) { r /= (i * i); p *= i; }
    }
    
    // Reduce p/denom
    var g = _gcd(p, denom); p /= g; denom /= g;
    
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

var _gcd = function(a, b) {
  a = Math.abs(a); b = Math.abs(b);
  while (b) { var t = b; b = a % b; a = t; }
  return a;
};

// --- Collect Wigner-d weights for structured computation ---

// --- Build structured amplitude for one vertex ---

/**
 * Build structured terms for all (ls, helicity) combinations of one vertex.
 * Returns:
 *   lsDict: { "l,s": { "la,lb,lc": [structuredTerms] } }
 *   children: [childSpinA, childSpinB]
 */
function _buildVertexStructTerms(vertex, thetaIdx, phiIdx) {
  var Ja = toSpin(vertex.Ja);
  var Jb = toSpin(vertex.Jb);
  var Jc = toSpin(vertex.Jc);
  var las = getHelicities(Ja), lbs = getHelicities(Jb), lcs = getHelicities(Jc);
  var lsPairs = getLSCombinations(Ja, Jb, Jc);
  
  var lsDict = {};
  
  for (var i = 0; i < lsPairs.length; i++) {
    var ls = lsPairs[i], lsKey = ls.l + ',' + ls.s;
    lsDict[lsKey] = {};
    
    for (var ia = 0; ia < las.length; ia++) {
      for (var ib = 0; ib < lbs.length; ib++) {
        for (var ic = 0; ic < lcs.length; ic++) {
          var la = las[ia].value, lb = lbs[ib].value, lc = lcs[ic].value;
          var terms = _getOneHelicityStruct(Ja.value, Jb.value, Jc.value, la, lb, lc, thetaIdx, phiIdx, ls.l, ls.s);
          if (terms && terms.length > 0) {
            lsDict[lsKey][la + ',' + lb + ',' + lc] = terms;
          }
        }
      }
    }
    if (Object.keys(lsDict[lsKey]).length === 0) delete lsDict[lsKey];
  }
  
  return {
    lsDict: lsDict,
    children: [Jb.value, Jc.value],
    childSpins: [Jb, Jc]
  };
}

function _getOneHelicityStruct(Ja, Jb, Jc, la, lb, lc, thetaIdx, phiIdx, l, s) {
  var delta = lb - lc;
  if (Math.abs(delta) > Ja + 1e-10) return null;
  
  var cg1 = computeCG(formatSpin(toSpin(Jb)), formatSpin(toSpin(lb)),
                      formatSpin(toSpin(Jc)), formatSpin(toSpin(-lc)),
                      formatSpin(toSpin(s)), formatSpin(toSpin(delta)));
  if (cg1.error || (Math.abs(cg1.decimal) < 1e-14 && cg1.symbolic === '0')) return null;
  
  var cg2 = computeCG(String(l), '0',
                      formatSpin(toSpin(s)), formatSpin(toSpin(delta)),
                      formatSpin(toSpin(Ja)), formatSpin(toSpin(delta)));
  if (cg2.error || (Math.abs(cg2.decimal) < 1e-14 && cg2.symbolic === '0')) return null;
  
  // Build coefficient string: sqrt((2l+1)/(2Ja+1)) * cg1 * cg2
  var sqrtPart = 'sqrt(' + (2*l + 1) + '/' + (Math.round(2*Ja) + 1) + ')';
  var coeffStr = sqrtPart + '*(' + cg1.symbolic + ')*(' + cg2.symbolic + ')';
  
  // Get exact Wigner-d weight strings
  var wdWeights = _getExactWignerDWeights(Ja, la, delta);
  if (wdWeights.length === 0) return null;
  
  var terms = [];
  for (var i = 0; i < wdWeights.length; i++) {
    var wt = wdWeights[i];
    var termCoeffStr = coeffStr + '*(' + wt.weightStr + ')';
    
    if (Math.abs(la) < 1e-10) {
      terms.push(_st(termCoeffStr, thetaIdx, wt.sinPow, wt.cosPow, phiIdx, '1', 0));
    } else {
      var absLA = Math.abs(la);
      var sinSign = la < 0 ? -1 : 1;
      // For cos: coefficient stays as-is. For sin: multiply by sinSign
      terms.push(_st(termCoeffStr, thetaIdx, wt.sinPow, wt.cosPow, phiIdx, 'cos', absLA));
      var sinStr = sinSign < 0 ? '-(' + termCoeffStr + ')' : termCoeffStr;
      terms.push(_sti(sinStr, thetaIdx, wt.sinPow, wt.cosPow, phiIdx, 'sin', absLA));
    }
  }
  return terms;
}

// --- Cascade combination (structured) ---

/**
 * Combine structured amplitudes across the decay chain.
 * Follows same helicity-sum logic as combineRecursive but with structured terms.
 * Returns: { "fullH": { "fullLs": [mergedStructuredTerms] } }
 */
function _combineStructured(ampStructs, node, idx) {
  var struct = ampStructs[idx];
  var result = {};
  
  var child0Decays = node.children && node.children[0] && node.children[0].children;
  var child1Decays = node.children && node.children[1] && node.children[1].children;
  
  for (var lsKey in struct.lsDict) {
    var hDict = struct.lsDict[lsKey];
    
    for (var hKey in hDict) {
      var hparts = hKey.split(',').map(parseFloat);
      var la = hparts[0], lb = hparts[1], lc = hparts[2];
      var terms = hDict[hKey]; // array of structured terms
      
      if (child0Decays && child1Decays) {
        var child0 = _combineStructured(ampStructs, node.children[0], idx + 1);
        var child0Size = countDecayVertices(node.children[0]);
        var child1 = _combineStructured(ampStructs, node.children[1], idx + 1 + child0Size);
        
        for (var h0 in child0) {
          var h0parts = h0.split(',').map(parseFloat);
          if (Math.abs(h0parts[0] - lb) > 1e-10) continue;
          
          for (var h1 in child1) {
            var h1parts = h1.split(',').map(parseFloat);
            if (Math.abs(h1parts[0] - lc) > 1e-10) continue;
            
            for (var ls0 in child0[h0]) {
              for (var ls1 in child1[h1]) {
                var fullLs = lsKey + ';' + ls0 + ';' + ls1;
                var fullH = la + ',' + h0parts.slice(1).join(',') + ',' + h1parts.slice(1).join(',');
                
                // Complex multiply: terms * child0[h0][ls0] * child1[h1][ls1]
                var prod = _threeWayStructuredMul(terms, child0[h0][ls0], child1[h1][ls1]);
                
                if (!(fullH in result)) result[fullH] = {};
                if (!(fullLs in result[fullH])) {
                  result[fullH][fullLs] = prod;
                } else {
                  result[fullH][fullLs] = result[fullH][fullLs].concat(prod);
                }
              }
            }
          }
        }
      } else if (child0Decays) {
        var child0 = _combineStructured(ampStructs, node.children[0], idx + 1);
        for (var h0 in child0) {
          var h0parts = h0.split(',').map(parseFloat);
          if (Math.abs(h0parts[0] - lb) > 1e-10) continue;
          
          for (var ls0 in child0[h0]) {
            var fullLs = lsKey + ';' + ls0;
            var fullH = la + ',' + h0parts.slice(1).join(',') + ',' + lc;
            var prod = _twoWayStructuredMul(terms, child0[h0][ls0]);
            
            if (!(fullH in result)) result[fullH] = {};
            if (!(fullLs in result[fullH])) {
              result[fullH][fullLs] = prod;
            } else {
              result[fullH][fullLs] = result[fullH][fullLs].concat(prod);
            }
          }
        }
      } else if (child1Decays) {
        var child0Size = countDecayVertices(node.children[0]);
        var child1 = _combineStructured(ampStructs, node.children[1], idx + 1 + child0Size);
        for (var h1 in child1) {
          var h1parts = h1.split(',').map(parseFloat);
          if (Math.abs(h1parts[0] - lc) > 1e-10) continue;
          
          for (var ls1 in child1[h1]) {
            var fullLs = lsKey + ';' + ls1;
            var fullH = la + ',' + lb + ',' + h1parts.slice(1).join(',');
            var prod = _twoWayStructuredMul(terms, child1[h1][ls1]);
            
            if (!(fullH in result)) result[fullH] = {};
            if (!(fullLs in result[fullH])) {
              result[fullH][fullLs] = prod;
            } else {
              result[fullH][fullLs] = result[fullH][fullLs].concat(prod);
            }
          }
        }
      } else {
        // Neither child decays - direct output
        var fullH = la + ',' + lb + ',' + lc;
        if (!(fullH in result)) result[fullH] = {};
        if (!(lsKey in result[fullH])) {
          result[fullH][lsKey] = terms.slice(); // clone
        } else {
          result[fullH][lsKey] = result[fullH][lsKey].concat(terms);
        }
      }
    }
  }
  
  return result;
}

/**
 * Multiply 3 structured term arrays: (A) * (B) * (C)
 * Complex multiplication tracking isImag flag
 */
function _threeWayStructuredMul(termsA, termsB, termsC) {
  var result = [];
  for (var a = 0; a < termsA.length; a++) {
    for (var b = 0; b < termsB.length; b++) {
      for (var c = 0; c < termsC.length; c++) {
        var t = _mulThree(termsA[a], termsB[b], termsC[c]);
        if (t && t.length) Array.prototype.push.apply(result, t);
      }
    }
  }
  return result;
}

function _twoWayStructuredMul(termsA, termsB) {
  var result = [];
  for (var a = 0; a < termsA.length; a++) {
    for (var b = 0; b < termsB.length; b++) {
      var t = _mulTwo(termsA[a], termsB[b]);
      if (t && t.length) Array.prototype.push.apply(result, t);
    }
  }
  return result;
}

function _mulTwo(ta, tb) {
  var theta = _mergeThetaPair(ta, tb);
  var phi = _mergePhiPair(ta, tb);
  
  var result = [];
  // (ar + i*ai)*(br + i*bi) = (ar*br - ai*bi) + i*(ar*bi + ai*br)
  // Coefficient: concatenate as Algebrite product strings
  
  if (!ta.im && !tb.im) {
    result.push({ s: '(' + ta.s + ')*(' + tb.s + ')', _tbatch: theta._tbatch, _pbatch: phi._pbatch, im: false });
  }
  if (ta.im && tb.im) {
    result.push({ s: '-(' + ta.s + ')*(' + tb.s + ')', _tbatch: theta._tbatch, _pbatch: phi._pbatch, im: false });
  }
  if (!ta.im && tb.im) {
    result.push({ s: '(' + ta.s + ')*(' + tb.s + ')', _tbatch: theta._tbatch, _pbatch: phi._pbatch, im: true });
  }
  if (ta.im && !tb.im) {
    result.push({ s: '(' + ta.s + ')*(' + tb.s + ')', _tbatch: theta._tbatch, _pbatch: phi._pbatch, im: true });
  }
  return result;
}

function _mulThree(ta, tb, tc) {
  // Multiply two then multiply by third
  var ab = _mulTwo(ta, tb);
  var result = [];
  for (var i = 0; i < ab.length; i++) {
    var ac = _mulTwo(ab[i], tc);
    for (var j = 0; j < ac.length; j++) {
      result.push(ac[j]);
    }
  }
  return result;
}

function _collectTheta(term) {
  var batch = [];
  if (term._tbatch) {
    for (var i = 0; i < term._tbatch.length; i++) {
      batch.push({ idx: term._tbatch[i].idx, sp: term._tbatch[i].sp, cp: term._tbatch[i].cp });
    }
  } else if (term.ti !== undefined) {
    batch.push({ idx: term.ti, sp: term.sp || 0, cp: term.cp || 0 });
  }
  return batch;
}

function _collectPhi(term) {
  var batch = [];
  if (term._pbatch) {
    for (var i = 0; i < term._pbatch.length; i++) {
      batch.push({ idx: term._pbatch[i].idx, pf: term._pbatch[i].pf, pm: term._pbatch[i].pm || 0 });
    }
  } else if (term.pi !== undefined && term.pi >= 0 && term.pf && term.pf !== '1') {
    batch.push({ idx: term.pi, pf: term.pf, pm: term.pm || 0 });
  }
  return batch;
}

/**
 * Merge theta info: combine same-idx entries, keep different-idx separate.
 * Returns { _tbatch: [{idx, sp, cp}] }
 */
function _mergeThetaPair(ta, tb) {
  var ba = _collectTheta(ta);
  var bb = _collectTheta(tb);
  // Merge matching entries
  for (var j = 0; j < bb.length; j++) {
    var found = false;
    for (var i = 0; i < ba.length; i++) {
      if (ba[i].idx === bb[j].idx) {
        ba[i].sp += bb[j].sp;
        ba[i].cp += bb[j].cp;
        found = true;
        break;
      }
    }
    if (!found) ba.push(bb[j]);
  }
  return { _tbatch: ba };
}

/**
 * Merge phi info: combine matching idx+pfunc entries.
 * Returns { _pbatch: [{idx, pf, pm}] }
 */
function _mergePhiPair(ta, tb) {
  var pa = _collectPhi(ta);
  var pb = _collectPhi(tb);
  for (var j = 0; j < pb.length; j++) {
    var found = false;
    for (var i = 0; i < pa.length; i++) {
      if (pa[i].idx === pb[j].idx && pa[i].pf === pb[j].pf) {
        // Same phi variable and function → could merge multipliers
        found = true;
        break;
      }
    }
    if (!found) pa.push(pb[j]);
  }
  return { _pbatch: pa };
}



// --- Half-angle expansion and grouping ---

/**
 * Expand all theta half-angle factors in a structured term batch.
 * Returns: [{c, thetaBasis: [{idx, func, k}], phiBasis: [{idx, func}], im}]
 * where func is "cos"/"sin"/"1" and k is multiplier of θ/2.
 */
function _expandStructTerm(term) {
  var tbatch = _collectTheta(term);
  var pbatch = _collectPhi(term);
  
  var expanded = [{ s: term.s, thetaBasis: [], phiBasis: pbatch, im: term.im }];
  
  for (var j = 0; j < tbatch.length; j++) {
    var tt = tbatch[j];
    var halfExp = expandHalfAngleBasis(tt.sp, tt.cp);
    if (halfExp.length === 0) continue;
    
    var newList = [];
    for (var e = 0; e < expanded.length; e++) {
      for (var h = 0; h < halfExp.length; h++) {
        var he = halfExp[h];
        if (he.s === '0') continue;
        var base = expanded[e];
        var nb = base.thetaBasis.slice();
        if (he.func !== '1') {
          nb.push({ idx: tt.idx, func: he.func, k: he.k });
        }
        // String coefficient multiplication
        var newS;
        if (he.s === '1') newS = base.s;
        else if (he.s === '-1') newS = '-(' + base.s + ')';
        else newS = '(' + base.s + ')*(' + he.s + ')';
        
        newList.push({ s: newS, thetaBasis: nb, phiBasis: base.phiBasis, im: base.im });
      }
    }
    expanded = newList;
  }
  
  return expanded;
}

function _groupExpandedTerms(expandedList) {
  var groups = {};
  for (var i = 0; i < expandedList.length; i++) {
    var t = expandedList[i];
    var key = _basisKey(t.thetaBasis) + '|' + _basisKey(t.phiBasis) + '|' + (t.im ? 'I' : 'R');
    if (!groups[key]) {
      groups[key] = { coeffStrs: [], thetaBasis: t.thetaBasis, phiBasis: t.phiBasis, im: t.im };
    }
    groups[key].coeffStrs.push(t.s);
  }
  var result = [];
  for (var k in groups) {
    result.push(groups[k]);
  }
  return result;
}

function _basisKey(basis) {
  var sorted = basis.slice().sort(function(a, b) {
    return a.idx - b.idx || ((a.func||a.pf||'') + '').localeCompare((b.func||b.pf||'') + '');
  });
  return sorted.map(function(x) {
    var f = x.func || x.pf || '';
    return x.idx + f + (x.k || '') + (x.pm || '');
  }).join(',');
}

// --- TOP-LEVEL: structured formula computation ---

/**
 * Compute simplified angular formula using structured terms.
 * Returns same format as getAngleFormula but with trig-simplified {real, imag}.
 */
function getAngleFormulaSimplified(decayTree) {
  if (typeof Algebrite === 'undefined') {
    return { error: 'Algebrite library is not loaded.' };
  }
  if (!decayTree.children) {
    return { error: 'Decay tree must have at least one decay.' };
  }
  
  var vertices = getVertices(decayTree);
  var nDecays = vertices.length;
  
  // Build structured terms for each vertex
  var structs = [];
  for (var i = 0; i < vertices.length; i++) {
    var sv = _buildVertexStructTerms(vertices[i], i, i);
    structs.push(sv);
  }
  
  // Cascade combination
  var rawStructured = _combineStructured(structs, decayTree, 0);
  
  // For each (helicity, LS) combination: expand half-angles, group, reconstruct
  var result = {};
  // Track phi and theta names
  var phiNames = [], thetaNames = [];
  for (var i = 0; i < nDecays; i++) {
    phiNames[i] = '\\phi_' + i;
    thetaNames[i] = '\\theta_' + i;
  }
  
  var algebriteThetaNames = [];
  for (var i = 0; i < nDecays; i++) {
    algebriteThetaNames[i] = 'theta_' + i;
  }
  
  for (var fullH in rawStructured) {
    var lsDict = rawStructured[fullH];
    result[fullH] = {};
    
    for (var fullLs in lsDict) {
      var terms = lsDict[fullLs];
      
      // Expand half-angles for each term
      var allExpanded = [];
      for (var t = 0; t < terms.length; t++) {
        var expanded = _expandStructTerm(terms[t]);
        allExpanded = allExpanded.concat(expanded);
      }
      
      // Group by basis
      var grouped = _groupExpandedTerms(allExpanded);
      
      // Reconstruct simplified expressions using Algebrite
      try {
        var realStrs = [], imagStrs = [];
        
        for (var g = 0; g < grouped.length; g++) {
          var grp = grouped[g];
          
          // Build coefficient expression from summed strings
          var coeffRaw;
          if (grp.coeffStrs.length === 1) {
            coeffRaw = grp.coeffStrs[0];
          } else {
            coeffRaw = '(' + grp.coeffStrs.join(')+(') + ')';
          }
          // Simplify the coefficient
          var coeffExpr = Algebrite.run('simplify(' + coeffRaw + ')').trim();
          if (coeffExpr === '0') continue;
          
          // Build trig factor product
          var trigParts = [];
          
          for (var ti = 0; ti < grp.thetaBasis.length; ti++) {
            var tb = grp.thetaBasis[ti];
            if (tb.func === '1' || tb.k === 0) continue;
            var tn = algebriteThetaNames[tb.idx];
            var arg;
            if (tb.k === 1) { arg = tn + '/2'; }
            else if (tb.k % 2 === 0) { var m = tb.k/2; arg = (m===1?'':m+'*') + tn; }
            else { arg = '(' + tb.k + '/2)*' + tn; }
            trigParts.push(tb.func + '(' + arg + ')');
          }
          
          for (var pi = 0; pi < grp.phiBasis.length; pi++) {
            var pb = grp.phiBasis[pi];
            var pn = 'phi_' + pb.idx;
            if (pb.pm && pb.pm !== 1) { trigParts.push(pb.pf + '(' + pb.pm + '*' + pn + ')'); }
            else { trigParts.push(pb.pf + '(' + pn + ')'); }
          }
          
          // Build complete term (coefficient already includes all numeric factors)
          // For imag: exclude *i from individual terms, factor it out later
          var termExpr = coeffExpr;
          if (trigParts.length > 0) {
            termExpr = termExpr + '*' + trigParts.join('*');
          }
          
          if (grp.im) imagStrs.push(termExpr);
          else realStrs.push(termExpr);
        }
        
        // Build Algebrite expressions, factoring i out of imag
        var realAlg = realStrs.length === 0 ? '0' : realStrs.join('+').replace(/\+\+/g,'+').replace(/\+\-/g,'-');
        if (realAlg[0] === '+') realAlg = realAlg.substring(1);
        
        var imagAlg;
        if (imagStrs.length === 0) {
          imagAlg = '0';
        } else {
          // No i prefix — imag is the real-valued coefficient of i
          imagAlg = imagStrs.join('+').replace(/\+\+/g,'+').replace(/\+\-/g,'-');
          if (imagAlg[0] === '+') imagAlg = imagAlg.substring(1);
        }
        
        var expanded = Algebrite.run('expand(' + realAlg + ')').trim();
        var realSimplified = Algebrite.run('simplify(' + expanded + ')').trim();
        var expImag = Algebrite.run('expand(' + imagAlg + ')').trim();
        var imagSimplified = Algebrite.run('simplify(' + expImag + ')').trim();
        
        result[fullH][fullLs] = {
          real: realSimplified,
          imag: imagSimplified
        };
      } catch (e) {
        result[fullH][fullLs] = {
          real: '0',
          imag: '0',
          _raw: grouped
        };
      }
    }
  }
  
  // Build angle list
  var angles = [];
  for (var i = 0; i < nDecays; i++) {
    angles.push('phi_' + i + ', theta_' + i);
  }
  
  return {
    helicities: result,
    nDecays: nDecays,
    angles: angles
  };
}


// Export
if (typeof module !== 'undefined' && module.exports) {
  module.exports = { 
    getAngleFormula: getAngleFormula, 
    getAngleFormulaSimplified: getAngleFormulaSimplified,
    resultToLatex: resultToLatex, 
    toSpin: toSpin,
    complex: complex,
    complexAdd: complexAdd,
    complexMul: complexMul
  };
}
