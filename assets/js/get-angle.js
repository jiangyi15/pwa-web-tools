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
 */
function complexSimplify(c) {
  try {
    return {
      real: Algebrite.run('simplify(' + c.real + ')').trim(),
      imag: Algebrite.run('simplify(' + c.imag + ')').trim()
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
  var r = c.real;
  var i = c.imag;
  
  try {
    r = Algebrite.run('printlatex(' + r + ')').trim();
    r = r.replace(/(?<!\\)cos/g, '\\cos').replace(/(?<!\\)sin/g, '\\sin');
    r = r.replace(/(?<!\\)phi/g, '\\phi').replace(/(?<!\\)theta/g, '\\theta');
  } catch(e) {}
  
  // Check if imag is 0
  if (i === '0' || i === '') {
    return r;
  }
  
  try {
    i = Algebrite.run('printlatex(' + i + ')').trim();
    i = i.replace(/(?<!\\)cos/g, '\\cos').replace(/(?<!\\)sin/g, '\\sin');
    i = i.replace(/(?<!\\)phi/g, '\\phi').replace(/(?<!\\)theta/g, '\\theta');
  } catch(e) {}
  
  if (i === '0') return r;
  if (r === '0' || r === '') return i + '\\,i';
  return '(' + r + ' + ' + i + '\\,i)';
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
    return toSpin(parseFloat(val));
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
  var spin2 = Math.abs(s.num);
  var values = [];
  for (var h2 = -spin2; h2 <= spin2; h2 += 2) {
    values.push({ num: h2, den: 2, value: h2 / 2 });
  }
  return values;
}

function getLSCombinations(Ja, Jb, Jc) {
  var a = toSpin(Ja);
  var b = toSpin(Jb);
  var c = toSpin(Jc);
  
  var Ja2 = Math.abs(a.num);
  var Jb2 = Math.abs(b.num);
  var Jc2 = Math.abs(c.num);
  
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
                
                if (!(fullH in result)) result[fullH] = {};
                result[fullH][fullLs] = combined;
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
            
            if (!(fullH in result)) result[fullH] = {};
            result[fullH][fullLs] = combined;
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
            
            if (!(fullH in result)) result[fullH] = {};
            result[fullH][fullLs] = combined;
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

// Export
if (typeof module !== 'undefined' && module.exports) {
  module.exports = { 
    getAngleFormula: getAngleFormula, 
    resultToLatex: resultToLatex, 
    toSpin: toSpin,
    complex: complex,
    complexAdd: complexAdd,
    complexMul: complexMul
  };
}
