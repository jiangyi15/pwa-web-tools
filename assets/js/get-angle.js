/**
 * Get Angle Formula - Helicity Amplitude Angular Formula Calculator
 * 
 * Uses a structured term pipeline built from CG coefficients
 * (via computeCGExact) and Wigner-d half-angle Fourier expansions,
 * combining coefficients via Surd/SurdSum exact arithmetic.
 * No CAS dependency — all arithmetic is pure BigInt + integer gcd.
 */



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
  
  var cg1 = computeCGExact(formatSpin(toSpin(Jb)), formatSpin(toSpin(lb)),
                           formatSpin(toSpin(Jc)), formatSpin(toSpin(-lc)),
                           formatSpin(toSpin(s)), formatSpin(toSpin(delta)));
  if (cg1.isZero()) return null;
  
  var cg2 = computeCGExact(String(l), '0',
                           formatSpin(toSpin(s)), formatSpin(toSpin(delta)),
                           formatSpin(toSpin(Ja)), formatSpin(toSpin(delta)));
  if (cg2.isZero()) return null;
  
  // Coefficient: sqrt((2l+1)/(2Ja+1)) * cg1 * cg2  as exact Surd
  var sqrtPart = Surd.div(Surd.fromRadicand(2 * l + 1),
                          Surd.fromRadicand(Math.round(2 * Ja) + 1));
  var coeff = Surd.mul(Surd.mul(sqrtPart, cg1), cg2);
  
  // Get exact Wigner-d weight strings
  var wdWeights = _getExactWignerDWeights(Ja, la, delta);
  if (wdWeights.length === 0) return null;
  
  var terms = [];
  for (var i = 0; i < wdWeights.length; i++) {
    var wt = wdWeights[i];
    // Pre-combine coefficient with Wigner-d weight as exact Surd
    var termSurd = Surd.mul(coeff, Surd.parse(wt.weightStr));
    var termCoeffStr = termSurd.toString();
    
    if (Math.abs(la) < 1e-10) {
      terms.push(_st(termCoeffStr, thetaIdx, wt.sinPow, wt.cosPow, phiIdx, '1', 0));
    } else {
      var absLA = Math.abs(la);
      var sinSign = la < 0 ? -1 : 1;
      terms.push(_st(termCoeffStr, thetaIdx, wt.sinPow, wt.cosPow, phiIdx, 'cos', absLA));
      // For sin: multiply coefficient by sinSign
      var sinSurd = Surd.scale(termSurd, sinSign);
      var sinStr = sinSurd.toString();
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
  
  // Multiply coefficients as exact Surds (both are clean Surd strings)
  var sA = Surd.parse(ta.s);
  var sB = Surd.parse(tb.s);
  
  var result = [];
  // (ar + i*ai)*(br + i*bi) = (ar*br - ai*bi) + i*(ar*bi + ai*br)
  
  if (!ta.im && !tb.im) {
    var s = Surd.mul(sA, sB);
    result.push({ s: s.toString(), _tbatch: theta._tbatch, _pbatch: phi._pbatch, im: false });
  }
  if (ta.im && tb.im) {
    // -ai*bi
    var s = Surd.scale(Surd.mul(sA, sB), -1);
    result.push({ s: s.toString(), _tbatch: theta._tbatch, _pbatch: phi._pbatch, im: false });
  }
  if (!ta.im && tb.im) {
    var s = Surd.mul(sA, sB);
    result.push({ s: s.toString(), _tbatch: theta._tbatch, _pbatch: phi._pbatch, im: true });
  }
  if (ta.im && !tb.im) {
    var s = Surd.mul(sA, sB);
    result.push({ s: s.toString(), _tbatch: theta._tbatch, _pbatch: phi._pbatch, im: true });
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
        // Multiply coefficient as exact Surd (he.s is a rational string like "1/4")
        var newS;
        if (he.s === '1') newS = base.s;
        else if (he.s === '-1') newS = Surd.scale(Surd.parse(base.s), -1).toString();
        else newS = Surd.mul(Surd.parse(base.s), Surd.parse(he.s)).toString();
        
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

function _pmToLatex(pm) {
  // pm is a phi multiplier: integer (1,2,3…) or half-integer (0.5, 1.5, …)
  if (Number.isInteger(pm)) return String(pm);
  var n = Math.round(pm * 2);
  return '\\frac{' + n + '}{2}';
}

// --- TOP-LEVEL: structured formula computation ---

/**
 * Compute simplified angular formula using structured terms.
 * Returns same format as getAngleFormula but with trig-simplified {real, imag}.
 */
function getAngleFormulaSimplified(decayTree) {
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
  
  // Detect J=0 root decay: both child phi terms share the same exponent,
  // so we can fix φ₁=0 and rename φ₂→χ (one fewer azimuthal variable).
  var phiCombine = null;
  if (nDecays >= 3 &&
      toSpin(vertices[0].Ja).value === 0 &&
      decayTree.children && decayTree.children[0] && decayTree.children[1] &&
      decayTree.children[0].children && decayTree.children[1].children) {
    var c0size = countDecayVertices(decayTree.children[0]);
    phiCombine = { removeIdx: 1, renameIdx: 1 + c0size };
  }
  
  // LaTeX names for phi and theta (built directly from structure)
  var phiLatexNames = [], thetaLatexNames = [];
  for (var i = 0; i < nDecays; i++) {
    if (phiCombine && i === phiCombine.removeIdx) {
      phiLatexNames[i] = '0';           // fixed to zero (not used after filtering)
    } else if (phiCombine && i === phiCombine.renameIdx) {
      // Wrapping parens so pm multiplier works on the whole sum, e.g. cos(2(φ₁+φ₂))
      phiLatexNames[i] = '(\\phi_{' + phiCombine.removeIdx + '}+\\phi_{' + phiCombine.renameIdx + '})';
    } else {
      phiLatexNames[i] = '\\phi_{' + i + '}';
    }
    thetaLatexNames[i] = '\\theta_{' + i + '}';
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
      
      // Apply J=0 phi substitution: φ_removeIdx = 0
      //   sin(l·0) = 0 → term vanishes
      //   cos(l·0) = 1 → phi factor drops
      if (phiCombine) {
        var filtered = [];
        for (var gi = 0; gi < grouped.length; gi++) {
          var grp = grouped[gi];
          var newPhi = [], skip = false;
          for (var pi = 0; pi < grp.phiBasis.length; pi++) {
            var pb = grp.phiBasis[pi];
            if (pb.idx === phiCombine.removeIdx) {
              if (pb.pf === 'sin') { skip = true; break; }
            } else {
              newPhi.push(pb);
            }
          }
          if (!skip) {
            filtered.push({
              coeffStrs: grp.coeffStrs,
              thetaBasis: grp.thetaBasis,
              phiBasis: newPhi,
              im: grp.im
            });
          }
        }
        grouped = filtered;
      }
      
      // Reconstruct expressions using SurdSum (no Algebrite)
      try {
        var realStrs = [], imagStrs = [];
        
        for (var g = 0; g < grouped.length; g++) {
          var grp = grouped[g];
          
          // Combine coefficients with SurdSum
          var sSum = new SurdSum();
          for (var ci = 0; ci < grp.coeffStrs.length; ci++) {
            sSum.add(Surd.parse(grp.coeffStrs[ci]));
          }
          if (sSum.isZero()) continue;
          // Use direct LaTeX from SurdSum internal representation
          var coeffStr = sSum.toLatex();
          // Wrap multi-group sums in parens so the trig product applies to all
          if (coeffStr.indexOf(' + ') !== -1 || coeffStr.indexOf(' - ') !== -1) {
            coeffStr = '(' + coeffStr + ')';
          }
          
          // Build trig factor product
          var trigParts = [];
          
          for (var ti = 0; ti < grp.thetaBasis.length; ti++) {
            var tb = grp.thetaBasis[ti];
            if (tb.func === '1' || tb.k === 0) continue;
            var tn = thetaLatexNames[tb.idx];
            var arg;
            if (tb.k === 1) {
              arg = '\\frac{' + tn + '}{2}';
            } else if (tb.k % 2 === 0) {
              var m = tb.k / 2;
              arg = (m === 1 ? '' : String(m)) + tn;
            } else {
              arg = '\\frac{' + tb.k + '}{2}' + tn;
            }
            trigParts.push('\\' + tb.func + '(' + arg + ')');
          }
          
          for (var pi = 0; pi < grp.phiBasis.length; pi++) {
            var pb = grp.phiBasis[pi];
            var pn = phiLatexNames[pb.idx];
            if (pb.pm && pb.pm !== 1) {
              trigParts.push('\\' + pb.pf + '(' + _pmToLatex(pb.pm) + pn + ')');
            } else {
              // Strip wrapping parens for pm=1 (combined J=0 name needs them for pm≠1)
              var name = pn;
              if (name.charAt(0) === '(' && name.charAt(name.length - 1) === ')') {
                name = name.substring(1, name.length - 1);
              }
              trigParts.push('\\' + pb.pf + '(' + name + ')');
            }
          }
          
          // Build complete term (all parts already LaTeX)
          var termExpr = coeffStr;
          if (trigParts.length > 0) {
            termExpr = termExpr + '\\!\\cdot\\!' + trigParts.join('\\!\\cdot\\!');
          }
          
          if (grp.im) imagStrs.push(termExpr);
          else realStrs.push(termExpr);
        }
        
        // Each group has a unique trig basis (from _groupExpandedTerms),
        // so no like-term combining needed across groups. Just join.
        var realExpr = realStrs.length === 0 ? '0' : realStrs.join('+').replace(/\+\-/g,'-');
        if (realExpr[0] === '+') realExpr = realExpr.substring(1);
        
        var imagExpr;
        if (imagStrs.length === 0) {
          imagExpr = '0';
        } else {
          imagExpr = imagStrs.join('+').replace(/\+\-/g,'-');
          if (imagExpr[0] === '+') imagExpr = imagExpr.substring(1);
        }
        
        result[fullH][fullLs] = {
          real: realExpr,
          imag: imagExpr
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
  
  // Build angle list (LaTeX names)
  var angles = [];
  for (var i = 0; i < nDecays; i++) {
    if (phiCombine && i === phiCombine.removeIdx) {
      angles.push('\\theta_' + i);             // only theta (phi fixed to 0)
    } else if (phiCombine && i === phiCombine.renameIdx) {
      angles.push('(\\phi_{' + phiCombine.removeIdx + '}+\\phi_{' + phiCombine.renameIdx + '}), \\theta_' + i);
    } else {
      angles.push('\\phi_' + i + ', \\theta_' + i);
    }
  }
  
  return {
    helicities: result,
    nDecays: nDecays,
    angles: angles,
  };
}


// Export
if (typeof module !== 'undefined' && module.exports) {
  module.exports = { 
    getAngleFormulaSimplified: getAngleFormulaSimplified,
    toSpin: toSpin
  };
}
