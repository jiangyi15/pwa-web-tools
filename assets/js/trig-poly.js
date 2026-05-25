/**
 * TrigPoly — trigonometric polynomial in half-angle Fourier basis.
 *
 * A trigonometric polynomial in variables identified by name strings
 * ("theta_0", "phi_1", "chi", etc.) expanded in the canonical basis:
 *   cos(k·name/2)   and   sin(k·name/2)    (k ∈ ℤ)
 *
 * Construction starts in "power form" (sin^sp·cos^cp for theta, raw
 * cos/sin factors for phi) to keep cascade multiplication compact.
 * Call .expand() to convert everything to the Fourier basis and apply
 * product-to-sum identities — the result has at most one factor per
 * variable name per term.
 *
 * Dependencies: Surd, SurdSum (surd.js), expandHalfAngleBasis (wigner-d.js)
 */

// ════════════════════════════════════════════════════════════════════
// TrigPoly constructor
// ════════════════════════════════════════════════════════════════════

function TrigPoly() {
  // Power form terms (before .expand()):
  //   { coeff: Surd, im: boolean,
  //     theta: [{name, sp, cp}],    // sin^sp·cos^cp per theta variable
  //     phi:   [{name, func, k}] }  // func(k·name), func ∈ {cos, sin}
  this._powerTerms = [];

  // Fourier form (after .expand()):
  //   null  —or—
  //   { key -> SurdSum }  where key is a sorted concatenation of
  //   "name:func|k" pairs (each name appears at most once),
  //   or "1" for the constant term.
  this._fourier = null;
}

// ════════════════════════════════════════════════════════════════════
// TERM CREATION
// ════════════════════════════════════════════════════════════════════

/**
 * Add a term in power form.
 *
 * @param {Surd} coeff       Exact coefficient
 * @param {boolean} im       True if this term is the imaginary part
 * @param {Array} thetaEntries  [{name, sp, cp}, ...] — sin^sp·cos^cp per variable
 * @param {Array} phiEntries    [{name, func, k}, ...] — func(k·name)
 * @return {TrigPoly} this
 */
TrigPoly.prototype.addPowerTerm = function (coeff, im, thetaEntries, phiEntries) {
  if (coeff.isZero()) return this;
  this._fourier = null;

  // Merge theta entries by name (same-name power exponents add)
  var thetaMap = {};
  for (var i = 0; i < thetaEntries.length; i++) {
    var t = thetaEntries[i];
    var key = t.name;
    if (thetaMap[key]) {
      thetaMap[key].sp += t.sp || 0;
      thetaMap[key].cp += t.cp || 0;
    } else {
      thetaMap[key] = { name: key, sp: t.sp || 0, cp: t.cp || 0 };
    }
  }
  // Remove entries with sp=cp=0
  var thetaList = [];
  for (var k in thetaMap) {
    var th = thetaMap[k];
    if (th.sp !== 0 || th.cp !== 0) thetaList.push(th);
  }

  this._powerTerms.push({
    coeff: coeff,
    im: !!im,
    theta: thetaList,
    phi: phiEntries ? phiEntries.slice() : []
  });
  return this;
};

/**
 * Add a single Fourier-factor term directly (use after expand).
 * @private
 */
TrigPoly.prototype._addFourierTerm = function (coeff, factors) {
  // factors = [{name, func, k}, ...] — each name at most once, sorted
  if (coeff.isZero()) return;
  if (!this._fourier) this._fourier = {};
  var key = _tpMakeKey(factors);
  _tpAddToMap(this._fourier, key, coeff);
};

// ════════════════════════════════════════════════════════════════════
// ADDITION
// ════════════════════════════════════════════════════════════════════

/**
 * Add another TrigPoly into this one, returning `this` for chaining.
 * Both operands must be in the same form (both power or both expanded).
 */
TrigPoly.prototype.add = function (other) {
  if (this._fourier !== null && other._fourier !== null) {
    // Both expanded — merge Fourier maps
    for (var key in other._fourier) {
      var sum = other._fourier[key];
      if (!sum.isEmpty()) {
        if (!this._fourier[key]) this._fourier[key] = new SurdSum();
        var terms = sum.terms();
        for (var i = 0; i < terms.length; i++) {
          this._fourier[key].add(terms[i]);
        }
      }
    }
  } else if (this._fourier === null && other._fourier === null) {
    // Both in power form — merge term lists
    for (var i = 0; i < other._powerTerms.length; i++) {
      var t = other._powerTerms[i];
      this._powerTerms.push({
        coeff: t.coeff,
        im: t.im,
        theta: t.theta.slice(),
        phi: t.phi.slice()
      });
    }
  } else {
    throw new Error('TrigPoly.add: cannot mix expanded and power-form operands');
  }
  return this;
};

// ════════════════════════════════════════════════════════════════════
// MULTIPLICATION
// ════════════════════════════════════════════════════════════════════

/**
 * Multiply two TrigPolys — returns a new TrigPoly.
 *
 * When multiplying cascade amplitudes (power form):
 *   Θ factors merge (sp/cp add), Φ factors concatenate,
 *   coeffs multiply, im = XOR.
 *
 * When both are expanded (Fourier form):
 *   Convolve the two Fourier maps with product-to-sum identities.
 *
 * @param {TrigPoly} other
 * @param {Object} [opts]
 * @param {boolean} [opts.conjugateSecond]  Conjugate `other` for |T|²
 * @return {TrigPoly}
 */
TrigPoly.prototype.mul = function (other, opts) {
  opts = opts || {};
  var conjSecond = !!opts.conjugateSecond;
  var result = new TrigPoly();

  if (this._fourier !== null && other._fourier !== null) {
    // Both expanded — convolve Fourier maps
    result._fourier = {};
    for (var k1 in this._fourier) {
      var v1 = this._fourier[k1];
      if (v1.isEmpty()) continue;
      for (var k2 in other._fourier) {
        var v2 = other._fourier[k2];
        if (v2.isEmpty()) continue;
        _tpConvolveKeys(result._fourier, k1, v1, k2, v2, conjSecond);
      }
    }
    return result;
  }

  // Power-form multiplication
  for (var i = 0; i < this._powerTerms.length; i++) {
    for (var j = 0; j < other._powerTerms.length; j++) {
      var ta = this._powerTerms[i];
      var tb = other._powerTerms[j];

      // Complex multiplication: (ar + i·ai) × (br + i·bi)
      //   real = ar·br - ai·bi
      //   imag = ar·bi + ai·br
      // With conjugateSecond: (ar + i·ai) × (br - i·bi)
      //   real = ar·br + ai·bi
      //   imag = ai·br - ar·bi

      var ar = ta.coeff, ai = ta.coeff;
      var br = tb.coeff, bi = tb.coeff;

      if (conjSecond) {
        // T1 × conj(T2) = (re1 + i·im1) × (re2 - i·im2)
        //   real = re1·re2 + im1·im2
        //   imag = im1·re2 - re1·im2

        if (!ta.im && !tb.im) {
          // real = re1·re2
          var coeff = Surd.mul(ar, br);
          if (!coeff.isZero()) result._addPowerTermFromMul(coeff, false, ta, tb);
        }
        if (ta.im && tb.im) {
          // real = (+1)·im1·im2  (not -im1·im2 because conj flips sign on im2)
          var coeff = Surd.mul(ai, bi);
          if (!coeff.isZero()) result._addPowerTermFromMul(coeff, false, ta, tb);
        }
        if (!ta.im && tb.im) {
          // imag = -re1·im2
          var coeff = Surd.scale(Surd.mul(ar, bi), -1);
          if (!coeff.isZero()) result._addPowerTermFromMul(coeff, true, ta, tb);
        }
        if (ta.im && !tb.im) {
          // imag = im1·re2
          var coeff = Surd.mul(ai, br);
          if (!coeff.isZero()) result._addPowerTermFromMul(coeff, true, ta, tb);
        }
      } else {
        // Standard: real = ar·br - ai·bi, imag = ar·bi + ai·br

        if (!ta.im && !tb.im) {
          // real(ar·br)
          var coeff = Surd.mul(ar, br);
          if (!coeff.isZero()) result._addPowerTermFromMul(coeff, false, ta, tb, 1, 1);
        }
        if (ta.im && tb.im) {
          // real(-ai·bi)
          var coeff = Surd.scale(Surd.mul(ai, bi), -1);
          if (!coeff.isZero()) result._addPowerTermFromMul(coeff, false, ta, tb, 1, 1);
        }
        if (!ta.im && tb.im) {
          // imag(ar·bi)
          var coeff = Surd.mul(ar, bi);
          if (!coeff.isZero()) result._addPowerTermFromMul(coeff, true, ta, tb, 1, 1);
        }
        if (ta.im && !tb.im) {
          // imag(ai·br)
          var coeff = Surd.mul(ai, br);
          if (!coeff.isZero()) result._addPowerTermFromMul(coeff, true, ta, tb, 1, 1);
        }
      }
    }
  }

  return result;
};

/**
 * Internal helper: build a power term from merging two input terms.
 */
TrigPoly.prototype._addPowerTermFromMul = function (coeff, im, ta, tb) {
  // Merge theta: same-name sp/cp add
  var thetaMap = {};
  for (var i = 0; i < ta.theta.length; i++) {
    var t = ta.theta[i];
    thetaMap[t.name] = { name: t.name, sp: t.sp, cp: t.cp };
  }
  for (var i = 0; i < tb.theta.length; i++) {
    var t = tb.theta[i];
    if (thetaMap[t.name]) {
      thetaMap[t.name].sp += t.sp;
      thetaMap[t.name].cp += t.cp;
    } else {
      thetaMap[t.name] = { name: t.name, sp: t.sp, cp: t.cp };
    }
  }
  var thetaList = [];
  for (var k in thetaMap) {
    var th = thetaMap[k];
    if (th.sp !== 0 || th.cp !== 0) thetaList.push(th);
  }

  // Phi: concatenate
  var phiList = [];
  for (var i = 0; i < ta.phi.length; i++) phiList.push(ta.phi[i]);
  for (var i = 0; i < tb.phi.length; i++) phiList.push(tb.phi[i]);

  this._powerTerms.push({
    coeff: coeff,
    im: im,
    theta: thetaList,
    phi: phiList
  });
};

// ════════════════════════════════════════════════════════════════════
// EXPANSION: power form → Fourier basis
// ════════════════════════════════════════════════════════════════════

/**
 * Expand the polynomial to the Fourier basis.
 *
 * For each power term:
 *   1. Theta factors: sin^sp·cos^cp → Σ c_m · func(m·θ/2) via expandHalfAngleBasis
 *   2. Phi factors:   product-to-sum identity per name-group
 *
 * After expansion, each term has at most one factor per variable name,
 * and all factors are in {cos(k·name/2), sin(k·name/2)} form.
 *
 * @return {TrigPoly} this (for chaining)
 */
TrigPoly.prototype.expand = function () {
  if (this._fourier !== null) return this; // already expanded
  if (this._powerTerms.length === 0) {
    this._fourier = {};
    return this;
  }

  this._fourier = {};

  for (var ti = 0; ti < this._powerTerms.length; ti++) {
    var pt = this._powerTerms[ti];
    this._expandOneTerm(pt);
  }

  // Free power form memory
  this._powerTerms = [];
  return this;
};

/**
 * Expand a single power-form term into the Fourier map.
 */
TrigPoly.prototype._expandOneTerm = function (pt) {
  var coeff = pt.coeff;

  // ── Step 1: Expand theta power factors via expandHalfAngleBasis ──
  // thetaExpansions starts as [{factors: [], surd: Surd.ONE}]
  var thetaExpansions = [{ factors: [], surd: Surd.ONE }];

  for (var ti = 0; ti < pt.theta.length; ti++) {
    var tt = pt.theta[ti];
    var halfExp = expandHalfAngleBasis(tt.sp, tt.cp);
    if (halfExp.length === 0) continue;

    var newList = [];
    for (var ei = 0; ei < thetaExpansions.length; ei++) {
      for (var hi = 0; hi < halfExp.length; hi++) {
        var he = halfExp[hi];
        if (he.s === '0') continue;

        var factors = thetaExpansions[ei].factors.slice();
        if (he.func !== '1') {
          factors.push({ name: tt.name, func: he.func, k: he.k });
        }
        var f = Surd.mul(thetaExpansions[ei].surd, Surd.parse(he.s));
        if (f.isZero()) continue;

        newList.push({ factors: factors, surd: f });
      }
    }
    thetaExpansions = newList;
  }

  // ── Step 2: Expand phi product-to-sum per name-group ──
  // Group phi factors by name
  var phiByName = {};
  for (var i = 0; i < pt.phi.length; i++) {
    var p = pt.phi[i];
    if (!phiByName[p.name]) phiByName[p.name] = [];
    phiByName[p.name].push(p);
  }

  // phiExpansions starts as [{factors: [], surd: Surd.ONE}]
  var phiExpansions = [{ factors: [], surd: Surd.ONE }];

  // Process phi groups in sorted name order for deterministic output
  var phiNameList = Object.keys(phiByName).sort();
  for (var ni = 0; ni < phiNameList.length; ni++) {
    var name = phiNameList[ni];
    var phis = phiByName[name];

    // Start with the first factor's options
    // A single phi factor is already in Fourier form
    // Two+ of the same name need product-to-sum
    var products = [{ pf: phis[0].func, pm: phis[0].k, factorSurd: Surd.ONE }];

    for (var pi = 1; pi < phis.length; pi++) {
      var newProducts = [];
      for (var ppi = 0; ppi < products.length; ppi++) {
        var expanded = _tpPhiProduct(products[ppi].pf, products[ppi].pm,
                                     phis[pi].func, phis[pi].k);
        for (var ei = 0; ei < expanded.length; ei++) {
          var f = Surd.mul(products[ppi].factorSurd,
                           _tpParseFracSurd(expanded[ei].factor));
          if (!f.isZero()) {
            newProducts.push({
              pf: expanded[ei].pf,
              pm: expanded[ei].pm,
              factorSurd: f
            });
          }
        }
      }
      products = newProducts;
    }

    // Cross phiExpansions × products
    var newPhiExp = [];
    for (var ei = 0; ei < phiExpansions.length; ei++) {
      for (var pj = 0; pj < products.length; pj++) {
        var factors = phiExpansions[ei].factors.slice();
        var p = products[pj];
        if (p.pf !== '1') {
          factors.push({ name: name, func: p.pf, k: p.pm });
        }
        var f = Surd.mul(phiExpansions[ei].surd, p.factorSurd);
        if (!f.isZero()) {
          newPhiExp.push({ factors: factors, surd: f });
        }
      }
    }
    phiExpansions = newPhiExp;
  }

  // ── Step 3: Combine theta × phi expansions ──
  for (var tei = 0; tei < thetaExpansions.length; tei++) {
    for (var pei = 0; pei < phiExpansions.length; pei++) {
      var tExp = thetaExpansions[tei];
      var pExp = phiExpansions[pei];
      if (tExp.surd.isZero() || pExp.surd.isZero()) continue;

      var totalCoeff = Surd.mul(coeff, Surd.mul(tExp.surd, pExp.surd));
      if (totalCoeff.isZero()) continue;

      var allFactors = tExp.factors.concat(pExp.factors);
      if (allFactors.length === 0) {
        _tpAddToMap(this._fourier, '1', totalCoeff);
      } else {
        // Remove factors with k=0 (sin(0)=0 kills term, cos(0)=1 drops)
        var filtered = [];
        var kill = false;
        for (var fi = 0; fi < allFactors.length; fi++) {
          var f = allFactors[fi];
          if (f.k === 0) {
            if (f.func === 'sin') { kill = true; break; }
            // cos(0) = 1 — skip
          } else {
            filtered.push(f);
          }
        }
        if (kill) continue;

        filtered.sort(_tpFactorSort);
        var key = _tpMakeKey(filtered);
        _tpAddToMap(this._fourier, key, totalCoeff);
      }
    }
  }

  // If nothing was added (all zero), ensure map is empty but exists
  if (!this._fourier) this._fourier = {};
};

// ════════════════════════════════════════════════════════════════════
// SUBSTITUTION
// ════════════════════════════════════════════════════════════════════

/**
 * Substitute variable names.
 *
 * Map format: { oldName: newNameOrNull }
 *   - If newNameOrNull is a string, rename the variable.
 *   - If newNameOrNull is null/undefined, eliminate the variable
 *     (required cos(0)=1, sin(0)=0 → term killed).
 *
 * This is used for J=0 phiCombine:
 *   substitute({ phi_1: null, phi_2: 'chi' })
 *
 * @param {Object} nameMap  oldName → newName (string) or null
 * @return {TrigPoly} this
 */
TrigPoly.prototype.substitute = function (nameMap) {
  // Ensure expanded
  if (this._fourier === null) this.expand();

  var newFourier = {};

  for (var key in this._fourier) {
    var sum = this._fourier[key];
    if (sum.isEmpty()) continue;

    if (key === '1') {
      // Constant might be affected if a factor with this name is non-constant
      // But '1' has no factors, so it stays unchanged
      newFourier['1'] = sum;
      continue;
    }

    // Parse key into factors, apply substitution
    var factors = _tpParseKey(key);
    var newFactors = [];
    var kill = false;

    for (var fi = 0; fi < factors.length; fi++) {
      var f = factors[fi];
      var newName = nameMap[f.name];

      if (newName === undefined) {
        // Name not in map — keep as-is
        newFactors.push(f);
      } else if (newName === null) {
        // Eliminate: cos(k·0/2)=1 (drop), sin(k·0/2)=0 (kill)
        if (f.func === 'sin') { kill = true; break; }
        // cos: drop the factor (constant value 1)
      } else {
        // Rename
        newFactors.push({ name: newName, func: f.func, k: f.k });
      }
    }

    if (kill) continue;

    newFactors.sort(_tpFactorSort);
    var newKey = _tpMakeKey(newFactors);
    if (!newFourier[newKey]) newFourier[newKey] = new SurdSum();
    var terms = sum.terms();
    for (var ti = 0; ti < terms.length; ti++) {
      newFourier[newKey].add(terms[ti]);
    }
  }

  this._fourier = newFourier;
  return this;
};

// ════════════════════════════════════════════════════════════════════
// OUTPUT
// ════════════════════════════════════════════════════════════════════

/**
 * Return the canonical Fourier map: { key: SurdSum }.
 * Keys are sorted factor lists: "name:func|k,name2:func2|k2".
 * The constant term has key "1".
 *
 * TrigPoly must be expanded before calling this.
 */
TrigPoly.prototype.toFourierMap = function () {
  if (this._fourier === null) this.expand();
  return this._fourier;
};

/**
 * Check if the polynomial is zero.
 */
TrigPoly.prototype.isZero = function () {
  if (this._fourier !== null) {
    for (var k in this._fourier) {
      if (!this._fourier[k].isEmpty()) return false;
    }
    return true;
  }
  return this._powerTerms.length === 0;
};

/**
 * Render the expanded polynomial to a LaTeX string.
 *
 * @param {Object} [opts]
 * @param {Object} [opts.phiCombine]  J=0 combine info for display: {fixIdx, chiIdx}
 * @return {string} LaTeX expression
 */
TrigPoly.prototype.toLatex = function (opts) {
  if (this._fourier === null) this.expand();
  opts = opts || {};

  var keys = Object.keys(this._fourier).sort(_tpSortKey);
  if (keys.length === 0) return '0';

  var parts = [];
  for (var ki = 0; ki < keys.length; ki++) {
    var key = keys[ki];
    var sum = this._fourier[key];
    if (sum.isEmpty()) continue;

    var coeffLatex = sum.toLatex();
    if (coeffLatex === '0') continue;

    var trigLatex = _tpKeyToLatex(key, opts);
    var term = coeffLatex;
    if (trigLatex) {
      term += '\\,' + trigLatex;
    }
    parts.push(term);
  }

  if (parts.length === 0) return '0';
  return parts.join(' + ').replace(/\+ -/g, '- ');
};

// ════════════════════════════════════════════════════════════════════
// PRIVATE HELPERS
// ════════════════════════════════════════════════════════════════════

/**
 * Product-to-sum for two phi factors (same variable).
 * Returns [{pf, pm, factor}] similar to _expandPhiProduct in the old code.
 *
 * cos(a)·cos(b) = ½[cos(a+b) + cos(a-b)]
 * sin(a)·sin(b) = ½[cos(a-b) - cos(a+b)]
 * cos(a)·sin(b) = ½[sin(a+b) - sin(a-b)]
 * sin(a)·cos(b) = ½[sin(a+b) + sin(a-b)]
 */
function _tpPhiProduct(pf1, pm1, pf2, pm2) {
  if (pf1 === '1') return [{ pf: pf2, pm: pm2 || 0, factor: '1' }];
  if (pf2 === '1') return [{ pf: pf1, pm: pm1 || 0, factor: '1' }];

  if (pf1 === 'cos' && pf2 === 'cos') {
    return [
      { pf: 'cos', pm: (pm1 || 0) + (pm2 || 0), factor: '1/2' },
      { pf: 'cos', pm: Math.abs((pm1 || 0) - (pm2 || 0)), factor: '1/2' }
    ];
  }
  if (pf1 === 'sin' && pf2 === 'sin') {
    return [
      { pf: 'cos', pm: Math.abs((pm1 || 0) - (pm2 || 0)), factor: '1/2' },
      { pf: 'cos', pm: (pm1 || 0) + (pm2 || 0), factor: '-1/2' }
    ];
  }
  if (pf1 === 'cos' && pf2 === 'sin') {
    // cos(a)·sin(b) = ½[sin(a+b) - sin(a-b)]
    // When a < b: sin(a-b) = -sin(b-a) → -½·(-sin(b-a)) = +½·sin(b-a)
    var pmA = pm1 || 0, pmB = pm2 || 0;
    return [
      { pf: 'sin', pm: pmA + pmB, factor: '1/2' },
      { pf: 'sin', pm: Math.abs(pmA - pmB), factor: (pmA >= pmB) ? '-1/2' : '1/2' }
    ];
  }
  if (pf1 === 'sin' && pf2 === 'cos') {
    // sin(a)·cos(b) = ½[sin(a+b) + sin(a-b)]
    // When a < b: sin(a-b) = -sin(b-a) → +½·(-sin(b-a)) = -½·sin(b-a)
    var pmA = pm1 || 0, pmB = pm2 || 0;
    return [
      { pf: 'sin', pm: pmA + pmB, factor: '1/2' },
      { pf: 'sin', pm: Math.abs(pmA - pmB), factor: (pmA >= pmB) ? '1/2' : '-1/2' }
    ];
  }
  return [];
}

/**
 * Convolve two Fourier keys (product-to-sum across all factors).
 * Used when both TrigPolys are in expanded form.
 */
function _tpConvolveKeys(resultMap, k1, v1, k2, v2, conjSecond) {
  // Parse keys into factor lists
  var factors1 = (k1 === '1') ? [] : _tpParseKey(k1);
  var factors2 = (k2 === '1') ? [] : _tpParseKey(k2);

  // Group factors by name
  var groups = {};
  for (var i = 0; i < factors1.length; i++) {
    var f = factors1[i];
    if (!groups[f.name]) groups[f.name] = [];
    groups[f.name].push({ func: f.func, k: f.k, src: 1 });
  }
  for (var i = 0; i < factors2.length; i++) {
    var f = factors2[i];
    if (!groups[f.name]) groups[f.name] = [];
    groups[f.name].push({ func: f.func, k: f.k, src: 2 });
  }

  var nameList = Object.keys(groups).sort();

  // Cross product of v1 × v2 SurdSums
  var v1terms = v1.terms();
  var v2terms = v2.terms();

  for (var vi = 0; vi < v1terms.length; vi++) {
    for (var vj = 0; vj < v2terms.length; vj++) {
      var coeff = Surd.mul(v1terms[vi], v2terms[vj]);
      if (coeff.isZero()) continue;

      // For each name group, apply product-to-sum
      var expansions = [{ factors: [], surd: Surd.ONE }];

      for (var ni = 0; ni < nameList.length; ni++) {
        var name = nameList[ni];
        var grp = groups[name];

        // Product within this name group
        var prods = [{ pf: grp[0].func, pm: grp[0].k, factorSurd: Surd.ONE }];
        for (var gi = 1; gi < grp.length; gi++) {
          var newProds = [];
          for (var ppi = 0; ppi < prods.length; ppi++) {
            var expanded = _tpPhiProduct(prods[ppi].pf, prods[ppi].pm,
                                         grp[gi].func, grp[gi].k);
            for (var ei = 0; ei < expanded.length; ei++) {
              var f = Surd.mul(prods[ppi].factorSurd,
                               _tpParseFracSurd(expanded[ei].factor));
              if (!f.isZero()) {
                newProds.push({ pf: expanded[ei].pf, pm: expanded[ei].pm, factorSurd: f });
              }
            }
          }
          prods = newProds;
        }

        // Cross with current expansions
        var newExp = [];
        for (var ei = 0; ei < expansions.length; ei++) {
          for (var pj = 0; pj < prods.length; pj++) {
            var factors = expansions[ei].factors.slice();
            if (prods[pj].pf !== '1') {
              factors.push({ name: name, func: prods[pj].pf, k: prods[pj].pm });
            }
            var f = Surd.mul(expansions[ei].surd, prods[pj].factorSurd);
            if (!f.isZero()) {
              newExp.push({ factors: factors, surd: f });
            }
          }
        }
        expansions = newExp;
      }

      // Add all resulting terms to the map
      for (var ei = 0; ei < expansions.length; ei++) {
        var totalCoeff = Surd.mul(coeff, expansions[ei].surd);
        if (totalCoeff.isZero()) continue;

        var factors = expansions[ei].factors;
        // Filter k=0 factors
        var filtered = [], kill = false;
        for (var fi = 0; fi < factors.length; fi++) {
          var f = factors[fi];
          if (f.k === 0) {
            if (f.func === 'sin') { kill = true; break; }
          } else {
            filtered.push(f);
          }
        }
        if (kill) continue;
        filtered.sort(_tpFactorSort);
        var key = _tpMakeKey(filtered);
        _tpAddToMap(resultMap, key, totalCoeff);
      }
    }
  }
}

/**
 * Parse a key string back to factor array.
 * Inverse of _tpMakeKey.
 */
function _tpParseKey(key) {
  if (key === '1') return [];
  var parts = key.split(',');
  var factors = [];
  for (var i = 0; i < parts.length; i++) {
    var m = parts[i].match(/([\w\d_]+):(\w+)\|(-?[\d.]+)/);
    if (m) {
      factors.push({ name: m[1], func: m[2], k: parseFloat(m[3]) });
    }
  }
  return factors;
}

/**
 * Build a canonical key string from sorted factors.
 */
function _tpMakeKey(factors) {
  if (!factors || factors.length === 0) return '1';
  var parts = [];
  for (var i = 0; i < factors.length; i++) {
    var f = factors[i];
    parts.push(f.name + ':' + f.func + '|' + f.k);
  }
  return parts.join(',');
}

/**
 * Sort comparator for factors: by name, then k, then func.
 */
function _tpFactorSort(a, b) {
  if (a.name !== b.name) return a.name < b.name ? -1 : 1;
  if (a.k !== b.k) return a.k - b.k;
  return a.func < b.func ? -1 : (a.func > b.func ? 1 : 0);
}

/**
 * Sort comparator for Fourier keys: by max |k|, then by length.
 */
function _tpSortKey(a, b) {
  if (a === '1') return -1;
  if (b === '1') return 1;
  var maxA = _tpMaxK(a);
  var maxB = _tpMaxK(b);
  if (maxA !== maxB) return maxA - maxB;
  if (a.length !== b.length) return a.length - b.length;
  return a < b ? -1 : (a > b ? 1 : 0);
}

function _tpMaxK(key) {
  var maxK = 0;
  var re = /\|(-?[\d.]+)/g;
  var m;
  while ((m = re.exec(key)) !== null) {
    maxK = Math.max(maxK, Math.abs(parseFloat(m[1])));
  }
  return maxK;
}

/**
 * Add a Surd to a map value (creating SurdSum if needed).
 */
function _tpAddToMap(map, key, surd) {
  if (surd.isZero()) return;
  if (!map[key]) map[key] = new SurdSum();
  map[key].add(surd);
  // Clean up if the result is zero
  if (map[key].isEmpty()) delete map[key];
}

/**
 * Parse a fraction string like "1/2", "-1/2" into a Surd.
 */
function _tpParseFracSurd(s) {
  try { return Surd.parse(s); } catch (e) { return Surd.ONE; }
}

// ════════════════════════════════════════════════════════════════════
// LaTeX RENDERING
// ════════════════════════════════════════════════════════════════════

/**
 * Render a single Fourier key to LaTeX.
 *
 * @param {string} key
 * @param {Object} [opts]
 * @param {Object} [opts.phiCombine]  J=0 combine info: {fixIdx, chiIdx}
 * @return {string} LaTeX trig product, empty for constant
 */
function _tpKeyToLatex(key, opts) {
  if (key === '1' || key === '') return '';
  opts = opts || {};
  var phiCombine = opts.phiCombine;

  var factors = _tpParseKey(key);
  var parts = [];

  for (var i = 0; i < factors.length; i++) {
    var f = factors[i];
    var name = f.name;

    // Determine the LaTeX name
    var latexName;
    if (name.indexOf('theta_') === 0) {
      var idx = name.substring(6);
      latexName = '\\theta_{' + idx + '}';
    } else if (name.indexOf('phi_') === 0) {
      var idx = name.substring(4);
      if (phiCombine && parseInt(idx, 10) < 0) {
        // Negative idx → χ variable → display as φ_fixIdx + φ_chiIdx
        latexName = '\\phi_{' + phiCombine.fixIdx + '}+\\phi_{' + Math.abs(parseInt(idx, 10)) + '}';
        if (f.k !== 1) latexName = '(' + latexName + ')';
      } else {
        latexName = '\\phi_{' + idx + '}';
      }
    } else if (name === 'chi') {
      if (phiCombine) {
        latexName = '\\phi_{' + phiCombine.fixIdx + '}+\\phi_{' + phiCombine.chiIdx + '}';
        if (f.k !== 1) latexName = '(' + latexName + ')';
      } else {
        latexName = '\\chi';
      }
    } else {
      latexName = name;
    }

    // Build argument: func(k·name/2)
    // For phi variables, we use full-angle: func(k·name)
    // For theta variables, we use half-angle: func(k·name/2)
    var isPhi = (name.indexOf('phi_') === 0 || name === 'chi');
    var k = f.k;
    var arg;
    if (isPhi) {
      // Phi: full-angle form cos(k·φ) or sin(k·φ)
      // Handle half-integer k (from half-integer helicities)
      var kInt = Math.round(k * 2);
      if (Math.abs(kInt - k * 2) < 1e-10 && kInt % 2 === 1) {
        arg = '\\frac{' + kInt + '}{2}' + latexName;
        if (kInt === 1) arg = '\\frac{' + latexName + '}{2}';
      } else if (k === 1) {
        arg = latexName;
      } else {
        arg = String(k) + latexName;
      }
    } else {
      // Theta: half-angle form cos(k·θ/2) or sin(k·θ/2)
      if (k === 1) {
        arg = '\\frac{' + latexName + '}{2}';
      } else if (k % 2 === 0) {
        var n = k / 2;
        arg = (n === 1 ? '' : String(n)) + latexName;
      } else {
        arg = '\\frac{' + k + '}{2}' + latexName;
      }
    }
    parts.push('\\' + f.func + '(' + arg + ')');
  }

  return parts.join('\\,');
}

// ════════════════════════════════════════════════════════════════════
// CONVENIENCE CONSTRUCTORS
// ════════════════════════════════════════════════════════════════════

/**
 * Create a TrigPoly from a power-form term (like _st in get-angle.js).
 *
 * @param {Surd} coeff       Exact coefficient
 * @param {string} thetaName   Variable name for theta (e.g. "theta_0")
 * @param {number} sinPow      Exponent of sin(theta/2)
 * @param {number} cosPow      Exponent of cos(theta/2)
 * @param {string} phiName     Variable name for phi (e.g. "phi_1")
 * @param {string} phiFunc     "cos"|"sin"|"1"
 * @param {number} phiK        Multiplier for phi (m in cos(m·φ))
 * @param {boolean} isImag     Is this the imaginary amplitude?
 * @return {TrigPoly}
 */
TrigPoly.fromPowerTerm = function (coeff, thetaName, sinPow, cosPow, phiName, phiFunc, phiK, isImag) {
  var tp = new TrigPoly();
  var thetaEntries = [{ name: thetaName, sp: sinPow || 0, cp: cosPow || 0 }];
  var phiEntries = [];
  if (phiFunc && phiFunc !== '1') {
    phiEntries.push({ name: phiName, func: phiFunc, k: phiK || 0 });
  }
  tp.addPowerTerm(coeff, !!isImag, thetaEntries, phiEntries);
  return tp;
};

/**
 * Deep-clone a TrigPoly.
 */
TrigPoly.prototype.clone = function () {
  var tp = new TrigPoly();
  if (this._fourier !== null) {
    tp._fourier = {};
    for (var key in this._fourier) {
      if (!this._fourier[key].isEmpty()) {
        tp._fourier[key] = new SurdSum();
        var terms = this._fourier[key].terms();
        for (var i = 0; i < terms.length; i++) {
          tp._fourier[key].add(terms[i]);
        }
      }
    }
  } else {
    for (var i = 0; i < this._powerTerms.length; i++) {
      var t = this._powerTerms[i];
      tp._powerTerms.push({
        coeff: t.coeff,
        im: t.im,
        theta: t.theta.slice(),
        phi: t.phi.slice()
      });
    }
  }
  return tp;
};

// ════════════════════════════════════════════════════════════════════
// MODULE EXPORT
// ════════════════════════════════════════════════════════════════════

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { TrigPoly: TrigPoly };
}
