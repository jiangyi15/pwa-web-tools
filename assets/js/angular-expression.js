/**
 * AngularExpression — computation and LaTeX rendering for angular distributions.
 * Pure functions, no DOM dependencies.
 *
 * Uses TrigPoly internally for trigonometric polynomial arithmetic.
 * Depends on: surd.js, cg.js, wigner-d.js, trig-poly.js, get-angle.js
 */
var AngularExpression = (function() {

  // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  // PUBLIC API
  // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  var api = {};

  // ── Main computation: I(θ,φ) = Σ|T|² for a single decay tree ──
  api.compute = function(tree, options) {
    if (!tree.children) return { error: 'Root must decay.' };

    var vertices = getVertices(tree);
    var nDecays = vertices.length;

    // Detect J=0 root decay
    var phiCombine = null;
    if (nDecays >= 3 &&
        toSpin(vertices[0].Ja).value === 0 &&
        tree.children && tree.children[0] && tree.children[1] &&
        tree.children[0].children && tree.children[1].children) {
      var c0size = countDecayVertices(tree.children[0]);
      phiCombine = { fixIdx: 1, chiIdx: 1 + c0size };
    }

    // Build structured vertex terms
    var thetaIdx = 0, phiIdx = 0;
    var vertexTerms = [];
    for (var vi = 0; vi < vertices.length; vi++) {
      vertexTerms.push(_buildVertexStructTerms(vertices[vi], thetaIdx, phiIdx));
      thetaIdx++;
      phiIdx++;
    }

    // Combine across cascade
    var combined = _combineStructured(vertexTerms, tree, 0);

    // Collect all helicity keys
    var allHK = Object.keys(combined);

    // ── Filter helicities if specified ──
    if (options && (options.rootHelicities || options.finalHelicities)) {
      var rootFilter = options.rootHelicities ? _normalizeHelicityList(options.rootHelicities) : null;
      var finalFilter = options.finalHelicities ? _normalizeHelicityList(options.finalHelicities) : null;

      allHK = allHK.filter(function(hk) {
        var parts = hk.split(',').map(function(s) { return parseFloat(s); });

        // Filter root helicity (position 0)
        if (rootFilter && _helIdx(parts[0], rootFilter) < 0) return false;

        // Filter all final/root particle helicities (all positions)
        if (finalFilter) {
          for (var i = 0; i < parts.length; i++) {
            if (_helIdx(parts[i], finalFilter) < 0) return false;
          }
        }

        return true;
      });
    }

    if (allHK.length === 0) return { error: 'No valid helicity combinations.' };

    var allLS = {};
    for (var hi = 0; hi < allHK.length; hi++) {
      var lsKeys = Object.keys(combined[allHK[hi]]);
      for (var li = 0; li < lsKeys.length; li++) {
        allLS[lsKeys[li]] = true;
      }
    }
    var lsKeyList = Object.keys(allLS).sort();

    // For each LS key, accumulate |T|² contributions from all helicities
    var lsMap = {};
    for (var li = 0; li < lsKeyList.length; li++) {
      lsMap[lsKeyList[li]] = { fourier: {} };
    }

    for (var hi = 0; hi < allHK.length; hi++) {
      var hk = allHK[hi];
      var lsDict = combined[hk];

      for (var li = 0; li < lsKeyList.length; li++) {
        var lk = lsKeyList[li];
        var terms = lsDict[lk];
        if (!terms || terms.length === 0) continue;

        // Build TrigPoly from cascade-combined power-form terms
        var tp = new TrigPoly();
        for (var ti = 0; ti < terms.length; ti++) {
          var t = terms[ti];
          var coeff = _tpParseSurd(t.s);
          var thetaEntries = _collectThetaEntries(t);
          var phiEntries = _collectPhiEntries(t);
          tp.addPowerTerm(coeff, !!t.im, thetaEntries, phiEntries);
        }

        // |T|² = T × conj(T)
        var tpSq = tp.mul(tp, { conjugateSecond: true });

        // Expand to Fourier basis
        tpSq.expand();

        // Apply J=0 φ₁→0, φ₂→χ substitution
        if (phiCombine) {
          var subMap = {};
          subMap['phi_' + phiCombine.fixIdx] = null;
          subMap['phi_' + phiCombine.chiIdx] = 'chi';
          tpSq.substitute(subMap);
        }

        lsMap[lk].fourier = tpSq.toFourierMap();
      }
    }

    // ── Interference terms ──
    var showInterference = options && options.showInterference;
    var interfMap = {};
    var interfImagMap = {};

    if (showInterference && lsKeyList.length >= 2) {
      for (var li1 = 0; li1 < lsKeyList.length; li1++) {
        for (var li2 = li1 + 1; li2 < lsKeyList.length; li2++) {
          var key = lsKeyList[li1] + '\u00d7' + lsKeyList[li2];
          interfMap[key] = { fourier: {} };
          interfImagMap[key] = { fourier: {} };
        }
      }

      for (var hi = 0; hi < allHK.length; hi++) {
        var hk = allHK[hi];
        var lsDict = combined[hk];

        for (var li1 = 0; li1 < lsKeyList.length; li1++) {
          var lk1 = lsKeyList[li1];
          var terms1 = lsDict[lk1];
          if (!terms1 || terms1.length === 0) continue;

          for (var li2 = li1 + 1; li2 < lsKeyList.length; li2++) {
            var lk2 = lsKeyList[li2];
            var terms2 = lsDict[lk2];
            if (!terms2 || terms2.length === 0) continue;

            var pairKey = lk1 + '\u00d7' + lk2;

            // Build TrigPolys
            var tp1 = _termsToTrigPoly(terms1);
            var tp2 = _termsToTrigPoly(terms2);

            // 2 Re[T1·T2*]
            var tpRe = tp1.mul(tp2, { conjugateSecond: true });
            // Multiply by 2 (the 2 Re factor)
            tpRe = _tpScale(tpRe, 2);
            tpRe.expand();
            if (phiCombine) {
              var subMap = {};
              subMap['phi_' + phiCombine.fixIdx] = null;
              subMap['phi_' + phiCombine.chiIdx] = 'chi';
              tpRe.substitute(subMap);
            }
            interfMap[pairKey].fourier = tpRe.toFourierMap();

            // 2 Im[T1·T2*] = 2 · (im1·re2 - re1·im2)
            var imagResult = _tpImagCrossMul(terms1, terms2);
            if (!imagResult.isZero()) {
              imagResult = _tpScale(imagResult, 2);
              imagResult.expand();
              _tpApplyPhiCombine(imagResult, phiCombine);
              interfImagMap[pairKey].fourier = imagResult.toFourierMap();
            }
          }
        }
      }
    }

    return { lsMap: lsMap, nDecays: nDecays, phiCombine: phiCombine,
             interfMap: interfMap, interfImagMap: interfImagMap,
             combinedRaw: combined, lsKeyList: lsKeyList };
  };

  // ── Cross-chain interference ──
  api.computeCrossChain = function(combined1, lsKeys1, combined2, lsKeys2, phiCombine) {
    var hk1List = Object.keys(combined1);
    var hk2List = Object.keys(combined2);

    var interfMap = {};
    var interfImagMap = {};
    for (var i1 = 0; i1 < lsKeys1.length; i1++) {
      for (var i2 = 0; i2 < lsKeys2.length; i2++) {
        var key = lsKeys1[i1] + '\u00d7' + lsKeys2[i2];
        interfMap[key] = { fourier: {} };
        interfImagMap[key] = { fourier: {} };
      }
    }

    for (var h1 = 0; h1 < hk1List.length; h1++) {
      var lsDict1 = combined1[hk1List[h1]];
      for (var h2 = 0; h2 < hk2List.length; h2++) {
        var lsDict2 = combined2[hk2List[h2]];

        for (var i1 = 0; i1 < lsKeys1.length; i1++) {
          var lk1 = lsKeys1[i1];
          var terms1 = lsDict1[lk1];
          if (!terms1 || terms1.length === 0) continue;

          for (var i2 = 0; i2 < lsKeys2.length; i2++) {
            var lk2 = lsKeys2[i2];
            var terms2 = lsDict2[lk2];
            if (!terms2 || terms2.length === 0) continue;

            var pairKey = lk1 + '\u00d7' + lk2;

            // Build TrigPolys
            var tp1 = _termsToTrigPoly(terms1);
            var tp2 = _termsToTrigPoly(terms2);

            // 2 Re[T1·T2*]
            var tpRe = tp1.mul(tp2, { conjugateSecond: true });
            tpRe = _tpScale(tpRe, 2);
            tpRe.expand();
            _tpApplyPhiCombine(tpRe, phiCombine);
            interfMap[pairKey].fourier = tpRe.toFourierMap();

            // 2 Im[T1·T2*]
            var imagResult = _tpImagCrossMul(terms1, terms2);
            if (!imagResult.isZero()) {
              imagResult = _tpScale(imagResult, 2);
              imagResult.expand();
              _tpApplyPhiCombine(imagResult, phiCombine);
              interfImagMap[pairKey].fourier = imagResult.toFourierMap();
            }
          }
        }
      }
    }

    return { interfMap: interfMap, interfImagMap: interfImagMap };
  };

  // ── Render a Fourier map to LaTeX string ──
  api.renderFourier = function(fourierMap, phiCombine) {
    var keys = Object.keys(fourierMap).sort(_tpSortKey);

    if (keys.length === 0) return '';

    var lines = '';
    var first = true;
    for (var ki = 0; ki < keys.length; ki++) {
      var key = keys[ki];
      var sum = fourierMap[key];
      if (sum.isEmpty()) continue;

      var sumStr = sum.toLatex();
      if (sumStr === '0') continue;

      var fourierStr = _tpKeyToLatex(key, { phiCombine: phiCombine });

      if (!first) {
        if (sumStr.charAt(0) === '-') {
          lines += ' \\!-\\! ';
          sumStr = sumStr.substring(1);
        } else {
          lines += ' \\!+\\! ';
        }
      }
      lines += sumStr;
      if (fourierStr) lines += '\\,' + fourierStr;
      first = false;
    }

    return lines;
  };

  // ── Render a single Fourier key to LaTeX ──
  api.fourierKeyToLatex = function(key, phiCombine) {
    return _tpKeyToLatex(key, { phiCombine: phiCombine });
  };

  // ── Collect all unique non-constant basis keys ──
  api.collectBasisKeys = function(results) {
    var basis = {};
    for (var ri = 0; ri < results.length; ri++) {
      var r = results[ri];
      if (!r.result) continue;
      var maps = [r.result.lsMap, r.result.interfMap, r.result.interfImagMap];
      for (var mi = 0; mi < maps.length; mi++) {
        var m = maps[mi];
        if (!m) continue;
        var mKeys = Object.keys(m);
        for (var ki = 0; ki < mKeys.length; ki++) {
          var fMap = m[mKeys[ki]].fourier;
          var fKeys = Object.keys(fMap);
          for (var fi = 0; fi < fKeys.length; fi++) {
            var fk = fKeys[fi];
            if (fk === '1') continue;
            if (fMap[fk].isEmpty()) continue;
            basis[fk] = true;
          }
        }
      }
    }
    return Object.keys(basis).sort(function(a, b) {
      if (a.length !== b.length) return a.length - b.length;
      return a.localeCompare(b);
    });
  };

  // ── Find first non-null phiCombine from a set of results ──
  api.findPhiCombine = function(results) {
    for (var ri = 0; ri < results.length; ri++) {
      if (results[ri].result && results[ri].result.phiCombine) return results[ri].result.phiCombine;
    }
    return null;
  };

  // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  // INTERNAL HELPERS
  // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  /**
   * Collect theta entries from a term object.
   * Handles all formats:
   *   - Cascade-combined: { _tbatch: [{idx, sp, cp}] }
   *   - Raw _st:          { ti, sp, cp }
   *   - Old format:        { THETA: [...], theta: [...] }
   * Returns [{name, sp, cp}] where name = "theta_{idx}".
   */
  function _collectThetaEntries(term) {
    var map = {};
    // _tbatch format (from _mulTwo/_mulThree cascade)
    if (term._tbatch) {
      for (var i = 0; i < term._tbatch.length; i++) {
        var e = term._tbatch[i];
        if (e.sp !== 0 || e.cp !== 0) {
          map[e.idx] = { sp: (map[e.idx] ? map[e.idx].sp : 0) + (e.sp || 0),
                         cp: (map[e.idx] ? map[e.idx].cp : 0) + (e.cp || 0) };
        }
      }
    }
    // _st format: ti/sp/cp
    if (term.ti !== undefined) {
      var idx = term.ti;
      map[idx] = { sp: (map[idx] ? map[idx].sp : 0) + (term.sp || 0),
                   cp: (map[idx] ? map[idx].cp : 0) + (term.cp || 0) };
    }
    // THETA format
    if (term.THETA) {
      for (var i = 0; i < term.THETA.length; i++) {
        var e = term.THETA[i];
        if (e.idx !== undefined && (e.sp || e.cp)) {
          map[e.idx] = { sp: (map[e.idx] ? map[e.idx].sp : 0) + (e.sp || 0),
                         cp: (map[e.idx] ? map[e.idx].cp : 0) + (e.cp || 0) };
        }
      }
    }
    // theta format (lowercase)
    if (term.theta) {
      for (var i = 0; i < term.theta.length; i++) {
        var e = term.theta[i];
        if (e.idx !== undefined) {
          map[e.idx] = { sp: (map[e.idx] ? map[e.idx].sp : 0) + (e.sp || 0),
                         cp: (map[e.idx] ? map[e.idx].cp : 0) + (e.cp || 0) };
        }
      }
    }
    // Convert to name-based array
    var result = [];
    for (var idx in map) {
      if (map[idx].sp !== 0 || map[idx].cp !== 0) {
        result.push({ name: 'theta_' + idx, sp: map[idx].sp, cp: map[idx].cp });
      }
    }
    return result;
  }

  /**
   * Collect phi entries from a term object.
   * Handles all formats:
   *   - Cascade-combined: { _pbatch: [{idx, pf, pm}] }
   *   - Raw _st:          { pi, pf, pm }
   *   - Old format:        { PHI: [...], phi: [...] }
   * Returns [{name, func, k}] where name = "phi_{idx}".
   */
  function _collectPhiEntries(term) {
    var entries = [];
    // _pbatch format
    if (term._pbatch) {
      for (var i = 0; i < term._pbatch.length; i++) {
        var e = term._pbatch[i];
        if (e.pf && e.pf !== '1') {
          entries.push({ name: 'phi_' + e.idx, func: e.pf, k: (e.pm || 0) });
        }
      }
    }
    // _st format: pi/pf/pm
    if (term.pi !== undefined && term.pf && term.pf !== '1') {
      entries.push({ name: 'phi_' + term.pi, func: term.pf, k: (term.pm || 0) });
    }
    // PHI format
    if (term.PHI) {
      for (var i = 0; i < term.PHI.length; i++) {
        var e = term.PHI[i];
        if (e.idx !== undefined && e.pf && e.pf !== '1') {
          entries.push({ name: 'phi_' + e.idx, func: e.pf, k: (e.pm || 0) });
        }
      }
    }
    // phi format (lowercase)
    if (term.phi) {
      for (var i = 0; i < term.phi.length; i++) {
        var e = term.phi[i];
        if (e.idx !== undefined && e.pf && e.pf !== '1') {
          entries.push({ name: 'phi_' + e.idx, func: e.pf, k: (e.pm || 0) });
        }
      }
    }
    return entries;
  }

  /**
   * Build a TrigPoly from a list of cascade-combined terms.
   */
  function _termsToTrigPoly(terms) {
    var tp = new TrigPoly();
    for (var i = 0; i < terms.length; i++) {
      var t = terms[i];
      var coeff = _tpParseSurd(t.s);
      var thetaEntries = _collectThetaEntries(t);
      var phiEntries = _collectPhiEntries(t);
      tp.addPowerTerm(coeff, !!t.im, thetaEntries, phiEntries);
    }
    return tp;
  }

  /**
   * Compute the imaginary part of T1 × conj(T2):
   *   Im[T1·conj(T2)] = im1·re2 - re1·im2
   * Returns a TrigPoly (in power form) for Im[T1·T2*].
   * Caller must expand() and substitute().
   */
  function _tpImagCrossMul(terms1, terms2) {
    var tp = new TrigPoly();
    for (var i = 0; i < terms1.length; i++) {
      for (var j = 0; j < terms2.length; j++) {
        var ta = terms1[i], tb = terms2[j];
        var coeffA = _tpParseSurd(ta.s);
        var coeffB = _tpParseSurd(tb.s);
        var thetaA = _collectThetaEntries(ta);
        var phiA = _collectPhiEntries(ta);
        var thetaB = _collectThetaEntries(tb);
        var phiB = _collectPhiEntries(tb);

        // im1·re2: result is imag
        if (ta.im && !tb.im) {
          var coeff = Surd.mul(coeffA, coeffB);
          if (!coeff.isZero()) {
            tp.addPowerTerm(coeff, true, _tpMergeTheta(thetaA, thetaB),
                            phiA.concat(phiB));
          }
        }
        // -re1·im2: result is imag
        if (!ta.im && tb.im) {
          var coeff = Surd.scale(Surd.mul(coeffA, coeffB), -1);
          if (!coeff.isZero()) {
            tp.addPowerTerm(coeff, true, _tpMergeTheta(thetaA, thetaB),
                            phiA.concat(phiB));
          }
        }
      }
    }
    return tp;
  }

  /**
   * Merge two theta entry arrays by name.
   * Returns [{name, sp, cp}].
   */
  function _tpMergeTheta(a, b) {
    var map = {};
    for (var i = 0; i < a.length; i++) {
      map[a[i].name] = { name: a[i].name, sp: a[i].sp, cp: a[i].cp };
    }
    for (var i = 0; i < b.length; i++) {
      if (map[b[i].name]) {
        map[b[i].name].sp += b[i].sp;
        map[b[i].name].cp += b[i].cp;
      } else {
        map[b[i].name] = { name: b[i].name, sp: b[i].sp, cp: b[i].cp };
      }
    }
    var result = [];
    for (var k in map) {
      if (map[k].sp !== 0 || map[k].cp !== 0) result.push(map[k]);
    }
    return result;
  }

  /**
   * Scale all coefficients in a TrigPoly (power form) by a rational factor.
   */
  function _tpScale(tp, num, den) {
    if (den === undefined) den = 1;
    // Work with power form: iterate terms, scale coefficients
    if (tp._fourier !== null) {
      throw new Error('_tpScale: not implemented for expanded form');
    }
    var result = new TrigPoly();
    for (var i = 0; i < tp._powerTerms.length; i++) {
      var pt = tp._powerTerms[i];
      var newCoeff = Surd.scale(pt.coeff, num, den);
      if (!newCoeff.isZero()) {
        result.addPowerTerm(newCoeff, pt.im, pt.theta, pt.phi);
      }
    }
    return result;
  }

  /**
   * Apply phiCombine substitution to a TrigPoly.
   */
  function _tpApplyPhiCombine(tp, phiCombine) {
    if (!phiCombine) return;
    var subMap = {};
    subMap['phi_' + phiCombine.fixIdx] = null;
    subMap['phi_' + phiCombine.chiIdx] = 'chi';
    tp.substitute(subMap);
  }

  /**
   * Parse a coefficient string to Surd, with safe fallback.
   */
  function _tpParseSurd(s) {
    try { return Surd.parse(s); } catch (e) { return Surd.ZERO; }
  }

  /**
   * Normalize a helicity list from string or array form.
   * Accepts: "1/2, -1/2", [0.5, -0.5], "1, 0, -1"
   * Returns: array of numbers (multiples of 0.5).
   */
  function _normalizeHelicityList(list) {
    if (typeof list === 'string') {
      var parts = list.replace(/\s/g, '').split(',');
      var result = [];
      for (var i = 0; i < parts.length; i++) {
        var p = parts[i];
        if (p === '') continue;
        var slashIdx = p.indexOf('/');
        var val;
        if (slashIdx >= 0) {
          var num = parseInt(p.substring(0, slashIdx), 10);
          var den = parseInt(p.substring(slashIdx + 1), 10);
          val = den ? num / den : num;
        } else {
          val = parseFloat(p);
        }
        if (!isNaN(val) && isFinite(val)) result.push(val);
      }
      return result;
    }
    if (Array.isArray(list)) return list.slice();
    return [];
  }

  /**
   * Find index of helicity value in list (with tolerance).
   */
  function _helIdx(val, list) {
    var v = Math.round(val * 1e10) / 1e10;
    for (var i = 0; i < list.length; i++) {
      if (Math.abs(v - list[i]) < 1e-10) return i;
    }
    return -1;
  }

  return api;
})();
