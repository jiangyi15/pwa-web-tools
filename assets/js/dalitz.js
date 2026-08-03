/**
 * Dalitz Plot Calculator — isobar-model intensity for a resonance decay
 * R → 2 + 3 in a single partial wave, modulated by a Breit–Wigner.
 *
 *   I(s12, s23) = |BW(m)|² · (q/mR)^{2L} · F_L(q·d)² · W(cosθ_R)
 *
 * where m = m_{23} = √s23 is the resonance invariant mass, θ_R the helicity
 * angle of the R → 2 + 3 decay, q(m) the two-body momentum in the R rest
 * frame, and F_L the Blatt–Weisskopf barrier factor (z = q·d).
 *
 * The angular part W(cosθ_R) is the θ_R marginal of the exact angular
 * distribution of the full chain A → R + 1, R → 2 + 3, computed with the
 * SAME machinery as the Angular Projection tool:
 *   - AngularExpression.compute() builds the exact helicity amplitude
 *     (BigInt CG coefficients + exact Wigner-d half-angle weights) in the
 *     TrigPoly Fourier basis;
 *   - projection.js integrates over the unobserved parent angles and the
 *     azimuths exactly via the analytic integral tables.
 * LS waves are defined per vertex exactly like the Angular Projection tool
 * ("l1,s1;l2,s2" — parent vertex A → R + 1; resonance vertex R → 2 + 3),
 * and the parent-vertex production helicity weights are included naturally.
 *
 * The parent decay A → R + 1 enters through the mass window
 * m ∈ [m2 + m3, mA − m1] and the production helicity structure; W(cosθ_R)
 * is m-independent, so I factorizes into |BW(m)|²·(q/mR)^{2L}·F_L² times
 * W(cosθ_R).
 *
 * Dependencies: surd.js, cg.js, wigner-d.js, trig-poly.js, get-angle.js,
 *               angular-expression.js, projection.js
 */

// ============================================================================
// KINEMATICS
// ============================================================================

/** Källén triangle function λ(x, y, z). */
function _dalitzKallen(x, y, z) {
  return x * x + y * y + z * z - 2 * x * y - 2 * x * z - 2 * y * z;
}

/**
 * Two-body decay momentum: |p| = sqrt(λ(s, m1², m2²)) / (2·sqrt(s)).
 * Returns NaN when the invariant mass s is below the two-body threshold.
 * (A tiny negative λ from float roundoff at the exact threshold is clamped.)
 */
function dalitzTwoBodyMomentum(s, m1Sq, m2Sq) {
  var lam = _dalitzKallen(s, m1Sq, m2Sq);
  if (lam < 0) {
    if (lam > -1e-9) lam = 0;
    else return NaN;
  }
  return Math.sqrt(lam) / (2 * Math.sqrt(s));
}

/** Blatt–Weisskopf barrier factor F_L(z), z = q·d (d = barrier radius). */
function dalitzBarrier(L, z) {
  var z2 = z * z;
  switch (L) {
    case 0: return 1;
    case 1: return 1 / Math.sqrt(1 + z2);
    case 2: return 1 / Math.sqrt(9 + 3 * z2 + z2 * z2);
    case 3: return 1 / Math.sqrt(225 + 45 * z2 + 6 * z2 * z2 + z2 * z2 * z2);
    case 4: return 1 / Math.sqrt(11025 + 1575 * z2 + 135 * z2 * z2 + 10 * z2 * z2 * z2 + z2 * z2 * z2 * z2);
    default: return 1;
  }
}

// ============================================================================
// DALITZ VARIABLES — mapping between (s12, s23) and (m_BC, cos θ_R)
// ============================================================================

/**
 * cos of the helicity angle θ_R of B in the BC rest frame, as a function of
 * the Dalitz variables (s12, s23) for the chain P → A + R(→ B + C):
 *
 *   cosθ_R = (s12 − m_A² − m_B² − 2 E_A* E_B*) / (2 |p_A*| q)
 *
 * θ_R is measured from the resonance flight direction (the quantization
 * axis, opposite to A), with E_A*, E_B*, |p_A*|, q in the BC rest frame.
 * NaN outside the physical Dalitz region.
 */
function dalitzCosTheta(mP, mA, mB, mC, s12, s23) {
  var m = Math.sqrt(s23);
  var eA = (mP * mP - s23 - mA * mA) / (2 * m);
  var eB = (s23 + mB * mB - mC * mC) / (2 * m);
  var pA = dalitzTwoBodyMomentum(s23, mP * mP, mA * mA);   // |p_A*| in BC frame
  var q = dalitzTwoBodyMomentum(s23, mB * mB, mC * mC);
  if (!isFinite(pA) || !isFinite(q) || pA < 1e-12 || q < 1e-12) return NaN;
  var ct = (s12 - mA * mA - mB * mB - 2 * eA * eB) / (2 * pA * q);
  // Outside the physical Dalitz region (|cosθ| > 1): not a valid point.
  if (Math.abs(ct) > 1 + 1e-9) return NaN;
  return Math.max(-1, Math.min(1, ct));
}

/**
 * s12 = m²(1,2) as a function of (s23, cos θ_R) — the inverse of
 * dalitzCosTheta:
 *
 *   s12 = m1² + m2² + 2 E1* E2* + 2 |p1*| q cosθ_R
 *
 * with E1*, E2*, |p1*|, q in the 2+3 rest frame. NaN when s23 is outside the
 * physical region.
 */
function dalitzS12(mP, mA, mB, mC, s23, cosT) {
  var m = Math.sqrt(s23);
  var e1 = (mP * mP - s23 - mA * mA) / (2 * m);
  var e2 = (s23 + mB * mB - mC * mC) / (2 * m);
  var p1 = dalitzTwoBodyMomentum(s23, mP * mP, mA * mA);
  var q = dalitzTwoBodyMomentum(s23, mB * mB, mC * mC);
  if (!isFinite(p1) || !isFinite(q)) return NaN;
  return mA * mA + mB * mB + 2 * e1 * e2 + 2 * p1 * q * cosT;
}

/**
 * Allowed s12 = m_AB² interval at fixed s23 = m_BC² (the Dalitz boundary at
 * that slice). Returns [s12min, s12max] or null when s23 is outside the
 * physical region.
 */
function dalitzS12Range(mP, mA, mB, mC, s23) {
  var m = Math.sqrt(s23);
  var eA = (mP * mP - s23 - mA * mA) / (2 * m);
  var eB = (s23 + mB * mB - mC * mC) / (2 * m);
  var pA = dalitzTwoBodyMomentum(s23, mP * mP, mA * mA);   // |p_A*| in BC frame
  var q = dalitzTwoBodyMomentum(s23, mB * mB, mC * mC);
  if (!isFinite(pA) || !isFinite(q)) return null;
  var em = eA + eB;
  return [em * em - (pA + q) * (pA + q), em * em - (pA - q) * (pA - q)];
}

/**
 * Global plotting bounds of the Dalitz region in (s12, s23).
 * Returns { s12Min, s12Max, s23Min, s23Max } (GeV²).
 */
function dalitzPlotBounds(mP, mA, mB, mC) {
  var s23Min = Math.pow(mB + mC, 2);
  var s23Max = Math.pow(mP - mA, 2);
  var s12Min = Infinity, s12Max = -Infinity;
  for (var i = 0; i <= 200; i++) {
    var s23 = s23Min + (s23Max - s23Min) * i / 200;
    var r = dalitzS12Range(mP, mA, mB, mC, s23);
    if (!r) continue;
    if (r[0] < s12Min) s12Min = r[0];
    if (r[1] > s12Max) s12Max = r[1];
  }
  return { s12Min: s12Min, s12Max: s12Max, s23Min: s23Min, s23Max: s23Max };
}

// ============================================================================
// ANGULAR DISTRIBUTION — reused from the Angular Projection pipeline
// ============================================================================

/**
 * Available (L, S) waves of the full chain A → R + 1, R → 2 + 3, in the
 * SAME per-vertex LS format as the Angular Projection tool: each key is
 * "l1,s1;l2,s2" (parent vertex A → R + 1; resonance vertex R → 2 + 3).
 * Computed by building the exact angular distribution (AngularExpression)
 * and returning the LS keys that have non-zero amplitude.
 *
 * @param {number} JP    parent spin (A)
 * @param {number} JR    resonance spin (R)
 * @param {number} JB, JC  resonance daughter spins (2, 3)
 * @param {number} JA    spectator spin (1)
 * @returns {Array<string>} lsKeys like "1,1;1,0" or "0,1;1,0"
 */
function dalitzWaves(JP, JR, JB, JC, JA) {
  var tree = {
    j: String(JP),
    children: [
      { j: String(JR), children: [{ j: String(JB) }, { j: String(JC) }] },
      { j: String(JA) }
    ]
  };
  var res = AngularExpression.compute(tree, { showInterference: false });
  if (res.error) return [];
  return (res.lsKeyList || []).slice();
}

// ============================================================================
// DALITZ SESSION
// ============================================================================

// Cache of the exact angular part (the θ_R marginal projection). It depends
// ONLY on the spins and the selected LS wave — NOT on masses or width — so
// it is computed once and reused whenever the resonance mass/width (or the
// BW/barrier options) change.
var _dalitzAngularCache = {};

function _dalitzAngularCacheKey(config, JP, JA, lsKey) {
  return [JP, config.JR, config.JB, config.JC, JA, lsKey,
          JSON.stringify(config.helicityFilters || null)].join('|');
}

/**
 * Build a Dalitz session — precomputes everything that does not depend on
 * the plot point.
 *
 * @param {Object} config {
 *   mP, mA,          // parent + spectator (particle 1) masses — mass window
 *   mR, GammaR, JR,  // resonance (R → 2 + 3)
 *   mB, mC, JB, JC,  // resonance daughters (particles 2, 3)
 *   JP, JA,          // parent spin (A) and spectator spin (1)
 *   lsKey,           // single partial wave "l1,s1;l2,s2" (from dalitzWaves)
 *   widthMode,       // 'running' | 'constant'
 *   barrier,         // boolean: include Blatt–Weisskopf barrier
 *   radius,          // barrier radius d (GeV^-1)
 *   helicityFilters  // optional {path: "λ-list"} for the full chain
 * }
 * @returns {{eval: Function, angular: Function, mProfile: Function,
 *            waves: Array<string>, qR: number, mMin: number, mMax: number,
 *            error: string}}
 */
function buildDalitzSession(config) {
  var mP = config.mP, mA = config.mA;
  var mR = config.mR, gR = config.GammaR, JR = config.JR;
  var mB = config.mB, mC = config.mC, JB = config.JB, JC = config.JC;
  var lsKey = config.lsKey;
  var JP = (config.JP === undefined) ? 0 : config.JP;
  var JA = (config.JA === undefined) ? 0 : config.JA;

  // Orbital momentum of the RESONANCE vertex (last "l" of "l1,s1;l2,s2"),
  // used for the momentum/barrier factors and the running-width exponent.
  var lKeyParts = lsKey.split(';');
  var L = parseInt(lKeyParts[lKeyParts.length - 1].split(',')[0], 10);

  var mMin = mB + mC;
  var mMax = mP - mA;
  if (mMax <= mMin + 1e-9) {
    return { error: 'Phase space is empty: mP − mA ≤ mB + mC.' };
  }

  var qR = dalitzTwoBodyMomentum(mR * mR, mB * mB, mC * mC);
  if (isNaN(qR)) {
    // Constant width never uses q_R, so a below-threshold m_R is fine (no
    // NaN). Only the running width's (q/q_R)^{2L+1} factor diverges — reject
    // it when clearly below threshold (a hair below, from float roundoff at
    // the boundary, is clamped to q_R = 0).
    if (config.widthMode === 'running' && mR < mB + mC - 1e-4) {
      return { error: 'Resonance mass is below the B + C threshold.' };
    }
    qR = 0;
  }

  // Angular part — the θ_R marginal of the exact full-chain angular
  // distribution A → R + 1, R → 2 + 3, in the selected combined LS wave.
  // Uses the SAME machinery as the Angular Projection tool (which defines
  // LS per vertex: "l1,s1;l2,s2") — computeProjections integrates over the
  // unobserved parent angles and azimuths exactly (integral tables).
  // The parent-vertex production helicity weights are included naturally.
  //
  // This depends only on the spins + wave, so it is CACHED: changing the
  // resonance mass/width (or barrier options) does not recompute the exact
  // angular formula — only the cheap BW × momentum factors below.
  var angKey = _dalitzAngularCacheKey(config, JP, JA, lsKey);
  var ang = _dalitzAngularCache[angKey];
  if (!ang) {
    var tree = {
      j: String(JP),
      children: [
        { j: String(JR), children: [{ j: String(JB) }, { j: String(JC) }] },
        { j: String(JA) }
      ]
    };
    var res = AngularExpression.compute(tree, {
      showInterference: false,
      helicityFilters: config.helicityFilters || null
    });
    ang = { error: res.error || null, thetaProj: null, waves: (res.lsKeyList || []).slice() };
    if (!res.error) {
      var lsMap = res.lsMap || {};
      var fourier = (lsMap[lsKey] && lsMap[lsKey].fourier) ? lsMap[lsKey].fourier : null;
      if (fourier) {
        // R → 2 + 3 is the second decay vertex → its angle is theta_1.
        ang.thetaProj = computeProjections(fourier, 2, null)['theta_1'];
      } else {
        ang.error = 'Partial wave ' + lsKey + ' is not available for these spins.';
      }
    }
    _dalitzAngularCache[angKey] = ang;
  }
  if (ang.error) return { error: ang.error };
  var thetaProj = ang.thetaProj;

  function angular(cosT) {
    var theta = Math.acos(Math.max(-1, Math.min(1, cosT)));
    return evaluateProjection(thetaProj, theta);
  }

  // |BW(m)|² — running (mass-dependent) or constant width.
  // With barriers ON the running width carries the L-wave threshold factor
  // (q/qR)^{2L+1} AND the Blatt–Weisskopf ratio [F_L(q)/F_L(qR)]², so it is
  // consistent with the (q/mR)^{2L}·F_L² factor in the amplitude. With
  // barriers OFF the width is Γ(m) = ΓR·(mR/m), independent of L.
  function bw2(m, q) {
    var mR2 = mR * mR, mSq = m * m;
    var G = gR;
    if (config.widthMode === 'running') {
      G = gR * (mR / m);
      if (config.barrier) {
        // (q/qR)^{2L+1} diverges when the resonance sits at (or below) the
        // threshold (qR→0), which zeroes the whole plot and makes the width
        // irrelevant. In that case keep the simple mass-dependent width.
        if (qR > 1e-9) {
          G *= Math.pow(q / qR, 2 * L + 1) *
               Math.pow(dalitzBarrier(L, q * config.radius) / dalitzBarrier(L, qR * config.radius), 2);
        }
      }
    }
    var den = mR2 - mSq;
    return 1 / (den * den + mR * mR * G * G);
  }

  // Momentum/barrier factor. With barriers OFF there is no L-dependent form
  // factor at all — a pure Breit-Wigner (J-independent line shape).
  function formFactor(m, q) {
    if (!config.barrier) return 1;
    var f = dalitzBarrier(L, q * config.radius);
    return Math.pow(q / mR, 2 * L) * f * f;
  }

  // Mass profile f(m) = |BW(m)|² · (q/mR)^{2L} · F_L(q·d)²
  function mProfile(m) {
    if (m < mMin - 1e-9 || m > mMax + 1e-9) return 0;
    var q = dalitzTwoBodyMomentum(m * m, mB * mB, mC * mC);
    if (!isFinite(q)) return 0;
    return formFactor(m, q) * bw2(m, q);
  }

  // Full intensity at a plot point
  function evalPoint(m, cosT) {
    if (m < mMin - 1e-9 || m > mMax + 1e-9) return 0;
    var q = dalitzTwoBodyMomentum(m * m, mB * mB, mC * mC);
    if (!isFinite(q)) return 0;
    return angular(cosT) * formFactor(m, q) * bw2(m, q);
  }

  return {
    eval: evalPoint,
    angular: angular,
    mProfile: mProfile,
    waves: ang.waves.slice(),
    qR: qR,
    mMin: mMin,
    mMax: mMax,
    config: config
  };
}

// ============================================================================
// NODE/BROWSER EXPORT
// ============================================================================

if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    buildDalitzSession: buildDalitzSession,
    dalitzWaves: dalitzWaves,
    dalitzTwoBodyMomentum: dalitzTwoBodyMomentum,
    dalitzBarrier: dalitzBarrier,
    dalitzCosTheta: dalitzCosTheta,
    dalitzS12: dalitzS12,
    dalitzS12Range: dalitzS12Range,
    dalitzPlotBounds: dalitzPlotBounds
  };
}
