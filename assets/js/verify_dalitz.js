// Verification for assets/js/dalitz.js
// Usage: node verify_dalitz.js  (run from assets/js)
'use strict';
const fs = require('fs');
const vm = require('vm');
for (const f of ['surd.js', 'cg.js', 'wigner-d.js', 'trig-poly.js', 'get-angle.js', 'angular-expression.js', 'projection.js', 'dalitz.js']) {
  vm.runInThisContext(fs.readFileSync(f, 'utf8'), { filename: f });
}
const { Surd } = global;

let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  \u2713 ' + name); }
  else { fail++; console.log('  \u2717 ' + name + (detail ? '  [' + detail + ']' : '')); }
}
function close(a, b, tol) { return Math.abs(a - b) < (tol === undefined ? 1e-9 : tol); }

const BASE = {
  mP: 1.87, mA: 0.14, mB: 0.14, mC: 0.14,
  mR: 0.78, GammaR: 0.15, JP: 0,
  widthMode: 'constant', barrier: false
};
function sess(over) {
  return global.buildDalitzSession(Object.assign({}, BASE, over));
}

console.log('=== kinematics ===');
{
  const q = global.dalitzTwoBodyMomentum(Math.pow(0.1396 + 0.1396, 2), 0.1396 * 0.1396, 0.1396 * 0.1396);
  check('two-body momentum ≈ 0 at threshold', close(q, 0, 1e-6), q);
  const q2 = global.dalitzTwoBodyMomentum(Math.pow(0.1396 + 0.1396 - 0.05, 2), 0.1396 * 0.1396, 0.1396 * 0.1396);
  check('below threshold -> NaN', isNaN(q2));
  const qR = global.dalitzTwoBodyMomentum(0.8955 * 0.8955, 0.4937 * 0.4937, 0.1396 * 0.1396);
  check('K*(892) q_R finite', isFinite(qR) && qR > 0.1, qR);
}

console.log('=== LS waves (projection-tool per-vertex format "l1,s1;l2,s2") ===');
{
  const w0 = global.dalitzWaves(0, 1, 0, 0, 0);
  check('parent-0, JR=1, spin-0: ["1,1;1,0"] (one wave)', w0.length === 1 && w0[0] === '1,1;1,0', JSON.stringify(w0));
  const w1 = global.dalitzWaves(1, 1, 0, 0, 0);
  check('parent-1, JR=1: three combined waves (0,1)(1,1)(2,1) × (1,0)', w1.length === 3 && w1[0] === '0,1;1,0' && w1[1] === '1,1;1,0' && w1[2] === '2,1;1,0', JSON.stringify(w1));
  const w00 = global.dalitzWaves(0, 0, 0, 0, 0);
  check('parent-0, JR=0: ["0,0;0,0"]', w00.length === 1 && w00[0] === '0,0;0,0', JSON.stringify(w00));
  const w12 = global.dalitzWaves(0, 0.5, 0.5, 0, 0.5);
  check('JR=1/2, JB=1/2, JC=0, JA=1/2: 4 waves', w12.length === 4 && w12[0] === '0,0;0,0.5' && w12[3] === '1,1;1,0.5', JSON.stringify(w12));
  const w32 = global.dalitzWaves(0, 1.5, 0.5, 0, 0.5);
  check('Δ-like (JR=3/2, N+π, JA=1/2): 4 waves', w32.length === 4 && w32[0] === '1,1;1,0.5', JSON.stringify(w32));
  const wBad = global.dalitzWaves(0, 1.5, 0.5, 0.5, 0.5);
  check('JR=3/2, JB=JC=1/2: no feasible wave (½⊗½ → integer S only)', wBad.length === 0, JSON.stringify(wBad));
  const wM = global.dalitzWaves(1, 1, 0.5, 0.5, 0);
  check('parent-1, JR=1, JB=JC=1/2: 3×4 = 12 waves', wM.length === 12, JSON.stringify(wM.length));
}

console.log('=== spin-0 daughters: angular part = Legendre^2 ===');
{
  function legendre(L, x) {
    if (L === 0) return 1;
    if (L === 1) return x;
    if (L === 2) return 0.5 * (3 * x * x - 1);
    if (L === 3) return 0.5 * (5 * x * x * x - 3 * x);
    return NaN;
  }
  for (const JR of [0, 1, 2]) {
    const lsKey = String(JR) + ',' + JR + ';' + JR + ',0';
    const s = sess({ JR: JR, JB: 0, JC: 0, lsKey: lsKey });
    check('session builds for JR=' + JR + ' (' + lsKey + ')', !s.error, s.error);
    if (s.error) continue;
    let ratioOk = true, ratio = 0;
    for (const x of [-0.9, -0.5, -0.1, 0.3, 0.7, 0.95]) {
      const w = s.angular(x);
      const P = legendre(JR, x);
      const r = w / (P * P);
      if (ratio === 0) ratio = r;
      else if (!close(r, ratio, 1e-6)) ratioOk = false;
    }
    check('W(cosθ) ∝ P_' + JR + '(cosθ)² (spin-0 daughters)', ratioOk);
  }
}

console.log('=== angular part vs direct CG + exact d-matrix formula ===');
{
  function cg(j1, m1, j2, m2, j3, m3) {
    const sd = global.computeCGExact(String(j1), String(m1), String(j2), String(m2), String(j3), String(m3));
    return sd.s * sd.p * Math.sqrt(sd.r) / sd.q;
  }
  function dNumeric(J, m1, m2, th) {
    const ws = global._getExactWignerDWeights(J, m1, m2);
    const sh = Math.sin(th / 2), ch = Math.cos(th / 2);
    let v = 0;
    for (const w of ws) {
      const sd = Surd.parse(w.weightStr);
      v += sd.s * sd.p * Math.sqrt(sd.r) / sd.q * Math.pow(sh, w.sinPow) * Math.pow(ch, w.cosPow);
    }
    return v;
  }
  // Independent R→2+3 vertex evaluation: W(θ) = Σ_{λR∈set} Σ_{λ2,λ3}
  // |√((2L+1)/(2JR+1))·⟨J2 λ2;J3 −λ3|S δ⟩·⟨L 0;S δ|JR δ⟩·d^{JR}_{λR,δ}(θ)|²
  function directW(JR, JB, JC, L, S, JA, cosT) {
    const th = Math.acos(cosT);
    const norm = Math.sqrt((2 * L + 1) / (2 * JR + 1));
    let total = 0;
    for (let lR = -JA; lR <= JA; lR += 1) {   // λR = λ1 (scalar parent)
      for (let lB = -JB; lB <= JB; lB += 1) {
        for (let lC = -JC; lC <= JC; lC += 1) {
          const d = lB - lC;
          const cgR = norm * cg(JB, lB, JC, -lC, S, d) * cg(L, 0, S, d, JR, d);
          if (cgR === 0) continue;
          total += cgR * cgR * Math.pow(dNumeric(JR, lR, d, th), 2);
        }
      }
    }
    return total;
  }
  // JP=0, JR=1, JB=JC=1/2, wave "1,1;1,1" (parent S-wave? parent (1,1); res (L=1,S=1))
  {
    const s = sess({ JR: 1, JB: 0.5, JC: 0.5, lsKey: '1,1;1,1' });
    check('session builds (JR=1, spin-1/2 daughters, P-wave)', !s.error, s.error);
    if (!s.error) {
      const ref = 0.2;
      let ok = true;
      for (const x of [-0.8, -0.3, 0.6, 0.9]) {
        const a = s.angular(x) / s.angular(ref);
        const b = directW(1, 0.5, 0.5, 1, 1, 0, x) / directW(1, 0.5, 0.5, 1, 1, 0, ref);
        if (!close(a, b, 1e-9)) { ok = false; console.log('    shape mismatch at cosθ=' + x, a, b); }
      }
      check('W(cosθ) ∝ Σ|cgR·d|² matches direct formula', ok);
    }
  }
  // JP=0, JR=1/2, JB=1/2, JC=0, JA=1/2, S-wave "0,0;0,0.5"
  {
    const s = sess({ JR: 0.5, JB: 0.5, JC: 0, JA: 0.5, lsKey: '0,0;0,0.5' });
    check('session builds (JR=1/2 S-wave)', !s.error, s.error);
    if (!s.error) {
      const ref = 0.1;
      let ok = true;
      for (const x of [-0.7, 0.8]) {
        const a = s.angular(x) / s.angular(ref);
        const b = directW(0.5, 0.5, 0, 0, 0.5, 0.5, x) / directW(0.5, 0.5, 0, 0, 0.5, 0.5, ref);
        if (!close(a, b, 1e-9)) { ok = false; console.log('    shape mismatch at cosθ=' + x, a, b); }
      }
      check('S-wave W(cosθ) matches direct formula', ok);
    }
  }
}

console.log('=== BW behavior ===');
{
  // constant width + barrier OFF → pure BW, peaked at mR: |BW(mR)|² = 1/(mR·Γ)²
  const s = sess({ JR: 1, JB: 0, JC: 0, lsKey: '1,1;1,0' });
  const peak = s.mProfile(0.78);
  check('constant BW peak (barrier off) = |BW(mR)|² = 1/(mR·Γ)²', close(peak * 0.78 * 0.78 * 0.15 * 0.15, 1, 1e-9), peak);
  // constant width + barrier off is J-INDEPENDENT (same line shape for all J)
  function spectrumPeak(J, wm, barrier) {
    const cfg = { mP: 1.86965, mA: 0.13957, mB: 0.13957, mC: 0.13957,
      mR: 0.77526, GammaR: 0.1478, JP: 0, JR: J, JB: 0, JC: 0,
      lsKey: J + ',' + J + ';' + J + ',0', widthMode: wm, barrier: barrier, radius: 3.0 };
    const s = global.buildDalitzSession(cfg);
    let mPeak = 0, vMax = -1;
    for (let m = 0.4; m <= 1.4; m += 0.002) {
      const v = s.mProfile(m);
      if (v > vMax) { vMax = v; mPeak = m; }
    }
    return mPeak;
  }
  const p0 = spectrumPeak(0, 'constant', false);
  const p5 = spectrumPeak(5, 'constant', false);
  check('constant+off: peak independent of J (J=0 vs J=5)', close(p0, p5, 1e-9) && close(p0, 0.77526, 0.01), p0 + ' vs ' + p5);
  // running width (barrier off): simple mass-dependent width, still L-independent
  const rp0 = spectrumPeak(0, 'running', false);
  const rp5 = spectrumPeak(5, 'running', false);
  check('running+off: peak independent of J', close(rp0, rp5, 1e-9), rp0 + ' vs ' + rp5);
  const sr = sess({ JR: 1, JB: 0, JC: 0, lsKey: '1,1;1,0', widthMode: 'running' });
  let mPeak = 0.78, vMax = 0;
  for (let m = 0.5; m <= 1.2; m += 0.002) {
    const v = sr.mProfile(m);
    if (v > vMax) { vMax = v; mPeak = m; }
  }
  check('running-width peak near mR', Math.abs(mPeak - 0.78) < 0.03, mPeak);
}

console.log('=== phase space & factorizability ===');
{
  const s = sess({ JR: 1, JB: 0, JC: 0, lsKey: '1,1;1,0', widthMode: 'running', barrier: true, radius: 3.0 });
  check('mMin = mB+mC', close(s.mMin, 0.28, 1e-9));
  check('mMax = mP−mA', close(s.mMax, 1.73, 1e-9));
  check('zero below mMin', s.eval(0.20, 0) === 0);
  check('zero above mMax', s.eval(1.80, 0) === 0);
  const f1 = s.mProfile(0.9), f2 = s.mProfile(1.1);
  const g05 = s.angular(0.5), g07 = s.angular(0.7);
  const r1 = s.eval(0.9, 0.5) / (f1 * g05);
  const r2 = s.eval(1.1, 0.7) / (f2 * g07);
  check('I(m,cosθ) = f(m)·g(cosθ) factorizes', close(r1, r2, 1e-9), r1 + ' vs ' + r2);
}

console.log('=== errors ===');
{
  const e1 = sess({ mP: 0.3, mA: 0.14, mB: 0.14, mC: 0.14, JR: 0, JB: 0, JC: 0, lsKey: '0,0;0,0' });
  check('empty phase space reported', !!e1.error, e1.error);
  const e2 = sess({ mR: 0.1, JR: 1, JB: 0, JC: 0, lsKey: '1,1;1,0', widthMode: 'running' });
  check('running width below threshold reported', !!e2.error, e2.error);
  // constant width below threshold is fine (no q_R dependence → no NaN)
  const e2b = sess({ mR: 0.1, JR: 1, JB: 0, JC: 0, lsKey: '1,1;1,0', widthMode: 'constant' });
  check('constant width below threshold OK (no NaN)', !e2b.error && isFinite(e2b.eval(0.9, 0.2)), e2b.error);
  const e3 = sess({ JR: 1, JB: 0, JC: 0, lsKey: '0,0;0,0' });
  check('infeasible wave reported', !!e3.error, e3.error);
}

console.log('=== Dalitz variables (s12, s23) mapping ===');
{
  const mP = 1.87, mA = 0.14, mB = 0.14, mC = 0.14;
  function checkPoint(s23) {
    const p = global.dalitzTwoBodyMomentum(mP * mP, mA * mA, s23);
    const ER = (mP * mP + s23 - mA * mA) / (2 * mP);
    const EA = (mP * mP + mA * mA - s23) / (2 * mP);
    const pR = { x: 0, y: 0, z: p }, pA = { x: 0, y: 0, z: -p };
    const g = ER / Math.sqrt(s23), bg = p / Math.sqrt(s23);
    const EAstar = g * EA + bg * p;
    const pAstarZ = -g * p - bg * EA;
    const q = global.dalitzTwoBodyMomentum(s23, mB * mB, mC * mC);
    const cosT = 0.3;
    const sinT = Math.sqrt(1 - cosT * cosT);
    const pB = { x: q * sinT, y: 0, z: q * cosT };
    const EBs = (s23 + mB * mB - mC * mC) / (2 * Math.sqrt(s23));
    const s12 = mA * mA + mB * mB + 2 * EAstar * EBs - 2 * (pAstarZ * pB.z);
    const ct = global.dalitzCosTheta(mP, mA, mB, mC, s12, s23);
    return Math.abs(ct - cosT) < 1e-9;
  }
  let ok = true;
  for (const s23 of [0.4, 0.6, 1.0, 1.4]) {
    if (!checkPoint(s23)) { ok = false; console.log('  mismatch at s23=' + s23); }
  }
  check('dalitzCosTheta inverts the 4-momentum construction (tf-pwa convention)', ok);
  const r = global.dalitzS12Range(mP, mA, mB, mC, 0.6);
  const c1 = global.dalitzCosTheta(mP, mA, mB, mC, r[0], 0.6);
  const c2 = global.dalitzCosTheta(mP, mA, mB, mC, r[1], 0.6);
  check('s12 range ends give |cosθ| = 1', close(Math.abs(c1), 1) && close(Math.abs(c2), 1), c1 + ', ' + c2);
  check('s12 outside range -> NaN', isNaN(global.dalitzCosTheta(mP, mA, mB, mC, r[0] - 0.1, 0.6)));
  const b = global.dalitzPlotBounds(mP, mA, mB, mC);
  check('plot bounds sane (s23Max = (mP−mA)²)', close(b.s23Max, Math.pow(1.73, 2), 1e-9) && b.s12Min < b.s12Max);
  // inverse: dalitzS12(s23, cosθ) -> dalitzCosTheta recovers cosθ
  let invOk = true;
  for (const s23 of [0.3, 0.6, 1.0, 1.8]) for (const ct of [-0.9, -0.2, 0.5, 0.95]) {
    const s12 = global.dalitzS12(mP, mA, mB, mC, s23, ct);
    const ct2 = global.dalitzCosTheta(mP, mA, mB, mC, s12, s23);
    if (!close(ct2, ct, 1e-9)) { invOk = false; console.log('  inverse fail', s23, ct, ct2); }
  }
  check('dalitzS12 inverts dalitzCosTheta', invOk);
  // sum rule: s12 + s13 + s23 = S
  const S = mP * mP + mA * mA + mB * mB + mC * mC;
  const s12t = global.dalitzS12(mP, mA, mB, mC, 0.6, 0.3);
  check('s12 + s13 + s23 = Σm²', close(s12t + (S - s12t - 0.6) + 0.6, S, 1e-9));
}

console.log('=== realistic example: D+ → ρ(770)π, ρ → ππ ===');
{
  const s = global.buildDalitzSession({
    mP: 1.86965, mA: 0.13957, mB: 0.13957, mC: 0.13957,
    mR: 0.77526, GammaR: 0.1478, JP: 0, JR: 1, JB: 0, JC: 0,
    lsKey: '1,1;1,0', widthMode: 'running', barrier: true, radius: 3.0
  });
  check('D+→ρπ session builds', !s.error, s.error);
  let vals = [];
  for (let m = 0.4; m <= 1.2; m += 0.005) vals.push([m, s.eval(m, 0.5)]);
  let mMax = 0, vMax = -1;
  for (const [m, v] of vals) if (v > vMax) { vMax = v; mMax = m; }
  check('ρ band peaks near 0.775', Math.abs(mMax - 0.775) < 0.06, mMax);
  const w0 = s.angular(0), w05 = s.angular(0.5), w1 = s.angular(1);
  check('angular ≈ cos²θ: W(0)≈0', close(w0, 0, 1e-9), w0);
  check('angular ≈ cos²θ: W(0.5)/W(1.0) = 1/4', close(w05 / w1, 0.25, 1e-9), w05 / w1);
}

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
