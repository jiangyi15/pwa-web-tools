// Verification script for the 4 critical fixes.
// Loads the libraries exactly like the browser does (shared globals via
// vm.runInThisContext, in the documented script-tag order).
// Usage: node verify_fixes.js  (run from assets/js)
'use strict';

const fs = require('fs');
const vm = require('vm');

for (const f of ['surd.js', 'cg.js', 'wigner-d.js', 'trig-poly.js', 'get-angle.js', 'angular-expression.js']) {
  vm.runInThisContext(fs.readFileSync(f, 'utf8'), { filename: f });
}

const { Surd } = global;

let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  \u2713 ' + name); }
  else { fail++; console.log('  \u2717 ' + name + (detail ? '  [' + detail + ']' : '')); }
}
function close(a, b, tol) { return Math.abs(a - b) < (tol === undefined ? 1e-9 : tol); }

console.log('=== Fix 1: decimal parsing ===');
{
  const s = global.toSpin('0.5');
  check('toSpin("0.5") -> 1/2', s.num === 1 && s.den === 2 && close(s.value, 0.5), JSON.stringify(s));
  const s2 = global.toSpin('1.5');
  check('toSpin("1.5") -> 3/2', s2.num === 3 && s2.den === 2, JSON.stringify(s2));
  const s3 = global.toSpin('-0.5');
  check('toSpin("-0.5") -> -1/2', s3.num === -1 && s3.den === 2, JSON.stringify(s3));
  let threw = false;
  try { global.toSpin('1.5junk'); } catch (e) { threw = true; }
  check('toSpin("1.5junk") throws', threw);
  const s4 = global.toSpin('3');
  check('toSpin("3") -> 3/1', s4.num === 3 && s4.den === 1, JSON.stringify(s4));

  const p = global.parseQuantumNumber('0.5');
  check('parseQuantumNumber("0.5") -> num 1 den 2', p.numerator === 1 && p.denominator === 2, JSON.stringify(p));
  const p2 = global.parseQuantumNumber('1.5');
  check('parseQuantumNumber("1.5") -> num 3 den 2', p2.numerator === 3 && p2.denominator === 2, JSON.stringify(p2));
  const p3 = global.parseQuantumNumber('-0.5');
  check('parseQuantumNumber("-0.5") -> -1/2', p3.numerator === -1 && p3.denominator === 2, JSON.stringify(p3));
  threw = false;
  try { global.parseQuantumNumber('2.5junk'); } catch (e) { threw = true; }
  check('parseQuantumNumber("2.5junk") throws', threw);

  const r = global.computeCG('1.5', '0.5', '0.5', '-0.5', '1', '0');
  check('computeCG(1.5,0.5,0.5,-0.5,1,0) = sqrt(2)/2', close(r.decimal, Math.SQRT1_2), JSON.stringify(r));
}

console.log('=== Fix 2: (j,m) parity validation ===');
{
  const r1 = global.computeCGExact('1/2', '1/2', '1/2', '0', '1', '1/2');
  check('CG <1/2 1/2; 1/2 0|1 1/2> = 0 (was 1)', r1.isZero());
  const r2 = global.computeCGExact('1', '0', '1/2', '0', '3/2', '0');
  check('CG <1 0; 1/2 0|3/2 0> = 0 (was 2sqrt(3)/3)', r2.isZero());
  const rules = global.checkSelectionRules('1/2', '1/2', '1/2', '0', '1', '1/2');
  check('checkSelectionRules rejects m=0 for j=1/2', rules.valid === false, rules.message);
  const ok1 = global.computeCGExact('1/2', '1/2', '1/2', '-1/2', '1', '0');
  check('valid CG: <1/2 1/2; 1/2 -1/2|1 0> = 1/sqrt(2)',
        !ok1.isZero() && close(ok1.s * ok1.p * Math.sqrt(ok1.r) / ok1.q, Math.SQRT1_2));
  const ok2 = global.computeCGExact('3/2', '1/2', '1', '1', '5/2', '3/2');
  check('valid CG: <3/2 1/2; 1 1|5/2 3/2> = sqrt(3/5) (sympy-verified)',
        !ok2.isZero() && close(ok2.s * ok2.p * Math.sqrt(ok2.r) / ok2.q, Math.sqrt(3 / 5)));
  const ok3 = global.computeCGExact('1', '1', '1', '0', '1', '1');
  check('valid CG: <1 1; 1 0|1 1> = 1/sqrt(2)',
        !ok3.isZero() && close(ok3.s * ok3.p * Math.sqrt(ok3.r) / ok3.q, Math.SQRT1_2));
}

console.log('=== Fix 3: exact Wigner-d weights (BigInt) ===');
{
  const w = global._getExactWignerDWeights(8, 8, 8);
  check('d^8_{8,8}: one weight, wStr "1"', w.length === 1 && w[0].weightStr === '1', JSON.stringify(w));
  const w2 = global._getExactWignerDWeights(9, 5, 9);
  const hasCorrect = w2.some(x => x.weightStr === '6*sqrt(85)');
  check('d^9_{5,9}: contains weight "6*sqrt(85)"', hasCorrect, JSON.stringify(w2));
  const w3 = global._getExactWignerDWeights(9.5, 9.5, 0.5);
  const hasSqrt = w3.some(x => x.weightStr === 'sqrt(92378)' || x.weightStr === '-sqrt(92378)');
  check('d^{9.5}_{9.5,0.5}: weight magnitude sqrt(92378)', hasSqrt, JSON.stringify(w3));
  let allParse = true;
  const w4 = global._getExactWignerDWeights(10, 10, 0);
  for (const x of w4) { try { Surd.parse(x.weightStr); } catch (e) { allParse = false; } }
  check('d^10_{10,0}: all weights parse as Surd (no corrupt radicand)', allParse && w4.length > 0);

  function floatD(j, m1, m2, beta) {
    let total = 0;
    const twoJ = Math.round(2 * j);
    const jpm1 = j + m1, jmm1 = j - m1, jpm2 = j + m2, jmm2 = j - m2;
    function fact(n) { let r = 1; for (let i = 2; i <= n; i++) r *= i; return r; }
    const num = fact(jpm1) * fact(jmm1) * fact(jpm2) * fact(jmm2);
    for (let l = 0; l <= twoJ; l++) {
      const k = (l + m2 - m1) / 2;
      if (Math.abs(k - Math.round(k)) > 1e-10) continue;
      const kk = Math.round(k);
      if (kk < Math.max(0, m2 - m1) || kk > Math.min(jmm1, jpm2)) continue;
      const den = fact(jmm1 - kk) * fact(jpm2 - kk) * fact(m1 - m2 + kk) * fact(kk);
      const sign = ((m1 - m2 + kk) % 2 === 0) ? 1 : -1;
      total += sign * Math.sqrt(num) / den * Math.pow(Math.sin(beta / 2), l) * Math.pow(Math.cos(beta / 2), twoJ - l);
    }
    return total;
  }
  let weightOk = true;
  for (const [j, m1, m2] of [[1, 1, 0], [2, 1, -1], [3, 2, 1], [5, 2, -3], [6, 3, 0]]) {
    const wts = global._getExactWignerDWeights(j, m1, m2);
    const beta = 0.9;
    let sum = 0;
    for (const x of wts) {
      const s = Surd.parse(x.weightStr);
      const val = (x.weightStr[0] === '-' ? -1 : 1) * s.p * Math.sqrt(s.r) / s.q;
      sum += val * Math.pow(Math.sin(beta / 2), x.sinPow) * Math.pow(Math.cos(beta / 2), x.cosPow);
    }
    if (!close(sum, floatD(j, m1, m2, beta), 1e-8)) { weightOk = false; console.log('    mismatch j=' + j, sum, floatD(j, m1, m2, beta)); }
  }
  check('weights reproduce float d-matrix for j<=6', weightOk);
  const t0 = Date.now();
  global._getExactWignerDWeights(10, 5, -3);
  check('d^10_{5,-3} fast (<300ms)', (Date.now() - t0) < 300);
}

console.log('=== Fix 4: interference contamination ===');
{
  const evalMap = (map, th, ph) => {
    let tot = 0;
    for (const key in map) {
      if (map[key].isEmpty()) continue;
      let c = 0;
      for (const t of map[key].terms()) c += t.s * t.p * Math.sqrt(t.r) / t.q;
      let fac = 1;
      if (key !== '1') {
        for (const pt of key.split(',')) {
          const m = pt.match(/([\w\d_]+):(\w+)\|(-?[\d.]+)/);
          const ang = m[1].indexOf('theta') === 0 ? th : ph;
          const k = parseFloat(m[3]);
          fac *= (m[2] === 'cos' ? Math.cos(k * ang / 2) : Math.sin(k * ang / 2));
        }
      }
      tot += c * fac;
    }
    return tot;
  };

  // T1 = cos(theta/2)*cos(phi/2) + i*sin(phi/2)
  // T2 = cos(theta/2)           + i*sin(phi/2)*cos(phi/2)
  const t1 = new global.TrigPoly();
  t1.addPowerTerm(Surd.ONE, false, [{ name: 'theta_0', sp: 0, cp: 1 }], [{ name: 'phi_1', func: 'cos', k: 2 }]);
  t1.addPowerTerm(Surd.ONE, true, [], [{ name: 'phi_1', func: 'sin', k: 2 }]);
  const t2 = new global.TrigPoly();
  t2.addPowerTerm(Surd.ONE, false, [{ name: 'theta_0', sp: 0, cp: 1 }], []);
  t2.addPowerTerm(Surd.ONE, true, [], [{ name: 'phi_1', func: 'sin', k: 2 }, { name: 'phi_1', func: 'cos', k: 2 }]);

  const th = 1.1, ph = 0.7;
  const re1 = Math.cos(th / 2) * Math.cos(ph), im1 = Math.sin(ph);
  const re2 = Math.cos(th / 2), im2 = Math.sin(ph) * Math.cos(ph);
  const trueRe = re1 * re2 + im1 * im2;
  const trueIm = im1 * re2 - re1 * im2;

  const prod = t1.mul(t2, { conjugateSecond: true });
  const reMap = prod.clone(); reMap.expand('real');
  const imMap = prod.clone(); imMap.expand('imag');
  const allMap = prod.clone(); allMap.expand();

  check('expand("real") = Re[T1 T2*]', close(evalMap(reMap._fourier, th, ph), trueRe, 1e-12), evalMap(reMap._fourier, th, ph) + ' vs ' + trueRe);
  check('expand("imag") = Im[T1 T2*]', close(evalMap(imMap._fourier, th, ph), trueIm, 1e-12), evalMap(imMap._fourier, th, ph) + ' vs ' + trueIm);
  const allVal = evalMap(allMap._fourier, th, ph);
  check('expand() (old) = Re+Im; differs from Re (the original bug)', close(allVal, trueRe + trueIm, 1e-12) && Math.abs(allVal - trueRe) > 1e-9);
  check('|T|² path: expand("real") = re^2+im^2', (() => {
    const sq = t1.mul(t1, { conjugateSecond: true });
    sq.expand('real');
    return close(evalMap(sq._fourier, th, ph), re1 * re1 + im1 * im1, 1e-12);
  })());

  // Full-pipeline end-to-end: J=0 -> (1 -> 0+0) + (1 -> 0+0)
  const tree = {
    j: '0',
    children: [
      { j: '1', children: [{ j: '0' }, { j: '0' }] },
      { j: '1', children: [{ j: '0' }, { j: '0' }] }
    ]
  };
  const res = global.AngularExpression.compute(tree, { showInterference: true });
  check('compute runs without error', !res.error, res.error);
  if (!res.error) {
    // 3 vertices: root (0→1+1) + two children (1→0+0)
    const theta = [1.2, 0.8, 0.5], phi = [0.3, 1.5, 0.7];
    const chiVal = phi[2];  // phi_2 is renamed to 'chi' by phiCombine

    const evalTerms = (terms, fixIdx, chiIdx) => {
      let re = 0, im = 0;
      for (const t of terms) {
        const s = Surd.parse(t.s);
        let v = s.s * s.p * Math.sqrt(s.r) / s.q;
        for (const e of (t._tbatch || [])) {
          const thv = theta[e.idx];
          if (e.sp) v *= Math.pow(Math.sin(thv / 2), e.sp);
          if (e.cp) v *= Math.pow(Math.cos(thv / 2), e.cp);
        }
        for (const p of (t._pbatch || [])) {
          // apply the same phiCombine substitution as the fourier path:
          // phi_fixIdx := 0 (eliminated), phi_chiIdx := chi
          let ang;
          if (p.idx === fixIdx) ang = 0;
          else if (p.idx === chiIdx) ang = chiVal;
          else ang = phi[p.idx];
          const k = p.pm || 0;
          if (p.pf === 'cos') v *= Math.cos(k * ang / 2);
          else if (p.pf === 'sin') v *= Math.sin(k * ang / 2);
        }
        if (t.im) im += v; else re += v;
      }
      return { re, im };
    };

    const fixIdx = res.phiCombine ? res.phiCombine.fixIdx : -1;
    const chiIdx = res.phiCombine ? res.phiCombine.chiIdx : -1;

    const allH = Object.keys(res.combinedRaw);
    const lsKeys = res.lsKeyList;

    let lsTotal = 0, pairTotal = 0;
    for (const lk of lsKeys) {
      let re = 0, im = 0;
      for (const hk of allH) {
        const terms = res.combinedRaw[hk] && res.combinedRaw[hk][lk];
        if (!terms) continue;
        const c = evalTerms(terms, fixIdx, chiIdx);
        re += c.re; im += c.im;
      }
      lsTotal += re * re + im * im;
    }
    for (let i = 0; i < lsKeys.length; i++) {
      for (let j = i + 1; j < lsKeys.length; j++) {
        let re1 = 0, im1 = 0, re2 = 0, im2 = 0;
        for (const hk of allH) {
          const t1 = res.combinedRaw[hk] && res.combinedRaw[hk][lsKeys[i]];
          const t2 = res.combinedRaw[hk] && res.combinedRaw[hk][lsKeys[j]];
          if (t1) { const c = evalTerms(t1, fixIdx, chiIdx); re1 += c.re; im1 += c.im; }
          if (t2) { const c = evalTerms(t2, fixIdx, chiIdx); re2 += c.re; im2 += c.im; }
        }
        pairTotal += 2 * (re1 * re2 + im1 * im2);
      }
    }
    const Itrue = lsTotal + pairTotal;

    const evalFMap = (fmap) => {
      let tot = 0;
      for (const key in fmap) {
        const sum = fmap[key];
        if (sum.isEmpty()) continue;
        let c = 0;
        for (const t of sum.terms()) c += t.s * t.p * Math.sqrt(t.r) / t.q;
        let fac = 1;
        if (key !== '1') {
          for (const pt of key.split(',')) {
            const m = pt.match(/([\w\d_]+):(\w+)\|(-?[\d.]+)/);
            const name = m[1];
            let ang;
            if (name.indexOf('theta_') === 0) ang = theta[parseInt(name.substring(6))];
            else if (name === 'chi') ang = chiVal;
            else ang = phi[parseInt(name.substring(4))];
            const k = parseFloat(m[3]);
            fac *= (m[2] === 'cos' ? Math.cos(k * ang / 2) : Math.sin(k * ang / 2));
          }
        }
        tot += c * fac;
      }
      return tot;
    };

    let lsLib = 0, ifLib = 0;
    for (const lk of lsKeys) lsLib += evalFMap(res.lsMap[lk].fourier);
    for (const pk in res.interfMap) ifLib += evalFMap(res.interfMap[pk].fourier);

    check('lsMap == sum |T_ls|^2', close(lsLib, lsTotal, 1e-9), lsLib + ' vs ' + lsTotal);
    check('interfMap == sum 2Re (fixed)', close(ifLib, pairTotal, 1e-9), ifLib + ' vs ' + pairTotal);
    check('lsMap + interfMap == |sum T|^2', close(lsLib + ifLib, Itrue, 1e-9), (lsLib + ifLib) + ' vs ' + Itrue);

    let pairOk = true;
    for (let i = 0; i < lsKeys.length && pairOk; i++) {
      for (let j = i + 1; j < lsKeys.length && pairOk; j++) {
        const pk = lsKeys[i] + '\u00d7' + lsKeys[j];
        const lib = evalFMap(res.interfMap[pk].fourier);
        let re1 = 0, im1 = 0, re2 = 0, im2 = 0;
        for (const hk of allH) {
          const t1 = res.combinedRaw[hk] && res.combinedRaw[hk][lsKeys[i]];
          const t2 = res.combinedRaw[hk] && res.combinedRaw[hk][lsKeys[j]];
          if (t1) { const c = evalTerms(t1, fixIdx, chiIdx); re1 += c.re; im1 += c.im; }
          if (t2) { const c = evalTerms(t2, fixIdx, chiIdx); re2 += c.re; im2 += c.im; }
        }
        const num = 2 * (re1 * re2 + im1 * im2);
        if (!close(lib, num, 1e-9)) { pairOk = false; console.log('    pair mismatch', pk, lib, num); }
      }
    }
    check('each interfMap pair == numeric 2Re', pairOk);
    check('interfImagMap still computed (2Im path intact)', Object.keys(res.interfImagMap).length > 0);
  }
}

console.log('=== regression: sanity checks ===');
{
  const r1 = global.runSanityChecksExact();
  check('runSanityChecksExact all pass', r1.passed === 10 && r1.failed === 0);
}

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
