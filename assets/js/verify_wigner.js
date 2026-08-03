// Verification: Algebrite removal from wigner-d.js + page wiring.
// Usage: node verify_wigner.js  (run from assets/js)
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

console.log('=== no Algebrite references in wigner-d.js ===');
{
  const src = fs.readFileSync('wigner-d.js', 'utf8');
  check('no "Algebrite.run" calls', !/Algebrite\.run/.test(src));
  check('no "typeof Algebrite" gates', !/typeof Algebrite/.test(src));
}

console.log('=== runWignerDSanityChecks ===');
{
  const r = global.runWignerDSanityChecks();
  check('wigner-d sanity suite passes', r.passed === 12 && r.failed === 0, JSON.stringify(r));
}

console.log('=== numeric cross-check vs independent float d-matrix ===');
{
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
  let n = 0, bad = 0;
  for (const [j, m1, m2] of [[1, 1, 0], [1, 0, 0], [1.5, 1.5, -0.5], [2, 1, -1], [3, 2, 2], [3.5, 0.5, -2.5], [5, 2, -3], [6, 3, 0], [7.5, 7.5, -7.5]]) {
    for (const beta of [0.3, 0.9, Math.PI / 2, 2.1]) {
      const r = global.computeWignerD(String(j), String(m1), String(m2), String(beta));
      if (r.error) { bad++; if (bad < 4) console.log('    error', j, m1, m2, r.error); continue; }
      n++;
      if (!close(r.decimal, floatD(j, m1, m2, beta), 1e-9)) {
        bad++;
        if (bad < 4) console.log('    numeric mismatch', j, m1, m2, beta, r.decimal, floatD(j, m1, m2, beta));
      }
    }
  }
  check('numeric d-elements match independent float formula (' + n + ' cases)', bad === 0);
}

console.log('=== symbolic closed forms ===');
{
  const cases = [
    ['1/2', '1/2', '1/2', 'cos(beta/2)'],
    ['1/2', '1/2', '-1/2', '-sin(beta/2)'],
    ['1', '1', '1', '1/2+1/2*cos(beta)'],
    ['1', '1', '0', '-sqrt(2)/2*sin(beta)'],
    ['1', '1', '-1', '1/2-1/2*cos(beta)'],
    ['1', '0', '0', 'cos(beta)'],
    ['2', '2', '2', '3/8+1/2*cos(beta)+1/8*cos(2*beta)'],  // cos^4(beta/2)
  ];
  let ok = true;
  for (const [j, m1, m2, exp] of cases) {
    const r = global.computeWignerD(j, m1, m2, '');
    if (r.symbolic !== exp) { ok = false; console.log('    symbolic mismatch', j, m1, m2, 'got', r.symbolic, 'want', exp); }
  }
  check('symbolic forms match known closed forms', ok);
}

console.log('=== numeric beta parsing ===');
{
  const r1 = global.computeWignerD('1', '1', '1', 'pi/2');
  check('beta="pi/2": d^1_{1,1} = cos^2(pi/4) = 0.5', close(r1.decimal, 0.5), r1.decimal);
  const r2 = global.computeWignerD('1/2', '1/2', '1/2', 'pi');
  check('beta="pi": d^{1/2}_{1/2,1/2} = cos(pi/2) = 0', close(r2.decimal, 0, 1e-10), r2.decimal);
  const r3 = global.computeWignerD('1', '0', '0', '0.5');
  check('beta="0.5": d^1_{0,0} = cos(0.5)', close(r3.decimal, Math.cos(0.5)), r3.decimal);
  const r4 = global.computeWignerD('1', '0', '0', '2*pi/3');
  check('beta="2*pi/3": d^1_{0,0} = cos(2pi/3)', close(r4.decimal, Math.cos(2 * Math.PI / 3)), r4.decimal);
  const r5 = global.computeWignerD('1', '0', '0', 'alpha');
  check('beta="alpha" (junk) errors loudly', !!r5.error, JSON.stringify(r5));
  const r6 = global.computeWignerD('1', '0', '0', 'pi/2; drop table');
  check('beta="pi/2; drop table" rejected', !!r6.error, JSON.stringify(r6));
}

console.log('=== symbolic (empty beta) returns symbolic + latex, decimal NaN ===');
{
  const r = global.computeWignerD('1', '1', '0', '');
  check('symbolic decimal is NaN (not shown in UI)', isNaN(r.decimal));
  check('latex has trig + \\beta', /\\sin|\\cos/.test(r.latex) && /\\beta/.test(r.latex), r.latex);
  check('latex renders sqrt for d^1_{1,0}', /\\sqrt\{2\}/.test(r.latex), r.latex);
}

console.log('=== computeWignerDSimplified ===');
{
  const r = global.computeWignerDSimplified('1', '1', '0');
  check('simplified d^1_{1,0} = -sqrt(2)/2*sin(beta)', r.symbolic === '-sqrt(2)/2*sin(beta)', JSON.stringify(r));
  check('simplified groups non-empty', r.groups && r.groups.length > 0);
  check('simplified latex has \\sin', /\\sin/.test(r.latex), r.latex);
  const r2 = global.computeWignerDSimplified('3/2', '3/2', '-1/2');
  // numeric cross-check of the simplified symbolic at several beta
  let numOk = true;
  for (const b of [0.4, 1.3]) {
    // evaluate the fourier groups numerically
    const wts = global._getExactWignerDWeights(1.5, 1.5, -0.5);
    let sum = 0;
    for (const w of wts) {
      const s = Surd.parse(w.weightStr);
      const v = (w.weightStr[0] === '-' ? -1 : 1) * s.p * Math.sqrt(s.r) / s.q;
      sum += v * Math.pow(Math.sin(b / 2), w.sinPow) * Math.pow(Math.cos(b / 2), w.cosPow);
    }
    const d = global.computeWignerD('3/2', '3/2', '-1/2', String(b));
    if (!close(d.decimal, sum, 1e-9)) { numOk = false; console.log('    mismatch at', b, d.decimal, sum); }
  }
  check('simplified-able half-integer element consistent', numOk);
  const r3 = global.computeWignerDSimplified('1/2', '1/2', '0');  // invalid parity → zero
  check('invalid parity (m=0 for j=1/2) → zero', r3.symbolic === '0' && r3.groups.length === 0, JSON.stringify(r3));
}

console.log('=== computeWignerDMatrix ===');
{
  const r = global.computeWignerDMatrix('1', 'pi/2');
  check('matrix 3x3 for j=1', r.matrix.length === 3 && r.matrix[0].length === 3);
  const m = r.matrix;
  check('matrix d^1_{0,0}(pi/2) = cos(pi/2) = 0', close(m[1][1].decimal, 0, 1e-10));
  check('matrix d^1_{1,1}(pi/2) = 1/2', close(m[0][0].decimal, 0.5, 1e-9));
  check('matrix symmetric d^1_{1,-1}(pi/2)', close(m[0][2].decimal, m[2][0].decimal, 1e-9));
  // orthonormality: sum_m d^{1}_{m,m1}(b) * d^{1}_{m,m2}(b) = delta
  const b = 0.7;
  const mm = global.computeWignerDMatrix('1', '0.7');
  const dv = mm.matrix.map(row => row.map(c => c.decimal));
  let orth = 0;
  for (let m1 = 0; m1 < 3; m1++) for (let m2 = 0; m2 < 3; m2++) {
    let s = 0;
    for (let k = 0; k < 3; k++) s += dv[k][m1] * dv[k][m2];
    if (m1 === m2) orth += Math.abs(s - 1);
    else orth += Math.abs(s);
  }
  check('matrix orthogonality (j=1)', orth < 1e-8, orth);
  const half = global.computeWignerDMatrix('3/2', '');
  check('half-integer matrix symbolic works (4x4)', !half.error && half.matrix.length === 4, half.error);
  // symbolic matrix cell has latex
  check('symbolic matrix cell has latex', /\\cos/.test(half.matrix[0][0].latex), half.matrix[0][0].latex);
}

console.log('=== validation / errors ===');
{
  const r1 = global.computeWignerD('1/2', '1/2', '0', '');
  check('m=0 for j=1/2 → parity error', !!r1.error, JSON.stringify(r1));
  const r2 = global.computeWignerD('1', '2', '0', '');
  check('|m|>j → error', !!r2.error);
  const r3 = global.computeWignerD('bad', '0', '0', '');
  check('bad j → error', !!r3.error);
}

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
