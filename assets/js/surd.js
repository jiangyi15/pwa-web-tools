/**
 * Surd — exact representation of sign * p * sqrt(r) / q
 * 
 * All quantum-number coefficients in the helicity amplitude pipeline
 * (Wigner-d weights, CG coefficients, half-angle expansion) are of 
 * this form, where p, q, r are non-negative integers:
 * 
 *   r = 1  →  pure rational:  sign * p / q
 *   r > 1  →  quadratic surd: sign * p * sqrt(r) / q
 * 
 * Construction normalises automatically:
 *   - perfect squares extracted from r
 *   - fraction p/q reduced by gcd
 *   - sign absorbed into p (p always non-negative)
 *   - p = 0 → ZERO
 */

function Surd(sign, p, q, r) {
  // ---------------------------------------------------------------
  // Normalise
  // ---------------------------------------------------------------
  if (q === undefined || q === 0) q = 1;
  if (r === undefined) r = 1;

  // Zero
  if (p === 0) {
    this.s = 1; this.p = 0; this.q = 1; this.r = 1;
    return;
  }

  // Sign
  if (q < 0) { sign = -sign; q = -q; }
  this.s = (sign < 0 || sign === undefined) ? -1 : 1;

  // Extract perfect squares from radicand
  var sq = 1, rem = r;
  for (var i = 2; i * i <= rem; i++) {
    while (rem % (i * i) === 0) { sq *= i; rem /= (i * i); }
  }
  p *= sq;
  r = rem;

  // Reduce fraction
  var g = Surd._gcd(p, q);
  this.p = p / g;
  this.q = q / g;
  this.r = r;
}

// ---- Internal helpers (defined before static constants!) --------

Surd._gcd = function (a, b) {
  a = Math.abs(a); b = Math.abs(b);
  while (b) { var t = b; b = a % b; a = t; }
  return a;
};

// ---- Static constants -------------------------------------------

Surd.ZERO = Object.freeze(new Surd(1, 0, 1, 1));
Surd.ONE  = Object.freeze(new Surd(1, 1, 1, 1));

// ---- Query methods -----------------------------------------------

Surd.prototype.isZero = function () { return this.p === 0; };
Surd.prototype.isOne  = function () {
  return this.s === 1 && this.p === 1 && this.q === 1 && this.r === 1;
};
Surd.prototype.isRational = function () { return this.r === 1; };
Surd.prototype.isPureSqrt  = function () {
  return this.p === 1 && this.q === 1 && this.r > 1;
};

Surd.prototype.clone = function () {
  return new Surd(this.s, this.p, this.q, this.r);
};

Surd.prototype.equals = function (other) {
  if (!(other instanceof Surd)) return false;
  return this.s === other.s && this.p === other.p &&
         this.q === other.q && this.r === other.r;
};

Surd.prototype.hashKey = function () {
  // Two Surds can be added/subtracted only when they have
  // the same sqrt(r)/q part.
  return this.r + '/' + this.q;
};

// ---- Arithmetic ---------------------------------------------------

/**
 * Multiply two surds.
 */
Surd.mul = function (a, b) {
  if (a.isZero() || b.isZero()) return Surd.ZERO;
  return new Surd(
    a.s * b.s,
    a.p * b.p,
    a.q * b.q,
    a.r * b.r
  );
};

/**
 * Multiply a surd by a rational num/den.
 */
Surd.scale = function (a, num, den) {
  if (a.isZero() || num === 0) return Surd.ZERO;
  if (den === undefined || den === 0) den = 1;
  return new Surd(a.s * (num < 0 ? -1 : 1), a.p * Math.abs(num), a.q * den, a.r);
};

/**
 * Divide a by b.
 */
Surd.div = function (a, b) {
  if (a.isZero()) return Surd.ZERO;
  // a / b = (s1*p1*sqrt(r1)/q1) / (s2*p2*sqrt(r2)/q2)
  //       = (s1/s2) * (p1*q2) * sqrt(r1/r2) / (q1*p2)
  // sqrt(r1/r2) = sqrt(r1*r2) / r2
  // Actually more carefully:
  // a/b = [s1*p1*sqrt(r1)/q1] / [s2*p2*sqrt(r2)/q2]
  //      = (s1/s2) * (p1*q2) * sqrt(r1) / (q1*p2) / sqrt(r2)
  //      = (s1/s2) * (p1*q2) * sqrt(r1*r2) / (q1*p2) / r2
  // Wait, let me be more careful:
  // sqrt(r1)/sqrt(r2) = sqrt(r1/r2) = sqrt(r1*r2)/r2
  // Actually: sqrt(r1)/sqrt(r2) = sqrt(r1*r2)/r2 when rationalizing
  // More precisely: 1/sqrt(r2) = sqrt(r2)/r2
  // So sqrt(r1)/sqrt(r2) = sqrt(r1)*sqrt(r2)/r2 = sqrt(r1*r2)/r2
  // 
  // So a/b = (s1/s2) * (p1*q2) * sqrt(r1*r2) / (r2 * q1 * p2)
  
  return new Surd(
    a.s * b.s,
    a.p * b.q,
    a.q * b.p * b.r,   // denominator includes r2 from rationalising
    a.r * b.r
  );
};

/**
 * Square a surd (rational result).
 */
Surd.sq = function (a) {
  if (a.isZero()) return Surd.ZERO;
  return new Surd(1, a.p * a.p * a.r, a.q * a.q, 1);
};

/**
 * sqrt of a surd (only if perfect square).
 * Returns null if not a perfect square.
 */
Surd.sqrt = function (a) {
  if (a.isZero()) return Surd.ZERO;
  if (a.r !== 1) return null; // sqrt of sqrt not supported
  // Need a.p and a.q to be perfect squares
  var sp = Math.round(Math.sqrt(a.p));
  var sq = Math.round(Math.sqrt(a.q));
  if (sp * sp !== a.p || sq * sq !== a.q) return null;
  return new Surd(a.s > 0 ? 1 : -1, sp, sq, 1);
};

// ---- String output -----------------------------------------------

/**
 * Algebrite-compatible string.
 * Examples: "1", "-1/2", "sqrt(3)/3", "-3*sqrt(5)/7"
 */
Surd.prototype.toString = function () {
  if (this.isZero()) return '0';
  var prefix = this.s < 0 ? '-' : '';
  
  if (this.r === 1) {
    // Pure rational
    if (this.q === 1) return prefix + this.p;
    if (this.p === 1) return prefix + '1/' + this.q;
    return prefix + this.p + '/' + this.q;
  }
  
  // Has sqrt factor
  var sqrtStr = 'sqrt(' + this.r + ')';
  if (this.p === 1) {
    // sqrt(r) / q  or  just sqrt(r)
    if (this.q === 1) return prefix + sqrtStr;
    return prefix + sqrtStr + '/' + this.q;
  }
  // p * sqrt(r)  or  p * sqrt(r) / q
  if (this.q === 1) return prefix + this.p + '*' + sqrtStr;
  return prefix + this.p + '*' + sqrtStr + '/' + this.q;
};

/**
 * LaTeX string.
 * Examples: "1", "-\frac{1}{2}", "\frac{\sqrt{3}}{3}", "-\frac{3\sqrt{5}}{7}"
 */
Surd.prototype.toLatex = function () {
  if (this.isZero()) return '0';
  
  var signLatex = this.s < 0 ? '-' : '';
  var numStr, denStr;
  
  if (this.r === 1) {
    // Pure rational
    if (this.q === 1) return signLatex + this.p;
    numStr = String(this.p);
    denStr = String(this.q);
  } else {
    var sqrtStr = '\\sqrt{' + this.r + '}';
    if (this.p === 1) numStr = sqrtStr;
    else numStr = this.p + sqrtStr;
    
    if (this.q === 1) return signLatex + numStr;
    denStr = String(this.q);
  }
  
  return signLatex + '\\frac{' + numStr + '}{' + denStr + '}';
};

// ---- Sum of surds (ordered collection) ---------------------------

/**
 * SurdSum — ordered sum of Surd terms.
 * Groups like terms by radicand r only, combining rational
 * coefficients via exact common-denominator arithmetic.
 */
function SurdSum() {
  this._groups = {};   // r -> {num: int, den: int>0}  (signed rational coeff)
  this._order = [];    // r insertion order
}

SurdSum.prototype.add = function (surd) {
  if (surd.isZero()) return this;
  var r = surd.r;
  // term = surd.s * surd.p / surd.q  (rational coefficient for sqrt(r))
  var tnum = surd.s * surd.p;
  var tden = surd.q;

  if (this._groups[r] === undefined) {
    this._groups[r] = { num: tnum, den: tden };
    this._order.push(r);
  } else {
    var g = this._groups[r];
    // tnum/tden + g.num/g.den
    var newNum = tnum * g.den + g.num * tden;
    var newDen = tden * g.den;
    if (newNum === 0) {
      delete this._groups[r];
      var idx = this._order.indexOf(r);
      if (idx !== -1) this._order.splice(idx, 1);
    } else {
      var gcd = Surd._gcd(Math.abs(newNum), newDen);
      g.num = newNum / gcd;
      g.den = newDen / gcd;
    }
  }
  return this;
};

SurdSum.prototype.addAll = function (surds) {
  for (var i = 0; i < surds.length; i++) this.add(surds[i]);
  return this;
};

SurdSum.prototype.terms = function () {
  var result = [];
  for (var i = 0; i < this._order.length; i++) {
    var r = this._order[i];
    var g = this._groups[r];
    // Build Surd from reduced coefficient g.num/g.den with sqrt(r)
    var sign = g.num < 0 ? -1 : 1;
    result.push(new Surd(sign, Math.abs(g.num), g.den, r));
  }
  return result;
};

SurdSum.prototype.isEmpty = function () {
  return this._order.length === 0;
};

SurdSum.prototype.isZero = function () {
  return this._order.length === 0;
};

/**
 * Helper: format a single group as Algebrite string.
 * (r = radicand, num = signed integer numerator, den = positive denominator)
 */
SurdSum._groupToString = function (r, num, den) {
  num = Number(num); den = Number(den);
  if (num === 0) return null;
  var sign = num < 0 ? '-' : '';
  var absNum = Math.abs(num);

  if (r === 1) {
    if (den === 1) return sign + absNum;
    if (absNum === 1) return sign + '1/' + den;
    return sign + absNum + '/' + den;
  }
  var sqrtStr = 'sqrt(' + r + ')';
  if (den === 1) {
    if (absNum === 1) return sign + sqrtStr;
    return sign + absNum + '*' + sqrtStr;
  }
  if (absNum === 1) return sign + sqrtStr + '/' + den;
  return sign + absNum + '*' + sqrtStr + '/' + den;
};

/**
 * Helper: format a single group as LaTeX string.
 */
SurdSum._groupToLatex = function (r, num, den) {
  num = Number(num); den = Number(den);
  if (num === 0) return null;
  var sign = num < 0 ? '-' : '';
  var absNum = Math.abs(num);

  if (r === 1) {
    if (den === 1) return sign + absNum;
    return sign + '\\frac{' + absNum + '}{' + den + '}';
  }
  var sqrtStr = '\\sqrt{' + r + '}';
  if (den === 1) {
    if (absNum === 1) return sign + sqrtStr;
    return sign + absNum + sqrtStr;
  }
  if (absNum === 1) return sign + '\\frac{' + sqrtStr + '}{' + den + '}';
  return sign + '\\frac{' + absNum + sqrtStr + '}{' + den + '}';
};

/**
 * Algebrite-compatible string of the sum.
 */
SurdSum.prototype.toString = function () {
  if (this.isEmpty()) return '0';
  var parts = [];
  for (var i = 0; i < this._order.length; i++) {
    var r = this._order[i];
    var g = this._groups[r];
    var s = SurdSum._groupToString(r, g.num, g.den);
    if (s !== null) parts.push(s);
  }
  if (parts.length === 0) return '0';
  var joined = parts.join('+').replace(/\+-/g, '-');
  // If first character is '+', strip it (shouldn't happen with proper sign handling)
  return joined[0] === '+' ? joined.substring(1) : joined;
};

/**
 * LaTeX string of the sum.
 */
SurdSum.prototype.toLatex = function () {
  if (this.isEmpty()) return '0';
  var parts = [];
  for (var i = 0; i < this._order.length; i++) {
    var r = this._order[i];
    var g = this._groups[r];
    var s = SurdSum._groupToLatex(r, g.num, g.den);
    if (s !== null) parts.push(s);
  }
  if (parts.length === 0) return '0';
  return parts.join(' + ').replace(/\+ -/g, '- ');
};

// ---- Convenience constructors ------------------------------------

/** Create from a rational numerator/denominator. */
Surd.rational = function (num, den) {
  return new Surd(num < 0 ? -1 : 1, Math.abs(num), den, 1);
};

/** Create from a given radicand r (sqrt(r) with p=q=1). */
Surd.fromRadicand = function (r) {
  return new Surd(1, 1, 1, r);
};

// ---- Node/browser export -----------------------------------------

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { Surd: Surd, SurdSum: SurdSum };
}

