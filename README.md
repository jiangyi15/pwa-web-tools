# Partial Wave Analysis Tools

A polished, modern static site with interactive calculators for quantum angular momentum coupling. Built with pure HTML, CSS, and JavaScript — no frameworks, no build step, no server required.

## Features

- **Clebsch-Gordan Coefficient Calculator**: Compute ⟨j₁ m₁ j₂ m₂ | J M⟩ coupling coefficients using the Racah formula
- **Symbolic Computation**: Powered by [Algebrite](http://algebrite.org/) for exact symbolic results
- **Half-Integer Support**: Accepts fractions like "1/2", "3/2" and decimals like "0.5", "1.5"
- **Symbolic Output**: Results displayed in exact form (e.g., `1/sqrt(2)`, `sqrt(2/3)`) plus decimal approximation
- **Selection Rule Validation**: Automatic checking of quantum number constraints
- **Offline-First**: Works completely offline via `file://` protocol
- **Responsive Design**: Mobile-first, works beautifully on all screen sizes

## Tools

| Tool | Status | Description |
|------|--------|-------------|
| **Clebsch-Gordan Calculator** | ✅ Available | Angular momentum coupling coefficients |

## Setup

Third-party libraries (Algebrite, KaTeX + fonts) are **not** committed to the
repo. Fetch them once with:

```bash
# Algebrite (symbolic math)
mkdir -p assets/vendor
curl -L -o assets/vendor/algebrite.bundle-for-browser.js \
  https://cdn.jsdelivr.net/npm/algebrite@1.4.0/dist/algebrite.bundle-for-browser.js

# KaTeX (LaTeX rendering)
mkdir -p assets/vendor/katex/fonts
curl -L -o assets/vendor/katex/katex.min.js  https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js
curl -L -o assets/vendor/katex/katex.min.css https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css
for f in KaTeX_AMS-Regular KaTeX_Caligraphic-Bold KaTeX_Caligraphic-Regular \
         KaTeX_Fraktur-Bold KaTeX_Fraktur-Regular KaTeX_Main-Bold KaTeX_Main-BoldItalic \
         KaTeX_Main-Italic KaTeX_Main-Regular KaTeX_Math-BoldItalic KaTeX_Math-Italic \
         KaTeX_SansSerif-Bold KaTeX_SansSerif-Italic KaTeX_SansSerif-Regular \
         KaTeX_Script-Regular KaTeX_Size1-Regular KaTeX_Size2-Regular KaTeX_Size3-Regular \
         KaTeX_Size4-Regular KaTeX_Typewriter-Regular; do
  curl -sL -o "assets/vendor/katex/fonts/$f.woff2" \
    "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/fonts/$f.woff2"
done
```

## Usage

### Local Development

Simply open `index.html` in your browser — no server needed:

```bash
# Clone or download the repo
cd physics-toolkit

# Open in browser
open index.html  # macOS
xdg-open index.html  # Linux
start index.html  # Windows
```

Or use a simple HTTP server:

```bash
# Python 3
python -m http.server 8000

# Node.js (npx)
npx serve .

# Then open http://localhost:8000
```

### GitHub Pages Deployment

1. Push this repository to GitHub
2. Go to **Settings** → **Pages**
3. Under **Source**, select **Deploy from a branch**
4. Select the **main** branch and **/ (root)** folder
5. Click **Save**
6. Your site will be live at `https://<username>.github.io/<repo-name>/`

The `.nojekyll` file ensures GitHub Pages serves the `assets/` directory without Jekyll processing.

## Clebsch-Gordan Calculator

The calculator implements the Racah formula:

```
⟨j₁ m₁ j₂ m₂ | J M⟩ = δ(m₁+m₂, M) · √((2J+1) · Δ(j₁,j₂,J)) ·
  √((j₁+m₁)!(j₁-m₁)!(j₂+m₂)!(j₂-m₂)!(J+M)!(J-M)!) ·
  Σₖ (-1)ᵏ / [k! (j₁+j₂-J-k)! (j₁-m₁-k)! (j₂+m₂-k)! (J-j₂+m₁+k)! (J-j₁-m₂+k)!]
```

where Δ(a,b,c) = (a+b-c)!(a-b+c)!(-a+b+c)! / (a+b+c+1)! is the triangle coefficient.

### Selection Rules

The coefficient is zero if any of these conditions are violated:

- **Projection conservation**: m₁ + m₂ = M
- **Triangle inequality**: |j₁ − j₂| ≤ J ≤ j₁ + j₂
- **Magnitude constraints**: |m₁| ≤ j₁, |m₂| ≤ j₂, |M| ≤ J
- **Quantum consistency**: j±m must be non-negative integers

### Input Formats

The calculator accepts several input formats:

| Format | Examples | Notes |
|--------|----------|-------|
| Fraction | `1/2`, `3/2`, `-1/2` | Exact half-integers |
| Decimal | `0.5`, `1.5`, `-0.5` | Converted to fractions internally |
| Integer | `0`, `1`, `2`, `-1` | Standard integer values |

### Examples

Try these well-known coefficients:

- ⟨½ ½ ½ −½ | 1 0⟩ = 1/√2 ≈ 0.707...
- ⟨1 0 1 0 | 2 0⟩ = √(2/3) ≈ 0.816...
- ⟨1 1 1 −1 | 0 0⟩ = 1/√3 ≈ 0.577...
- ⟨1 0 1 0 | 0 0⟩ = −1/√3 ≈ −0.577...

## Technical Details

### Symbolic Computation with Algebrite

The calculator uses [Algebrite](http://algebrite.org/) — a pure JavaScript computer algebra system — for exact symbolic computation:

1. Parse inputs to exact fractions (integers or half-integers)
2. Build the Racah formula expression as an Algebrite-compatible string
3. Compute the sum over k for all valid terms
4. Simplify the result symbolically using Algebrite's `simplify()`
5. Display the exact symbolic form and decimal approximation via `float()`

This approach provides mathematically exact results (e.g., `1/sqrt(2)`) rather than floating-point approximations, while handling arbitrarily complex expressions.

### Console Tests

Open the browser console on the calculator page to see automatic sanity checks:

```javascript
runSanityChecks();  // Re-run tests
```

## Project Structure

```
physics-toolkit/
├── index.html              # Landing page
├── tools/
│   └── cg-calculator.html  # CG coefficient calculator
├── assets/
│   ├── css/
│   │   └── style.css       # All styles (dark theme)
│   ├── js/
│   │   ├── cg.js           # CG coefficient math (Algebrite)
│   │   └── main.js         # UI helpers
│   └── vendor/
│       └── algebrite.bundle-for-browser.js  # Algebrite CAS library
├── .nojekyll               # Disable Jekyll on GitHub Pages
└── README.md               # This file
```

## Browser Support

- Chrome/Chromium 90+
- Firefox 90+
- Safari 14+
- Edge 90+

Requires BigInt support (ES2020+).

## License

MIT License — feel free to use, modify, and distribute.

## Acknowledgments

- **[Algebrite](http://algebrite.org/)** — Symbolic computation powered by this excellent pure-JavaScript computer algebra system. Algebrite is released under the MIT License.

## Contributing

Contributions welcome! Ideas for future tools:

- Wigner 3-j, 6-j, 9-j symbol calculators
- Racah coefficients
- Angular momentum eigenvalue problems
- Spherical tensor operators

---

Built with care for the physics community. ⚛️
