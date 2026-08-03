> **AI-Generated Project** — This project was developed through an interactive AI-assisted workflow using Opencode (Claude-based agent), iterating on physics tooling for hadron spectroscopy.
> 
> **Live Site**: [https://jiangyi15.github.io/pwa-web-tools](https://jiangyi15.github.io/pwa-web-tools)

# Partial Wave Analysis Tools

A polished, modern static site with interactive calculators for quantum angular momentum coupling and particle decay kinematics. Built with pure HTML, CSS, and JavaScript — no frameworks, no build step, no server required.

## Tools

| Tool | Description |
|------|-------------|
| **Angular Formula Calculator** | Compute angular formulas for particle decay chains using helicity formalism. Supports cascade decays with exact symbolic computation of Wigner D-functions and CG coefficients. |
| **Angular Distribution Calculator** | Compute the full angular distribution I(θ,φ) = \|T\|² for a decay chain, decomposed into per-LS intensities and 2·Re[T<sub>ls₁</sub>·T<sub>ls₂</sub>*] interference terms (real and imaginary parts) with exact symbolic output. |
| **Angular Projection Tool** | 1D projections of the angular distribution onto each angle variable (θ or φ), unit-area normalized to a proper probability density, with formulas and plots. |
| **Decay Angles 3D Viewer** | Interactive 3D visualization of particle decay angles in the helicity formalism (Three.js). |
| **LS → Helicity Converter** | Express helicity amplitude coefficients H<sub>λ<sub>B</sub>,λ<sub>C</sub></sub> in terms of LS coupling amplitudes g<sub>l,s</sub> with exact Clebsch-Gordan coefficients. Includes parity selection (P<sub>A</sub>·P<sub>B</sub>·P<sub>C</sub> = (−1)<sup>l</sup>). |
| **Clebsch-Gordan Coefficient Calculator** | Compute ⟨j₁ m₁ j₂ m₂ | J M⟩ coupling coefficients using the Racah formula with exact BigInt/surd arithmetic. |
| **Wigner d-Matrix Calculator** | Compute Wigner small d-matrix elements d<sup>j</sup><sub>m₁,m₂</sub>(β) using exact BigInt/surd computation. Supports single elements and full matrix display. |

## Setup

No third-party CAS or vendored libraries needed — all computation uses pure JavaScript with BigInt/surd arithmetic (the former Algebrite dependency has been fully removed).

Two CDN libraries load at runtime:

- **MathJax** — LaTeX rendering (all pages)
- **Three.js** — 3D scene in the Decay Angles 3D viewer (that page only)

```bash
# No setup required — CDN libraries load from jsDelivr at runtime
# Works offline upon first cache
```

## Usage

### Local Development

Simply open `index.html` in your browser — no server needed:

```bash
# Clone or download the repo
cd pwa-web-tools

# Open in browser
open index.html  # macOS
xdg-open index.html  # Linux
start index.html  # Windows
```

Or use a simple HTTP server (required for MathJax CDN to load):

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

## Technical Details

### Symbolic Computation with Surd.js

All calculators use a custom `Surd` arithmetic library for exact symbolic computation — no computer algebra system required:

- **Surd** — exact representation of `sign × p × √r / q` where p, q, r are non-negative integers
- **SurdSum** — sum of Surd terms for representing combined coefficients, grouped by radicand with exact common-denominator arithmetic
- Pure BigInt arithmetic with gcd reduction — no floating-point error
- CG coefficients computed via the Racah formula using exact BigInt factorials (`cg.js`)
- Wigner d-matrix elements computed via exact BigInt weights and half-angle Fourier expansion with exact binomial coefficients (`wigner-d.js`)
- Trigonometric polynomial arithmetic for the amplitude pipeline (`trig-poly.js`) and angular-distribution computation (`angular-expression.js`)

### Angular Projection Normalization

The **Angular Projection** tool computes 1D projections of the angular distribution I(θ,φ) = |T|² onto each individual angle variable (θ or φ). The projection is the marginal distribution:

- **θ projection** at cos θ: I(cos θ) = ∫ I(θ, φ) dφ (integrated over φ, plotted vs cos θ)
- **φ projection** at φ: I(φ) = ∫ I(θ, φ) d(cos θ) (integrated over cos θ, plotted vs φ)

Each projection is then normalized by the **computed total integral** ∫ I(θ,φ) d(cos θ) dφ, so that every projection integrates to 1 — i.e. it becomes the probability density of that angle variable:

- ∫ I<sub>proj</sub>(cos θ) d(cos θ) = 1
- ∫ I<sub>proj</sub>(φ) dφ = 1

For a flat |T|² = 1 the marginal equals the measure of the perpendicular variables divided by the total angular measure 2<sup>N<sub>θ</sub></sup> × (2π)<sup>N<sub>φ</sub></sup>, giving the expected flat densities:

- I<sub>proj</sub>(cos θ) = 1/2  (cos θ range = 2)
- I<sub>proj</sub>(φ) = 1/(2π)  (φ range = 2π)

### Helicity Amplitude Pipeline

The Angular Distribution and Angular Projection tools combine:

1. **CG coefficients** from `cg.js` — exact Clebsch-Gordan coefficients
2. **Wigner d-matrix** from `wigner-d.js` — exact BigInt weights with half-angle expansions of d<sup>j</sup><sub>m₁,m₂</sub>(β)
3. **LS coupling** from `get-angle.js` — conversion between LS and helicity bases with parity selection
4. **Trigonometric polynomial arithmetic** from `trig-poly.js` — exact multiplication/expansion in the {cos(kθ/2), sin(kφ/2)} Fourier basis, with real/imaginary parts kept separate (2·Re vs 2·Im interference terms)
5. **Angular distributions** from `angular-expression.js` — per-LS intensity maps, interference maps, and LaTeX rendering

All coefficients are combined via `Surd` arithmetic, producing LaTeX output directly from the internal SurdSum representation — no string parsing, no regex post-processing.

## Project Structure

```
pwa-web-tools/
├── index.html                  # Landing page with tool cards
├── .nojekyll                   # Disable Jekyll on GitHub Pages
├── assets/
│   ├── css/
│   │   └── style.css           # All styles (dark theme, responsive)
│   └── js/
│       ├── surd.js             # Exact surd arithmetic (Surd, SurdSum)
│       ├── cg.js               # Clebsch-Gordan coefficients (BigInt Racah)
│       ├── wigner-d.js         # Wigner d-matrix computation (BigInt)
│       ├── trig-poly.js        # Trigonometric polynomial arithmetic
│       ├── get-angle.js        # Helicity amplitude formulas (LS coupling)
│       ├── angular-expression.js  # Angular distributions + interference maps
│       └── main.js             # Shared UI bootstrap/utilities
├── tools/
│   ├── get-angle-calculator.html   # Angular Formula Calculator
│   ├── angular-distribution.html   # Angular Distribution Calculator
│   ├── angular-projection.html     # Angular Projection Tool
│   ├── decay-angles-3d.html        # Decay Angles 3D Viewer
│   ├── ls-to-helicity.html         # LS → Helicity Converter
│   ├── cg-calculator.html          # CG Coefficient Calculator
│   └── wigner-d-calculator.html    # Wigner d-Matrix Calculator
└── README.md                   # This file
```

## Browser Support

- Chrome/Chromium 90+
- Firefox 90+
- Safari 14+
- Edge 90+

Requires BigInt support (ES2020+). The **Decay Angles 3D Viewer** additionally uses an ES module import map for Three.js, which requires Chrome 89+ / Firefox 108+ / Safari 16.4+ / Edge 89+.

## License

MIT License — feel free to use, modify, and distribute.

## Contributing

Contributions welcome!

Built with care for the hadron spectroscopy community. ⚛️
