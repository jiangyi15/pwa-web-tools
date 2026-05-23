# Partial Wave Analysis Tools

A polished, modern static site with interactive calculators for quantum angular momentum coupling and particle decay kinematics. Built with pure HTML, CSS, and JavaScript — no frameworks, no build step, no server required.

## Tools

| Tool | Description |
|------|-------------|
| **Angular Formula Calculator** | Compute angular formulas for particle decay chains using helicity formalism. Supports cascade decays with exact symbolic computation of Wigner D-functions and CG coefficients. |
| **LS → Helicity Converter** | Express helicity amplitude coefficients H<sub>λ<sub>B</sub>,λ<sub>C</sub></sub> in terms of LS coupling amplitudes g<sub>l,s</sub> with exact Clebsch-Gordan coefficients. Includes parity selection (P<sub>A</sub>·P<sub>B</sub>·P<sub>C</sub> = (−1)<sup>l</sup>). |
| **Clebsch-Gordan Coefficient Calculator** | Compute ⟨j₁ m₁ j₂ m₂ | J M⟩ coupling coefficients using the Racah formula with exact symbolic surd arithmetic. |
| **Wigner d-Matrix Calculator** | Compute Wigner small d-matrix elements d<sup>j</sup><sub>m₁,m₂</sub>(β) using exact symbolic computation. Supports single elements and full matrix display. |

## Setup

No third-party vendor dependencies needed. All computation uses pure JavaScript with BigInt arithmetic.

For LaTeX rendering, MathJax is loaded from CDN:

```bash
# No setup required — MathJax loads from jsDelivr CDN at runtime
# Works offline upon first cache
```

## Usage

### Local Development

Simply open `index.html` in your browser — no server needed:

```bash
# Clone or download the repo
cd pwa-tools

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

All calculators use a custom `Surd` arithmetic library for exact symbolic computation:

- **Surd** — exact representation of `sign × p × √r / q` where p, q, r are non-negative integers
- **SurdSum** — sum of Surd terms for representing combined coefficients
- Pure BigInt arithmetic with gcd reduction — no floating-point error
- CG coefficients computed via the Racah formula using exact integer factorials
- Wigner d-matrix elements computed via half-angle Fourier expansion with exact binomial coefficients

### Helicity Amplitude Pipeline

The Angular Formula calculator combines:

1. **CG coefficients** from `cg.js` — exact Clebsch-Gordan coefficients
2. **Wigner d-matrix** from `wigner-d.js` — half-angle expansions of d<sup>j</sup><sub>m₁,m₂</sub>(β)
3. **LS coupling** from `get-angle.js` — conversion between LS and helicity bases with parity selection

All coefficients are combined via `Surd` arithmetic, producing LaTeX output directly from the internal SurdSum representation — no string parsing, no regex post-processing.

## Project Structure

```
pwa-tools/
├── index.html                  # Landing page with tool cards
├── .nojekyll                   # Disable Jekyll on GitHub Pages
├── assets/
│   ├── css/
│   │   └── style.css           # All styles (dark theme, responsive)
│   └── js/
│       ├── surd.js             # Exact surd arithmetic (Surd, SurdSum)
│       ├── cg.js               # Clebsch-Gordan coefficients
│       ├── wigner-d.js         # Wigner d-matrix computation
│       └── get-angle.js        # Helicity amplitude formulas
├── tools/
│   ├── get-angle-calculator.html   # Angular Formula Calculator
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

Requires BigInt support (ES2020+).

## License

MIT License — feel free to use, modify, and distribute.

## Contributing

Contributions welcome! Ideas for future tools:

- Wigner 3-j, 6-j, 9-j symbol calculators
- Racah coefficients
- Angular momentum eigenvalue problems
- Spherical tensor operators
- Decay width / Dalitz plot analysis

---

Built with care for the hadron spectroscopy community. ⚛️
