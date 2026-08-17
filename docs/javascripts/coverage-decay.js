/* Interactive terminal coverage-decay explorer.
 *
 * Self-contained: evaluates the analytic model c(x) from the coverage-decay
 * guide page in closed form (normal CDF via the Abramowitz–Stegun 7.1.26 erf
 * approximation, |error| < 1.5e-7) and renders an inline SVG, in the same
 * hand-written-SVG idiom as the teloclip extend HTML report. Colours come
 * from CSS custom properties defined in stylesheets/coverage-decay.css, so
 * the light/dark theme toggle needs no JavaScript here.
 *
 * The maths mirrors teloclip/model/decay.py exactly; if you change one,
 * change the other (tests spot-check parity between the two).
 */
(function () {
  'use strict';

  // --- model -------------------------------------------------------------

  function erf(z) {
    // Abramowitz & Stegun 7.1.26, odd symmetry.
    var sign = z < 0 ? -1 : 1;
    var az = Math.abs(z);
    var t = 1 / (1 + 0.3275911 * az);
    var poly =
      t *
      (0.254829592 +
        t *
          (-0.284496736 +
            t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    return sign * (1 - poly * Math.exp(-az * az));
  }

  function phi(z) {
    return 0.5 * (1 + erf(z / Math.SQRT2));
  }

  // Lognormal helpers, parameterised by real-space mean m and sd s.
  function logParams(m, s) {
    var sigma2 = Math.log(1 + (s * s) / (m * m));
    return { mu: Math.log(m) - sigma2 / 2, sigma: Math.sqrt(sigma2) };
  }

  function survival(u, p) {
    if (u <= 0) return 1;
    return 1 - phi((Math.log(u) - p.mu) / p.sigma);
  }

  function partialExpectation(a, p, mean) {
    // E[min(a, l)] for lognormal l.
    if (a <= 0) return 0;
    var z = (Math.log(a) - p.mu) / p.sigma;
    return mean * phi(z - p.sigma) + a * (1 - phi(z));
  }

  // Relative coverage c(x) for the current parameter set. Returns a function
  // of x so per-frame slider updates precompute the constants once.
  function makeCurve(params) {
    if (params.dist === 'fixed') {
      var L = params.mean;
      if (params.cutoff > L) return null; // size selection removes everything
      return function (x) {
        return Math.min(x + 1, L) / L;
      };
    }
    var p = logParams(params.mean, params.sd);
    var mean = params.mean;
    var Lmin = params.cutoff;
    if (Lmin <= 0) {
      return function (x) {
        return partialExpectation(x + 1, p, mean) / mean;
      };
    }
    var survMin = survival(Lmin, p);
    if (survMin < 1e-9) return null;
    var capMin = partialExpectation(Lmin, p, mean);
    var truncMean = Lmin + (mean - capMin) / survMin;
    return function (x) {
      var a = x + 1;
      var numer =
        a <= Lmin
          ? a
          : Lmin + (partialExpectation(a, p, mean) - capMin) / survMin;
      return numer / truncMean;
    };
  }

  // --- rendering ---------------------------------------------------------

  var W = 880;
  var H = 300;
  var PAD_L = 52;
  var PAD_R = 16;
  var PAD_T = 18;
  var PAD_B = 44;
  var N_POINTS = 240;
  var SVG_NS = 'http://www.w3.org/2000/svg';

  function el(tag, attrs, text) {
    var node = document.createElementNS(SVG_NS, tag);
    for (var k in attrs) node.setAttribute(k, attrs[k]);
    if (text !== undefined) node.textContent = text;
    return node;
  }

  function niceStep(range, target) {
    var raw = range / target;
    var mag = Math.pow(10, Math.floor(Math.log10(raw)));
    var candidates = [1, 2, 5, 10];
    for (var i = 0; i < candidates.length; i++) {
      if (candidates[i] * mag >= raw) return candidates[i] * mag;
    }
    return 10 * mag;
  }

  function drawChart(svg, params) {
    while (svg.firstChild) svg.removeChild(svg.firstChild);

    var plotW = W - PAD_L - PAD_R;
    var plotH = H - PAD_T - PAD_B;
    var xMax = params.window;
    var yMax = params.depth * 1.08;
    var curve = makeCurve(params);

    function sx(x) {
      return PAD_L + (x / xMax) * plotW;
    }
    function sy(y) {
      return PAD_T + plotH - (y / yMax) * plotH;
    }

    // Gridlines and axis labels.
    var xStep = niceStep(xMax, 6);
    for (var gx = 0; gx <= xMax + 1; gx += xStep) {
      svg.appendChild(
        el('line', {
          x1: sx(gx), y1: PAD_T, x2: sx(gx), y2: PAD_T + plotH,
          class: 'cd-grid',
        })
      );
      svg.appendChild(
        el(
          'text',
          { x: sx(gx), y: H - PAD_B + 16, class: 'cd-tick', 'text-anchor': 'middle' },
          (gx / 1000).toFixed(0)
        )
      );
    }
    var yStep = niceStep(yMax, 5);
    for (var gy = 0; gy <= yMax; gy += yStep) {
      svg.appendChild(
        el('line', {
          x1: PAD_L, y1: sy(gy), x2: PAD_L + plotW, y2: sy(gy),
          class: 'cd-grid',
        })
      );
      svg.appendChild(
        el(
          'text',
          { x: PAD_L - 8, y: sy(gy) + 3, class: 'cd-tick', 'text-anchor': 'end' },
          gy.toFixed(0)
        )
      );
    }

    // Axes.
    svg.appendChild(
      el('line', { x1: PAD_L, y1: PAD_T + plotH, x2: PAD_L + plotW, y2: PAD_T + plotH, class: 'cd-axis' })
    );
    svg.appendChild(
      el('line', { x1: PAD_L, y1: PAD_T, x2: PAD_L, y2: PAD_T + plotH, class: 'cd-axis' })
    );
    svg.appendChild(
      el(
        'text',
        { x: PAD_L + plotW / 2, y: H - 8, class: 'cd-label', 'text-anchor': 'middle' },
        'Distance from contig end (kb)'
      )
    );
    var yLabel = el(
      'text',
      { x: 14, y: PAD_T + plotH / 2, class: 'cd-label', 'text-anchor': 'middle' },
      'Expected coverage (reads)'
    );
    yLabel.setAttribute('transform', 'rotate(-90 14 ' + (PAD_T + plotH / 2) + ')');
    svg.appendChild(yLabel);

    if (!curve) {
      svg.appendChild(
        el(
          'text',
          { x: PAD_L + plotW / 2, y: PAD_T + plotH / 2, class: 'cd-warning', 'text-anchor': 'middle' },
          'Cutoff removes every fragment — no coverage anywhere'
        )
      );
      return null;
    }

    // Interior depth reference line.
    svg.appendChild(
      el('line', {
        x1: PAD_L, y1: sy(params.depth), x2: PAD_L + plotW, y2: sy(params.depth),
        class: 'cd-ref',
      })
    );
    svg.appendChild(
      el(
        'text',
        { x: PAD_L + plotW - 4, y: sy(params.depth) - 5, class: 'cd-ref-label', 'text-anchor': 'end' },
        'interior depth ' + params.depth + '×'
      )
    );

    // The decay curve.
    var pts = [];
    for (var i = 0; i <= N_POINTS; i++) {
      var x = (xMax * i) / N_POINTS;
      pts.push(sx(x).toFixed(1) + ',' + sy(params.depth * curve(x)).toFixed(1));
    }
    svg.appendChild(
      el('polyline', { points: pts.join(' '), class: 'cd-curve', fill: 'none' })
    );
    return curve;
  }

  // --- controls ----------------------------------------------------------

  // Defaults mirror the PacBio HiFi-like library used in the static figures:
  // lognormal 15 +/- 3 kb with a strict 10 kb lower-bound cutoff.
  var SLIDERS = [
    { key: 'depth', label: 'Interior depth', min: 5, max: 100, step: 1, value: 30, unit: '×' },
    { key: 'mean', label: 'Mean fragment length', min: 1000, max: 30000, step: 500, value: 15000, unit: 'kb', kb: true },
    { key: 'sd', label: 'Fragment length sd', min: 500, max: 15000, step: 250, value: 3000, unit: 'kb', kb: true },
    { key: 'cutoff', label: 'Lower-bound cutoff (discard shorter)', min: 0, max: 25000, step: 500, value: 10000, unit: 'kb', kb: true },
    { key: 'window', label: 'x-axis range', min: 10000, max: 100000, step: 5000, value: 35000, unit: 'kb', kb: true },
  ];

  function fmt(spec, value) {
    return spec.kb ? (value / 1000).toFixed(1) + ' ' + spec.unit : value + ' ' + spec.unit;
  }

  function build(root) {
    var params = { dist: 'lognormal' };
    SLIDERS.forEach(function (s) {
      params[s.key] = s.value;
    });

    var controls = document.createElement('div');
    controls.className = 'cd-controls';

    // Distribution toggle.
    var distWrap = document.createElement('label');
    distWrap.className = 'cd-control';
    var distName = document.createElement('span');
    distName.className = 'cd-control-name';
    distName.textContent = 'Fragment lengths';
    var distSel = document.createElement('select');
    ['lognormal', 'fixed'].forEach(function (opt) {
      var o = document.createElement('option');
      o.value = opt;
      o.textContent = opt === 'lognormal' ? 'lognormal (mean, sd)' : 'fixed length';
      distSel.appendChild(o);
    });
    distWrap.appendChild(distName);
    distWrap.appendChild(distSel);
    controls.appendChild(distWrap);

    var inputs = {};
    var readouts = {};
    SLIDERS.forEach(function (spec) {
      var wrap = document.createElement('label');
      wrap.className = 'cd-control';
      var name = document.createElement('span');
      name.className = 'cd-control-name';
      name.textContent = spec.label;
      var value = document.createElement('span');
      value.className = 'cd-control-value';
      value.textContent = fmt(spec, spec.value);
      var input = document.createElement('input');
      input.type = 'range';
      input.min = spec.min;
      input.max = spec.max;
      input.step = spec.step;
      input.value = spec.value;
      wrap.appendChild(name);
      wrap.appendChild(input);
      wrap.appendChild(value);
      controls.appendChild(wrap);
      inputs[spec.key] = input;
      readouts[spec.key] = value;
    });

    var svg = el('svg', {
      viewBox: '0 0 ' + W + ' ' + H,
      class: 'cd-chart',
      role: 'img',
      'aria-label':
        'Expected coverage against distance from the contig end for the chosen parameters',
    });

    var summary = document.createElement('p');
    summary.className = 'cd-summary';

    function update() {
      params.dist = distSel.value;
      SLIDERS.forEach(function (spec) {
        params[spec.key] = Number(inputs[spec.key].value);
        readouts[spec.key].textContent = fmt(spec, params[spec.key]);
      });
      // sd is meaningless for a fixed length.
      inputs.sd.disabled = params.dist === 'fixed';
      inputs.sd.closest('.cd-control').classList.toggle('cd-disabled', params.dist === 'fixed');

      var curve = drawChart(svg, params);
      if (curve) {
        var atEnd = params.depth * curve(0);
        var at1kb = params.depth * curve(1000);
        summary.textContent =
          'Expected reads spanning the outermost base: ' + atEnd.toFixed(2) +
          ' • covering 1 kb in: ' + at1kb.toFixed(1) +
          ' • interior: ' + params.depth + '. ' +
          'Overhang candidates can only come from that first number.';
      } else {
        summary.textContent =
          'The size-selection cutoff is above every fragment length: nothing survives.';
      }
    }

    distSel.addEventListener('input', update);
    SLIDERS.forEach(function (spec) {
      inputs[spec.key].addEventListener('input', update);
    });

    root.appendChild(controls);
    root.appendChild(svg);
    root.appendChild(summary);
    update();
  }

  function init() {
    var root = document.getElementById('coverage-decay-app');
    if (!root || root.dataset.cdInit) return;
    root.dataset.cdInit = '1';
    build(root);
  }

  // The site uses instant navigation; Material exposes document$ for exactly
  // this re-init problem. Fall back to plain load events elsewhere.
  if (window.document$ && window.document$.subscribe) {
    window.document$.subscribe(init);
  } else {
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', init);
    } else {
      init();
    }
  }
})();
