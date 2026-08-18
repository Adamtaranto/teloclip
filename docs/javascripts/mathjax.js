/* MathJax 3 configuration for pymdownx.arithmatex (generic mode).
 *
 * Arithmatex wraps maths in elements with the `arithmatex` class, so MathJax
 * is told to process exactly those and nothing else. The site uses instant
 * navigation, so re-typeset on each page load event when the Material
 * `document$` observable is available; otherwise the initial typeset on
 * script load covers the classic full-page navigation case.
 */
window.MathJax = {
  tex: {
    inlineMath: [['\\(', '\\)']],
    displayMath: [['\\[', '\\]']],
    processEscapes: true,
    processEnvironments: true,
  },
  options: {
    ignoreHtmlClass: '.*|',
    processHtmlClass: 'arithmatex',
  },
};

if (window.document$ && window.document$.subscribe) {
  window.document$.subscribe(function () {
    if (window.MathJax.typesetPromise) {
      window.MathJax.typesetClear();
      window.MathJax.typesetPromise();
    }
  });
}
