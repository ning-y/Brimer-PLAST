(function () {
  'use strict';

  // ── State ─────────────────────────────────────────────────
  let genomePath = null;
  let gtfPath = null;

  const GENE_PLACEHOLDER  = 'e.g. LTBP4';
  const TRANS_PLACEHOLDER = 'e.g. ENST00000204005';

  // Active running requests: requestId -> rowIndex
  const runningMap = new Map();

  // ── DOM refs ──────────────────────────────────────────────
  const genomeBtn = document.getElementById('genome-btn');
  const genomePathEl = document.getElementById('genome-path');
  const gtfBtn = document.getElementById('gtf-btn');
  const gtfPathEl = document.getElementById('gtf-path');
  const targetsBody = document.getElementById('targets-body');
  const globalPdfArea = document.getElementById('global-pdf-area');
  const globalPdfLink = document.getElementById('global-pdf-link');

  // ── File pickers ─────────────────────────────────────────
  genomeBtn.addEventListener('click', async () => {
    const path = await window.api.selectFasta();
    if (path) { genomePath = path; genomePathEl.textContent = path; }
  });
  gtfBtn.addEventListener('click', async () => {
    const path = await window.api.selectGtf();
    if (path) { gtfPath = path; gtfPathEl.textContent = path; }
  });

  // ── Create a row ─────────────────────────────────────────
  function createRow() {
    const tr = document.createElement('tr');
    tr.className = 'target-row';

    // Target type dropdown
    const tdTarget = document.createElement('td');
    const sel = document.createElement('select');
    sel.className = 'target-type';
    const optG = document.createElement('option'); optG.value = 'gene'; optG.textContent = 'Gene';
    const optT = document.createElement('option'); optT.value = 'transcript'; optT.textContent = 'Transcript';
    sel.appendChild(optG);
    sel.appendChild(optT);
    tdTarget.appendChild(sel);

    // Name input
    const tdName = document.createElement('td');
    const inp = document.createElement('input');
    inp.type = 'text';
    inp.className = 'target-name';
    inp.placeholder = GENE_PLACEHOLDER;
    tdName.appendChild(inp);

    // Optional params (collapsible)
    const tdParams = document.createElement('td');
    tdParams.className = 'td-params';
    const toggleBtn = document.createElement('span');
    toggleBtn.className = 'params-toggle';
    toggleBtn.textContent = '▸';
    toggleBtn.title = 'Primer parameters';
    tdParams.appendChild(toggleBtn);

    const panel = document.createElement('div');
    panel.className = 'params-panel';
    panel.innerHTML = `
      <div class="params-grid">
        <label>Product min <input class="pp" data-key="product-min" type="number" value="80"></label>
        <label>Product max <input class="pp" data-key="product-max" type="number" value="200"></label>
        <label>Opt Tm (°C) <input class="pp" data-key="opt-tm" type="number" value="60" step="0.1"></label>
        <label>Min Tm (°C) <input class="pp" data-key="min-tm" type="number" value="57" step="0.1"></label>
        <label>Max Tm (°C) <input class="pp" data-key="max-tm" type="number" value="63" step="0.1"></label>
        <label>Opt length <input class="pp" data-key="opt-size" type="number" value="20"></label>
        <label>Min length <input class="pp" data-key="min-size" type="number" value="18"></label>
        <label>Max length <input class="pp" data-key="max-size" type="number" value="25"></label>
        <label>Min GC% <input class="pp" data-key="min-gc" type="number" value="40"></label>
        <label>Max GC% <input class="pp" data-key="max-gc" type="number" value="60"></label>
        <label>Pairs/chain <input class="pp" data-key="num-return" type="number" value="50"></label>
        <label>Max amplicon <input class="pp" data-key="max-amplicon" type="number" value="2000"></label>
      </div>
    `;
    tdParams.appendChild(panel);

    toggleBtn.addEventListener('click', () => {
      const open = panel.classList.toggle('open');
      toggleBtn.textContent = open ? '▾' : '▸';
      // Adjust table position if panel grows beyond viewport
      if (open) {
        const rect = panel.getBoundingClientRect();
        const rightEdge = rect.right;
        const vw = window.innerWidth;
        if (rightEdge > vw) {
          panel.style.left = 'auto';
          panel.style.right = '0';
        }
      }
    });
    // Close panel when clicking outside
    document.addEventListener('click', (e) => {
      if (!tdParams.contains(e.target)) {
        panel.classList.remove('open');
        toggleBtn.textContent = '▸';
      }
    });

    // Results cell
    const tdResults = document.createElement('td');
    tdResults.className = 'result-cell';
    const runBtn = document.createElement('button');
    runBtn.className = 'run-btn';
    runBtn.textContent = 'Run';
    tdResults.appendChild(runBtn);

    const statusEl = document.createElement('div');
    statusEl.className = 'status-text';
    statusEl.style.display = 'none';
    tdResults.appendChild(statusEl);

    const resultLink = document.createElement('span');
    resultLink.className = 'result-link';
    resultLink.style.display = 'none';
    tdResults.appendChild(resultLink);

    const errorEl = document.createElement('span');
    errorEl.className = 'error-msg';
    errorEl.style.display = 'none';
    errorEl.innerHTML = '<span class="error-tooltip"></span>';
    tdResults.appendChild(errorEl);

    // Assemble row
    tr.appendChild(tdTarget);
    tr.appendChild(tdName);
    tr.appendChild(tdParams);
    tr.appendChild(tdResults);
    targetsBody.appendChild(tr);

    // ── Placeholder switching ───────────────────────
    sel.addEventListener('change', () => {
      inp.placeholder = sel.value === 'gene' ? GENE_PLACEHOLDER : TRANS_PLACEHOLDER;
    });

    // ── Tab navigation ─────────────────────────────
    // Collect tabbable elements in this row
    function getTabbables() {
      const tabbables = [sel, inp];
      // If params panel is open, add the input fields
      if (panel.classList.contains('open')) {
        panel.querySelectorAll('.pp').forEach(el => tabbables.push(el));
      }
      tabbables.push(runBtn);
      return tabbables;
    }

    tr.addEventListener('keydown', (e) => {
      if (e.key !== 'Tab') return;
      e.preventDefault();
      const tabbables = getTabbables();
      const idx = tabbables.indexOf(e.target);
      if (idx === -1) return;

      if (e.shiftKey) {
        // Shift+Tab: previous
        if (idx > 0) {
          tabbables[idx - 1].focus();
        } else {
          // Go to previous row's last field
          const prevRow = tr.previousElementSibling;
          if (prevRow) {
            const prevTabbables = getRowTabbables(prevRow);
            if (prevTabbables.length > 0) prevTabbables[prevTabbables.length - 1].focus();
            else prevRow.querySelector('.target-type')?.focus();
          }
        }
      } else {
        // Tab: next
        if (idx < tabbables.length - 1) {
          tabbables[idx + 1].focus();
        } else {
          // Last field: go to next row or add new row
          const nextRow = tr.nextElementSibling;
          if (nextRow) {
            const nextTabbables = getRowTabbables(nextRow);
            if (nextTabbables.length > 0) nextTabbables[0].focus();
            else nextRow.querySelector('.target-type')?.focus();
          } else {
            // Add a new row and focus its first field
            const newRow = createRow();
            newRow.querySelector('.target-type').focus();
          }
        }
      }
    });

    // ── Run button ────────────────────────────────
    runBtn.addEventListener('click', async () => {
      if (runBtn.disabled) return;
      if (!genomePath) { showRowError(tr, 'Select a genome FASTA first.'); return; }
      if (!gtfPath) { showRowError(tr, 'Select a GTF annotation first.'); return; }
      const name = inp.value.trim();
      if (!name) { showRowError(tr, 'Enter a target name.'); return; }

      clearRowState(tr);
      runBtn.disabled = true;
      runBtn.textContent = 'Running…';
      statusEl.style.display = 'block';
      statusEl.textContent = 'Starting…';

      const targetType = sel.value;

      // Read optional params from panel
      const primerArgs = {};
      panel.querySelectorAll('.pp').forEach(el => {
        const key = el.dataset.key;
        const val = parseFloat(el.value);
        if (isNaN(val)) return;
        switch (key) {
          case 'product-min': break; // handled below
          case 'product-max': break;
          case 'num-return': primerArgs.PRIMER_NUM_RETURN = val; break;
          case 'opt-tm':     primerArgs.PRIMER_OPT_TM = val; break;
          case 'min-tm':     primerArgs.PRIMER_MIN_TM = val; break;
          case 'max-tm':     primerArgs.PRIMER_MAX_TM = val; break;
          case 'opt-size':   primerArgs.PRIMER_OPT_SIZE = val; break;
          case 'min-size':   primerArgs.PRIMER_MIN_SIZE = val; break;
          case 'max-size':   primerArgs.PRIMER_MAX_SIZE = val; break;
          case 'min-gc':     primerArgs.PRIMER_MIN_GC = val; break;
          case 'max-gc':     primerArgs.PRIMER_MAX_GC = val; break;
          case 'max-amplicon': break; // handled below
        }
      });
      const productMin = parseFloat(panel.querySelector('[data-key="product-min"]').value) || 80;
      const productMax = parseFloat(panel.querySelector('[data-key="product-max"]').value) || 200;
      primerArgs.PRIMER_PRODUCT_SIZE_RANGE = `${productMin}-${productMax}`;
      const maxAmplicon = parseFloat(panel.querySelector('[data-key="max-amplicon"]').value) || 2000;

      try {
        const { promise, requestId } = await window.api.runPipeline({
          genome: genomePath,
          annotations: gtfPath,
          target_key: name,
          target_type: targetType,
          primer_args: primerArgs,
          max_amplicon: maxAmplicon,
          pdf_output_dir: null,
        });

        runningMap.set(requestId, tr);

        const result = await promise;
        runningMap.delete(requestId);

        runBtn.disabled = false;
        runBtn.textContent = 'Run';
        statusEl.style.display = 'none';

        const count = result.filtered_pairs ? result.filtered_pairs.length : 0;
        resultLink.style.display = 'inline';
        resultLink.textContent = `${count} primer${count !== 1 ? 's' : ''} found`;

        if (result.pdf_path) {
          resultLink.dataset.pdfPath = result.pdf_path;
          resultLink.title = 'Click to open PDF report';
          resultLink.style.cursor = 'pointer';
          resultLink.addEventListener('click', () => {
            window.api.openPdf(result.pdf_path);
          });
        } else {
          resultLink.dataset.pdfPath = '';
          resultLink.title = '';
          resultLink.style.cursor = 'default';
        }
      } catch (err) {
        runningMap.delete(requestId);
        runBtn.disabled = false;
        runBtn.textContent = 'Run';
        statusEl.style.display = 'none';
        showRowError(tr, err.message);
      }
    });

    return tr;
  }

  // ── Helpers ──────────────────────────────────────────────
  function getRowTabbables(row) {
    const sel = row.querySelector('.target-type');
    const inp = row.querySelector('.target-name');
    const panel = row.querySelector('.params-panel');
    const runBtn = row.querySelector('.run-btn');
    const tabbables = [sel, inp].filter(Boolean);
    if (panel && panel.classList.contains('open')) {
      panel.querySelectorAll('.pp').forEach(el => tabbables.push(el));
    }
    if (runBtn) tabbables.push(runBtn);
    return tabbables;
  }

  function showRowError(tr, msg) {
    const errorEl = tr.querySelector('.error-msg');
    const tooltip = errorEl.querySelector('.error-tooltip');
    if (tooltip) tooltip.textContent = msg;
    errorEl.style.display = 'inline';
    // Auto-hide after 8 seconds
    setTimeout(() => { errorEl.style.display = 'none'; }, 8000);
  }

  function clearRowState(tr) {
    const statusEl = tr.querySelector('.status-text');
    const resultLink = tr.querySelector('.result-link');
    const errorEl = tr.querySelector('.error-msg');
    if (statusEl) { statusEl.style.display = 'none'; statusEl.textContent = ''; }
    if (resultLink) { resultLink.style.display = 'none'; resultLink.textContent = ''; resultLink.dataset.pdfPath = ''; }
    if (errorEl) { errorEl.style.display = 'none'; }
  }

  // ── Global progress/error listeners ──────────────────────
  window.api.onProgress((data) => {
    const row = runningMap.get(data._requestId);
    if (!row) return;
    const statusEl = row.querySelector('.status-text');
    if (statusEl) statusEl.textContent = data.message || '';
  });

  window.api.onError((data) => {
    const row = runningMap.get(data._requestId);
    if (!row) return;
    runningMap.delete(data._requestId);
    const runBtn = row.querySelector('.run-btn');
    const statusEl = row.querySelector('.status-text');
    if (runBtn) { runBtn.disabled = false; runBtn.textContent = 'Run'; }
    if (statusEl) statusEl.style.display = 'none';
    showRowError(row, data.message || 'Unknown error');
  });

  // ── Theme toggle ────────────────────────────────────────
  const htmlEl = document.documentElement;
  const themeToggle = document.getElementById('theme-toggle');

  function setTheme(theme) {
    htmlEl.setAttribute('data-theme', theme);
    themeToggle.textContent = theme === 'dark' ? '🌙' : '☀️';
    try { localStorage.setItem('brimer-theme', theme); } catch (e) { /* noop */ }
  }

  (function initTheme() {
    const saved = (() => { try { return localStorage.getItem('brimer-theme'); } catch (e) { return null; } })();
    if (saved === 'light' || saved === 'dark') {
      setTheme(saved);
    } else {
      const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
      setTheme(prefersDark ? 'dark' : 'light');
    }
  })();

  themeToggle.addEventListener('click', () => {
    const current = htmlEl.getAttribute('data-theme') || 'light';
    setTheme(current === 'dark' ? 'light' : 'dark');
  });

  // Listen for system preference changes
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', (e) => {
    const saved = (() => { try { return localStorage.getItem('brimer-theme'); } catch (e) { return null; } })();
    if (!saved) {
      setTheme(e.matches ? 'dark' : 'light');
    }
  });

  // ── Initial row ──────────────────────────────────────────
  createRow().querySelector('.target-type').focus();

  // ── Global PDF link ──────────────────────────────────────
  globalPdfLink.addEventListener('click', (e) => {
    e.preventDefault();
    const href = globalPdfLink.getAttribute('href');
    if (href && href !== '#') window.api.openPdf(href);
  });

})();