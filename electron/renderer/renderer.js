(function () {
  'use strict';

  try {
    // ── State ─────────────────────────────────────────────────
  let genomePath = (() => { try { return localStorage.getItem('brimer-genome-path'); } catch (e) { return null; } })();
  let gtfPath = (() => { try { return localStorage.getItem('brimer-gtf-path'); } catch (e) { return null; } })();

  const GENE_PLACEHOLDER  = 'e.g. LTBP4';
  const TRANS_PLACEHOLDER = 'e.g. ENST00000204005';

  // Active running requests: requestId -> rowIndex
  const runningMap = new Map();

  // Hardcoded defaults for highlight comparison
  const DEFAULTS = {
    'product-min': 80,
    'product-max': 200,
    'opt-tm': 60,
    'min-tm': 57,
    'max-tm': 63,
    'opt-size': 20,
    'min-size': 18,
    'max-size': 25,
    'min-gc': 40,
    'max-gc': 60,
    'num-return': 50,
    'max-amplicon': 2000,
  };

  // ── DOM refs ──────────────────────────────────────────────
  const genomeBtn = document.getElementById('genome-btn');
  const genomePathEl = document.getElementById('genome-path');
  const gtfBtn = document.getElementById('gtf-btn');
  const gtfPathEl = document.getElementById('gtf-path');
  const targetsBody = document.getElementById('targets-body');

  // Load saved paths into UI
  if (genomePath) genomePathEl.textContent = genomePath;
  if (gtfPath) gtfPathEl.textContent = gtfPath;

  // ── File pickers ─────────────────────────────────────────
  genomeBtn.addEventListener('click', async () => {
    const path = await window.api.selectFasta();
    if (path) { 
      genomePath = path; 
      genomePathEl.textContent = path;
      try { localStorage.setItem('brimer-genome-path', path); } catch (e) { /* ignore */ }
    }
  });
  gtfBtn.addEventListener('click', async () => {
    const path = await window.api.selectGtf();
    if (path) { 
      gtfPath = path; 
      gtfPathEl.textContent = path;
      try { localStorage.setItem('brimer-gtf-path', path); } catch (e) { /* ignore */ }
    }
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
    inp.className = 'target-name is-gene';
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
        <label>Product min <input class="pp" data-key="product-min" type="number" value="${DEFAULTS['product-min']}"></label>
        <label>Product max <input class="pp" data-key="product-max" type="number" value="${DEFAULTS['product-max']}"></label>
        <label>Opt Tm (°C) <input class="pp" data-key="opt-tm" type="number" value="${DEFAULTS['opt-tm']}" step="0.1"></label>
        <label>Min Tm (°C) <input class="pp" data-key="min-tm" type="number" value="${DEFAULTS['min-tm']}" step="0.1"></label>
        <label>Max Tm (°C) <input class="pp" data-key="max-tm" type="number" value="${DEFAULTS['max-tm']}" step="0.1"></label>
        <label>Opt length <input class="pp" data-key="opt-size" type="number" value="${DEFAULTS['opt-size']}"></label>
        <label>Min length <input class="pp" data-key="min-size" type="number" value="${DEFAULTS['min-size']}"></label>
        <label>Max length <input class="pp" data-key="max-size" type="number" value="${DEFAULTS['max-size']}"></label>
        <label>Min GC% <input class="pp" data-key="min-gc" type="number" value="${DEFAULTS['min-gc']}"></label>
        <label>Max GC% <input class="pp" data-key="max-gc" type="number" value="${DEFAULTS['max-gc']}"></label>
        <label>Pairs/chain <input class="pp" data-key="num-return" type="number" value="${DEFAULTS['num-return']}"></label>
        <label>Max amplicon <input class="pp" data-key="max-amplicon" type="number" value="${DEFAULTS['max-amplicon']}"></label>
      </div>
    `;
    tdParams.appendChild(panel);

    const inputs = panel.querySelectorAll('.pp');
    function updateHighlights() {
      let anyChanged = false;
      inputs.forEach(inpEl => {
        const key = inpEl.dataset.key;
        const val = parseFloat(inpEl.value);
        const isChanged = !isNaN(val) && val !== DEFAULTS[key];
        inpEl.classList.toggle('highlighted', isChanged);
        if (isChanged) anyChanged = true;
      });
      toggleBtn.classList.toggle('highlighted', anyChanged);
    }
    inputs.forEach(inpEl => inpEl.addEventListener('input', updateHighlights));

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

    const debugLink = document.createElement('a');
    debugLink.className = 'debug-archive-link';
    debugLink.style.display = 'none';
    debugLink.textContent = 'Something went wrong. Save this debug archive for troubleshooting.';
    debugLink.href = '#';
    tdResults.appendChild(debugLink);

    const cancelBtn = document.createElement('button');
    cancelBtn.className = 'cancel-btn';
    cancelBtn.textContent = '✕';
    cancelBtn.title = 'Cancel queued job';
    cancelBtn.style.display = 'none';
    tdResults.appendChild(cancelBtn);

    // Assemble row
    tr.appendChild(tdTarget);
    tr.appendChild(tdName);
    tr.appendChild(tdParams);
    tr.appendChild(tdResults);
    targetsBody.appendChild(tr);

    // ── Placeholder switching ───────────────────────
    sel.addEventListener('change', () => {
      const isGene = sel.value === 'gene';
      inp.placeholder = isGene ? GENE_PLACEHOLDER : TRANS_PLACEHOLDER;
      inp.classList.toggle('is-gene', isGene);
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
          // Last field: go to next row
          const nextRow = tr.nextElementSibling;
          if (nextRow) {
            const nextTabbables = getRowTabbables(nextRow);
            if (nextTabbables.length > 0) nextTabbables[0].focus();
            else nextRow.querySelector('.target-type')?.focus();
          }
        }
      }
    });

    // Shared across Run and Cancel handlers
    let requestId;

    // ── Run button ────────────────────────────────
    runBtn.addEventListener('click', async () => {
      if (runBtn.disabled) return;
      if (!genomePath) { showRowError(tr, 'Select a genome FASTA first.'); return; }
      if (!gtfPath) { showRowError(tr, 'Select a GTF annotation first.'); return; }
      const name = inp.value.trim();
      if (!name) { showRowError(tr, 'Enter a target name.'); return; }

      // Lock row inputs forever
      sel.disabled = true;
      inp.readOnly = true;
      panel.querySelectorAll('input').forEach(el => el.disabled = true);

      // Spawn new row if this is the last one
      if (!tr.nextElementSibling) {
        createRow();
      }

      clearRowState(tr);
      runBtn.style.display = 'none';
      cancelBtn.style.display = 'none';
      statusEl.style.display = '';  // let CSS display:inline take effect
      statusEl.textContent = ''; // Set by queue-status event

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
        // Don't await runPipeline() — extract requestId synchronously so
        // data-request-id and runningMap are set before the event loop
        // processes the first queue-status IPC event from main.
        const runResult = window.api.runPipeline({
          genome: genomePath,
          annotations: gtfPath,
          target_key: name,
          target_type: targetType,
          primer_args: primerArgs,
          max_amplicon: maxAmplicon,
          pdf_output_dir: null,
        });

        requestId = runResult.requestId;
        tr.dataset.requestId = requestId;
        const promise = runResult.promise;

        runningMap.set(requestId, tr);

        const result = await promise;
        runningMap.delete(requestId);
        cancelBtn.style.display = 'none';

        if (result && result.cancelled) {
          // Cancelled from queue — row stays locked, shows cancelled
          statusEl.textContent = 'Cancelled';
          return;
        }

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
        const alreadyHandled = requestId ? !runningMap.delete(requestId) : false;
        if (requestId) tr.dataset.requestId = requestId;
        cancelBtn.style.display = 'none';
        statusEl.style.display = 'none';
        if (!alreadyHandled) {
          showRowError(tr, err.message);
        }
      }
    });

    // ── Cancel button ─────────────────────────────
    cancelBtn.addEventListener('click', () => {
      window.api.cancelQueuedJob(requestId).catch(() => {});
      // The run-pipeline promise will reject with '__CANCELLED__' and the
      // catch handler above resets the row.
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

  function showRowError(tr, msg, debugZip) {
    const errorEl = tr.querySelector('.error-msg');
    const tooltip = errorEl.querySelector('.error-tooltip');
    if (tooltip && msg !== null) tooltip.textContent = msg;
    if (msg !== null) errorEl.style.display = 'inline';

    // Persistent debug archive link (never auto-hides)
    const debugLink = tr.querySelector('.debug-archive-link');
    if (debugZip && debugLink) {
      debugLink.dataset.zipPath = debugZip;
      debugLink.style.display = 'inline';
      // Replace to remove old click listeners
      const newLink = debugLink.cloneNode(true);
      debugLink.parentNode.replaceChild(newLink, debugLink);
      newLink.addEventListener('click', function (e) {
        e.preventDefault();
        window.api.openDebugZip(this.dataset.zipPath);
      });
    }
  }

  function clearRowState(tr) {
    const statusEl = tr.querySelector('.status-text');
    const resultLink = tr.querySelector('.result-link');
    const errorEl = tr.querySelector('.error-msg');
    const debugLink = tr.querySelector('.debug-archive-link');
    const cancelBtn = tr.querySelector('.cancel-btn');
    delete tr.dataset.requestId;
    if (statusEl) { statusEl.style.display = 'none'; statusEl.textContent = ''; }
    if (resultLink) { resultLink.style.display = 'none'; resultLink.textContent = ''; resultLink.dataset.pdfPath = ''; }
    if (errorEl) { errorEl.style.display = 'none'; }
    if (debugLink) { debugLink.style.display = 'none'; debugLink.dataset.zipPath = ''; }
    if (cancelBtn) { cancelBtn.style.display = 'none'; }
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
    row.dataset.requestId = data._requestId;
    const statusEl = row.querySelector('.status-text');
    if (statusEl) statusEl.style.display = 'none';
    showRowError(row, data.message || 'Unknown error', data.debug_zip || null);
  });

  // Delayed debug ZIP notification (e.g. from crash fallback).
  // Uses data-request-id on the <tr> so it works even after runningMap
  // entry was deleted by onError.
  window.api.onDebugZip((data) => {
    if (!data._requestId) return;
    const row = document.querySelector('tr[data-request-id="' + data._requestId + '"]');
    if (!row) return;
    showRowError(row, null, data.debug_zip || null);
  });

  // ── Queue status updates ─────────────────────────
  window.api.onQueueUpdate((entries) => {
    entries.forEach(function (queueEntry) {
      var row = document.querySelector('tr[data-request-id="' + queueEntry.requestId + '"]');
      if (!row) return;
      var statusEl = row.querySelector('.status-text');
      var cancelBtn = row.querySelector('.cancel-btn');
      if (!statusEl) return;

      if (queueEntry.status === 'running') {
        if (cancelBtn) cancelBtn.style.display = 'none';
        // Don't overwrite progress messages; only set if still empty
        if (statusEl.textContent === '' || statusEl.textContent.indexOf('Queued') === 0) {
          statusEl.textContent = 'Starting…';
        }
      } else if (queueEntry.status === 'queued') {
        statusEl.textContent = 'Queued (' + queueEntry.position + ' of ' + queueEntry.total + ')';
        if (cancelBtn) cancelBtn.style.display = 'inline';
      }
    });
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

  // ── Version badge ────────────────────────────────────────
  window.api.getVersion().then(function (v) {
    var el = document.getElementById('version-text');
    if (el) el.textContent = 'v' + v;
  }).catch(function () {
    // silent — non-critical
  });

  // ── Initial row ──────────────────────────────────────────
  createRow().querySelector('.target-type').focus();
  } catch (e) {
    if (typeof window.showFatalError === 'function') {
      window.showFatalError(e.message, e.stack);
    } else {
      document.body.textContent = 'Fatal error: ' + (e.message || e);
    }
  }

})();