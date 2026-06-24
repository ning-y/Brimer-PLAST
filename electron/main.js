const { app, BrowserWindow, ipcMain, dialog, shell } = require('electron');
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const readline = require('readline');
const archiver = require('archiver');

// ── CLI mode detection ────────────────────────────────────────────
const IS_CLI_MODE = process.argv.includes('--cli');

function parseCliArgs() {
  const args = {};
  for (const arg of process.argv) {
    const m = arg.match(/^--([^=]+)=(.*)$/);
    if (m) args[m[1]] = m[2];
  }
  return args;
}

function buildPipelineParams(cliArgs) {
  const genome = cliArgs.genome;
  const annotations = cliArgs.annotations;
  const targetGene = cliArgs['target-gene'];
  const targetTranscript = cliArgs['target-transcript'];

  if (!genome) throw new Error('--genome=<path> is required');
  if (!annotations) throw new Error('--annotations=<path> is required');
  if (!targetGene && !targetTranscript) {
    throw new Error('--target-gene=<name> or --target-transcript=<id> is required');
  }

  const targetType = targetGene ? 'gene' : 'transcript';
  const targetKey = targetGene || targetTranscript;

  const primerArgs = {};
  if (cliArgs['num-return'])       primerArgs.PRIMER_NUM_RETURN = parseInt(cliArgs['num-return'], 10);
  if (cliArgs['min-tm'])           primerArgs.PRIMER_MIN_TM = parseFloat(cliArgs['min-tm']);
  if (cliArgs['max-tm'])           primerArgs.PRIMER_MAX_TM = parseFloat(cliArgs['max-tm']);
  if (cliArgs['opt-tm'])           primerArgs.PRIMER_OPT_TM = parseFloat(cliArgs['opt-tm']);
  if (cliArgs['min-size'])         primerArgs.PRIMER_MIN_SIZE = parseInt(cliArgs['min-size'], 10);
  if (cliArgs['max-size'])         primerArgs.PRIMER_MAX_SIZE = parseInt(cliArgs['max-size'], 10);
  if (cliArgs['opt-size'])         primerArgs.PRIMER_OPT_SIZE = parseInt(cliArgs['opt-size'], 10);
  if (cliArgs['min-gc'])           primerArgs.PRIMER_MIN_GC = parseFloat(cliArgs['min-gc']);
  if (cliArgs['max-gc'])           primerArgs.PRIMER_MAX_GC = parseFloat(cliArgs['max-gc']);
  if (cliArgs['product-min'] || cliArgs['product-max']) {
    const pmin = parseInt(cliArgs['product-min'] || '80', 10);
    const pmax = parseInt(cliArgs['product-max'] || '200', 10);
    primerArgs.PRIMER_PRODUCT_SIZE_RANGE = `${pmin}-${pmax}`;
  }

  return {
    genome,
    annotations,
    target_key: targetKey,
    target_type: targetType,
    primer_args: primerArgs,
    max_amplicon: parseInt(cliArgs['max-amplicon'] || '2000', 10),
    tnblast_timeout: parseInt(cliArgs['tntblast-timeout'] || '3600', 10),
    pdf_output_dir: cliArgs['pdf-dir'] || null,
  };
}

async function runCliPipeline(params) {
  const bin = findSidecar();
  let proc;
  if (typeof bin === 'string') {
    proc = spawn(bin, [], { stdio: ['pipe', 'pipe', 'pipe'] });
  } else {
    proc = spawn(bin.command, bin.args, { stdio: ['pipe', 'pipe', 'pipe'] });
  }

  const sidecarId = nextId++;

  return new Promise((resolve, reject) => {
    const rl = readline.createInterface({ input: proc.stdout });
    rl.on('line', (line) => {
      if (!line.trim()) return;
      try {
        const msg = JSON.parse(line);
        if (msg.id !== sidecarId) return;
        if (msg.status === 'progress') {
          console.error(`[${msg.pct}%] ${msg.message}`);
        } else if (msg.status === 'ok') {
          resolve(msg.result);
        } else if (msg.status === 'error') {
          reject(new Error(msg.message));
        }
      } catch (e) {
        console.error('Failed to parse sidecar output:', line, e);
      }
    });

    proc.stderr.on('data', (chunk) => {
      console.error('[sidecar]', chunk.toString().trimEnd());
    });

    proc.on('close', (code) => {
      reject(new Error(`Pipeline process exited unexpectedly (code ${code}).`));
    });

    const request = JSON.stringify({
      id: sidecarId, command: 'run_pipeline', params,
    }) + '\n';
    proc.stdin.write(request);
    proc.stdin.end();
  });
}
const pendingRequests = new Map();
let nextId = 1;

// ── Serial job queue ────────────────────────────────────────────────
let sidecarProcess = null;
let isProcessing = false;
const jobQueue = [];

function findSidecar() {
  const prodPathExe = path.join(process.resourcesPath, 'pybrimer.exe');
  const prodPath = path.join(process.resourcesPath, 'pybrimer');

  let binPath = null;
  if (process.platform === 'win32') {
    if (fs.existsSync(prodPathExe)) binPath = prodPathExe;
  } else {
    if (fs.existsSync(prodPath)) binPath = prodPath;
  }

  if (binPath) {
    try {
      fs.chmodSync(binPath, 0o755);
    } catch (_) { /* best-effort; may fail on read-only fs */ }
    return binPath;
  }

  const devPath = path.join(__dirname, '..', 'sidecar.py');
  if (fs.existsSync(devPath)) return { command: 'python3', args: [devPath] };
  return { command: 'python', args: [devPath] };
}

function spawnSidecar() {
  const bin = findSidecar();

  // Build tnBLAST path relative to app resources
  const tntblastName = process.platform === 'win32' ? 'tntblast.exe' : 'tntblast';
  const tntblastPath = path.join(process.resourcesPath, tntblastName);
  const env = { ...process.env };
  if (fs.existsSync(tntblastPath)) {
    env.TNTBLAST_PATH = tntblastPath;
  }

  const spawnOpts = { stdio: ['pipe', 'pipe', 'pipe'], env };

  let proc;
  if (typeof bin === 'string') {
    proc = spawn(bin, [], spawnOpts);
  } else {
    proc = spawn(bin.command, bin.args, spawnOpts);
  }

  const rl = readline.createInterface({ input: proc.stdout });
  rl.on('line', (line) => {
    if (!line.trim()) return;
    try {
      const msg = JSON.parse(line);
      handleSidecarMessage(msg);
    } catch (e) {
      console.error('Failed to parse sidecar output:', line, e);
    }
  });

  proc.stderr.on('data', (chunk) => {
    const str = chunk.toString();
    console.error('[sidecar]', str.trimEnd());
    // Accumulate stderr for all entries using this process
    for (const [, entry] of pendingRequests) {
      if (entry.process === proc) {
        entry.stderr += str;
      }
    }
  });

  proc.on('close', (code) => {
    sidecarProcess = null;
    // Reject the in-flight request (if any)
    for (const [sid, entry] of pendingRequests) {
      if (entry.process === proc) {
        const errMsg = `Pipeline process exited unexpectedly (code ${code}).`;
        entry.event.sender.send('pipeline-error', {
          status: 'error',
          message: errMsg,
          _requestId: entry.requestId,
        });
        entry.reject(new Error(errMsg));
        pendingRequests.delete(sid);

        // Fire-and-forget: create fallback ZIP and notify renderer when ready
        if (entry.debugDir) {
          createFallbackDebugZipAsync(entry.debugDir, entry.requestId, entry.stderr || '')
            .then((zipPath) => {
              if (zipPath) {
                entry.event.sender.send('pipeline-debug-zip', {
                  _requestId: entry.requestId,
                  debug_zip: zipPath,
                });
              }
            })
            .catch((e) => console.error('Fallback ZIP creation failed:', e));
        }
      }
    }
    isProcessing = false;
    // Reject all queued jobs — the sidecar is dead
    const errMsg = `Pipeline process exited unexpectedly (code ${code}).`;
    for (const entry of jobQueue) {
      entry.event.sender.send('pipeline-error', {
        status: 'error',
        message: errMsg,
        _requestId: entry.requestId,
      });
      entry.reject(new Error(errMsg));
    }
    jobQueue.length = 0;
    broadcastQueueStatus();
  });

  return proc;
}

function handleSidecarMessage(msg) {
  const entry = pendingRequests.get(msg.id);
  if (!entry) {
    console.warn(`Response for unknown sidecar request ${msg.id}`);
    return;
  }

  switch (msg.status) {
    case 'progress':
      entry.event.sender.send('pipeline-progress', {
        ...msg,
        _requestId: entry.requestId,
      });
      break;
    case 'ok':
      pendingRequests.delete(msg.id);
      isProcessing = false;
      entry.resolve(msg.result);
      processNextInQueue();
      break;
    case 'error':
      pendingRequests.delete(msg.id);
      isProcessing = false;
      entry.event.sender.send('pipeline-error', {
        status: 'error',
        message: msg.message,
        debug_zip: msg.debug_zip || null,
        _requestId: entry.requestId,
      });
      entry.reject(new Error(msg.message));
      processNextInQueue();
      break;
  }
}

function cancelAllForEvent(event) {
  for (const [sid, entry] of pendingRequests) {
    if (entry.event.sender === event.sender) {
      entry.process.kill();
      pendingRequests.delete(sid);
    }
  }
  // Also drain the queue for this sender
  const remaining = [];
  for (const entry of jobQueue) {
    if (entry.event.sender === event.sender) {
      entry.reject(new Error('__CANCELLED__'));
    } else {
      remaining.push(entry);
    }
  }
  jobQueue.length = 0;
  jobQueue.push(...remaining);
  broadcastQueueStatus();
}

function killAllSidecars() {
  for (const [, entry] of pendingRequests) {
    entry.process.kill();
  }
  pendingRequests.clear();
  // Drain the queue
  for (const entry of jobQueue) {
    entry.reject(new Error('__CANCELLED__'));
  }
  jobQueue.length = 0;
  if (sidecarProcess) {
    sidecarProcess.kill();
    sidecarProcess = null;
  }
  isProcessing = false;
}

// ── Serial job queue management ─────────────────────────────────

function processNextInQueue() {
  if (isProcessing) return;
  if (jobQueue.length === 0) {
    killSidecarIfIdle();
    return;
  }

  isProcessing = true;

  if (!sidecarProcess) {
    sidecarProcess = spawnSidecar();
  }

  const entry = jobQueue.shift();
  const sidecarId = nextId++;

  pendingRequests.set(sidecarId, {
    resolve: entry.resolve,
    reject: entry.reject,
    event: entry.event,
    process: sidecarProcess,
    requestId: entry.requestId,
    debugDir: entry.debugDir,
    stderr: '',
  });

  const request = JSON.stringify({ id: sidecarId, command: 'run_pipeline', params: entry.params }) + '\n';
  sidecarProcess.stdin.write(request);
  broadcastQueueStatus();
}

function broadcastQueueStatus() {
  const entries = [];
  // Running job (no position — not in the queue)
  for (const [sid, entry] of pendingRequests) {
    entries.push({ requestId: entry.requestId, status: 'running' });
  }
  // Queued jobs — position is 1-based within the queue
  let pos = 1;
  const total = jobQueue.length;
  for (const entry of jobQueue) {
    entries.push({ requestId: entry.requestId, status: 'queued', position: pos, total: total });
    pos++;
  }
  BrowserWindow.getAllWindows().forEach(win => {
    win.webContents.send('queue-status', entries);
  });
}

function killSidecarIfIdle() {
  if (sidecarProcess && pendingRequests.size === 0 && jobQueue.length === 0) {
    sidecarProcess.kill();
    sidecarProcess = null;
  }
}

// ── Custom CLI argument parsing ────────────────────────────────
function parseTntblastTimeout() {
  for (const arg of process.argv) {
    if (arg.startsWith('--tntblast-timeout=')) {
      const val = parseInt(arg.split('=')[1], 10);
      if (!isNaN(val) && val > 0) return val;
    }
  }
  return 3600;  // default 60 min
}

// ── Window ─────────────────────────────────────────────────────────
function createWindow() {
  const win = new BrowserWindow({
    width: 960,
    height: 720,
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      nodeIntegration: false,
    },
  });

  win.loadFile(path.join(__dirname, 'renderer', 'index.html'));

  // Kill all sidecar processes when the window is closed (catches macOS where
  // window-all-closed does not quit the app).
  win.on('closed', killAllSidecars);
}

// ── IPC handlers ───────────────────────────────────────────────────

app.whenReady().then(() => {
  if (IS_CLI_MODE) {
    try {
      const cliArgs = parseCliArgs();
      const params = buildPipelineParams(cliArgs);
      runCliPipeline(params).then((result) => {
        console.log(JSON.stringify(result, null, 2));
        app.exit(0);
      }).catch((err) => {
        console.error(err.message);
        app.exit(1);
      });
    } catch (err) {
      console.error(err.message);
      console.error('Usage: Brimer-PLAST.AppImage --cli --genome=<path> --annotations=<path> --target-gene=<name>');
      app.exit(1);
    }
    return;
  }

  createWindow();

  ipcMain.handle('select-file', async (_event, opts) => {
    const result = await dialog.showOpenDialog({
      title: opts.title,
      properties: ['openFile'],
      filters: opts.filters,
    });
    if (result.canceled || result.filePaths.length === 0) return null;
    return result.filePaths[0];
  });

  ipcMain.handle('run-pipeline', async (event, params) => {
    const { _requestId, pdf_output_dir, ...rest } = params;
    const reportsDir = path.join(app.getPath('userData'), 'reports');
    const debugDir = path.join(app.getPath('userData'), 'debug-logs');
    const fullParams = {
      ...rest,
      pdf_output_dir: reportsDir,
      debug_dir: debugDir,
      tnblast_timeout: parseTntblastTimeout(),
    };

    return new Promise((resolve, reject) => {
      jobQueue.push({ resolve, reject, event, params: fullParams, requestId: _requestId, debugDir });
      processNextInQueue();
      broadcastQueueStatus();
    });
  });

  ipcMain.handle('cancel-queued-job', (_event, requestId) => {
    const idx = jobQueue.findIndex(e => e.requestId === requestId);
    if (idx === -1) return false;
    const [entry] = jobQueue.splice(idx, 1);
    entry.resolve({ cancelled: true });
    broadcastQueueStatus();
    return true;
  });

  ipcMain.handle('get-app-title', async () => {
    try {
      const { execSync } = require('child_process');
      const out = execSync('python3 -c "from brimer_plast import APP_TITLE; print(APP_TITLE)"', {
        encoding: 'utf-8',
        timeout: 10000,
      });
      return out.trim();
    } catch (_) {
      return 'Brimer-PLAST by Wang Linfa Lab';
    }
  });

  ipcMain.handle('get-version', async () => {
    try {
      // Dev path: ask Python directly (gives full dev/local suffixes)
      const { execSync } = require('child_process');
      const out = execSync('python3 -c "from brimer_plast import __version__; print(__version__)"', {
        encoding: 'utf-8',
        timeout: 10000,
      });
      const v = out.trim();
      if (v) return v;
    } catch (_) {}
    // Production path: use version embedded in Electron package.json
    try {
      const v = app.getVersion();
      if (v && v !== '0.0.0') return v;
    } catch (_) {}
    return '';
  });

  ipcMain.handle('open-pdf', async (_event, filePath) => {
    shell.openPath(filePath);
  });

  // ── Debug ZIP (opening from renderer) ─────────────────────────────────
  ipcMain.handle('open-debug-zip', async (_event, filePath) => {
    shell.openPath(filePath);
  });

  app.on('window-all-closed', () => {
    killAllSidecars();
    if (process.platform !== 'darwin') app.quit();
  });

  app.on('before-quit', killAllSidecars);
});


// ── Fallback debug ZIP (async, Node.js side) ────────────────────────────────
// Used when the sidecar crashes without writing an error response.  Collects
// whatever files survive in the debug directory and bundles them into a ZIP.

async function createFallbackDebugZipAsync(debugDir, requestId, stderrText) {
  if (!debugDir) return null;
  try {
    const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
    const zipName = `debug_crash_${requestId}_${timestamp}.zip`;
    const zipPath = path.join(debugDir, zipName);

    await new Promise((resolve, reject) => {
      const output = fs.createWriteStream(zipPath);
      const archive = archiver('zip', { zlib: { level: 9 } });

      output.on('close', resolve);
      archive.on('error', reject);

      archive.pipe(output);

      // JSONL log (surviving partial log from sidecar)
      const logPath = path.join(debugDir, `pipeline_${requestId}.jsonl`);
      if (fs.existsSync(logPath)) {
        archive.file(logPath, { name: 'pipeline_log.jsonl' });
      }

      // Sidecar stderr
      archive.append(stderrText || '(no stderr captured)', { name: 'sidecar_stderr.txt' });

      // System info
      const sysInfo = [
        `OS: ${process.platform} ${process.arch}`,
        `Node.js: ${process.version}`,
        `Electron: ${process.versions.electron}`,
      ].join('\n') + '\n';
      archive.append(sysInfo, { name: 'system_info.txt' });

      // Any surviving debug artifacts
      for (const fname of ['assays.txt', 'tntblast_genome.txt', 'tntblast_transcriptome.txt']) {
        const fpath = path.join(debugDir, fname);
        if (fs.existsSync(fpath)) {
          archive.file(fpath, { name: fname });
        }
      }

      archive.finalize();
    });

    return zipPath;
  } catch (e) {
    console.error('Failed to create fallback debug ZIP:', e);
    return null;
  }
}