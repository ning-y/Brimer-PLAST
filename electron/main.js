const { app, BrowserWindow, ipcMain, dialog, shell } = require('electron');
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const readline = require('readline');

// ── Sidecar process ────────────────────────────────────────────────
let sidecarProcess = null;
let pendingRequest = null;  // { resolve, reject, event } for the current IPC

function findSidecar() {
  // Production: bundled PyInstaller binary in resources/
  const prodPath = path.join(process.resourcesPath, 'pybrimer');
  if (fs.existsSync(prodPath)) return prodPath;
  if (fs.existsSync(prodPath + '.exe')) return prodPath + '.exe';

  // Development: python sidecar.py from repo root
  // Assume we're in electron/ so repo root is ..
  const devPath = path.join(__dirname, '..', 'sidecar.py');
  if (fs.existsSync(devPath)) return { command: 'python3', args: [devPath] };

  // Fallback: try 'python' instead of 'python3'
  return { command: 'python', args: [devPath] };
}

function startSidecar() {
  const bin = findSidecar();

  if (typeof bin === 'string') {
    // PyInstaller binary
    sidecarProcess = spawn(bin, [], { stdio: ['pipe', 'pipe', 'pipe'] });
  } else {
    // Python script
    sidecarProcess = spawn(bin.command, bin.args, { stdio: ['pipe', 'pipe', 'pipe'] });
  }

  // Parse JSON lines from sidecar stdout
  const rl = readline.createInterface({ input: sidecarProcess.stdout });
  rl.on('line', (line) => {
    if (!line.trim()) return;
    try {
      const msg = JSON.parse(line);
      handleSidecarMessage(msg);
    } catch (e) {
      console.error('Failed to parse sidecar output:', line, e);
    }
  });

  // Log stderr for debugging
  sidecarProcess.stderr.on('data', (chunk) => {
    console.error('[sidecar]', chunk.toString().trimEnd());
  });

  // Handle crash
  sidecarProcess.on('close', (code) => {
    console.log(`Sidecar exited with code ${code}`);
    sidecarProcess = null;

    // Notify renderer if a request was pending
    if (pendingRequest) {
      const msg = `Pipeline process exited unexpectedly (code ${code}). Please restart the app.`;
      pendingRequest.event.sender.send('pipeline-error', msg);
      pendingRequest.reject(new Error(msg));
      pendingRequest = null;
    }
  });
}

function stopSidecar() {
  if (sidecarProcess) {
    sidecarProcess.kill();
    sidecarProcess = null;
  }
}

function handleSidecarMessage(msg) {
  if (!pendingRequest) return;

  const { event } = pendingRequest;

  switch (msg.status) {
    case 'progress':
      event.sender.send('pipeline-progress', msg);
      break;
    case 'ok':
      pendingRequest.resolve(msg.result);
      pendingRequest = null;
      break;
    case 'error':
      event.sender.send('pipeline-error', msg.message);
      pendingRequest.reject(new Error(msg.message));
      pendingRequest = null;
      break;
  }
}

// ── Window ─────────────────────────────────────────────────────────
function createWindow() {
  const win = new BrowserWindow({
    width: 900,
    height: 700,
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      nodeIntegration: false,
    },
  });

  win.loadFile(path.join(__dirname, 'renderer', 'index.html'));
}

// ── IPC handlers ───────────────────────────────────────────────────

app.whenReady().then(() => {
  startSidecar();
  createWindow();

  // File dialogs
  ipcMain.handle('select-file', async (_event, opts) => {
    const result = await dialog.showOpenDialog({
      title: opts.title,
      properties: ['openFile'],
      filters: opts.filters,
    });
    if (result.canceled || result.filePaths.length === 0) return null;
    return result.filePaths[0];
  });

  // Run pipeline
  ipcMain.handle('run-pipeline', async (event, params) => {
    // Ensure the sidecar is running
    if (!sidecarProcess) {
      throw new Error('Pipeline process is not running. Please restart the app.');
    }

    // Inject PDF output directory (app data dir)
    const reportsDir = path.join(app.getPath('userData'), 'reports');
    const fullParams = { ...params, pdf_output_dir: reportsDir };

    return new Promise((resolve, reject) => {
      const id = Date.now();
      pendingRequest = { resolve, reject, event };

      const request = JSON.stringify({ id, command: 'run_pipeline', params: fullParams }) + '\n';
      sidecarProcess.stdin.write(request);
    });
  });

  // Open PDF in system viewer
  ipcMain.handle('open-pdf', async (_event, filePath) => {
    shell.openPath(filePath);
  });

  app.on('window-all-closed', () => {
    stopSidecar();
    if (process.platform !== 'darwin') app.quit();
  });

  app.on('before-quit', () => {
    stopSidecar();
  });
});