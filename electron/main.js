const { app, BrowserWindow, ipcMain, dialog, shell } = require('electron');
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const readline = require('readline');

// ── Sidecar management ─────────────────────────────────────────────
// Map<sidecarId, { resolve, reject, event, process, requestId }>
const pendingRequests = new Map();
let nextId = 1;

function findSidecar() {
  const prodPath = path.join(process.resourcesPath, 'pybrimer');
  if (fs.existsSync(prodPath)) return prodPath;
  if (fs.existsSync(prodPath + '.exe')) return prodPath + '.exe';
  const devPath = path.join(__dirname, '..', 'sidecar.py');
  if (fs.existsSync(devPath)) return { command: 'python3', args: [devPath] };
  return { command: 'python', args: [devPath] };
}

function spawnSidecar() {
  const bin = findSidecar();
  let proc;
  if (typeof bin === 'string') {
    proc = spawn(bin, [], { stdio: ['pipe', 'pipe', 'pipe'] });
  } else {
    proc = spawn(bin.command, bin.args, { stdio: ['pipe', 'pipe', 'pipe'] });
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
    console.error('[sidecar]', chunk.toString().trimEnd());
  });

  proc.on('close', (code) => {
    // Find and reject all entries using this process
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
      }
    }
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
      entry.process.kill();
      pendingRequests.delete(msg.id);
      entry.resolve(msg.result);
      break;
    case 'error':
      entry.process.kill();
      pendingRequests.delete(msg.id);
      entry.event.sender.send('pipeline-error', {
        status: 'error',
        message: msg.message,
        _requestId: entry.requestId,
      });
      entry.reject(new Error(msg.message));
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
}

// ── IPC handlers ───────────────────────────────────────────────────

app.whenReady().then(() => {
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
    const fullParams = { ...rest, pdf_output_dir: reportsDir };

    return new Promise((resolve, reject) => {
      const sidecarId = nextId++;
      const proc = spawnSidecar();

      pendingRequests.set(sidecarId, { resolve, reject, event, process: proc, requestId: _requestId });

      const request = JSON.stringify({ id: sidecarId, command: 'run_pipeline', params: fullParams }) + '\n';
      proc.stdin.write(request);
      proc.stdin.end();
    });
  });

  ipcMain.handle('open-pdf', async (_event, filePath) => {
    shell.openPath(filePath);
  });

  app.on('window-all-closed', () => {
    if (process.platform !== 'darwin') app.quit();
  });

  app.on('before-quit', () => {
    for (const [, entry] of pendingRequests) {
      entry.process.kill();
    }
    pendingRequests.clear();
  });
});