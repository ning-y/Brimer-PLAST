const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('api', {
  selectFasta: () =>
    ipcRenderer.invoke('select-file', {
      title: 'Select genome FASTA',
      filters: [
        { name: 'FASTA', extensions: ['fna', 'fa', 'fasta', 'fna.gz', 'fa.gz'] },
      ],
    }),
  selectGtf: () =>
    ipcRenderer.invoke('select-file', {
      title: 'Select GTF annotation',
      filters: [{ name: 'GTF', extensions: ['gtf'] }],
    }),

  // Returns { promise, requestId }
  runPipeline: (params) => {
    const requestId = crypto.randomUUID();
    const promise = ipcRenderer.invoke('run-pipeline', { ...params, _requestId: requestId });
    return { promise, requestId };
  },

  // Progress listener — data includes _requestId to route to the right row
  onProgress: (callback) => {
    ipcRenderer.on('pipeline-progress', (_event, data) => callback(data));
  },
  onError: (callback) => {
    ipcRenderer.on('pipeline-error', (_event, data) => callback(data));
  },
  onDebugZip: (callback) => {
    ipcRenderer.on('pipeline-debug-zip', (_event, data) => callback(data));
  },

  openPdf: (filePath) => ipcRenderer.invoke('open-pdf', filePath),
  openDebugZip: (filePath) => ipcRenderer.invoke('open-debug-zip', filePath),
});