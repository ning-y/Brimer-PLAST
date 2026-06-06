const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('api', {
  // File pickers
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

  // Pipeline
  runPipeline: (params) => ipcRenderer.invoke('run-pipeline', params),

  // Progress and error listeners
  onProgress: (callback) => {
    ipcRenderer.on('pipeline-progress', (_event, data) => callback(data));
  },
  onError: (callback) => {
    ipcRenderer.on('pipeline-error', (_event, msg) => callback(msg));
  },

  // Open PDF in system viewer
  openPdf: (filePath) => ipcRenderer.invoke('open-pdf', filePath),
});