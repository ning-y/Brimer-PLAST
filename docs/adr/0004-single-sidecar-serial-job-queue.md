# Single sidecar with serial job queue

The Electron app spawns one persistent Python sidecar process per session and
feeds pipeline jobs through it serially (FIFO). Jobs submitted while a previous
job is running join a queue and are displayed as "Queued (N of M)". The queue
prevents resource contention because tnBLAST uses all CPU threads greedily —
running multiple instances concurrently would cause thread thrashing.

Concurrent submission (the previous behavior) was removed because tnBLAST does
not limit its own thread count. A single serialized sidecar guarantees at most
one tnBLAST invocation at a time.

## Key design points

- **Single sidecar, serial jobs.** One `sidecar.py` process stays alive for the
  lifetime of the app session. Jobs are written to its stdin as JSON lines.
  `handle_run_pipeline` is synchronous, so the second command sits in the
  stdin buffer until the first finishes. The sidecar required no code changes.
- **Per-row ✕ cancel (queued only).** A queued-but-not-yet-started job can be
  removed from the queue by clicking ✕. The row's inputs unlock. A cancelled
  row can be re-submitted, which places it at the end of the queue.
- **No cancel of a running job.** Killing the sidecar mid-pipeline is unreliable
  (tnBLAST may leave temp files, primer3 state is opaque). The existing
  debug-archive-on-crash mechanism is sufficient for the running case.
- **Live position updates.** When a job finishes or errors, the remaining queued
  rows' status text is re-rendered (e.g. "Queued (3 of 3)" → "Queued (2 of 2)").
- **Sidecar dies when queue empties.** Once all jobs are done and the queue is
  empty, the sidecar is killed. A later Run press spawns a fresh one.

## Considered Options

- **Concurrent sidecar-per-job (previous behavior).** Rejected because multiple
  tnBLAST instances would compete for CPU threads, degrading throughput.
- **Kill-and-restart for cancellation.** Rejected as unreliable — safer to let
  a running job finish and only allow cancellation before it starts.