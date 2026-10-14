// Promise-based wrapper around the HDF5 Web Worker.

export class Hdf5Client {
  constructor() {
    // carry the cache-busting ?v= query of this module over to the worker script
    this.worker = new Worker(new URL(`./hdf5-worker.js${new URL(import.meta.url).search}`, import.meta.url));
    this.pending = new Map();
    this.nextId = 1;
    this.worker.onmessage = (event) => this.onMessage(event.data);
    this.worker.onerror = (event) => {
      const message = event.message || "HDF5 worker failed to start";
      this.pending.forEach(({reject}) => reject(new Error(message)));
      this.pending.clear();
    };
  }

  onMessage(msg) {
    const job = this.pending.get(msg.requestId);
    if (!job) return;
    if (msg.type === "progress") {
      job.onProgress?.(msg);
      return;
    }
    this.pending.delete(msg.requestId);
    if (msg.type === "error") job.reject(new Error(msg.message));
    else job.resolve(msg.payload);
  }

  request(message, {transfer = [], onProgress = null} = {}) {
    const requestId = this.nextId;
    this.nextId += 1;
    return new Promise((resolve, reject) => {
      this.pending.set(requestId, {resolve, reject, onProgress});
      this.worker.postMessage({...message, requestId}, transfer);
    });
  }

  // source: {file: File} | {url: string} | {buffer: ArrayBuffer}
  open(source, onProgress) {
    const transfer = source.buffer ? [source.buffer] : [];
    return this.request({type: "open", ...source}, {transfer, onProgress});
  }

  spectrum(id, bin) {
    return this.request({type: "spectrum", id, bin});
  }

  close(id) {
    return this.request({type: "close", id});
  }
}
