require("dotenv").config();
const express = require("express");
const path = require("path");
const https = require("https");
const http  = require("http");
const fs    = require("fs");
const { execSync } = require("child_process");
const Database = require("better-sqlite3");

const PROCESSING_SERVICE_URL = process.env.PROCESSING_SERVICE_URL || "http://127.0.0.1:5200";
const NAS_WATCH_PATH = process.env.NAS_WATCH_PATH || "/mnt/nas/input/pyl/astro/input";

function callProcessingService(job_id, fits_dir, target, selected_files = null, target_filter = null) {
  return new Promise((resolve, reject) => {
    const payload = { job_id, fits_dir, target };
    if (selected_files && selected_files.length) payload.selected_files = selected_files;
    if (target_filter) payload.target_filter = target_filter;
    // When the fits_dir is a SNAPSHOT folder, pass the target name as object_filter
    // so only frames matching that OBJECT header are processed (not all mixed targets).
    if (fits_dir && path.basename(fits_dir).toUpperCase() === "SNAPSHOT" && target) {
      payload.object_filter = target;
    }
    const body = JSON.stringify(payload);
    const url  = new URL(PROCESSING_SERVICE_URL + "/process");
    const lib  = url.protocol === "https:" ? https : http;
    const req  = lib.request(
      { hostname: url.hostname, port: url.port || (url.protocol === "https:" ? 443 : 80),
        path: url.pathname, method: "POST",
        headers: { "Content-Type": "application/json", "Content-Length": Buffer.byteLength(body) } },
      (res) => {
        let data = "";
        res.on("data", (c) => (data += c));
        res.on("end", () => {
          let parsed;
          try { parsed = JSON.parse(data); } catch { parsed = { raw: data }; }
          // 409 = job already active in processing service → treat as idempotent success
          if (res.statusCode === 409) {
            resolve({ success: true, alreadyActive: true, ...parsed });
          } else if (res.statusCode >= 400) {
            resolve({ success: false, error: parsed.error || `HTTP ${res.statusCode}`, ...parsed });
          } else {
            resolve(parsed);
          }
        });
      }
    );
    req.on("error", reject);
    req.write(body);
    req.end();
  });
}

const app = express();

const PORT = Number(process.env.PORT || 3000);
const DEFAULT_NINA_HOST = process.env.NINA_HOST || "192.168.1.174";
const DEFAULT_NINA_PORT = process.env.NINA_PORT || "1888";
const DEFAULT_NINA_PROTOCOL = process.env.NINA_PROTOCOL || "http";
const DEFAULT_OCS_HOST = process.env.OCS_HOST || "192.168.1.220";
const TNS_BOT_ID = process.env.TNS_BOT_ID || "";
const TNS_BOT_NAME = process.env.TNS_BOT_NAME || "";
const TNS_API_KEY = process.env.TNS_API_KEY || "";
const DEVICE_TYPES = ["mount", "camera", "filterwheel", "focuser"];
const STATUS_DEVICE_MAP = {
  mount: ["Mount", "Connected"],
  camera: ["Camera", "Connected"],
  filterwheel: ["FilterWheel", "Connected"],
  focuser: ["Focuser", "Connected"],
  dome: ["Dome", "Connected"],
  guider: ["Guider", "Connected"],
  rotator: ["Rotator", "Connected"],
  safetymonitor: ["SafetyMonitor", "Connected"],
  switch: ["Switch", "Connected"],
  weather: ["WeatherData", "Connected"],
  flatdevice: ["FlatDevice", "Connected"],
};

app.use(express.json());
app.use(express.static(path.join(__dirname, "public")));

// ── Processed data browser ────────────────────────────────────────────────────
const DATA_DIR = path.join(__dirname, "..", "data");

// Serve files for download: GET /data/2026-03-19/SN%202026fvx/G/120s/res.fit
app.use("/data", express.static(DATA_DIR, { dotfiles: "deny" }));

// List all result files: GET /api/data/results
app.get("/api/data/results", (req, res) => {
  const results = [];
  function walk(dir, rel = "") {
    let entries;
    try { entries = fs.readdirSync(dir, { withFileTypes: true }); } catch { return; }
    for (const e of entries) {
      const relPath = rel ? `${rel}/${e.name}` : e.name;
      if (e.isDirectory()) { walk(path.join(dir, e.name), relPath); }
      else if (e.name === "res.fit") {
        const stat = fs.statSync(path.join(dir, e.name));
        const parts = relPath.split("/");
        const dirPath = path.dirname(relPath);
        const previewRel = `${dirPath}/res_preview.png`;
        const previewExists = fs.existsSync(path.join(DATA_DIR, previewRel));
        const encode = s => encodeURIComponent(s);
        results.push({
          path:        relPath,
          url:         `/data/${relPath.split("/").map(encode).join("/")}`,
          preview_url: previewExists ? `/data/${previewRel.split("/").map(encode).join("/")}` : null,
          date:        parts[0] || "",
          target:      parts[1] || "",
          filter:      parts[2] || "",
          exp:         parts[3] || "",
          size_mb:     (stat.size / 1e6).toFixed(1),
          mtime:       stat.mtime.toISOString(),
        });
      }
    }
  }
  walk(DATA_DIR);
  results.sort((a, b) => b.mtime.localeCompare(a.mtime));
  res.json(results);
});

// ── Persistent DB (targets + pipeline jobs) ───────────────────────────────────

const DB_FILE = process.env.DB_FILE
  ? (path.isAbsolute(process.env.DB_FILE) ? process.env.DB_FILE : path.join(__dirname, process.env.DB_FILE))
  : path.join(__dirname, "..", "nightmanager.db");
const db = new Database(DB_FILE);
db.exec(`
  CREATE TABLE IF NOT EXISTS targets (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    name       TEXT    NOT NULL,
    ra         TEXT,
    dec        TEXT,
    ra_deg     REAL,
    dec_deg    REAL,
    type       TEXT,
    source     TEXT    DEFAULT 'manual',
    notes      TEXT,
    created_at TEXT    DEFAULT (datetime('now'))
  );
  CREATE TABLE IF NOT EXISTS pipeline_jobs (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    target          TEXT,
    filter          TEXT,
    exposure        TEXT,
    fits_dir        TEXT,
    status          TEXT DEFAULT 'queued',
    stdweb_task_id  TEXT,
    stdweb_url      TEXT,
    error           TEXT,
    created_at      TEXT DEFAULT (datetime('now')),
    updated_at      TEXT DEFAULT (datetime('now'))
  );
  CREATE TABLE IF NOT EXISTS pipeline_results (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    job_id          INTEGER NOT NULL REFERENCES pipeline_jobs(id) ON DELETE CASCADE,
    target          TEXT,
    filter          TEXT,
    exposure        TEXT,
    n_frames        INTEGER,
    status          TEXT DEFAULT 'processing',
    stdweb_task_id  TEXT,
    stdweb_url      TEXT,
    error           TEXT,
    created_at      TEXT DEFAULT (datetime('now')),
    updated_at      TEXT DEFAULT (datetime('now'))
  );
  CREATE TABLE IF NOT EXISTS seq_queue (
    position   INTEGER PRIMARY KEY,
    name       TEXT    NOT NULL,
    ra         TEXT,
    dec        TEXT,
    ra_deg     REAL,
    dec_deg    REAL,
    done       INTEGER DEFAULT 0,
    added_at   INTEGER
  );
`);

// Migrations
try { db.prepare("ALTER TABLE pipeline_results ADD COLUMN target TEXT").run(); } catch { /* already exists */ }
try { db.prepare("ALTER TABLE pipeline_results ADD COLUMN frame_previews TEXT").run(); } catch { /* already exists */ }
try { db.prepare("ALTER TABLE pipeline_jobs ADD COLUMN target_filter TEXT").run(); } catch { /* already exists */ }
try { db.prepare("ALTER TABLE pipeline_results ADD COLUMN stdweb_state TEXT").run(); } catch { /* already exists */ }
try { db.prepare("ALTER TABLE pipeline_results ADD COLUMN obs_date TEXT").run(); } catch { /* already exists */ }
try { db.prepare("ALTER TABLE pipeline_results ADD COLUMN mjd REAL").run(); } catch { /* already exists */ }
try { db.prepare("ALTER TABLE pipeline_results ADD COLUMN mag_ap REAL").run(); } catch { /* already exists */ }
try { db.prepare("ALTER TABLE pipeline_results ADD COLUMN magerr_ap REAL").run(); } catch { /* already exists */ }
try { db.prepare("ALTER TABLE pipeline_results ADD COLUMN mag_sub REAL").run(); } catch { /* already exists */ }
try { db.prepare("ALTER TABLE pipeline_results ADD COLUMN magerr_sub REAL").run(); } catch { /* already exists */ }
try { db.prepare("ALTER TABLE pipeline_results ADD COLUMN mag_sub_ul REAL").run(); } catch { /* already exists */ }
db.exec(`
  CREATE TABLE IF NOT EXISTS ocs_history (
    id        INTEGER PRIMARY KEY AUTOINCREMENT,
    ts        TEXT NOT NULL DEFAULT (datetime('now')),
    roof      TEXT,
    safe      TEXT,
    rain      TEXT,
    temp      REAL,
    humidity  REAL,
    pressure  REAL,
    sky       REAL
  );
`);

// ── OCS history poller ────────────────────────────────────────────────────────

/** Strip HTML tags and trim whitespace from an OCS field value. */
function ocsClean(v) {
  return (v || "").replace(/<[^>]*>/g, "").replace(/&[a-z]+;/gi, "").trim();
}

async function pollAndStoreOcs() {
  try {
    const host = DEFAULT_OCS_HOST;
    const { ok, fields } = await ocsStatus(host);
    if (!ok) return;
    const temp     = parseFloat(fields.wea_temp);
    const humidity = parseFloat(fields.wea_humd);
    const pressure = parseFloat(fields.wea_pres);
    const sky      = parseFloat(fields.wea_sq || "");   // e.g. " 5.4 mpsas"
    db.prepare(
      `INSERT INTO ocs_history (roof, safe, rain, temp, humidity, pressure, sky) VALUES (?, ?, ?, ?, ?, ?, ?)`
    ).run(
      ocsClean(fields.roof_sta)  || null,
      ocsClean(fields.stat_safe) || null,
      ocsClean(fields.wea_rain)  || null,
      isNaN(temp)     ? null : temp,
      isNaN(humidity) ? null : humidity,
      isNaN(pressure) ? null : pressure,
      isNaN(sky)      ? null : sky,
    );
  } catch (e) {
    // OCS unreachable — skip silently (don't fill logs with noise)
    if (process.env.DEBUG) console.error("[OCS poller]", e.message);
  }
}

setInterval(pollAndStoreOcs, 10 * 60 * 1000);  // every 10 min
setTimeout(pollAndStoreOcs, 8000);              // initial reading after server startup

// ── Queue persistence helpers ─────────────────────────────────────────────────

const _saveQueue = db.transaction(() => {
  db.prepare("DELETE FROM seq_queue").run();
  const ins = db.prepare(
    "INSERT INTO seq_queue (position, name, ra, dec, ra_deg, dec_deg, done, added_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
  );
  seqState.queue.forEach((t, i) => {
    ins.run(i, t.name, t.ra || "", t.dec || "", t.raDeg, t.decDeg, t.done ? 1 : 0, t.addedAt || Date.now());
  });
});

function saveQueue() {
  try { _saveQueue(); } catch (e) { console.error("saveQueue error:", e.message); }
}

function loadQueue() {
  try {
    const rows = db.prepare("SELECT * FROM seq_queue ORDER BY position").all();
    seqState.queue = rows.map(r => ({
      name:     r.name,
      ra:       r.ra,
      dec:      r.dec,
      raDeg:    r.ra_deg,
      decDeg:   r.dec_deg,
      done:     r.done === 1,
      addedAt:  r.added_at,
    }));
    if (seqState.queue.length) {
      console.log(`Restored ${seqState.queue.length} target(s) from seq_queue`);
    }
  } catch (e) {
    console.error("loadQueue error:", e.message);
  }
}

function normalizeTargetConfig(payload = {}) {
  const host = String(payload.host || DEFAULT_NINA_HOST).trim();
  const protocol = String(payload.protocol || DEFAULT_NINA_PROTOCOL).trim().toLowerCase();
  const port = String(payload.port || DEFAULT_NINA_PORT).trim();

  if (!host) {
    throw new Error("Host is required");
  }

  if (!/^\d+$/.test(port)) {
    throw new Error("Port must be numeric");
  }

  if (!["http", "https"].includes(protocol)) {
    throw new Error("Protocol must be http or https");
  }

  return { host, port, protocol };
}

function buildBaseUrl(target) {
  return `${target.protocol}://${target.host}:${target.port}`;
}

function buildUrl(target, endpointPath, query = {}) {
  const url = new URL(`${buildBaseUrl(target)}${endpointPath}`);
  for (const [key, value] of Object.entries(query)) {
    if (value === undefined || value === null || value === "") {
      continue;
    }
    url.searchParams.set(key, String(value));
  }
  return url.toString();
}

async function callNina(target, endpointPath, query = {}) {
  const url = buildUrl(target, endpointPath, query);
  const controller = new AbortController();
  const timeout = setTimeout(() => controller.abort(), 10000);

  try {
    const response = await fetch(url, {
      method: "GET",
      signal: controller.signal,
      headers: {
        Accept: "application/json",
      },
    });

    const text = await response.text();
    let body;
    try {
      body = text ? JSON.parse(text) : null;
    } catch {
      body = text;
    }

    return {
      ok: response.ok,
      status: response.status,
      body,
      url,
    };
  } finally {
    clearTimeout(timeout);
  }
}

function deviceStatusFromEquipmentInfo(equipmentInfo) {
  const statuses = {};
  for (const [device, path] of Object.entries(STATUS_DEVICE_MAP)) {
    statuses[device] = Boolean(equipmentInfo?.[path[0]]?.[path[1]]);
  }
  return statuses;
}

function isNinaApiSuccess(result) {
  if (!result.ok) {
    return false;
  }

  if (result.body && typeof result.body === "object" && "Success" in result.body) {
    return Boolean(result.body.Success);
  }

  return true;
}

function normalizeNinaError(error) {
  return error.name === "AbortError" ? "Connection timed out" : error.message;
}

async function getEquipmentInfo(target) {
  const result = await callNina(target, "/v2/api/equipment/info");
  const equipmentInfo = typeof result.body === "object" ? result.body?.Response || null : null;
  const statuses = equipmentInfo ? deviceStatusFromEquipmentInfo(equipmentInfo) : null;
  return { result, equipmentInfo, statuses };
}

app.get("/api/config/defaults", (req, res) => {
  res.json({
    host: DEFAULT_NINA_HOST,
    port: DEFAULT_NINA_PORT,
    protocol: DEFAULT_NINA_PROTOCOL,
    ocsHost: DEFAULT_OCS_HOST,
    tns: {
      botId: TNS_BOT_ID,
      botName: TNS_BOT_NAME,
      hasApiKey: Boolean(TNS_API_KEY),
    },
  });
});

app.post("/api/nina/test", async (req, res) => {
  let target;
  try {
    target = normalizeTargetConfig(req.body);
  } catch (error) {
    return res.status(400).json({ success: false, error: error.message });
  }

  try {
    const { result, equipmentInfo, statuses } = await getEquipmentInfo(target);
    const success = isNinaApiSuccess(result);
    return res.status(success ? 200 : 502).json({
      success,
      target,
      status: result.status,
      url: result.url,
      devices: statuses,
      equipment: equipmentInfo,
      raw: result.body,
    });
  } catch (error) {
    return res.status(502).json({
      success: false,
      target,
      error: normalizeNinaError(error),
    });
  }
});

app.post("/api/nina/devices/status", async (req, res) => {
  let target;
  try {
    target = normalizeTargetConfig(req.body);
  } catch (error) {
    return res.status(400).json({ success: false, error: error.message });
  }

  try {
    const { result, equipmentInfo, statuses } = await getEquipmentInfo(target);
    const success = isNinaApiSuccess(result);
    return res.status(success ? 200 : 502).json({
      success,
      target,
      status: result.status,
      devices: statuses,
      equipment: equipmentInfo,
      raw: result.body,
    });
  } catch (error) {
    return res.status(502).json({
      success: false,
      target,
      error: normalizeNinaError(error),
    });
  }
});

app.post("/api/nina/devices/connect", async (req, res) => {
  let target;
  try {
    target = normalizeTargetConfig(req.body);
  } catch (error) {
    return res.status(400).json({ success: false, error: error.message });
  }

  const requestedDevice = req.body?.device;
  const selectedDevices = requestedDevice === "all"
    ? DEVICE_TYPES
    : DEVICE_TYPES.includes(requestedDevice)
      ? [requestedDevice]
      : null;

  if (!selectedDevices) {
    return res.status(400).json({
      success: false,
      error: "Device must be one of: all, mount, camera, filterwheel, focuser",
    });
  }

  try {
    const connectionResults = await Promise.all(
      selectedDevices.map(async (device) => {
        const result = await callNina(target, `/v2/api/equipment/${device}/connect`);
        return {
          device,
          success: isNinaApiSuccess(result),
          status: result.status,
          response: result.body,
        };
      }),
    );

    const { result: infoResult, statuses } = await getEquipmentInfo(target);

    const allSuccess = connectionResults.every((result) => result.success);
    return res.status(allSuccess ? 200 : 502).json({
      success: allSuccess,
      target,
      results: connectionResults,
      devices: statuses,
      infoStatus: infoResult.status,
    });
  } catch (error) {
    return res.status(502).json({
      success: false,
      target,
      error: normalizeNinaError(error),
    });
  }
});

app.post("/api/nina/actions/mount", async (req, res) => {
  let target;
  try {
    target = normalizeTargetConfig(req.body);
  } catch (error) {
    return res.status(400).json({ success: false, error: error.message });
  }

  const command = String(req.body?.command || "").trim().toLowerCase();
  const endpointByCommand = {
    park: "/v2/api/equipment/mount/park",
    unpark: "/v2/api/equipment/mount/unpark",
    home: "/v2/api/equipment/mount/home",
  };
  const endpoint = endpointByCommand[command];

  if (!endpoint) {
    return res.status(400).json({
      success: false,
      error: "Command must be one of: park, unpark, home",
    });
  }

  try {
    const actionResult = await callNina(target, endpoint);
    const actionSuccess = isNinaApiSuccess(actionResult);
    const { statuses } = await getEquipmentInfo(target);
    return res.status(actionSuccess ? 200 : 502).json({
      success: actionSuccess,
      target,
      command,
      result: actionResult.body,
      devices: statuses,
    });
  } catch (error) {
    return res.status(502).json({
      success: false,
      target,
      error: normalizeNinaError(error),
    });
  }
});

app.post("/api/nina/actions/camera/capture", async (req, res) => {
  let target;
  try {
    target = normalizeTargetConfig(req.body);
  } catch (error) {
    return res.status(400).json({ success: false, error: error.message });
  }

  const duration = Number(req.body?.duration);
  const gain = Number(req.body?.gain);
  const savePath = req.body?.savePath ? String(req.body.savePath).trim() : "";
  const targetName = req.body?.targetName ? String(req.body.targetName).trim() : "";

  if (!Number.isFinite(duration) || duration <= 0) {
    return res.status(400).json({
      success: false,
      error: "duration must be a positive number (seconds)",
    });
  }

  if (!Number.isFinite(gain) || gain < 0) {
    return res.status(400).json({
      success: false,
      error: "gain must be a number >= 0",
    });
  }

  try {
    let savePathResult = null;
    if (savePath) {
      savePathResult = await callNina(target, "/v2/api/profile/change-value", {
        settingpath: "ImageFileSettings-FilePath",
        newValue: savePath,
      });
    }

    const captureResult = await callNina(target, "/v2/api/equipment/camera/capture", {
      duration,
      gain,
      save: true,
      ...(targetName && { targetName }),
    });

    const captureSuccess = isNinaApiSuccess(captureResult)
      && (!savePathResult || isNinaApiSuccess(savePathResult));
    const { statuses } = await getEquipmentInfo(target);

    return res.status(captureSuccess ? 200 : 502).json({
      success: captureSuccess,
      target,
      capture: captureResult.body,
      savePath,
      savePathUpdate: savePathResult ? savePathResult.body : null,
      devices: statuses,
    });
  } catch (error) {
    return res.status(502).json({
      success: false,
      target,
      error: normalizeNinaError(error),
    });
  }
});

app.post("/api/nina/actions/filterwheel/change", async (req, res) => {
  let target;
  try {
    target = normalizeTargetConfig(req.body);
  } catch (error) {
    return res.status(400).json({ success: false, error: error.message });
  }

  const filterId = Number(req.body?.filterId);
  if (!Number.isInteger(filterId) || filterId < 0) {
    return res.status(400).json({
      success: false,
      error: "filterId must be an integer >= 0",
    });
  }

  try {
    const actionResult = await callNina(target, "/v2/api/equipment/filterwheel/change-filter", {
      filterId,
    });
    const actionSuccess = isNinaApiSuccess(actionResult);
    const { statuses } = await getEquipmentInfo(target);
    return res.status(actionSuccess ? 200 : 502).json({
      success: actionSuccess,
      target,
      filterId,
      result: actionResult.body,
      devices: statuses,
    });
  } catch (error) {
    return res.status(502).json({
      success: false,
      target,
      error: normalizeNinaError(error),
    });
  }
});

app.post("/api/nina/actions/focuser/move-relative", async (req, res) => {
  let target;
  try {
    target = normalizeTargetConfig(req.body);
  } catch (error) {
    return res.status(400).json({ success: false, error: error.message });
  }

  const direction = String(req.body?.direction || "").trim().toLowerCase();
  const steps = Number(req.body?.steps);
  if (!["in", "out"].includes(direction)) {
    return res.status(400).json({
      success: false,
      error: "direction must be 'in' or 'out'",
    });
  }
  if (!Number.isInteger(steps) || steps <= 0) {
    return res.status(400).json({
      success: false,
      error: "steps must be a positive integer",
    });
  }

  try {
    const { equipmentInfo, statuses } = await getEquipmentInfo(target);
    const currentPosition = Number(equipmentInfo?.Focuser?.Position);
    if (!Number.isFinite(currentPosition)) {
      return res.status(502).json({
        success: false,
        target,
        error: "Could not read focuser position",
      });
    }

    const signedDelta = direction === "out" ? steps : -steps;
    const targetPosition = Math.max(0, Math.round(currentPosition + signedDelta));

    const moveResult = await callNina(target, "/v2/api/equipment/focuser/move", {
      position: targetPosition,
    });
    const actionSuccess = isNinaApiSuccess(moveResult);
    const updatedInfo = await getEquipmentInfo(target);

    return res.status(actionSuccess ? 200 : 502).json({
      success: actionSuccess,
      target,
      direction,
      steps,
      fromPosition: currentPosition,
      toPosition: targetPosition,
      result: moveResult.body,
      devices: updatedInfo.statuses || statuses,
      focuserPosition: updatedInfo.equipmentInfo?.Focuser?.Position,
    });
  } catch (error) {
    return res.status(502).json({
      success: false,
      target,
      error: normalizeNinaError(error),
    });
  }
});

// ── Sequence Engine ───────────────────────────────────────────────────────────

const seqState = {
  running: false,
  aborted: false,
  queue: [],
  currentTargetIdx: -1,
  currentTarget: null,
  currentStep: null,
  progress: null,      // { filter, frame, frames }
  log: [],
  lastAutofocusTime: null,
  error: null,
  manualMode: false,
  waitingForStep: null, // non-null when paused waiting for /api/sequence/next
  ninaConfig: null,    // set when sequence starts, used by abort to reach NINA
};

loadQueue();

// Resolvers for the manual-step confirmation promise
let _confirmResolve = null;
let _confirmReject  = null;

const SEQ_MAX_LOG = 300;

function seqLog(msg, level = "info") {
  const ts = new Date().toLocaleTimeString("en-GB", { hour12: false });
  const entry = { ts, msg, level };
  seqState.log.push(entry);
  if (seqState.log.length > SEQ_MAX_LOG) seqState.log.shift();
  console.log(`[SEQ ${level.toUpperCase()}] ${ts} ${msg}`);
}

function checkAbort() {
  if (seqState.aborted) throw new Error("__ABORTED__");
}

/**
 * In manual mode: pause the sequence and wait for the user to press "Next".
 * In auto mode: returns immediately.
 */
async function waitForConfirmation(stepDesc) {
  if (!seqState.manualMode) return;
  checkAbort();
  seqState.waitingForStep = stepDesc;
  seqLog(`⏸  Paused — press "Next Step" to continue`, "warn");
  await new Promise((resolve, reject) => {
    _confirmResolve = resolve;
    _confirmReject  = reject;
  });
  seqState.waitingForStep = null;
  seqLog(`▶  Confirmed: ${stepDesc}`);
  checkAbort();
}

/** Like callNina but with a custom timeout — needed for long exposures */
async function callNinaLong(target, endpointPath, query = {}, timeoutMs = 30000) {
  const url = buildUrl(target, endpointPath, query);
  const controller = new AbortController();
  const timer = setTimeout(() => controller.abort(), timeoutMs);
  try {
    const response = await fetch(url, {
      method: "GET",
      signal: controller.signal,
      headers: { Accept: "application/json" },
    });
    const text = await response.text();
    let body;
    try { body = text ? JSON.parse(text) : null; } catch { body = text; }
    return { ok: response.ok, status: response.status, body, url };
  } finally {
    clearTimeout(timer);
  }
}

/** Poll equipment info until condFn returns truthy, or timeout */
async function pollEquipment(target, condFn, timeoutMs, intervalMs = 5000) {
  const deadline = Date.now() + timeoutMs;
  while (Date.now() < deadline) {
    checkAbort();
    try {
      const { equipmentInfo } = await getEquipmentInfo(target);
      if (condFn(equipmentInfo)) return equipmentInfo;
    } catch { /* ignore transient errors */ }
    await new Promise(r => setTimeout(r, intervalMs));
  }
  throw new Error("Timeout waiting for condition");
}

// ── Sequence step helpers ─────────────────────────────────────────────────────

async function stepCoolCamera(target, targetTempC = -5) {
  seqState.currentStep = "Checking camera temperature...";
  seqLog("Checking camera temperature...");
  const { equipmentInfo } = await getEquipmentInfo(target);
  const camera = equipmentInfo?.Camera;

  if (!camera?.Connected) {
    seqLog("Camera not connected — skipping cooling", "warn");
    return;
  }

  const currentTemp = Number(camera.Temperature);
  if (currentTemp <= targetTempC + 0.5) {
    seqLog(`Camera already at ${currentTemp.toFixed(1)}°C ✓`);
    return;
  }

  seqLog(`Cooling camera to ${targetTempC}°C (currently ${currentTemp.toFixed(1)}°C)...`);
  // Endpoint: /v2/api/equipment/camera/cool?temperature=X&minutes=Y
  await callNinaLong(target, "/v2/api/equipment/camera/cool", {
    temperature: targetTempC,
    minutes: 10,
  }, 15000);

  // Poll until within 0.5°C of target
  await pollEquipment(
    target,
    (info) => {
      const t = Number(info?.Camera?.Temperature);
      seqState.currentStep = `Cooling: ${t.toFixed(1)}°C → ${targetTempC}°C`;
      return t <= targetTempC + 0.5;
    },
    20 * 60 * 1000,
    10000,
  );
  seqLog(`Camera cooled to ${targetTempC}°C ✓`);
}

async function stepUnparkMount(target) {
  seqState.currentStep = "Preparing mount...";
  const { equipmentInfo } = await getEquipmentInfo(target);
  const mount = equipmentInfo?.Mount;

  if (!mount?.Connected) {
    throw new Error("Mount not connected");
  }

  // Log the actual mount state for transparency
  seqLog(
    `  Mount: AtPark=${mount.AtPark}  AtHome=${mount.AtHome}  ` +
    `Tracking=${mount.TrackingEnabled ? mount.TrackingMode : "OFF"}  ` +
    `RA=${mount.RightAscensionString ?? "?"}  Dec=${mount.DeclinationString ?? "?"}`
  );

  // Unpark only if formally parked — then poll until the mount confirms AtPark=false
  if (mount.AtPark === true) {
    seqLog("Unparking mount...");
    await callNinaLong(target, "/v2/api/equipment/mount/unpark", {}, 60000);

    // Wait for the mount to physically clear the parked flag (up to 30 s)
    const UNPARK_TIMEOUT  = 30000;
    const UNPARK_INTERVAL = 1500;
    const deadline = Date.now() + UNPARK_TIMEOUT;
    let unparked = false;
    while (Date.now() < deadline) {
      await new Promise(r => setTimeout(r, UNPARK_INTERVAL));
      const { equipmentInfo: ei } = await getEquipmentInfo(target);
      if (ei?.Mount?.AtPark === false) { unparked = true; break; }
    }
    if (!unparked) seqLog("  Unpark confirmation timed out — continuing anyway", "warn");
    else           seqLog("Mount unparked ✓");
  } else if (mount.AtHome === true) {
    seqLog("Mount at home position — ready to slew ✓");
  } else {
    seqLog("Mount ready ✓");
  }

  // Re-read mount state after potential unpark before touching tracking
  const { equipmentInfo: freshEI } = await getEquipmentInfo(target);
  const freshMount = freshEI?.Mount ?? mount;

  // Enable sidereal tracking if off
  if (!freshMount.TrackingEnabled) {
    seqLog("Enabling sidereal tracking...");
    const tr = await callNinaLong(target, "/v2/api/equipment/mount/tracking", { on: true }, 15000);
    if (!isNinaApiSuccess(tr)) {
      seqLog(`  Tracking enable warning: ${tr.body?.Error || tr.status}`, "warn");
    } else {
      seqLog("Sidereal tracking enabled ✓");
    }
  } else {
    seqLog(`Tracking already on (${freshMount.TrackingMode}) ✓`);
  }
}

// ── Server-side spherical astronomy (mirrors public/astronomy.js) ─────────────

const _D2R = Math.PI / 180;
const _R2D = 180 / Math.PI;

/** Angular separation in arcseconds between two points (ra/dec in degrees) */
function _angularSepArcsec(ra1Deg, dec1Deg, ra2Deg, dec2Deg) {
  const dRa  = (ra2Deg  - ra1Deg)  * _D2R;
  const dDec = (dec2Deg - dec1Deg) * _D2R;
  const a = Math.sin(dDec / 2) ** 2
    + Math.cos(dec1Deg * _D2R) * Math.cos(dec2Deg * _D2R) * Math.sin(dRa / 2) ** 2;
  return 2 * Math.asin(Math.sqrt(Math.min(1, a))) * _R2D * 3600;
}

function _toJD(date) {
  return date.getTime() / 86400000 + 2440587.5;
}

function _gmstDeg(jd) {
  const T = (jd - 2451545.0) / 36525.0;
  const g =
    280.46061837 +
    360.98564736629 * (jd - 2451545.0) +
    0.000387933 * T * T -
    (T * T * T) / 38710000;
  return ((g % 360) + 360) % 360;
}

function serverComputeAltAz(raDeg, decDeg, latDeg, lonDeg, date = new Date()) {
  const jd  = _toJD(date);
  const lst = (_gmstDeg(jd) + lonDeg + 360) % 360;
  const ha  = (lst - raDeg + 360) % 360;
  const haR  = ha  * _D2R;
  const decR = decDeg * _D2R;
  const latR = latDeg * _D2R;

  const sinAlt =
    Math.sin(decR) * Math.sin(latR) +
    Math.cos(decR) * Math.cos(latR) * Math.cos(haR);
  const altR = Math.asin(Math.max(-1, Math.min(1, sinAlt)));

  const cosAz =
    (Math.sin(decR) - Math.sin(altR) * Math.sin(latR)) /
    (Math.cos(altR) * Math.cos(latR));
  let azR = Math.acos(Math.max(-1, Math.min(1, cosAz)));
  if (Math.sin(haR) > 0) azR = 2 * Math.PI - azR;

  return { alt: altR * _R2D, az: azR * _R2D };
}

/** Abort if target is below the horizon (< MIN_ALT_DEG) */
const MIN_ALT_DEG = 5;

async function checkTargetAboveHorizon(target, raDeg, decDeg, name) {
  const { equipmentInfo } = await getEquipmentInfo(target);
  const mount = equipmentInfo?.Mount;
  const lat = mount?.SiteLatitude;
  const lon = mount?.SiteLongitude;

  if (lat == null || lon == null || lat === 0 && lon === 0) {
    seqLog(`  Horizon check skipped — no observer location in mount (lat=${lat} lon=${lon})`, "warn");
    return;
  }

  const { alt, az } = serverComputeAltAz(raDeg, decDeg, lat, lon, new Date());
  const altStr = alt.toFixed(1);
  const azStr  = az.toFixed(1);

  if (alt < MIN_ALT_DEG) {
    throw new Error(
      `Target ${name} is below the horizon! ` +
      `Alt ${altStr}° Az ${azStr}° from lat=${lat.toFixed(2)}° lon=${lon.toFixed(2)}°. ` +
      `Minimum altitude is ${MIN_ALT_DEG}°.`
    );
  }

  seqLog(`  Horizon check: Alt ${altStr}°  Az ${azStr}°  (observer lat ${lat.toFixed(2)}° lon ${lon.toFixed(2)}°) ✓`);
}

/**
 * Slew to target, optionally with NINA's built-in plate-solve-and-center.
 *
 * NINA slew endpoint: GET /v2/api/equipment/mount/slew
 *   ra            – RA in degrees (Angle.ByDegree)
 *   dec           – Dec in degrees
 *   waitForResult – true → blocks until done; false → returns immediately
 *   center        – true → use NINA Center instruction (capture → solve → sync → re-slew)
 */
async function stepSlewToTarget(target, raDeg, decDeg, name, { center = false } = {}) {
  const raHours = raDeg / 15;
  seqState.currentStep = `Slewing${center ? " + centering" : ""} to ${name}...`;
  seqLog(`Slewing${center ? " + centering" : ""} to ${name} (RA ${raHours.toFixed(4)}h  Dec ${decDeg >= 0 ? "+" : ""}${decDeg.toFixed(4)}°)...`);

  // Tell NINA the target name so it writes OBJECT into every FITS header
  await callNinaLong(target, "/v2/api/equipment/mount/target", {
    name, ra: raDeg, dec: decDeg,
  }, 10000).catch(() => {/* endpoint may not exist on older NINA — non-fatal */});

  // Helper: plain slew (used as fallback when centering fails)
  async function plainSlew() {
    seqState.currentStep = `Slewing to ${name}...`;
    const r = await callNinaLong(target, "/v2/api/equipment/mount/slew", {
      ra: raDeg, dec: decDeg, waitForResult: true, name,
    }, 3 * 60 * 1000);
    if (!isNinaApiSuccess(r)) {
      const msg = r.body?.Error || r.body?.Response || JSON.stringify(r.body) || "(no details)";
      throw new Error(`Slew failed: ${msg}`);
    }
  }

  if (center) {
    // NINA runs the full solve→sync→re-slew loop internally.
    // Its HTTP response reflects only the first motor-slew, not the whole centering
    // cycle, so Success=false / "Slew failed" is normal mid-centering noise — we
    // must NOT fall back; instead we just kick it off and poll until the mount stops.
    try {
      await callNinaLong(target, "/v2/api/equipment/mount/slew", {
        ra: raDeg, dec: decDeg, waitForResult: true, center: true, name,
      }, 10 * 60 * 1000);
    } catch (fetchErr) {
      // Network/timeout error means we never even reached NINA — fall back
      seqLog(`  Centering request error: ${fetchErr.message} — falling back to plain slew`, "warn");
      await plainSlew();
      seqLog("  Plain slew complete (no centering)");
    }

    // Poll until NINA finishes all centering iterations, logging RA/Dec error each pass
    seqState.currentStep = `Centering ${name} — waiting for convergence...`;
    const SETTLE_TIMEOUT  = 8 * 60 * 1000;
    const SETTLE_INTERVAL = 2500;
    const deadline = Date.now() + SETTLE_TIMEOUT;
    let lastErrArcsec = null;

    while (Date.now() < deadline) {
      checkAbort();
      await new Promise(r => setTimeout(r, SETTLE_INTERVAL));
      const { equipmentInfo: ei } = await getEquipmentInfo(target);
      const m = ei?.Mount;
      if (!m) continue;

      // Compute angular separation (arcsec) between current mount pos and target
      if (m.RightAscension != null && m.Declination != null) {
        const errArcsec = _angularSepArcsec(
          m.RightAscension * 15, m.Declination,   // NINA returns RA in hours
          raDeg, decDeg,
        );
        const errRounded = Math.round(errArcsec);
        // Log when error changes by more than 10" (avoids spam while slewing)
        if (lastErrArcsec === null || Math.abs(errRounded - lastErrArcsec) > 10) {
          seqLog(`  Centering: offset ${errRounded}"  (Slewing=${m.Slewing})`);
          lastErrArcsec = errRounded;
          seqState.currentStep = `Centering ${name} — offset ${errRounded}"`;
        }
      }

      if (m.Slewing !== true) break;
    }

    const finalErrTxt = lastErrArcsec !== null ? `  final offset ${lastErrArcsec}"` : "";
    seqLog(`Slew + centering complete ✓${finalErrTxt}`);
  } else {
    await plainSlew();
  }

  // Log final mount position for verification
  const { equipmentInfo: final } = await getEquipmentInfo(target);
  const m = final?.Mount;
  const posStr = `RA ${m?.RightAscensionString ?? "?"}  Dec ${m?.DeclinationString ?? "?"}  Tracking: ${m?.TrackingEnabled ? m?.TrackingMode : "OFF"}`;
  seqLog(`Slew complete ✓  Position: ${posStr}`);
  checkAbort();
}

async function stepStartGuiding(target) {
  seqState.currentStep = "Starting PHD2 guiding...";
  seqLog("Starting PHD2 guiding...");

  // Endpoint: /v2/api/equipment/guider/start
  const result = await callNinaLong(target, "/v2/api/equipment/guider/start", {}, 2 * 60 * 1000);
  if (!isNinaApiSuccess(result)) {
    seqLog(`Guiding start note: ${JSON.stringify(result.body)}`, "warn");
    // Continue — guiding may already be running
  }

  seqState.currentStep = "Waiting for guiding to settle...";
  seqLog("Waiting for guiding to settle...");

  // Poll until State === "Guiding"
  await pollEquipment(
    target,
    (info) => {
      const s = info?.Guider?.State;
      return s === "Guiding" || s === "LostLock";
    },
    5 * 60 * 1000,
    3000,
  );
  seqLog("Guiding active ✓");
}

async function stepAutofocus(target) {
  const now = Date.now();
  if (seqState.lastAutofocusTime && (now - seqState.lastAutofocusTime) < 3600000) {
    const minsAgo = Math.round((now - seqState.lastAutofocusTime) / 60000);
    seqLog(`Skipping autofocus (last run ${minsAgo} min ago) ✓`);
    return;
  }

  seqState.currentStep = "Running autofocus...";
  seqLog("Running autofocus...");

  // Record starting position so we can detect movement
  const { equipmentInfo: infoBefore } = await getEquipmentInfo(target);
  const startPos = infoBefore?.Focuser?.Position;
  seqLog(`  Focuser start position: ${startPos ?? "unknown"}`);

  // Trigger AF — non-blocking, returns "Autofocus started" immediately
  const result = await callNinaLong(target, "/v2/api/equipment/focuser/auto-focus", {}, 15000);

  if (!isNinaApiSuccess(result)) {
    seqLog(`  AF trigger failed: ${result.body?.Error || result.status} — skipping`, "warn");
    return;
  }
  seqLog(`  ${result.body?.Response ?? "AF started"} — monitoring position...`);

  // Give NINA time to start the AF routine (first capture + analysis)
  await new Promise(r => setTimeout(r, 15000));
  checkAbort();

  // Poll focuser position until it stabilises after movement, or timeout
  const AF_DEADLINE   = Date.now() + 10 * 60 * 1000; // 10 min hard limit
  const STABLE_NEEDED = 4;   // consecutive polls with no change → done  (4 × 5s = 20s)
  const NO_MOVE_LIMIT = 24;  // polls before giving up on movement detection (2 min)

  let lastPos     = startPos;
  let movementSeen = false;
  let stableCount  = 0;
  let pollCount    = 0;

  while (Date.now() < AF_DEADLINE) {
    checkAbort();
    await new Promise(r => setTimeout(r, 5000));
    pollCount++;

    try {
      const { equipmentInfo } = await getEquipmentInfo(target);
      const pos      = equipmentInfo?.Focuser?.Position;
      const isMoving = equipmentInfo?.Focuser?.IsMoving;

      seqState.currentStep = `Autofocus: pos ${pos ?? "?"}${isMoving ? " ●" : ""}`;

      if (pos !== lastPos) {
        movementSeen = true;
        stableCount  = 0;
        seqLog(`  Focuser moving: ${lastPos} → ${pos}`);
        lastPos = pos;
      } else {
        stableCount++;
      }

      // Seen movement + position stable for STABLE_NEEDED polls → AF done
      if (movementSeen && stableCount >= STABLE_NEEDED) {
        seqState.lastAutofocusTime = Date.now();
        seqLog(`  Autofocus complete — position: ${startPos} → ${pos} ✓`);
        return;
      }

      // No movement at all after 2 min → AF may have been instant or skipped
      if (!movementSeen && pollCount >= NO_MOVE_LIMIT) {
        seqLog(`  No focuser movement detected after ${Math.round(NO_MOVE_LIMIT * 5 / 60)} min — assuming AF complete`);
        seqState.lastAutofocusTime = Date.now();
        return;
      }
    } catch { /* ignore transient errors */ }
  }

  seqLog("  AF wait timed out (10 min) — proceeding", "warn");
  seqState.lastAutofocusTime = Date.now();
}

/** Wait the exposure duration with periodic abort checks every 5 s */
async function waitExposure(durationSec) {
  const expireAt = Date.now() + durationSec * 1000;
  while (Date.now() < expireAt) {
    checkAbort();
    const remaining = expireAt - Date.now();
    await new Promise(r => setTimeout(r, Math.min(5000, Math.max(0, remaining))));
  }
}

/** Poll until Camera.IsExposing is false (download + save complete) */
async function waitForCameraIdle(target, timeoutMs = 120000) {
  return pollEquipment(
    target,
    (info) => info?.Camera?.IsExposing === false,
    timeoutMs,
    3000,
  );
}

/**
 * Abort immediately if NINA is already capturing (e.g. user manually started
 * imaging in NINA). Call before the sequence starts and before each capture.
 */
async function assertCameraIdle(target, context = "") {
  const { equipmentInfo } = await getEquipmentInfo(target);
  const cam = equipmentInfo?.Camera;
  if (!cam?.Connected) throw new Error("Camera not connected");
  if (cam?.IsExposing === true) {
    const where = context ? ` (${context})` : "";
    throw new Error(
      `Camera is already exposing${where} — stop the ongoing capture in NINA before running the sequence.`
    );
  }
}

function findFilterByName(filters, name) {
  const n = name.toLowerCase();
  return (
    filters.find(f => f.Name === name) ||                                      // exact
    filters.find(f => (f.Name || "").toLowerCase() === n) ||                   // case-insensitive exact
    filters.find(f => (f.Name || "").toLowerCase().startsWith(n)) ||           // prefix: "G" → "GBP"
    filters.find(f => (f.Name || "").toLowerCase().endsWith(n)) ||             // suffix: "RP" → "GRP"
    filters.find(f => (f.Name || "").toLowerCase().includes(n)) ||             // substring
    null
  );
}


/**
 * Snapshot all .fit/.fits filenames currently present in every SNAPSHOT directory.
 * Call this BEFORE triggering a capture; pass the result to patchFitsFilterHeaders.
 * Using a before/after set difference avoids any dependency on clock timestamps,
 * which is important because the NAS (Windows) clock can differ from the Linux host.
 */
function snapshotNasFiles() {
  const existing = new Set();
  try {
    for (const entry of fs.readdirSync(NAS_WATCH_PATH, { withFileTypes: true })) {
      if (!entry.isDirectory()) continue;
      const snapDir = path.join(NAS_WATCH_PATH, entry.name, "SNAPSHOT");
      if (!fs.existsSync(snapDir)) continue;
      try {
        for (const name of fs.readdirSync(snapDir)) {
          if (/\.(fit|fits)$/i.test(name)) existing.add(path.join(snapDir, name));
        }
      } catch { /* ignore unreadable dirs */ }
    }
  } catch (e) {
    console.error(`[patchFitsFilter] snapshotNasFiles failed: ${e.message}`);
  }
  return existing;
}

/**
 * After NINA saves a frame, find the new FITS file (by set-difference with beforeFiles)
 * and patch its FILTER header to filterName.
 *
 * Using before/after set difference avoids mtime clock-skew issues between Linux host
 * and the Windows machine where NINA saves files to the NAS.
 *
 * Retries up to ~15 s to handle NFS write-flush delays.
 * Returns the number of files patched.
 */
async function patchFitsFilterHeaders(beforeFiles, filterName) {
  const sleep = ms => new Promise(r => setTimeout(r, ms));
  const patched = new Set();

  // Collect all .../DATE/SNAPSHOT/ directories
  const snapshotDirs = [];
  try {
    for (const entry of fs.readdirSync(NAS_WATCH_PATH, { withFileTypes: true })) {
      if (!entry.isDirectory()) continue;
      const snapDir = path.join(NAS_WATCH_PATH, entry.name, "SNAPSHOT");
      if (fs.existsSync(snapDir)) snapshotDirs.push(snapDir);
    }
  } catch (e) {
    console.error(`[patchFitsFilter] cannot read NAS dir: ${e.message}`);
    return 0;
  }

  // Retry loop — NFS directory listings can be stale for several seconds
  for (let attempt = 0; attempt < 6; attempt++) {
    if (attempt > 0) await sleep(3000);

    for (const snapshotDir of snapshotDirs) {
      let names;
      try { names = fs.readdirSync(snapshotDir); } catch { continue; }

      for (const name of names) {
        if (!/\.(fit|fits)$/i.test(name)) continue;
        const fullPath = path.join(snapshotDir, name);
        if (beforeFiles.has(fullPath)) continue; // existed before capture — skip
        if (patched.has(fullPath)) continue;     // already patched — skip

        // Wait for file size to stabilise (NINA may still be flushing to NAS)
        let prevSize = -1;
        for (let s = 0; s < 8; s++) {
          const sz = (() => { try { return fs.statSync(fullPath).size; } catch { return 0; } })();
          if (sz > 0 && sz === prevSize) break;
          prevSize = sz;
          await sleep(500);
        }

        // Patch FILTER header via astropy (setval works in-place, single line)
        try {
          const cmd = `from astropy.io import fits; fits.setval(${JSON.stringify(fullPath)}, 'FILTER', value=${JSON.stringify(filterName)})`;
          execSync(`.venv/bin/python -c "${cmd.replace(/"/g, '\\"')}"`,
            { cwd: path.join(__dirname, ".."), timeout: 15000 });

          // Rename the file: replace the filter token in the filename
          // Format: YYYY-MM-DD_HH-MM-SS_FILTER_TEMP_DURATION_SEQ.ext
          const ext  = path.extname(name);
          const base = path.basename(name, ext);
          const parts = base.split("_");
          // parts[2] is the filter token written by NINA (e.g. "RP")
          let newName = name;
          if (parts.length >= 3 && parts[2] !== filterName) {
            parts[2] = filterName;
            newName = parts.join("_") + ext;
            const newPath = path.join(snapshotDir, newName);
            fs.renameSync(fullPath, newPath);
            console.log(`[patchFitsFilter] renamed: ${name} → ${newName}`);
          } else {
            console.log(`[patchFitsFilter] header patched: ${name} → FILTER=${filterName}`);
          }
          patched.add(path.join(snapshotDir, newName));
        } catch (e) {
          console.error(`[patchFitsFilter] failed for ${name}: ${e.message}`);
        }
      }
    }

    if (patched.size > 0) break;
  }

  if (patched.size === 0) {
    console.warn(`[patchFitsFilter] no new files found after retries (filter=${filterName})`);
  }
  return patched.size;
}

async function stepCaptureFilter(target, targetName, filterName, count, duration, gain, filters) {
  checkAbort();
  const filter = findFilterByName(filters, filterName);
  if (!filter) {
    const available = filters.map(f => f.Name).join(", ");
    seqLog(`Filter "${filterName}" not found (available: ${available || "none"}) — skipping`, "warn");
    return;
  }

  seqLog(`Changing to filter: ${filter.Name} (ID ${filter.Id})...`);
  const changeResult = await callNinaLong(
    target, "/v2/api/equipment/filterwheel/change-filter", { filterName: filter.Name }, 30000,
  );
  if (!isNinaApiSuccess(changeResult)) {
    throw new Error(`Filter change to "${filter.Name}" failed`);
  }
  seqLog(`Filter → ${filter.Name} ✓`);
  checkAbort();

  let saved = 0;
  for (let i = 1; i <= count; i++) {
    checkAbort();
    seqState.currentStep = `${filterName}: frame ${i}/${count}`;
    seqState.progress = { filter: filterName, frame: i, frames: count };

    seqLog(`${filterName} frame ${i}/${count}: triggering ${duration}s exposure (gain ${gain})`);

    const capturePayload = { duration, gain, save: true, targetName, filter: filterName };
    const beforeFiles = snapshotNasFiles(); // snapshot before capture for set-difference check

    // ninaAPI capture is non-blocking — it starts the exposure and returns immediately
    let captureResult;
    try {
      captureResult = await callNinaLong(
        target, "/v2/api/equipment/camera/capture", capturePayload, 30000,
      );
    } catch (fetchErr) {
      // Network-level error (fetch failed, ECONNREFUSED, ETIMEDOUT) — retry once after 15s
      seqLog(`${filterName} ${i}/${count}: network error (${fetchErr.message}) — retrying in 15s...`, "warn");
      await new Promise(r => setTimeout(r, 15000));
      try {
        captureResult = await callNinaLong(
          target, "/v2/api/equipment/camera/capture", capturePayload, 30000,
        );
      } catch (retryErr) {
        seqLog(`${filterName} ${i}/${count}: retry also failed (${retryErr.message}) — skipping frame`, "warn");
        continue;
      }
    }

    if (!isNinaApiSuccess(captureResult)) {
      seqLog(`${filterName} ${i}/${count}: failed to start — ${captureResult.body?.Error || captureResult.status}`, "warn");
      // If already exposing, wait for it to finish then retry once
      if (captureResult.status === 409) {
        seqLog(`Camera busy — waiting for current exposure to finish...`);
        try { await waitForCameraIdle(target, (duration + 120) * 1000); } catch { /* timeout */ }
        const retry = await callNinaLong(
          target, "/v2/api/equipment/camera/capture", capturePayload, 30000,
        );
        if (!isNinaApiSuccess(retry)) {
          seqLog(`${filterName} ${i}/${count}: retry failed — skipping frame`, "warn");
          continue;
        }
      } else {
        continue;
      }
    }

    // Wait for the exposure to run (abort-checkable sleep in 5 s increments)
    seqState.currentStep = `${filterName}: exposing ${i}/${count} (${duration}s)`;
    await waitExposure(duration);
    checkAbort();

    // Wait for download + save to complete (Camera.IsExposing → false)
    seqState.currentStep = `${filterName}: downloading ${i}/${count}`;
    try {
      await waitForCameraIdle(target, 120 * 1000);
    } catch {
      // If camera doesn't report state, just add a fixed buffer
      await new Promise(r => setTimeout(r, 20000));
    }

    // Patch any new FITS files with the correct FILTER header (await — it retries on NFS delay)
    await patchFitsFilterHeaders(beforeFiles, filter.Name);

    saved++;
    seqLog(`${filterName} ${i}/${count} saved ✓`);
  }
  seqLog(`${filterName} complete (${saved}/${count}) ✓`);
}

// ── Main sequence runner ──────────────────────────────────────────────────────

async function runSequence(ninaConfig, seqConfig) {
  const { duration, gain, count, filters: filterNames, solveEnabled } = seqConfig;

  seqState.running    = true;
  seqState.aborted    = false;
  seqState.error      = null;
  seqState.manualMode = seqConfig.manualMode === true;
  seqState.ninaConfig = ninaConfig;

  seqLog(`=== Sequence started ${seqState.manualMode ? "(MANUAL step mode)" : "(auto)"} ===`);

  /**
   * Home the mount using FindHome (physical sensor search).
   * The NINA endpoint is non-blocking — we poll until AtHome=true.
   * If AtHome is already true the driver skips FindHome, so we force
   * an unpark first to clear the flag, then home.
   * Never throws — always logs and returns.
   */
  async function safeHome() {
    seqState.currentStep = "Homing mount...";
    try {
      // FindHome is skipped when AtHome=true, and blocked when AtPark=true.
      // Unpark in both cases to let FindHome execute physically.
      const { equipmentInfo: preInfo } = await getEquipmentInfo(ninaConfig);
      const mountPre = preInfo?.Mount;
      if (mountPre?.AtPark === true || mountPre?.AtHome === true) {
        const reason = mountPre.AtPark ? "parked" : "already at home flag";
        seqLog(`  Mount ${reason} — unparking before FindHome...`);
        await callNinaLong(ninaConfig, "/v2/api/equipment/mount/unpark", {}, 30000);
        // Poll until actually unparked (up to 20 s)
        const deadline2 = Date.now() + 20000;
        while (Date.now() < deadline2) {
          await new Promise(r => setTimeout(r, 1500));
          const { equipmentInfo: ei2 } = await getEquipmentInfo(ninaConfig);
          if (ei2?.Mount?.AtPark === false) break;
        }
      }

      // Trigger the (non-blocking) FindHome
      seqLog("  Sending FindHome command...");
      const result = await callNinaLong(ninaConfig, "/v2/api/equipment/mount/home", {}, 15000);
      const resp = result.body?.Response ?? "";
      seqLog(`  Home command accepted: ${resp}`);

      // Poll until the mount is no longer slewing AND AtHome is true
      const HOME_TIMEOUT = 5 * 60 * 1000;
      const POLL_INTERVAL = 3000;
      const deadline = Date.now() + HOME_TIMEOUT;

      while (Date.now() < deadline) {
        await new Promise(r => setTimeout(r, POLL_INTERVAL));
        const { equipmentInfo } = await getEquipmentInfo(ninaConfig);
        const m = equipmentInfo?.Mount;
        seqState.currentStep = `Homing: slewing=${m?.Slewing}  atHome=${m?.AtHome}`;
        if (m?.AtHome === true && m?.Slewing !== true) {
          seqLog("Mount homed ✓ (FindHome complete)");
          return;
        }
      }

      seqLog("Home timed out — mount may not have reached home sensors", "warn");
    } catch (e) {
      seqLog(`Home failed: ${e.message}`, "warn");
    }
  }

  try {
    // ── Pre-flight: ensure NINA is not already capturing ─────────────────────
    await assertCameraIdle(ninaConfig, "sequence start");

    // ── Step 1: Cool camera ──────────────────────────────────────────────────
    try {
      await stepCoolCamera(ninaConfig, -5);
    } catch (err) {
      if (err.message === "__ABORTED__") throw err;
      seqLog(`✗ Camera cooling failed: ${err.message} — aborting sequence`, "error");
      await safeHome();
      return;                      // stop sequence, mount is safe
    }

    // ── Step 1b: Filter checklist — verify all requested filters exist ──────────
    let verifiedFilters = [];
    try {
      const { equipmentInfo: fwInfo } = await getEquipmentInfo(ninaConfig);
      const availableFilters = fwInfo?.FilterWheel?.AvailableFilters || [];
      const availableNames = availableFilters.map(f => f.Name).filter(Boolean);

      seqLog("── Filter checklist ────────────────────────────────");
      seqLog(`  Filterwheel reports: ${availableNames.join(", ") || "(none)"}`);
      seqLog(`  Requested filters  : ${filterNames.join(", ")}`);
      seqLog("");

      const missing = [];
      for (const name of filterNames) {
        const found = findFilterByName(availableFilters, name);
        if (found) {
          seqLog(`  ✓  ${name.padEnd(6)} → matched as "${found.Name}" (ID ${found.Id})`);
          verifiedFilters.push(name);
        } else {
          seqLog(`  ✗  ${name.padEnd(6)} → NOT FOUND in filterwheel — will be skipped`);
          missing.push(name);
        }
      }

      if (missing.length) {
        seqLog(`  ⚠️  ${missing.length} filter(s) missing — proceeding with: ${verifiedFilters.join(", ") || "none"}`);
      } else {
        seqLog(`  All ${verifiedFilters.length} filter(s) OK ✓`);
      }
      seqLog("────────────────────────────────────────────────────");

      if (!verifiedFilters.length) {
        seqLog("✗ No valid filters found — aborting sequence", "error");
        await safeHome();
        return;
      }

      await waitForConfirmation("Filter checklist done — proceed to imaging?");
      checkAbort();
    } catch (err) {
      if (err.message === "__ABORTED__") throw err;
      seqLog(`Filter checklist failed: ${err.message} — continuing with requested filters`, "warn");
      verifiedFilters = [...filterNames];
    }

    // ── Step 1c: Switch to G filter (clear dark/lum filter before imaging) ──────
    try {
      const { equipmentInfo: fwInfo } = await getEquipmentInfo(ninaConfig);
      const filters = fwInfo?.FilterWheel?.AvailableFilters || [];
      const gFilter = findFilterByName(filters, "G");
      if (gFilter) {
        seqLog(`Switching to G filter (ID ${gFilter.Id}) to clear any dark filter...`);
        const r = await callNinaLong(ninaConfig, "/v2/api/equipment/filterwheel/change-filter",
          { filterName: gFilter.Name }, 30000);
        if (isNinaApiSuccess(r)) seqLog("Filter → G ✓");
        else seqLog("Filter change to G returned non-success (continuing)", "warn");
      } else {
        const available = filters.map(f => f.Name).join(", ");
        seqLog(`G filter not found in filterwheel (available: ${available || "none"}) — skipping pre-filter switch`, "warn");
      }
    } catch (err) {
      if (err.message === "__ABORTED__") throw err;
      seqLog(`Filter switch warning: ${err.message} — continuing anyway`, "warn");
    }

    await waitForConfirmation("Camera cooled — proceed to first target?");
    checkAbort();

    const pending = seqState.queue.filter(t => !t.done);
    const total = pending.length;
    let done = 0;

    for (let i = 0; i < seqState.queue.length; i++) {
      const t = seqState.queue[i];
      if (t.done) continue;
      checkAbort();

      done++;
      seqState.currentTargetIdx = i;
      seqState.currentTarget = `${t.name} (${done}/${total})`;
      seqLog(`--- Target ${done}/${total}: ${t.name} ---`);

      try {
        // Fetch available filters before each target
        const { equipmentInfo } = await getEquipmentInfo(ninaConfig);
        const filters = equipmentInfo?.FilterWheel?.AvailableFilters || [];
        if (!filters.length) seqLog("No filters found in filterwheel", "warn");

        // ── Step 2: Unpark + enable tracking ────────────────────────────────
        await stepUnparkMount(ninaConfig);
        const { equipmentInfo: mInfo } = await getEquipmentInfo(ninaConfig);
        const mPos = mInfo?.Mount;
        const mDesc = mPos
          ? `RA ${mPos.RightAscensionString ?? "?"}  Dec ${mPos.DeclinationString ?? "?"}  Tracking: ${mPos.TrackingEnabled ? mPos.TrackingMode : "OFF"}`
          : "unknown position";
        await waitForConfirmation(`Mount ready (${mDesc}) — proceed to slew to ${t.name}?`);
        checkAbort();

        // ── Step 3: Horizon check + Slew (+ optional plate-solve centering) ─
        await checkTargetAboveHorizon(ninaConfig, t.raDeg, t.decDeg, t.name);
        const useCenter = solveEnabled !== false;
        await stepSlewToTarget(ninaConfig, t.raDeg, t.decDeg, t.name, { center: useCenter });
        if (!useCenter) seqLog("Plate solve centering disabled — skipping");
        const slewDoneMsg = useCenter
          ? `Slew + centering complete — proceed to start guiding?`
          : `Slew complete — proceed to start guiding?`;
        await waitForConfirmation(slewDoneMsg);
        checkAbort();

        // ── Step 4: Guiding ──────────────────────────────────────────────────
        await stepStartGuiding(ninaConfig);
        await waitForConfirmation(`Guiding active — proceed to autofocus?`);
        checkAbort();

        // ── Step 5: Autofocus ────────────────────────────────────────────────
        await stepAutofocus(ninaConfig);
        await waitForConfirmation(`Autofocus done — proceed to captures?`);
        checkAbort();

        // ── Step 6: Capture each filter ──────────────────────────────────────
        for (const filterName of verifiedFilters) {
          checkAbort();
          await assertCameraIdle(ninaConfig, `before ${filterName} captures`);
          await stepCaptureFilter(ninaConfig, t.name, filterName, count, duration, gain, filters);
          await waitForConfirmation(`${filterName} captures complete — proceed to next filter?`);
          checkAbort();
        }

        t.done = true;
        saveQueue();
        seqLog(`${t.name} — all captures complete ✓`);
        await waitForConfirmation(total > 1 ? `${t.name} done — proceed to next target?` : `${t.name} done — home mount?`);
        checkAbort();

      } catch (err) {
        // User abort propagates up to the outer handler — stop the whole sequence
        if (err.message === "__ABORTED__") throw err;

        const isTransient = /fetch failed|ECONNREFUSED|ETIMEDOUT|ENOTFOUND/i.test(err.message);
        const MAX_TARGET_RETRIES = 2;
        let targetSucceeded = false;

        if (isTransient) {
          for (let attempt = 1; attempt <= MAX_TARGET_RETRIES && !targetSucceeded; attempt++) {
            seqLog(`✗ ${t.name} transient error (${err.message}) — retry ${attempt}/${MAX_TARGET_RETRIES} in 30s...`, "warn");
            await new Promise(r => setTimeout(r, 30000));
            try {
              // Re-run the full target block after transient error
              const { equipmentInfo: reInfo } = await getEquipmentInfo(ninaConfig);
              const reFilters = reInfo?.FilterWheel?.AvailableFilters || [];
              await stepUnparkMount(ninaConfig);
              await checkTargetAboveHorizon(ninaConfig, t.raDeg, t.decDeg, t.name);
              const useCenter = solveEnabled !== false;
              await stepSlewToTarget(ninaConfig, t.raDeg, t.decDeg, t.name, { center: useCenter });
              await stepStartGuiding(ninaConfig);
              await stepAutofocus(ninaConfig);
              for (const filterName of verifiedFilters) {
                checkAbort();
                await assertCameraIdle(ninaConfig, `before ${filterName} captures`);
                await stepCaptureFilter(ninaConfig, t.name, filterName, count, duration, gain, reFilters);
              }
              t.done = true;
              saveQueue();
              seqLog(`${t.name} — recovered and complete ✓`);
              targetSucceeded = true;
            } catch (retryErr) {
              if (retryErr.message === "__ABORTED__") throw retryErr;
              err = retryErr; // update err for next iteration / final log
            }
          }
        }

        if (!targetSucceeded) {
          t.error = err.message;
          seqLog(`✗ ${t.name} failed: ${err.message}`, "error");
          seqLog(`  Homing mount before next target...`);
          await safeHome();
          seqLog(`  Continuing to next target`);
        }
      }
    }

    // ── Final: Home mount ──────────────────────────────────────────────────
    seqLog("All targets complete. Homing mount...");
    await safeHome();
    seqLog("=== Sequence complete ===");

  } catch (err) {
    if (err.message === "__ABORTED__") {
      seqLog("=== Sequence aborted by user ===", "warn");
    } else {
      seqLog(`=== Sequence error: ${err.message} ===`, "error");
      seqState.error = err.message;
    }
  } finally {
    seqState.running = false;
    seqState.ninaConfig = null;
    seqState.currentTarget = null;
    seqState.currentStep = null;
    seqState.progress = null;
    seqState.currentTargetIdx = -1;
    seqState.waitingForStep = null;
    // Release any pending confirmation so the promise doesn't leak
    if (_confirmResolve) { _confirmResolve(); _confirmResolve = null; _confirmReject = null; }
  }
}

// ── Sequence routes ───────────────────────────────────────────────────────────

app.get("/api/sequence/state", (req, res) => {
  res.json({
    running: seqState.running,
    manualMode: seqState.manualMode,
    waitingForStep: seqState.waitingForStep,
    queue: seqState.queue,
    currentTargetIdx: seqState.currentTargetIdx,
    currentTarget: seqState.currentTarget,
    currentStep: seqState.currentStep,
    progress: seqState.progress,
    log: seqState.log,
    lastAutofocusTime: seqState.lastAutofocusTime,
    error: seqState.error,
  });
});

app.post("/api/sequence/queue", (req, res) => {
  const name = String(req.body?.name || "").trim();
  const raDeg = Number(req.body?.raDeg);
  const decDeg = Number(req.body?.decDeg);
  if (!name || !Number.isFinite(raDeg) || !Number.isFinite(decDeg)) {
    return res.status(400).json({ success: false, error: "name, raDeg, decDeg required" });
  }
  seqState.queue.push({
    name,
    raDeg,
    decDeg,
    ra: String(req.body?.ra || ""),
    dec: String(req.body?.dec || ""),
    done: false,
    addedAt: Date.now(),
  });
  saveQueue();
  return res.json({ success: true, queue: seqState.queue });
});

app.delete("/api/sequence/queue/:index", (req, res) => {
  const idx = Number(req.params.index);
  if (!Number.isInteger(idx) || idx < 0 || idx >= seqState.queue.length) {
    return res.status(400).json({ success: false, error: "Invalid queue index" });
  }
  if (seqState.running && idx === seqState.currentTargetIdx) {
    return res.status(409).json({ success: false, error: "Cannot remove currently active target" });
  }
  seqState.queue.splice(idx, 1);
  saveQueue();
  return res.json({ success: true, queue: seqState.queue });
});

app.post("/api/sequence/run", (req, res) => {
  if (seqState.running) {
    return res.status(409).json({ success: false, error: "Sequence already running" });
  }
  const pending = seqState.queue.filter(t => !t.done);
  if (!pending.length) {
    return res.status(400).json({ success: false, error: "No pending targets in queue" });
  }

  let ninaConfig;
  try { ninaConfig = normalizeTargetConfig(req.body); }
  catch (err) { return res.status(400).json({ success: false, error: err.message }); }

  const seqConfig = {
    duration:       Number(req.body?.duration)       || 120,
    gain:           Number(req.body?.gain)           || 10,
    count:          Number(req.body?.count)          || 10,
    filters:        Array.isArray(req.body?.filters) ? req.body.filters : ["G", "BP", "RP"],  // real filterwheel names
    solveEnabled:   req.body?.solveEnabled !== false,
    solveExp:       Number(req.body?.solveExp)       || 5,
    solveThreshold: Number(req.body?.solveThreshold) || 60,
    manualMode:     req.body?.manualMode === true,
  };

  // Fire-and-forget — do not await
  runSequence(ninaConfig, seqConfig).catch(err => {
    console.error("Sequence runner uncaught error:", err);
  });

  return res.json({ success: true, message: "Sequence started", pendingTargets: pending.length });
});

app.post("/api/sequence/abort", (req, res) => {
  if (!seqState.running) {
    return res.status(400).json({ success: false, error: "No sequence running" });
  }
  seqState.aborted = true;

  // If paused waiting for confirmation, reject the promise so it throws __ABORTED__
  if (_confirmReject) {
    _confirmReject(new Error("__ABORTED__"));
    _confirmResolve = null;
    _confirmReject  = null;
  }

  // Tell NINA to stop the current exposure and guider immediately (fire-and-forget)
  if (seqState.ninaConfig) {
    const cfg = seqState.ninaConfig;
    callNinaLong(cfg, "/v2/api/equipment/camera/abort", {}, 8000)
      .then(() => seqLog("Camera exposure aborted ✓"))
      .catch(() => {/* non-fatal */});
    callNinaLong(cfg, "/v2/api/equipment/guider/stop", {}, 8000)
      .then(() => seqLog("Guider stopped ✓"))
      .catch(() => {/* non-fatal */});
  }

  return res.json({ success: true, message: "Abort signal sent — stopping camera and guider" });
});

app.post("/api/sequence/next", (req, res) => {
  if (!seqState.waitingForStep || !_confirmResolve) {
    return res.status(400).json({ success: false, error: "Sequence is not waiting for confirmation" });
  }
  const step = seqState.waitingForStep;
  _confirmResolve();
  _confirmResolve = null;
  _confirmReject  = null;
  return res.json({ success: true, confirmedStep: step });
});

app.post("/api/sequence/manual-mode", (req, res) => {
  const enabled = req.body?.enabled === true;
  seqState.manualMode = enabled;
  seqLog(`Manual step mode ${enabled ? "ON" : "OFF"}`, "info");
  // If switching OFF while paused, auto-resolve so the sequence continues
  if (!enabled && seqState.waitingForStep && _confirmResolve) {
    seqLog(`▶ Auto-continuing (manual mode turned off)`, "info");
    _confirmResolve();
    _confirmResolve = null;
    _confirmReject  = null;
    seqState.waitingForStep = null;
  }
  return res.json({ success: true, manualMode: enabled });
});

app.post("/api/sequence/clear", (req, res) => {
  if (seqState.running) {
    return res.status(409).json({ success: false, error: "Cannot clear while sequence is running" });
  }
  seqState.queue = [];
  seqState.log = [];
  seqState.error = null;
  seqState.currentTargetIdx = -1;
  saveQueue();
  return res.json({ success: true });
});

app.post("/api/sequence/reset-af", (req, res) => {
  seqState.lastAutofocusTime = null;
  return res.json({ success: true, message: "Autofocus timer reset — AF will run on next sequence" });
});

app.post("/api/sequence/restart", (req, res) => {
  if (seqState.running) {
    return res.status(409).json({ success: false, error: "Cannot restart while sequence is running" });
  }
  seqState.queue.forEach(t => { t.done = false; delete t.error; });
  seqState.log = [];
  seqState.error = null;
  seqState.currentTargetIdx = -1;
  saveQueue();
  return res.json({ success: true, queue: seqState.queue });
});

app.post("/api/sequence/reorder", (req, res) => {
  if (seqState.running) {
    return res.status(409).json({ success: false, error: "Cannot reorder while sequence is running" });
  }
  const order = req.body?.order; // array of target names in the desired order
  if (!Array.isArray(order)) {
    return res.status(400).json({ success: false, error: "order array of names required" });
  }
  const nameIndex = new Map(order.map((name, i) => [name, i]));
  seqState.queue.sort((a, b) => {
    const ia = nameIndex.has(a.name) ? nameIndex.get(a.name) : Infinity;
    const ib = nameIndex.has(b.name) ? nameIndex.get(b.name) : Infinity;
    return ia - ib;
  });
  saveQueue();
  return res.json({ success: true, queue: seqState.queue });
});

// ── TNS (Transient Name Server) ───────────────────────────────────────────────

function parseTnsName(raw) {
  // Accept "AT 2026fuh", "SN 2026fuh", "2026fuh" → returns "2026fuh"
  return String(raw || "").trim().replace(/^(AT|SN|FRB)\s+/i, "").trim();
}

app.post("/api/target/tns", async (req, res) => {
  const botId = String(req.body?.botId || TNS_BOT_ID).trim();
  const botName = String(req.body?.botName || TNS_BOT_NAME).trim();
  const apiKey = String(req.body?.apiKey || TNS_API_KEY).trim();
  const rawName = String(req.body?.name || "").trim();

  if (!botId || !botName || !apiKey) {
    return res.status(400).json({
      success: false,
      error: "TNS credentials required: botId, botName, apiKey",
    });
  }

  const objname = parseTnsName(rawName);
  if (!objname) {
    return res.status(400).json({ success: false, error: "Target name is required" });
  }

  const tnsMarker = JSON.stringify({ tns_id: Number(botId), type: "bot", name: botName });
  const queryData = JSON.stringify({ objname, photometry: "0", spectra: "0" });

  const body = new URLSearchParams({ api_key: apiKey, data: queryData });

  const controller = new AbortController();
  const timeout = setTimeout(() => controller.abort(), 15000);

  try {
    const response = await fetch("https://www.wis-tns.org/api/get/object", {
      method: "POST",
      headers: {
        "User-Agent": `tns_marker${tnsMarker}`,
        "Content-Type": "application/x-www-form-urlencoded",
      },
      body: body.toString(),
      signal: controller.signal,
    });

    const json = await response.json();

    // TNS API returns { id_code, id_message, data: { objname, ... } }
    const obj = json?.data;

    if (!response.ok || !obj?.objname) {
      return res.status(response.ok ? 404 : response.status).json({
        success: false,
        error: json?.id_message || `TNS returned status ${response.status}`,
        raw: json,
      });
    }

    const prefix = obj.name_prefix || "AT";
    return res.json({
      success: true,
      target: {
        name: obj.objname,
        prefix,
        fullName: `${prefix} ${obj.objname}`,
        ra: obj.ra ?? null,
        dec: obj.dec ?? null,
        raDeg: obj.radeg ?? null,
        decDeg: obj.decdeg ?? null,
        redshift: obj.redshift ?? null,
        type: obj.object_type?.name ?? null,
        discoveryDate: obj.discoverydate ?? null,
        discoveryMag: obj.discoverymag ?? null,
        discoveryFilter: obj.discmagfilter?.name ?? null,
        hostName: obj.hostname ?? null,
        hostRedshift: obj.host_redshift ?? null,
        reporters: obj.discoverer ?? obj.reporter ?? null,
        reportingGroup: obj.reporting_group?.group_name ?? null,
        sourceGroup: obj.discovery_data_source?.group_name ?? null,
        internalNames: obj.internal_names ?? null,
        tnsUrl: `https://www.wis-tns.org/object/${obj.objname}`,
      },
      raw: obj,
    });
  } catch (error) {
    return res.status(502).json({
      success: false,
      error: error.name === "AbortError" ? "TNS request timed out" : error.message,
    });
  } finally {
    clearTimeout(timeout);
  }
});

// ── AstroColibri ─────────────────────────────────────────────────────────────
// GET /api/target/astrocolibri?name=2025mq
// Queries https://astro-colibri.science/event?trigger_id=TNS{name}
// No authentication needed.

app.get("/api/target/astrocolibri", async (req, res) => {
  const rawName = String(req.query.name || "").trim();
  const objname = parseTnsName(rawName);  // strip "AT ", "SN " prefix
  if (!objname) {
    return res.status(400).json({ success: false, error: "Target name is required" });
  }

  const triggerId = `TNS${objname}`;
  const url = `https://astro-colibri.science/event?trigger_id=${encodeURIComponent(triggerId)}`;
  const controller = new AbortController();
  const timeout = setTimeout(() => controller.abort(), 12000);

  try {
    const response = await fetch(url, { signal: controller.signal });
    const json = await response.json();

    // Response is either an object (found) or [{error: "event not found", ...}]
    const data = Array.isArray(json) ? json[0] : json;

    if (data?.error || !data?.ra) {
      return res.status(404).json({
        success: false,
        error: data?.error || "Not found on AstroColibri",
        source: "astrocolibri",
      });
    }

    // Normalise ra/dec to sexagesimal strings for display
    const raHMS  = raDecToHMS(data.ra);
    const decDMS = raDecToDMS(data.dec);

    return res.json({
      success: true,
      source: "astrocolibri",
      target: {
        name:        objname,
        prefix:      data.source_name?.split(" ")[0] || "AT",
        fullName:    data.source_name || `AT ${objname}`,
        ra:          raHMS,
        dec:         decDMS,
        raDeg:       data.ra,
        decDeg:      data.dec,
        redshift:    data.redshift > 0 ? data.redshift : null,
        type:        data.classification || null,
        discoveryDate: data.timestamp
          ? new Date(data.timestamp).toISOString().slice(0, 10) : null,
        constellation: data.constellation || null,
        hostName:    null,
        reporters:   null,
        astrocolibriId: data.astro_colibri_id || null,
        colibriUrl:  `https://astro-colibri.science/?trigger_id=${triggerId}`,
        tnsUrl:      `https://www.wis-tns.org/object/${objname}`,
      },
    });
  } catch (err) {
    return res.status(502).json({
      success: false,
      error: err.name === "AbortError" ? "AstroColibri request timed out" : err.message,
      source: "astrocolibri",
    });
  } finally {
    clearTimeout(timeout);
  }
});

// Helpers: decimal degrees → sexagesimal strings
function raDecToHMS(ra) {
  const h = Math.floor(ra / 15);
  const rem = (ra / 15 - h) * 60;
  const m = Math.floor(rem);
  const s = ((rem - m) * 60).toFixed(2);
  return `${String(h).padStart(2,"0")}:${String(m).padStart(2,"0")}:${String(s).padStart(5,"0")}`;
}
function raDecToDMS(dec) {
  const sign = dec < 0 ? "-" : "+";
  const abs  = Math.abs(dec);
  const d    = Math.floor(abs);
  const rem  = (abs - d) * 60;
  const m    = Math.floor(rem);
  const s    = ((rem - m) * 60).toFixed(1);
  return `${sign}${String(d).padStart(2,"0")}:${String(m).padStart(2,"0")}:${String(s).padStart(4,"0")}`;
}

// ── OCS (Observatory Control System) ─────────────────────────────────────────

function validateOcsHost(host) {
  if (!host || typeof host !== "string" || !host.trim()) {
    throw new Error("OCS host is required");
  }
  return host.trim();
}

async function ocsGet(host, params = {}) {
  const url = new URL(`http://${host}/index-ajax-get.txt`);
  for (const [k, v] of Object.entries(params)) {
    url.searchParams.set(k, v);
  }
  url.searchParams.set("x", Date.now());

  const controller = new AbortController();
  // OCS command endpoints accept the request and return an empty body slowly;
  // we abort after 3s — the command has already been received by then.
  const timeout = setTimeout(() => controller.abort(), 3000);
  try {
    const response = await fetch(url.toString(), { signal: controller.signal });
    let text = "";
    try { text = await response.text(); } catch { /* empty body is fine */ }
    return { ok: true, status: response.status, text, url: url.toString() };
  } catch (error) {
    if (error.name === "AbortError") {
      // Command was sent, OCS just didn't respond in time — that's expected
      return { ok: true, status: 202, text: "", url: url.toString() };
    }
    throw error;
  } finally {
    clearTimeout(timeout);
  }
}

async function ocsStatus(host) {
  const controller = new AbortController();
  const timeout = setTimeout(() => controller.abort(), 8000);
  try {
    const response = await fetch(`http://${host}/index.txt?x=${Date.now()}`, {
      signal: controller.signal,
    });
    const text = await response.text();
    const fields = {};
    for (const line of text.split("\n")) {
      const sep = line.indexOf("|");
      if (sep !== -1) {
        fields[line.slice(0, sep).trim()] = line.slice(sep + 1).trim();
      }
    }
    return { ok: response.ok, status: response.status, fields };
  } finally {
    clearTimeout(timeout);
  }
}

app.post("/api/ocs/status", async (req, res) => {
  let host;
  try {
    host = validateOcsHost(req.body?.ocsHost || DEFAULT_OCS_HOST);
  } catch (error) {
    return res.status(400).json({ success: false, error: error.message });
  }

  try {
    const { ok, status, fields } = await ocsStatus(host);
    return res.status(ok ? 200 : 502).json({
      success: ok,
      host,
      status,
      roof: {
        status: fields.roof_sta || "Unknown",
        error: (fields.roof_err || "").replace(/<[^>]*>/g, "").trim(),
        safe: fields.stat_safe || "Unknown",
        rain: fields.wea_rain || "Unknown",
        temp: fields.wea_temp || "Invalid",
        humidity: fields.wea_humd || "Invalid",
        pressure: fields.wea_pres || "Invalid",
      },
      raw: fields,
    });
  } catch (error) {
    return res.status(502).json({
      success: false,
      host,
      error: error.name === "AbortError" ? "Connection timed out" : error.message,
    });
  }
});

app.post("/api/ocs/roof", async (req, res) => {
  let host;
  try {
    host = validateOcsHost(req.body?.ocsHost || DEFAULT_OCS_HOST);
  } catch (error) {
    return res.status(400).json({ success: false, error: error.message });
  }

  const command = String(req.body?.command || "").trim().toLowerCase();
  if (!["open", "close", "stop"].includes(command)) {
    return res.status(400).json({
      success: false,
      error: "Command must be one of: open, close, stop",
    });
  }

  try {
    const result = await ocsGet(host, { roof: command });
    const { fields } = await ocsStatus(host);
    return res.status(result.ok ? 200 : 502).json({
      success: result.ok,
      host,
      command,
      roof: {
        status: fields.roof_sta || "Unknown",
        error: (fields.roof_err || "").replace(/<[^>]*>/g, "").trim(),
        safe: fields.stat_safe || "Unknown",
      },
    });
  } catch (error) {
    return res.status(502).json({
      success: false,
      host,
      error: error.name === "AbortError" ? "Connection timed out" : error.message,
    });
  }
});

app.get("/api/ocs/history", (req, res) => {
  const hours = Math.min(Number(req.query.hours) || 48, 168); // max 7 days
  const rows = db.prepare(
    `SELECT id, ts, roof, safe, rain, temp, humidity, pressure, sky
     FROM ocs_history
     WHERE ts >= datetime('now', ? || ' hours')
     ORDER BY ts ASC`
  ).all(`-${hours}`);
  res.json({ success: true, history: rows, hours });
});

app.post("/api/ocs/poll-now", async (req, res) => {
  try {
    const host = validateOcsHost(req.body?.ocsHost || DEFAULT_OCS_HOST);
    const { ok, fields } = await ocsStatus(host);
    if (!ok) return res.status(502).json({ success: false, error: "OCS returned non-OK status" });

    const temp     = parseFloat(fields.wea_temp);
    const humidity = parseFloat(fields.wea_humd);
    const pressure = parseFloat(fields.wea_pres);
    const sky      = parseFloat(fields.wea_sq || "");

    const result = db.prepare(
      `INSERT INTO ocs_history (roof, safe, rain, temp, humidity, pressure, sky) VALUES (?, ?, ?, ?, ?, ?, ?)`
    ).run(
      ocsClean(fields.roof_sta)  || null,
      ocsClean(fields.stat_safe) || null,
      ocsClean(fields.wea_rain)  || null,
      isNaN(temp)     ? null : temp,
      isNaN(humidity) ? null : humidity,
      isNaN(pressure) ? null : pressure,
      isNaN(sky)      ? null : sky,
    );

    const row = db.prepare("SELECT * FROM ocs_history WHERE id = ?").get(result.lastInsertRowid);
    res.json({ success: true, row });
  } catch (error) {
    res.status(502).json({
      success: false,
      error: error.name === "AbortError" ? "OCS connection timed out" : error.message,
    });
  }
});

// ── OCS backfill: scrape 60-min history from OCS embedded charts ──────────────

app.post("/api/ocs/backfill", async (req, res) => {
  try {
    const host = validateOcsHost(req.body?.ocsHost || DEFAULT_OCS_HOST);
    const now = Date.now();

    // Fetch all three time windows + sky page in parallel
    const [weatherRecent, weather24h, weather48h, skyRecent, sky24h, sky48h] = await Promise.all([
      fetch(`http://${host}/weatherpage.htm?x=${now}`,             { signal: AbortSignal.timeout(8000) }).then(r => r.text()),
      fetch(`http://${host}/weatherpage.htm?chart=last24&x=${now}`,{ signal: AbortSignal.timeout(8000) }).then(r => r.text()),
      fetch(`http://${host}/weatherpage.htm?chart=last48&x=${now}`,{ signal: AbortSignal.timeout(8000) }).then(r => r.text()),
      fetch(`http://${host}/skypage.htm?x=${now}`,                 { signal: AbortSignal.timeout(8000) }).then(r => r.text()),
      fetch(`http://${host}/skypage.htm?chart=last24&x=${now}`,    { signal: AbortSignal.timeout(8000) }).then(r => r.text()),
      fetch(`http://${host}/skypage.htm?chart=last48&x=${now}`,    { signal: AbortSignal.timeout(8000) }).then(r => r.text()),
    ]);

    // Parse {x: hoursAgo_or_minutesAgo, y: value} from an embedded Chart.js dataset.
    // Returns array of { msAgo, value } — x unit auto-detected by label suffix.
    function parseDataset(html, labelFragment) {
      const re = new RegExp(
        `label:\\s*'[^']*${labelFragment}[^']*'[\\s\\S]*?data:\\s*(\\[[^\\]]*\\])`,
        "i"
      );
      const m = html.match(re);
      if (!m) return [];
      // Detect unit from the label text
      const labelMatch = html.match(new RegExp(`label:\\s*'([^']*${labelFragment}[^']*)'`, "i"));
      const isHours = labelMatch && /hours/i.test(labelMatch[1]);
      return [...m[1].matchAll(/\{x:([\d.]+),y:([\d.-]+)\}/g)].map(p => ({
        msAgo: parseFloat(p[1]) * (isHours ? 3600000 : 60000),
        value: parseFloat(p[2]),
      }));
    }

    // Collect all points from all windows
    const allPts = {
      temp:     [...parseDataset(weatherRecent, "Temperature"), ...parseDataset(weather24h, "Temperature"), ...parseDataset(weather48h, "Temperature")],
      pressure: [...parseDataset(weatherRecent, "Pressure"),    ...parseDataset(weather24h, "Pressure"),    ...parseDataset(weather48h, "Pressure")],
      humidity: [...parseDataset(weatherRecent, "Humidity"),    ...parseDataset(weather24h, "Humidity"),    ...parseDataset(weather48h, "Humidity")],
      sky:      [...parseDataset(skyRecent,     "Sky Quality"), ...parseDataset(sky24h,     "Sky Quality"), ...parseDataset(sky48h,     "Sky Quality")],
    };

    // Build map keyed by minute-rounded timestamp string
    const byTs = new Map();
    const roundToMin = ms => {
      const d = new Date(now - ms);
      d.setSeconds(0, 0);
      return d.toISOString().replace("T", " ").slice(0, 19);
    };
    for (const [key, pts] of Object.entries(allPts)) {
      for (const { msAgo, value } of pts) {
        const ts = roundToMin(msAgo);
        if (!byTs.has(ts)) byTs.set(ts, {});
        // Prefer finer-resolution (recent) data — only fill if not yet set
        if (byTs.get(ts)[key] === undefined) byTs.get(ts)[key] = value;
      }
    }

    db.exec(`CREATE UNIQUE INDEX IF NOT EXISTS ocs_history_ts_unique ON ocs_history(ts)`);
    const stmt = db.prepare(
      `INSERT OR IGNORE INTO ocs_history (ts, temp, humidity, pressure, sky) VALUES (?, ?, ?, ?, ?)`
    );

    let inserted = 0;
    db.transaction(() => {
      for (const [ts, vals] of byTs) {
        const r = stmt.run(ts, vals.temp ?? null, vals.humidity ?? null, vals.pressure ?? null, vals.sky ?? null);
        if (r.changes) inserted++;
      }
    })();

    res.json({ success: true, inserted, total: byTs.size, windows: ["recent (60 min)", "last24h", "last48h"] });
  } catch (error) {
    res.status(502).json({
      success: false,
      error: error.name === "AbortError" || error.name === "TimeoutError"
        ? "OCS connection timed out"
        : error.message,
    });
  }
});

// ── Saved Targets routes ──────────────────────────────────────────────────────

app.get("/api/targets", (req, res) => {
  const targets = db.prepare("SELECT * FROM targets ORDER BY created_at DESC").all();
  res.json({ success: true, targets });
});

app.post("/api/targets", (req, res) => {
  const { name, ra, dec, ra_deg, dec_deg, type, source, notes } = req.body;
  if (!name) return res.status(400).json({ success: false, error: "name is required" });
  const stmt = db.prepare(
    "INSERT INTO targets (name, ra, dec, ra_deg, dec_deg, type, source, notes) VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
  );
  const result = stmt.run(
    name,
    ra   || null,
    dec  || null,
    ra_deg  != null ? Number(ra_deg)  : null,
    dec_deg != null ? Number(dec_deg) : null,
    type   || null,
    source || "manual",
    notes  || null
  );
  res.json({ success: true, id: result.lastInsertRowid });
});

app.delete("/api/targets/:id", (req, res) => {
  db.prepare("DELETE FROM targets WHERE id = ?").run(req.params.id);
  res.json({ success: true });
});

// ── Pipeline jobs routes ──────────────────────────────────────────────────────

app.get("/api/pipeline/jobs", (req, res) => {
  const jobs = db.prepare(
    "SELECT * FROM pipeline_jobs ORDER BY created_at DESC LIMIT 100"
  ).all();
  if (!jobs.length) return res.json({ success: true, jobs: [] });
  const ids = jobs.map(j => j.id);
  const placeholders = ids.map(() => "?").join(",");
  const results = db.prepare(
    `SELECT * FROM pipeline_results WHERE job_id IN (${placeholders}) ORDER BY created_at ASC`
  ).all(...ids);
  const resultsByJob = {};
  for (const r of results) {
    if (!resultsByJob[r.job_id]) resultsByJob[r.job_id] = [];
    resultsByJob[r.job_id].push(r);
  }
  for (const j of jobs) j.results = resultsByJob[j.id] || [];
  res.json({ success: true, jobs });
});

app.delete("/api/pipeline/results/:id", (req, res) => {
  db.prepare("DELETE FROM pipeline_results WHERE id = ?").run(req.params.id);
  res.json({ success: true });
});

app.post("/api/pipeline/trigger", async (req, res) => {
  const { fits_dir, target, filter, exposure, selected_files, target_filter } = req.body;
  if (!fits_dir) return res.status(400).json({ success: false, error: "fits_dir is required" });

  const result = db.prepare(
    "INSERT INTO pipeline_jobs (target, filter, exposure, fits_dir, status, target_filter) VALUES (?, ?, ?, ?, ?, ?)"
  ).run(target || null, filter || null, exposure || null, fits_dir, "queued", target_filter || null);

  const job_id = result.lastInsertRowid;

  callProcessingService(job_id, fits_dir, target || "Unknown", selected_files || null, target_filter || null)
    .then((svc) => {
      if (!svc.success) {
        db.prepare("UPDATE pipeline_jobs SET status='error', error=? WHERE id=?")
          .run(svc.error || "Service call failed", job_id);
      }
    })
    .catch((err) => {
      console.error("Processing service unreachable:", err.message);
      db.prepare("UPDATE pipeline_jobs SET status='error', error=? WHERE id=?")
        .run("Processing service unreachable — is it running?", job_id);
    });

  res.json({ success: true, id: job_id });
});

app.post("/api/pipeline/rerun/:id", async (req, res) => {
  const oldJob = db.prepare("SELECT * FROM pipeline_jobs WHERE id=?").get(req.params.id);
  if (!oldJob) return res.status(404).json({ success: false, error: "Job not found" });

  const result = db.prepare(
    "INSERT INTO pipeline_jobs (target, filter, exposure, fits_dir, status, target_filter) VALUES (?, ?, ?, ?, ?, ?)"
  ).run(oldJob.target, oldJob.filter, oldJob.exposure, oldJob.fits_dir, "queued", oldJob.target_filter || null);

  const job_id = result.lastInsertRowid;

  callProcessingService(job_id, oldJob.fits_dir, oldJob.target || "Unknown", null, oldJob.target_filter || null)
    .then((svc) => {
      if (!svc.success)
        db.prepare("UPDATE pipeline_jobs SET status='error', error=? WHERE id=?")
          .run(svc.error || "Service call failed", job_id);
    })
    .catch((err) => {
      db.prepare("UPDATE pipeline_jobs SET status='error', error=? WHERE id=?")
        .run("Processing service unreachable — is it running?", job_id);
    });

  res.json({ success: true, id: job_id });
});

// ── Retry a single pipeline step ─────────────────────────────────────────────
// POST /api/pipeline/results/:id/retry-step   body: { step: "calibrate"|"solve"|"upload"|"inspect"|"photometry"|"subtraction" }
//
// Local steps (calibrate/solve/upload) call the processing service /resume
// endpoint so only the failed step and everything after it is re-run.
// STDWeb-only steps (inspect/photometry/subtraction) additionally offer the
// option to re-trigger directly via STDWeb API (the poller then catches up).
const PROC_URL = process.env.PROC_URL || "http://127.0.0.1:5200";

app.post("/api/pipeline/results/:id/retry-step", async (req, res) => {
  const result = db.prepare("SELECT * FROM pipeline_results WHERE id=?").get(req.params.id);
  if (!result) return res.status(404).json({ success: false, error: "Result not found" });

  const step = req.body?.step;
  const VALID_STEPS = ["calibrate", "solve", "upload", "inspect", "photometry", "subtraction"];
  if (!VALID_STEPS.includes(step)) {
    return res.status(400).json({ success: false, error: `Unknown step: ${step}. Valid: ${VALID_STEPS.join(", ")}` });
  }

  try {
    // All steps — delegate to processing service /resume
    const r = await fetch(`${PROC_URL}/resume`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ result_id: result.id, from_step: step }),
    });
    const data = await r.json();
    if (!data.success) return res.status(400).json(data);

    // Optimistically clear error so UI shows running immediately
    db.prepare("UPDATE pipeline_results SET status='processing', error=NULL, updated_at=datetime('now') WHERE id=?")
      .run(result.id);

    res.json({ success: true, step, message: "Resuming from step" });
  } catch (e) {
    res.status(500).json({ success: false, error: `Processing service unreachable: ${e.message}` });
  }
});

// ── NAS browser ───────────────────────────────────────────────────────────────

function hasFits(dir) {
  try {
    return fs.readdirSync(dir).some((f) => /\.(fit|fits)$/i.test(f));
  } catch { return false; }
}

function listSubdirs(dir) {
  try {
    return fs.readdirSync(dir, { withFileTypes: true })
      .filter((d) => d.isDirectory())
      .map((d) => d.name)
      .sort();
  } catch { return []; }
}

/** Read FITS OBJECT keyword from first header block without astropy */
function readFitsObject(filePath) {
  try {
    const BLOCK = 2880;
    const fd = fs.openSync(filePath, "r");
    const buf = Buffer.alloc(BLOCK);
    let end = false;
    let object = null;
    while (!end) {
      const bytesRead = fs.readSync(fd, buf, 0, BLOCK, null);
      if (bytesRead === 0) break;
      for (let i = 0; i < bytesRead; i += 80) {
        const card = buf.slice(i, i + 80).toString("ascii");
        if (card.startsWith("END ") || card.trim() === "END") { end = true; break; }
        if (card.startsWith("OBJECT  =") || card.startsWith("OBJECT =")) {
          const val = card.slice(10).split("/")[0].trim().replace(/^'+|'+$/g, "").trim();
          if (val) object = val;
        }
      }
    }
    fs.closeSync(fd);
    return object;
  } catch { return null; }
}

/** Strip common transient prefixes for deduplication */
function canonicalTarget(name) {
  return name.replace(/^(AT|SN|TCP|CSS|MASTER|PSN|PNV)\s*/i, "").trim().toUpperCase();
}

/** Given a list of raw target names, merge those with the same canonical form.
 *  Prefers "AT" prefix over "SN"; falls back to the first seen. */
function deduplicateTargets(names) {
  const map = new Map(); // canonical → preferred raw name
  for (const n of names) {
    const key = canonicalTarget(n);
    if (!map.has(key)) { map.set(key, n); continue; }
    const prev = map.get(key);
    // Prefer AT over SN, otherwise keep first
    if (/^SN\s*/i.test(prev) && /^AT\s*/i.test(n)) map.set(key, n);
  }
  return [...map.values()];
}

/** For a SNAPSHOT folder, scan FITS files and return unique real target names */
function snapshotTargets(folderPath) {
  try {
    const SKIP_OBJECTS = new Set(["snapshot", "dark", "flat", "bias", "test_target",
                                   "test_target_v2", "test", "unknown", "none", ""]);
    const files = fs.readdirSync(folderPath).filter((f) => /\.(fit|fits)$/i.test(f));
    const seen = new Set();
    for (const f of files) {
      const obj = readFitsObject(path.join(folderPath, f));
      if (obj && !SKIP_OBJECTS.has(obj.toLowerCase())) seen.add(obj);
    }
    return deduplicateTargets([...seen]);
  } catch { return []; }
}

// GET /api/pipeline/nas-dates
// Step 1 — fast: returns only the list of date folders.
// Does NOT read any FITS headers (no SNAPSHOT scanning).
app.get("/api/pipeline/nas-dates", (req, res) => {
  try {
    const dates = fs.readdirSync(NAS_WATCH_PATH, { withFileTypes: true })
      .filter((d) => d.isDirectory() && /^\d{4}-\d{2}-\d{2}$/.test(d.name))
      .map((d) => d.name)
      .sort()
      .reverse();
    res.json({ success: true, dates, base: NAS_WATCH_PATH });
  } catch (err) {
    res.json({ success: true, dates: [], base: NAS_WATCH_PATH, error: err.message });
  }
});

// GET /api/pipeline/nas-dates/:date
// Step 2 — on demand: scan one night's folder and return its targets.
// Reads FITS headers only for that night's SNAPSHOT (to discover real target names).
app.get("/api/pipeline/nas-dates/:date", (req, res) => {
  const { date } = req.params;
  if (!/^\d{4}-\d{2}-\d{2}$/.test(date)) {
    return res.status(400).json({ success: false, error: "Invalid date format" });
  }
  try {
    const datePath = path.join(NAS_WATCH_PATH, date);
    if (!fs.existsSync(datePath)) {
      return res.json({ success: true, date, targets: [] });
    }

    const subs = listSubdirs(datePath);
    const targets = [];

    for (const sub of subs) {
      const subPath = path.join(datePath, sub);
      const directFits = hasFits(subPath);
      const deepFits = !directFits && listSubdirs(subPath).some((s) =>
        hasFits(path.join(subPath, s))
      );
      if (!directFits && !deepFits) continue;

      if (sub.toUpperCase() === "SNAPSHOT" && directFits) {
        const objects = snapshotTargets(subPath);
        if (objects.length > 0) {
          for (const obj of objects) {
            targets.push({ name: obj, path: subPath, hasFits: true, snapshot: true });
          }
        } else {
          targets.push({ name: sub, path: subPath, hasFits: true });
        }
      } else {
        targets.push({ name: sub, path: subPath, hasFits: true });
      }
    }

    res.json({ success: true, date, targets });
  } catch (err) {
    res.status(500).json({ success: false, error: err.message });
  }
});

app.get("/api/pipeline/job/:id/log", async (req, res) => {
  try {
    const url = `${PROCESSING_SERVICE_URL}/jobs/${req.params.id}/log`;
    http.get(url, (r) => {
      let data = "";
      r.on("data", (c) => (data += c));
      r.on("end", () => {
        try { res.json(JSON.parse(data)); } catch { res.json({ success: false, log: data }); }
      });
    }).on("error", (e) => res.json({ success: false, error: e.message }));
  } catch (err) {
    res.json({ success: false, error: err.message });
  }
});

app.delete("/api/pipeline/jobs/:id", (req, res) => {
  db.prepare("DELETE FROM pipeline_jobs WHERE id = ?").run(req.params.id);
  res.json({ success: true });
});

// ── Integration test: filter-patch end-to-end ─────────────────────────────────
// POST /api/test/filter-patch
// Exercises the real production path:
//   stepCaptureFilter (change filter → capture → waitForCameraIdle → patchFitsFilterHeaders)
// then reads back the FITS header to confirm it was corrected.
// Runs 1-second exposures for each of the two test filters (G and BP by default).
app.post("/api/test/filter-patch", async (req, res) => {
  const testFilters = (req.body?.filters || "G,BP,RP").split(",").map(s => s.trim()).filter(Boolean);
  const duration    = Number(req.body?.duration ?? 5);
  const gain        = Number(req.body?.gain ?? 10);
  const count       = Math.max(1, Math.min(10, Number(req.body?.count ?? 2)));

  let ninaConfig;
  try { ninaConfig = normalizeTargetConfig(req.body || {}); }
  catch (e) { return res.status(400).json({ error: e.message }); }

  const log = [];
  const results = [];

  const tlog = (...args) => {
    const line = args.join(" ");
    log.push(line);
    console.log("[filter-patch-test]", line);
  };

  try {
    // Get available filters from NINA
    const { equipmentInfo } = await getEquipmentInfo(ninaConfig);
    const filters = equipmentInfo?.FilterWheel?.AvailableFilters || [];
    const available = filters.map(f => f.Name).join(", ");
    tlog(`Available filters: ${available}`);
    if (!filters.length) return res.status(503).json({ error: "No filters available from NINA", log });

    // Temporarily enable sequence state so stepCaptureFilter can run
    const prevAborted  = seqState.aborted;
    const prevRunning  = seqState.running;
    seqState.aborted   = false;
    seqState.running   = true;

    try {
      for (const filterName of testFilters) {
        tlog(`\n--- Testing filter: ${filterName} ---`);

        const filter = findFilterByName(filters, filterName);
        if (!filter) {
          tlog(`  SKIP: '${filterName}' not found in filterwheel`);
          results.push({ filter: filterName, status: "SKIP", reason: "not in filterwheel" });
          continue;
        }

        // Snapshot existing files BEFORE capture (avoids NAS/Windows clock-skew issue)
        const beforeFiles = snapshotNasFiles();

        // Run the real stepCaptureFilter — this triggers capture + patch
        tlog(`  Calling stepCaptureFilter (${count} × ${duration}s, gain ${gain}) ...`);
        await stepCaptureFilter(ninaConfig, "FILTER_TEST", filterName, count, duration, gain, filters);

        // Find new FITS files by set-difference
        const newFiles = [];
        for (const entry of fs.readdirSync(NAS_WATCH_PATH, { withFileTypes: true })) {
          if (!entry.isDirectory()) continue;
          const snapDir = path.join(NAS_WATCH_PATH, entry.name, "SNAPSHOT");
          if (!fs.existsSync(snapDir)) continue;
          for (const name of fs.readdirSync(snapDir)) {
            if (!/\.(fit|fits)$/i.test(name)) continue;
            const fullPath = path.join(snapDir, name);
            if (beforeFiles.has(fullPath)) continue; // pre-existing — skip
            // Read the FITS FILTER header
            let header_filter = "?";
            try {
              const readScript = `from astropy.io import fits; h = fits.getheader(${JSON.stringify(fullPath)}, memmap=False); print(h.get('FILTER','?'))`;
              header_filter = execSync(
                `.venv/bin/python -c "${readScript.replace(/"/g, '\\"')}"`,
                { cwd: path.join(__dirname, ".."), timeout: 10000 }
              ).toString().trim();
            } catch (e) {
              header_filter = `ERROR: ${e.message.slice(0, 80)}`;
            }
            const pass = header_filter === filter.Name;
            tlog(`  File: ${name}  FILTER header='${header_filter}'  expected='${filter.Name}'  → ${pass ? "PASS ✓" : "FAIL ✗"}`);
            newFiles.push({ name, header_filter, expected: filter.Name, pass });
          }
        }

        if (newFiles.length === 0) {
          tlog(`  WARNING: No new files found — camera may not have saved`);
          results.push({ filter: filterName, status: "NO_FILE", files: [] });
        } else {
          const allPass = newFiles.every(f => f.pass);
          results.push({ filter: filterName, status: allPass ? "PASS" : "FAIL", files: newFiles });
        }
      }
    } finally {
      seqState.aborted  = prevAborted;
      seqState.running  = prevRunning;
    }

    const overall = results.every(r => r.status === "PASS" || r.status === "SKIP")
      ? "PASS" : "FAIL";
    tlog(`\n=== Overall result: ${overall} ===`);
    res.json({ overall, results, log });

  } catch (e) {
    tlog(`ERROR: ${e.message}`);
    res.status(500).json({ error: e.message, log });
  }
});

// ── STDWeb photometry proxy ───────────────────────────────────────────────────
// GET /api/stdweb/task/:task_id/photometry
// Returns photometry measurements for a task. Serves from DB when already
// stored; fetches from STDWeb live otherwise and persists the result.
app.get("/api/stdweb/task/:task_id/photometry", async (req, res) => {
  const { task_id } = req.params;
  const STDWEB_URL   = process.env.STDWEB_URL   || "http://86.253.141.183:7000";
  const STDWEB_TOKEN = process.env.STDWEB_TOKEN  || "1e296ddd6738af45467b7bc6558c00a9524447ab";
  const headers = { "Authorization": `Token ${STDWEB_TOKEN}` };

  try {
    // Serve from DB only when photometry AND subtraction data are both present.
    // If mag_sub and mag_sub_ul are both NULL the subtraction hasn't been
    // stored yet — fetch live so we get the full result and persist it.
    const stored = db.prepare(
      `SELECT r.id, r.target, r.filter, r.obs_date, r.mjd,
              r.mag_ap, r.magerr_ap, r.mag_sub, r.magerr_sub, r.mag_sub_ul, r.status
       FROM pipeline_results r
       WHERE r.stdweb_task_id = ?
         AND r.mjd IS NOT NULL
         AND (r.mag_sub IS NOT NULL OR r.mag_sub_ul IS NOT NULL)
       LIMIT 1`
    ).get(task_id);

    if (stored) {
      return res.json({
        task_id: parseInt(task_id),
        target: stored.target,
        filter: stored.filter,
        obs_date: stored.obs_date,
        mjd:    stored.mjd,
        direct: stored.mag_ap    != null ? { mag: stored.mag_ap,    magerr: stored.magerr_ap }    : null,
        sub:    stored.mag_sub   != null ? { mag: stored.mag_sub,   magerr: stored.magerr_sub }   : null,
        sub_ul: stored.mag_sub_ul != null ? { ul: stored.mag_sub_ul } : null,
        from_db: true,
      });
    }

    // Fetch live from STDWeb
    const taskJson = await fetch(`${STDWEB_URL}/api/tasks/${task_id}/`, { headers })
      .then(r => r.json()).catch(() => ({}));
    const filter = taskJson?.config?.filter || "?";

    const html = await fetch(`${STDWEB_URL}/tasks/${task_id}`, { headers })
      .then(r => r.text()).catch(() => "");

    const mjdMatch    = html.match(/MJD is ([\d.]+)/);
    const directMatch = html.match(/Primary target magnitude is \w+ = ([\d.]+) \+\/- ([\d.]+)/);
    const subMatch    = html.match(/Target magnitude is \w+ = ([\d.]+) \+\/- ([\d.]+)/);
    const subUlMatch  = !subMatch && html.match(/Target magnitude upper limit is \w+ > ([\d.]+)/);

    const mjd     = mjdMatch    ? parseFloat(mjdMatch[1])    : null;
    const direct  = directMatch ? { mag: parseFloat(directMatch[1]), magerr: parseFloat(directMatch[2]) } : null;
    const sub     = subMatch    ? { mag: parseFloat(subMatch[1]),    magerr: parseFloat(subMatch[2]) }    : null;
    const sub_ul  = subUlMatch  ? { ul: parseFloat(subUlMatch[1]) } : null;

    // Persist to DB if we got useful data
    if (mjd != null) {
      const obsDate = mjd ? (() => {
        const jd = mjd + 2400000.5;
        const d = new Date((jd - 2440587.5) * 86400000);
        return d.toISOString().slice(0, 10);
      })() : null;

      try {
        db.prepare(
          `UPDATE pipeline_results
           SET mjd=?, obs_date=COALESCE(obs_date,?),
               mag_ap=?, magerr_ap=?, mag_sub=?, magerr_sub=?, mag_sub_ul=?,
               updated_at=datetime('now')
           WHERE stdweb_task_id=?`
        ).run(mjd, obsDate,
              direct?.mag ?? null, direct?.magerr ?? null,
              sub?.mag    ?? null, sub?.magerr    ?? null,
              sub_ul?.ul  ?? null,
              task_id);
      } catch { /* non-fatal */ }
    }

    res.json({
      task_id: parseInt(task_id),
      target: taskJson?.config?.target || null,
      filter,
      mjd,
      direct,
      sub,
      sub_ul,
      from_db: false,
    });
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

// ── Pipeline history search ───────────────────────────────────────────────────
// GET /api/pipeline/history?target=&filter=&date_from=&date_to=&limit=
app.get("/api/pipeline/history", (req, res) => {
  try {
    const { target, filter, date_from, date_to, limit = 200 } = req.query;
    let sql = `
      SELECT r.id, r.job_id, r.target, r.filter, r.exposure, r.n_frames,
             r.obs_date, r.mjd, r.mag_ap, r.magerr_ap, r.mag_sub, r.magerr_sub,
             r.mag_sub_ul, r.stdweb_task_id, r.stdweb_url, r.status,
             j.fits_dir, j.created_at
      FROM pipeline_results r
      LEFT JOIN pipeline_jobs j ON j.id = r.job_id
      WHERE r.status = 'done'`;
    const params = [];
    if (target) { sql += ` AND lower(r.target) LIKE lower(?)`; params.push(`%${target}%`); }
    if (filter) { sql += ` AND r.filter = ?`; params.push(filter); }
    if (date_from) { sql += ` AND r.obs_date >= ?`; params.push(date_from); }
    if (date_to)   { sql += ` AND r.obs_date <= ?`; params.push(date_to); }
    sql += ` ORDER BY r.obs_date DESC, r.id DESC LIMIT ?`;
    params.push(parseInt(limit));

    const rows = db.prepare(sql).all(...params);
    res.json({ success: true, results: rows });
  } catch (err) {
    res.status(500).json({ success: false, error: err.message });
  }
});

app.use((req, res) => {
  res.sendFile(path.join(__dirname, "public", "index.html"));
});

app.listen(PORT, () => {
  console.log(`NINA control server listening on http://localhost:${PORT}`);
});
