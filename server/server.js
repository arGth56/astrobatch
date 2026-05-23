require("dotenv").config();
const express = require("express");
const path = require("path");
const https = require("https");
const http  = require("http");
const fs    = require("fs");
const { execSync, spawn } = require("child_process");
const Database = require("better-sqlite3");
const nodemailer = require("nodemailer");
const alerts = require("./alerts");

const TOO_MOSAIC_MAX_SIDE = Math.max(1, parseInt(process.env.TOO_MOSAIC_MAX_SIDE || "8", 10));
const PROCESSING_SERVICE_URL = process.env.PROCESSING_SERVICE_URL || "http://127.0.0.1:5200";
const NAS_WATCH_PATH = process.env.NAS_WATCH_PATH || "/mnt/nas/input/pyl/astro/input";
const PYTHON_BIN_CANDIDATE = process.env.PYTHON_BIN || path.join(__dirname, "..", ".venv", "bin", "python");
const PYTHON_BIN = fs.existsSync(PYTHON_BIN_CANDIDATE) ? PYTHON_BIN_CANDIDATE : "python3";

/** ToO-Fermi-801108913 → "Fermi"; non-ToO → null (TNS/manual). */
function parseTooBroker(name) {
  const m = String(name || "").trim().match(/^ToO-([^-]+)-/i);
  return m ? m[1] : null;
}

/** Alerts default to transient detection (no STDWeb target); TNS/SN always use target photometry. */
function stdwebUseTargetForObject(objectName) {
  const broker = parseTooBroker(objectName);
  if (!broker) return 1;
  const strat = getStrategy(broker);
  return strat.stdweb_use_target ? 1 : 0;
}

function callProcessingService(job_id, fits_dir, target, selected_files = null, target_filter = null, manual_selection = false, force_fresh = false, stdweb_use_target = null) {
  return new Promise((resolve, reject) => {
    const payload = { job_id, fits_dir, target };
    if (selected_files && selected_files.length) payload.selected_files = selected_files;
    if (target_filter) payload.target_filter = target_filter;
    if (manual_selection) payload.manual_selection = true;
    if (force_fresh) payload.force_fresh = true;
    if (stdweb_use_target !== null && stdweb_use_target !== undefined) {
      payload.stdweb_use_target = !!stdweb_use_target;
    }
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
const STDWEB_TOKEN = process.env.STDWEB_TOKEN || "";
// Runtime helper — checks DB first so the UI save takes effect without restart
const getStdwebToken = () => getSetting("secret_stdweb_token") || process.env.STDWEB_TOKEN || "";
const DEVICE_TYPES = ["mount", "camera", "filterwheel", "focuser", "flatdevice"];
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

function warnMissingSecrets() {
  const warnings = [];
  if (!STDWEB_TOKEN) {
    warnings.push("STDWEB_TOKEN is missing (STDWeb health/photometry proxy will be disabled).");
  }
  if (!TNS_API_KEY || !TNS_BOT_ID || !TNS_BOT_NAME) {
    warnings.push("TNS credentials are incomplete (TNS lookup requires TNS_BOT_ID, TNS_BOT_NAME, TNS_API_KEY).");
  }
  if (warnings.length) {
    console.warn("[startup] Missing optional secrets:");
    for (const w of warnings) console.warn(`[startup] - ${w}`);
  }
}
warnMissingSecrets();

// ── Email notifications ───────────────────────────────────────────────────────
function getEmailConfig() {
  return {
    host:     getSetting("email_host")     || process.env.EMAIL_HOST     || "",
    port:     parseInt(getSetting("email_port") || process.env.EMAIL_PORT || "587", 10),
    secure:   (getSetting("email_secure")  || process.env.EMAIL_SECURE   || "false") === "true",
    user:     getSetting("email_user")     || process.env.EMAIL_USER     || "",
    pass:     getSetting("secret_email_pass") || process.env.EMAIL_PASS  || "",
    from:     getSetting("email_from")     || process.env.EMAIL_FROM     || "",
    to:       getSetting("email_to")       || process.env.EMAIL_TO       || "",
  };
}

async function sendAlertEmail(subject, body) {
  const cfg = getEmailConfig();
  if (!cfg.host || !cfg.to || !cfg.user || !cfg.pass) {
    seqLog("[email] Not configured — skipping notification", "warn");
    return;
  }
  try {
    const transporter = nodemailer.createTransport({
      host: cfg.host, port: cfg.port, secure: cfg.secure,
      auth: { user: cfg.user, pass: cfg.pass },
    });
    await transporter.sendMail({
      from: cfg.from || cfg.user,
      to:   cfg.to,
      subject,
      text: body,
    });
    seqLog(`[email] Notification sent to ${cfg.to}`);
  } catch (e) {
    seqLog(`[email] Failed to send: ${e.message}`, "error");
  }
}

app.use(express.json());
app.use(express.static(path.join(__dirname, "public")));

// ── Processed data browser ────────────────────────────────────────────────────
const DATA_DIR  = path.join(__dirname, "..", "data");
const NAS_OUTPUT = process.env.NAS_OUTPUT || "/mnt/nas/input/pyl/astro/output";

// Serve files for download: GET /data/2026-03-19/SN%202026fvx/G/120s/res.fit
app.use("/data",       express.static(DATA_DIR,   { dotfiles: "deny" }));
// Serve NAS output stacks: GET /nas-output/2026-03-31/sn_2026fvx/G/res_preview.png
app.use("/nas-output", express.static(NAS_OUTPUT, { dotfiles: "deny" }));

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
    use_color       INTEGER DEFAULT 1,
    refine_wcs      INTEGER DEFAULT 1,
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
try { db.prepare("ALTER TABLE pipeline_jobs ADD COLUMN use_color INTEGER DEFAULT 0").run(); } catch { /* already exists */ }
try { db.prepare("ALTER TABLE pipeline_jobs ADD COLUMN refine_wcs INTEGER DEFAULT 1").run(); } catch { /* already exists */ }
try { db.prepare("ALTER TABLE pipeline_jobs ADD COLUMN stdweb_use_target INTEGER").run(); } catch { /* already exists */ }
try { db.prepare("ALTER TABLE pipeline_results ADD COLUMN stdweb_state TEXT").run(); } catch { /* already exists */ }
try { db.prepare("ALTER TABLE pipeline_results ADD COLUMN obs_date TEXT").run(); } catch { /* already exists */ }
// seq_queue source column: "manual" (user-added) or "alert" (pushed from alert broker)
try { db.prepare("ALTER TABLE seq_queue ADD COLUMN source TEXT DEFAULT 'manual'").run(); } catch { /* already exists */ }
// seq_queue cycles column
try { db.prepare("ALTER TABLE seq_queue ADD COLUMN cycles INTEGER DEFAULT 1").run(); } catch { /* already exists */ }
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
    sky       REAL,
    ir_sky    REAL
  );
`);
// Migration: add ir_sky column if DB was created before this version
try { db.exec(`ALTER TABLE ocs_history ADD COLUMN ir_sky REAL`); } catch (_) {}

// ── Rain forecast log table ──────────────────────────────────────────────────
db.exec(`
  CREATE TABLE IF NOT EXISTS rain_forecast_log (
    id           INTEGER PRIMARY KEY AUTOINCREMENT,
    ts           TEXT NOT NULL DEFAULT (datetime('now')),
    warn_level   TEXT,
    max_prob_30  REAL,
    max_prob_60  REAL,
    slots_json   TEXT
  );
  CREATE INDEX IF NOT EXISTS idx_rfl_ts ON rain_forecast_log(ts DESC);
`);

// ── Autofocus log table ──────────────────────────────────────────────────────
db.exec(`
  CREATE TABLE IF NOT EXISTS autofocus_log (
    id        INTEGER PRIMARY KEY AUTOINCREMENT,
    filter    TEXT    NOT NULL,
    position  INTEGER NOT NULL,
    ts        TEXT    NOT NULL DEFAULT (datetime('now'))
  );
  CREATE INDEX IF NOT EXISTS idx_af_filter_ts ON autofocus_log(filter, ts DESC);
`);

// ── Key-value settings table (persistent user preferences) ───────────────────
db.exec(`
  CREATE TABLE IF NOT EXISTS settings (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
  );
`);
const getSetting = (key, fallback = null) => {
  const row = db.prepare("SELECT value FROM settings WHERE key=?").get(key);
  return row ? row.value : fallback;
};
const setSetting = (key, value) => {
  db.prepare("INSERT OR REPLACE INTO settings (key, value) VALUES (?,?)").run(key, String(value));
};

// ── Per-broker alert strategy table ──────────────────────────────────────────
db.exec(`
  CREATE TABLE IF NOT EXISTS alert_strategies (
    broker       TEXT PRIMARY KEY,
    enabled      INTEGER NOT NULL DEFAULT 1,
    mode         TEXT    NOT NULL DEFAULT 'too',
    exposure     REAL    NOT NULL DEFAULT 30,
    gain         INTEGER NOT NULL DEFAULT 10,
    filter_cycle TEXT    NOT NULL DEFAULT 'G,RP,BP',
    rapid_count  INTEGER NOT NULL DEFAULT 15,
    do_af        INTEGER NOT NULL DEFAULT 0,
    do_guiding   INTEGER NOT NULL DEFAULT 0,
    do_center    INTEGER NOT NULL DEFAULT 0,
    min_alt      REAL    NOT NULL DEFAULT 25,
    max_err_deg  REAL    NOT NULL DEFAULT 1.0,
    notes        TEXT    DEFAULT '',
    notify_email INTEGER NOT NULL DEFAULT 1
  );
`);

const ALERT_STRATEGY_SEEDS = [
  { broker: "Swift",    mode: "too",    exposure: 30,  filter_cycle: "G,RP,BP", rapid_count: 15, max_err_deg: 1.0,  notes: "Fast BAT localization — rapid color cycling" },
  { broker: "Fermi",    mode: "ignore", exposure: 30,  filter_cycle: "G,RP,BP", rapid_count: 15, max_err_deg: 0.1,  notes: "Error box too large; only LAT-refined useful" },
  { broker: "EP",       mode: "too",    exposure: 120, filter_cycle: "G,BP,RP", rapid_count: 10, max_err_deg: 1.0,  notes: "WXT ~4-8 arcmin — fits FOV" },
  { broker: "IceCube",  mode: "queue",  exposure: 120, filter_cycle: "G,BP,RP", rapid_count: 10, max_err_deg: 1.0,  notes: "Queue for manual review" },
  { broker: "LVK",      mode: "ignore", exposure: 120, filter_cycle: "G,BP,RP", rapid_count: 10, max_err_deg: 1.0,  notes: "Cannot tile GW skymaps" },
  { broker: "INTEGRAL", mode: "ignore", exposure: 30,  filter_cycle: "G,RP,BP", rapid_count: 15, max_err_deg: 1.0,  notes: "Science ops ended Feb 2025 — re-entry 2029" },
  { broker: "AGILE",    mode: "ignore", exposure: 30,  filter_cycle: "G,RP,BP", rapid_count: 15, max_err_deg: 1.0,  notes: "Decommissioned Feb 2024 — re-entered atmosphere" },
  { broker: "GECAM",    mode: "too",    exposure: 30,  filter_cycle: "G,RP,BP", rapid_count: 15, max_err_deg: 1.0,  notes: "Good localization — rapid response" },
  { broker: "MAXI",     mode: "queue",  exposure: 120, filter_cycle: "G,BP,RP", rapid_count: 10, max_err_deg: 1.0,  notes: "Queue for review" },
  { broker: "IPN",      mode: "queue",  exposure: 120, filter_cycle: "G,BP,RP", rapid_count: 10, max_err_deg: 1.0,  notes: "Queue for review" },
  { broker: "SN-nu",    mode: "queue",  exposure: 120, filter_cycle: "G,BP,RP", rapid_count: 10, max_err_deg: 1.0,  notes: "Queue for review" },
  { broker: "SVOM",     mode: "too",    exposure: 30,  filter_cycle: "G,RP,BP", rapid_count: 15, max_err_deg: 1.0,  notes: "ECLAIRs/MXT GRB — rapid response like Swift" },
  { broker: "GCN",      mode: "queue",  exposure: 120, filter_cycle: "G,BP,RP", rapid_count: 10, max_err_deg: 1.0,  notes: "Fallback — queue for review" },
];

const _stratSeedStmt = db.prepare(`
  INSERT OR IGNORE INTO alert_strategies
    (broker, mode, exposure, filter_cycle, rapid_count, max_err_deg, notes)
  VALUES (@broker, @mode, @exposure, @filter_cycle, @rapid_count, @max_err_deg, @notes)
`);
for (const s of ALERT_STRATEGY_SEEDS) _stratSeedStmt.run(s);

// Migration: add notify_email column if it doesn't exist yet
try {
  db.prepare("ALTER TABLE alert_strategies ADD COLUMN notify_email INTEGER NOT NULL DEFAULT 1").run();
} catch { /* column already exists */ }
// STDWeb: 0 = transient detection (default for alerts); 1 = fill target field (TNS-style)
try {
  db.prepare("ALTER TABLE alert_strategies ADD COLUMN stdweb_use_target INTEGER NOT NULL DEFAULT 0").run();
} catch { /* column already exists */ }

const STRATEGY_DEFAULTS = {
  enabled: 1, mode: "too", exposure: 30, gain: 10, filter_cycle: "G,RP,BP",
  rapid_count: 15, do_af: 0, do_guiding: 0, do_center: 0, min_alt: 25, max_err_deg: 1.0, notes: "",
  stdweb_use_target: 0,
};

function getStrategy(broker) {
  const row = db.prepare("SELECT * FROM alert_strategies WHERE broker = ?").get(broker);
  return row || { broker, ...STRATEGY_DEFAULTS };
}

function getAllStrategies() {
  return db.prepare("SELECT * FROM alert_strategies ORDER BY broker").all();
}

// ── OCS history poller ────────────────────────────────────────────────────────

/** Strip HTML tags and trim whitespace from an OCS field value. */
function ocsClean(v) {
  return (v || "").replace(/<[^>]*>/g, "").replace(/&[a-z]+;/gi, "").trim();
}

/** Parse a numeric value from an OCS field string (handles " 25.7 °C", "nan", "--"). */
function ocsNum(v) {
  if (!v) return null;
  const n = parseFloat(ocsClean(v));
  return isNaN(n) ? null : n;
}

async function pollAndStoreOcs() {
  try {
    const host = DEFAULT_OCS_HOST;
    const { ok, fields } = await ocsStatus(host);
    if (!ok) return;
    db.prepare(
      `INSERT INTO ocs_history (roof, safe, rain, temp, humidity, pressure, sky, ir_sky)
       VALUES (?, ?, ?, ?, ?, ?, ?, ?)`
    ).run(
      ocsClean(fields.roof_sta)  || null,
      ocsClean(fields.stat_safe) || null,
      ocsClean(fields.wea_rain)  || null,
      ocsNum(fields.wea_temp),
      ocsNum(fields.wea_humd),
      ocsNum(fields.wea_pres),
      ocsNum(fields.wea_sq),
      ocsNum(fields.wea_irsky),
    );
  } catch (e) {
    if (process.env.DEBUG) console.error("[OCS poller]", e.message);
  }
}

setInterval(pollAndStoreOcs, 10 * 60 * 1000);  // store to DB every 10 min
setTimeout(pollAndStoreOcs, 8000);              // initial reading after server startup

// ── Rain forecast periodic logger ────────────────────────────────────────────
const storeRainForecast = db.prepare(
  `INSERT INTO rain_forecast_log (warn_level, max_prob_30, max_prob_60, slots_json)
   VALUES (?, ?, ?, ?)`
);
async function pollAndStoreRainForecast() {
  try {
    const fc = await fetchRainForecast();
    storeRainForecast.run(fc.warnLevel, fc.maxProb30, fc.maxProb60, JSON.stringify(fc.slots));
    if (process.env.DEBUG) console.log(`[rain-forecast] stored: ${fc.warnLevel} (${fc.maxProb30}%/30m)`);
  } catch (e) {
    if (process.env.DEBUG) console.error("[rain-forecast]", e.message);
  }
}
setInterval(pollAndStoreRainForecast, 15 * 60 * 1000); // every 15 min (matches forecast resolution)
setTimeout(pollAndStoreRainForecast, 12000);            // initial store after startup

// ── Queue persistence helpers ─────────────────────────────────────────────────

const _saveQueue = db.transaction(() => {
  db.prepare("DELETE FROM seq_queue").run();
  const ins = db.prepare(
    "INSERT INTO seq_queue (position, name, ra, dec, ra_deg, dec_deg, done, added_at, source, cycles) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
  );
  seqState.queue.forEach((t, i) => {
    ins.run(i, t.name, t.ra || "", t.dec || "", t.raDeg, t.decDeg, t.done ? 1 : 0, t.addedAt || Date.now(), t.source || "manual", t.cycles ?? 1);
  });
});

function saveQueue() {
  try { _saveQueue(); } catch (e) { console.error("saveQueue error:", e.message); }
}

function loadQueue() {
  try {
    const rows = db.prepare("SELECT * FROM seq_queue ORDER BY position").all();
    const seen = new Set();
    seqState.queue = rows
      .filter(r => {
        if (seen.has(r.name)) return false;  // deduplicate by name
        seen.add(r.name);
        return true;
      })
      .map(r => ({
        name:     r.name,
        ra:       r.ra,
        dec:      r.dec,
        raDeg:    r.ra_deg,
        decDeg:   r.dec_deg,
        done:     r.done === 1,
        addedAt:  r.added_at,
        source:   r.source || "manual",
        cycles:   r.cycles ?? 1,
        filters:  r.filters ? JSON.parse(r.filters) : null,
        count:    r.count    ?? null,
        duration: r.duration ?? null,
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
    colibriUid: getSetting("colibri_uid", "") || process.env.COLIBRI_UID || "",
  });
});

// GET  /api/settings/:key
// POST /api/settings  { key, value }
const SETTING_DEFAULTS = {
  minAlt:               "20",
  zenithLimit:          "70",
  meridianGap:          "10",
  minStars:             "10",
  maxCloudRecoveryWaitMin: "12",
  frameCheckEnabled:    "false",
  frameCheckThreshold:  "5",
  morning_park_hour:    "8",
  daily_reset_hour:     "12",
  alert_notify_email:   "true",
  // Inscribed-circle FOV (°) on sensor — used for ToO mosaic tile spacing
  telescope_fov_deg:    "0.86",
  // Focuser temperature model: pos ≈ slope × (FOCTEMP − 15) + ref15
  af_slope_G:           "10.83",
  af_slope_BP:          "9.73",
  af_slope_RP:          "12.08",
  af_ref15_G:           "4195",
  af_ref15_BP:          "4202",
  af_ref15_RP:          "4202",
};

app.get("/api/settings/:key", (req, res) => {
  const key = req.params.key;
  const val = getSetting(key) ?? SETTING_DEFAULTS[key] ?? null;
  if (val === null) return res.status(404).json({ success: false, error: "not found" });
  res.json({ success: true, key, value: val });
});

app.post("/api/settings", (req, res) => {
  const { key, value } = req.body || {};
  if (!key) return res.status(400).json({ success: false, error: "key required" });
  if (key === "telescope_fov_deg") {
    const n = parseFloat(value);
    if (!Number.isFinite(n) || n <= 0.05 || n > 30) {
      return res.status(400).json({
        success: false,
        error: "telescope_fov_deg must be a number between 0.05 and 30 (degrees)",
      });
    }
    setSetting(key, String(n));
    return res.json({ success: true });
  }
  setSetting(key, value ?? "");
  res.json({ success: true });
});

/** User-configured inscribed-circle FOV (°); env TOO_FOV_DEG overrides only if DB unset. */
function getTooFovDeg() {
  const fromDb = parseFloat(getSetting("telescope_fov_deg", ""));
  if (Number.isFinite(fromDb) && fromDb > 0.05 && fromDb <= 30) return fromDb;
  const fromEnv = parseFloat(process.env.TOO_FOV_DEG || "");
  if (Number.isFinite(fromEnv) && fromEnv > 0.05) return fromEnv;
  return parseFloat(SETTING_DEFAULTS.telescope_fov_deg || "0.86");
}

/** Single pointing when err ≤ threshold; default threshold = FOV/2 (2×err fits in one field). */
function getTooMosaicThresholdDeg(fovDeg = getTooFovDeg()) {
  if (process.env.TOO_MOSAIC_THRESHOLD_DEG != null && process.env.TOO_MOSAIC_THRESHOLD_DEG !== "") {
    const env = parseFloat(process.env.TOO_MOSAIC_THRESHOLD_DEG);
    if (Number.isFinite(env) && env > 0) return env;
  }
  return fovDeg / 2;
}

// ── Secrets ──────────────────────────────────────────────────────────────────
// GET  /api/secrets/:key  — returns { set: bool } only, never the value
// POST /api/secrets       — { key, value } stores in DB (survives restarts, no file edit needed)
const ALLOWED_SECRET_KEYS = new Set(["secret_stdweb_token", "secret_email_pass"]);

app.get("/api/secrets/:key", (req, res) => {
  const key = req.params.key;
  if (!ALLOWED_SECRET_KEYS.has(key))
    return res.status(400).json({ success: false, error: "unknown secret key" });
  const val = getSetting(key) || "";
  res.json({ success: true, key, set: val.length > 0 });
});

app.post("/api/secrets", (req, res) => {
  const { key, value } = req.body || {};
  if (!key) return res.status(400).json({ success: false, error: "key required" });
  if (!ALLOWED_SECRET_KEYS.has(key))
    return res.status(400).json({ success: false, error: "unknown secret key" });
  setSetting(key, value ?? "");
  res.json({ success: true });
});

// POST /api/email/test — send a test email with current config
app.post("/api/email/test", async (req, res) => {
  const cfg = getEmailConfig();
  if (!cfg.host || !cfg.to || !cfg.user || !cfg.pass)
    return res.status(400).json({ success: false, error: "Email not fully configured (host, user, password, to required)" });
  try {
    const transporter = nodemailer.createTransport({
      host: cfg.host, port: cfg.port, secure: cfg.secure,
      auth: { user: cfg.user, pass: cfg.pass },
    });
    await transporter.sendMail({
      from: cfg.from || cfg.user,
      to:   cfg.to,
      subject: "🔭 AstroBatch — test notification",
      text:  `This is a test email from AstroBatch.\n\nSent at: ${new Date().toISOString()}\n`,
    });
    res.json({ success: true, to: cfg.to });
  } catch (e) {
    res.status(500).json({ success: false, error: e.message });
  }
});

function runQualityReplay(date, minStars, trailElong, trailSizeFactor, maxEvents) {
  const script = `
import json, sys, math
from pathlib import Path
import numpy as np
from astropy.io import fits
import sep

date = sys.argv[1]
snap = Path(sys.argv[2])
min_stars = int(sys.argv[3])
trail_elong = float(sys.argv[4])
trail_size_factor = float(sys.argv[5])
max_events = int(sys.argv[6])

if not snap.exists():
    print(json.dumps({"success": False, "error": f"Snapshot folder not found: {snap}"}))
    sys.exit(0)

rows = []
for fp in sorted(snap.glob("*.fit*")):
    try:
        data = fits.getdata(fp, 0)
        h = fits.getheader(fp, 0)
    except Exception:
        continue
    arr = np.asarray(data, dtype=np.float32)
    if arr.ndim > 2:
        arr = arr[0]
    arr = np.ascontiguousarray(np.nan_to_num(arr, nan=np.nanmedian(arr)))
    try:
        b = sep.Background(arr)
        objs = sep.extract(arr - b, 4.5 * b.globalrms, minarea=5)
        if len(objs):
            a = np.array(objs["a"])
            bb = np.array(objs["b"])
            flux = np.array(objs["flux"])
            good = (a > 0.8) & (bb > 0.6) & (flux > 0)
            a = a[good]
            bb = bb[good]
            stars = int(len(a))
            elong = float(np.median(a / np.maximum(bb, 1e-3))) if stars else None
            size = float(np.median(np.sqrt(a * bb))) if stars else None
        else:
            stars, elong, size = 0, None, None
    except Exception:
        stars, elong, size = 0, None, None
    dt = h.get("DATE-OBS") or h.get("DATEOBS") or h.get("DATE") or ""
    obj = str(h.get("OBJECT", ""))
    rows.append((dt, fp.name, obj, stars, elong, size))

rows.sort(key=lambda r: r[0])
recent_sizes = []
in_cloud = False
events = []

for dt, name, obj, stars, elong, size in rows:
    if stars < min_stars and not in_cloud:
        events.append({
            "ts": dt,
            "type": "CLOUD_RECOVERY_START",
            "image": name,
            "target": obj,
            "reason": f"stars={stars}<min{min_stars}"
        })
        in_cloud = True
    elif stars >= min_stars and in_cloud:
        in_cloud = False

    if stars >= min_stars:
        med_size = float(np.median(recent_sizes)) if recent_sizes else None
        size_trail = (
            med_size is not None and size is not None and
            size > med_size * trail_size_factor
        )
        elong_trail = (elong is not None and elong >= trail_elong)
        if elong_trail or size_trail:
            bits = []
            if elong_trail:
                bits.append(f"elong={elong:.2f}")
            if size_trail:
                bits.append(f"size_jump={size:.2f}>{med_size * trail_size_factor:.2f}")
            events.append({
                "ts": dt,
                "type": "GUIDING_RECOVERY_START",
                "image": name,
                "target": obj,
                "reason": ", ".join(bits) if bits else "shape anomaly"
            })
        elif size is not None and not math.isnan(size) and size > 0:
            recent_sizes.append(size)
            if len(recent_sizes) > 8:
                recent_sizes.pop(0)

if max_events > 0:
    events = events[:max_events]

print(json.dumps({
    "success": True,
    "date": date,
    "snapshotPath": str(snap),
    "frameCount": len(rows),
    "activationCount": len(events),
    "activations": events
}))
`.trim();

  return new Promise((resolve, reject) => {
    const args = [
      "-c", script, date,
      path.join(NAS_WATCH_PATH, date, "SNAPSHOT"),
      String(minStars), String(trailElong), String(trailSizeFactor), String(maxEvents),
    ];
    const proc = spawn(PYTHON_BIN, args, { stdio: ["ignore", "pipe", "pipe"] });
    let out = "";
    let err = "";
    proc.stdout.on("data", (d) => { out += d.toString(); });
    proc.stderr.on("data", (d) => { err += d.toString(); });
    proc.on("error", reject);
    proc.on("close", (code) => {
      if (code !== 0) {
        return reject(new Error(err.trim() || `Replay exited with code ${code}`));
      }
      try {
        resolve(JSON.parse(out));
      } catch (e) {
        reject(new Error(`Invalid replay output: ${e.message}`));
      }
    });
  });
}

app.get("/api/quality/replay", async (req, res) => {
  try {
    const date = String(req.query.date || "").trim();
    if (!/^\d{4}-\d{2}-\d{2}$/.test(date)) {
      return res.status(400).json({ success: false, error: "date must be YYYY-MM-DD", received: date });
    }
    const minStarsRaw = Number(req.query.minStars ?? getSetting("minStars", SETTING_DEFAULTS.minStars));
    const minStars = Number.isFinite(minStarsRaw) ? Math.max(0, Math.round(minStarsRaw)) : 10;
    const maxEventsRaw = Number(req.query.maxEvents ?? 200);
    const maxEvents = Number.isFinite(maxEventsRaw) ? Math.max(1, Math.min(1000, Math.round(maxEventsRaw))) : 200;
    const result = await runQualityReplay(date, minStars, 1.85, 2.0, maxEvents);
    if (!result?.success) return res.status(404).json(result);
    res.json(result);
  } catch (e) {
    res.status(500).json({ success: false, error: e.message });
  }
});

// ── Alerts API routes (registered early for Express 5 compatibility) ─────────

app.get("/api/alerts/config", (req, res) => {
  res.json({
    success: true,
    alertReady: seqState.alertReady,
    tooRunning: seqState.tooRunning,
    tooFovDeg: getTooFovDeg(),
    tooMosaicThresholdDeg: getTooMosaicThresholdDeg(),
    tooMosaicMaxSide: TOO_MOSAIC_MAX_SIDE,
    alertPrepRunning: seqState.alertPrepRunning,
    alertPrepStep: seqState.alertPrepStep,
  });
});

app.post("/api/alerts/config", (req, res) => {
  if (req.body?.alertReady !== undefined) {
    seqState.alertReady = Boolean(req.body.alertReady);
    setSetting("alertReady", seqState.alertReady ? "true" : "false");
    seqLog(`Alert-Ready mode ${seqState.alertReady ? "ENABLED" : "DISABLED"}`, seqState.alertReady ? "warn" : "info");
  }
  res.json({ success: true, alertReady: seqState.alertReady });
});

app.post("/api/alerts/prep", async (req, res) => {
  if (seqState.alertPrepRunning) {
    return res.status(409).json({ success: false, error: "alert prep already running" });
  }
  if (seqState.running) {
    return res.status(409).json({ success: false, error: "a sequence is already running — abort it first" });
  }
  const ninaConfig = normalizeTargetConfig(req.body);
  const seqConfig = {
    filters: req.body?.filters || ["L", "G", "BP", "RP"],
    gain:    parseInt(req.body?.gain) || 30,
  };
  res.json({ success: true, message: "alert prep started" });
  runAlertPrep(ninaConfig, seqConfig).catch(e => {
    console.error("[alert-prep] uncaught:", e);
  });
});

app.post("/api/alerts/prep/abort", (req, res) => {
  if (!seqState.alertPrepRunning) {
    return res.json({ success: true, message: "not running" });
  }
  seqState.aborted = true;
  res.json({ success: true, message: "abort requested" });
});

// ── Inject a simulated alert (dev/test) ──────────────────────────────────────
// POST /api/alerts/inject  { broker, trigger_id, ra, dec, err_deg, classification }
// Runs through the exact same strategy/altitude/moon checks as a real GCN alert.
app.post("/api/alerts/inject", async (req, res) => {
  const { broker = "SVOM", trigger_id = "TEST-" + Date.now(),
          ra, dec, err_deg = 0.1, classification = "GRB (simulated)" } = req.body || {};

  if (!Number.isFinite(Number(ra)) || !Number.isFinite(Number(dec))) {
    return res.status(400).json({ success: false, error: "ra and dec are required" });
  }
  const raDeg  = Number(ra);
  const decDeg = Number(dec);
  const errDeg = Number(err_deg);

  const strategy = getStrategy(broker);
  if (!strategy) return res.status(400).json({ success: false, error: `No strategy for broker: ${broker}` });

  const site = resolveObserverSite();
  const { alt } = serverComputeAltAz(raDeg, decDeg, site.lat, site.lon, new Date());
  const moon  = alerts.moonRaDec(new Date());
  const moonSep = alerts.angularSep(raDeg, decDeg, moon.ra, moon.dec);

  const minAltDeg    = strategy.min_alt      ?? 25;
  const maxErrDeg    = strategy.max_err_deg  ?? 1;
  const minMoonSep   = 20;

  let action = "rejected", action_reason = null;
  const rejectReasons = [];
  if (!Number.isFinite(alt)) rejectReasons.push("no-coords");
  else {
    if (alt < minAltDeg)     rejectReasons.push(`alt-too-low:${alt.toFixed(1)}deg`);
    if (moonSep < minMoonSep) rejectReasons.push(`moon-too-close:${moonSep.toFixed(1)}deg`);
    if (errDeg > maxErrDeg)  rejectReasons.push(`err-too-large:${errDeg.toFixed(2)}deg(limit=${maxErrDeg})`);
  }

  const ok = rejectReasons.length === 0;
  if (ok) {
    try {
      const alertObj = { broker, trigger_id, ra: raDeg, dec: decDeg,
                         err_deg: errDeg, alt, moonSep, strategy };
      const mode = strategy.mode || "too";
      _gcnPushToQueue({ name: `${broker}-${trigger_id}`, raDeg, decDeg,
                        alert: alertObj, mode });
      saveQueue();
      action = mode === "too" ? "queued" : "queued-only";
      action_reason = `strategy:${mode} (injected)`;
    } catch (e) {
      action = "rejected"; action_reason = `enqueue-error:${e.message}`;
    }
  } else {
    action_reason = rejectReasons.join(", ");
  }

  db.prepare(`
    INSERT INTO alerts (received_at, event_time, broker, topic, trigger_id, classification,
      ra, dec, err_deg, alt_now, moon_sep, action, action_reason, colibri_id, raw)
    VALUES (datetime('now'), datetime('now'), ?, 'injected', ?, ?,
      ?, ?, ?, ?, ?, ?, ?, null, '(injected via API)')
  `).run(broker, String(trigger_id), classification,
         raDeg, decDeg, errDeg,
         Number.isFinite(alt) ? alt : null,
         Number.isFinite(moonSep) ? moonSep : null,
         action, action_reason);

  seqLog(`🧪 Injected alert: ${broker} ${trigger_id} RA=${raDeg.toFixed(3)}° Dec=${decDeg.toFixed(3)}° alt=${alt?.toFixed(1)}° → ${action}${action_reason ? " ("+action_reason+")" : ""}`, ok ? "warn" : "info");

  res.json({ success: true, action, action_reason,
             alt: alt?.toFixed(1), moonSep: moonSep?.toFixed(1),
             strategy: { mode: strategy.mode, enabled: strategy.enabled } });
});

// ── Alert strategy routes ────────────────────────────────────────────────────

app.get("/api/alert-strategies", (req, res) => {
  res.json({ success: true, strategies: getAllStrategies() });
});

app.put("/api/alert-strategies/:broker", (req, res) => {
  const broker = req.params.broker;
  const existing = getStrategy(broker);
  if (!existing) return res.status(404).json({ success: false, error: "unknown broker" });

  const fields = {};
  const allowed = {
    enabled:      v => { const n = parseInt(v); return (n === 0 || n === 1) ? n : undefined; },
    mode:         v => ["too", "queue", "ignore"].includes(v) ? v : undefined,
    exposure:     v => { const n = parseFloat(v); return (n > 0 && n <= 600) ? n : undefined; },
    gain:         v => { const n = parseInt(v);   return (n >= 0 && n <= 100) ? n : undefined; },
    filter_cycle: v => (typeof v === "string" && v.trim()) ? v.trim() : undefined,
    rapid_count:  v => { const n = parseInt(v);   return (n >= 1 && n <= 200) ? n : undefined; },
    do_af:        v => { const n = parseInt(v); return (n === 0 || n === 1) ? n : undefined; },
    do_guiding:   v => { const n = parseInt(v); return (n === 0 || n === 1) ? n : undefined; },
    do_center:    v => { const n = parseInt(v); return (n === 0 || n === 1) ? n : undefined; },
    min_alt:      v => { const n = parseFloat(v); return (n >= 0 && n <= 90) ? n : undefined; },
    max_err_deg:  v => { const n = parseFloat(v); return (n > 0 && n <= 180) ? n : undefined; },
    notes:        v => (typeof v === "string") ? v.slice(0, 500) : undefined,
    stdweb_use_target: v => { const n = parseInt(v); return (n === 0 || n === 1) ? n : undefined; },
  };

  for (const [key, validate] of Object.entries(allowed)) {
    if (req.body?.[key] !== undefined) {
      const clean = validate(req.body[key]);
      if (clean !== undefined) fields[key] = clean;
    }
  }

  if (Object.keys(fields).length === 0) {
    return res.status(400).json({ success: false, error: "no valid fields to update" });
  }

  const sets = Object.keys(fields).map(k => `${k} = @${k}`).join(", ");
  db.prepare(`UPDATE alert_strategies SET ${sets} WHERE broker = @broker`).run({ ...fields, broker });
  res.json({ success: true, strategy: getStrategy(broker) });
});

app.post("/api/alert-strategies/:broker/reset", (req, res) => {
  const broker = req.params.broker;
  const seed = ALERT_STRATEGY_SEEDS.find(s => s.broker === broker);
  if (!seed) return res.status(404).json({ success: false, error: "unknown broker" });
  db.prepare("DELETE FROM alert_strategies WHERE broker = ?").run(broker);
  _stratSeedStmt.run(seed);
  res.json({ success: true, strategy: getStrategy(broker) });
});

// ── Alerts list/detail routes ────────────────────────────────────────────────

app.get("/api/alerts", (req, res) => {
  const limit  = Math.min(500, Math.max(1, parseInt(req.query.limit)  || 100));
  const offset = Math.max(0, parseInt(req.query.offset) || 0);
  const rows   = alerts.listAlerts(db, { limit, offset });
  res.json({ success: true, alerts: rows });
});

app.get("/api/alerts/:id", (req, res) => {
  const row = alerts.getAlert(db, parseInt(req.params.id));
  if (!row) return res.status(404).json({ success: false, error: "not found" });
  res.json({ success: true, alert: row });
});

app.post("/api/alerts/:id/queue", (req, res) => {
  const row = alerts.getAlert(db, parseInt(req.params.id));
  if (!row) return res.status(404).json({ success: false, error: "not found" });
  if (!Number.isFinite(row.ra) || !Number.isFinite(row.dec)) {
    return res.status(400).json({ success: false, error: "alert has no coordinates" });
  }
  seqState.queue.push({
    name:   `${row.broker}-${row.trigger_id || row.id}`,
    raDeg:  row.ra, decDeg: row.dec, ra: "", dec: "",
    done: false, addedAt: Date.now(),
  });
  saveQueue();
  db.prepare("UPDATE alerts SET action = 'manual', action_reason = 'queued-by-operator' WHERE id = ?").run(row.id);
  res.json({ success: true });
});

app.post("/api/alerts/:id/ignore", (req, res) => {
  db.prepare("UPDATE alerts SET action = 'ignored', action_reason = ? WHERE id = ?")
    .run(String(req.body?.reason || "operator"), parseInt(req.params.id));
  res.json({ success: true });
});

// ── Sky Watchdog ─────────────────────────────────────────────────────────────
// Config persisted in settings DB; runtime state in memory.

const watchdog = {
  // config (loaded from DB on startup, saved on POST /api/watchdog)
  enabled:      false,
  skyTempLimit: -4,    // °C  — park if IR sky > this
  sqLimit:      16,    // mag/arcsec² — park if SQ < this (0 = disabled)
  retentionMin: 10,    // minutes of clear sky before resuming

  // runtime state
  state:        "off",   // "off" | "clear" | "bad" | "recovering"
  clearSince:   null,    // Date when sky turned good again
  badSince:     null,    // Date when sky first went bad (this episode)
  lastCheck:    null,    // Date of last OCS check
  lastMsg:      "",
  parked:       false,   // did WE park the mount?
  coverClosed:  false,   // did WE close the dust cover during parking?
  _interval:    null,
};

// Load persisted watchdog config
(function loadWatchdogConfig() {
  try {
    const raw = getSetting("watchdog_config");
    if (raw) {
      const cfg = JSON.parse(raw);
      if (cfg.enabled      !== undefined) watchdog.enabled      = !!cfg.enabled;
      if (cfg.skyTempLimit !== undefined) watchdog.skyTempLimit = Number(cfg.skyTempLimit);
      if (cfg.sqLimit      !== undefined) watchdog.sqLimit      = Number(cfg.sqLimit);
      if (cfg.retentionMin !== undefined) watchdog.retentionMin = Number(cfg.retentionMin);
      if (cfg.morningParkHour !== undefined) setSetting("morning_park_hour", String(Number(cfg.morningParkHour)));
    }
  } catch (_) {}
  watchdog.state = watchdog.enabled ? "clear" : "off";
})();

async function runWatchdogCheck() {
  if (!watchdog.enabled) { watchdog.state = "off"; return; }

  watchdog.lastCheck = new Date();

  let fields;
  try {
    const { ok, fields: f } = await ocsStatus(DEFAULT_OCS_HOST);
    if (!ok) {
      watchdog.lastMsg = "OCS unreachable";
      if (watchdog.state === "off") watchdog.state = "clear"; // enabled but can't check
      return;
    }
    fields = f;
  } catch (e) {
    watchdog.lastMsg = `OCS error: ${e.message}`;
    if (watchdog.state === "off") watchdog.state = "clear";
    return;
  }

  const irSky = ocsNum(fields.wea_irsky);
  const sq    = ocsNum(fields.wea_sq);
  const rain  = ocsClean(fields.wea_rain).toLowerCase();

  // Evaluate sky quality
  const rainBad    = rain && rain !== "no" && rain !== "0" && rain !== "clear" && rain !== "dry";
  const irBad      = irSky !== null && irSky > watchdog.skyTempLimit;
  const sqBad      = watchdog.sqLimit > 0 && sq !== null && sq < watchdog.sqLimit;
  const conditionsBad = rainBad || irBad || sqBad;

  const reasons = [];
  if (rainBad) reasons.push(`rain=${rain}`);
  if (irBad)   reasons.push(`IR=${irSky}°C > ${watchdog.skyTempLimit}°C`);
  if (sqBad)   reasons.push(`SQ=${sq} < ${watchdog.sqLimit}`);
  watchdog.lastMsg = conditionsBad ? `Bad: ${reasons.join(", ")}` : `Clear (IR=${irSky ?? "?"}°C, SQ=${sq ?? "?"})`;

  // Rain radar forecast check — non-blocking, runs in background
  fetchRainForecast().then(fc => {
    watchdog.rainForecast = fc;
    if (fc.warnLevel === "imminent" && !conditionsBad) {
      console.warn(`[watchdog] Rain radar: IMMINENT (${fc.maxProb30}% in next 30 min) — prepare to close`);
      if (seqState.running) seqLog(`⚠ Rain radar: ${fc.maxProb30}% chance in next 30 min`, "warn");
    }
  }).catch(() => {});  // forecast failure never blocks the watchdog

  if (conditionsBad) {
    watchdog.clearSince = null;
    if (watchdog.state !== "bad") {
      watchdog.state    = "bad";
      watchdog.badSince = watchdog.badSince ?? new Date(); // keep first onset time
      console.log(`[watchdog] Conditions bad — ${watchdog.lastMsg}`);

      // ── Park mount (always, not just when sequence is running) ─────────────
      if (!watchdog.parked) {
        console.log("[watchdog] Parking mount due to bad conditions...");
        try {
          const cfg = seqState.ninaConfig ?? { host: DEFAULT_NINA_HOST, port: DEFAULT_NINA_PORT, protocol: DEFAULT_NINA_PROTOCOL };
          await fetch(`http://${cfg.host}:${cfg.port}/v2/api/equipment/mount/park`, {
            method: "GET", signal: AbortSignal.timeout(15000),
          }).catch(() => {});
          watchdog.parked = true;
          console.log("[watchdog] Mount parked.");
        } catch (e) {
          console.error("[watchdog] Park failed:", e.message);
        }
      }

      // ── Rain: command OCS roof close with retries (backup to OCS firmware) ─
      if (rainBad) {
        console.log("[watchdog] Rain detected — commanding OCS roof close (with retries)...");
        (async () => {
          for (let attempt = 1; attempt <= 3; attempt++) {
            try {
              await ocsGet(DEFAULT_OCS_HOST, { roof: "close" });
              console.log(`[watchdog] OCS roof close command sent ✓ (attempt ${attempt})`);
              break;
            } catch (e) {
              console.error(`[watchdog] OCS roof close attempt ${attempt} failed: ${e.message}`);
              if (attempt < 3) await new Promise(r => setTimeout(r, 10_000));
            }
          }
        })();
      }

      // ── Close dust cover (always — protect optics regardless of sequence state) ──
      const cfg = seqState.ninaConfig ?? { host: DEFAULT_NINA_HOST, port: DEFAULT_NINA_PORT, protocol: DEFAULT_NINA_PROTOCOL };
      try {
        const fd = await getCoverState(cfg).catch(() => null);
        const cs = normCoverState(fd?.CoverState);
        if (cs === 0) {
          console.log("[watchdog] Dust cover not present — skipping close.");
        } else if (cs === 1) {
          console.log("[watchdog] Dust cover already closed ✓");
          watchdog.coverClosed = true;
        } else {
          console.log(`[watchdog] Closing dust cover (current state: ${fd?.CoverState ?? "unknown"})...`);
          await callNinaLong(cfg, "/v2/api/equipment/flatdevice/set-cover", { closed: true }, 30000);
          const deadline = Date.now() + 30000;
          let closed = false;
          while (Date.now() < deadline) {
            await new Promise(r => setTimeout(r, 1000));
            const st = normCoverState((await getCoverState(cfg).catch(() => null))?.CoverState);
            if (st === 1) { closed = true; break; }
            if (st === 101) break;
          }
          watchdog.coverClosed = closed;
          console.log(closed ? "[watchdog] Dust cover closed ✓" : "[watchdog] Dust cover close confirmation timed out");
        }
      } catch (e) {
        console.error("[watchdog] Dust cover close failed:", e.message);
      }
    }
  } else {
    // Conditions good
    watchdog.badSince = null; // reset bad-duration counter
    if (watchdog.state === "bad" || watchdog.state === "recovering") {
      if (!watchdog.clearSince) watchdog.clearSince = new Date();
      const clearedMs = Date.now() - watchdog.clearSince.getTime();
      const neededMs  = watchdog.retentionMin * 60 * 1000;
      if (clearedMs >= neededMs) {
        // Retention passed — unpark if we parked it
        watchdog.state = "clear";
        if (watchdog.parked) {
          console.log("[watchdog] Retention passed — checking roof before unparking...");
          try {
            const { ok, fields } = await ocsStatus(DEFAULT_OCS_HOST);
            if (!ok || !fields) {
              console.warn("[watchdog] OCS unreachable — staying parked until OCS confirms roof open");
              watchdog.state = "bad";
              return;
            }
            const roofStatus = ocsClean(fields?.roof_sta ?? "").toLowerCase();
            if (roofStatus && !roofStatus.includes("open")) {
              console.warn(`[watchdog] Roof is closed (${ocsClean(fields.roof_sta)}) — skipping unpark, staying parked`);
              watchdog.state = "bad";
              return;
            }
          } catch {
            console.warn("[watchdog] OCS error during roof check — staying parked");
            watchdog.state = "bad";
            return;
          }
          console.log("[watchdog] Roof open — unparking mount...");
          try {
            const cfg = seqState.ninaConfig;
            if (cfg) {
              await fetch(`http://${cfg.host}:${cfg.port}/v2/api/equipment/mount/unpark`, {
                method: "GET", signal: AbortSignal.timeout(15000),
              }).catch(() => {});
            }
            watchdog.parked = false;
            console.log("[watchdog] Mount unparked — sequence will retry on next frame check.");
          } catch (e) {
            console.error("[watchdog] Unpark failed:", e.message);
          }

          // Re-open dust cover if we closed it during parking
          if (watchdog.coverClosed) {
            try {
              const cfg = seqState.ninaConfig;
              if (cfg) {
                // Safety: only open if OCS roof is confirmed open
                let roofOpen = true;
                try {
                  const { ok, fields } = await ocsStatus(DEFAULT_OCS_HOST);
                  if (ok) {
                    const roofRaw = ocsClean(fields?.roof_sta ?? "").toLowerCase();
                    if (roofRaw && !roofRaw.includes("open")) roofOpen = false;
                  }
                } catch { /* OCS unreachable — assume open */ }

                if (!roofOpen) {
                  console.log("[watchdog] Roof not confirmed open — keeping dust cover closed.");
                } else {
                  console.log("[watchdog] Opening dust cover after recovery...");
                  await callNinaLong(cfg, "/v2/api/equipment/flatdevice/set-cover", { closed: false }, 30000);
                  const deadline = Date.now() + 30000;
                  let opened = false;
                  while (Date.now() < deadline) {
                    await new Promise(r => setTimeout(r, 1000));
                    const st = normCoverState((await getCoverState(cfg).catch(() => null))?.CoverState);
                    if (st === 3) { opened = true; break; }
                    if (st === 101) break;
                  }
                  watchdog.coverClosed = !opened;
                  console.log(opened ? "[watchdog] Dust cover open ✓" : "[watchdog] Dust cover open confirmation timed out");
                }
              }
            } catch (e) {
              console.error("[watchdog] Dust cover re-open failed:", e.message);
            }
          }
        }
      } else {
        watchdog.state = "recovering";
        const remainMin = Math.ceil((neededMs - clearedMs) / 60000);
        watchdog.lastMsg += ` — recovering, ${remainMin} min left`;
      }
    } else {
      watchdog.state = "clear";
    }
  }
}

// Start / stop watchdog interval whenever enabled changes
function applyWatchdogInterval() {
  if (watchdog._interval) { clearInterval(watchdog._interval); watchdog._interval = null; }
  if (watchdog.enabled) {
    watchdog._interval = setInterval(runWatchdogCheck, 10 * 1000);  // safety check every 10s
    runWatchdogCheck(); // immediate first check
  }
}
applyWatchdogInterval();

// ── Mount park monitor ────────────────────────────────────────────────────────
// Polls mount AtPark state every 15 s. If the mount transitions to parked for
// ANY reason (watchdog, NINA sequence, manual), automatically closes the dust
// cover to protect the optics. Independent of the watchdog.

const mountMonitor = {
  lastAtPark: null,   // last known AtPark value (null = unknown)
  lastCheck:  null,   // Date of last poll
  _interval:  null,
};

async function runMountParkCheck() {
  mountMonitor.lastCheck = new Date();
  try {
    const cfg = seqState.ninaConfig ?? {
      host: DEFAULT_NINA_HOST, port: DEFAULT_NINA_PORT, protocol: DEFAULT_NINA_PROTOCOL,
    };
    const { equipmentInfo } = await getEquipmentInfo(cfg);
    const mount = equipmentInfo?.Mount;
    if (!mount?.Connected) { mountMonitor.lastAtPark = null; return; }

    const atPark = mount.AtPark === true;

    // Transition: not parked → parked
    if (atPark && mountMonitor.lastAtPark === false) {
      console.log("[mount-monitor] Mount just parked — checking dust cover...");
      try {
        const fd = await getCoverState(cfg).catch(() => null);
        const cs = normCoverState(fd?.CoverState);
        if (cs === 0) {
          console.log("[mount-monitor] Dust cover not present — skipping.");
        } else if (cs === 1) {
          console.log("[mount-monitor] Dust cover already closed ✓");
        } else {
          console.log(`[mount-monitor] Closing dust cover (state: ${fd?.CoverState ?? "unknown"})...`);
          await callNinaLong(cfg, "/v2/api/equipment/flatdevice/set-cover", { closed: true }, 30000);
          const deadline = Date.now() + 30000;
          let closed = false;
          while (Date.now() < deadline) {
            await new Promise(r => setTimeout(r, 1000));
            const st = normCoverState((await getCoverState(cfg).catch(() => null))?.CoverState);
            if (st === 1) { closed = true; break; }
            if (st === 101) break;
          }
          console.log(closed
            ? "[mount-monitor] Dust cover closed ✓"
            : "[mount-monitor] Dust cover close timed out — check manually!");
        }
      } catch (e) {
        console.error("[mount-monitor] Dust cover close failed:", e.message);
      }
    }

    mountMonitor.lastAtPark = atPark;
  } catch { /* NINA unreachable — skip silently */ }
}

mountMonitor._interval = setInterval(runMountParkCheck, 5 * 1000);
setTimeout(runMountParkCheck, 5000); // initial check shortly after startup

app.get("/api/watchdog", (req, res) => {
  res.json({
    enabled:        watchdog.enabled,
    skyTempLimit:   watchdog.skyTempLimit,
    sqLimit:        watchdog.sqLimit,
    retentionMin:   watchdog.retentionMin,
    morningParkHour: parseInt(getSetting("morning_park_hour") ?? "8", 10),
    state:        watchdog.state,
    lastCheck:    watchdog.lastCheck,
    lastMsg:      watchdog.lastMsg,
    parked:         watchdog.parked,
    coverClosed:    watchdog.coverClosed,
    badSince:       watchdog.badSince,
    safetyOverride: watchdog.safetyOverride,
    rainForecast:   watchdog.rainForecast || null,
    mountMonitor: {
      lastAtPark:   mountMonitor.lastAtPark,
      lastCheck:    mountMonitor.lastCheck ?? null,
      intervalMs:   5000,
    },
  });
});

app.post("/api/watchdog/safety-override", (req, res) => {
  watchdog.safetyOverride = !watchdog.safetyOverride;
  console.log(`[watchdog] Safety override ${watchdog.safetyOverride ? "ENABLED" : "disabled"}`);
  res.json({ success: true, safetyOverride: watchdog.safetyOverride });
});

app.post("/api/watchdog", (req, res) => {
  const b = req.body || {};
  watchdog.enabled      = !!b.enabled;
  watchdog.skyTempLimit = Number(b.skyTempLimit) ?? -4;
  watchdog.sqLimit      = Number(b.sqLimit)      ?? 16;
  watchdog.retentionMin = Math.max(1, Number(b.retentionMin) || 10);
  if (b.morningParkHour !== undefined) setSetting("morning_park_hour", String(Math.max(0, Math.min(12, Number(b.morningParkHour) || 8))));
  // Set state synchronously so the response reflects the new state immediately
  if (!watchdog.enabled) {
    watchdog.state  = "off";
    watchdog.parked = false;
  } else if (watchdog.state === "off") {
    watchdog.state = "clear"; // will be refined on next runWatchdogCheck
  }
  setSetting("watchdog_config", JSON.stringify({
    enabled: watchdog.enabled, skyTempLimit: watchdog.skyTempLimit,
    sqLimit: watchdog.sqLimit, retentionMin: watchdog.retentionMin,
  }));
  applyWatchdogInterval();
  res.json({
    success: true,
    watchdog: {
      enabled:         watchdog.enabled,
      skyTempLimit:    watchdog.skyTempLimit,
      sqLimit:         watchdog.sqLimit,
      retentionMin:    watchdog.retentionMin,
      morningParkHour: parseInt(getSetting("morning_park_hour") ?? "8", 10),
      state:           watchdog.state,
      lastCheck:       watchdog.lastCheck,
      lastMsg:         watchdog.lastMsg,
      parked:          watchdog.parked,
      badSince:        watchdog.badSince,
      safetyOverride:  watchdog.safetyOverride,
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
      error: "Device must be one of: all, mount, camera, filterwheel, focuser, flatdevice",
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

    const { result: infoResult, equipmentInfo, statuses } = await getEquipmentInfo(target);

    const allSuccess = connectionResults.every((result) => result.success);
    return res.status(allSuccess ? 200 : 502).json({
      success: allSuccess,
      target,
      results: connectionResults,
      devices: statuses,
      equipment: equipmentInfo,
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
    // Auto-connect filterwheel if needed (same logic as sequence)
    const { equipmentInfo: fwPre } = await getEquipmentInfo(target);
    if (!fwPre?.FilterWheel?.Connected) {
      const cr = await callNinaLong(target, "/v2/api/equipment/filterwheel/connect", {}, 15000);
      if (!isNinaApiSuccess(cr)) {
        return res.status(502).json({ success: false, error: "Filterwheel failed to connect" });
      }
    }

    // Send the change command with a longer timeout (wheel can take up to ~15s)
    const changeResult = await callNinaLong(
      target, "/v2/api/equipment/filterwheel/change-filter", { filterId }, 30000,
    );
    if (!isNinaApiSuccess(changeResult)) {
      return res.status(502).json({
        success: false,
        error: `NINA rejected filter change to ID ${filterId}: ${changeResult.body?.Error || changeResult.status}`,
      });
    }

    // Poll until the wheel reports the expected position (up to 15s)
    const POLL_MS = 500;
    const DEADLINE = Date.now() + 15000;
    let confirmed = false;
    let actualFilter = null;
    while (Date.now() < DEADLINE) {
      await new Promise(r => setTimeout(r, POLL_MS));
      try {
        const { equipmentInfo } = await getEquipmentInfo(target);
        actualFilter = equipmentInfo?.FilterWheel?.SelectedFilter;
        if (actualFilter && Number(actualFilter.Id) === filterId) { confirmed = true; break; }
      } catch { /* ignore transient poll errors */ }
    }

    const { statuses } = await getEquipmentInfo(target);
    return res.status(confirmed ? 200 : 502).json({
      success: confirmed,
      filterId,
      confirmed,
      actualFilter,
      devices: statuses,
      ...(!confirmed && { error: `Wheel did not reach slot ${filterId} within 15 s` }),
    });
  } catch (error) {
    return res.status(502).json({
      success: false,
      target,
      error: normalizeNinaError(error),
    });
  }
});

// ── Dust cover (flat device) manual actions ───────────────────────────────────

async function flatDeviceAction(req, res, action) {
  let target;
  try { target = normalizeTargetConfig(req.body); }
  catch (e) { return res.status(400).json({ success: false, error: e.message }); }
  try {
    // Safety gate: refuse to open the cover if the OCS roof is not open
    if (action === "open") {
      try {
        const { ok, fields } = await ocsStatus(DEFAULT_OCS_HOST);
        if (ok) {
          const roofRaw = ocsClean(fields?.roof_sta ?? "").toLowerCase();
          if (roofRaw && !roofRaw.includes("open")) {
            return res.status(403).json({
              success: false,
              error: `Roof is not open (OCS: "${ocsClean(fields.roof_sta)}") — cover open blocked for safety`,
              roofStatus: ocsClean(fields.roof_sta),
            });
          }
        }
      } catch { /* OCS unreachable — allow but let sequence/user decide */ }
    }
    await ensureFlatDeviceConnected(target);
    const r = await callNinaLong(
      target, "/v2/api/equipment/flatdevice/set-cover", { closed: action === "close" }, 30000,
    );
    if (!isNinaApiSuccess(r)) {
      return res.status(502).json({ success: false, error: r.body?.Error || `NINA returned ${r.status}` });
    }
    // Poll for final state
    const want = action === "open" ? 3 : 1;
    const deadline = Date.now() + 30000;
    let coverState = null;
    let coverStateNum = 100;
    while (Date.now() < deadline) {
      await new Promise(x => setTimeout(x, 1000));
      const fd = await getCoverState(target);
      coverState = fd?.CoverState;
      coverStateNum = normCoverState(coverState);
      if (coverStateNum === want) break;
      if (coverStateNum === 101) break; // error
    }
    const stateNames = { 0: "NotPresent", 1: "Closed", 2: "Moving", 3: "Open", 100: "Unknown", 101: "Error" };
    const confirmed = coverStateNum === want;
    return res.status(confirmed ? 200 : 502).json({
      success: confirmed,
      action,
      coverState,
      coverStateLabel: stateNames[coverStateNum] ?? String(coverState ?? coverStateNum),
      ...(!confirmed && { error: `Cover did not reach ${action === "open" ? "Open" : "Closed"} state within 30 s` }),
    });
  } catch (e) {
    return res.status(502).json({ success: false, error: normalizeNinaError(e) });
  }
}

app.post("/api/nina/actions/flatdevice/open",    (req, res) => flatDeviceAction(req, res, "open"));
app.post("/api/nina/actions/flatdevice/close",   (req, res) => flatDeviceAction(req, res, "close"));
app.post("/api/nina/actions/flatdevice/connect", async (req, res) => {
  let target;
  try { target = normalizeTargetConfig(req.body); }
  catch (e) { return res.status(400).json({ success: false, error: e.message }); }
  try {
    const cr = await callNinaLong(target, "/v2/api/equipment/flatdevice/connect", {}, 15000);
    const { equipmentInfo } = await getEquipmentInfo(target);
    const ok = isNinaApiSuccess(cr);
    return res.status(ok ? 200 : 502).json({
      success: ok,
      flatDevice: equipmentInfo?.FlatDevice ?? null,
      devices: equipmentInfo ? deviceStatusFromEquipmentInfo(equipmentInfo) : null,
      error: ok ? undefined : (cr.body?.Error || cr.body?.error || `NINA connect failed (${cr.status})`),
    });
  } catch (e) {
    return res.status(502).json({ success: false, error: normalizeNinaError(e) });
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
  lastAutofocusTime: (() => { const v = getSetting("lastAutofocusTime"); return v ? Number(v) : null; })(),
  error: null,
  manualMode: false,
  waitingForStep: null, // non-null when paused waiting for /api/sequence/next
  ninaConfig: null,    // set when sequence starts, used by abort to reach NINA
  // ── ToO (Target-of-Opportunity) interrupt state ──────────────────────────
  alertReady: getSetting("alertReady", "false") === "true",
  tooInterrupt: null,     // set to { alert } when a ToO fires; checked by checkAbort
  tooRunning: false,      // true while the ToO sub-sequence is executing
  tooResumeState: null,   // saved state of the interrupted target for resume
};

loadQueue();

// ── Daily noon reset ──────────────────────────────────────────────────────────
function dailyReset({ manual = false } = {}) {
  const before = seqState.queue.length;
  // Remove alert-sourced targets (they belonged to last night)
  const alertRemoved = seqState.queue.filter(t => t.source === "alert");
  seqState.queue = seqState.queue.filter(t => t.source !== "alert");
  // Reset done=false on all remaining manual targets so they run again tonight
  seqState.queue.forEach(t => { t.done = false; });
  saveQueue();
  const msg = `Daily reset: removed ${alertRemoved.length} alert target(s), reset ${seqState.queue.length} manual target(s). (${manual ? "manual trigger" : "scheduled noon"})`;
  console.log("[daily-reset]", msg);
  if (alertRemoved.length || seqState.queue.length) seqLog(msg, "info");
  return { removed: alertRemoved.map(t => t.name), kept: seqState.queue.map(t => t.name) };
}

function scheduleNoonReset() {
  const resetHourLocal = parseInt(getSetting("daily_reset_hour", "12"), 10);
  const now = new Date();
  // Compute next occurrence of resetHourLocal:00 local time (server TZ)
  const next = new Date(now);
  next.setHours(resetHourLocal, 0, 0, 0);
  if (next <= now) next.setDate(next.getDate() + 1); // already passed today → tomorrow
  const msUntil = next - now;
  console.log(`[daily-reset] Scheduled for ${next.toLocaleString()} (in ${Math.round(msUntil / 60000)} min)`);
  setTimeout(() => {
    dailyReset();
    setInterval(dailyReset, 24 * 60 * 60 * 1000); // every 24 h after first fire
  }, msUntil);
}
scheduleNoonReset();

// ── Per-filter autofocus cache (DB-backed) ────────────────────────────────────
const FILTER_AF_CACHE_MS = 2 * 60 * 60 * 1000; // 2 hours

// ── Temperature-based focuser position prediction ──────────────────────────
// Linear regression from May 2026 autofocus runs (current camera configuration).
// Formula: pos ≈ slope × FOCTEMP + intercept  (FOCTEMP in °C from focuser thermistor)
// R²≈0.60 for combined, ±20 steps typical error vs ~±100 steps depth-of-focus.
const AF_TEMP_REGRESSION = {
  G:   { slope: 10.83, intercept: 4033 },
  BP:  { slope:  9.73, intercept: 4056 },
  RP:  { slope: 12.08, intercept: 4021 },
  ALL: { slope: 11.09, intercept: 4034 }, // combined / fallback
};

/** Return predicted focuser position from thermistor temperature (configurable model).
 *  Formula: pos = slope × (FOCTEMP − 15) + ref15
 *  Parameters are read live from DB settings so the UI can update them without restart. */
function predictFocuserPos(filterName, foctemp) {
  if (!Number.isFinite(foctemp)) return null;
  const knownFilters = ["G", "BP", "RP"];
  const f = knownFilters.includes(filterName) ? filterName : null;
  const reg = AF_TEMP_REGRESSION[f] ?? AF_TEMP_REGRESSION.ALL;
  // Live DB read (SETTING_DEFAULTS provide the May 2026 values as fallback)
  const slope = parseFloat(f ? getSetting(`af_slope_${f}`) : null) || reg.slope;
  const ref15 = parseFloat(f ? getSetting(`af_ref15_${f}`) : null) || Math.round(reg.slope * 15 + reg.intercept);
  return Math.round(slope * (foctemp - 15) + ref15);
}

/** Read current focuser thermistor temperature (°C) from NINA live equipment info. */
async function getFocuserTemp(target) {
  try {
    const { equipmentInfo } = await getEquipmentInfo(target);
    const temp = equipmentInfo?.Focuser?.Temperature;
    return Number.isFinite(temp) ? temp : null;
  } catch {
    return null;
  }
}

const _afInsertStmt = db.prepare(
  `INSERT INTO autofocus_log (filter, position, ts) VALUES (?, ?, datetime('now'))`
);
const _afLatestStmt = db.prepare(
  `SELECT position, ts FROM autofocus_log WHERE filter = ? ORDER BY ts DESC LIMIT 1`
);

function saveAfResult(filterName, pos) {
  if (!Number.isFinite(pos)) return;
  _afInsertStmt.run(filterName, pos);
}

function getLatestAfResult(filterName) {
  const row = _afLatestStmt.get(filterName);
  if (!row) return null;
  return { pos: row.position, ts: new Date(row.ts + "Z").getTime() };
}

function getAfResultIfFresh(filterName, maxAgeMs = FILTER_AF_CACHE_MS) {
  const entry = getLatestAfResult(filterName);
  if (!entry) return null;
  if ((Date.now() - entry.ts) > maxAgeMs) return null;
  return entry;
}

async function stepMoveFocuserAbsolute(target, pos) {
  await callNinaLong(target, "/v2/api/equipment/focuser/move", { position: pos }, 30000);
  const deadline = Date.now() + 20000;
  while (Date.now() < deadline) {
    await new Promise(r => setTimeout(r, 500));
    const { equipmentInfo } = await getEquipmentInfo(target);
    if (!equipmentInfo?.Focuser?.IsMoving) break;
  }
}

// Restore the last N log lines from tonight's log file so the UI shows
// recent history immediately after a server restart.
(function restoreLogFromFile() {
  try {
    const SEQ_LOG_DIR_BOOT = path.join(__dirname, "logs");
    const d = new Date();
    if (d.getHours() < 12) d.setDate(d.getDate() - 1);
    const ymd = d.toISOString().slice(0, 10);
    const logPath = path.join(SEQ_LOG_DIR_BOOT, `seq-${ymd}.log`);
    if (!fs.existsSync(logPath)) return;
    const lines = fs.readFileSync(logPath, "utf8").trim().split("\n").filter(Boolean);
    const RESTORE = 200;
    for (const line of lines.slice(-RESTORE)) {
      // Format: 2026-04-07T00:21:30.416Z  [INFO]  message
      const m = line.match(/^(\S+)\s+\[(\w+)\]\s+(.*)/);
      if (!m) continue;
      const [, iso, lvl, msg] = m;
      const ts = new Date(iso).toLocaleTimeString("en-GB", { hour12: false });
      seqState.log.push({ ts, msg: msg.trim(), level: lvl.toLowerCase() });
    }
    if (seqState.log.length) console.log(`Restored ${seqState.log.length} log line(s) from ${logPath}`);
  } catch { /* non-fatal */ }
})();

// Resolvers for the manual-step confirmation promise
let _confirmResolve = null;
let _confirmReject  = null;

const SEQ_MAX_LOG = 500;

// ── Persistent sequence log file ─────────────────────────────────────────────
// One file per night (YYYY-MM-DD), written to server/logs/seq-YYYY-MM-DD.log
const SEQ_LOG_DIR = path.join(__dirname, "logs");
fs.mkdirSync(SEQ_LOG_DIR, { recursive: true });

function seqLogFile() {
  const d = new Date();
  // Use the local date, but treat nights as starting at noon (so 01:00 still
  // belongs to "last night"). If hour < 12, roll back one day.
  const midnight = new Date(d);
  if (d.getHours() < 12) midnight.setDate(midnight.getDate() - 1);
  const ymd = midnight.toISOString().slice(0, 10);
  return path.join(SEQ_LOG_DIR, `seq-${ymd}.log`);
}

function seqLog(msg, level = "info") {
  const now = new Date();
  const ts  = now.toLocaleTimeString("en-GB", { hour12: false });
  const entry = { ts, msg, level };
  seqState.log.push(entry);
  if (seqState.log.length > SEQ_MAX_LOG) seqState.log.shift();
  const consoleLine = `[SEQ ${level.toUpperCase()}] ${ts} ${msg}`;
  console.log(consoleLine);
  // Append to nightly log file (non-blocking, errors silently ignored)
  const iso  = now.toISOString();
  const line = `${iso}  [${level.toUpperCase().padEnd(4)}]  ${msg}\n`;
  fs.appendFile(seqLogFile(), line, () => {});
}

function checkAbort() {
  if (seqState.aborted) throw new Error("__ABORTED__");
  if (seqState.tooInterrupt && !seqState.tooRunning) {
    const err = new Error("__ALERT_INTERRUPT__");
    err.code = "__ALERT_INTERRUPT__";
    err.alert = seqState.tooInterrupt;
    throw err;
  }
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
  await callNinaLong(target, "/v2/api/equipment/camera/cool", {
    temperature: targetTempC,
    minutes: 0,
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

/**
 * Verify the OCS roof is open before any mount movement.
 * Throws if roof is confirmed closed OR if OCS is unreachable (fail-safe).
 * @param {string} context  – short label shown in the error/log (e.g. "unpark", "slew")
 */
async function assertRoofOpen(context = "mount move") {
  let ok, fields;
  try {
    ({ ok, fields } = await ocsStatus(DEFAULT_OCS_HOST));
  } catch (e) {
    const err = new Error(`OCS unreachable — cannot verify roof before ${context}, aborting for safety`);
    err.code = "__ROOF_CLOSED__";
    seqLog(`✗ ${err.message}`, "error");
    throw err;
  }
  if (!ok) {
    const err = new Error(`OCS unreachable — cannot verify roof before ${context}, aborting for safety`);
    err.code = "__ROOF_CLOSED__";
    seqLog(`✗ ${err.message}`, "error");
    throw err;
  }
  const roofStatus = ocsClean(fields?.roof_sta ?? "").toLowerCase();
  if (roofStatus && !roofStatus.includes("open")) {
    const err = new Error(`Roof is closed (OCS: "${ocsClean(fields.roof_sta)}") — aborting ${context} for safety`);
    err.code = "__ROOF_CLOSED__";
    seqLog(`✗ ${err.message}`, "error");
    throw err;
  }
}

async function stepUnparkMount(target) {
  seqState.currentStep = "Preparing mount...";
  await assertRoofOpen("unpark");
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

/**
 * Precess J2000 coordinates to JNOW (approximate, IAU 1976).
 * Returns { ra, dec } in degrees for the current epoch.
 */
function _j2000ToJnow(raDeg, decDeg, date = new Date()) {
  const T = (date.getFullYear() + (date.getMonth() + 1) / 12 - 2000) / 100;
  const d2r = _D2R, r2d = _R2D;
  const zetaA  = ((2306.2181 + (1.39656 - 0.000139 * T) * T) * T
                + (0.30188 - 0.000345 * T) * T * T + 0.017998 * T * T * T) * d2r / 3600;
  const zA     = ((2306.2181 + (1.39656 - 0.000139 * T) * T) * T
                + (1.09468  + 0.000066 * T) * T * T + 0.018203 * T * T * T) * d2r / 3600;
  const thetaA = ((2004.3109 - (0.85330 + 0.000217 * T) * T) * T
                - (0.42665  + 0.000217 * T) * T * T - 0.041775 * T * T * T) * d2r / 3600;
  const ra0  = raDeg  * d2r;
  const dec0 = decDeg * d2r;
  const A = Math.cos(dec0) * Math.sin(ra0 + zetaA);
  const B = Math.cos(thetaA) * Math.cos(dec0) * Math.cos(ra0 + zetaA) - Math.sin(thetaA) * Math.sin(dec0);
  const C = Math.sin(thetaA) * Math.cos(dec0) * Math.cos(ra0 + zetaA) + Math.cos(thetaA) * Math.sin(dec0);
  return {
    ra:  ((Math.atan2(A, B) * r2d) + zA * r2d + 360) % 360,
    dec: Math.asin(C) * r2d,
  };
}

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
  const haRaw = (lst - raDeg + 360) % 360;                  // 0–360°
  const haSigned = haRaw > 180 ? haRaw - 360 : haRaw;       // −180 … +180° (+ = West)
  const haR  = haRaw * _D2R;
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

  return { alt: altR * _R2D, az: azR * _R2D, ha: haSigned };
}

/**
 * Mount safety check: called before every frame exposure.
 * Throws an error with code "__LIMIT_REACHED__" if:
 *   - mount tracking is OFF
 *   - current altitude < minAlt  (too low)
 *   - current altitude > zenithLimit
 *   - |current HA| < meridianGap  (too close to meridian)
 *
 * @param {object} target       – NINA config for getEquipmentInfo
 * @param {object} limits       – { minAlt, meridianGap, zenithLimit }
 */
async function checkMountSafetyForFrame(target, { minAlt = 20, meridianGap = 10, zenithLimit = 70 } = {}) {
  // ── Morning cutoff: end sequence if past configured local hour ─────────────
  // Checked here (not just in waitForWatchdogClear) so it fires even on clear nights.
  {
    const morningParkHour = parseInt(getSetting("morning_park_hour") ?? "8", 10);
    const localHour = new Date().getHours();
    if (localHour >= morningParkHour && localHour < 14) {
      seqLog(`🌅 Morning cutoff (${morningParkHour}:00 local) — ending sequence and parking`, "warn");
      // Park mount
      try {
        await fetch(`http://${target.host}:${target.port}/v2/api/equipment/mount/park`, {
          method: "GET", signal: AbortSignal.timeout(15000),
        }).catch(() => {});
      } catch { /* non-fatal */ }
      seqState.aborted = true;
      checkAbort(); // throws __ABORTED__
    }
  }

  // Watchdog takes priority: if sky is bad, wait here — do NOT throw __LIMIT_REACHED__
  // because that would consume a retry slot; instead block until the watchdog clears.
  if (watchdog.enabled && (watchdog.state === "bad" || watchdog.state === "recovering")) {
    seqLog("  Safety check: watchdog reports bad sky — holding frame until sky clears", "warn");
    while (watchdog.state === "bad" || watchdog.state === "recovering") {
      checkAbort();
      const morningParkHour = parseInt(getSetting("morning_park_hour") ?? "8", 10);
      const localHour = new Date().getHours();
      if (localHour >= morningParkHour && localHour < 14) {
        seqLog(`🌅 Morning cutoff (${morningParkHour}:00 local) — ending sequence, mount already parked`, "warn");
        seqState.aborted = true;
        checkAbort();
      }
      await new Promise(r => setTimeout(r, 60_000));
    }
    seqLog("  Safety check: watchdog cleared — proceeding with frame");
  }

  let mount;
  {
    let lastErr;
    for (let attempt = 1; attempt <= 3; attempt++) {
      try {
        const { equipmentInfo } = await getEquipmentInfo(target);
        mount = equipmentInfo?.Mount;
        lastErr = null;
        break;
      } catch (e) {
        lastErr = e;
        if (attempt < 3) await new Promise(r => setTimeout(r, 5000));
      }
    }
    if (lastErr) {
      seqLog("  Safety check: could not reach NINA after 3 attempts — skipping check", "warn");
      return; // non-fatal: proceed with exposure if NINA is unreachable
    }
  }

  if (!mount) {
    seqLog("  Safety check: no mount info — skipping check", "warn");
    return;
  }

  // a) Tracking check — but only if watchdog did NOT park us (it will unpark autonomously)
  if (mount.TrackingEnabled === false && !watchdog.parked) {
    const err = new Error("Mount tracking is OFF — pausing target");
    err.code = "__LIMIT_REACHED__";
    throw err;
  }
  // If watchdog parked it but state is now clear, simply return — arbiter loop will re-slew
  if (mount.TrackingEnabled === false && watchdog.parked) {
    seqLog("  Safety check: mount parked by watchdog but sky just cleared — arbiter will re-slew", "warn");
    return;
  }

  // b) Altitude / HA checks (require live mount position and observer site)
  const lat    = mount.SiteLatitude;
  const lon    = mount.SiteLongitude;
  const raHrs  = mount.RightAscension;  // NINA returns RA in hours
  const decDeg = mount.Declination;

  if (lat == null || lon == null || raHrs == null || decDeg == null) {
    seqLog("  Safety check: missing mount position data — skipping alt/HA check", "warn");
    return;
  }

  const raDeg = raHrs * 15; // hours → degrees
  const { alt, ha } = serverComputeAltAz(raDeg, decDeg, lat, lon, new Date());

  // c) Minimum altitude
  if (alt < minAlt) {
    const err = new Error(
      `Target too low: ${alt.toFixed(1)}° < min altitude ${minAlt}° — pausing target`
    );
    err.code = "__LIMIT_REACHED__";
    throw err;
  }

  // d) Zenith limit
  if (alt > zenithLimit) {
    const err = new Error(
      `Target too high: ${alt.toFixed(1)}° > zenith limit ${zenithLimit}° — pausing target`
    );
    err.code = "__LIMIT_REACHED__";
    throw err;
  }

  // e) Meridian gap: |HA| < meridianGap means the target is within the forbidden zone
  //    around the meridian (HA=0). Positive HA = West, negative = East.
  if (Math.abs(ha) < meridianGap) {
    const side = ha >= 0 ? "W" : "E";
    const err = new Error(
      `Target within meridian gap: HA ${ha.toFixed(1)}°${side} < ${meridianGap}° — pausing target`
    );
    err.code = "__LIMIT_REACHED__";
    throw err;
  }

  seqLog(
    `  Safety ✓  Alt ${alt.toFixed(1)}° (${minAlt}°–${zenithLimit}°)  |HA| ${Math.abs(ha).toFixed(1)}° (gap ${meridianGap}°)  Tracking: ON`,
  );
}

/**
 * Among a list of candidate targets, return the one currently at the highest
 * altitude that is inside all safety limits. Returns null if none qualify.
 *
 * Scoring: primary = altitude (highest first, i.e. best sky conditions).
 * Tie-break: "urgency" — prefer a target that is closest to leaving the valid
 * window (setting or entering the forbidden zone), so we don't miss it.
 *
 * @param {Array}  candidates  – target objects with { raDeg, decDeg, name }
 * @param {number} lat         – observer latitude (degrees)
 * @param {number} lon         – observer longitude (degrees)
 * @param {object} limits      – { minAlt, zenithLimit, meridianGap }
 * @returns {{ target, alt, ha, reason }|null}
 */
function pickBestTarget(candidates, lat, lon, { minAlt = 20, zenithLimit = 70, meridianGap = 10 } = {}) {
  const now = new Date();

  const scored = candidates
    .filter(t => t.raDeg != null && t.decDeg != null)
    .map(t => {
      const { alt, ha } = serverComputeAltAz(t.raDeg, t.decDeg, lat, lon, now);
      const absHA = Math.abs(ha);
      const valid = alt >= minAlt && alt <= zenithLimit && absHA >= meridianGap;

      // Estimate minutes until the target leaves the valid window:
      //   either sets below minAlt, rises above zenithLimit, or enters meridian gap.
      // We sample 5-min steps forward (up to 4 h) to find the first violation.
      let minutesUntilInvalid = Infinity;
      if (valid) {
        for (let m = 5; m <= 240; m += 5) {
          const future = new Date(now.getTime() + m * 60000);
          const { alt: fAlt, ha: fHA } = serverComputeAltAz(t.raDeg, t.decDeg, lat, lon, future);
          if (fAlt < minAlt || fAlt > zenithLimit || Math.abs(fHA) < meridianGap) {
            minutesUntilInvalid = m;
            break;
          }
        }
      }

      return { target: t, alt, ha, absHA, valid, minutesUntilInvalid };
    })
    .filter(s => s.valid);

  if (!scored.length) return null;

  // Sort: highest altitude first; among ties (< 2°), prefer the most urgent (sets soonest)
  scored.sort((a, b) => {
    if (Math.abs(b.alt - a.alt) > 2) return b.alt - a.alt;       // altitude wins
    return a.minutesUntilInvalid - b.minutesUntilInvalid;         // urgency tie-break
  });

  const best = scored[0];
  return {
    target:            best.target,
    alt:               best.alt,
    ha:                best.ha,
    minutesUntilValid: 0,
    minutesUntilInvalid: best.minutesUntilInvalid,
    allScores:         scored.map(s => ({
      name: s.target.name,
      alt:  s.alt.toFixed(1),
      ha:   s.ha.toFixed(1),
      minsLeft: s.minutesUntilInvalid === Infinity ? ">4h" : `${s.minutesUntilInvalid}min`,
    })),
  };
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
/**
 * Block until the watchdog reports clear sky (state !== "bad" / "recovering").
 * Called from stepSlewWithCloudRetry, checkMountSafetyForFrame, and runSequence.
 */
async function waitForWatchdogClear() {
  if (watchdog.state !== "bad" && watchdog.state !== "recovering") return;
  seqLog("⛅ Watchdog: bad conditions — sequence suspended until sky clears", "warn");
  while (watchdog.state === "bad" || watchdog.state === "recovering") {
    checkAbort();
    // Morning cutoff: abort once local time reaches the configured park hour
    // (mount is already parked by the watchdog loop, this just ends the sequence)
    const morningParkHour = parseInt(getSetting("morning_park_hour") ?? "8", 10);
    const localHour = new Date().getHours();
    if (localHour >= morningParkHour && localHour < 14) {
      seqLog(`🌅 Watchdog: past morning cutoff (${morningParkHour}:00 local) — ending sequence, mount already parked`, "warn");
      seqState.aborted = true;
      checkAbort(); // throws __ABORTED__
    }
    const badMin = watchdog.badSince ? Math.round((Date.now() - watchdog.badSince.getTime()) / 60_000) : "?";
    seqState.currentStep = `Watchdog: waiting for clear sky — bad for ${badMin} min (${watchdog.lastMsg || watchdog.state})`;
    await new Promise(r => setTimeout(r, 60_000));
  }
  seqLog("☀ Watchdog: sky clear — resuming", "info");
}

/**
 * Slew + center with automatic cloud-pause retry.
 *
 * NINA's centering endpoint runs capture→solve→sync→re-slew internally.
 * If the solve FAILS (no stars visible = clouds) NINA returns a failure
 * response. We detect that, wait cloudWaitMs, then retry.
 * Large offset after a successful solve is NOT a cloud indicator — it just
 * means the mount was off and NINA corrected it.
 *
 * After maxRetries consecutive solve failures the error propagates so the
 * arbiter can pause this target and move to the next one.
 */
async function stepSlewWithCloudRetry(
  target, raDeg, decDeg, name,
  { center = false, cloudWaitMs = 2 * 60 * 1000, maxRetries = 4 } = {}
) {
  if (!center) {
    // No plate solve requested — plain slew, no retry logic needed
    await stepSlewToTarget(target, raDeg, decDeg, name, { center: false });
    return;
  }

  // Watchdog gate before slewing
  await waitForWatchdogClear();
  checkAbort();

  // Delegate entirely to stepSlewToTarget — it handles the full slew+solve
  // cycle and polls until the mount converges. NINA is the authority on
  // whether the plate solve succeeded; we don't second-guess it here.
  // Cloud detection happens per-frame via the star count gate.
  await stepSlewToTarget(target, raDeg, decDeg, name, { center });
}

async function stepSlewToTarget(target, raDeg, decDeg, name, { center = false } = {}) {
  await assertRoofOpen("slew");
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
    // Fire the SlewAndCenter command. NINA will: goto → capture → plate solve →
    // sync → re-slew → repeat until within its own tolerance.
    let centeringOk = false;
    try {
      const centerResult = await callNinaLong(target, "/v2/api/equipment/mount/slew", {
        ra: raDeg, dec: decDeg, waitForResult: true, center: true, name,
      }, 10 * 60 * 1000);
      centeringOk = isNinaApiSuccess(centerResult);
      // Don't warn yet — we'll check everExposed after polling to distinguish
      // a real plate-solve failure from NINA's spurious "Slew failed" on plain gotos.
    } catch (fetchErr) {
      seqLog(`  Centering request error: ${fetchErr.message} — falling back to plain slew`, "warn");
      await plainSlew();
      seqLog("  Plain slew complete (no centering)");
    }

    // Poll until NINA's full SlewAndCenter cycle completes.
    // NINA iterates: goto → capture plate-solve image (IsExposing=true) →
    // solve → correction-slew (Slewing=true) → repeat until within tolerance.
    // Done when BOTH mount.Slewing=false AND camera.IsExposing=false for SETTLE_STABLE ms.
    // NOTE: we do NOT verify the final offset here because NINA reports coordinates
    // in JNOW while our targets are stored as J2000. The precession for 2026 is
    // ~10-22 arcmin depending on Dec — exactly the "offsets" we used to misdiagnose
    // as plate-solve failures. Trust NINA's own centering loop.
    seqState.currentStep = `Centering ${name} — waiting for convergence...`;
    const SETTLE_TIMEOUT  = 10 * 60 * 1000;
    const POLL_INTERVAL   = 2500;
    const SETTLE_STABLE   = 10000;  // 10s both-idle → centering done
    const deadline  = Date.now() + SETTLE_TIMEOUT;
    let stableMs    = 0;
    let lastLogRA   = null;

    while (Date.now() < deadline) {
      checkAbort();
      await new Promise(r => setTimeout(r, POLL_INTERVAL));

      const { equipmentInfo: ei } = await getEquipmentInfo(target);
      const m = ei?.Mount;
      if (!m) continue;

      // Log position changes (informational only — JNOW coords, not comparable to J2000 target)
      const raLog = m.RightAscensionString;
      if (raLog !== lastLogRA) {
        seqLog(`  Centering: mount at ${raLog} / ${m.DeclinationString}  (Slewing=${m.Slewing}, Exposing=${ei?.Camera?.IsExposing ?? "?"})`);
        seqState.currentStep = `Centering ${name} — Slewing=${m.Slewing} Exposing=${ei?.Camera?.IsExposing ?? "?"}`;
        lastLogRA = raLog;
      }

      const mountBusy  = m.Slewing === true;
      const cameraBusy = ei?.Camera?.IsExposing === true;
      if (mountBusy || cameraBusy) {
        stableMs = 0;  // still in slew or plate-solve image capture
      } else {
        stableMs += POLL_INTERVAL;
        if (stableMs >= SETTLE_STABLE) break;  // both idle for 10s → done
      }
    }

    // NINA's centering loop is complete. Compare mount JNOW position to
    // target converted to JNOW — this is the definitive centering quality check.
    const { equipmentInfo: eiFinal } = await getEquipmentInfo(target).catch(() => ({ equipmentInfo: null }));
    const mFinal = eiFinal?.Mount;
    const finalPos = mFinal ? `${mFinal.RightAscensionString} / ${mFinal.DeclinationString}` : "unknown";

    if (mFinal?.RightAscension != null && mFinal?.Declination != null) {
      const targetJnow = _j2000ToJnow(raDeg, decDeg);
      const mountRaDeg = mFinal.RightAscension * 15;  // NINA RA is in hours
      const offsetArcsec = _angularSepArcsec(mountRaDeg, mFinal.Declination, targetJnow.ra, targetJnow.dec);
      const MAX_CENTER_OFFSET = 120; // arcsec — good centering
      if (offsetArcsec <= MAX_CENTER_OFFSET) {
        seqLog(`Slew + centering complete ✓  offset ${Math.round(offsetArcsec)}"  mount at ${finalPos} (JNOW)`);
      } else {
        seqLog(`  ⚠ Centering offset ${Math.round(offsetArcsec)}" (limit ${MAX_CENTER_OFFSET}") — plate solve may have failed (clouds?)`, "warn");
        seqLog(`    Watchdog will suspend imaging if sky is too cloudy`, "warn");
        seqLog(`Slew + centering complete ⚠  offset ${Math.round(offsetArcsec)}"  mount at ${finalPos} (JNOW)`);
      }
    } else {
      seqLog(`Slew + centering complete ✓  mount at ${finalPos} (JNOW)`);
    }

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

async function stepStartGuiding(target, opts = {}) {
  const MAX_ATTEMPTS   = 3;
  const SETTLE_TIMEOUT = 3 * 60 * 1000;  // 3 min per attempt to reach "Guiding"
  const RETRY_PAUSE    = 15_000;          // 15 s between stop → restart
  const allowLostLock  = opts.allowLostLock === true;

  for (let attempt = 1; attempt <= MAX_ATTEMPTS; attempt++) {
    seqState.currentStep = attempt === 1
      ? "Starting PHD2 guiding..."
      : `PHD2 guiding — retry ${attempt}/${MAX_ATTEMPTS}...`;
    if (attempt === 1) seqLog("Starting PHD2 guiding...");
    else               seqLog(`PHD2 guiding retry ${attempt}/${MAX_ATTEMPTS} — stopping then restarting...`, "warn");

    // Stop guiding on retries so PHD2 re-calibrates / re-acquires star
    if (attempt > 1) {
      await callNinaLong(target, "/v2/api/equipment/guider/stop", {}, 30_000).catch(() => {});
      await new Promise(r => setTimeout(r, RETRY_PAUSE));
      checkAbort();
    }

    const result = await callNinaLong(target, "/v2/api/equipment/guider/start", {}, 2 * 60 * 1000);
    if (!isNinaApiSuccess(result)) {
      seqLog(`  Guiding start note: ${JSON.stringify(result.body)}`, "warn");
    }

    seqState.currentStep = `Waiting for PHD2 to settle (attempt ${attempt}/${MAX_ATTEMPTS})...`;
    seqLog(`Waiting for guiding to settle...`);

    try {
      await pollEquipment(
        target,
        (info) => {
          const s = info?.Guider?.State;
          return s === "Guiding" || (allowLostLock && s === "LostLock");
        },
        SETTLE_TIMEOUT,
        3000,
      );
      seqLog("Guiding active ✓");
      return true; // success
    } catch {
      if (attempt < MAX_ATTEMPTS) {
        seqLog(`  PHD2 did not settle within ${SETTLE_TIMEOUT / 60000} min — retrying...`, "warn");
      } else {
        seqLog(`⚠ PHD2 failed to settle after ${MAX_ATTEMPTS} attempts — continuing without guiding`, "warn");
      }
    }
  }
  return false;
}

// ── Dust Cover (Alnitak Flip-Flat / NINA FlatDevice) ─────────────────────────
// CoverState: 1=Closed  2=Moving  3=Open  0=NotPresent  100=Unknown  101=Error

async function getCoverState(target) {
  const { equipmentInfo } = await getEquipmentInfo(target);
  return equipmentInfo?.FlatDevice ?? null;
}

async function ensureFlatDeviceConnected(target) {
  const fd = await getCoverState(target);
  if (!fd?.Connected) {
    seqLog("Flat device not connected — connecting...");
    const cr = await callNinaLong(target, "/v2/api/equipment/flatdevice/connect", {}, 15000);
    if (!isNinaApiSuccess(cr)) throw new Error("Flat device failed to connect");
    seqLog("Flat device connected ✓");
  }
}

/** Normalise CoverState: NINA may return "Open"/"Closed"/… string or ASCOM int */
function normCoverState(raw) {
  if (typeof raw === "number") return raw;
  // NINA serialises the enum to its name
  const map = { notpresent: 0, closed: 1, moving: 2, open: 3, unknown: 100, error: 101 };
  return map[String(raw).toLowerCase().trim()] ?? 100;
}

function isCoverDisabled() {
  return getSetting("coverDisabled", "false") === "true";
}

async function stepOpenCover(target) {
  if (isCoverDisabled()) {
    seqLog("Dust cover disabled (maintenance mode) — skipping open");
    return;
  }
  seqState.currentStep = "Opening dust cover...";
  seqLog("Opening dust cover...");
  await ensureFlatDeviceConnected(target);

  const r = await callNinaLong(target, "/v2/api/equipment/flatdevice/set-cover", { closed: false }, 30000);
  if (!isNinaApiSuccess(r)) throw new Error(`Cover open command rejected: ${r.body?.Error || r.status}`);

  // Poll until CoverState = Open (3), up to 30 s
  const deadline = Date.now() + 30000;
  while (Date.now() < deadline) {
    await new Promise(res => setTimeout(res, 1000));
    const fd = await getCoverState(target);
    const cs = normCoverState(fd?.CoverState);
    if (cs === 3) { seqLog("Dust cover open ✓"); return; }
    if (cs === 101) throw new Error("Cover reported error state");
  }
  seqLog("Cover open confirmation timed out — continuing anyway", "warn");
}

async function stepCloseCover(target, { strict = true } = {}) {
  if (isCoverDisabled()) {
    seqLog("Dust cover disabled (maintenance mode) — skipping close");
    return true;
  }
  seqState.currentStep = "Closing dust cover...";
  seqLog("Closing dust cover...");
  await ensureFlatDeviceConnected(target);

  const r = await callNinaLong(target, "/v2/api/equipment/flatdevice/set-cover", { closed: true }, 30000);
  if (!isNinaApiSuccess(r)) throw new Error(`Cover close command rejected: ${r.body?.Error || r.status}`);

  // Poll until CoverState = Closed (1), up to 30 s
  const deadline = Date.now() + 30000;
  let confirmed = false;
  while (Date.now() < deadline) {
    await new Promise(res => setTimeout(res, 1000));
    const fd = await getCoverState(target);
    const cs = normCoverState(fd?.CoverState);
    if (cs === 1) { confirmed = true; break; }
    if (cs === 101) {
      if (strict) throw new Error("Cover reported error state during close");
      seqLog("Cover error state during close", "warn");
      break;
    }
  }
  if (confirmed) {
    seqLog("Dust cover closed ✓");
  } else {
    const msg = "Cover close confirmation timed out — verify manually before leaving!";
    if (strict) seqLog(msg, "error"); else seqLog(msg, "warn");
  }
  return confirmed;
}

async function stepParkMount(target) {
  seqState.currentStep = "Parking mount...";
  seqLog("Parking mount...");
  await callNinaLong(target, "/v2/api/equipment/mount/park", {}, 60000);
  // Poll up to 30s for AtPark=true
  const deadline = Date.now() + 30000;
  while (Date.now() < deadline) {
    await new Promise(r => setTimeout(r, 2000));
    const { equipmentInfo } = await getEquipmentInfo(target);
    if (equipmentInfo?.Mount?.AtPark === true) { seqLog("Mount parked ✓"); return; }
  }
  seqLog("Park confirmation timed out — verify manually", "warn");
}

/**
 * Run autofocus for a specific filter.
 * - Changes the filterwheel to filterName first (so NINA AF uses the right filter).
 * - If a fresh DB entry exists (< FILTER_AF_CACHE_MS old), just moves the focuser
 *   to that position instead of running a full AF sequence.
 * - On a full run, saves the result to autofocus_log DB table for future reuse.
 * availableFilters: the filterwheel filter list from NINA equipment info.
 */
async function stepAutofocusForFilter(target, filterName, availableFilters) {
  // ── 1. Change to this filter before AF so NINA uses the right bandpass ────
  const filter = findFilterByName(availableFilters, filterName);
  if (filter) {
    await changeFilterVerified(target, filter);
    checkAbort();
  }

  // ── 2. Check DB cache ──────────────────────────────────────────────────────
  const cached = getAfResultIfFresh(filterName);
  if (cached) {
    const ageMin = Math.round((Date.now() - cached.ts) / 60000);
    seqLog(`AF [${filterName}]: DB pos ${cached.pos} (${ageMin} min ago) — moving focuser ✓`);
    seqState.currentStep = `AF [${filterName}]: moving to cached position ${cached.pos}...`;
    await stepMoveFocuserAbsolute(target, cached.pos);
    return false; // cache hit — guiding was not disturbed
  }

  // ── 3. Full AF run ────────────────────────────────────────────────────────
  seqState.currentStep = `Running autofocus [${filterName}]...`;
  seqLog(`Running autofocus [${filterName}]...`);

  const { equipmentInfo: infoBefore } = await getEquipmentInfo(target);
  const startPos = infoBefore?.Focuser?.Position;
  seqLog(`  Focuser start position: ${startPos ?? "unknown"}`);

  const result = await callNinaLong(target, "/v2/api/equipment/focuser/auto-focus", {}, 15000);
  if (!isNinaApiSuccess(result)) {
    seqLog(`  AF trigger failed: ${result.body?.Error || result.status} — skipping`, "warn");
    return;
  }
  seqLog(`  ${result.body?.Response ?? "AF started"} — monitoring position...`);

  await new Promise(r => setTimeout(r, 5000)); // brief pause to let AF motion begin
  checkAbort();

  const AF_DEADLINE   = Date.now() + 12 * 60 * 1000;
  const STABLE_NEEDED = 6;     // 6 × 5 s = 30 s stable → AF done
  const MIN_MOVES     = 3;
  const MIN_AF_MS     = 60 * 1000; // min 60s before declaring done
  const afStartMs     = Date.now();
  const NO_MOVE_LIMIT = 30;

  let lastPos     = startPos;
  let moveCount   = 0;
  let stableCount = 0;
  let pollCount   = 0;
  let lastLogPos  = startPos;
  let finalPos    = startPos;

  while (Date.now() < AF_DEADLINE) {
    checkAbort();
    await new Promise(r => setTimeout(r, 5000));
    pollCount++;

    try {
      const { equipmentInfo } = await getEquipmentInfo(target);
      const pos        = equipmentInfo?.Focuser?.Position;
      const isMoving   = equipmentInfo?.Focuser?.IsMoving;
      const isExposing = equipmentInfo?.Camera?.IsExposing;

      seqState.currentStep = `Autofocus [${filterName}]: pos ${pos ?? "?"}${isMoving ? " ●" : ""}${isExposing ? " 📷" : ""}`;
      finalPos = pos ?? finalPos;

      if (pos !== lastPos) {
        moveCount++;
        stableCount = 0;
        if (lastLogPos == null || Math.abs((pos ?? 0) - (lastLogPos ?? 0)) >= 5) {
          seqLog(`  Focuser move #${moveCount}: ${lastPos} → ${pos}`);
          lastLogPos = pos;
        }
        lastPos = pos;
      } else if (!isMoving && !isExposing) {
        stableCount++;
      } else {
        stableCount = 0;
      }

      const elapsed = Date.now() - afStartMs;
      if (moveCount >= MIN_MOVES && stableCount >= STABLE_NEEDED && elapsed >= MIN_AF_MS) {
        seqLog(`  Autofocus [${filterName}] complete — position: ${startPos} → ${pos} (${moveCount} moves) ✓`);
        saveAfResult(filterName, pos);
        return true; // full AF ran
      }

      if (moveCount === 0 && pollCount >= NO_MOVE_LIMIT) {
        seqLog(`  No focuser movement — assuming AF skipped/complete [${filterName}]`);
        if (Number.isFinite(finalPos)) saveAfResult(filterName, finalPos);
        return true; // full AF ran (or was skipped by NINA)
      }
    } catch { /* ignore transient errors */ }
  }

  seqLog(`  AF [${filterName}] timed out (12 min) — proceeding`, "warn");
  if (Number.isFinite(finalPos)) saveAfResult(filterName, finalPos);
  return true; // full AF ran (timed out)
}

/** Legacy wrapper — kept for any direct call sites outside the filter loop */
async function stepAutofocus(target, availableFilters = []) {
  // Run AF without a specific filter (uses whatever filter is currently in place)
  const { equipmentInfo } = await getEquipmentInfo(target).catch(() => ({ equipmentInfo: null }));
  const currentFilter = equipmentInfo?.FilterWheel?.SelectedFilter?.Name ?? "__current__";
  await stepAutofocusForFilter(target, currentFilter, availableFilters);
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

/**
 * Fetch NINA's post-capture star statistics for the last image.
 * Returns stats object (stars, hfr, eccentricity, roundness, elongation)
 * or null if the endpoint is unavailable.
 */
async function fetchLastImageStats(target) {
  const toNum = (v) => {
    const n = Number(v);
    return Number.isFinite(n) ? n : null;
  };
  try {
    const res = await callNinaLong(
      target, "/v2/api/equipment/camera/capture/statistics", {}, 30000,
    );
    const r = res.body?.Response;
    if (!r) return null;
    const stars = toNum(r.Stars ?? r.DetectedStars);
    if (stars === null) return null;
    const hfr = toNum(r.HFR);
    const eccentricity = toNum(
      r.Eccentricity
      ?? r.MedianEccentricity
      ?? r.StarEccentricity
      ?? r.AverageEccentricity
    );
    const roundness = toNum(
      r.Roundness
      ?? r.MedianRoundness
      ?? r.StarRoundness
      ?? r.AverageRoundness
    );
    const directElong = toNum(
      r.Elongation
      ?? r.MedianElongation
      ?? r.StarElongation
      ?? r.AverageElongation
    );
    const derivedElong = (
      eccentricity !== null && eccentricity >= 0 && eccentricity < 0.995
    ) ? (1 / Math.sqrt(1 - (eccentricity * eccentricity))) : null;
    const elongation = directElong ?? derivedElong;

    return { stars, hfr, eccentricity, roundness, elongation };
  } catch {
    return null;
  }
}

function medianNumber(values) {
  const nums = values.filter(v => Number.isFinite(v)).sort((a, b) => a - b);
  if (!nums.length) return null;
  const mid = Math.floor(nums.length / 2);
  return nums.length % 2 === 0
    ? (nums[mid - 1] + nums[mid]) / 2
    : nums[mid];
}

/**
 * Issue a filter change and verify the wheel physically reached the requested
 * position by polling NINA until the reported slot matches.
 *
 * Strategy:
 *   - Poll every POLL_MS for up to MOVE_TIMEOUT_MS (covers worst-case long
 *     moves across the full wheel, empirically up to ~10 s).
 *   - If the position is confirmed, return immediately (fast moves ~1 s finish
 *     quickly without waiting for the full timeout).
 *   - If still wrong after the timeout, retry the command once then poll again.
 *   - Throw if still wrong after MAX_TRIES, aborting the sequence rather than
 *     imaging through the wrong filter.
 */
async function changeFilterVerified(target, filter) {
  const POLL_MS         = 500;   // how often to query position while waiting
  const MOVE_TIMEOUT_MS = 15000; // max time to wait for wheel to reach position
  const MAX_TRIES       = 2;

  for (let attempt = 1; attempt <= MAX_TRIES; attempt++) {
    seqLog(`Changing to filter: ${filter.Name} (ID ${filter.Id})${attempt > 1 ? ` (retry ${attempt})` : ""}...`);
    const changeResult = await callNinaLong(
      target, "/v2/api/equipment/filterwheel/change-filter", { filterId: filter.Id }, 30000,
    );
    if (!isNinaApiSuccess(changeResult)) {
      throw new Error(`Filter change to "${filter.Name}" (ID ${filter.Id}) failed`);
    }

    // Poll until the wheel reports the target position or we time out
    const deadline = Date.now() + MOVE_TIMEOUT_MS;
    let actual = null;
    let confirmed = false;
    while (Date.now() < deadline) {
      await new Promise(r => setTimeout(r, POLL_MS));
      try {
        const { equipmentInfo: fwPost } = await getEquipmentInfo(target);
        actual = fwPost?.FilterWheel?.SelectedFilter;
        if (actual && Number(actual.Id) === Number(filter.Id)) {
          confirmed = true;
          break;
        }
      } catch { /* transient query error — keep polling */ }
    }

    if (confirmed) {
      const elapsed = MOVE_TIMEOUT_MS - Math.max(0, deadline - Date.now());
      seqLog(`Filter → ${filter.Name} ✓  (confirmed slot ${actual.Id} in ~${Math.round(elapsed / 100) * 100} ms)`);
      return;
    }

    const actualName = actual ? `${actual.Name} (ID ${actual.Id})` : "unknown";
    seqLog(`Filter wheel did not reach ${filter.Name} (ID ${filter.Id}) within ${MOVE_TIMEOUT_MS / 1000}s — got ${actualName}`, "warn");

    if (attempt === MAX_TRIES) {
      throw new Error(
        `Filter wheel stuck: requested ${filter.Name} (ID ${filter.Id}) but wheel reports ${actualName} after ${MAX_TRIES} attempts`
      );
    }
    seqLog("Retrying filter change...", "warn");
  }
}

async function stepCaptureFilter(target, targetName, filterName, count, duration, gain, filters, safetyLimits = {}, targetRaDeg = null, targetDecDeg = null) {
  checkAbort();
  const filter = findFilterByName(filters, filterName);
  if (!filter) {
    const available = filters.map(f => f.Name).join(", ");
    seqLog(`Filter "${filterName}" not found (available: ${available || "none"}) — skipping`, "warn");
    return;
  }

  // Ensure filterwheel is connected before changing
  const { equipmentInfo: fwPre } = await getEquipmentInfo(target);
  if (!fwPre?.FilterWheel?.Connected) {
    seqLog("Filterwheel not connected — reconnecting...");
    const cr = await callNinaLong(target, "/v2/api/equipment/filterwheel/connect", {}, 15000);
    if (!isNinaApiSuccess(cr)) throw new Error("Filterwheel failed to connect");
    seqLog("Filterwheel connected ✓");
  }

  await changeFilterVerified(target, filter);
  checkAbort();

  let saved = 0;
  const recentGoodHfr = [];
  for (let i = 1; i <= count; i++) {
    checkAbort();

    // Mount safety check before every frame (tracking, zenith limit, meridian gap)
    await checkMountSafetyForFrame(target, safetyLimits);

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
      await new Promise(r => setTimeout(r, 20000));
    }

    // ── Step A: Star count — fast cloud gate ─────────────────────────────
    // NINA already ran star detection; we just fetch the result.
    // If too few stars → clouds → discard frame, wait, retake.
    // Only proceed to plate-solve drift check when stars look OK.
    const MIN_STARS          = safetyLimits.minStars ?? 10;
    const CLOUD_WAIT_MS      = 90_000;          // used by re-centering retries below
    const CLOUD_PROBE_EXP_S  = 5;               // probe frame for cloud recovery
    const CLOUD_PROBE_GAP_MS = 20_000;          // pause between probes
    const cloudRecoveryWaitMinRaw = Number(safetyLimits.maxCloudRecoveryWaitMin ?? 12);
    const cloudRecoveryWaitMin = Number.isFinite(cloudRecoveryWaitMinRaw)
      ? Math.min(120, Math.max(1, cloudRecoveryWaitMinRaw))
      : 12;
    const CLOUD_PROBE_MAX_MS = cloudRecoveryWaitMin * 60_000; // max wait for stars to come back
    const MAX_RETAKES        = 3;               // max cloud retakes before giving up on this frame
    const TRAIL_ELONG_LIMIT  = 1.85;
    const TRAIL_ECC_LIMIT    = 0.72;
    const TRAIL_HFR_FACTOR   = 2.0;
    const MAX_TRAIL_RECOVERY = 2;

    /**
     * Retake the current frame.
     * Updates `beforeFiles` in the outer scope and returns the new stats,
     * or null if the retake itself failed.
     */
    const doRetake = async () => {
      await checkMountSafetyForFrame(target, safetyLimits);
      const retakeBefore = snapshotNasFiles();
      try {
        const r2 = await callNinaLong(target,
          "/v2/api/equipment/camera/capture",
          { duration, gain, save: true, targetName, filter: filterName },
          30000,
        );
        if (!isNinaApiSuccess(r2)) { seqLog(`Retake failed — skipping`, "warn"); return null; }
      } catch (e) { seqLog(`Retake network error: ${e.message} — skipping`, "warn"); return null; }
      await waitExposure(duration);
      checkAbort();
      try { await waitForCameraIdle(target, 120_000); } catch { await new Promise(r => setTimeout(r, 20_000)); }
      await patchFitsFilterHeaders(retakeBefore, filter.Name);
      return await fetchLastImageStats(target);
    };

    const waitForCloudRecovery = async () => {
      const deadline = Date.now() + CLOUD_PROBE_MAX_MS;
      let probeNo = 0;
      while (Date.now() < deadline) {
        checkAbort();
        probeNo++;
        seqState.currentStep = `⛅ Cloud wait: 5s probe ${probeNo}`;
        await checkMountSafetyForFrame(target, safetyLimits);

        try {
          const probeResult = await callNinaLong(
            target,
            "/v2/api/equipment/camera/capture",
            { duration: CLOUD_PROBE_EXP_S, gain, save: false, filter: filterName },
            30000,
          );
          if (!isNinaApiSuccess(probeResult)) {
            seqLog(`Cloud probe rejected: ${probeResult.body?.Error || probeResult.status}`, "warn");
          }
        } catch (e) {
          seqLog(`Cloud probe network error: ${e.message} — retrying`, "warn");
        }

        await waitExposure(CLOUD_PROBE_EXP_S);
        checkAbort();
        try { await waitForCameraIdle(target, 120_000); } catch { await new Promise(r => setTimeout(r, 20_000)); }

        const probeStats = await fetchLastImageStats(target);
        if (probeStats?.stars != null) {
          if (probeStats.stars >= MIN_STARS) {
            seqLog(`☀ Cloud recovery: ${probeStats.stars} stars on 5s probe — resuming science frame`);
            return true;
          }
          seqLog(`⛅ Cloud probe ${probeNo}: only ${probeStats.stars} stars (need ${MIN_STARS})`, "warn");
        } else {
          seqLog(`⛅ Cloud probe ${probeNo}: no star stats returned`, "warn");
        }

        const gapDeadline = Date.now() + CLOUD_PROBE_GAP_MS;
        while (Date.now() < gapDeadline) {
          checkAbort();
          await new Promise(r => setTimeout(r, 5000));
        }
      }
      return false;
    };

    const isLikelyTrailing = (stats) => {
      if (!stats || stats.stars == null || stats.stars < MIN_STARS) return false;
      const baselineHfr = medianNumber(recentGoodHfr);
      const hfrTrail = (
        baselineHfr !== null
        && Number.isFinite(stats.hfr)
        && stats.hfr > baselineHfr * TRAIL_HFR_FACTOR
      );
      const eccTrail = Number.isFinite(stats.eccentricity) && stats.eccentricity >= TRAIL_ECC_LIMIT;
      const elongTrail = Number.isFinite(stats.elongation) && stats.elongation >= TRAIL_ELONG_LIMIT;
      return hfrTrail || eccTrail || elongTrail;
    };

    let imgStats  = await fetchLastImageStats(target);
    let retakes   = 0;
    let frameGood = true;   // set false only if we exhaust retakes

    if (imgStats !== null && MIN_STARS > 0) {
      while (imgStats !== null && imgStats.stars < MIN_STARS && retakes < MAX_RETAKES) {
        retakes++;
        seqLog(
          `⛅ ${filterName} ${i}/${count}: ${imgStats.stars} stars < min ${MIN_STARS} ` +
          `— cloud passage detected. Probing with 5s frames before retake ` +
          `(${retakes}/${MAX_RETAKES})...`,
          "warn",
        );
        seqState.currentStep = `⛅ Cloud — probing sky (${retakes}/${MAX_RETAKES})`;
        const recovered = await waitForCloudRecovery();
        if (!recovered) {
          seqLog(`⛅ Cloud did not clear within ${cloudRecoveryWaitMin} min`, "warn");
          break;
        }
        checkAbort();
        imgStats = await doRetake();
        if (imgStats === null) break; // retake itself failed; frame already patched above
      }

      if (imgStats !== null && imgStats.stars < MIN_STARS) {
        seqLog(`⛅ ${filterName} ${i}/${count}: still only ${imgStats.stars} stars after ${retakes} retakes — frame kept but flagged`, "warn");
        frameGood = false;
      } else if (imgStats !== null) {
        const hfrStr = imgStats.hfr > 0 ? `  HFR=${imgStats.hfr.toFixed(2)}` : "";
        seqLog(`  Stars: ${imgStats.stars}${hfrStr} ✓`);
      }
    }

    // ── Step A2: Trailing gate — restart guiding + retake ─────────────────
    let trailRecovery = 0;
    while (isLikelyTrailing(imgStats) && trailRecovery < MAX_TRAIL_RECOVERY) {
      trailRecovery++;
      const statsBits = [];
      if (Number.isFinite(imgStats?.hfr)) statsBits.push(`HFR=${imgStats.hfr.toFixed(2)}`);
      if (Number.isFinite(imgStats?.eccentricity)) statsBits.push(`Ecc=${imgStats.eccentricity.toFixed(2)}`);
      if (Number.isFinite(imgStats?.elongation)) statsBits.push(`Elong=${imgStats.elongation.toFixed(2)}`);
      seqLog(
        `⚠ ${filterName} ${i}/${count}: likely trailing (${statsBits.join(", ") || "shape anomaly"}) — restarting guiding + retake (${trailRecovery}/${MAX_TRAIL_RECOVERY})`,
        "warn",
      );
      const guidingOk = await stepStartGuiding(target);
      if (!guidingOk) {
        seqLog("  Guiding did not settle after restart attempt — retaking anyway", "warn");
      }
      imgStats = await doRetake();
      if (imgStats === null) break;
    }
    if (isLikelyTrailing(imgStats)) {
      seqLog(`⚠ ${filterName} ${i}/${count}: frame still looks trailed after guiding recovery`, "warn");
      frameGood = false;
    }
    if (imgStats !== null && imgStats.stars >= MIN_STARS && Number.isFinite(imgStats.hfr) && imgStats.hfr > 0) {
      recentGoodHfr.push(imgStats.hfr);
      if (recentGoodHfr.length > 8) recentGoodHfr.shift();
    }

    // ── Step B: Drift check — plate-solve only when stars are present ─────
    // Skipped when: frameCheckEnabled is false, stats say too few stars, or stats unavailable.
    const doFrameCheck   = safetyLimits.frameCheckEnabled === true;
    const driftThreshold = (safetyLimits.frameCheckThresholdArcmin ?? 5) * 60; // arcmin → arcsec
    const starsOkForSolve = imgStats === null           // stats not available → attempt anyway
                          || imgStats.stars >= MIN_STARS;

    if (doFrameCheck && starsOkForSolve) {
      seqState.currentStep = `${filterName}: drift check ${i}/${count}`;
      try {
        // Use NINA's solve-and-sync via the centering endpoint with waitForResult
        const solveResult = await callNinaLong(
          target, "/v2/api/equipment/camera/capture",
          { duration: safetyLimits.solveExp ?? 5, gain, solve: true, save: false },
          60_000,
        );
        const ps = solveResult.body?.Response?.PlateSolveResult ?? solveResult.body?.PlateSolveResult;
        if (ps?.Success === true || ps?.success === true) {
          // Compare solved coords to target coords to compute drift
          const solvedRaDeg  = (ps.Coordinates?.RA  ?? ps.RA  ?? null);
          const solvedDecDeg = (ps.Coordinates?.Dec ?? ps.Dec ?? null);
          if (solvedRaDeg !== null && solvedDecDeg !== null) {
            const driftArcsec = _angularSepArcsec(solvedRaDeg, solvedDecDeg, targetRaDeg, targetDecDeg);
            if (driftArcsec > driftThreshold) {
              seqLog(
                `⚠ Drift: ${Math.round(driftArcsec)}" > ${driftThreshold}" — re-centering...`,
                "warn",
              );
              await stepSlewWithCloudRetry(target, targetRaDeg, targetDecDeg, targetName,
                { center: true, cloudWaitMs: CLOUD_WAIT_MS, maxRetries: 2 });
            } else {
              seqLog(`  Drift: ${Math.round(driftArcsec)}" ✓`);
            }
          }
        } else if (!starsOkForSolve || imgStats?.stars > 0) {
          seqLog(`  Drift solve returned no result (${ps?.Reason ?? "unknown"}) — continuing`, "warn");
        }
      } catch (e) {
        seqLog(`  Drift check error: ${e.message} — continuing`, "warn");
      }
    } else if (doFrameCheck && !starsOkForSolve) {
      seqLog(`  Drift check skipped — too few stars for a reliable solve`);
    }

    // Patch any new FITS files with the correct FILTER header
    await patchFitsFilterHeaders(beforeFiles, filter.Name);

    saved++;
    seqLog(`${filterName} ${i}/${count} saved ✓${frameGood ? "" : " (⛅ flagged)"}`);
  }
  seqLog(`${filterName} complete (${saved}/${count}) ✓`);
}

// ── Main sequence runner ──────────────────────────────────────────────────────

async function runSequence(ninaConfig, seqConfig) {
  const {
    duration, gain, count, filters: filterNames, solveEnabled,
    minAlt = 20, meridianGap = 10, zenithLimit = 70, minStars = 10, maxCloudRecoveryWaitMin = 12,
    frameCheckEnabled = false, frameCheckThresholdArcmin = 5, solveExp = 5,
  } = seqConfig;
  function liveSafetyLimits() {
    return {
      minAlt:           Number(getSetting("minAlt",          minAlt)),
      meridianGap:      Number(getSetting("meridianGap",     meridianGap)),
      zenithLimit:      Number(getSetting("zenithLimit",     zenithLimit)),
      minStars:         Number(getSetting("minStars",        minStars)),
      maxCloudRecoveryWaitMin: Number(getSetting("maxCloudRecoveryWaitMin", maxCloudRecoveryWaitMin)),
      frameCheckEnabled:        getSetting("frameCheckEnabled", String(frameCheckEnabled)) === "true",
      frameCheckThresholdArcmin: Number(getSetting("frameCheckThresholdArcmin", frameCheckThresholdArcmin)),
      solveExp:         Number(getSetting("solveExp",        solveExp)),
    };
  }
  const safetyLimits = new Proxy({}, {
    get: (_, prop) => liveSafetyLimits()[prop],
    ownKeys: () => Object.keys(liveSafetyLimits()),
    getOwnPropertyDescriptor: (_, prop) => {
      const v = liveSafetyLimits()[prop];
      return v !== undefined ? { configurable: true, enumerable: true, value: v } : undefined;
    },
  });

  // Maximum time to wait for any target to enter a valid sky zone before giving up
  const MAX_WAIT_NO_TARGET_MS = 45 * 60 * 1000; // 45 min
  const WAIT_POLL_INTERVAL_MS =  5 * 60 * 1000; //  5 min between re-checks

  seqState.running    = true;
  seqState.aborted    = false;
  seqState.error      = null;
  seqState.manualMode = seqConfig.manualMode === true;
  seqState.ninaConfig = ninaConfig;
  // Reset stale daytime badSince so morning-cutoff timer starts fresh for tonight
  watchdog.badSince   = null;

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
    await assertRoofOpen("home");
    try {
      const { equipmentInfo: preInfo } = await getEquipmentInfo(ninaConfig);
      const mountPre = preInfo?.Mount;

      if (mountPre?.AtHome === true && mountPre?.AtPark !== true) {
        seqLog("Mount already at home ✓");
        return;
      }

      if (mountPre?.AtPark === true) {
        seqLog(`  Mount parked — unparking before FindHome...`);
        await callNinaLong(ninaConfig, "/v2/api/equipment/mount/unpark", {}, 30000);
        const deadline2 = Date.now() + 20000;
        while (Date.now() < deadline2) {
          await new Promise(r => setTimeout(r, 1500));
          const { equipmentInfo: ei2 } = await getEquipmentInfo(ninaConfig);
          if (ei2?.Mount?.AtPark === false) break;
        }
      }

      seqLog("  Sending FindHome command...");
      const result = await callNinaLong(ninaConfig, "/v2/api/equipment/mount/home", {}, 15000);
      const resp = result.body?.Response ?? "";
      seqLog(`  Home command accepted: ${resp}`);

      if (typeof resp === "string" && /already.?hom/i.test(resp)) {
        seqLog("Mount homed ✓ (already at home position)");
        return;
      }

      const HOME_TIMEOUT = 5 * 60 * 1000;
      const POLL_INTERVAL = 3000;
      const deadline = Date.now() + HOME_TIMEOUT;

      while (Date.now() < deadline) {
        try { checkAbort(); } catch { seqLog("  Home interrupted by abort — stopping", "warn"); return; }
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

    // ── Step 0a: Wait for watchdog to clear before moving anything ────────────
    if (watchdog.enabled && (watchdog.state === "bad" || watchdog.state === "recovering")) {
      seqLog("⛅ Watchdog: conditions bad at sequence start — waiting before moving mount", "warn");
      await waitForWatchdogClear();
    }

    // ── Step 0b: Check roof, unpark, go home ─────────────────────────────────
    try {
      // Verify OCS roof is open before moving anything
      const { ok: ocsOk, fields } = await ocsStatus(DEFAULT_OCS_HOST).catch(() => ({ ok: false, fields: {} }));
      const roofStatus = ocsClean(fields?.roof_sta ?? "").toLowerCase();
      const isClosed   = ocsOk && roofStatus && !roofStatus.includes("open");
      if (isClosed) {
        seqLog(`✗ Roof appears closed (OCS: "${ocsClean(fields.roof_sta)}") — aborting sequence for safety`, "error");
        return;
      }
      if (!ocsOk) {
        seqLog("OCS unreachable — cannot verify roof status, proceeding with caution", "warn");
      } else {
        seqLog(`Roof status: ${ocsClean(fields.roof_sta)} ✓`);
      }

      // Unpark + go home
      seqLog("── Mount pre-flight ─────────────────────────────────");
      await stepUnparkMount(ninaConfig);
      await safeHome();
      seqLog("Mount at home position ✓");
    } catch (err) {
      if (err.message === "__ABORTED__") throw err;
      seqLog(`✗ Mount pre-flight failed: ${err.message} — aborting sequence`, "error");
      return;
    }

    // ── Step 0b: Open dust cover ──────────────────────────────────────────────
    try {
      await stepOpenCover(ninaConfig);
    } catch (err) {
      if (err.message === "__ABORTED__") throw err;
      seqLog(`✗ Dust cover failed to open: ${err.message} — aborting sequence`, "error");
      return;
    }

    // ── Step 1: Cool camera ──────────────────────────────────────────────────
    try {
      await stepCoolCamera(ninaConfig, -5);
    } catch (err) {
      if (err.message === "__ABORTED__") throw err;
      seqLog(`✗ Camera cooling failed: ${err.message} — aborting sequence`, "error");
      await stepCloseCover(ninaConfig, { strict: false });
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

    // ── Step 1c: Ensure filterwheel is connected, then clear dark filter ──────────
    try {
      // Connect filterwheel if not already connected
      const { equipmentInfo: fwCheck } = await getEquipmentInfo(ninaConfig);
      if (!fwCheck?.FilterWheel?.Connected) {
        seqLog("Filterwheel not connected — connecting...");
        const cr = await callNinaLong(ninaConfig, "/v2/api/equipment/filterwheel/connect", {}, 15000);
        if (isNinaApiSuccess(cr)) seqLog("Filterwheel connected ✓");
        else seqLog("Filterwheel connect returned non-success — continuing", "warn");
      }
      const { equipmentInfo: fwInfo } = await getEquipmentInfo(ninaConfig);
      const filters = fwInfo?.FilterWheel?.AvailableFilters || [];
      const gFilter = findFilterByName(filters, "G");
      if (gFilter) {
        seqLog(`Switching to G filter (ID ${gFilter.Id}) to clear any dark filter...`);
        try {
          await changeFilterVerified(ninaConfig, gFilter);
        } catch (fwErr) {
          seqLog(`Startup filter switch warning: ${fwErr.message} — continuing anyway`, "warn");
        }
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

    // ── Dynamic target arbiter ────────────────────────────────────────────────
    //
    // Instead of a fixed order, at every decision point the sequencer picks the
    // *best available* target: highest altitude inside all safety limits.
    // Tie-break: most urgent (closest to leaving the valid window).
    // If no target is currently valid, we wait WAIT_POLL_INTERVAL_MS and retry.
    // The sequence ends when all targets are done or the wait budget runs out.

    const remaining = new Set(seqState.queue.filter(t => !t.done));
    const total = remaining.size;
    let done = 0;
    let totalWaitedMs = 0;

    // waitForWatchdogClear is defined at top level (above runSequence) so it
    // can also be called from stepSlewWithCloudRetry and checkMountSafetyForFrame.

    // Helper: get observer lat/lon from mount via NINA (cached per iteration)
    async function getObserverLatLon() {
      try {
        const { equipmentInfo } = await getEquipmentInfo(ninaConfig);
        const m = equipmentInfo?.Mount;
        if (m?.SiteLatitude != null && m?.SiteLongitude != null &&
            !(m.SiteLatitude === 0 && m.SiteLongitude === 0)) {
          return { lat: m.SiteLatitude, lon: m.SiteLongitude };
        }
      } catch { /* fall through */ }
      return null;
    }

    while (remaining.size > 0) {
      checkAbort();

      // ── Watchdog gate: suspend if sky is bad ──────────────────────────────
      await waitForWatchdogClear();
      checkAbort();

      // ── Pick best target ──────────────────────────────────────────────────
      const obs = await getObserverLatLon();
      let chosen = null;

      if (obs) {
        const result = pickBestTarget([...remaining], obs.lat, obs.lon, safetyLimits);
        if (result) {
          chosen = result.target;
          // Log the full ranking so the user can see why this target was chosen
          const ranking = result.allScores
            .map((s, idx) => `${idx === 0 ? "→" : " "}${s.name} alt${s.alt}° HA${s.ha}° [${s.minsLeft}]`)
            .join("  ");
          seqLog(`  Arbiter picks: ${ranking}`);
        }
      } else {
        // No observer location — fall back to queue order
        seqLog("  Arbiter: no observer location from mount — falling back to queue order", "warn");
        chosen = [...remaining][0];
      }

      // ── No valid target right now — wait and retry ────────────────────────
      if (!chosen) {
        if (totalWaitedMs >= MAX_WAIT_NO_TARGET_MS) {
          seqLog(`No valid target for ${MAX_WAIT_NO_TARGET_MS / 60000} min — ending sequence`, "warn");
          break;
        }
        seqLog(
          `⏳ No target in valid sky zone — waiting ${WAIT_POLL_INTERVAL_MS / 60000} min ` +
          `(${Math.round((MAX_WAIT_NO_TARGET_MS - totalWaitedMs) / 60000)} min budget left)...`,
          "warn",
        );
        await new Promise(r => setTimeout(r, WAIT_POLL_INTERVAL_MS));
        totalWaitedMs += WAIT_POLL_INTERVAL_MS;
        continue;
      }

      // Reset wait budget whenever we find a valid target
      totalWaitedMs = 0;

      const t = chosen;
      done++;
      seqState.currentTargetIdx = seqState.queue.indexOf(t);
      seqState.currentTarget = `${t.name} (${done}/${total})`;
      seqLog(`--- Target ${done}/${total}: ${t.name} ---`);

      // Declare per-target override vars here so both the try AND catch/retry
      // blocks can access them (const inside try is not visible in catch).
      let targetCount      = count;
      let targetDuration   = duration;
      let targetFilters    = verifiedFilters;
      let completedFilters = new Set();   // tracks filters done before any retry

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

        // Step 3 intentionally removed — AF now runs per-filter after slew (Step 6)

        // ── Step 4: Horizon check + Slew (+ cloud-retry plate-solve centering) ─
        await checkTargetAboveHorizon(ninaConfig, t.raDeg, t.decDeg, t.name);
        const useCenter = solveEnabled !== false;
        await stepSlewWithCloudRetry(ninaConfig, t.raDeg, t.decDeg, t.name, {
          center: useCenter, cloudWaitMs: 2 * 60 * 1000, maxRetries: 4,
        });
        if (!useCenter) seqLog("Plate solve centering disabled — skipping");
        const slewDoneMsg = useCenter
          ? `Slew + centering complete — proceed to start guiding?`
          : `Slew complete — proceed to start guiding?`;
        await waitForConfirmation(slewDoneMsg);
        checkAbort();

        // ── Step 5: Guiding ──────────────────────────────────────────────────
        const guidingOk = await stepStartGuiding(ninaConfig);
        if (!guidingOk) seqLog("⚠ Guiding failed to settle — continuing with caution", "warn");
        await waitForConfirmation(guidingOk ? `Guiding active — proceed to captures?` : `Guiding did not settle — continue anyway?`);
        checkAbort();

        // ── Step 6: Capture each filter ──────────────────────────────────────
        // Apply per-target overrides (vars declared above the try so catch can see them)
        targetCount    = (t.count    != null && t.count    > 0) ? t.count    : count;
        targetDuration = (t.duration != null && t.duration > 0) ? t.duration : duration;
        // Per-target filter override: intersect with verifiedFilters so we never use an absent filter
        targetFilters  = (Array.isArray(t.filters) && t.filters.length > 0)
          ? t.filters.filter(f => verifiedFilters.includes(f))
          : verifiedFilters;
        if (t.count    != null) seqLog(`  ↳ per-target frames  : ${targetCount} (global: ${count})`);
        if (t.duration != null) seqLog(`  ↳ per-target duration: ${targetDuration}s (global: ${duration}s)`);
        if (t.filters  != null) seqLog(`  ↳ per-target filters : ${targetFilters.join(", ")} (global: ${verifiedFilters.join(", ")})`);
        completedFilters = new Set();   // reset for this target (declared above try)
        for (const filterName of targetFilters) {
          checkAbort();
          // AF per filter: changes filter wheel then runs AF or uses 2h cache
          const fullAfRan = await stepAutofocusForFilter(ninaConfig, filterName, filters);
          // If a full AF ran, PHD2 guiding may have been disturbed — restart it
          if (fullAfRan) {
            seqLog(`AF [${filterName}] complete — restarting guiding before captures`);
            const guidingOkAf = await stepStartGuiding(ninaConfig);
            if (!guidingOkAf) seqLog("⚠ Guiding did not settle after AF — continuing with caution", "warn");
          }
          await waitForConfirmation(`Autofocus [${filterName}] done — proceed to captures?`);
          checkAbort();
          await assertCameraIdle(ninaConfig, `before ${filterName} captures`);
          // stepCaptureFilter will confirm/re-set the filter (no-op since already set by AF step)
          await stepCaptureFilter(ninaConfig, t.name, filterName, targetCount, targetDuration, gain, filters, safetyLimits, t.raDeg, t.decDeg);
          completedFilters.add(filterName);
          await waitForConfirmation(`${filterName} captures complete — proceed to next filter?`);
          checkAbort();
        }

        t.done = true;
        t.paused = false;
        remaining.delete(t);
        saveQueue();
        seqLog(`${t.name} — all captures complete ✓`);
        await waitForConfirmation(remaining.size > 0 ? `${t.name} done — proceed to next target?` : `${t.name} done — home mount?`);
        checkAbort();

      } catch (err) {
        if (err.message === "__ABORTED__") throw err;

        // ── ToO alert interrupt ───────────────────────────────────────────
        if (err.code === "__ALERT_INTERRUPT__") {
          const alertData = err.alert;
          seqLog(`🚨 ToO INTERRUPT: ${alertData.broker} ${alertData.trigger_id} — pausing ${t.name}`, "warn");

          // Save resume state
          seqState.tooResumeState = {
            target: t,
            completedFilters: [...completedFilters],
            targetFilters,
            targetCount,
            targetDuration,
          };

          // Abort the current exposure
          try {
            await callNinaLong(ninaConfig, "/v2/api/equipment/camera/abort", {}, 8000);
            seqLog("  Camera exposure aborted");
          } catch { /* best effort */ }

          // Run the ToO sub-sequence
          try {
            await runTooSequence(ninaConfig, alertData, seqConfig);
          } catch (tooErr) {
            seqLog(`  ToO sub-sequence error: ${tooErr.message}`, "error");
          }

          // Clear interrupt flag
          seqState.tooInterrupt = null;
          seqState.tooResumeState = null;

          // Resume: re-slew to the paused target and continue with remaining filters
          seqLog(`↩ Resuming ${t.name} after ToO...`);
          try {
            await stepUnparkMount(ninaConfig);
            await checkTargetAboveHorizon(ninaConfig, t.raDeg, t.decDeg, t.name);
            const useCenter = solveEnabled !== false;
            await stepSlewWithCloudRetry(ninaConfig, t.raDeg, t.decDeg, t.name, {
              center: useCenter, cloudWaitMs: 2 * 60 * 1000, maxRetries: 4,
            });
            const guidingOkResume = await stepStartGuiding(ninaConfig);
            if (!guidingOkResume) seqLog("⚠ Guiding failed to settle on resume — continuing with caution", "warn");

            const { equipmentInfo: reInfo } = await getEquipmentInfo(ninaConfig);
            const reFilters = reInfo?.FilterWheel?.AvailableFilters || [];
            const remainingFilters = targetFilters.filter(f => !completedFilters.has(f));
            if (remainingFilters.length > 0) {
              seqLog(`  Resuming filters: [${remainingFilters.join(", ")}] (done: [${[...completedFilters].join(", ")}])`);
              for (const filterName of remainingFilters) {
                checkAbort();
                const fullAfRanResume = await stepAutofocusForFilter(ninaConfig, filterName, reFilters);
                if (fullAfRanResume) {
                  seqLog(`AF [${filterName}] complete — restarting guiding before captures (resume)`);
                  const gOk = await stepStartGuiding(ninaConfig);
                  if (!gOk) seqLog("⚠ Guiding did not settle after AF (resume) — continuing with caution", "warn");
                }
                checkAbort();
                await assertCameraIdle(ninaConfig, `before ${filterName} captures (resume)`);
                await stepCaptureFilter(ninaConfig, t.name, filterName, targetCount, targetDuration, gain, reFilters, safetyLimits, t.raDeg, t.decDeg);
                completedFilters.add(filterName);
              }
            }
            t.done = true;
            remaining.delete(t);
            saveQueue();
            seqLog(`${t.name} — resumed and complete ✓`);
          } catch (resumeErr) {
            if (resumeErr.message === "__ABORTED__") throw resumeErr;
            seqLog(`  Resume of ${t.name} failed: ${resumeErr.message}`, "error");
            t.paused = true;
            await safeHome();
          }
          continue;
        }

        // ── Safety limit hit mid-capture ──────────────────────────────────
        if (err.code === "__LIMIT_REACHED__") {
          seqLog(`⚠ ${t.name} paused (limit): ${err.message}`, "warn");
          t.paused = true;
          seqLog(`  Homing mount — arbiter will re-evaluate target order`);
          await safeHome();
          continue;
        }

        // ── Plate solve failed (clouds) ────────────────────────────────────
        // The target stays in `remaining`. Home the mount and let the arbiter
        // try a different target; this one will be retried when the sky clears.
        if (err.code === "__CLOUD_SOLVE_FAIL__") {
          seqLog(`⛅ ${t.name} paused — plate solve failed after all retries (clouds). Moving to next target.`, "warn");
          t.paused = true;
          t.cloudPaused = true;
          await safeHome();
          continue;
        }

        // ── Transient network error ──────────────────────────────────────────
        const isTransient = /fetch failed|ECONNREFUSED|ETIMEDOUT|ENOTFOUND/i.test(err.message);
        const MAX_TARGET_RETRIES = 2;
        let targetSucceeded = false;

        if (isTransient) {
          for (let attempt = 1; attempt <= MAX_TARGET_RETRIES && !targetSucceeded; attempt++) {
            seqLog(`✗ ${t.name} transient error (${err.message}) — retry ${attempt}/${MAX_TARGET_RETRIES} in 30s...`, "warn");
            await new Promise(r => setTimeout(r, 30000));
            try {
              const { equipmentInfo: reInfo } = await getEquipmentInfo(ninaConfig);
              const reFilters = reInfo?.FilterWheel?.AvailableFilters || [];
              await stepUnparkMount(ninaConfig);
              await stepAutofocus(ninaConfig);
              await checkTargetAboveHorizon(ninaConfig, t.raDeg, t.decDeg, t.name);
              const useCenter = solveEnabled !== false;
              await stepSlewWithCloudRetry(ninaConfig, t.raDeg, t.decDeg, t.name, {
                center: useCenter, cloudWaitMs: 2 * 60 * 1000, maxRetries: 4,
              });
              const guidingOkRetry = await stepStartGuiding(ninaConfig);
              if (!guidingOkRetry) seqLog("⚠ Guiding failed to settle after retry slew — continuing with caution", "warn");
              const remainingFilters = targetFilters.filter(f => !completedFilters.has(f));
              if (completedFilters.size > 0) seqLog(`  ↳ retry: skipping already-done filters [${[...completedFilters].join(", ")}], running: [${remainingFilters.join(", ")}]`);
              for (const filterName of remainingFilters) {
                checkAbort();
                const fullAfRanRetry = await stepAutofocusForFilter(ninaConfig, filterName, reFilters);
                if (fullAfRanRetry) {
                  seqLog(`AF [${filterName}] complete — restarting guiding before captures (retry)`);
                  const gOk = await stepStartGuiding(ninaConfig);
                  if (!gOk) seqLog("⚠ Guiding did not settle after AF (retry) — continuing with caution", "warn");
                }
                checkAbort();
                await assertCameraIdle(ninaConfig, `before ${filterName} captures`);
                await stepCaptureFilter(ninaConfig, t.name, filterName, targetCount, targetDuration, gain, reFilters, safetyLimits, t.raDeg, t.decDeg);
                completedFilters.add(filterName);
              }
              t.done = true;
              remaining.delete(t);
              saveQueue();
              seqLog(`${t.name} — recovered and complete ✓`);
              targetSucceeded = true;
            } catch (retryErr) {
              if (retryErr.message === "__ABORTED__") throw retryErr;
              err = retryErr;
            }
          }
        }

        if (!targetSucceeded) {
          t.error = err.message;
          remaining.delete(t); // hard failure — remove from pool
          seqLog(`✗ ${t.name} failed: ${err.message}`, "error");
          seqLog(`  Homing mount before next target...`);
          await safeHome();
          seqLog(`  Continuing to next target`);
        }
      }
    }

    // ── Final: Home → close cover → park ─────────────────────────────────────
    seqLog("All targets complete. Homing mount...");
    await safeHome();
    seqLog("Mount at home ✓");
    try {
      const confirmed = await stepCloseCover(ninaConfig, { strict: false });
      if (!confirmed) seqLog("⚠ Verify dust cover is closed before leaving!", "warn");
    } catch (err) {
      seqLog(`⚠ Cover close error: ${err.message} — verify manually`, "warn");
    }
    try {
      await stepParkMount(ninaConfig);
    } catch (err) {
      seqLog(`⚠ Park failed: ${err.message} — verify mount is safe`, "warn");
    }
    seqLog("=== Sequence complete ===");

  } catch (err) {
    if (err.message === "__ABORTED__") {
      seqLog("=== Sequence aborted by user — stopping immediately ===", "warn");
      // User abort: no teardown, leave mount/cover as-is
    } else {
      seqLog(`=== Sequence error: ${err.message} ===`, "error");
      seqState.error = err.message;
      // Error teardown only (not user abort)
      try { await safeHome(); } catch { /* ignore */ }
      try { await stepCloseCover(ninaConfig, { strict: false }); } catch { /* ignore */ }
      try { await stepParkMount(ninaConfig); } catch { /* ignore */ }
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

// ── Alert Prep: "Ready Obs for Alert" ────────────────────────────────────────
// Pre-flight that gets the telescope into standby: roof check, unpark, home,
// open cover, cool camera, AF every filter, then park at home with tracking off.
// After this the telescope is ready to react in seconds: just slew and shoot.

async function runAlertPrep(ninaConfig, seqConfig) {
  if (seqState.alertPrepRunning) {
    seqLog("Alert prep already running — ignoring duplicate request", "warn");
    return;
  }
  seqState.alertPrepRunning = true;
  seqState.alertPrepStep = "Starting...";
  seqState.ninaConfig = ninaConfig;
  const { filters: filterNames = [], gain } = seqConfig;

  seqLog("═══ ALERT PREP: Ready Obs for Alert ═══");

  function prepStep(msg) {
    seqState.alertPrepStep = msg;
    seqLog(`  ${msg}`);
  }

  try {
    // ── 1. Roof check (manual — just warn) ─────────────────────────────────
    prepStep("Checking roof status...");
    try {
      const { ok: ocsOk, fields } = await ocsStatus(DEFAULT_OCS_HOST).catch(() => ({ ok: false, fields: {} }));
      const roofStatus = ocsClean(fields?.roof_sta ?? "").toLowerCase();
      const isClosed = ocsOk && roofStatus && !roofStatus.includes("open");
      if (isClosed) {
        seqLog(`✗ Roof appears CLOSED (OCS: "${ocsClean(fields.roof_sta)}") — aborting alert prep for safety`, "error");
        seqState.alertPrepStep = "ABORTED: Roof closed";
        return;
      }
      if (!ocsOk) {
        seqLog("  OCS unreachable — cannot verify roof. Proceeding with caution.", "warn");
      } else {
        seqLog(`  Roof: ${ocsClean(fields.roof_sta)} ✓`);
      }
    } catch (e) {
      seqLog(`  Roof check error: ${e.message} — proceeding with caution`, "warn");
    }
    if (seqState.aborted) return;

    // ── 2. Unpark ──────────────────────────────────────────────────────────
    prepStep("Unparking mount...");
    await stepUnparkMount(ninaConfig);
    if (seqState.aborted) return;

    // ── 3. Home ────────────────────────────────────────────────────────────
    prepStep("Homing mount...");
    // Inline safeHome logic (safeHome is scoped inside runSequence)
    try {
      const { equipmentInfo: preInfo } = await getEquipmentInfo(ninaConfig);
      const mountPre = preInfo?.Mount;
      if (mountPre?.AtPark === true || mountPre?.AtHome === true) {
        const reason = mountPre.AtPark ? "parked" : "already at home flag";
        seqLog(`    Mount ${reason} — unparking before FindHome...`);
        await callNinaLong(ninaConfig, "/v2/api/equipment/mount/unpark", {}, 30000);
        const deadline = Date.now() + 20000;
        while (Date.now() < deadline) {
          await new Promise(r => setTimeout(r, 1500));
          const { equipmentInfo: ei2 } = await getEquipmentInfo(ninaConfig);
          if (ei2?.Mount?.AtPark === false) break;
        }
      }
      seqLog("    Sending FindHome command...");
      const result = await callNinaLong(ninaConfig, "/v2/api/equipment/mount/home", {}, 15000);
      const resp = result.body?.Response ?? "";
      if (typeof resp === "string" && /already.?hom/i.test(resp)) {
        seqLog("    Mount homed ✓ (already at home position)");
      } else {
        const HOME_TIMEOUT = 5 * 60 * 1000;
        const deadline = Date.now() + HOME_TIMEOUT;
        while (Date.now() < deadline) {
          await new Promise(r => setTimeout(r, 3000));
          if (seqState.aborted) return;
          const { equipmentInfo } = await getEquipmentInfo(ninaConfig);
          const m = equipmentInfo?.Mount;
          seqState.alertPrepStep = `Homing: slewing=${m?.Slewing} atHome=${m?.AtHome}`;
          if (m?.AtHome === true && m?.Slewing !== true) {
            seqLog("    Mount homed ✓");
            break;
          }
        }
      }
    } catch (e) {
      seqLog(`    Home failed: ${e.message}`, "warn");
    }
    if (seqState.aborted) return;

    // ── 4. Open dust cover ─────────────────────────────────────────────────
    prepStep("Opening dust cover...");
    try {
      await stepOpenCover(ninaConfig);
    } catch (e) {
      seqLog(`    Cover open failed: ${e.message} — continuing`, "warn");
    }
    if (seqState.aborted) return;

    // ── 5. Cool camera ─────────────────────────────────────────────────────
    prepStep("Cooling camera...");
    try {
      await stepCoolCamera(ninaConfig, -5);
    } catch (e) {
      seqLog(`    Camera cooling failed: ${e.message} — continuing`, "warn");
    }
    if (seqState.aborted) return;

    // ── 6. AF on each filter (results are cached for 2 hours) ──────────────
    const { equipmentInfo: fwInfo } = await getEquipmentInfo(ninaConfig);
    const availableFilters = fwInfo?.FilterWheel?.AvailableFilters || [];
    const availableNames = availableFilters.map(f => f.Name).filter(Boolean);
    seqLog(`  Filterwheel: ${availableNames.join(", ") || "(none)"}`);

    // Need to enable tracking + slew to a bright area for AF to work
    prepStep("Enabling tracking for autofocus...");
    const { equipmentInfo: mInfo } = await getEquipmentInfo(ninaConfig);
    if (!mInfo?.Mount?.TrackingEnabled) {
      await callNinaLong(ninaConfig, "/v2/api/equipment/mount/tracking", { on: true }, 15000);
    }
    if (seqState.aborted) return;

    // Slew to a good AF target: current meridian, observatory latitude (near zenith)
    const site = resolveObserverSite();
    const jd = _toJD(new Date());
    const lstDeg = (_gmstDeg(jd) + site.lon + 360) % 360;
    const afRaDeg  = lstDeg;
    const afDecDeg = Math.min(80, Math.max(-80, site.lat));
    prepStep(`Slewing to AF position (RA=${(afRaDeg/15).toFixed(2)}h Dec=${afDecDeg.toFixed(1)}°)...`);
    await stepSlewToTarget(ninaConfig, afRaDeg, afDecDeg, "AlertPrep-AF", { center: false });
    await new Promise(r => setTimeout(r, 3000));
    if (seqState.aborted) return;

    const validFilters = filterNames.filter(f => findFilterByName(availableFilters, f));
    if (!validFilters.length && availableNames.length) {
      seqLog("  No configured filters match filterwheel — using all available", "warn");
      validFilters.push(...availableNames);
    }

    for (let i = 0; i < validFilters.length; i++) {
      if (seqState.aborted) return;
      const filterName = validFilters[i];
      prepStep(`Autofocus [${filterName}] (${i + 1}/${validFilters.length})...`);
      try {
        await stepAutofocusForFilter(ninaConfig, filterName, availableFilters);
        seqLog(`  AF [${filterName}] complete ✓`);
      } catch (e) {
        seqLog(`  AF [${filterName}] failed: ${e.message} — continuing`, "warn");
      }
    }
    if (seqState.aborted) return;

    // ── 7. Go home, disable tracking, wait ─────────────────────────────────
    prepStep("Returning to home position...");
    try {
      await callNinaLong(ninaConfig, "/v2/api/equipment/mount/home", {}, 15000);
      const deadline = Date.now() + 3 * 60 * 1000;
      while (Date.now() < deadline) {
        await new Promise(r => setTimeout(r, 3000));
        const { equipmentInfo } = await getEquipmentInfo(ninaConfig);
        if (equipmentInfo?.Mount?.AtHome === true && !equipmentInfo?.Mount?.Slewing) break;
      }
    } catch (e) {
      seqLog(`    Home failed: ${e.message}`, "warn");
    }

    prepStep("Disabling tracking...");
    try {
      await callNinaLong(ninaConfig, "/v2/api/equipment/mount/tracking", { on: false }, 15000);
      seqLog("  Tracking disabled ✓");
    } catch (e) {
      seqLog(`  Tracking disable failed: ${e.message}`, "warn");
    }

    // Auto-enable alert-ready mode
    seqState.alertReady = true;
    setSetting("alertReady", "true");

    seqLog("═══ ALERT PREP COMPLETE — Observatory is ready ═══");
    seqLog("  Cover: open  |  Camera: cooled  |  AF: cached  |  Mount: home, tracking off");
    seqLog("  Alert-Ready mode: ARMED — waiting for alerts...", "warn");
    seqState.alertPrepStep = "READY — Waiting for alerts";

  } catch (e) {
    if (e.message === "__ABORTED__") {
      seqLog("═══ ALERT PREP ABORTED ═══", "warn");
      seqState.alertPrepStep = "Aborted";
    } else {
      seqLog(`═══ ALERT PREP ERROR: ${e.message} ═══`, "error");
      seqState.alertPrepStep = `Error: ${e.message}`;
    }
  } finally {
    seqState.alertPrepRunning = false;
    seqState.aborted = false;
    seqState.ninaConfig = ninaConfig;
  }
}

// ── Target-of-Opportunity (ToO) sub-sequence ─────────────────────────────────
// Runs when an alert interrupt fires. Uses the same camera settings as the main
// sequence. Implements a mosaic strategy based on the error box vs FOV.
//
// Mosaic rules (see computeMosaicGrid):
//   err_deg ≤ FOV/2  → single center pointing (localization fits in one inscribed circle)
//   err_deg > FOV/2  → n×n grid with side = ceil(2×err/FOV), tile spacing = FOV
//                      center-out spiral; one filter cycle per tile per epoch
//
/** Center-out spiral visitation order for an n×n tile grid */
function orderTilesSpiral(tiles, side) {
  const cr = (side - 1) / 2;
  const cc = (side - 1) / 2;
  return [...tiles].sort((a, b) => {
    const drA = a.row - cr;
    const dcA = a.col - cc;
    const drB = b.row - cr;
    const dcB = b.col - cc;
    const da = drA * drA + dcA * dcA;
    const db = drB * drB + dcB * dcB;
    if (da !== db) return da - db;
    return Math.atan2(drA, dcA) - Math.atan2(drB, dcB);
  });
}

function computeMosaicGrid(raDeg, decDeg, errDeg, fovDeg) {
  const threshold = getTooMosaicThresholdDeg(fovDeg);

  if (!Number.isFinite(errDeg) || errDeg <= threshold) {
    return {
      grid: [{ ra: raDeg, dec: decDeg, label: "center" }],
      mode: "single",
      side: 1,
      requestedSide: 1,
      capped: false,
    };
  }

  const requestedSide = Math.max(1, Math.ceil((2 * errDeg) / fovDeg));
  const side = Math.min(requestedSide, TOO_MOSAIC_MAX_SIDE);
  const capped = requestedSide > side;
  const cosDec = Math.cos(decDeg * Math.PI / 180) || 1;
  const step = fovDeg;

  const tiles = [];
  for (let row = 0; row < side; row++) {
    for (let col = 0; col < side; col++) {
      const dra = ((col - (side - 1) / 2) * step) / cosDec;
      const ddec = (row - (side - 1) / 2) * step;
      tiles.push({
        ra: raDeg + dra,
        dec: decDeg + ddec,
        row,
        col,
        label: `T${row + 1}-${col + 1}`,
      });
    }
  }

  return {
    grid: orderTilesSpiral(tiles, side),
    mode: "spiral",
    side,
    requestedSide,
    capped,
  };
}

async function runTooSequence(ninaConfig, alertData, seqConfig) {
  seqState.tooRunning = true;
  const { duration, gain, count, filters: filterNames } = seqConfig;
  const { ra, dec, err_deg, broker, trigger_id, strategy } = alertData;
  const tooName = `ToO-${broker}-${trigger_id || "alert"}`;

  // ── Resolve strategy params (fallback to seqConfig if no strategy) ──────
  const strat = strategy || getStrategy(broker);
  const rapidExposure  = strat.exposure     || duration;
  const rapidGain      = strat.gain         ?? gain;
  const rapidCount     = Math.min(count, strat.rapid_count || 15);
  const filterCycle    = (strat.filter_cycle || "G,RP,BP").split(",").map(s => s.trim()).filter(Boolean);
  const doAf           = Boolean(strat.do_af);
  const doGuiding      = Boolean(strat.do_guiding);
  const doCenter       = Boolean(strat.do_center);

  seqLog(`🚨 ═══ ToO SUB-SEQUENCE: ${tooName} ═══`);
  seqLog(`  RA=${ra.toFixed(4)}° Dec=${dec.toFixed(4)}° err=${(err_deg || 0).toFixed(3)}° FOV=${getTooFovDeg()}°`);
  seqLog(`  Strategy: mode=${strat.mode} exp=${rapidExposure}s gain=${rapidGain} filters=[${filterCycle}] rapid=${rapidCount} AF=${doAf} guide=${doGuiding} center=${doCenter}`);

  try {
    const fovDeg = getTooFovDeg();
    const mosaicThreshold = getTooMosaicThresholdDeg(fovDeg);
    const mosaic = computeMosaicGrid(ra, dec, err_deg || 0, fovDeg);
    const grid = mosaic.grid;
    const useSpiral = mosaic.mode === "spiral";
    const centerSpanDeg = ((mosaic.side - 1) * fovDeg).toFixed(2);
    const outerExtentDeg = (mosaic.side * fovDeg).toFixed(2);
    let mosaicMsg = `  Mosaic: ${grid.length} tile(s)`;
    if (useSpiral) {
      mosaicMsg += ` ${mosaic.side}×${mosaic.side} spiral (err=${(err_deg || 0).toFixed(2)}°, FOV=${fovDeg}°, single if err≤${mosaicThreshold.toFixed(2)}°`;
      mosaicMsg += ` → ~${centerSpanDeg}° between outer centers, ~${outerExtentDeg}° full width)`;
      if (mosaic.capped) {
        mosaicMsg += ` [capped from ${mosaic.requestedSide}×${mosaic.requestedSide}]`;
      }
      mosaicMsg += ` — ${grid.map(g => g.label).join(" → ")}`;
    } else {
      mosaicMsg += " — center";
    }
    seqLog(mosaicMsg);

    const { equipmentInfo } = await getEquipmentInfo(ninaConfig);
    const availFilters = equipmentInfo?.FilterWheel?.AvailableFilters || [];

    // ── Phase 1 fast-path: first filter from cycle + focuser restore ────────
    const firstFilterName = filterCycle[0] || "G";
    const firstFilter = findFilterByName(availFilters, firstFilterName);
    if (firstFilter) {
      await changeFilterVerified(ninaConfig, firstFilter);
      if (doAf) {
        seqLog(`  Phase 1 AF: running autofocus for ${firstFilterName}`);
        await stepAutofocusForFilter(ninaConfig, firstFilterName, availFilters);
      } else {
        // Fix 1: prefer recent AF result (<2h); only fall back to temp model if AF is stale/absent
        const lastAf = getAfResultIfFresh(firstFilterName); // returns null if >2h old
        if (lastAf) {
          const ageMin = Math.round((Date.now() - lastAf.ts) / 60000);
          seqLog(`  Fast focus: using recent AF pos ${lastAf.pos} for ${firstFilterName} (${ageMin} min ago)`);
          await stepMoveFocuserAbsolute(ninaConfig, lastAf.pos);
        } else {
          // No recent AF — fall back to temperature model
          const foctemp   = await getFocuserTemp(ninaConfig);
          const predicted = predictFocuserPos(firstFilterName, foctemp);
          const staleAf   = getLatestAfResult(firstFilterName);

          if (predicted != null) {
            const drift = staleAf ? predicted - staleAf.pos : null;
            const driftStr = drift != null
              ? ` (last AF pos ${staleAf.pos} ${Math.round((Date.now()-staleAf.ts)/60000)}min ago, Δ${drift >= 0 ? "+" : ""}${drift} steps)`
              : " (no prior AF in DB)";
            seqLog(`  🌡 Fast focus (no recent AF): FOCTEMP=${foctemp?.toFixed(1)}°C → predicted pos ${predicted} [${firstFilterName}]${driftStr}`);
            if (drift != null && Math.abs(drift) > 100) {
              seqLog(`  ⚠ Drift ${Math.abs(drift)} steps > 100 from last AF — check camera position`, "warn");
            }
            await stepMoveFocuserAbsolute(ninaConfig, predicted);
          } else if (staleAf) {
            const ageMin = Math.round((Date.now() - staleAf.ts) / 60000);
            seqLog(`  Fast focus: FOCTEMP unavailable — using stale AF pos ${staleAf.pos} for ${firstFilterName} (${ageMin} min ago)`, "warn");
            await stepMoveFocuserAbsolute(ninaConfig, staleAf.pos);
          } else {
            seqLog(`  No AF data and no FOCTEMP for ${firstFilterName} — proceeding with current focuser position`, "warn");
          }
        }
      }
    } else {
      seqLog(`  ${firstFilterName} filter not found — using first available`, "warn");
    }
    const fastFilter = firstFilter ? firstFilterName : (filterCycle.find(f => findFilterByName(availFilters, f)) || filterNames[0] || null);
    if (!fastFilter) {
      seqLog("  No valid filters — aborting ToO", "error");
      return;
    }

    await stepUnparkMount(ninaConfig);

    // Fix 5: update status text to show ToO target while interrupt is running
    seqState.currentTarget = `🚨 ToO: ${tooName}`;

    // ── Phase 1: rapid filter-cycle imaging (time-critical) ─────────────────
    const validCycleFilters = filterCycle.filter(f => findFilterByName(availFilters, f));
    const cycleLen    = validCycleFilters.length || 1;
    // rapid_count = number of complete filter cycles (e.g. SVOM=1 → G+RP+BP once, Swift=15 → 15×G+RP+BP)
    const rapidCycles = Math.max(1, strat.rapid_count || 1);
    const rapidCount  = rapidCycles * cycleLen;   // total frames in Phase 1
    const totalEpochs = rapidCycles;              // one epoch per complete filter cycle
    seqLog(`  Phase 1: ${rapidCycles} × [${validCycleFilters.join(",")}] = ${rapidCount} frames → ${totalEpochs} epoch(s) (AF=${doAf}, guide=${doGuiding}, center=${doCenter})`);

    const captureRapidFrame = async (pt, epochNum, epochTag, cycleFilter, frameIdx) => {
      const epochName = `${tooName}-${pt.label}-${epochTag}`;
      const cFilter = findFilterByName(availFilters, cycleFilter);
      if (cFilter) await changeFilterVerified(ninaConfig, cFilter);
      seqState.currentStep = `ToO ${pt.label} ${epochTag} ${cycleFilter}: ${frameIdx}/${rapidCount} [RAPID]`;
      seqState.progress = { filter: cycleFilter, frame: frameIdx, frames: rapidCount, epoch: epochNum, totalEpochs };
      seqLog(`  ${pt.label} ${epochTag} ${cycleFilter} ${frameIdx}/${rapidCount}: ${rapidExposure}s`);
      await callNinaLong(ninaConfig, "/v2/api/equipment/camera/capture",
        { duration: rapidExposure, gain: rapidGain, save: true, targetName: epochName, filter: cycleFilter }, 30000);
      await waitExposure(rapidExposure);
      try { await waitForCameraIdle(ninaConfig, 120000); } catch { await new Promise(r => setTimeout(r, 20000)); }
    };

    if (useSpiral) {
      // Spiral: one filter cycle per tile per epoch (epochs interleaved across the mosaic)
      for (let ep = 1; ep <= rapidCycles; ep++) {
        if (seqState.aborted) break;
        const epochTag = `E${String(ep).padStart(2, "0")}`;
        seqLog(`  ── Epoch ${epochTag} — spiral ${grid.length} tiles × [${validCycleFilters.join(",")}] ──`);
        for (let ti = 0; ti < grid.length; ti++) {
          const pt = grid[ti];
          if (seqState.aborted) break;
          seqLog(`  → Slewing to ${tooName} [${pt.label}] (${ti + 1}/${grid.length})`);
          await stepSlewToTarget(ninaConfig, pt.ra, pt.dec, `${tooName}-${pt.label}`, { center: doCenter });
          if (!doCenter) await new Promise(r => setTimeout(r, 2000));
          if (doGuiding) {
            seqLog(`  Starting PHD2 guiding for ToO...`);
            const gOk = await stepStartGuiding(ninaConfig);
            if (!gOk) seqLog("  ⚠ Guiding did not settle for ToO — proceeding without guiding", "warn");
          }
          for (let fi = 0; fi < cycleLen; fi++) {
            if (seqState.aborted) break;
            const cycleFilter = validCycleFilters[fi] || fastFilter;
            const frameIdx = (ep - 1) * grid.length * cycleLen + ti * cycleLen + fi + 1;
            await captureRapidFrame(pt, ep, epochTag, cycleFilter, frameIdx);
          }
        }
      }
    } else {
      // Single pointing: all epochs at center before moving on
      for (const pt of grid) {
        if (seqState.aborted) break;
        seqLog(`  → Slewing to ${tooName} [${pt.label}]`);
        await stepSlewToTarget(ninaConfig, pt.ra, pt.dec, `${tooName}-${pt.label}`, { center: doCenter });
        if (!doCenter) await new Promise(r => setTimeout(r, 2000));
        if (doGuiding) {
          seqLog(`  Starting PHD2 guiding for ToO...`);
          const gOk = await stepStartGuiding(ninaConfig);
          if (!gOk) seqLog("  ⚠ Guiding did not settle for ToO — proceeding without guiding", "warn");
        }
        let currentEpoch = 0;
        for (let i = 1; i <= rapidCount; i++) {
          if (seqState.aborted) break;
          const cycleFilter = validCycleFilters[(i - 1) % cycleLen] || fastFilter;
          const epochNum = Math.ceil(i / cycleLen);
          const epochTag = `E${String(epochNum).padStart(2, "0")}`;
          if (epochNum !== currentEpoch) {
            seqLog(`  ── Epoch ${epochTag} (frames ${(epochNum - 1) * cycleLen + 1}–${Math.min(epochNum * cycleLen, rapidCount)}) ──`);
            currentEpoch = epochNum;
          }
          await captureRapidFrame(pt, epochNum, epochTag, cycleFilter, i);
        }
      }
    }

    // ── Phase 2: multi-filter deep follow-up (remaining cycles, with optional AF) ─
    const tooFilters = filterCycle.filter(f => findFilterByName(availFilters, f));
    // Fix 6: remaining cycles = total strategy cycles minus the 1 rapid cycle already done
    const remainingCount = Math.max(0, rapidCycles - 1);
    if (remainingCount > 0 && tooFilters.length > 0 && !seqState.aborted) {
      // Phase 2 epoch offset: continue numbering after Phase 1
      const p2EpochOffset = totalEpochs;
      seqLog(`  Phase 2: multi-filter follow-up — ${remainingCount} more epoch(s) × [${tooFilters.join(", ")}] (epochs E${String(p2EpochOffset + 1).padStart(2,"0")}+)`);

      const runPhase2AtTile = async (pt, ep, epochNum, epochTag, epochName, ti, tileCount) => {
        seqLog(`  → Slewing to ${tooName} [${pt.label}] (Phase 2${tileCount ? ` ${ti + 1}/${tileCount}` : ""})`);
        await stepSlewToTarget(ninaConfig, pt.ra, pt.dec, `${tooName}-${pt.label}`, { center: false });
        await new Promise(r => setTimeout(r, 2000));
        if (doGuiding) {
          const gOk2 = await stepStartGuiding(ninaConfig);
          if (!gOk2) seqLog("  ⚠ Guiding did not settle (Phase 2) — proceeding", "warn");
        }
        seqLog(`  ── Epoch ${epochTag} (Phase 2, pass ${ep}/${remainingCount}) ──`);
        for (const filterName of tooFilters) {
          if (seqState.aborted) break;
          const filter = findFilterByName(availFilters, filterName);
          if (!filter) continue;
          if (ep === 1 && doAf) {
            const fullAfRan = await stepAutofocusForFilter(ninaConfig, filterName, availFilters);
            if (fullAfRan && doGuiding) {
              seqLog(`  AF [${filterName}] complete — restarting guiding (Phase 2)`);
              await stepStartGuiding(ninaConfig);
            }
          }
          await changeFilterVerified(ninaConfig, filter);
          seqState.currentStep = `ToO ${pt.label} ${epochTag} ${filterName}: ${ep}/${remainingCount}`;
          seqState.progress = { filter: filterName, frame: ep, frames: remainingCount, epoch: epochNum };
          seqLog(`  ${pt.label} ${epochTag} ${filterName} ${ep}/${remainingCount}: ${duration}s`);
          await callNinaLong(ninaConfig, "/v2/api/equipment/camera/capture",
            { duration, gain, save: true, targetName: epochName, filter: filterName }, 30000);
          await waitExposure(duration);
          try { await waitForCameraIdle(ninaConfig, 120000); } catch { await new Promise(r => setTimeout(r, 20000)); }
        }
      };

      if (useSpiral) {
        for (let ep = 1; ep <= remainingCount; ep++) {
          if (seqState.aborted) break;
          const epochNum = p2EpochOffset + ep;
          const epochTag = `E${String(epochNum).padStart(2, "0")}`;
          for (let ti = 0; ti < grid.length; ti++) {
            if (seqState.aborted) break;
            const pt = grid[ti];
            const epochName = `${tooName}-${pt.label}-${epochTag}`;
            await runPhase2AtTile(pt, ep, epochNum, epochTag, epochName, ti, grid.length);
          }
        }
      } else {
        for (const pt of grid) {
          if (seqState.aborted) break;
          for (let ep = 1; ep <= remainingCount; ep++) {
            if (seqState.aborted) break;
            const epochNum = p2EpochOffset + ep;
            const epochTag = `E${String(epochNum).padStart(2, "0")}`;
            const epochName = `${tooName}-${pt.label}-${epochTag}`;
            await runPhase2AtTile(pt, ep, epochNum, epochTag, epochName, null, null);
          }
        }
      }
    }

    seqLog(`🚨 ═══ ToO ${tooName} COMPLETE ═══`);

    try {
      db.prepare("UPDATE alerts SET action = 'too-observed', action_reason = ? WHERE trigger_id = ?")
        .run(`${grid.length}-pointings`, String(trigger_id));
    } catch { /* non-fatal */ }

  } finally {
    seqState.tooRunning = false;
    seqState.tooInterrupt = null;
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
    alertReady: seqState.alertReady,
    tooRunning: seqState.tooRunning,
    alertPrepRunning: seqState.alertPrepRunning,
    alertPrepStep: seqState.alertPrepStep,
  });
});

// Return list of available nightly log files
app.get("/api/sequence/logs", (req, res) => {
  try {
    const files = fs.readdirSync(SEQ_LOG_DIR)
      .filter(f => f.startsWith("seq-") && f.endsWith(".log"))
      .sort().reverse();
    res.json(files);
  } catch { res.json([]); }
});

// Return content of a specific nightly log (default: tonight)
app.get("/api/sequence/logs/:date", (req, res) => {
  const date = req.params.date.replace(/[^0-9-]/g, ""); // sanitise
  const file = path.join(SEQ_LOG_DIR, `seq-${date}.log`);
  if (!fs.existsSync(file)) return res.status(404).json({ error: "Log not found" });
  res.setHeader("Content-Type", "text/plain; charset=utf-8");
  res.sendFile(file);
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

// PATCH a single queue item (e.g. per-target frame count override)
app.patch("/api/sequence/queue/:index", (req, res) => {
  const idx = Number(req.params.index);
  if (!Number.isInteger(idx) || idx < 0 || idx >= seqState.queue.length) {
    return res.status(400).json({ success: false, error: "Invalid queue index" });
  }
  const item = seqState.queue[idx];
  // Per-target overrides: count, duration, filters, waitTime
  if (req.body?.count !== undefined) {
    const c = Number(req.body.count);
    item.count = (Number.isFinite(c) && c > 0) ? c : null;
  }
  if (req.body?.duration !== undefined) {
    const d = Number(req.body.duration);
    item.duration = (Number.isFinite(d) && d > 0) ? d : null;
  }
  if (req.body?.filters !== undefined) {
    const f = req.body.filters;
    if (Array.isArray(f) && f.length > 0) {
      item.filters = f.map(s => String(s).trim()).filter(Boolean);
    } else if (typeof f === "string" && f.trim()) {
      item.filters = f.split(",").map(s => s.trim()).filter(Boolean);
    } else {
      item.filters = null;
    }
  }
  if (req.body?.waitTime !== undefined) {
    item.waitTime = String(req.body.waitTime);
  }
  saveQueue();
  return res.json({ success: true, item });
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

// ── Daily reset endpoint ──────────────────────────────────────────────────────
app.post("/api/sequence/daily-reset", (req, res) => {
  if (seqState.running) {
    return res.status(409).json({ success: false, error: "Cannot reset while sequence is running" });
  }
  const result = dailyReset({ manual: true });
  res.json({ success: true, ...result });
});

app.get("/api/sequence/daily-reset/info", (req, res) => {
  const hour = parseInt(getSetting("daily_reset_hour", "12"), 10);
  const alertTargets  = seqState.queue.filter(t => t.source === "alert").map(t => t.name);
  const manualTargets = seqState.queue.filter(t => t.source !== "alert").map(t => t.name);
  res.json({ success: true, resetHour: hour, alertTargets, manualTargets });
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
    solveEnabled:              req.body?.solveEnabled !== false,
    solveExp:                  Number(req.body?.solveExp)                  || 5,
    solveThreshold:            Number(req.body?.solveThreshold)            || 60,
    manualMode:                req.body?.manualMode === true,
    minAlt:                    Number(req.body?.minAlt)                    || 20,
    zenithLimit:               Number(req.body?.zenithLimit)               || 70,
    meridianGap:               Number(req.body?.meridianGap)               || 10,
    minStars:                  Number(req.body?.minStars)                  || 10,
    frameCheckEnabled:         req.body?.frameCheckEnabled === true,
    frameCheckThresholdArcmin: Number(req.body?.frameCheckThresholdArcmin) || 5,
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
  setSetting("lastAutofocusTime", "");
  return res.json({ success: true, message: "Autofocus timer reset — AF will run on next sequence" });
});

// Returns the temperature-predicted focus position for the current FOCTEMP,
// alongside the last stored AF result per filter. Useful for debugging the model.
app.get("/api/autofocus/predict", async (req, res) => {
  const ninaTarget = seqState.ninaConfig || { host: DEFAULT_NINA_HOST, port: DEFAULT_NINA_PORT, protocol: DEFAULT_NINA_PROTOCOL };
  const foctemp = await getFocuserTemp(ninaTarget);
  const predictions = {};
  for (const f of ["G", "BP", "RP", "ALL"]) {
    const reg = AF_TEMP_REGRESSION[f];
    const lastAf = getLatestAfResult(f === "ALL" ? "__none__" : f);
    predictions[f] = {
      slope: reg.slope,
      intercept: reg.intercept,
      predicted: foctemp != null ? predictFocuserPos(f, foctemp) : null,
      lastStoredPos: lastAf?.pos ?? null,
      lastStoredTs:  lastAf?.ts  ?? null,
      drift: (foctemp != null && lastAf) ? predictFocuserPos(f, foctemp) - lastAf.pos : null,
    };
  }
  res.json({ success: true, foctemp, predictions });
});

app.get("/api/autofocus/history", (req, res) => {
  const limit = Math.min(200, Math.max(1, parseInt(req.query.limit) || 50));
  const filter = req.query.filter;
  let rows;
  if (filter) {
    rows = db.prepare(
      `SELECT id, filter, position, ts FROM autofocus_log WHERE filter = ? ORDER BY ts DESC LIMIT ?`
    ).all(filter, limit);
  } else {
    rows = db.prepare(
      `SELECT id, filter, position, ts FROM autofocus_log ORDER BY ts DESC LIMIT ?`
    ).all(limit);
  }
  const latest = {};
  for (const r of rows) {
    if (!latest[r.filter]) latest[r.filter] = { position: r.position, ts: r.ts };
  }
  res.json({ success: true, latest, history: rows });
});

app.post("/api/sequence/restart", (req, res) => {
  // If a sequence is running, abort it first (same as clicking Abort)
  if (seqState.running) {
    seqState.aborted = true;
    if (_confirmReject) {
      _confirmReject(new Error("__ABORTED__"));
      _confirmResolve = null;
      _confirmReject  = null;
    }
    if (seqState.ninaConfig) {
      const cfg = seqState.ninaConfig;
      callNinaLong(cfg, "/v2/api/equipment/camera/abort", {}, 8000).catch(() => {});
      callNinaLong(cfg, "/v2/api/equipment/guider/stop",  {}, 8000).catch(() => {});
    }
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
        type:        (data.type || data.event_type || data.classification || null),
        discoveryDate: data.timestamp
          ? new Date(data.timestamp).toISOString().slice(0, 10) : null,
        constellation: data.constellation || null,
        hostName:    null,
        reporters:   null,
        astrocolibriId: data.astro_colibri_id || null,
        colibriUrl:  `https://astro-colibri.com/?trigger_id=${triggerId}`,
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

// ── Open-Meteo weather proxy ──────────────────────────────────────────────────
//
// GET /api/weather/openmeteo?lat=47.5&lon=2.3
//
// Fetches the next 48 h of hourly weather from Open-Meteo (free, no API key)
// and returns a normalised { points: [...] } array in the format expected by
// the night-plan weather widget.

app.get("/api/weather/openmeteo", async (req, res) => {
  const lat = parseFloat(req.query.lat);
  const lon = parseFloat(req.query.lon);
  if (isNaN(lat) || isNaN(lon)) {
    return res.status(400).json({ success: false, error: "lat and lon are required" });
  }

  const params = new URLSearchParams({
    latitude:           lat,
    longitude:          lon,
    hourly:             "cloudcover,precipitation,relativehumidity_2m,windspeed_10m,winddirection_10m,temperature_2m",
    wind_speed_unit:    "ms",
    precipitation_unit: "mm",
    timezone:           "UTC",
    forecast_days:      2,
  });

  try {
    const r = await fetch(`https://api.open-meteo.com/v1/forecast?${params}`, {
      signal: AbortSignal.timeout(10000),
    });
    if (!r.ok) {
      const txt = await r.text().catch(() => "");
      return res.status(502).json({ success: false, error: `Open-Meteo HTTP ${r.status}`, detail: txt.slice(0, 200) });
    }

    const data  = await r.json();
    const h     = data.hourly;
    const times = h.time || [];

    const points = times.map((t, i) => {
      const precip = h.precipitation?.[i]   ?? 0;
      const temp   = h.temperature_2m?.[i]  ?? 10;
      const precipType = precip > 0
        ? (temp < 0 ? "snow" : "rain")
        : "";

      return {
        time:       t + ":00", // Open-Meteo gives "YYYY-MM-DDTHH:MM" without seconds
        cloudCover: h.cloudcover?.[i]            ?? 0,
        precip,
        precipType,
        humidity:   h.relativehumidity_2m?.[i]   ?? 0,
        windSpeed:  h.windspeed_10m?.[i]          ?? 0,
        windDir:    h.winddirection_10m?.[i]      ?? 0,
        temp,
      };
    });

    res.json({ success: true, points, lat, lon });
  } catch (err) {
    res.status(502).json({
      success: false,
      error: err.name === "AbortError" ? "Open-Meteo request timed out" : err.message,
    });
  }
});

// ── Plan Tonight: query Astro-COLIBRI /latest_transients + compute visibility ─
//
// GET /api/tonight?uid=...&days=7&lat=48.5&lon=2.3&minAlt=20
//
// Calls the Astro-COLIBRI /latest_transients endpoint, then for every returned
// event computes tonight's visibility window using our existing serverComputeAltAz
// helper.  Only events with maxAlt ≥ minAlt are returned, sorted by maxAlt desc.

app.get("/api/tonight", async (req, res) => {
  const uid    = String(req.query.uid || getSetting("colibri_uid", "") || process.env.COLIBRI_UID || "").trim();
  const days   = Math.min(30, Math.max(1, parseInt(req.query.days)  || 7));
  const lat    = parseFloat(req.query.lat);
  const lon    = parseFloat(req.query.lon);
  const minAlt = parseFloat(req.query.minAlt) || 20;

  if (!uid) {
    return res.status(400).json({ success: false, error: "Astro-COLIBRI user ID (uid) is required" });
  }

  // ── Time range: last `days` days up to now ──────────────────────────────
  const now     = new Date();
  const minDate = new Date(now.getTime() - days * 86400000);
  // Pass an explicit filter so GRBs, SNe, AT objects and other high-energy
  // transients are always included regardless of the user's saved app preferences.
  const body    = {
    uid,
    time_range: {
      min: minDate.toISOString().slice(0, 19),
      max: now.toISOString().slice(0, 19),
    },
  };

  let events;
  try {
    const controller = new AbortController();
    const timeout = setTimeout(() => controller.abort(), 20000);
    const r = await fetch("https://astro-colibri.science/latest_transients", {
      method:  "POST",
      headers: { "Content-Type": "application/json" },
      body:    JSON.stringify(body),
      signal:  controller.signal,
    });
    clearTimeout(timeout);

    if (!r.ok) {
      const txt = await r.text().catch(() => "");
      return res.status(502).json({
        success: false,
        error: `Astro-COLIBRI returned HTTP ${r.status}`,
        detail: txt.slice(0, 200),
      });
    }

    const data = await r.json();

    // The API returns {"message": "..."} for auth/input errors (bad UID, rate limit, etc.)
    if (data?.message && !data?.voevents) {
      return res.status(400).json({
        success: false,
        error: `Astro-COLIBRI: ${data.message}`,
      });
    }

    // Real response key is "voevents"
    events = Array.isArray(data?.voevents)      ? data.voevents
           : Array.isArray(data)                ? data
           : Array.isArray(data?.events)        ? data.events
           : Array.isArray(data?.transients)    ? data.transients
           : Array.isArray(data?.results)       ? data.results
           : [];
  } catch (err) {
    return res.status(502).json({
      success: false,
      error: err.name === "AbortError" ? "Astro-COLIBRI request timed out" : err.message,
    });
  }

  // ── Compute tonight's visibility for each event ─────────────────────────
  const hasObs = !isNaN(lat) && !isNaN(lon);

  // Approximate astronomical night: dusk ~ 21:00 local, dawn ~ 05:00 local
  // We use simple fixed UTC offsets from today; accurate enough for filtering.
  const today       = new Date(now);
  today.setUTCHours(18, 0, 0, 0);
  const nightStartUtc = today;
  const nightEndUtc   = new Date(today.getTime() + 12 * 3600000); // +12h

  const enriched = events
    .filter(e => e.ra != null && e.dec != null)
    .map(e => {
      const name   = e.source_name || e.name || e.trigger_id || "Unknown";
      // Use type first; classification on SNe contains subtype ("Ia","II") not the category
      const type   = (e.type || e.event_type || e.classification || "—").toLowerCase();
      // transient_flux is the optical magnitude for SNe; last_mag is often null
      const lastMag = e.transient_flux ?? e.last_mag ?? e.mag ?? null;
      // timestamp is in milliseconds
      const discDate = e.timestamp
        ? new Date(e.timestamp > 1e12 ? e.timestamp : e.timestamp * 1000).toISOString().slice(0, 16).replace("T", " ")
        : (e.time ? String(e.time).slice(0, 16).replace("T", " ") : "—");
      // Build source URL depending on event type
      const triggerId = e.trigger_id || "";
      const evType = (e.type || e.classification || "").toLowerCase();
      const isGrb  = evType === "grb" || /^GRB/i.test(triggerId) || /^GRB/i.test(name);
      const isTns  = /^TNS/i.test(triggerId) || /^(AT|SN)\s/i.test(name);
      const tnsId  = triggerId.replace(/^TNS/i, "") || name.replace(/^(AT|SN)\s*/i, "");
      // Prefer the embedded GCN/ATel URL when present; fall back by type
      const sourceUrl = e.url
        ? e.url
        : isGrb
          ? `https://gcn.nasa.gov/circulars?query=${encodeURIComponent(name)}`
          : isTns
            ? `https://www.wis-tns.org/object/${tnsId}`
            : `https://astro-colibri.com/?trigger_id=${triggerId}`;
      const tnsUrl = sourceUrl; // kept as tnsUrl for backward compat with frontend
      // AstroColibri deep-link: astro-colibri.com is the Flutter web app (not .science
      // which is just the landing page). The app reads ?trigger_id= for routing.
      const acId   = e.astro_colibri_id || "";
      const colUrl = `https://astro-colibri.com/?trigger_id=${encodeURIComponent(triggerId || acId)}`;

      let maxAlt = null, transitTime = null, windowStart = null, windowEnd = null, windowMin = 0;

      if (hasObs) {
        const STEP_MS = 5 * 60 * 1000;
        const steps   = Math.round(12 * 3600000 / STEP_MS);
        let   best    = -Infinity;
        let   inWin   = false;

        for (let s = 0; s <= steps; s++) {
          const t   = new Date(nightStartUtc.getTime() + s * STEP_MS);
          const { alt } = serverComputeAltAz(e.ra, e.dec, lat, lon, t);

          if (alt > best) { best = alt; transitTime = t; }
          if (alt >= minAlt) {
            if (!inWin) { windowStart = t; inWin = true; }
            windowEnd = t;
          } else {
            inWin = false;
          }
        }
        maxAlt    = best;
        windowMin = (windowStart && windowEnd)
          ? Math.round((windowEnd - windowStart) / 60000) : 0;
      }

      // err is in degrees; keep raw value for adaptive display
      const errArcsec = e.err > 0 ? e.err : null;

      // Age since discovery in fractional hours
      const discMs = e.timestamp
        ? (e.timestamp > 1e12 ? e.timestamp : e.timestamp * 1000)
        : (e.time ? new Date(e.time).getTime() : null);
      const ageHours = discMs ? (now.getTime() - discMs) / 3600000 : null;

      return {
        name, type, ra: e.ra, dec: e.dec,
        lastMag:   lastMag != null ? Number(lastMag).toFixed(1) : "—",
        discDate,
        ageHours:  ageHours != null ? Math.round(ageHours * 10) / 10 : null,
        redshift:  e.redshift > 0 ? e.redshift : null,
        errDeg: errArcsec,
        observatory: e.observatory || null,
        gcnUrl:    e.url || null,
        maxAlt:    maxAlt != null ? Math.round(maxAlt) : null,
        transitUtc: transitTime ? transitTime.toISOString().slice(11, 16) : "—",
        windowStart: windowStart ? windowStart.toISOString().slice(11, 16) : "—",
        windowEnd:   windowEnd   ? windowEnd.toISOString().slice(11, 16)   : "—",
        windowMin,
        tnsUrl, colUrl, acId,
      };
    })
    .filter(e => !hasObs || e.maxAlt >= minAlt)
    .sort((a, b) => (b.maxAlt ?? 0) - (a.maxAlt ?? 0));

  res.json({ success: true, count: enriched.length, events: enriched, days, uid });
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

  // Safety gate for close: verify mount is parked AND dust cover is closed
  if (command === "close") {
    const ninaTarget = {
      host:     DEFAULT_NINA_HOST,
      port:     DEFAULT_NINA_PORT,
      protocol: "http",
    };
    const blocks = [];
    try {
      const { equipmentInfo } = await getEquipmentInfo(ninaTarget);
      const mount = equipmentInfo?.Mount;
      const fd    = equipmentInfo?.FlatDevice;

      if (mount?.Connected && mount.AtPark !== true) {
        blocks.push(`mount is not parked (AtPark=${mount.AtPark})`);
      }
      if (fd?.Connected && normCoverState(fd.CoverState) !== 1) {
        const label = { 0: "Not Present", 2: "Moving", 3: "Open", 100: "Unknown", 101: "Error" }[normCoverState(fd.CoverState)] ?? String(fd.CoverState);
        blocks.push(`dust cover is not closed (${label})`);
      }
    } catch {
      // NINA unreachable — warn but don't hard-block (manual emergency close must still work)
    }

    if (blocks.length) {
      return res.status(403).json({
        success: false,
        error: `Roof close blocked: ${blocks.join(" and ")}`,
        blocks,
      });
    }
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

// ── Rain radar forecast (Open-Meteo, free, no key) ───────────────────────────
// GET /api/weather/forecast  → next 2h in 15-min slots + rain warning level
const OBS_LAT = parseFloat(process.env.OBS_LAT || "47.75");
const OBS_LON = parseFloat(process.env.OBS_LON || "-2.83");

let _rainForecastCache = null;
let _rainForecastTs    = 0;
const RAIN_FORECAST_TTL_MS = 10 * 60 * 1000; // refresh every 10 min

async function fetchRainForecast() {
  if (_rainForecastCache && Date.now() - _rainForecastTs < RAIN_FORECAST_TTL_MS) {
    return _rainForecastCache;
  }
  const url = `https://api.open-meteo.com/v1/forecast?latitude=${OBS_LAT}&longitude=${OBS_LON}` +
    `&minutely_15=precipitation,precipitation_probability&forecast_days=1&timezone=UTC`;
  const r = await fetch(url, { signal: AbortSignal.timeout(8000) });
  if (!r.ok) throw new Error(`Open-Meteo HTTP ${r.status}`);
  const d = await r.json();
  const times = d.minutely_15?.time || [];
  const prob  = d.minutely_15?.precipitation_probability || [];
  const prec  = d.minutely_15?.precipitation || [];
  const now = new Date().toISOString().slice(0, 16);
  const startIdx = Math.max(0, times.findIndex(t => t >= now));
  // Next 8 slots = 2 hours
  const slots = [];
  for (let i = startIdx; i < Math.min(startIdx + 8, times.length); i++) {
    slots.push({ time: times[i], prob: prob[i] ?? null, prec: prec[i] ?? null });
  }
  // Warning level: max prob in next 30 min (2 slots)
  const next30Prob = slots.slice(0, 2).map(s => s.prob ?? 0);
  const next60Prob = slots.slice(0, 4).map(s => s.prob ?? 0);
  const maxProb30 = Math.max(...next30Prob);
  const maxProb60 = Math.max(...next60Prob);
  const warnLevel = maxProb30 >= 70 ? "imminent" : maxProb60 >= 50 ? "likely" : maxProb60 >= 25 ? "possible" : "clear";
  const result = { slots, maxProb30, maxProb60, warnLevel, fetchedAt: new Date().toISOString() };
  _rainForecastCache = result;
  _rainForecastTs    = Date.now();
  return result;
}

app.get("/api/weather/forecast", async (req, res) => {
  try {
    const forecast = await fetchRainForecast();
    res.json({ success: true, ...forecast });
  } catch (e) {
    res.status(500).json({ success: false, error: e.message });
  }
});

app.get("/api/weather/forecast/history", (req, res) => {
  const hours = Math.min(Number(req.query.hours) || 48, 168);
  const rows = db.prepare(
    `SELECT id, ts, warn_level, max_prob_30, max_prob_60, slots_json
     FROM rain_forecast_log
     WHERE ts >= datetime('now', ? || ' hours')
     ORDER BY ts ASC`
  ).all(`-${hours}`);
  res.json({ success: true, history: rows.map(r => ({ ...r, slots: JSON.parse(r.slots_json || "[]") })), hours });
});

app.get("/api/ocs/history", (req, res) => {
  const hours = Math.min(Number(req.query.hours) || 48, 168); // max 7 days
  const rows = db.prepare(
    `SELECT id, ts, roof, safe, rain, temp, humidity, pressure, sky, ir_sky
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

    const result = db.prepare(
      `INSERT INTO ocs_history (roof, safe, rain, temp, humidity, pressure, sky, ir_sky)
       VALUES (?, ?, ?, ?, ?, ?, ?, ?)`
    ).run(
      ocsClean(fields.roof_sta)  || null,
      ocsClean(fields.stat_safe) || null,
      ocsClean(fields.wea_rain)  || null,
      ocsNum(fields.wea_temp),
      ocsNum(fields.wea_humd),
      ocsNum(fields.wea_pres),
      ocsNum(fields.wea_sq),
      ocsNum(fields.wea_irsky),
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

// ── OCS IR thermal camera ──────────────────────────────────────────────────────

app.get("/api/ocs/ir-image", async (req, res) => {
  try {
    const host = validateOcsHost(req.query?.ocsHost || DEFAULT_OCS_HOST);
    const response = await fetch(`http://${host}/ir-image.txt?x=${Date.now()}`, {
      signal: AbortSignal.timeout(5000),
    });
    if (!response.ok) return res.status(502).json({ success: false, error: `OCS returned ${response.status}` });

    const text = await response.text();
    const pixels = [];  // [{col, row, temp}]
    let avg = null;

    for (const line of text.trim().split("\n")) {
      const parts = line.trim().split(",");
      if (parts[0] === "avg") {
        avg = parseFloat(parts[1]);
      } else if (parts.length === 3) {
        const col = parseInt(parts[0], 10);
        const row = parseInt(parts[1], 10);
        const temp = parseFloat(parts[2]);
        if (!isNaN(col) && !isNaN(row) && !isNaN(temp)) {
          pixels.push({ col, row, temp });
        }
      }
    }

    res.json({ success: true, pixels, avg, cols: 16, rows: 4 });
  } catch (error) {
    res.status(502).json({
      success: false,
      error: error.name === "AbortError" || error.name === "TimeoutError"
        ? "OCS IR camera timed out"
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

// ── Alerts (GCN + AstroColibri) ──────────────────────────────────────────────
alerts.initAlertsTable(db);

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
  for (const j of jobs) {
    j.results = resultsByJob[j.id] || [];
    // Enrich each result with preview_url / stack_url if the files exist on disk
    for (const r of j.results) {
      const slug   = (r.target || "").toLowerCase().replace(/\s+/g, "_");
      const enc    = s => encodeURIComponent(s);
      // Check DATA_DIR first (with exposure subdir), then NAS_OUTPUT (flat)
      const candidates = [
        {
          stackPath:   path.join(DATA_DIR, r.obs_date||"", slug, r.filter||"", `${r.exposure||""}s`, "res.fit"),
          previewPath: path.join(DATA_DIR, r.obs_date||"", slug, r.filter||"", `${r.exposure||""}s`, "res_preview.png"),
          stackUrl:    `/data/${[r.obs_date, slug, r.filter, `${r.exposure}s`, "res.fit"].map(enc).join("/")}`,
          previewUrl:  `/data/${[r.obs_date, slug, r.filter, `${r.exposure}s`, "res_preview.png"].map(enc).join("/")}`,
        },
        {
          stackPath:   path.join(NAS_OUTPUT, r.obs_date||"", slug, r.filter||"", "res.fit"),
          previewPath: path.join(NAS_OUTPUT, r.obs_date||"", slug, r.filter||"", "res_preview.png"),
          stackUrl:    `/nas-output/${[r.obs_date, slug, r.filter, "res.fit"].map(enc).join("/")}`,
          previewUrl:  `/nas-output/${[r.obs_date, slug, r.filter, "res_preview.png"].map(enc).join("/")}`,
        },
      ];
      for (const c of candidates) {
        if (!r.stack_url   && fs.existsSync(c.stackPath))   r.stack_url   = c.stackUrl;
        if (!r.preview_url && fs.existsSync(c.previewPath)) r.preview_url = c.previewUrl;
        if (r.stack_url && r.preview_url) break;
      }
    }
  }
  res.json({ success: true, jobs });
});

app.delete("/api/pipeline/results/:id", (req, res) => {
  db.prepare("DELETE FROM pipeline_results WHERE id = ?").run(req.params.id);
  res.json({ success: true });
});

app.get("/api/pipeline/result/:id/download", (req, res) => {
  const row = db.prepare(
    "SELECT r.*, j.fits_dir FROM pipeline_results r JOIN pipeline_jobs j ON j.id = r.job_id WHERE r.id = ?"
  ).get(req.params.id);
  if (!row) return res.status(404).json({ error: "Result not found" });

  // date_str: use obs_date if set, otherwise extract from fits_dir (last numeric-date-like segment)
  let dateStr = row.obs_date;
  if (!dateStr && row.fits_dir) {
    const parts = row.fits_dir.split("/");
    const datePart = parts.find(p => /^\d{4}-\d{2}-\d{2}$/.test(p));
    dateStr = datePart || parts[parts.length - 2] || "unknown";
  }
  dateStr = dateStr || "unknown";

  const target = (row.target || "unknown").toLowerCase().replace(/\s+/g, "_");
  const filter = (row.filter || "unknown").toUpperCase();

  const fitsPath = path.join(NAS_OUTPUT, dateStr, target, filter, "res.fit");

  if (!fs.existsSync(fitsPath)) {
    return res.status(404).json({ error: `FITS file not found: ${fitsPath}` });
  }

  const filename = `${dateStr}_${target}_${filter}.fits`;
  res.setHeader("Content-Disposition", `attachment; filename="${filename}"`);
  res.setHeader("Content-Type", "application/fits");
  res.sendFile(path.resolve(fitsPath));
});

app.post("/api/pipeline/trigger", async (req, res) => {
  const { fits_dir, target, filter, exposure, selected_files, target_filter, manual_selection, force_fresh, use_color, refine_wcs } = req.body;
  if (!fits_dir) return res.status(400).json({ success: false, error: "fits_dir is required" });

  const useColorVal    = use_color    === false ? 0 : 1;
  const refineWcsVal   = refine_wcs   === false ? 0 : 1;
  const filterName     = target_filter || target || null;
  const stdwebTargetVal = req.body?.stdweb_use_target !== undefined
    ? (req.body.stdweb_use_target ? 1 : 0)
    : stdwebUseTargetForObject(filterName);

  const result = db.prepare(
    "INSERT INTO pipeline_jobs (target, filter, exposure, fits_dir, status, target_filter, use_color, refine_wcs, stdweb_use_target) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)"
  ).run(target || null, filter || null, exposure || null, fits_dir, "queued", target_filter || null, useColorVal, refineWcsVal, stdwebTargetVal);

  const job_id = result.lastInsertRowid;

  callProcessingService(job_id, fits_dir, target || "Unknown", selected_files || null, target_filter || null, !!manual_selection, !!force_fresh, !!stdwebTargetVal)
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

  const stdwebTargetVal = oldJob.stdweb_use_target != null
    ? oldJob.stdweb_use_target
    : stdwebUseTargetForObject(oldJob.target_filter || oldJob.target);

  const result = db.prepare(
    "INSERT INTO pipeline_jobs (target, filter, exposure, fits_dir, status, target_filter, use_color, refine_wcs, stdweb_use_target) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)"
  ).run(oldJob.target, oldJob.filter, oldJob.exposure, oldJob.fits_dir, "queued", oldJob.target_filter || null, oldJob.use_color ? 1 : 0, oldJob.refine_wcs != null ? oldJob.refine_wcs : 1, stdwebTargetVal);

  const job_id = result.lastInsertRowid;

  callProcessingService(job_id, oldJob.fits_dir, oldJob.target || "Unknown", null, oldJob.target_filter || null, false, false, !!stdwebTargetVal)
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

/** For a SNAPSHOT folder, scan FITS files and return unique real target names.
 *  Stops early once the target set has been stable for STABLE_STREAK files in
 *  a row, avoiding full scans of large nights (200+ frames) over slow NAS. */
function snapshotTargets(folderPath) {
  try {
    const SKIP_OBJECTS = new Set(["snapshot", "dark", "flat", "bias", "test_target",
                                   "test_target_v2", "test", "unknown", "none", ""]);
    const MAX_FILES    = 200;  // hard cap — covers long nights with many targets
    const STABLE_STREAK = 30; // only stop early if we haven't found a NEW target in 30 files
                               // AND we already have at least 2 distinct targets (avoids
                               // breaking out early when the first target has many frames)
    const files = fs.readdirSync(folderPath).filter((f) => /\.(fit|fits)$/i.test(f));
    const seen = new Set();
    let streak = 0;
    for (const f of files.slice(0, MAX_FILES)) {
      const obj = readFitsObject(path.join(folderPath, f));
      if (obj && !SKIP_OBJECTS.has(obj.toLowerCase()) && !seen.has(obj)) {
        seen.add(obj);
        streak = 0;
      } else {
        // Only allow early exit once we have ≥2 distinct targets, to avoid stopping
        // prematurely when the first target fills many consecutive files.
        if (++streak >= STABLE_STREAK && seen.size >= 2) break;
      }
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

app.delete("/api/pipeline/jobs/:id", async (req, res) => {
  const jobId = req.params.id;
  const PROC_URL = process.env.PROCESSING_URL || "http://localhost:5200";

  // Tell the processing service to cancel the running/queued job so it
  // kills any active Siril subprocess and moves to the next job.
  try {
    await fetch(`${PROC_URL}/cancel/${jobId}`, { method: "POST", signal: AbortSignal.timeout(5000) });
  } catch { /* non-fatal — maybe already finished or service unreachable */ }

  // Delete the job (pipeline_results cascade-deleted via FK ON DELETE CASCADE)
  db.prepare("DELETE FROM pipeline_jobs WHERE id = ?").run(jobId);

  res.json({ success: true });
});

// ── Manual frame selection endpoints ─────────────────────────────────────────
// GET  /api/pipeline/jobs/:id/selection  → frame list + stats
// POST /api/pipeline/jobs/:id/selection  → confirm selection

app.get("/api/pipeline/jobs/:id/selection", async (req, res) => {
  try {
    const r = await fetch(`${PROCESSING_SERVICE_URL}/jobs/${req.params.id}/selection`);
    const text = await r.text();
    let data;
    try { data = JSON.parse(text); } catch { data = { success: false, error: `Service error ${r.status}` }; }
    res.status(r.ok ? 200 : r.status).json(data);
  } catch (err) {
    res.status(502).json({ success: false, error: err.message });
  }
});

app.post("/api/pipeline/jobs/:id/selection", async (req, res) => {
  try {
    const r = await fetch(`${PROCESSING_SERVICE_URL}/jobs/${req.params.id}/selection`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(req.body),
    });
    const text = await r.text();
    let data;
    try { data = JSON.parse(text); } catch { data = { success: false, error: `Service error ${r.status}` }; }
    res.status(r.ok ? 200 : r.status).json(data);
  } catch (err) {
    res.status(502).json({ success: false, error: err.message });
  }
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

// ── STDWeb health check ───────────────────────────────────────────────────────
// GET /api/stdweb/health  — returns { reachable, url, status? }
app.get("/api/stdweb/health", async (req, res) => {
  const STDWEB_URL = process.env.STDWEB_URL || "http://86.253.141.183:7000";
  const tok = getStdwebToken();
  if (!tok) {
    return res.status(400).json({ reachable: false, error: "STDWEB_TOKEN is not configured", url: STDWEB_URL });
  }
  try {
    const r = await fetch(`${STDWEB_URL}/api/tasks/`, {
      headers: { "Authorization": `Token ${tok}` },
      signal: AbortSignal.timeout(5000),
    });
    res.json({ reachable: r.ok || r.status < 500, httpStatus: r.status, url: STDWEB_URL });
  } catch (err) {
    res.json({ reachable: false, error: err.message, url: STDWEB_URL });
  }
});

// ── STDWeb photometry proxy ───────────────────────────────────────────────────
// GET /api/stdweb/task/:task_id/photometry[?refresh=1]
// Returns photometry measurements for a task. Serves from DB when already
// fully populated; fetches from STDWeb live otherwise and persists the result.
// Pass ?refresh=1 to force a live re-fetch (e.g. user just changed photometry
// params on STDWeb and wants the new numbers).
app.get("/api/stdweb/task/:task_id/photometry", async (req, res) => {
  const { task_id } = req.params;
  const refresh = req.query.refresh === "1" || req.query.refresh === "true";
  const STDWEB_URL = process.env.STDWEB_URL || "http://86.253.141.183:7000";
  const tok = getStdwebToken();
  if (!tok) {
    return res.status(400).json({ success: false, error: "STDWEB_TOKEN is not configured" });
  }
  const headers = { "Authorization": `Token ${tok}` };

  try {
    // Serve from DB only when photometry AND subtraction data are both present
    // AND the caller didn't ask for a forced refresh.
    if (!refresh) {
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
    }

    // Fetch live from STDWeb
    const taskJson = await fetch(`${STDWEB_URL}/api/tasks/${task_id}/`, { headers })
      .then(r => r.json()).catch(() => ({}));
    const filter = taskJson?.config?.filter || "?";

    const html = await fetch(`${STDWEB_URL}/tasks/${task_id}`, { headers })
      .then(r => r.text()).catch(() => "");

    // Match both plain photometry and color-term photometry, e.g.
    //   "Primary target magnitude is BPmag = 14.26 +/- 0.01"
    //   "Primary target magnitude is BPmag - 0.03 (BPmag - RPmag) = 14.26 +/- 0.01"
    // The "Target magnitude is ..." line is anchored to a line start so we don't
    // accidentally match the "Primary target magnitude is ..." line above it.
    // The upper-limit regex tolerates both literal '>' and the HTML entity '&gt;'.
    const mjdMatch    = html.match(/MJD is ([\d.]+)/);
    const directMatch = html.match(/Primary target magnitude is .+? = ([\d.]+) \+\/- ([\d.]+)/);
    const subMatch    = html.match(/(?:^|\n)\s*Target magnitude is .+? = ([\d.]+) \+\/- ([\d.]+)/m);
    const subUlMatch  = !subMatch && html.match(/(?:^|\n)\s*Target magnitude upper limit is .+? (?:>|&gt;) ([\d.]+)/m);

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
        // Use COALESCE(new, existing) for each mag field so a partial parse
        // (e.g. subtraction not re-run yet) doesn't wipe previously stored
        // values. obs_date keeps its existing value if already set.
        db.prepare(
          `UPDATE pipeline_results
           SET mjd        = ?,
               obs_date   = COALESCE(obs_date,   ?),
               mag_ap     = COALESCE(?, mag_ap),
               magerr_ap  = COALESCE(?, magerr_ap),
               mag_sub    = COALESCE(?, mag_sub),
               magerr_sub = COALESCE(?, magerr_sub),
               mag_sub_ul = COALESCE(?, mag_sub_ul),
               updated_at = datetime('now')
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

// ── GCN-Kafka alert listener ─────────────────────────────────────────────────
function resolveObserverSite() {
  const envLat = parseFloat(process.env.OBS_LAT || "");
  const envLon = parseFloat(process.env.OBS_LON || "");
  if (Number.isFinite(envLat) && Number.isFinite(envLon)) return { lat: envLat, lon: envLon };
  const sLat = parseFloat(getSetting("obs_lat", ""));
  const sLon = parseFloat(getSetting("obs_lon", ""));
  if (Number.isFinite(sLat) && Number.isFinite(sLon)) return { lat: sLat, lon: sLon };
  return { lat: 47.7, lon: -3.0 };
}

// Stored reference so /api/alerts/inject can reuse the same callback
function _gcnPushToQueue({ name, raDeg, decDeg, alert, mode }) {
    const isToo = (mode === "too");
    const shouldNotify = getSetting("alert_notify_email", "true") !== "false";
    if (isToo && seqState.alertReady && seqState.running && !seqState.tooRunning && alert) {
      seqState.tooInterrupt = alert;
      seqLog(`🚨 ToO ALERT received: ${alert.broker} ${alert.trigger_id} — interrupt will fire at next checkAbort`, "warn");
      if (shouldNotify) sendAlertEmail(
        `🚨 ToO Alert: ${alert.broker} ${alert.trigger_id}`,
        `A Target-of-Opportunity alert has been received and will interrupt the current sequence.\n\n` +
        `Broker:     ${alert.broker}\n` +
        `Trigger ID: ${alert.trigger_id}\n` +
        `RA:         ${alert.ra?.toFixed?.(4) ?? "?"}°\n` +
        `Dec:        ${alert.dec?.toFixed?.(4) ?? "?"}°\n` +
        `Error:      ${alert.err_deg?.toFixed?.(3) ?? "?"}°\n` +
        `Time:       ${new Date().toISOString()}\n`
      );
      return;
    }
    seqState.queue.push({
      name, raDeg, decDeg, ra: "", dec: "",
      done: false, addedAt: Date.now(), source: "alert",
    });
    if (shouldNotify) sendAlertEmail(
      `📡 Alert queued: ${name}`,
      `A new alert has been added to the Night Plan queue.\n\n` +
      `Target:     ${name}\n` +
      `RA:         ${raDeg?.toFixed?.(4) ?? "?"}°\n` +
      `Dec:        ${decDeg?.toFixed?.(4) ?? "?"}°\n` +
      `Time:       ${new Date().toISOString()}\n`
    );
}

alerts.startGcnListener({
  db,
  getStrategy,
  pushToQueue: (opts) => {
    _gcnPushToQueue(opts);
  },
  saveQueue,
  computeAltAz:  serverComputeAltAz,
  getMountSite:  resolveObserverSite,
  minAltDeg:     parseFloat(process.env.ALERT_MIN_ALT_DEG || "25"),
  minMoonSepDeg: parseFloat(process.env.ALERT_MIN_MOON_SEP_DEG || "20"),
});


app.listen(PORT, () => {
  console.log(`NINA control server listening on http://localhost:${PORT}`);
});
