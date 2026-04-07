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
    sky       REAL,
    ir_sky    REAL
  );
`);
// Migration: add ir_sky column if DB was created before this version
try { db.exec(`ALTER TABLE ocs_history ADD COLUMN ir_sky REAL`); } catch (_) {}

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
    colibriUid: getSetting("colibri_uid", ""),
  });
});

// GET  /api/settings/:key
// POST /api/settings  { key, value }
const SETTING_DEFAULTS = {
  minAlt:               "20",
  zenithLimit:          "70",
  meridianGap:          "10",
  minStars:             "10",
  frameCheckEnabled:    "false",
  frameCheckThreshold:  "5",
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
  setSetting(key, value ?? "");
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
  maxBadMin:    90,    // abort sequence if sky has been bad for this long (dawn protection)

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
      if (cfg.maxBadMin    !== undefined) watchdog.maxBadMin    = Number(cfg.maxBadMin);
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

  if (conditionsBad) {
    watchdog.clearSince = null;
    if (watchdog.state !== "bad") {
      watchdog.state    = "bad";
      watchdog.badSince = watchdog.badSince ?? new Date(); // keep first onset time
      console.log(`[watchdog] Conditions bad — ${watchdog.lastMsg}`);
      // Park mount + close dust cover if sequence is running and not already parked
      if (seqState.running && !watchdog.parked) {
        console.log("[watchdog] Parking mount due to bad conditions...");
        try {
          const cfg = seqState.ninaConfig;
          if (cfg) {
            await fetch(`http://${cfg.host}:${cfg.port}/v2/api/equipment/mount/park`, {
              method: "GET", signal: AbortSignal.timeout(15000),
            }).catch(() => {});
          }
          watchdog.parked = true;
          console.log("[watchdog] Mount parked.");
        } catch (e) {
          console.error("[watchdog] Park failed:", e.message);
        }

        // Close dust cover immediately after parking — protect optics from rain
        try {
          const cfg = seqState.ninaConfig;
          if (cfg) {
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
              // Poll up to 30 s for confirmation
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
          }
        } catch (e) {
          console.error("[watchdog] Dust cover close failed:", e.message);
        }
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
          console.log("[watchdog] Retention passed — unparking mount...");
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
                    const roofRaw = ocsClean(fields?.stat_roof ?? "").toLowerCase();
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
    watchdog._interval = setInterval(runWatchdogCheck, 60 * 1000);
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
    enabled:      watchdog.enabled,
    skyTempLimit: watchdog.skyTempLimit,
    sqLimit:      watchdog.sqLimit,
    retentionMin: watchdog.retentionMin,
    maxBadMin:    watchdog.maxBadMin,
    state:        watchdog.state,
    lastCheck:    watchdog.lastCheck,
    lastMsg:      watchdog.lastMsg,
    parked:         watchdog.parked,
    coverClosed:    watchdog.coverClosed,
    badSince:       watchdog.badSince,
    safetyOverride: watchdog.safetyOverride,
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
  watchdog.maxBadMin    = Math.max(10, Number(b.maxBadMin)   || 90);
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
    maxBadMin: watchdog.maxBadMin,
  }));
  applyWatchdogInterval();
  res.json({
    success: true,
    watchdog: {
      enabled:        watchdog.enabled,
      skyTempLimit:   watchdog.skyTempLimit,
      sqLimit:        watchdog.sqLimit,
      retentionMin:   watchdog.retentionMin,
      maxBadMin:      watchdog.maxBadMin,
      state:          watchdog.state,
      lastCheck:      watchdog.lastCheck,
      lastMsg:        watchdog.lastMsg,
      parked:         watchdog.parked,
      badSince:       watchdog.badSince,
      safetyOverride: watchdog.safetyOverride,
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
          const roofRaw = ocsClean(fields?.stat_roof ?? "").toLowerCase();
          if (roofRaw && !roofRaw.includes("open")) {
            return res.status(403).json({
              success: false,
              error: `Roof is not open (OCS: "${ocsClean(fields.stat_roof)}") — cover open blocked for safety`,
              roofStatus: ocsClean(fields.stat_roof),
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
      coverStateLabel: stateNames[coverState] ?? String(coverState),
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
    return res.status(isNinaApiSuccess(cr) ? 200 : 502).json({
      success: isNinaApiSuccess(cr),
      flatDevice: equipmentInfo?.FlatDevice ?? null,
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
};

loadQueue();

// ── Per-filter autofocus cache ────────────────────────────────────────────────
// Stores the best known focuser position per filter, with timestamp.
// Entries expire after FILTER_AF_CACHE_MS; within that window the focuser
// is simply moved to the cached position instead of running a full AF run.
const FILTER_AF_CACHE_MS = 2 * 60 * 60 * 1000; // 2 hours
let filterAfCache = {};  // { filterName: { pos: number, ts: number } }

(function loadFilterAfCache() {
  try {
    const v = getSetting("filterAfCache");
    if (v) filterAfCache = JSON.parse(v);
  } catch { filterAfCache = {}; }
})();

function saveFilterAfCache() {
  try { setSetting("filterAfCache", JSON.stringify(filterAfCache)); } catch {}
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
  // Watchdog takes priority: if sky is bad, wait here — do NOT throw __LIMIT_REACHED__
  // because that would consume a retry slot; instead block until the watchdog clears.
  if (watchdog.enabled && (watchdog.state === "bad" || watchdog.state === "recovering")) {
    seqLog("  Safety check: watchdog reports bad sky — holding frame until sky clears", "warn");
    while (watchdog.state === "bad" || watchdog.state === "recovering") {
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
  const maxBadMs = (watchdog.maxBadMin ?? 90) * 60_000;
  while (watchdog.state === "bad" || watchdog.state === "recovering") {
    checkAbort();
    // If conditions have been bad longer than maxBadMin → dawn or persistent closure → abort
    if (watchdog.badSince && (Date.now() - watchdog.badSince.getTime()) > maxBadMs) {
      const badMin = Math.round((Date.now() - watchdog.badSince.getTime()) / 60_000);
      seqLog(`🌅 Watchdog: sky has been bad for ${badMin} min (limit ${watchdog.maxBadMin} min) — aborting sequence (likely dawn)`, "warn");
      seqState.aborted = true;
      checkAbort(); // throws __ABORTED__
    }
    const badMin = watchdog.badSince ? Math.round((Date.now() - watchdog.badSince.getTime()) / 60_000) : "?";
    seqState.currentStep = `Watchdog: waiting for clear sky — bad for ${badMin} min / max ${watchdog.maxBadMin} min (${watchdog.lastMsg || watchdog.state})`;
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
    // sync (writes correction into mount model) → re-slew → repeat until within
    // its own tolerance. We don't trust the HTTP response (NINA often returns
    // Success=false on intermediate slews); instead we poll mount state.
    try {
      await callNinaLong(target, "/v2/api/equipment/mount/slew", {
        ra: raDeg, dec: decDeg, waitForResult: true, center: true, name,
      }, 10 * 60 * 1000);
    } catch (fetchErr) {
      seqLog(`  Centering request error: ${fetchErr.message} — falling back to plain slew`, "warn");
      await plainSlew();
      seqLog("  Plain slew complete (no centering)");
    }

    // Poll until the mount has been continuously NOT slewing for SETTLE_STABLE ms.
    // This catches the full solve→sync→re-slew cycle: the mount may stop briefly
    // after the first goto, then start again after the sync. We only declare done
    // once it has been truly still for a sustained period.
    seqState.currentStep = `Centering ${name} — waiting for convergence...`;
    const SETTLE_TIMEOUT  = 10 * 60 * 1000;  // hard cap
    const POLL_INTERVAL   = 2500;
    const SETTLE_STABLE   = 5000;             // must be not-slewing for 5 s to be "done"
    const deadline  = Date.now() + SETTLE_TIMEOUT;
    let stableMs    = 0;
    let lastErrArcsec = null;

    while (Date.now() < deadline) {
      checkAbort();
      await new Promise(r => setTimeout(r, POLL_INTERVAL));

      const { equipmentInfo: ei } = await getEquipmentInfo(target);
      const m = ei?.Mount;
      if (!m) continue;

      if (m.RightAscension != null && m.Declination != null) {
        const errArcsec = _angularSepArcsec(
          m.RightAscension * 15, m.Declination,   // NINA returns RA in hours
          raDeg, decDeg,
        );
        const errRounded = Math.round(errArcsec);
        if (lastErrArcsec === null || Math.abs(errRounded - lastErrArcsec) > 10) {
          seqLog(`  Centering: offset ${errRounded}"  (Slewing=${m.Slewing})`);
          lastErrArcsec = errRounded;
          seqState.currentStep = `Centering ${name} — offset ${errRounded}"`;
        }
      }

      if (m.Slewing === true) {
        stableMs = 0;  // reset — still moving
      } else {
        stableMs += POLL_INTERVAL;
        if (stableMs >= SETTLE_STABLE) break;  // sustained stillness → done
      }
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
  const MAX_ATTEMPTS   = 3;
  const SETTLE_TIMEOUT = 3 * 60 * 1000;  // 3 min per attempt to reach "Guiding"
  const RETRY_PAUSE    = 15_000;          // 15 s between stop → restart

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
          return s === "Guiding" || s === "LostLock";
        },
        SETTLE_TIMEOUT,
        3000,
      );
      seqLog("Guiding active ✓");
      return; // success
    } catch {
      if (attempt < MAX_ATTEMPTS) {
        seqLog(`  PHD2 did not settle within ${SETTLE_TIMEOUT / 60000} min — retrying...`, "warn");
      } else {
        seqLog(`⚠ PHD2 failed to settle after ${MAX_ATTEMPTS} attempts — continuing without guiding`, "warn");
      }
    }
  }
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

async function stepOpenCover(target) {
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
 * - If a cached position exists that is < FILTER_AF_CACHE_MS old, just moves the
 *   focuser to that position instead of running a full AF sequence.
 * - On a full run, saves the result to filterAfCache for future reuse.
 * availableFilters: the filterwheel filter list from NINA equipment info.
 */
async function stepAutofocusForFilter(target, filterName, availableFilters) {
  // ── 1. Change to this filter before AF so NINA uses the right bandpass ────
  const filter = findFilterByName(availableFilters, filterName);
  if (filter) {
    await changeFilterVerified(target, filter);
    checkAbort();
  }

  // ── 2. Check per-filter cache ─────────────────────────────────────────────
  const cached = filterAfCache[filterName];
  const now    = Date.now();
  if (cached && Number.isFinite(cached.pos) && (now - cached.ts) < FILTER_AF_CACHE_MS) {
    const ageMin = Math.round((now - cached.ts) / 60000);
    seqLog(`AF [${filterName}]: cached pos ${cached.pos} (${ageMin} min ago) — moving focuser ✓`);
    seqState.currentStep = `AF [${filterName}]: moving to cached position ${cached.pos}...`;
    await stepMoveFocuserAbsolute(target, cached.pos);
    return;
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

  await new Promise(r => setTimeout(r, 12000));
  checkAbort();

  const AF_DEADLINE   = Date.now() + 12 * 60 * 1000;
  const STABLE_NEEDED = 12;    // 12 × 5 s = 60 s stable
  const MIN_MOVES     = 3;
  const MIN_AF_MS     = 90 * 1000;
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
        filterAfCache[filterName] = { pos, ts: Date.now() };
        saveFilterAfCache();
        return;
      }

      if (moveCount === 0 && pollCount >= NO_MOVE_LIMIT) {
        seqLog(`  No focuser movement — assuming AF skipped/complete [${filterName}]`);
        if (Number.isFinite(finalPos)) {
          filterAfCache[filterName] = { pos: finalPos, ts: Date.now() };
          saveFilterAfCache();
        }
        return;
      }
    } catch { /* ignore transient errors */ }
  }

  seqLog(`  AF [${filterName}] timed out (12 min) — proceeding`, "warn");
  if (Number.isFinite(finalPos)) {
    filterAfCache[filterName] = { pos: finalPos, ts: Date.now() };
    saveFilterAfCache();
  }
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
 * Returns { stars, hfr } or null if the endpoint is unavailable.
 */
async function fetchLastImageStats(target) {
  try {
    const res = await callNinaLong(
      target, "/v2/api/equipment/camera/capture/statistics", {}, 30000,
    );
    const r = res.body?.Response;
    if (!r) return null;
    const stars = r.Stars ?? r.DetectedStars ?? null;
    const hfr   = r.HFR ?? null;
    return (stars !== null) ? { stars: Number(stars), hfr: Number(hfr) } : null;
  } catch {
    return null;
  }
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
    const MIN_STARS      = safetyLimits.minStars ?? 10;
    const CLOUD_WAIT_MS  = 90_000;   // 90 s pause between retries
    const MAX_RETAKES    = 3;        // max cloud retakes before giving up on this frame

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

    let imgStats  = await fetchLastImageStats(target);
    let retakes   = 0;
    let frameGood = true;   // set false only if we exhaust retakes

    if (imgStats !== null && MIN_STARS > 0) {
      while (imgStats !== null && imgStats.stars < MIN_STARS && retakes < MAX_RETAKES) {
        retakes++;
        seqLog(
          `⛅ ${filterName} ${i}/${count}: ${imgStats.stars} stars < min ${MIN_STARS} ` +
          `— cloud passage detected. Waiting ${CLOUD_WAIT_MS / 1000}s then retaking ` +
          `(${retakes}/${MAX_RETAKES})...`,
          "warn",
        );
        seqState.currentStep = `⛅ Cloud — retake ${retakes}/${MAX_RETAKES} waiting...`;
        // Poll every 30 s — bail out early if watchdog says sky clear
        const frameDeadline = Date.now() + CLOUD_WAIT_MS;
        while (Date.now() < frameDeadline) {
          checkAbort();
          await new Promise(r => setTimeout(r, 30_000));
          if (watchdog.enabled && watchdog.state === "clear") {
            seqLog("☀ Sky cleared (watchdog) — retrying frame early");
            break;
          }
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
    minAlt = 20, meridianGap = 10, zenithLimit = 70, minStars = 10,
    frameCheckEnabled = false, frameCheckThresholdArcmin = 5, solveExp = 5,
  } = seqConfig;
  const safetyLimits = {
    minAlt, meridianGap, zenithLimit, minStars,
    frameCheckEnabled, frameCheckThresholdArcmin, solveExp,
  };

  // Maximum time to wait for any target to enter a valid sky zone before giving up
  const MAX_WAIT_NO_TARGET_MS = 45 * 60 * 1000; // 45 min
  const WAIT_POLL_INTERVAL_MS =  5 * 60 * 1000; //  5 min between re-checks

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

      // If NINA says the mount was already homed (skipped physical movement),
      // accept immediately — the unpark cleared AtHome but the mount didn't move.
      if (typeof resp === "string" && /already.?hom/i.test(resp)) {
        seqLog("Mount homed ✓ (already at home position)");
        return;
      }

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

    // ── Step 0: Check roof, unpark, go home ──────────────────────────────────
    try {
      // Verify OCS roof is open before moving anything
      const { ok: ocsOk, fields } = await ocsStatus(DEFAULT_OCS_HOST).catch(() => ({ ok: false, fields: {} }));
      const roofStatus = ocsClean(fields?.stat_roof ?? "").toLowerCase();
      const isClosed   = ocsOk && roofStatus && !roofStatus.includes("open");
      if (isClosed) {
        seqLog(`✗ Roof appears closed (OCS: "${ocsClean(fields.stat_roof)}") — aborting sequence for safety`, "error");
        return;
      }
      if (!ocsOk) {
        seqLog("OCS unreachable — cannot verify roof status, proceeding with caution", "warn");
      } else {
        seqLog(`Roof status: ${ocsClean(fields.stat_roof)} ✓`);
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
        await stepStartGuiding(ninaConfig);
        await waitForConfirmation(`Guiding active — proceed to captures?`);
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
          await stepAutofocusForFilter(ninaConfig, filterName, filters);
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
              await stepStartGuiding(ninaConfig);
              const remainingFilters = targetFilters.filter(f => !completedFilters.has(f));
              if (completedFilters.size > 0) seqLog(`  ↳ retry: skipping already-done filters [${[...completedFilters].join(", ")}], running: [${remainingFilters.join(", ")}]`);
              for (const filterName of remainingFilters) {
                checkAbort();
                await stepAutofocusForFilter(ninaConfig, filterName, reFilters);
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
  const uid    = String(req.query.uid  || getSetting("colibri_uid", "") || "").trim();
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

// ── STDWeb health check ───────────────────────────────────────────────────────
// GET /api/stdweb/health  — returns { reachable, url, status? }
app.get("/api/stdweb/health", async (req, res) => {
  const STDWEB_URL   = process.env.STDWEB_URL   || "http://86.253.141.183:7000";
  const STDWEB_TOKEN = process.env.STDWEB_TOKEN  || "1e296ddd6738af45467b7bc6558c00a9524447ab";
  try {
    const r = await fetch(`${STDWEB_URL}/api/tasks/`, {
      headers: { "Authorization": `Token ${STDWEB_TOKEN}` },
      signal: AbortSignal.timeout(5000),
    });
    res.json({ reachable: r.ok || r.status < 500, httpStatus: r.status, url: STDWEB_URL });
  } catch (err) {
    res.json({ reachable: false, error: err.message, url: STDWEB_URL });
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
