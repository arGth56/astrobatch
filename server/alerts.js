// server/alerts.js — GCN-Kafka listener + AstroColibri enrichment
"use strict";

const { Kafka } = require("gcn-kafka");

// ── Curated topic list ──────────────────────────────────────────────────────
const TOPICS = [
  // Swift
  "gcn.classic.text.SWIFT_BAT_GRB_POS_ACK",
  "gcn.classic.text.SWIFT_BAT_QL_POS",
  "gcn.classic.text.SWIFT_XRT_POSITION",
  "gcn.classic.text.SWIFT_XRT_CENTROID",
  "gcn.classic.text.SWIFT_UVOT_POS",
  // Fermi
  "gcn.classic.text.FERMI_GBM_FLT_POS",
  "gcn.classic.text.FERMI_GBM_GND_POS",
  "gcn.classic.text.FERMI_GBM_FIN_POS",
  "gcn.classic.text.FERMI_LAT_POS_INI",
  "gcn.classic.text.FERMI_LAT_POS_UPD",
  "gcn.classic.text.FERMI_LAT_GND",
  "gcn.classic.text.FERMI_LAT_OFFLINE",
  // IceCube / AMON
  "gcn.classic.text.ICECUBE_ASTROTRACK_GOLD",
  "gcn.classic.text.ICECUBE_ASTROTRACK_BRONZE",
  "gcn.classic.text.ICECUBE_CASCADE",
  "gcn.classic.text.AMON_ICECUBE_COINC",
  "gcn.classic.text.AMON_ICECUBE_EHE",
  "gcn.classic.text.AMON_ICECUBE_HESE",
  "gcn.classic.text.AMON_NU_EM_COINC",
  // LVK
  "gcn.classic.text.LVC_EARLY_WARNING",
  "gcn.classic.text.LVC_PRELIMINARY",
  "gcn.classic.text.LVC_INITIAL",
  "gcn.classic.text.LVC_UPDATE",
  "gcn.classic.text.LVC_COUNTERPART",
  "gcn.classic.text.LVC_RETRACTION",
  // INTEGRAL — science ops ended Feb 2025, re-entry 2029
  // "gcn.classic.text.INTEGRAL_WAKEUP",
  // "gcn.classic.text.INTEGRAL_REFINED",
  // "gcn.classic.text.INTEGRAL_OFFLINE",
  // "gcn.classic.text.INTEGRAL_WEAK",
  // AGILE — decommissioned Feb 2024 (re-entered atmosphere), topics kept for archive
  // "gcn.classic.text.AGILE_GRB_WAKEUP",
  // "gcn.classic.text.AGILE_GRB_GROUND",
  // "gcn.classic.text.AGILE_GRB_REFINED",
  // GECAM
  "gcn.classic.text.GECAM_FLT",
  "gcn.classic.text.GECAM_GND",
  // MAXI
  "gcn.classic.text.MAXI_UNKNOWN",
  "gcn.classic.text.MAXI_KNOWN",
  // IPN
  "gcn.classic.text.IPN_POS",
  // SN neutrino
  "gcn.classic.text.SK_SN",
  "gcn.classic.text.SNEWS",
  // Einstein Probe (JSON)
  "gcn.notices.einstein_probe.wxt.alert",
  // SVOM (VOEvent XML) — GRM omitted: spectrometer, never carries coordinates
  "gcn.notices.svom.voevent.eclairs",
  "gcn.notices.svom.voevent.mxt",
];

// ── DB table ────────────────────────────────────────────────────────────────
function initAlertsTable(db) {
  db.exec(`
    CREATE TABLE IF NOT EXISTS alerts (
      id              INTEGER PRIMARY KEY AUTOINCREMENT,
      received_at     TEXT DEFAULT (datetime('now')),
      event_time      TEXT,
      broker          TEXT,
      topic           TEXT,
      trigger_id      TEXT,
      classification  TEXT,
      ra              REAL,
      dec             REAL,
      err_deg         REAL,
      alt_now         REAL,
      moon_sep        REAL,
      action          TEXT,
      action_reason   TEXT,
      colibri_id      TEXT,
      raw             TEXT
    );
    CREATE INDEX IF NOT EXISTS idx_alerts_trigger ON alerts(trigger_id);
    CREATE INDEX IF NOT EXISTS idx_alerts_received ON alerts(received_at DESC);
  `);
}

// ── Parsers ─────────────────────────────────────────────────────────────────
function parseTextNotice(text) {
  const lines = String(text || "").split(/\r?\n/);
  const out = {};
  for (const ln of lines) {
    const m = ln.match(/^([A-Z0-9_\-]+):\s*(.*?)\s*$/i);
    if (m) out[m[1]] = m[2];
  }
  return out;
}

function extractFromText(kv) {
  const parseDeg = (s) => {
    if (!s) return null;
    const m = String(s).match(/(-?\d+(?:\.\d+)?)\s*d/);
    return m ? parseFloat(m[1]) : null;
  };
  const parseErr = (s) => {
    if (!s) return null;
    const m = String(s).match(/(-?\d+(?:\.\d+)?)\s*\[(arcmin|deg|arcsec)?\]?/i);
    if (!m) return null;
    const v = parseFloat(m[1]);
    const u = (m[2] || "arcmin").toLowerCase();
    return u === "deg" ? v : u === "arcsec" ? v / 3600 : v / 60;
  };
  const ra  = parseDeg(kv.GRB_RA || kv.SRC_RA || kv.RA);
  const dec = parseDeg(kv.GRB_DEC || kv.SRC_DEC || kv.DEC);
  const err = parseErr(kv.GRB_ERROR || kv.SRC_ERROR || kv.ERROR);
  const trigger_id =
    kv.TRIGGER_NUM || kv.TRIGGER_ID || kv.EVENT_ID || kv.NOTICE_DATE || null;
  const time = kv.GRB_DATE || kv.EVENT_DATE || kv.NOTICE_DATE || null;
  return { ra, dec, err, trigger_id, time };
}

function extractFromVoevent(xml) {
  const tag = (name) => {
    const m = xml.match(new RegExp(`<${name}[^>]*>([^<]*)</${name}>`, "i"));
    return m ? m[1].trim() : null;
  };
  const param = (name) => {
    const m = xml.match(new RegExp(`name="${name}"[^>]*value="([^"]*)"`, "i"));
    return m ? m[1] : null;
  };
  const c1m = xml.match(/<C1>([^<]+)<\/C1>/i);
  const c2m = xml.match(/<C2>([^<]+)<\/C2>/i);
  const errm = xml.match(/<Error2Radius>([^<]+)<\/Error2Radius>/i);
  const ra  = c1m ? parseFloat(c1m[1]) : null;
  const dec = c2m ? parseFloat(c2m[1]) : null;
  const err = errm ? parseFloat(errm[1]) : null;
  const trigger_id = param("Burst_Id") || tag("ivorn")?.split("#")[1] || null;
  const time = tag("ISOTime") || null;
  return {
    ra:  Number.isFinite(ra)  ? ra  : null,
    dec: Number.isFinite(dec) ? dec : null,
    err: Number.isFinite(err) ? err : null,
    trigger_id,
    time,
  };
}

function extractFromJson(j) {
  const ra  = Number.isFinite(j.ra)  ? j.ra  : null;
  const dec = Number.isFinite(j.dec) ? j.dec : null;
  const err =
    Number.isFinite(j.ra_err)      ? j.ra_err      :
    Number.isFinite(j.pos_err)     ? j.pos_err     :
    Number.isFinite(j.err_arcsec)  ? j.err_arcsec / 3600 :
    null;
  const trigger_id = j.trigger_id || j.id || j.event_id || null;
  const time = j.trigger_time || j.event_time || j.isotime || null;
  return { ra, dec, err, trigger_id, time };
}

function brokerFromTopic(topic) {
  if (topic.includes("SWIFT"))       return "Swift";
  if (topic.includes("FERMI"))       return "Fermi";
  if (topic.includes("ICECUBE") || topic.includes("AMON")) return "IceCube";
  if (topic.includes("LVC") || topic.includes("igwn"))     return "LVK";
  if (topic.includes("einstein_probe")) return "EP";
  if (topic.includes("svom"))        return "SVOM";
  if (topic.includes("INTEGRAL"))    return "INTEGRAL";
  if (topic.includes("AGILE"))       return "AGILE";
  if (topic.includes("GECAM"))       return "GECAM";
  if (topic.includes("MAXI"))        return "MAXI";
  if (topic.includes("IPN"))         return "IPN";
  if (topic.includes("SK_SN") || topic.includes("SNEWS")) return "SN-nu";
  return "GCN";
}

// ── Moon RA/Dec (low-precision) ─────────────────────────────────────────────
function moonRaDec(date = new Date()) {
  const jd = date.getTime() / 86400000 + 2440587.5;
  const T = (jd - 2451545.0) / 36525;
  const L = 218.316 + 481267.881 * T;
  const M = 134.963 + 477198.867 * T;
  const F =  93.272 + 483202.018 * T;
  const rad = Math.PI / 180;
  const lambda = L + 6.289 * Math.sin(M * rad);
  const beta   =       5.128 * Math.sin(F * rad);
  const eps = 23.43929;
  const sinLam = Math.sin(lambda * rad), cosLam = Math.cos(lambda * rad);
  const sinBet = Math.sin(beta   * rad), cosBet = Math.cos(beta   * rad);
  const sinEps = Math.sin(eps    * rad), cosEps = Math.cos(eps    * rad);
  const ra  = Math.atan2(sinLam * cosEps - (sinBet / cosBet) * sinEps, cosLam) / rad;
  const dec = Math.asin(sinBet * cosEps + cosBet * sinEps * sinLam) / rad;
  return { ra: ((ra % 360) + 360) % 360, dec };
}

function angularSep(ra1, dec1, ra2, dec2) {
  const rad = Math.PI / 180;
  const a1 = ra1 * rad, a2 = ra2 * rad;
  const d1 = dec1 * rad, d2 = dec2 * rad;
  const cosSep =
    Math.sin(d1) * Math.sin(d2) +
    Math.cos(d1) * Math.cos(d2) * Math.cos(a1 - a2);
  return Math.acos(Math.max(-1, Math.min(1, cosSep))) / rad;
}

// ── AstroColibri enrichment ─────────────────────────────────────────────────
async function enrichWithColibri(triggerId) {
  if (!triggerId) return null;
  try {
    const ctrl = new AbortController();
    const to = setTimeout(() => ctrl.abort(), 5000);
    const r = await fetch(
      `https://astro-colibri.science/event?trigger_id=${encodeURIComponent(triggerId)}`,
      { signal: ctrl.signal }
    );
    clearTimeout(to);
    if (!r.ok) return null;
    const j = await r.json();
    return j && j.astro_colibri_id ? j.astro_colibri_id : null;
  } catch (_) {
    return null;
  }
}

// ── Decision engine ─────────────────────────────────────────────────────────
function decide({ ra, dec, err, minAltDeg, minMoonSepDeg, maxErrDeg, getMountSite, computeAltAz }) {
  if (!Number.isFinite(ra) || !Number.isFinite(dec)) {
    return { ok: false, reason: "no-coords" };
  }
  const site = getMountSite();
  if (!site) {
    return { ok: true, alt: null, moonSep: null, note: "no-site-yet" };
  }
  const { alt } = computeAltAz(ra, dec, site.lat, site.lon, new Date());
  if (!Number.isFinite(alt) || alt < minAltDeg) {
    return { ok: false, reason: `alt-too-low:${(alt || 0).toFixed(1)}deg`, alt };
  }
  const moon = moonRaDec();
  const moonSep = angularSep(ra, dec, moon.ra, moon.dec);
  if (moonSep < minMoonSepDeg) {
    return { ok: false, reason: `moon-too-close:${moonSep.toFixed(1)}deg`, alt, moonSep };
  }
  const errLimit = Number.isFinite(maxErrDeg) ? maxErrDeg : 1;
  if (Number.isFinite(err) && err > errLimit) {
    return { ok: false, reason: `err-too-large:${err.toFixed(2)}deg(limit=${errLimit})`, alt, moonSep };
  }
  return { ok: true, alt, moonSep };
}

// ── Main entry point ────────────────────────────────────────────────────────
function startGcnListener(opts) {
  const {
    db, pushToQueue, saveQueue, computeAltAz, getMountSite,
    getStrategy, logger = console, minAltDeg = 25, minMoonSepDeg = 20,
  } = opts;

  const clientId     = process.env.GCN_CLIENT_ID;
  const clientSecret = process.env.GCN_CLIENT_SECRET;
  if (!clientId || !clientSecret) {
    logger.log("[alerts] GCN_CLIENT_ID / GCN_CLIENT_SECRET not set — listener disabled.");
    return { started: false };
  }

  initAlertsTable(db);

  const kafka = new Kafka({ client_id: clientId, client_secret: clientSecret });
  const consumer = kafka.consumer();

  const insertStmt = db.prepare(`
    INSERT INTO alerts
      (received_at, event_time, broker, topic, trigger_id, classification,
       ra, dec, err_deg, alt_now, moon_sep, action, action_reason, colibri_id, raw)
    VALUES (datetime('now'), @event_time, @broker, @topic, @trigger_id, @classification,
       @ra, @dec, @err_deg, @alt_now, @moon_sep, @action, @action_reason, @colibri_id, @raw)
  `);

  async function handleMessage(topic, value) {
    const raw = value?.toString?.() || "";
    const broker = brokerFromTopic(topic);
    let parsed = {};
    let fields = { ra: null, dec: null, err: null, trigger_id: null, time: null };
    const isVoevent = topic.includes("voevent");
    const isJson    = !isVoevent && topic.includes("gcn.notices.");

    try {
      if (isVoevent) {
        parsed = { raw_xml: true };
        fields = extractFromVoevent(raw);
      } else if (isJson) {
        parsed = JSON.parse(raw);
        fields = extractFromJson(parsed);
      } else {
        parsed = parseTextNotice(raw);
        fields = extractFromText(parsed);
      }
    } catch (e) {
      logger.warn(`[alerts] parse error topic=${topic}: ${e.message}`);
    }

    // ── Look up per-broker strategy ─────────────────────────────────────────
    const strategy = getStrategy ? getStrategy(broker) : null;
    if (strategy && (!strategy.enabled || strategy.mode === "ignore")) {
      logger.log(`[alerts] ${broker} ${fields.trigger_id || "?"} → dropped (strategy: ${strategy.mode}, enabled=${strategy.enabled})`);
      insertStmt.run({
        event_time: fields.time || null, broker, topic,
        trigger_id: fields.trigger_id ? String(fields.trigger_id) : null,
      classification: parsed.NOTICE_TYPE || parsed.classification || (parsed.raw_xml ? "VOEvent" : null),
      ra: Number.isFinite(fields.ra) ? fields.ra : null,
      dec: Number.isFinite(fields.dec) ? fields.dec : null,
      err_deg: Number.isFinite(fields.err) ? fields.err : null,
      alt_now: null, moon_sep: null,
      action: "ignored", action_reason: `strategy:${strategy.mode}`,
        colibri_id: null, raw,
      });
      return;
    }

    const stratMinAlt    = strategy?.min_alt      ?? minAltDeg;
    const stratMaxErrDeg = strategy?.max_err_deg  ?? 1;

    const decision = decide({
      ra: fields.ra, dec: fields.dec, err: fields.err,
      minAltDeg: stratMinAlt, minMoonSepDeg, maxErrDeg: stratMaxErrDeg,
      getMountSite, computeAltAz,
    });

    // "no-coords" means the notice has no position (e.g. GRM spectrometer) — distinct from "rejected by criteria"
    let action = decision.reason === "no-coords" ? "no-position" : "rejected";
    let action_reason = decision.reason || null;
    const stratMode = strategy?.mode || "too";

    if (decision.ok) {
      try {
        pushToQueue({
          name: `${broker}-${fields.trigger_id || "alert"}`,
          raDeg: fields.ra,
          decDeg: fields.dec,
          alert: {
            broker,
            trigger_id: fields.trigger_id,
            ra: fields.ra,
            dec: fields.dec,
            err_deg: fields.err,
            alt: decision.alt,
            moonSep: decision.moonSep,
            strategy,
          },
          mode: stratMode,
        });
        if (saveQueue) saveQueue();
        action = stratMode === "too" ? "queued" : "queued-only";
        action_reason = decision.note || `strategy:${stratMode}`;
      } catch (e) {
        action = "rejected";
        action_reason = `enqueue-error:${e.message}`;
      }
    }

    const colibri_id = decision.ok ? await enrichWithColibri(fields.trigger_id) : null;

    insertStmt.run({
      event_time:     fields.time || null,
      broker,
      topic,
      trigger_id:     fields.trigger_id ? String(fields.trigger_id) : null,
      classification: parsed.NOTICE_TYPE || parsed.classification || (parsed.raw_xml ? "VOEvent" : null),
      ra:             Number.isFinite(fields.ra)  ? fields.ra  : null,
      dec:            Number.isFinite(fields.dec) ? fields.dec : null,
      err_deg:        Number.isFinite(fields.err) ? fields.err : null,
      alt_now:        Number.isFinite(decision.alt)     ? decision.alt     : null,
      moon_sep:       Number.isFinite(decision.moonSep) ? decision.moonSep : null,
      action,
      action_reason,
      colibri_id,
      raw,
    });

    logger.log(
      `[alerts] ${broker} ${fields.trigger_id || "?"} ` +
      `RA=${fields.ra?.toFixed?.(3) ?? "?"} Dec=${fields.dec?.toFixed?.(3) ?? "?"} ` +
      `alt=${decision.alt?.toFixed?.(1) ?? "?"}° moonSep=${decision.moonSep?.toFixed?.(1) ?? "?"}° → ${action}${action_reason ? " ("+action_reason+")" : ""}`
    );
  }

  const INITIAL_BACKOFF_MS = 5_000;
  const MAX_BACKOFF_MS     = 5 * 60_000;

  async function connectAndRun() {
    try {
      await consumer.subscribe({ topics: TOPICS });
    } catch (err) {
      if (err?.type === "TOPIC_AUTHORIZATION_FAILED") {
        logger.warn("[alerts] some GCN topics are not accessible — continuing with the others.");
      } else {
        throw err;
      }
    }
    await consumer.run({
      eachMessage: async ({ topic, message }) => {
        await handleMessage(topic, message.value);
      },
    });
    logger.log(`[alerts] GCN-Kafka listener started — ${TOPICS.length} topics subscribed.`);
  }

  (async () => {
    let backoff = INITIAL_BACKOFF_MS;
    let attempt = 0;

    // eslint-disable-next-line no-constant-condition
    while (true) {
      try {
        await connectAndRun();
        backoff = INITIAL_BACKOFF_MS;
        attempt = 0;
        await new Promise(() => {});
      } catch (err) {
        attempt++;
        logger.error(`[alerts] consumer crashed (attempt ${attempt}), restarting in ${(backoff / 1000).toFixed(0)}s: ${err?.message || err}`);
        await new Promise((r) => setTimeout(r, backoff));
        backoff = Math.min(backoff * 2, MAX_BACKOFF_MS);

        try { await consumer.disconnect(); } catch (_) {}
      }
    }
  })();

  return { started: true, consumer };
}

function listAlerts(db, { limit = 100, offset = 0 } = {}) {
  return db
    .prepare(
      `SELECT id, received_at, event_time, broker, topic, trigger_id, classification,
              ra, dec, err_deg, alt_now, moon_sep, action, action_reason, colibri_id
       FROM alerts ORDER BY id DESC LIMIT ? OFFSET ?`
    )
    .all(limit, offset);
}

function getAlert(db, id) {
  return db.prepare(`SELECT * FROM alerts WHERE id = ?`).get(id);
}

module.exports = {
  TOPICS, initAlertsTable, startGcnListener, listAlerts, getAlert,
  parseTextNotice, extractFromText, extractFromJson, extractFromVoevent,
  angularSep, moonRaDec,
};
