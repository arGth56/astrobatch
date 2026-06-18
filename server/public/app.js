const form = document.getElementById("config-form");
const connectionStatus = document.getElementById("connection-status");
const logEl = document.getElementById("log");
const refreshButton = document.getElementById("refresh-status");
const connectButtons = Array.from(document.querySelectorAll("button[data-device]"));

const tabDevicesBtn   = document.getElementById("tab-devices");
const tabActionsBtn   = document.getElementById("tab-actions");
const tabDomeBtn      = document.getElementById("tab-dome");
const tabTargetBtn    = document.getElementById("tab-target");
const tabTodoBtn      = document.getElementById("tab-todo");
const tabAlertsBtn    = document.getElementById("tab-alerts");
const tabPipelineBtn  = document.getElementById("tab-pipeline");
const tabHistoryBtn   = document.getElementById("tab-history");
const tabSettingsBtn  = document.getElementById("tab-settings");
const panelDevices    = document.getElementById("panel-devices");
const panelActions    = document.getElementById("panel-actions");
const panelDome       = document.getElementById("panel-dome");
const panelTarget     = document.getElementById("panel-target");
const panelTodo       = document.getElementById("panel-todo");
const panelAlerts     = document.getElementById("panel-alerts");
const panelPipeline   = document.getElementById("panel-pipeline");
const panelHistory    = document.getElementById("panel-history");
const panelSettings   = document.getElementById("panel-settings");

const cameraCaptureForm = document.getElementById("camera-capture-form");
const filterwheelForm = document.getElementById("filterwheel-form");
const focuserInBtn = document.getElementById("focuser-in");
const focuserOutBtn = document.getElementById("focuser-out");
const mountParkBtn = document.getElementById("mount-park");
const mountUnparkBtn = document.getElementById("mount-unpark");
const mountHomeBtn = document.getElementById("mount-home");
const filterSelect = document.getElementById("filterwheel-filter-id");

// NINA-sourced devices
const deviceStatusEls = {
  mount: document.getElementById("status-mount"),
  camera: document.getElementById("status-camera"),
  filterwheel: document.getElementById("status-filterwheel"),
  focuser: document.getElementById("status-focuser"),
  guider: document.getElementById("status-guider"),
  rotator: document.getElementById("status-rotator"),
  flatdevice: document.getElementById("status-flatdevice"),
  switch: document.getElementById("status-switch"),
  weather: document.getElementById("status-weather"),
};

// OCS-sourced device rows (Dome + Safety Monitor)
const ocsDomeEl = document.getElementById("status-dome");
const ocsDomeDetailEl = document.getElementById("detail-dome");
const ocsSafetyEl = document.getElementById("status-safetymonitor");
const ocsSafetyDetailEl = document.getElementById("detail-safetymonitor");

let lastEquipment = null;
let ocsHost = "192.168.1.220";
const NINA_CONFIG_STORAGE_KEY = "ninaConfig";
let ocsConnected = false;
let stdwebConnected = false;
let observerLocation = { lat: null, lon: null, elevation: null };
let lastTargetCoords = { raDeg: null, decDeg: null };
let lastTnsTarget = null;

// Plan Tonight state
let _tonightEvents = [];
let _tonightSort   = { col: "maxAlt", dir: -1 }; // default: highest altitude first

// Sequence polling
let seqPollInterval    = null;
let seqLogRenderedCount = 0;

// Drag-and-drop state (module-level so re-renders during polling don't lose it)
let _seqDragSrc        = null;
let _seqDragging       = false;
let _seqDragInsertAfter = false;

// Track which per-target override panels are currently open (by queue index)
const _openOverridePanels = new Set();

function currentConfig() {
  return {
    host: document.getElementById("host").value.trim(),
    port: document.getElementById("port").value.trim(),
    protocol: document.getElementById("protocol").value,
  };
}

function saveNinaConfigToStorage() {
  try {
    localStorage.setItem(NINA_CONFIG_STORAGE_KEY, JSON.stringify(currentConfig()));
  } catch { /* private browsing / quota */ }
}

function loadNinaConfigFromStorage() {
  try {
    const raw = localStorage.getItem(NINA_CONFIG_STORAGE_KEY);
    return raw ? JSON.parse(raw) : null;
  } catch {
    return null;
  }
}

function applyNinaConfigFields(cfg) {
  if (!cfg?.host) return;
  const hostEl = document.getElementById("host");
  const portEl = document.getElementById("port");
  const protoEl = document.getElementById("protocol");
  if (hostEl) hostEl.value = cfg.host;
  if (portEl && cfg.port) portEl.value = cfg.port;
  if (protoEl && cfg.protocol) protoEl.value = cfg.protocol;
}

function applyNinaConfigFromStorage() {
  applyNinaConfigFields(loadNinaConfigFromStorage());
}

function setLog(value) {
  logEl.textContent = typeof value === "string" ? value : JSON.stringify(value, null, 2);
}

function setActiveTab(tab) {
  const tabs = { devices: false, actions: false, dome: false, target: false, todo: false, alerts: false, pipeline: false, history: false, settings: false };
  tabs[tab] = true;
  tabDevicesBtn.classList.toggle("active",  tabs.devices);
  tabActionsBtn.classList.toggle("active",  tabs.actions);
  tabDomeBtn.classList.toggle("active",     tabs.dome);
  tabTargetBtn.classList.toggle("active",   tabs.target);
  tabTodoBtn.classList.toggle("active",     tabs.todo);
  tabAlertsBtn.classList.toggle("active",   tabs.alerts);
  tabPipelineBtn.classList.toggle("active", tabs.pipeline);
  tabHistoryBtn.classList.toggle("active",  tabs.history);
  tabSettingsBtn.classList.toggle("active", tabs.settings);
  panelDevices.classList.toggle("active",   tabs.devices);
  panelActions.classList.toggle("active",   tabs.actions);
  panelDome.classList.toggle("active",      tabs.dome);
  panelTarget.classList.toggle("active",    tabs.target);
  panelTodo.classList.toggle("active",      tabs.todo);
  panelAlerts.classList.toggle("active",    tabs.alerts);
  panelPipeline.classList.toggle("active",  tabs.pipeline);
  panelHistory.classList.toggle("active",   tabs.history);
  panelSettings.classList.toggle("active",  tabs.settings);

  if (tab === "todo") {
    startSeqPolling();
    initNightPlan();
    refreshNightPlan();
  } else {
    stopSeqPolling();
  }
  if (tab === "actions") {
    // Refresh filter list when Actions tab opens
    refreshStatus().catch(() => {});
  }
  if (tab === "pipeline") {
    startPipelinePolling();
  } else {
    stopPipelinePolling();
  }
  if (tab !== "settings") {
    stopWatchdogPolling();
  }
}

function currentOcsHost() {
  // Primary source: top-level OCS host input; fallback to Dome tab input
  const top = document.getElementById("ocs-host-top");
  return (top?.value.trim()) || document.getElementById("dome-ocs-host").value.trim() || ocsHost;
}

function syncOcsHostInputs(host) {
  const top = document.getElementById("ocs-host-top");
  const dome = document.getElementById("dome-ocs-host");
  if (top) top.value = host;
  if (dome) dome.value = host;
  ocsHost = host;
}

function setOcsStatus(roof = {}, reachable = true) {
  // Dome tab detail view
  const roofEl = document.getElementById("dome-roof-status");
  const safeEl = document.getElementById("dome-safe");

  const roofStatus = (roof.status || "Unknown").replace(/<[^>]*>/g, "");
  roofEl.textContent = reachable ? roofStatus : "Unreachable";
  const isOpen = /open/i.test(roofStatus);
  const isClosed = /close/i.test(roofStatus);
  roofEl.className = reachable ? (isOpen ? "online" : isClosed ? "offline" : "") : "offline";

  const safe = roof.safe || "Unknown";
  safeEl.textContent = reachable ? safe : "Unreachable";
  safeEl.className = reachable && /safe/i.test(safe) ? "online" : "offline";

  document.getElementById("dome-rain").textContent = roof.rain || "–";
  document.getElementById("dome-temp").textContent = roof.temp || "–";
  document.getElementById("dome-humidity").textContent = roof.humidity || "–";
  document.getElementById("dome-pressure").textContent = roof.pressure || "–";

  // Devices tab — Dome row (roof open/closed)
  ocsDomeEl.textContent = reachable ? roofStatus : "Unreachable";
  ocsDomeEl.className = reachable ? (isOpen ? "online" : isClosed ? "offline" : "") : "offline";
  if (ocsDomeDetailEl) {
    ocsDomeDetailEl.textContent = roof.rain ? `Rain: ${roof.rain}` : "";
  }

  // Devices tab — Safety Monitor row
  ocsSafetyEl.textContent = reachable ? safe : "Unreachable";
  ocsSafetyEl.className = reachable && /safe/i.test(safe) ? "online" : "offline";
  if (ocsSafetyDetailEl) {
    const details = [];
    if (roof.temp && roof.temp !== "Invalid") details.push(`${roof.temp}°C`);
    if (roof.humidity && roof.humidity !== "Invalid") details.push(`${roof.humidity}% RH`);
    ocsSafetyDetailEl.textContent = details.join("  ");
  }

  // Update global OCS connection state
  ocsConnected = reachable;
  updateOcsIndicator();
  updateSeqRunButton();
}

function updateOcsIndicator() {
  const el = document.getElementById("ocs-connection-status");
  if (!el) return;
  const hostEl = document.getElementById("ocs-host-display");
  if (hostEl) hostEl.textContent = currentOcsHost();
  if (ocsConnected) {
    el.textContent = `OCS ✓ ${currentOcsHost()}`;
    el.className = "ocs-status-indicator online";
  } else {
    el.textContent = `OCS ✗ ${currentOcsHost()} — not reachable`;
    el.className = "ocs-status-indicator offline";
  }
}

function updateStdwebIndicator(url, reachable, detail) {
  const el     = document.getElementById("stdweb-connection-status");
  const hostEl = document.getElementById("stdweb-host-display");
  const shortUrl = url ? url.replace(/^https?:\/\//, "") : "—";
  if (hostEl) hostEl.textContent = shortUrl;
  if (!el) return;
  stdwebConnected = reachable;
  if (reachable) {
    el.textContent = `STDWeb ✓ ${shortUrl}`;
    el.className = "ocs-status-indicator online";
  } else {
    el.textContent = `STDWeb ✗ ${shortUrl}${detail ? ` — ${detail}` : " — not reachable"}`;
    el.className = "ocs-status-indicator offline";
  }
}

async function checkStdwebHealth() {
  const el = document.getElementById("stdweb-connection-status");
  if (el) { el.textContent = "STDWeb — checking…"; el.className = "ocs-status-indicator"; }
  try {
    const data = await fetch("/api/stdweb/health").then(r => r.json());
    updateStdwebIndicator(data.url, data.reachable, data.error || (data.httpStatus && !data.reachable ? `HTTP ${data.httpStatus}` : null));
  } catch (err) {
    updateStdwebIndicator(null, false, err.message);
  }
}

function updateSeqRunButton() {
  const runBtn = document.getElementById("seq-run");
  if (!runBtn) return;
  const seqRunning = runBtn.disabled && document.getElementById("seq-abort") && !document.getElementById("seq-abort").disabled;
  if (seqRunning) return; // sequence already running — don't touch it
  if (!ocsConnected) {
    runBtn.disabled = true;
    runBtn.title = "OCS not connected — cannot launch sequence";
  } else {
    runBtn.disabled = false;
    runBtn.title = "";
  }
}

async function refreshOcsStatus() {
  const host = currentOcsHost();
  setLog(`Fetching OCS status from ${host}...`);
  try {
    const payload = await postJson("/api/ocs/status", { ocsHost: host });
    setOcsStatus(payload.roof || {}, true);
    setLog(payload);
  } catch (error) {
    setOcsStatus({}, false);
    setLog({ success: false, action: "ocs-status", error: error.message });
  }
}

function currentTnsCredentials() {
  return {
    botId: document.getElementById("tns-bot-id").value.trim(),
    botName: document.getElementById("tns-bot-name").value.trim(),
    apiKey: document.getElementById("tns-api-key").value.trim(),
  };
}

function setTargetResult(target) {
  const el = (id, val) => {
    const e = document.getElementById(id);
    if (e) e.textContent = val ?? "–";
  };
  // Name + source badge
  const nameBadge = target._source
    ? ` <span style="font-size:11px;background:#1e3a5f;color:#93c5fd;padding:2px 6px;border-radius:4px;font-weight:normal;">${target._source}</span>`
    : "";
  document.getElementById("target-full-name").innerHTML =
    (target.fullName || target.name) + nameBadge;

  // Primary link: prefer TNS, add AstroColibri if available
  const link = document.getElementById("target-tns-link");
  link.href = target.tnsUrl || "#";
  link.textContent = "View on TNS ↗";
  link.style.display = target.tnsUrl ? "" : "none";

  // AstroColibri link (injected after TNS link if available)
  let acLink = document.getElementById("target-colibri-link");
  if (!acLink) {
    acLink = document.createElement("a");
    acLink.id = "target-colibri-link";
    acLink.target = "_blank";
    acLink.className = "tns-link";
    acLink.style.marginLeft = "10px";
    link.parentNode.insertBefore(acLink, link.nextSibling);
  }
  if (target.colibriUrl) {
    acLink.href = target.colibriUrl;
    acLink.textContent = "View on AstroColibri ↗";
    acLink.style.display = "";
  } else {
    acLink.style.display = "none";
  }

  el("t-type", target.type);
  el("t-ra", target.ra || null);
  el("t-dec", target.dec || null);
  el("t-redshift", target.redshift != null ? target.redshift : null);
  el("t-disc-date", target.discoveryDate);
  el("t-disc-mag", target.discoveryMag != null ? target.discoveryMag : null);
  el("t-disc-filter", target.discoveryFilter);
  el("t-host", target.hostName);
  el("t-host-z", target.hostRedshift != null ? target.hostRedshift : null);
  el("t-reporters", target.reporters);
  el("t-group", target.reportingGroup);
  el("t-internal", target.internalNames);
  el("t-constellation", target.constellation);

  // Store for "Add to ToDo"
  lastTnsTarget = target;
  const addBtn = document.getElementById("add-to-todo");
  if (addBtn) {
    addBtn.style.display = target.raDeg != null ? "inline-block" : "none";
  }

  document.getElementById("target-result").style.display = "block";

  // Draw sky charts if we have coordinates
  if (target.raDeg != null && target.decDeg != null) {
    lastTargetCoords = { raDeg: target.raDeg, decDeg: target.decDeg };
    document.getElementById("target-visibility").style.display = "block";
    updateObsLocationDisplay();
    drawVisibility();
  }
}

function updateObsLocationDisplay() {
  const textEl = document.getElementById("obs-location-text");
  if (!textEl) return;
  if (observerLocation.lat != null && observerLocation.lon != null) {
    const lat = observerLocation.lat.toFixed(4);
    const lon = observerLocation.lon.toFixed(4);
    const elev = observerLocation.elevation != null ? ` · ${Math.round(observerLocation.elevation)} m` : "";
    textEl.textContent = `${lat}°, ${lon}°${elev}`;
    // Fill the form fields as well
    const latInput = document.getElementById("obs-lat");
    const lonInput = document.getElementById("obs-lon");
    const elevInput = document.getElementById("obs-elev");
    if (latInput && !latInput.value) latInput.value = lat;
    if (lonInput && !lonInput.value) lonInput.value = lon;
    if (elevInput && !elevInput.value) elevInput.value = Math.round(observerLocation.elevation ?? 0);
  } else {
    textEl.textContent = "Unknown — enter manually";
  }
}

function drawVisibility() {
  const { raDeg, decDeg } = lastTargetCoords;
  if (raDeg == null || decDeg == null) return;

  // Use stored location or fall back to form inputs
  let lat = observerLocation.lat;
  let lon = observerLocation.lon;
  let elevation = observerLocation.elevation ?? 0;

  const latInput = document.getElementById("obs-lat");
  const lonInput = document.getElementById("obs-lon");
  const elevInput = document.getElementById("obs-elev");
  if (latInput?.value) lat = parseFloat(latInput.value);
  if (lonInput?.value) lon = parseFloat(lonInput.value);
  if (elevInput?.value) elevation = parseFloat(elevInput.value);

  if (lat == null || lon == null || isNaN(lat) || isNaN(lon)) {
    document.getElementById("obs-location-text").textContent =
      "Unknown — enter Lat/Lon above";
    return;
  }

  const now = new Date();
  const HOURS = 12;
  const STEP = 15; // minutes

  const targetPoints = Astronomy.computeCurve(raDeg, decDeg, lat, lon, now, HOURS, STEP);
  const sunPoints = Astronomy.computeSunCurve(lat, lon, now, HOURS, STEP);
  const { alt, az } = Astronomy.computeAltAz(raDeg, decDeg, lat, lon, now);

  // Sky dome
  const dome = document.getElementById("sky-dome");
  Astronomy.drawSkyDome(dome, alt, az);

  // Altitude chart — size canvas to its container width
  const chart = document.getElementById("alt-chart");
  const containerWidth = chart.parentElement.clientWidth;
  if (containerWidth > 0) chart.width = containerWidth;
  Astronomy.drawAltChart(chart, targetPoints, sunPoints);
}

function saveToHistory(target, source) {
  if (target.raDeg == null) return;
  postJson("/api/targets", {
    name:    target.fullName || target.name,
    ra:      target.ra      || null,
    dec:     target.dec     || null,
    ra_deg:  target.raDeg,
    dec_deg: target.decDeg,
    type:    target.type    || null,
    source,
  }).then(loadHistory).catch(() => {});
}

async function lookupTarget() {
  const name = document.getElementById("target-name").value.trim();
  if (!name) return;

  const { botId, botName, apiKey } = currentTnsCredentials();
  document.getElementById("target-result").style.display = "none";

  // ── Step 1: TNS ────────────────────────────────────────────────────────────
  setLog(`Looking up ${name} on TNS…`);
  try {
    const payload = await postJson("/api/target/tns", { name, botId, botName, apiKey });
    if (payload.success) {
      setTargetResult({ ...payload.target, _source: "TNS" });
      saveToHistory(payload.target, "tns");
      setLog(payload);
      return;
    }
    // TNS failed — try AstroColibri
    setLog({ ...payload, _note: "TNS did not find it — trying AstroColibri…" });
  } catch (err) {
    setLog({ success: false, error: err.message, _note: "TNS error — trying AstroColibri…" });
  }

  // ── Step 2: AstroColibri fallback ─────────────────────────────────────────
  setLog(`TNS: no result — searching AstroColibri for ${name}…`);
  try {
    const res = await fetch(`/api/target/astrocolibri?name=${encodeURIComponent(name)}`);
    const payload = await res.json();
    if (payload.success) {
      setTargetResult({ ...payload.target, _source: "AstroColibri" });
      saveToHistory(payload.target, "astrocolibri");
      setLog({ ...payload, _note: "Found on AstroColibri (TNS had no result)" });
    } else {
      setLog({ success: false, error: `Not found on TNS or AstroColibri: ${payload.error}` });
    }
  } catch (err) {
    setLog({ success: false, error: `AstroColibri error: ${err.message}` });
  }
}

async function roofCommand(command) {
  const host = currentOcsHost();
  setLog(`Sending roof command: ${command}...`);
  try {
    const resp = await fetch("/api/ocs/roof", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ ocsHost: host, command }),
    });
    const payload = await resp.json();
    if (resp.status === 403) {
      alert(`⛔ ${payload.error}`);
      setLog(payload);
      return;
    }
    setOcsStatus(payload.roof || {});
    setLog(payload);
  } catch (error) {
    setLog({ success: false, action: "roof", command, error: error.message });
  }
}

function setDeviceStatuses(devices = {}) {
  for (const [device, el] of Object.entries(deviceStatusEls)) {
    const connected = devices?.[device];
    if (typeof connected !== "boolean") {
      el.textContent = "Unknown";
      el.className = "offline";
      continue;
    }
    el.textContent = connected ? "Connected" : "Disconnected";
    el.className = connected ? "online" : "offline";
  }
}

function updateFilterOptions(equipment = null) {
  const raw = equipment?.FilterWheel?.AvailableFilters || [];
  // Keep only slots with a non-empty name; fallback to 4 numeric slots so the
  // dropdown is always usable even when NINA is offline.
  const named   = raw.filter(f => f.Name && f.Name.trim());
  const filters  = named.length ? named : raw.length ? raw : null;
  const selectedId = equipment?.FilterWheel?.SelectedFilter?.Id;

  filterSelect.innerHTML = "";

  if (filters) {
    for (const filter of filters) {
      const option = document.createElement("option");
      option.value = String(filter.Id);
      const label = filter.Name && filter.Name.trim() ? filter.Name : `Filter ${filter.Id}`;
      option.textContent = `${label} (slot ${filter.Id})`;
      if (selectedId === filter.Id) option.selected = true;
      filterSelect.appendChild(option);
    }
  } else {
    // NINA unreachable — offer 4 numeric slots as fallback
    for (let i = 0; i < 4; i++) {
      const option = document.createElement("option");
      option.value = String(i);
      option.textContent = `Slot ${i}`;
      filterSelect.appendChild(option);
    }
  }
}

async function postJson(url, data) {
  const response = await fetch(url, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(data),
  });

  let payload;
  try {
    payload = await response.json();
  } catch {
    payload = { success: false, error: "Invalid JSON response from server" };
  }

  if (!response.ok) {
    throw new Error(payload?.error || `Request failed with ${response.status}`);
  }
  return payload;
}

async function refreshStatus() {
  const [ninaResult, ocsResult] = await Promise.allSettled([
    postJson("/api/nina/devices/status", currentConfig()),
    postJson("/api/ocs/status", { ocsHost: currentOcsHost() }),
  ]);

  // Update NINA host display
  const ninaHostEl = document.getElementById("nina-host-display");
  if (ninaHostEl) {
    const cfg = currentConfig();
    ninaHostEl.textContent = `${cfg.host}:${cfg.port}`;
  }

  if (ninaResult.status === "fulfilled") {
    const payload = ninaResult.value;
    connectionStatus.textContent = payload.success
      ? `Connected to ${payload.target.protocol}://${payload.target.host}:${payload.target.port}`
      : `Connection failed (${payload.status || "unknown status"})`;
    connectionStatus.classList.toggle("online", !!payload.success);
    connectionStatus.classList.toggle("offline", !payload.success);
    setDeviceStatuses(payload.devices || {});
    lastEquipment = payload.equipment || null;
    updateFilterOptions(lastEquipment);
    updateFwCurrentBadge(lastEquipment);
    updateCoverBadge(lastEquipment?.FlatDevice ?? null, payload.devices?.flatdevice);

    // Store site coordinates when mount is connected and has real coordinates
    const mount = lastEquipment?.Mount;
    if (mount?.Connected && (mount.SiteLatitude !== 0 || mount.SiteLongitude !== 0)) {
      observerLocation = {
        lat: mount.SiteLatitude,
        lon: mount.SiteLongitude,
        elevation: mount.SiteElevation ?? 0,
      };
      updateObsLocationDisplay();
    }

    setLog(payload);
  } else {
    connectionStatus.textContent = `Connection failed: ${ninaResult.reason?.message}`;
    connectionStatus.classList.add("offline");
    connectionStatus.classList.remove("online");
    lastEquipment = null;
    setDeviceStatuses({});
    updateCoverBadge(null);
    setLog({ success: false, error: ninaResult.reason?.message });
  }

  if (ocsResult.status === "fulfilled") {
    setOcsStatus(ocsResult.value.roof || {}, true);
  } else {
    setOcsStatus({}, false);
  }
}

async function connectDevice(device) {
  setLog(`Connecting ${device}...`);
  try {
    const payload = await postJson("/api/nina/devices/connect", {
      ...currentConfig(),
      device,
    });
    setDeviceStatuses(payload.devices || {});
    await refreshStatus();
    setLog(payload);
  } catch (error) {
    setLog({ success: false, device, error: error.message });
  }
}

async function runMountCommand(command) {
  setLog(`Mount ${command}...`);
  try {
    const payload = await postJson("/api/nina/actions/mount", {
      ...currentConfig(),
      command,
    });
    setDeviceStatuses(payload.devices || {});
    await refreshStatus();
    setLog(payload);
  } catch (error) {
    setLog({ success: false, command, error: error.message });
  }
}

async function runCameraCapture() {
  const duration = Number(document.getElementById("camera-duration").value);
  const gain = Number(document.getElementById("camera-gain").value);
  const savePath = document.getElementById("camera-save-path").value.trim();
  setLog("Triggering capture...");
  try {
    const payload = await postJson("/api/nina/actions/camera/capture", {
      ...currentConfig(),
      duration,
      gain,
      savePath,
    });
    await refreshStatus();
    setLog(payload);
  } catch (error) {
    setLog({ success: false, action: "capture", error: error.message });
  }
}

function updateFwCurrentBadge(equipment) {
  const fw = equipment?.FilterWheel;
  const label = document.getElementById("fw-current-label");
  const slot  = document.getElementById("fw-current-slot");
  if (!label || !slot) return;
  if (fw?.Connected && fw.SelectedFilter) {
    label.textContent = fw.SelectedFilter.Name || `Filter ${fw.SelectedFilter.Id}`;
    label.style.color = "#a78bfa";
    slot.textContent  = `(slot ${fw.SelectedFilter.Id})`;
  } else {
    label.textContent = fw?.Connected ? "Unknown" : "Not connected";
    label.style.color = "#9ca3af";
    slot.textContent  = "";
  }
}

async function changeFilter() {
  const filterId    = Number(filterSelect.value);
  const filterName  = filterSelect.options[filterSelect.selectedIndex]?.text || String(filterId);
  const msgEl       = document.getElementById("fw-result-msg");

  // Capture before state
  const beforeFilter = lastEquipment?.FilterWheel?.SelectedFilter;
  const beforeName   = beforeFilter?.Name || (beforeFilter ? `Filter ${beforeFilter.Id}` : "Unknown");
  const beforeSlot   = beforeFilter?.Id ?? "?";

  if (msgEl) msgEl.innerHTML = `<span style="color:#6b7280">Sending change to ${filterName}…</span>`;
  setLog(`Changing filter to ${filterName} (id=${filterId})...`);

  try {
    const payload = await postJson("/api/nina/actions/filterwheel/change", {
      ...currentConfig(),
      filterId,
    });
    await refreshStatus();  // updates lastEquipment

    // Capture after state
    const afterFilter = lastEquipment?.FilterWheel?.SelectedFilter;
    const afterName   = afterFilter?.Name || (afterFilter ? `Filter ${afterFilter.Id}` : "Unknown");
    const afterSlot   = afterFilter?.Id ?? "?";

    const moved = afterSlot !== beforeSlot;
    if (msgEl) {
      if (moved) {
        msgEl.innerHTML =
          `<span style="color:#34d399">✔ Moved: <b>${beforeName}</b> (slot ${beforeSlot}) → <b>${afterName}</b> (slot ${afterSlot})</span>`;
      } else {
        msgEl.innerHTML =
          `<span style="color:#fbbf24">⚠ Filter wheel position unchanged — still <b>${afterName}</b> (slot ${afterSlot}). Was it already on this filter, or did it not move?</span>`;
      }
    }
    setLog(payload);
  } catch (error) {
    if (msgEl) msgEl.innerHTML = `<span style="color:#f87171">✘ Error: ${error.message}</span>`;
    setLog({ success: false, action: "filter-change", error: error.message });
  }
}

async function verifyFilterPosition() {
  const msgEl = document.getElementById("fw-result-msg");
  if (msgEl) msgEl.innerHTML = `<span style="color:#6b7280">Querying NINA…</span>`;
  await refreshStatus();
  const fw = lastEquipment?.FilterWheel;
  if (msgEl) {
    if (fw?.Connected && fw.SelectedFilter) {
      const filters = fw.AvailableFilters || [];
      const list = filters.map(f => `${f.Name}=slot${f.Id}`).join(", ");
      msgEl.innerHTML =
        `<span style="color:#60a5fa">Current: <b>${fw.SelectedFilter.Name}</b> (slot ${fw.SelectedFilter.Id})</span>` +
        (list ? `<br><span style="font-size:11px;color:#6b7280">Available: ${list}</span>` : "");
    } else {
      msgEl.innerHTML = `<span style="color:#f87171">Filter wheel not connected</span>`;
    }
  }
}

async function moveFocuser(direction) {
  const steps = Number(document.getElementById("focuser-steps").value);
  setLog(`Moving focuser ${direction} by ${steps}...`);
  try {
    const payload = await postJson("/api/nina/actions/focuser/move-relative", {
      ...currentConfig(),
      direction,
      steps,
    });
    await refreshStatus();
    setLog(payload);
  } catch (error) {
    setLog({ success: false, action: "focuser-move", direction, error: error.message });
  }
}

form.addEventListener("submit", async (event) => {
  event.preventDefault();
  saveNinaConfigToStorage();
  // Update NINA host display in connect-bar
  const ninaHostEl = document.getElementById("nina-host-display");
  if (ninaHostEl) {
    const host = document.getElementById("host").value;
    const port = document.getElementById("port").value;
    ninaHostEl.textContent = `${host}:${port}`;
  }
  await refreshStatus();
});

// ── Connect-bar Test buttons ──────────────────────────────────────────────────

document.getElementById("btn-nina-test")?.addEventListener("click", async () => {
  const btn = document.getElementById("btn-nina-test");
  const el  = document.getElementById("connection-status");
  btn.disabled = true;
  btn.textContent = "…";
  if (el) { el.textContent = "Testing…"; el.className = "status"; }
  try {
    await refreshStatus();
  } finally {
    btn.disabled = false;
    btn.textContent = "Test";
  }
});

document.getElementById("btn-ocs-test")?.addEventListener("click", async () => {
  const btn = document.getElementById("btn-ocs-test");
  btn.disabled = true;
  btn.textContent = "…";
  try {
    await refreshOcsStatus();
  } finally {
    btn.disabled = false;
    btn.textContent = "Test";
  }
});

document.getElementById("btn-stdweb-test")?.addEventListener("click", async () => {
  const btn = document.getElementById("btn-stdweb-test");
  btn.disabled = true;
  btn.textContent = "…";
  try {
    await checkStdwebHealth();
  } finally {
    btn.disabled = false;
    btn.textContent = "Test";
  }
});

// ── STDWeb token save ─────────────────────────────────────────────────────────
async function loadStdwebTokenStatus() {
  try {
    const r = await fetch("/api/secrets/secret_stdweb_token");
    if (!r.ok) return;
    const data = await r.json();
    const el = document.getElementById("stdweb-token-status");
    if (el) el.textContent = data.set ? "Token configured ✓" : "No token saved yet";
    if (el) el.style.color = data.set ? "#86efac" : "#94a3b8";
  } catch {}
}

document.getElementById("stdweb-token-save")?.addEventListener("click", async () => {
  const input = document.getElementById("stdweb-token-input");
  const status = document.getElementById("stdweb-token-status");
  const btn = document.getElementById("stdweb-token-save");
  const value = input?.value.trim();
  if (!value) { if (status) { status.textContent = "Enter a token first"; status.style.color = "#fca5a5"; } return; }
  btn.disabled = true;
  btn.textContent = "Saving…";
  try {
    const r = await fetch("/api/secrets", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ key: "secret_stdweb_token", value }),
    });
    if (!r.ok) throw new Error((await r.json()).error || r.statusText);
    if (input) input.value = "";
    if (status) { status.textContent = "Token saved ✓"; status.style.color = "#86efac"; }
    await checkStdwebHealth();
  } catch (e) {
    if (status) { status.textContent = "Error: " + e.message; status.style.color = "#fca5a5"; }
  } finally {
    btn.disabled = false;
    btn.textContent = "Save Token";
  }
});

tabDevicesBtn.addEventListener("click",  () => setActiveTab("devices"));
tabActionsBtn.addEventListener("click",  () => setActiveTab("actions"));
tabDomeBtn.addEventListener("click",     () => { setActiveTab("dome"); refreshOcsStatus(); loadOcsHistory(); loadIrImage(); loadRainForecast(); });
tabTargetBtn.addEventListener("click",   () => setActiveTab("target"));

// ── Email notification settings ───────────────────────────────────────────────
async function loadEmailConfig() {
  const keys = ["email_host","email_port","email_secure","email_user","email_from","email_to"];
  try {
    const vals = await Promise.all(keys.map(k => fetch(`/api/settings/${k}`).then(r => r.ok ? r.json() : null)));
    const map = {};
    keys.forEach((k, i) => { if (vals[i]?.success) map[k] = vals[i].value; });
    if (map.email_host)   document.getElementById("email-host").value   = map.email_host;
    if (map.email_port)   document.getElementById("email-port").value   = map.email_port;
    if (map.email_secure) document.getElementById("email-secure").checked = map.email_secure === "true";
    if (map.email_user)   document.getElementById("email-user").value   = map.email_user;
    if (map.email_from)   document.getElementById("email-from").value   = map.email_from;
    if (map.email_to)     document.getElementById("email-to").value     = map.email_to;
    // show masked pass status
    const passRes = await fetch("/api/secrets/secret_email_pass");
    const passData = passRes.ok ? await passRes.json() : null;
    const st = document.getElementById("email-status");
    if (st && passData?.set) { st.textContent = "Password saved ✓"; st.style.color = "#86efac"; }
  } catch {}
}

document.getElementById("email-save")?.addEventListener("click", async () => {
  const btn = document.getElementById("email-save");
  const st  = document.getElementById("email-status");
  btn.disabled = true; btn.textContent = "Saving…";
  try {
    const settings = {
      email_host:   document.getElementById("email-host").value.trim(),
      email_port:   document.getElementById("email-port").value,
      email_secure: document.getElementById("email-secure").checked ? "true" : "false",
      email_user:   document.getElementById("email-user").value.trim(),
      email_from:   document.getElementById("email-from").value.trim(),
      email_to:     document.getElementById("email-to").value.trim(),
    };
    await Promise.all(Object.entries(settings).map(([key, value]) =>
      fetch("/api/settings", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ key, value }) })
    ));
    const pass = document.getElementById("email-pass").value;
    if (pass) {
      await fetch("/api/secrets", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ key: "secret_email_pass", value: pass }) });
      document.getElementById("email-pass").value = "";
    }
    if (st) { st.textContent = "Saved ✓"; st.style.color = "#86efac"; }
  } catch (e) {
    if (st) { st.textContent = "Error: " + e.message; st.style.color = "#fca5a5"; }
  } finally {
    btn.disabled = false; btn.textContent = "Save";
  }
});

document.getElementById("email-test")?.addEventListener("click", async () => {
  const btn = document.getElementById("email-test");
  const st  = document.getElementById("email-status");
  btn.disabled = true; btn.textContent = "Sending…";
  try {
    const r = await fetch("/api/email/test", { method: "POST" });
    const d = await r.json();
    if (d.success) { st.textContent = `Test sent to ${d.to} ✓`; st.style.color = "#86efac"; }
    else           { st.textContent = `Failed: ${d.error}`;      st.style.color = "#fca5a5"; }
  } catch (e) {
    st.textContent = "Error: " + e.message; st.style.color = "#fca5a5";
  } finally {
    btn.disabled = false; btn.textContent = "Send Test Email";
  }
});
tabAlertsBtn.addEventListener("click",   () => { setActiveTab("alerts"); loadAlerts(); loadAlertConfig(); loadStrategies(); });
tabSettingsBtn.addEventListener("click", () => {
  setActiveTab("settings");
  loadWatchdogConfig();
  loadArbiterConfig();
  initQualityReplayDefaults();
  startWatchdogPolling();
  loadCoverDisabled();
  loadStdwebTokenStatus();
  loadEmailConfig();
  loadMpcSettings();
});

// ── MPC / SkyBoT settings ─────────────────────────────────────────────────────
async function loadMpcSettings() {
  try {
    const res = await fetch("/api/settings/mpc_mag_limit");
    const d   = res.ok ? await res.json() : null;
    const inp = document.getElementById("mpc-mag-limit");
    if (d?.success && d.value != null) {
      _mpcMagLimit = parseFloat(d.value) || 21;
      if (inp) inp.value = _mpcMagLimit;
    }
  } catch {}
}

document.getElementById("mpc-settings-save")?.addEventListener("click", async () => {
  const btn = document.getElementById("mpc-settings-save");
  const st  = document.getElementById("mpc-settings-status");
  const val = document.getElementById("mpc-mag-limit")?.value;
  btn.disabled = true; btn.textContent = "Saving…";
  try {
    await fetch("/api/settings", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ key: "mpc_mag_limit", value: val }),
    });
    _mpcMagLimit = parseFloat(val) || 21;
    if (st) { st.textContent = "Saved ✓"; st.style.color = "#86efac"; }
  } catch (e) {
    if (st) { st.textContent = "Error: " + e.message; st.style.color = "#fca5a5"; }
  } finally {
    btn.disabled = false; btn.textContent = "Save";
    setTimeout(() => { if (st) st.textContent = ""; }, 3000);
  }
});

async function loadCoverDisabled() {
  try {
    const res = await fetch("/api/settings/coverDisabled");
    const data = await res.json();
    document.getElementById("set-cover-disabled").checked = data.value === "true";
  } catch { /* ignore */ }
}

document.getElementById("set-cover-disabled").addEventListener("change", async (e) => {
  try {
    await postJson("/api/settings", { key: "coverDisabled", value: e.target.checked ? "true" : "false" });
  } catch { /* ignore */ }
});

document.getElementById("target-lookup-form").addEventListener("submit", async (event) => {
  event.preventDefault();
  await lookupTarget();
});

document.getElementById("dome-ocs-config-form").addEventListener("submit", async (event) => {
  event.preventDefault();
  syncOcsHostInputs(document.getElementById("dome-ocs-host").value.trim() || ocsHost);
  await refreshOcsStatus();
});

document.getElementById("ocs-connect-form")?.addEventListener("submit", async (event) => {
  event.preventDefault();
  syncOcsHostInputs(document.getElementById("ocs-host-top").value.trim() || ocsHost);
  await refreshOcsStatus();
});

document.getElementById("dome-open").addEventListener("click", async () => {
  if (confirm("Open the roof?")) await roofCommand("open");
});
document.getElementById("dome-close").addEventListener("click", async () => {
  const mount = lastEquipment?.Mount;
  const fd    = lastEquipment?.FlatDevice;

  // Build list of blocking conditions
  const blocks = [];
  if (mount?.Connected && mount.AtPark !== true)
    blocks.push(`mount is not parked (AtPark=${mount.AtPark}${mount.Slewing ? ", slewing" : ""})`);
  if (fd?.Connected && normCoverStateClient(fd.CoverState) !== 1)
    blocks.push(`dust cover is not closed (${COVER_STATE_LABELS[normCoverStateClient(fd.CoverState)] ?? fd.CoverState})`);

  if (blocks.length) {
    alert(`⛔ Cannot close roof:\n\n• ${blocks.join("\n• ")}\n\nFix these first.`);
    return;
  }

  const ninaUnknown = !mount?.Connected;
  const coverUnknown = !fd?.Connected;
  let warning = "";
  if (ninaUnknown)   warning += "⚠ NINA not connected — cannot verify mount park state.\n";
  if (coverUnknown)  warning += "⚠ Cover device not connected — cannot verify dust cover state.\n";

  const msg = warning
    ? `${warning}\nClose roof anyway?`
    : "Mount parked ✓  Cover closed ✓\n\nClose the roof?";

  if (!confirm(msg)) return;
  await roofCommand("close");
});
document.getElementById("dome-stop").addEventListener("click", async () => roofCommand("stop"));

refreshButton.addEventListener("click", async () => {
  await refreshStatus();
});

for (const button of connectButtons) {
  button.addEventListener("click", async () => {
    await connectDevice(button.dataset.device);
  });
}

mountParkBtn.addEventListener("click", async () => runMountCommand("park"));
mountUnparkBtn.addEventListener("click", async () => runMountCommand("unpark"));
mountHomeBtn.addEventListener("click", async () => runMountCommand("home"));

cameraCaptureForm.addEventListener("submit", async (event) => {
  event.preventDefault();
  await runCameraCapture();
});

filterwheelForm.addEventListener("submit", async (event) => {
  event.preventDefault();
  await changeFilter();
});

document.getElementById("fw-verify-btn")?.addEventListener("click", () => verifyFilterPosition());

focuserInBtn.addEventListener("click", async () => moveFocuser("in"));
focuserOutBtn.addEventListener("click", async () => moveFocuser("out"));

// ── Dust Cover ────────────────────────────────────────────────────────────────

// NINA returns CoverState as a string ("Open", "Closed", "Moving", …) or int
function normCoverStateClient(raw) {
  if (typeof raw === "number") return raw;
  const map = { notpresent: 0, closed: 1, moving: 2, open: 3, unknown: 100, error: 101 };
  return map[String(raw ?? "").toLowerCase().trim()] ?? 100;
}

const COVER_STATE_LABELS = { 0: "Not Present", 1: "Closed", 2: "Moving…", 3: "Open", 100: "Unknown", 101: "Error" };
const COVER_STATE_COLORS = { 0: "#6b7280", 1: "#4ade80", 2: "#fbbf24", 3: "#34d399", 100: "#9ca3af", 101: "#f87171" };

function updateCoverBadge(fd, flatdeviceConnected) {
  const stateEl  = document.getElementById("cover-state-label");
  const detailEl = document.getElementById("detail-flatdevice");
  const statusEl = deviceStatusEls.flatdevice;
  if (!fd) {
    if (stateEl) { stateEl.textContent = "—"; stateEl.style.color = ""; }
    if (detailEl) { detailEl.textContent = ""; detailEl.style.color = ""; }
    if (statusEl && flatdeviceConnected === undefined) {
      statusEl.textContent = "Unknown";
      statusEl.className = "offline";
    }
    return;
  }
  const stateNum = normCoverStateClient(fd.CoverState);
  const label    = COVER_STATE_LABELS[stateNum] ?? String(fd.CoverState ?? "Unknown");
  const color    = COVER_STATE_COLORS[stateNum] ?? "#9ca3af";
  if (stateEl)  { stateEl.textContent = label; stateEl.style.color = color; }
  if (detailEl) {
    const name = fd.DisplayName || fd.Name || "";
    detailEl.textContent = fd.Connected
      ? (name ? `${name} · ${label}` : label)
      : (name ? `${name} · not connected` : "Not connected");
    detailEl.style.color = fd.Connected ? color : "#6b7280";
  }
  if (statusEl) {
    if (fd.Connected) {
      statusEl.textContent = label;
      statusEl.className = "online";
    } else {
      statusEl.textContent = "Disconnected";
      statusEl.className = "offline";
    }
  }
}

async function coverAction(action) {
  const msgEl = document.getElementById("cover-result-msg");
  saveNinaConfigToStorage();
  if (msgEl) msgEl.innerHTML = `<span style="color:#6b7280">${action === "open" ? "Opening" : "Closing"} cover…</span>`;
  try {
    const resp = await fetch(`/api/nina/actions/flatdevice/${action}`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(currentConfig()),
    });
    const r = await resp.json();
    if (msgEl) {
      if (r.success) {
        msgEl.innerHTML = `<span style="color:#34d399">✔ Cover ${action === "open" ? "opened" : "closed"} ✓</span>`;
      } else if (resp.status === 403) {
        msgEl.innerHTML = `<span style="color:#f87171">⛔ ${r.error}</span>`;
      } else {
        msgEl.innerHTML = `<span style="color:#f87171">✘ ${r.error || "Failed"}</span>`;
      }
    }
    if (!r.success) setLog(r);
    await refreshStatus();
  } catch (err) {
    if (msgEl) msgEl.innerHTML = `<span style="color:#f87171">✘ ${err.message}</span>`;
    setLog({ success: false, action: `cover-${action}`, error: err.message });
  }
}

document.getElementById("cover-open-btn").addEventListener("click", () => coverAction("open"));
document.getElementById("cover-close-btn").addEventListener("click", () => coverAction("close"));
document.getElementById("cover-connect-btn").addEventListener("click", async () => {
  const msgEl = document.getElementById("cover-result-msg");
  if (msgEl) msgEl.innerHTML = `<span style="color:#6b7280">Connecting…</span>`;
  saveNinaConfigToStorage();
  try {
    const resp = await fetch("/api/nina/actions/flatdevice/connect", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(currentConfig()),
    });
    const r = await resp.json();
    if (msgEl) {
      msgEl.innerHTML = r.success
        ? `<span style="color:#34d399">✔ Connected</span>`
        : `<span style="color:#f87171">✘ ${r.error || "Failed"}</span>`;
    }
    if (!r.success) setLog(r);
    if (r.flatDevice) updateCoverBadge(r.flatDevice, r.devices?.flatdevice);
    await refreshStatus();
  } catch (err) {
    if (msgEl) msgEl.innerHTML = `<span style="color:#f87171">✘ ${err.message}</span>`;
    setLog({ success: false, action: "cover-connect", error: err.message });
  }
});

async function loadDefaults() {
  try {
    const response = await fetch("/api/config/defaults");
    const defaults = await response.json();
    document.getElementById("host").value = defaults.host || "192.168.1.174";
    document.getElementById("port").value = defaults.port || "1888";
    document.getElementById("protocol").value = defaults.protocol || "http";
    applyNinaConfigFromStorage();
    ocsHost = defaults.ocsHost || "192.168.1.220";
    syncOcsHostInputs(ocsHost);
    if (defaults.tns?.botId) document.getElementById("tns-bot-id").value = defaults.tns.botId;
    if (defaults.tns?.botName) document.getElementById("tns-bot-name").value = defaults.tns.botName;
    if (defaults.tns?.hasApiKey) document.getElementById("tns-api-key").placeholder = "Loaded from server env";
    // Restore Astro-COLIBRI UID from server DB (persists across browsers/sessions)
    if (defaults.colibriUid) {
      const el = document.getElementById("colibri-uid");
      if (el && !el.value) el.value = defaults.colibriUid;
    }
  } catch {
    document.getElementById("host").value = "192.168.1.174";
    document.getElementById("port").value = "1888";
    document.getElementById("protocol").value = "http";
    applyNinaConfigFromStorage();
    syncOcsHostInputs("192.168.1.220");
  }
  // Restore tonight search settings from localStorage
  ["tonight-days", "tonight-minalt"].forEach(id => {
    const saved = localStorage.getItem(id);
    const el    = document.getElementById(id);
    if (el && saved) el.value = saved;
  });
  // Update NINA host display in connect-bar
  const ninaHostEl = document.getElementById("nina-host-display");
  if (ninaHostEl) {
    const host = document.getElementById("host").value;
    const port = document.getElementById("port").value;
    ninaHostEl.textContent = `${host}:${port}`;
  }
  // Update OCS host display
  updateOcsIndicator();
  // Check STDWeb reachability
  checkStdwebHealth();
}

// ── Watchdog ──────────────────────────────────────────────────────────────────

async function loadWatchdogConfig() {
  try {
    const res  = await fetch("/api/watchdog");
    if (!res.ok) return;
    const cfg  = await res.json();
    document.getElementById("watchdog-enable").checked      = !!cfg.enabled;
    document.getElementById("watchdog-sky-temp").value      = cfg.skyTempLimit ?? -4;
    document.getElementById("watchdog-sq").value            = cfg.sqLimit      ?? 16;
    document.getElementById("watchdog-retention").value     = cfg.retentionMin ?? 10;
    const mph = document.getElementById("watchdog-morning-park");
    if (mph) mph.value = cfg.morningParkHour ?? 8;
    const iroof = document.getElementById("watchdog-ignore-roof");
    if (iroof) iroof.checked = !!cfg.ignoreRoof;
    updateWatchdogStatusPanel(cfg);
  } catch {}
}

function updateWatchdogStatusPanel(cfg) {
  const dot   = document.getElementById("watchdog-state-dot");
  const text  = document.getElementById("watchdog-state-text");
  const last  = document.getElementById("watchdog-last-check");
  const badge = document.getElementById("watchdog-status-badge");
  if (!dot) return;

  const colors = { off: "#6b7280", clear: "#22c55e", bad: "#ef4444", recovering: "#f59e0b" };
  const labels = { off: "Disabled", clear: "Clear ✓", bad: "⚠ Bad conditions — mount parked", recovering: "Recovering…" };

  const state = cfg.state || "off";
  dot.style.background = colors[state] || "#6b7280";
  text.textContent     = labels[state] || state;
  if (cfg.lastMsg && state !== "off") text.textContent += `  (${cfg.lastMsg})`;
  if (state === "bad" && cfg.badSince) {
    const badMin = Math.round((Date.now() - new Date(cfg.badSince).getTime()) / 60_000);
    text.textContent += `  — bad for ${badMin} min / max ${cfg.maxBadMin ?? 90} min`;
  }
  if (cfg.lastCheck) {
    const d = new Date(cfg.lastCheck);
    last.textContent = `Last check: ${d.toLocaleTimeString()}`;
  } else {
    last.textContent = cfg.enabled ? "Not checked yet" : "";
  }

  // Mount park monitor row
  const mmDot  = document.getElementById("mount-monitor-dot");
  const mmText = document.getElementById("mount-monitor-text");
  const mmLast = document.getElementById("mount-monitor-last");
  if (mmDot && mmText) {
    const mm = cfg.mountMonitor;
    const atPark = mm?.lastAtPark;
    if (atPark === null || atPark === undefined) {
      mmDot.style.background = "#6b7280";
      mmText.textContent = "Park monitor: NINA unreachable";
      mmText.style.color = "#6b7280";
    } else if (atPark === true) {
      mmDot.style.background = "#f59e0b";
      mmText.textContent = "Park monitor: Mount PARKED — cover closed";
      mmText.style.color = "#f59e0b";
    } else {
      mmDot.style.background = "#22c55e";
      mmText.textContent = "Park monitor: Mount tracking ✓";
      mmText.style.color = "#94a3b8";
    }
    if (mm?.lastCheck && mmLast) {
      mmLast.textContent = `Last: ${new Date(mm.lastCheck).toLocaleTimeString()}`;
    }
  }

  // Header badge
  if (badge) {
    if (state === "bad") {
      badge.textContent = "⛅ Watchdog: bad";
      badge.className   = "watchdog-badge watchdog-badge-bad";
      badge.style.display = "";
    } else if (state === "recovering") {
      badge.textContent = "⏳ Recovering";
      badge.className   = "watchdog-badge watchdog-badge-recovering";
      badge.style.display = "";
    } else if (state === "clear") {
      badge.textContent = "☀ Clear";
      badge.className   = "watchdog-badge watchdog-badge-clear";
      badge.style.display = "";
    } else {
      badge.style.display = "none";
    }
  }
}

document.getElementById("watchdog-save")?.addEventListener("click", async () => {
  const body = {
    enabled:      document.getElementById("watchdog-enable").checked,
    skyTempLimit: Number(document.getElementById("watchdog-sky-temp").value),
    sqLimit:      Number(document.getElementById("watchdog-sq").value),
    retentionMin: Number(document.getElementById("watchdog-retention").value),
    morningParkHour: Number(document.getElementById("watchdog-morning-park")?.value ?? 8),
    ignoreRoof:   document.getElementById("watchdog-ignore-roof")?.checked ?? false,
  };
  try {
    const res = await fetch("/api/watchdog", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    const data = await res.json();
    updateWatchdogStatusPanel(data.watchdog ?? data);
    // Flash the button to confirm
    const btn = document.getElementById("watchdog-save");
    btn.textContent = "Saved ✓";
    setTimeout(() => { btn.textContent = "Save Watchdog"; }, 2000);
  } catch {}
});

// Poll watchdog status every 30 s when settings tab is active
let _watchdogPollInterval = null;
function startWatchdogPolling() {
  if (_watchdogPollInterval) return;
  _watchdogPollInterval = setInterval(async () => {
    try {
      const res = await fetch("/api/watchdog");
      if (res.ok) updateWatchdogStatusPanel(await res.json());
    } catch {}
  }, 30_000);
}
function stopWatchdogPolling() {
  if (_watchdogPollInterval) { clearInterval(_watchdogPollInterval); _watchdogPollInterval = null; }
}

// Also poll when header badge is visible (any tab)
setInterval(async () => {
  try {
    const res = await fetch("/api/watchdog");
    if (res.ok) updateWatchdogStatusPanel(await res.json());
  } catch {}
}, 60_000);

// ── Astro-COLIBRI UID: save to server DB on blur ──────────────────────────────
document.getElementById("colibri-uid")?.addEventListener("change", async () => {
  const uid = document.getElementById("colibri-uid").value.trim();
  await fetch("/api/settings", {
    method:  "POST",
    headers: { "Content-Type": "application/json" },
    body:    JSON.stringify({ key: "colibri_uid", value: uid }),
  }).catch(() => {});
});
["tonight-days", "tonight-minalt"].forEach(id => {
  document.getElementById(id)?.addEventListener("change", () => {
    localStorage.setItem(id, document.getElementById(id).value);
  });
});

// ── Plan Tonight ──────────────────────────────────────────────────────────────

function formatErr(deg) {
  if (deg == null || deg <= 0) return "—";
  const arcsec = deg * 3600;
  if (arcsec < 60)   return arcsec.toFixed(0) + "″";
  if (arcsec < 3600) return (arcsec / 60).toFixed(1) + "′";
  return deg.toFixed(2) + "°";
}

function formatAge(hours) {
  if (hours == null || hours < 0) return "—";
  if (hours < 1)   return Math.round(hours * 60) + "m";
  if (hours < 24)  return hours.toFixed(1) + "h";
  return (hours / 24).toFixed(1) + "d";
}

function tonightSortBy(col) {
  if (_tonightSort.col === col) {
    _tonightSort.dir *= -1;
  } else {
    _tonightSort.col = col;
    // Default: descending for numbers, ascending for strings
    _tonightSort.dir = ["maxAlt", "windowMin"].includes(col) ? -1 : 1;
  }
  renderTonightTable();
}

function renderTonightTable() {
  const table      = document.getElementById("tonight-table");
  const tbody      = document.getElementById("tonight-tbody");
  const empty      = document.getElementById("tonight-empty");
  const filtersDiv = document.getElementById("tonight-filters");

  // Collect active type filters
  const checkboxes  = filtersDiv ? filtersDiv.querySelectorAll("input[type=checkbox]") : [];
  const activeTypes = new Set();
  checkboxes.forEach(cb => { if (cb.checked) activeTypes.add(cb.value); });

  let rows = _tonightEvents.filter(e =>
    activeTypes.size === 0 || activeTypes.has(e.type)
  );

  // Sort
  const { col, dir } = _tonightSort;
  rows.sort((a, b) => {
    let av = a[col], bv = b[col];
    // Treat dash / null / -1 (unknown mag) as worst
    const isNull = v => v == null || v === "—" || v === -1;
    if (isNull(av) && isNull(bv)) return 0;
    if (isNull(av)) return 1;
    if (isNull(bv)) return -1;
    if (typeof av === "number" && typeof bv === "number") return (av - bv) * dir;
    return String(av).localeCompare(String(bv)) * dir;
  });

  // Update sort indicators
  document.querySelectorAll("#tonight-table th[data-col]").forEach(th => {
    const ind = th.querySelector(".sort-ind");
    if (!ind) return;
    ind.textContent = th.dataset.col === col ? (dir === 1 ? " ▲" : " ▼") : " ⇅";
    th.classList.toggle("sorted", th.dataset.col === col);
  });

  tbody.innerHTML = "";
  if (!rows.length) {
    empty.textContent  = "No transients match the selected filters.";
    empty.style.display = "block";
    table.style.display = "none";
    return;
  }
  empty.style.display = "none";

  rows.forEach(e => {
    const winTxt = e.windowStart !== "—" ? `${e.windowStart}–${e.windowEnd} UTC` : "—";
    const durTxt = e.windowMin > 0 ? `${e.windowMin} min` : "—";
    const altTxt = e.maxAlt != null ? `${e.maxAlt}°` : "—";
    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td><a href="${e.tnsUrl}" target="_blank" style="color:#60a5fa;">${e.name}</a></td>
      <td style="white-space:nowrap;">
        <button class="btn-small tonight-add-btn"
                data-name="${e.name.replace(/"/g,'&quot;').replace(/'/g,'&#39;')}"
                data-ra="${e.ra}" data-dec="${e.dec}">+ Add</button>
        <a href="${e.colUrl}" target="_blank" class="btn-small" style="text-decoration:none;display:inline-block;" title="${e.acId || e.colUrl}">🔗 AC</a>
      </td>
      <td><span class="tonight-type-badge tonight-type-${e.type}">${e.type}</span></td>
      <td style="white-space:nowrap;">${formatAge(e.ageHours)}</td>
      <td style="white-space:nowrap;">${formatErr(e.errDeg)}</td>
      <td>${altTxt}</td>
      <td>${winTxt}</td>
      <td>${durTxt}</td>
      <td>${e.transitUtc !== "—" ? e.transitUtc + " UTC" : "—"}</td>
      <td>${e.lastMag !== "-1.0" && e.lastMag !== "-1" ? e.lastMag : "—"}</td>
    `;
    tbody.appendChild(tr);
  });

  table.style.display = "table";
}

async function planTonight() {
  const uid    = (document.getElementById("colibri-uid")?.value || "").trim();
  const days   = parseInt(document.getElementById("tonight-days")?.value)   || 7;
  const minAlt = parseInt(document.getElementById("tonight-minalt")?.value) || 20;

  if (!uid) {
    alert("Please enter your Astro-COLIBRI User ID first.\nFind it at astro-colibri.com/account");
    return;
  }

  const btn  = document.getElementById("plan-tonight-btn");
  const orig = btn.textContent;
  btn.disabled    = true;
  btn.textContent = "⏳ Loading…";

  const panel      = document.getElementById("tonight-panel");
  const header     = document.getElementById("tonight-header");
  const empty      = document.getElementById("tonight-empty");
  const table      = document.getElementById("tonight-table");
  const filtersDiv = document.getElementById("tonight-filters");

  panel.style.display  = "block";
  table.style.display  = "none";
  empty.style.display  = "none";
  header.textContent   = "Querying Astro-COLIBRI…";

  try {
    let lat = "", lon = "";
    if (typeof observerLocation !== "undefined" && observerLocation?.lat) {
      lat = observerLocation.lat;
      lon = observerLocation.lon;
    } else {
      try {
        const eq = await postJson("/api/nina/devices/status", currentConfig());
        const m  = eq?.equipmentInfo?.Mount;
        if (m?.SiteLatitude && m?.SiteLongitude &&
            !(m.SiteLatitude === 0 && m.SiteLongitude === 0)) {
          lat = m.SiteLatitude;
          lon = m.SiteLongitude;
        }
      } catch { /* no mount — proceed without location */ }
    }

    const params = new URLSearchParams({ uid, days, minAlt });
    if (lat) { params.set("lat", lat); params.set("lon", lon); }

    const data = await fetch(`/api/tonight?${params}`).then(r => r.json());

    if (!data.success) {
      header.textContent = `Error: ${data.error}${data.detail ? " — " + data.detail : ""}`;
      return;
    }

    const hasPos = lat !== "";
    const n      = data.events.length;
    header.textContent = `${n} transient${n !== 1 ? "s" : ""} in the last ${data.days} days` +
      (hasPos ? ` · visible tonight (alt ≥ ${minAlt}°)` : " · no observer location for altitude filter") +
      ` — ${new Date().toLocaleTimeString()}`;

    _tonightEvents = data.events || [];

    if (!_tonightEvents.length) {
      empty.textContent  = hasPos
        ? `No transients above ${minAlt}° altitude tonight. Try increasing "Days back" or lowering "Min altitude".`
        : "No transients found for this time range.";
      empty.style.display = "block";
      return;
    }

    // Build type filter checkboxes from unique types in this result set
    if (filtersDiv) {
      const typeCounts = {};
      _tonightEvents.forEach(e => { typeCounts[e.type] = (typeCounts[e.type] || 0) + 1; });
      const types = Object.keys(typeCounts).sort();
      filtersDiv.innerHTML =
        `<span class="tonight-filter-label">Filter:</span>` +
        types.map(t =>
          `<label class="tonight-filter-chip">` +
          `<input type="checkbox" value="${t}" checked onchange="renderTonightTable()">` +
          `<span class="tonight-type-badge tonight-type-${t}">${t}</span>` +
          `<span class="tonight-filter-count">${typeCounts[t]}</span>` +
          `</label>`
        ).join("") +
        `<button class="btn-small" style="margin-left:8px;" onclick="tonightSelectAllTypes(true)">All</button>` +
        `<button class="btn-small" onclick="tonightSelectAllTypes(false)">None</button>`;
      filtersDiv.style.display = "flex";
    }

    renderTonightTable();
  } catch (err) {
    header.textContent = `Error: ${err.message}`;
  } finally {
    btn.disabled    = false;
    btn.textContent = orig;
  }
}

function tonightSelectAllTypes(checked) {
  const filtersDiv = document.getElementById("tonight-filters");
  if (!filtersDiv) return;
  filtersDiv.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = checked; });
  renderTonightTable();
}

function addTonightTarget(name, raDeg, decDeg, btn) {
  const raH  = raDeg / 15;
  const raHH = Math.floor(raH);
  const raM  = Math.floor((raH - raHH) * 60);
  const raS  = (((raH - raHH) * 60 - raM) * 60).toFixed(1);
  const raStr = `${String(raHH).padStart(2,"0")}:${String(raM).padStart(2,"0")}:${String(raS).padStart(4,"0")}`;
  const sign  = decDeg < 0 ? "-" : "+";
  const absD  = Math.abs(decDeg);
  const dD    = Math.floor(absD);
  const dM    = Math.floor((absD - dD) * 60);
  const dS    = (((absD - dD) * 60 - dM) * 60).toFixed(0);
  const decStr= `${sign}${String(dD).padStart(2,"0")}:${String(dM).padStart(2,"0")}:${String(dS).padStart(2,"0")}`;

  if (btn) { btn.disabled = true; btn.textContent = "…"; }

  fetch("/api/sequence/queue", {
    method:  "POST",
    headers: { "Content-Type": "application/json" },
    body:    JSON.stringify({ name, ra: raStr, dec: decStr, raDeg, decDeg }),
  })
    .then(r => r.json())
    .then(d => {
      if (d.success) {
        if (btn) { btn.textContent = "✓ Added"; }
        if (typeof refreshSeqState === "function") refreshSeqState();
      } else {
        if (btn) { btn.disabled = false; btn.textContent = "+ Add"; }
        alert("Failed to add: " + (d.error || "unknown error"));
      }
    })
    .catch(err => {
      if (btn) { btn.disabled = false; btn.textContent = "+ Add"; }
      alert("Error: " + err.message);
    });
}

document.getElementById("plan-tonight-btn")?.addEventListener("click", planTonight);

document.getElementById("tonight-tbody")?.addEventListener("click", function(ev) {
  const btn = ev.target.closest(".tonight-add-btn");
  if (!btn) return;
  addTonightTarget(btn.dataset.name, parseFloat(btn.dataset.ra), parseFloat(btn.dataset.dec), btn);
});

// ── Sequence / ToDo tab ──────────────────────────────────────────────────────

function seqConfig() {
  const filterRaw = document.getElementById("seq-filters").value;
  return {
    duration:       Number(document.getElementById("seq-duration").value)       || 30,
    gain:           Number(document.getElementById("seq-gain").value)           || 10,
    count:          Number(document.getElementById("seq-count").value)          || 40,
    filters:        filterRaw.split(",").map(s => s.trim()).filter(Boolean),
    solveEnabled:              document.getElementById("seq-solve-enable").checked,
    solveExp:                  Number(document.getElementById("seq-solve-exp").value)                  || 5,
    solveThreshold:            Number(document.getElementById("seq-solve-threshold").value)            || 60,
    frameCheckEnabled:         document.getElementById("set-frame-check-enable")?.checked    ?? false,
    frameCheckThresholdArcmin: Number(document.getElementById("set-frame-check-threshold")?.value) || 5,
    manualMode:                document.getElementById("seq-manual-mode").checked,
    minAlt:                    Number(document.getElementById("set-minalt")?.value)          || 20,
    meridianGap:               Number(document.getElementById("set-meridian-gap")?.value)   || 10,
    zenithLimit:               Number(document.getElementById("set-zenith-limit")?.value)   || 70,
    minStars:                  Number(document.getElementById("set-min-stars")?.value)      ?? 10,
    maxCloudRecoveryWaitMin:   Number(document.getElementById("set-max-cloud-recovery-min")?.value) ?? 12,
  };
}

// ── Arbiter / Safety limits — load from DB, save to DB ───────────────────────
async function loadArbiterConfig() {
  try {
    const keys = [
      "minAlt", "zenithLimit", "meridianGap",
      "minStars", "maxCloudRecoveryWaitMin",
      "frameCheckEnabled", "frameCheckThreshold",
      "daily_reset_hour",
      "telescope_fov_deg",
      "af_slope_G", "af_slope_BP", "af_slope_RP",
      "af_ref15_G", "af_ref15_BP", "af_ref15_RP",
    ];
    const results = await Promise.all(keys.map(k =>
      fetch(`/api/settings/${k}`).then(r => r.ok ? r.json() : null).catch(() => null)
    ));
    const [minAlt, zenithLimit, meridianGap, minStars, maxCloudRecoveryWaitMin,
           frameCheckEnabled, frameCheckThreshold, dailyResetHour,
           telescopeFovDeg,
           afSlopeG, afSlopeBP, afSlopeRP,
           afRef15G, afRef15BP, afRef15RP] = results;
    if (minAlt          != null) document.getElementById("set-minalt").value                   = minAlt.value          ?? 20;
    if (zenithLimit     != null) document.getElementById("set-zenith-limit").value             = zenithLimit.value     ?? 70;
    if (meridianGap     != null) document.getElementById("set-meridian-gap").value             = meridianGap.value     ?? 10;
    if (minStars        != null) document.getElementById("set-min-stars").value                = minStars.value        ?? 10;
    if (maxCloudRecoveryWaitMin != null) document.getElementById("set-max-cloud-recovery-min").value = maxCloudRecoveryWaitMin.value ?? 12;
    if (frameCheckEnabled != null) document.getElementById("set-frame-check-enable").checked  = frameCheckEnabled.value === "true";
    if (frameCheckThreshold != null) document.getElementById("set-frame-check-threshold").value = frameCheckThreshold.value ?? 5;
    if (dailyResetHour  != null) { const el = document.getElementById("set-daily-reset-hour"); if (el) el.value = dailyResetHour.value ?? 12; }
    if (telescopeFovDeg != null) {
      const el = document.getElementById("set-telescope-fov-deg");
      if (el) el.value = telescopeFovDeg.value ?? "0.86";
    }
    updateTelescopeFovPreview();
    if (afSlopeG  != null) { const el = document.getElementById("set-af-slope-G");  if (el) el.value = afSlopeG.value  ?? 10.83; }
    if (afSlopeBP != null) { const el = document.getElementById("set-af-slope-BP"); if (el) el.value = afSlopeBP.value ?? 9.73; }
    if (afSlopeRP != null) { const el = document.getElementById("set-af-slope-RP"); if (el) el.value = afSlopeRP.value ?? 12.08; }
    if (afRef15G  != null) { const el = document.getElementById("set-af-ref15-G");  if (el) el.value = afRef15G.value  ?? 4195; }
    if (afRef15BP != null) { const el = document.getElementById("set-af-ref15-BP"); if (el) el.value = afRef15BP.value ?? 4202; }
    if (afRef15RP != null) { const el = document.getElementById("set-af-ref15-RP"); if (el) el.value = afRef15RP.value ?? 4202; }
    _updateAfModelPreview();
  } catch { /* non-fatal */ }
}

function updateTelescopeFovPreview() {
  const el = document.getElementById("telescope-fov-preview");
  const input = document.getElementById("set-telescope-fov-deg");
  if (!el || !input) return;
  const fov = parseFloat(input.value);
  if (!Number.isFinite(fov) || fov <= 0) {
    el.textContent = "";
    return;
  }
  const maxSide = 8;
  const capSpan = ((maxSide - 1) * fov).toFixed(2);
  const capWidth = (maxSide * fov).toFixed(2);
  const singleThreshold = (fov / 2).toFixed(2);
  const exampleErr = 1.0;
  const n = Math.min(maxSide, Math.max(1, Math.ceil((2 * exampleErr) / fov)));
  const exSpan = ((n - 1) * fov).toFixed(2);
  el.textContent =
    `Single pointing when err ≤ ${singleThreshold}° (FOV/2). ` +
    `Max mosaic (8×8 cap): ~${capSpan}° between outer tile centers, ~${capWidth}° full width. ` +
    `Example: err=1.0° → ${n}×${n} tiles (~${exSpan}° between outer centers).`;
}

document.getElementById("set-telescope-fov-deg")?.addEventListener("input", updateTelescopeFovPreview);

document.getElementById("telescope-fov-save")?.addEventListener("click", async () => {
  const input = document.getElementById("set-telescope-fov-deg");
  const status = document.getElementById("telescope-fov-status");
  const btn = document.getElementById("telescope-fov-save");
  const val = parseFloat(input?.value);
  if (!Number.isFinite(val) || val <= 0.05 || val > 30) {
    if (status) status.textContent = "Enter a value between 0.05 and 30°";
    return;
  }
  try {
    const resp = await fetch("/api/settings", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ key: "telescope_fov_deg", value: String(val) }),
    });
    const data = await resp.json();
    if (!resp.ok || !data.success) {
      if (status) status.textContent = data.error || "Save failed";
      return;
    }
    if (status) status.textContent = "Saved ✓";
    updateTelescopeFovPreview();
    const orig = btn.textContent;
    btn.textContent = "Saved ✓";
    btn.disabled = true;
    setTimeout(() => { btn.textContent = orig; btn.disabled = false; if (status) status.textContent = ""; }, 2000);
  } catch (e) {
    if (status) status.textContent = e.message;
  }
});

document.getElementById("arbiter-save")?.addEventListener("click", async () => {
  const settings = {
    minAlt:               document.getElementById("set-minalt").value,
    zenithLimit:          document.getElementById("set-zenith-limit").value,
    meridianGap:          document.getElementById("set-meridian-gap").value,
    frameCheckEnabled:    String(document.getElementById("set-frame-check-enable").checked),
    frameCheckThreshold:  document.getElementById("set-frame-check-threshold").value,
  };
  try {
    await Promise.all(Object.entries(settings).map(([key, value]) =>
      fetch("/api/settings", { method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ key, value }) })
    ));
    const btn = document.getElementById("arbiter-save");
    const orig = btn.textContent;
    btn.textContent = "Saved ✓"; btn.disabled = true;
    setTimeout(() => { btn.textContent = orig; btn.disabled = false; }, 2000);
  } catch (e) {
    alert("Failed to save arbiter settings: " + e.message);
  }
});

document.getElementById("quality-save")?.addEventListener("click", async () => {
  const settings = {
    minStars:               document.getElementById("set-min-stars").value,
    maxCloudRecoveryWaitMin: document.getElementById("set-max-cloud-recovery-min").value,
  };
  try {
    await Promise.all(Object.entries(settings).map(([key, value]) =>
      fetch("/api/settings", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ key, value }),
      })
    ));
    const btn = document.getElementById("quality-save");
    const orig = btn.textContent;
    btn.textContent = "Saved ✓";
    btn.disabled = true;
    setTimeout(() => { btn.textContent = orig; btn.disabled = false; }, 2000);
  } catch (e) {
    alert("Failed to save quality gates: " + e.message);
  }
});

// ── Daily reset hour save ─────────────────────────────────────────────────────
document.getElementById("daily-reset-hour-save")?.addEventListener("click", async () => {
  const val = document.getElementById("set-daily-reset-hour")?.value;
  try {
    await fetch("/api/settings", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ key: "daily_reset_hour", value: val }),
    });
    const btn = document.getElementById("daily-reset-hour-save");
    btn.textContent = "Saved ✓"; btn.disabled = true;
    setTimeout(() => { btn.textContent = "Save"; btn.disabled = false; }, 2000);
  } catch (e) { alert("Save failed: " + e.message); }
});

// ── Focuser temperature model ─────────────────────────────────────────────────
function _updateAfModelPreview() {
  const el = document.getElementById("af-model-preview");
  if (!el) return;
  const filters = ["G", "BP", "RP"];
  const parts = filters.map(f => {
    const slope = parseFloat(document.getElementById(`set-af-slope-${f}`)?.value);
    const ref15 = parseFloat(document.getElementById(`set-af-ref15-${f}`)?.value);
    if (!slope || !ref15) return null;
    const p10 = Math.round(slope * (10 - 15) + ref15);
    const p15 = Math.round(ref15);
    const p20 = Math.round(slope * (20 - 15) + ref15);
    return `${f}: ${p10}@10° / ${p15}@15° / ${p20}@20°`;
  }).filter(Boolean);
  el.textContent = parts.length ? "Preview — " + parts.join("  ·  ") : "";
}
["G","BP","RP"].forEach(f => {
  document.getElementById(`set-af-slope-${f}`)?.addEventListener("input", _updateAfModelPreview);
  document.getElementById(`set-af-ref15-${f}`)?.addEventListener("input", _updateAfModelPreview);
});

document.getElementById("af-model-save")?.addEventListener("click", async () => {
  const settings = {};
  ["G","BP","RP"].forEach(f => {
    settings[`af_slope_${f}`] = document.getElementById(`set-af-slope-${f}`)?.value;
    settings[`af_ref15_${f}`] = document.getElementById(`set-af-ref15-${f}`)?.value;
  });
  try {
    await Promise.all(Object.entries(settings).map(([key, value]) =>
      fetch("/api/settings", { method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ key, value }) })
    ));
    const btn = document.getElementById("af-model-save");
    btn.textContent = "Saved ✓"; btn.disabled = true;
    setTimeout(() => { btn.textContent = "Save"; btn.disabled = false; }, 2000);
  } catch (e) { alert("Save failed: " + e.message); }
});

// ── Daily reset button (Night Plan tab) ──────────────────────────────────────
document.getElementById("seq-daily-reset-btn")?.addEventListener("click", async () => {
  const info = await (await fetch("/api/sequence/daily-reset/info")).json();
  const alertNames  = info.alertTargets?.join(", ")  || "(none)";
  const manualNames = info.manualTargets?.join(", ") || "(none)";
  const ok = confirm(
    `Daily Night Plan Reset\n\n` +
    `Will REMOVE (alert targets): ${alertNames}\n` +
    `Will KEEP & RESET (monitored): ${manualNames}\n\n` +
    `Continue?`
  );
  if (!ok) return;
  const r = await fetch("/api/sequence/daily-reset", { method: "POST" });
  const d = await r.json();
  if (d.success) {
    const infoEl = document.getElementById("seq-daily-reset-info");
    if (infoEl) infoEl.textContent = `Reset done — removed: ${d.removed?.join(", ") || "none"}, kept: ${d.kept?.join(", ") || "none"}`;
    refreshSeqState();
  } else {
    alert("Reset failed: " + d.error);
  }
});

async function initQualityReplayDefaults() {
  const dateInput = document.getElementById("quality-replay-date");
  if (!dateInput) return;
  if (!dateInput.value) {
    dateInput.value = new Date().toISOString().slice(0, 10);
  }
  try {
    const res = await fetch("/api/pipeline/nas-dates");
    const data = await res.json();
    if (data?.success && Array.isArray(data.dates) && data.dates.length) {
      dateInput.value = data.dates[0];
    }
  } catch { /* non-fatal */ }
}

function renderQualityReplay(result) {
  const body = document.getElementById("quality-replay-body");
  if (!body) return;
  body.innerHTML = "";
  const activations = Array.isArray(result?.activations) ? result.activations : [];
  if (!activations.length) {
    const tr = document.createElement("tr");
    const td = document.createElement("td");
    td.colSpan = 5;
    td.style.padding = "8px";
    td.style.color = "#64748b";
    td.textContent = "No activation events for this night.";
    tr.appendChild(td);
    body.appendChild(tr);
    return;
  }
  activations.forEach((ev) => {
    const tr = document.createElement("tr");
    const cells = [
      ev.ts || "",
      ev.type || "",
      ev.target || "",
      ev.image || "",
      ev.reason || "",
    ];
    cells.forEach((txt, idx) => {
      const td = document.createElement("td");
      td.style.padding = "6px";
      td.style.borderBottom = "1px solid #1e293b";
      td.style.verticalAlign = "top";
      if (idx === 1) {
        td.style.color = String(txt).includes("CLOUD") ? "#fbbf24" : "#fca5a5";
      }
      td.textContent = txt;
      tr.appendChild(td);
    });
    body.appendChild(tr);
  });
}

document.getElementById("quality-replay-run")?.addEventListener("click", async () => {
  const date = document.getElementById("quality-replay-date")?.value;
  const maxEvents = Number(document.getElementById("quality-replay-max-events")?.value || 200);
  const status = document.getElementById("quality-replay-status");
  const btn = document.getElementById("quality-replay-run");
  if (!date) {
    if (status) status.textContent = "Pick a date first.";
    return;
  }
  if (btn) { btn.disabled = true; btn.textContent = "Running..."; }
  if (status) status.textContent = "Replaying night FITS (can take 10-60s)...";
  try {
    const q = new URLSearchParams({
      date,
      maxEvents: String(Math.max(10, Math.min(1000, Math.round(maxEvents) || 200))),
    });
    const res = await fetch(`/api/quality/replay?${q.toString()}`);
    const data = await res.json();
    if (!res.ok || !data?.success) throw new Error(data?.error || `HTTP ${res.status}`);
    renderQualityReplay(data);
    if (status) status.textContent = `${data.activationCount ?? 0} activations from ${data.frameCount ?? 0} frames (${data.date})`;
  } catch (e) {
    if (status) status.textContent = `Replay failed: ${e.message}`;
  } finally {
    if (btn) { btn.disabled = false; btn.textContent = "Run Replay"; }
  }
});

function renderSeqQueue(queue) {
  // Don't tear down the list while the user is mid-drag — it would destroy the
  // dragged element and reset the drag source, causing items to vanish.
  if (_seqDragging) return;

  const list  = document.getElementById("seq-queue-list");
  const empty = document.getElementById("seq-queue-empty");
  list.innerHTML = "";

  const nonWaitItems = queue.filter(t => !t.done && t.itemType !== "wait");
  if (!queue.length || (!nonWaitItems.length && !queue.some(t => t.itemType === "wait"))) {
    empty.style.display = queue.length === 0 ? "block" : "none";
  } else {
    empty.style.display = "none";
  }
  if (!queue.length) return;

  queue.forEach((t, i) => {
    const isWait   = t.itemType === "wait";
    const li       = document.createElement("li");
    li.className   = "seq-queue-item" + (t.done ? " done" : "") + (isWait ? " wait-item" : "");
    li.setAttribute("draggable", t.done ? "false" : "true");
    li.dataset.idx = i;

    if (isWait) {
      li.innerHTML = `
        <span class="seq-queue-num">${i + 1}</span>
        <span style="font-size:16px;flex-shrink:0;">⏳</span>
        <span class="seq-queue-name">Wait until</span>
        <input class="seq-queue-wait-input" type="time" value="${t.waitTime || "00:00"}" data-idx="${i}" />
        ${t.done
          ? `<span class="seq-queue-done-badge">✓ Done</span>`
          : `<button class="seq-queue-remove" data-idx="${i}" type="button">×</button>`}
      `;
    } else {
      const raH      = (t.raDeg || 0) / 15;
      const raHH     = Math.floor(raH);
      const raM      = String(Math.floor((raH - raHH) * 60)).padStart(2, "0");
      const decSign  = (t.decDeg || 0) >= 0 ? "+" : "";

      const gCount    = Number(document.getElementById("seq-count").value)    || 40;
      const gDuration = Number(document.getElementById("seq-duration").value) || 30;
      const gFilters  = (document.getElementById("seq-filters").value || "G,BP,RP")
                          .split(",").map(s => s.trim()).filter(Boolean);

      // Filter checkboxes — null means "use all global filters" (show all checked)
      // A subset array means explicit override (only those filters checked)
      const hasFilterOverride = Array.isArray(t.filters) && t.filters.length > 0;
      const activeFilters     = hasFilterOverride ? t.filters : gFilters;
      const filterCheckboxes  = gFilters.map(f => {
        const checked = activeFilters.includes(f) ? "checked" : "";
        return `<label class="seq-filter-cb-label">
          <input type="checkbox" class="seq-filter-cb" data-filter="${f}" data-idx="${i}" ${checked} />
          <span>${f}</span>
        </label>`;
      }).join("");

      // Badge summarising active overrides (filter badge only when a real subset is set)
      const overrideParts = [];
      if (hasFilterOverride) overrideParts.push(t.filters.join(","));
      if (t.count    != null) overrideParts.push(`${t.count}fr`);
      if (t.duration != null) overrideParts.push(`${t.duration}s`);
      const badge = overrideParts.length
        ? `<span class="seq-queue-override-badge" title="Per-target overrides">${overrideParts.join(" · ")}</span>`
        : "";

      li.innerHTML = `
        <div class="seq-queue-row">
          <span class="seq-queue-num">${i + 1}</span>
          <span class="seq-queue-name">${t.name}</span>
          ${t.source === "alert" ? `<span class="seq-source-badge alert-source" title="Added from alert broker — will be removed at daily reset">alert</span>` : `<span class="seq-source-badge manual-source" title="Manually added — kept across daily reset">monitored</span>`}
          <span class="seq-queue-coords">α&nbsp;${raHH}h${raM}m &nbsp;δ&nbsp;${decSign}${(t.decDeg || 0).toFixed(1)}°</span>
          ${badge}
          ${t.done
            ? `<span class="seq-queue-done-badge">✓ Done</span>`
            : `<button class="seq-queue-cfg-btn" data-idx="${i}" type="button" title="Per-target overrides">⚙</button>
               <button class="seq-queue-remove"  data-idx="${i}" type="button">×</button>`}
        </div>
        ${!t.done ? `
        <div class="seq-queue-override-panel" data-idx="${i}" style="display:none;">
          <div class="seq-override-row">
            <div class="seq-override-label">
              <span>Filters <small style="color:#334155;font-style:italic;">(unchecked = use all)</small></span>
              <div class="seq-filter-cb-group">${filterCheckboxes}</div>
            </div>
            <label class="seq-override-label">Frames / filter
              <input class="seq-override-input" type="number" min="1" step="1"
                     placeholder="${gCount}" value="${t.count ?? ""}"
                     data-field="count" data-idx="${i}" />
            </label>
            <label class="seq-override-label">Duration (s)
              <input class="seq-override-input" type="number" min="1" step="1"
                     placeholder="${gDuration}" value="${t.duration ?? ""}"
                     data-field="duration" data-idx="${i}" />
            </label>
          </div>
        </div>` : ""}
      `;
    }

    // ── Drag-and-drop ───────────────────────────────────────────────────────
    li.addEventListener("dragstart", e => {
      _seqDragSrc  = i;
      _seqDragging = true;
      e.dataTransfer.effectAllowed = "move";
      setTimeout(() => li.classList.add("dragging"), 0);
    });
    li.addEventListener("dragend", () => {
      li.classList.remove("dragging");
      // Only refresh when drag was CANCELLED (no drop fired).
      // If drop handled it, _seqDragSrc is already null and _seqDragging
      // will be cleared by the drop handler after its async work finishes.
      if (_seqDragSrc !== null) {
        _seqDragSrc  = null;
        _seqDragging = false;
        refreshSeqState();
      }
    });
    li.addEventListener("dragover", e => {
      e.preventDefault();
      const rowH = li.querySelector(".seq-queue-row")?.offsetHeight ?? li.offsetHeight;
      _seqDragInsertAfter = e.offsetY >= rowH / 2;
      li.classList.toggle("drag-insert-before", !_seqDragInsertAfter);
      li.classList.toggle("drag-insert-after",   _seqDragInsertAfter);
    });
    li.addEventListener("dragleave", () => {
      li.classList.remove("drag-insert-before", "drag-insert-after");
    });
    li.addEventListener("drop", async e => {
      e.preventDefault();
      li.classList.remove("drag-insert-before", "drag-insert-after");
      if (_seqDragSrc === null) return;

      const from = _seqDragSrc;
      // Clear src immediately so dragend (which fires next) knows drop was handled
      // and won't race with our async refresh below.
      _seqDragSrc = null;

      // gap = insertion point in original array (0 = before first, n = after last)
      const gap = _seqDragInsertAfter ? i + 1 : i;

      // No-op: drop lands right before or after the source item
      if (gap === from || gap === from + 1) {
        _seqDragging = false;
        await refreshSeqState();
        return;
      }

      // `to` = insertion index in the array AFTER `from` is removed
      const to = gap > from ? gap - 1 : gap;

      try {
        await fetch("/api/sequence/move", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ from, to }),
        });
      } finally {
        // Always unblock renders and refresh, even on error
        _seqDragging = false;
        await refreshSeqState();
      }
    });

    list.appendChild(li);
  });

  // Remove buttons
  for (const btn of list.querySelectorAll(".seq-queue-remove")) {
    btn.addEventListener("click", async () => {
      await fetch(`/api/sequence/queue/${btn.dataset.idx}`, { method: "DELETE" });
      await refreshSeqState();
      if (typeof refreshNightPlan === "function") refreshNightPlan();
    });
  }

  // Wait-time edits: update the wait item's time on change
  for (const inp of list.querySelectorAll(".seq-queue-wait-input")) {
    inp.addEventListener("change", async () => {
      const idx = parseInt(inp.dataset.idx);
      const newTime = inp.value; // HH:MM
      await fetch(`/api/sequence/queue/${idx}`, {
        method: "PATCH",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ waitTime: newTime }),
      });
      await refreshSeqState();
    });
  }

  // Restore panels that were open before re-render
  for (const idx of _openOverridePanels) {
    const panel = list.querySelector(`.seq-queue-override-panel[data-idx="${idx}"]`);
    if (panel) panel.style.display = "block";
  }

  // ⚙ Toggle override panel
  for (const btn of list.querySelectorAll(".seq-queue-cfg-btn")) {
    btn.addEventListener("click", e => {
      e.stopPropagation();
      const idx   = btn.dataset.idx;
      const panel = list.querySelector(`.seq-queue-override-panel[data-idx="${idx}"]`);
      if (!panel) return;
      const isOpen = panel.style.display !== "none";
      panel.style.display = isOpen ? "none" : "block";
      if (isOpen) _openOverridePanels.delete(idx);
      else        _openOverridePanels.add(idx);
    });
  }

  // Helper: update the override badge from current DOM state.
  // "all global filters checked" = no filter override (badge hidden for filters).
  function _updateOverrideBadge(li, idx) {
    const gFiltersNow = (document.getElementById("seq-filters")?.value || "G,BP,RP")
                          .split(",").map(s => s.trim()).filter(Boolean);
    const checked  = Array.from(li.querySelectorAll(`.seq-filter-cb[data-idx="${idx}"]:checked`))
                          .map(c => c.dataset.filter);
    const allChecked = gFiltersNow.every(f => checked.includes(f)) && checked.length === gFiltersNow.length;
    const countV    = li.querySelector(`[data-field="count"]`)?.value.trim()    || "";
    const durationV = li.querySelector(`[data-field="duration"]`)?.value.trim() || "";
    const parts = [];
    if (!allChecked && checked.length) parts.push(checked.join(","));
    if (countV)          parts.push(`${countV}fr`);
    if (durationV)       parts.push(`${durationV}s`);
    let badge = li.querySelector(".seq-queue-override-badge");
    if (parts.length) {
      if (!badge) {
        badge = document.createElement("span");
        badge.className = "seq-queue-override-badge";
        li.querySelector(".seq-queue-cfg-btn")?.before(badge);
      }
      badge.textContent = parts.join(" · ");
    } else if (badge) {
      badge.remove();
    }
  }

  // Filter checkboxes — instant save on click.
  // If all global filters are re-checked, save null to clear the override.
  for (const cb of list.querySelectorAll(".seq-filter-cb")) {
    cb.addEventListener("change", async () => {
      const idx = parseInt(cb.dataset.idx);
      const li  = cb.closest(".seq-queue-item");
      const gFiltersNow = (document.getElementById("seq-filters")?.value || "G,BP,RP")
                            .split(",").map(s => s.trim()).filter(Boolean);
      const checked = Array.from(li.querySelectorAll(`.seq-filter-cb[data-idx="${idx}"]:checked`))
                           .map(c => c.dataset.filter);
      const allChecked = gFiltersNow.every(f => checked.includes(f)) && checked.length === gFiltersNow.length;
      // Save null when all filters are checked (= use global default, no override needed)
      const filtersPayload = allChecked ? null : (checked.length ? checked : null);
      await fetch(`/api/sequence/queue/${idx}`, {
        method: "PATCH",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ filters: filtersPayload }),
      });
      _updateOverrideBadge(li, idx);
    });
  }

  // Frames / duration inputs — save on blur or Enter
  for (const inp of list.querySelectorAll(".seq-override-input")) {
    const save = async () => {
      const idx   = parseInt(inp.dataset.idx);
      const field = inp.dataset.field;
      const raw   = inp.value.trim();
      const value = raw === "" ? null : Math.max(1, Number(raw));
      await fetch(`/api/sequence/queue/${idx}`, {
        method: "PATCH",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ [field]: value }),
      });
      _updateOverrideBadge(inp.closest(".seq-queue-item"), idx);
    };
    inp.addEventListener("blur",    save);
    inp.addEventListener("keydown", e => { if (e.key === "Enter") { e.preventDefault(); inp.blur(); } });
  }
}

function renderSeqStatus(state) {
  const dot         = document.getElementById("seq-status-dot");
  const label       = document.getElementById("seq-status-label");
  const targetEl    = document.getElementById("seq-status-target");
  const stepEl      = document.getElementById("seq-status-step");
  const progWrap    = document.getElementById("seq-progress-wrap");
  const progBar     = document.getElementById("seq-progress-bar");
  const confirmRow  = document.getElementById("seq-confirm-row");
  const confirmStep = document.getElementById("seq-confirm-step");
  const runBtn      = document.getElementById("seq-run");
  const abortBtn    = document.getElementById("seq-abort");
  const nextBtn     = document.getElementById("seq-next");

  const paused = state.running && Boolean(state.waitingForStep);

  if (paused) {
    dot.className = "seq-dot seq-dot-paused";
    label.textContent = "Paused — waiting for confirmation";
    runBtn.disabled = true;
    abortBtn.disabled = false;
    nextBtn.disabled = false;
  } else if (state.running) {
    dot.className = "seq-dot seq-dot-running";
    label.textContent = "Running";
    runBtn.disabled = true;
    abortBtn.disabled = false;
    nextBtn.disabled = true;
  } else if (state.error) {
    dot.className = "seq-dot seq-dot-error";
    label.textContent = "Error";
    runBtn.disabled = !ocsConnected;
    runBtn.title = ocsConnected ? "" : "OCS not connected — cannot launch sequence";
    abortBtn.disabled = true;
    nextBtn.disabled = true;
  } else {
    const allDone = state.queue.length > 0 && state.queue.every(t => t.done);
    dot.className = "seq-dot " + (allDone ? "seq-dot-done" : "seq-dot-idle");
    label.textContent = allDone ? "All targets complete ✓" : "Idle";
    runBtn.disabled = !ocsConnected;
    runBtn.title = ocsConnected ? "" : "OCS not connected — cannot launch sequence";
    abortBtn.disabled = true;
    nextBtn.disabled = true;
  }

  // Confirmation bar
  if (state.waitingForStep) {
    confirmRow.style.display = "flex";
    confirmStep.textContent = state.waitingForStep;
  } else {
    confirmRow.style.display = "none";
  }

  targetEl.style.display = state.currentTarget ? "block" : "none";
  if (state.currentTarget) targetEl.textContent = `Target: ${state.currentTarget}`;

  stepEl.style.display = state.currentStep && !paused ? "block" : "none";
  if (state.currentStep && !paused) stepEl.textContent = state.currentStep;

  if (state.progress) {
    progWrap.style.display = "block";
    progBar.style.width = `${(state.progress.frame / state.progress.frames) * 100}%`;
  } else {
    progWrap.style.display = "none";
  }

  // Per-frame solve badge
  const frameSolveRow = document.getElementById("seq-frame-solve-row");
  if (frameSolveRow) {
    const fs = state.lastFrameSolve;
    if (!fs || !state.running) {
      frameSolveRow.style.display = "none";
    } else if (!fs.solved) {
      frameSolveRow.style.display = "block";
      frameSolveRow.innerHTML =
        `<span style="color:#f59e0b;">⚠ Frame check: solve failed</span>` +
        `<span style="color:#9ca3af;margin-left:6px;">${fs.filename || ""}</span>`;
    } else {
      const color  = fs.on_target ? "#22c55e" : "#ef4444";
      const icon   = fs.on_target ? "✓" : "✗";
      const label  = fs.on_target ? "on target" : "DRIFTED";
      const offset = fs.offset_arcmin != null ? `${fs.offset_arcmin}'` : "?";
      frameSolveRow.style.display = "block";
      frameSolveRow.innerHTML =
        `<span style="color:${color};font-weight:600;">${icon} Frame check: ${label}</span>` +
        `<span style="color:#9ca3af;margin-left:6px;">offset ${offset} — ${fs.filename || ""}</span>`;
    }
  }
}

function renderSeqLog(log) {
  if (log.length === seqLogRenderedCount) return;
  const el = document.getElementById("seq-log");
  // Append only new lines
  const newEntries = log.slice(seqLogRenderedCount);
  seqLogRenderedCount = log.length;
  for (const entry of newEntries) {
    const line = document.createElement("span");
    line.className = `seq-log-line ${entry.level || "info"}`;
    line.innerHTML = `<span class="seq-ts">[${entry.ts}]</span><span class="seq-msg">${entry.msg}</span>`;
    el.appendChild(line);
    el.appendChild(document.createTextNode("\n"));
  }
  el.scrollTop = el.scrollHeight;
}

async function refreshSeqState() {
  try {
    const res = await fetch("/api/sequence/state");
    const state = await res.json();
    // Don't re-render the queue while the user is typing in an override panel
    const editingOverride = document.activeElement?.closest(".seq-queue-override-panel");
    if (!editingOverride) renderSeqQueue(state.queue);
    renderSeqStatus(state);
    renderSeqLog(state.log);

    // Safety override modal — show/hide based on state
    if (state.waitingForSafetyOverride) {
      showSafetyModal(state.waitingForSafetyOverride);
    } else {
      hideSafetyModal();
    }
  } catch { /* ignore */ }
}

function startSeqPolling() {
  if (seqPollInterval) return;
  refreshSeqState();
  seqPollInterval = setInterval(refreshSeqState, 2000);
}

function stopSeqPolling() {
  if (!seqPollInterval) return;
  clearInterval(seqPollInterval);
  seqPollInterval = null;
}

async function addToQueue(target) {
  try {
    const result = await postJson("/api/sequence/queue", {
      name:   target.fullName || target.name,
      raDeg:  target.raDeg,
      decDeg: target.decDeg,
      ra:     target.ra,
      dec:    target.dec,
    });
    if (result.success) {
      setLog(`Added "${target.fullName || target.name}" to observation queue`);
      setActiveTab("todo");
      await refreshSeqState();
      if (typeof refreshNightPlan === "function") refreshNightPlan();
    } else {
      setLog({ success: false, error: result.error });
    }
  } catch (error) {
    setLog({ success: false, error: error.message });
  }
}

document.getElementById("add-to-todo").addEventListener("click", () => {
  if (lastTnsTarget) addToQueue(lastTnsTarget);
});

tabTodoBtn.addEventListener("click", () => setActiveTab("todo"));

document.getElementById("seq-run").addEventListener("click", async () => {
  if (!ocsConnected) {
    setLog({ success: false, error: `OCS not connected (${currentOcsHost()}) — cannot launch sequence. Check OCS connection in the Dome tab.` });
    return;
  }
  try {
    const result = await postJson("/api/sequence/run", { ...currentConfig(), ...seqConfig() });
    setLog(result);
    await refreshSeqState();
  } catch (error) {
    setLog({ success: false, error: error.message });
  }
});

document.getElementById("seq-abort").addEventListener("click", async () => {
  try {
    const result = await postJson("/api/sequence/abort", {});
    setLog(result);
  } catch (error) {
    setLog({ success: false, error: error.message });
  }
});

// ── Safety Override Modal ─────────────────────────────────────────────────────

function showSafetyModal(details) {
  const modal = document.getElementById("safety-modal");
  if (!modal) return;
  document.getElementById("safety-modal-message").textContent = details.message || "Safety condition not met.";
  document.getElementById("sm-ocs").textContent  = details.ocsReachable ? "✓ Reachable" : "✗ Unreachable";
  document.getElementById("sm-ocs").style.color  = details.ocsReachable ? "#4ade80" : "#f87171";
  document.getElementById("sm-safe").textContent = details.safe || "?";
  document.getElementById("sm-safe").style.color = /safe/i.test(details.safe || "") ? "#4ade80" : "#f87171";
  document.getElementById("sm-roof").textContent = details.roof || "?";
  document.getElementById("sm-roof").style.color = /open/i.test(details.roof || "") ? "#4ade80" : "#f87171";
  document.getElementById("sm-rain").textContent = details.rain || "?";
  modal.style.display = "flex";
}

function hideSafetyModal() {
  const modal = document.getElementById("safety-modal");
  if (modal) modal.style.display = "none";
}

document.getElementById("safety-override-btn").addEventListener("click", async () => {
  try {
    await postJson("/api/sequence/safety-override", {});
    hideSafetyModal();
  } catch (error) {
    setLog({ success: false, error: error.message });
  }
});

document.getElementById("safety-abort-btn").addEventListener("click", async () => {
  try {
    await postJson("/api/sequence/abort", {});
    hideSafetyModal();
  } catch (error) {
    setLog({ success: false, error: error.message });
  }
});

document.getElementById("seq-restart").addEventListener("click", async () => {
  try {
    const result = await postJson("/api/sequence/restart", {});
    seqLogRenderedCount = 0;
    document.getElementById("seq-log").innerHTML = "";
    setLog(result.success ? "Queue restarted — all targets marked pending." : result);
    await refreshSeqState();
  } catch (error) {
    setLog({ success: false, error: error.message });
  }
});

document.getElementById("seq-clear").addEventListener("click", async () => {
  try {
    const result = await postJson("/api/sequence/clear", {});
    seqLogRenderedCount = 0;
    document.getElementById("seq-log").innerHTML = "";
    await refreshSeqState();
    setLog(result);
  } catch (error) {
    setLog({ success: false, error: error.message });
  }
});

document.getElementById("seq-add-wait").addEventListener("click", async () => {
  const t = prompt("Wait until (HH:MM):", "23:30");
  if (!t || !t.match(/^\d{1,2}:\d{2}$/)) return;
  const hh = String(parseInt(t.split(":")[0])).padStart(2, "0");
  const mm = t.split(":")[1];
  await fetch("/api/sequence/queue/wait", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ waitTime: `${hh}:${mm}` }),
  });
  await refreshSeqState();
});

document.getElementById("seq-reset-af").addEventListener("click", async () => {
  try {
    const result = await postJson("/api/sequence/reset-af", {});
    setLog(result);
  } catch (error) {
    setLog({ success: false, error: error.message });
  }
});

document.getElementById("seq-next").addEventListener("click", async () => {
  try {
    const result = await postJson("/api/sequence/next", {});
    setLog(result);
    await refreshSeqState();
  } catch (error) {
    setLog({ success: false, error: error.message });
  }
});

document.getElementById("seq-manual-mode").addEventListener("change", async function () {
  try {
    await postJson("/api/sequence/manual-mode", { enabled: this.checked });
    await refreshSeqState();
  } catch (error) {
    setLog({ success: false, error: error.message });
  }
});

document.getElementById("toggle-obs-form").addEventListener("click", () => {
  const form = document.getElementById("obs-form");
  form.style.display = form.style.display === "none" ? "flex" : "none";
});

document.getElementById("apply-obs-location").addEventListener("click", () => {
  const lat = parseFloat(document.getElementById("obs-lat").value);
  const lon = parseFloat(document.getElementById("obs-lon").value);
  const elevation = parseFloat(document.getElementById("obs-elev").value) || 0;
  if (!isNaN(lat) && !isNaN(lon)) {
    observerLocation = { lat, lon, elevation };
    updateObsLocationDisplay();
    drawVisibility();
  }
});

// Redraw on window resize so the alt chart fills its container
window.addEventListener("resize", () => {
  if (lastTargetCoords.raDeg != null) drawVisibility();
});

tabPipelineBtn.addEventListener("click",  () => setActiveTab("pipeline"));
tabHistoryBtn.addEventListener("click",   () => { setActiveTab("history"); loadHistory(); });

// ── Target History ────────────────────────────────────────────────────────────

function renderHistory(targets) {
  const tbody = document.getElementById("history-tbody");
  const table = document.getElementById("history-table");
  const empty = document.getElementById("history-empty");
  const count = document.getElementById("history-count");
  tbody.innerHTML = "";

  count.textContent = targets.length ? `(${targets.length})` : "";

  if (!targets.length) {
    table.style.display = "none";
    empty.style.display = "block";
    return;
  }
  table.style.display = "";
  empty.style.display = "none";

  for (const t of targets) {
    const tr = document.createElement("tr");
    const raStr  = t.ra_deg  != null ? Number(t.ra_deg).toFixed(4)  : t.ra  || "–";
    const decStr = t.dec_deg != null ? Number(t.dec_deg).toFixed(4) : t.dec || "–";
    tr.innerHTML = `
      <td><strong>${t.name}</strong></td>
      <td>${t.type || "–"}</td>
      <td style="font-family:monospace;font-size:12px;">${raStr}</td>
      <td style="font-family:monospace;font-size:12px;">${decStr}</td>
      <td>
        <div class="buttons" style="margin:0;gap:5px;">
          <button class="btn-small btn-history-lookup" data-name="${t.name}" type="button">↻ Lookup</button>
          ${t.ra_deg != null ? `<button class="btn-small btn-history-queue" data-id="${t.id}" type="button">→ Queue</button>` : ""}
          <button class="btn-small btn-history-del" data-id="${t.id}" type="button" style="color:#fca5a5;border-color:#7f1d1d;">✕</button>
        </div>
      </td>
    `;
    tbody.appendChild(tr);
  }

  for (const btn of tbody.querySelectorAll(".btn-history-lookup")) {
    btn.addEventListener("click", () => {
      document.getElementById("target-name").value = btn.dataset.name;
      lookupTarget();
    });
  }

  for (const btn of tbody.querySelectorAll(".btn-history-queue")) {
    btn.addEventListener("click", () => {
      const t = targets.find(x => String(x.id) === btn.dataset.id);
      if (t) addToQueue({ name: t.name, fullName: t.name, raDeg: t.ra_deg, decDeg: t.dec_deg, ra: t.ra, dec: t.dec });
    });
  }

  for (const btn of tbody.querySelectorAll(".btn-history-del")) {
    btn.addEventListener("click", async () => {
      await fetch(`/api/targets/${btn.dataset.id}`, { method: "DELETE" });
      loadHistory();
    });
  }
}

async function loadHistory() {
  try {
    const res = await fetch("/api/targets");
    const data = await res.json();
    renderHistory(data.targets || []);
  } catch { /* ignore */ }
}

document.getElementById("history-clear").addEventListener("click", async () => {
  if (!confirm("Clear all history?")) return;
  try {
    const res = await fetch("/api/targets");
    const data = await res.json();
    await Promise.all((data.targets || []).map(t =>
      fetch(`/api/targets/${t.id}`, { method: "DELETE" })
    ));
    loadHistory();
  } catch { /* ignore */ }
});

// ── Pipeline ──────────────────────────────────────────────────────────────────

let pipelinePollInterval = null;
let activeLogJobId = null;
let logPollInterval = null;

const JOB_STATUS_CLASS = {
  queued:      "job-status-queued",
  scanning:    "job-status-running",
  splitting:   "job-status-running",
  calibrating: "job-status-running",
  solving:     "job-status-running",
  uploading:   "job-status-running",
  done:        "job-status-done",
  error:       "job-status-error",
};

const JOB_STATUS_STEPS = {
  queued:             "Queued…",
  scanning:           "Scanning FITS files…",
  splitting:          "Copying to work dir…",
  calibrating:        "Calibrating + stacking (Siril)…",
  solving:            "Plate solving…",
  uploading:          "Uploading to STDWeb…",
  selection_pending:  "⏸ Waiting for frame selection…",
  done:               "Complete",
  error:              "Error",
};

// ── NAS date / target picker ──────────────────────────────────────────────────

const _nasTargetCache = {};   // date → targets array, so re-selecting a date skips the scan
let _mpcMagLimit = 21;        // default; overwritten by loadMpcSettings()

async function loadNasDates() {
  const dateSel    = document.getElementById("pipeline-date-select");
  const targetSel  = document.getElementById("pipeline-target-select");
  const fitsInput  = document.getElementById("pipeline-fits-dir");
  const targetInput = document.getElementById("pipeline-target");

  // ── Step 1: load the date list (fast — no FITS reading) ──────────────────
  dateSel.disabled = true;
  dateSel.innerHTML = '<option value="">⏳ Loading nights…</option>';
  targetSel.style.display = "none";

  try {
    const res  = await fetch("/api/pipeline/nas-dates");
    const data = await res.json();
    const dates = data.dates || [];

    dateSel.innerHTML = '<option value="">— pick a night —</option>';
    for (const date of dates) {
      const opt = document.createElement("option");
      opt.value = date;
      opt.textContent = date;
      dateSel.appendChild(opt);
    }
    if (!dates.length) {
      dateSel.innerHTML = '<option value="">No data on NAS</option>';
    }
  } catch {
    dateSel.innerHTML = '<option value="">NAS unavailable</option>';
  }
  dateSel.disabled = false;

  function populateTargets(targets) {
    targetSel.innerHTML = '<option value="">— pick a target —</option>';
    if (!targets.length) {
      targetSel.innerHTML = '<option value="">No FITS data found</option>';
      targetSel.disabled = false;
      return;
    }
    const hasRealTargets = targets.some((t) => t.name !== "SNAPSHOT");
    for (const t of targets) {
      if (t.name === "SNAPSHOT" && hasRealTargets) continue;
      const opt = document.createElement("option");
      opt.value = t.path;
      opt.dataset.name = t.name;
      if (t.snapshot) opt.dataset.snapshot = "1";
      opt.textContent = t.name === "SNAPSHOT" ? "SNAPSHOT (unsorted)" : t.name;
      targetSel.appendChild(opt);
    }
    targetSel.disabled = false;

    // Auto-select if only one target
    if (targets.length === 1) {
      targetSel.value = targets[0].path;
      fitsInput.value = targets[0].path;
      if (targets[0].name !== "SNAPSHOT") targetInput.value = targets[0].name;
    }
  }

  // ── Step 2: when a night is chosen, scan that folder for targets ──────────
  dateSel.addEventListener("change", async () => {
    const date = dateSel.value;
    fitsInput.value = "";

    if (!date) {
      targetSel.style.display = "none";
      targetSel.disabled = false;
      return;
    }

    targetSel.style.display = "";

    // Use cached result if available (avoids re-scanning after queuing a job)
    if (_nasTargetCache[date]) {
      populateTargets(_nasTargetCache[date]);
      return;
    }

    targetSel.innerHTML = '<option value="">⏳ Scanning night…</option>';
    targetSel.disabled = true;

    let targets = [];
    try {
      const res  = await fetch(`/api/pipeline/nas-dates/${date}`);
      const data = await res.json();
      targets = data.targets || [];
      _nasTargetCache[date] = targets;   // cache for this session
    } catch {
      targetSel.innerHTML = '<option value="">Scan failed</option>';
      targetSel.disabled = false;
      return;
    }

    populateTargets(targets);
  });

  // ── Step 3: when a target is chosen, fill in the path/name fields ─────────
  targetSel.addEventListener("change", () => {
    if (!targetSel.value) return;
    fitsInput.value = targetSel.value;
    const selOpt = targetSel.selectedOptions[0];
    const name = selOpt?.dataset.name || "";
    if (name && name !== "SNAPSHOT") targetInput.value = name;
  });
}

// ── Job log modal ─────────────────────────────────────────────────────────────

function openLogModal(job_id) {
  activeLogJobId = job_id;
  document.getElementById("pipeline-log-modal").style.display = "block";
  fetchAndShowLog(job_id);
  logPollInterval = setInterval(() => fetchAndShowLog(activeLogJobId), 3000);
}

function closeLogModal() {
  document.getElementById("pipeline-log-modal").style.display = "none";
  clearInterval(logPollInterval);
  logPollInterval = null;
  activeLogJobId = null;
}

async function fetchAndShowLog(job_id) {
  try {
    const res  = await fetch(`/api/pipeline/job/${job_id}/log`);
    const data = await res.json();
    const pre  = document.getElementById("pipeline-log-content");
    const atBottom = pre.scrollHeight - pre.scrollTop - pre.clientHeight < 40;
    pre.textContent = data.log || data.error || "No log yet.";
    if (atBottom) pre.scrollTop = pre.scrollHeight;
  } catch { /* ignore */ }
}

document.getElementById("pipeline-log-close")
  .addEventListener("click", closeLogModal);

document.getElementById("pipeline-log-modal")
  .addEventListener("click", (e) => {
    if (e.target === e.currentTarget) closeLogModal();
  });

document.getElementById("pipeline-log-copy").addEventListener("click", () => {
  const text = document.getElementById("pipeline-log-content").textContent;
  const btn  = document.getElementById("pipeline-log-copy");
  const done = () => { btn.textContent = "✓ Copied"; setTimeout(() => { btn.textContent = "⎘ Copy"; }, 1500); };
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(done).catch(() => fallbackCopy(text, done));
  } else {
    fallbackCopy(text, done);
  }
});

// ── All-jobs photometry summary modal ────────────────────────────────────────

const _jobPhotModal      = document.getElementById("job-phot-modal");
const _jobPhotModalTitle = document.getElementById("job-phot-modal-title");
const _jobPhotModalBody  = document.getElementById("job-phot-modal-body");
const _jobPhotCopyBtn    = document.getElementById("job-phot-copy");
const _btnCopyAllMeasures = document.getElementById("btn-copy-all-measures");
_btnCopyAllMeasures.addEventListener("click", () => openAllPhotModal(_allJobs));

document.getElementById("job-phot-modal-close")
  .addEventListener("click", () => { _jobPhotModal.style.display = "none"; });
_jobPhotModal.addEventListener("click", (e) => {
  if (e.target === _jobPhotModal) _jobPhotModal.style.display = "none";
});

_jobPhotCopyBtn.addEventListener("click", () => {
  const text = _jobPhotCopyBtn.dataset.copyText || "";
  if (!text) return;
  const done = () => {
    _jobPhotCopyBtn.textContent = "✓ Copied";
    setTimeout(() => { _jobPhotCopyBtn.textContent = "⧉ Copy all"; }, 2000);
  };
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(done).catch(() => fallbackCopy(text, done));
  } else {
    fallbackCopy(text, done);
  }
});

async function openAllPhotModal(jobs) {
  // Collect all results across ALL jobs that have a STDWeb task
  const SKIP = new Set(["pending","running","uploading","error","uploaded"]);
  const photResults = [];
  for (const job of jobs) {
    for (const r of (job.results || [])) {
      if (r.stdweb_task_id && !SKIP.has(r.status)) {
        photResults.push({ ...r, jobTarget: job.target });
      }
    }
  }

  _jobPhotModalTitle.textContent = `📊 All measurements`;
  _jobPhotModalBody.innerHTML = `<div style="color:#9ca3af;font-size:12px;">Fetching…</div>`;
  _jobPhotCopyBtn.dataset.copyText = "";
  _jobPhotModal.style.display = "block";

  const allLines = [];
  const rows = [];

  for (const r of photResults) {
    const taskId = r.stdweb_task_id;

    if (!_photCache[taskId] || (!_photCache[taskId].sub && !_photCache[taskId].sub_ul)) {
      try {
        const resp = await fetch(`/api/stdweb/task/${taskId}/photometry`);
        _photCache[taskId] = await resp.json();
      } catch (e) {
        _photCache[taskId] = { error: e.message };
      }
    }

    const d      = _photCache[taskId];
    const mjd    = d.mjd    != null ? d.mjd    : "?";
    const filter = d.filter || r.filter || "?";
    const target = d.target || r.jobTarget || r.target || "";
    const tgt    = target ? target + ", " : "";
    const stack  = (r.n_frames && r.exposure) ? `, ${r.n_frames}×${r.exposure}s` : "";

    const lines = [];
    if (d.direct) lines.push(`${tgt}${mjd}, ${filter}${stack}, Aperture photometry, ${d.direct.mag.toFixed(2)}, ${d.direct.magerr.toFixed(2)}`);
    if (d.sub)    lines.push(`${tgt}${mjd}, ${filter}${stack}, Template substraction aperture photometry, ${d.sub.mag.toFixed(2)}, ${d.sub.magerr.toFixed(2)}`);
    else if (d.sub_ul) lines.push(`${tgt}${mjd}, ${filter}${stack}, Template substraction aperture photometry, >${d.sub_ul.ul.toFixed(2)}, –`);

    rows.push({ filter, taskId, lines, target, candidates: d.candidates || null, error: d.error || null });
    allLines.push(...lines);
  }

  if (!rows.length) {
    _jobPhotModalBody.innerHTML = `<span style="color:#9ca3af;">No measurements available yet.</span>`;
    return;
  }

  const html = rows.map(row => {
    const tgtLabel = row.target ? `<span style="color:#d1fae5;">${row.target}</span> · ` : "";
    if (row.error) {
      return `<div style="margin-bottom:12px;">
        <div style="font-size:11px;color:#6b7280;margin-bottom:4px;">${tgtLabel}${row.filter} · <a href="//${location.host.replace(/:\d+/,'')}:7000/tasks/${row.taskId}" target="_blank" style="color:#60a5fa;">#${row.taskId}</a></div>
        <div style="color:#f87171;font-size:11px;">Error: ${row.error}</div>
      </div>`;
    }
    if (!row.lines.length) {
      if (row.candidates && row.candidates.length) {
        const candLines = row.candidates.map((c, i) => {
          const filt = c.filter || row.filter;
          const magerr = c.magerr != null && !isNaN(c.magerr) ? ` ± ${c.magerr.toFixed(3)}` : "";
          const pos = (c.ra != null && c.dec != null) ? ` (RA ${c.ra.toFixed(5)}, Dec ${c.dec.toFixed(5)})` : "";
          return `Candidate ${i+1}${pos}: ${filt} = ${c.mag.toFixed(2)}${magerr}`;
        });
        const candEscaped = candLines.map(l => l.replace(/</g,"&lt;")).join("<br>");
        return `<div style="margin-bottom:14px;">
          <div style="font-size:11px;color:#6b7280;margin-bottom:4px;">${tgtLabel}${row.filter} · <a href="//${location.host.replace(/:\d+/,'')}:7000/tasks/${row.taskId}" target="_blank" style="color:#60a5fa;">#${row.taskId} ↗</a></div>
          <div style="color:#fbbf24;font-size:11px;margin-bottom:4px;">⚠ Target position unknown — ${row.candidates.length} transient candidate(s):</div>
          <pre style="margin:0;font-size:12px;color:#fde68a;font-family:monospace;line-height:1.7;">${candEscaped}</pre>
        </div>`;
      }
      return `<div style="margin-bottom:12px;">
        <div style="font-size:11px;color:#6b7280;margin-bottom:4px;">${tgtLabel}${row.filter} · <a href="//${location.host.replace(/:\d+/,'')}:7000/tasks/${row.taskId}" target="_blank" style="color:#60a5fa;">#${row.taskId}</a></div>
        <div style="color:#9ca3af;font-size:11px;">No magnitude found</div>
      </div>`;
    }
    const escaped = row.lines.map(l => l.replace(/</g,"&lt;")).join("<br>");
    return `<div style="margin-bottom:14px;">
      <div style="font-size:11px;color:#6b7280;margin-bottom:4px;">${tgtLabel}${row.filter} · <a href="//${location.host.replace(/:\d+/,'')}:7000/tasks/${row.taskId}" target="_blank" style="color:#60a5fa;">#${row.taskId} ↗</a></div>
      <pre style="margin:0;font-size:12px;color:#a7f3d0;font-family:monospace;line-height:1.7;">${escaped}</pre>
    </div>`;
  }).join("");

  _jobPhotModalBody.innerHTML = html;
  _jobPhotCopyBtn.dataset.copyText = allLines.join("\n");
}

function fallbackCopy(text, onSuccess) {
  const ta = document.createElement("textarea");
  ta.value = text;
  ta.style.cssText = "position:fixed;top:-9999px;left:-9999px;opacity:0;";
  document.body.appendChild(ta);
  ta.focus();
  ta.select();
  try { document.execCommand("copy"); onSuccess(); } catch { /* ignore */ }
  document.body.removeChild(ta);
}

// ── Render jobs ───────────────────────────────────────────────────────────────
// Photometry cache persists across re-renders so the data isn't re-fetched
// every 4-second poll and open rows stay open.
const _photCache = {};

// ── Steps definitions ─────────────────────────────────────────────────────────
const PIPELINE_STEPS = [
  { key: "calibrate",   label: "Calibrate & Stack",  statuses: ["calibrating"] },
  { key: "solve",       label: "Plate Solve",         statuses: ["solving"] },
  { key: "mpc_check",   label: "MPC Check",           statuses: [] },
  { key: "upload",      label: "Upload to STDWeb",    statuses: ["uploading", "uploaded"] },
  { key: "inspect",     label: "Inspect",             statuses: ["inspecting"] },
  { key: "photometry",  label: "Photometry",          statuses: ["photometry"] },
  { key: "subtraction", label: "Subtraction",         statuses: ["subtraction"] },
];

// Status order for deciding which steps are done vs pending
const STATUS_ORDER = [
  "queued","scanning","splitting","calibrating","solving",
  "uploading","uploaded","inspecting","photometry","subtraction","done","error"
];

// Map stdweb_state values to which step they belong to
const STDWEB_STATE_TO_STEP = {
  inspect_done: "inspect", inspect_failed: "inspect", inspect_error: "inspect",
  photometry_done: "photometry", photometry_failed: "photometry",
  subtraction_done: "subtraction", subtraction_failed: "subtraction",
};

function _getStepState(step, result) {
  if (result.status === "error") {
    // Determine which step caused the error
    let errorStep = null;
    if (result.stdweb_state) {
      // Map known terminal stdweb states
      errorStep = STDWEB_STATE_TO_STEP[result.stdweb_state] || null;
      // Also handle active states like "inspect_running" -> "inspect"
      if (!errorStep) {
        const key = Object.keys(STDWEB_STATE_TO_STEP).find(k =>
          result.stdweb_state.startsWith(k.split("_")[0]));
        if (key) errorStep = STDWEB_STATE_TO_STEP[key];
      }
    }
    if (!errorStep && result.error) {
      // Guess from error text
      if (/siril|calibrat|stack|registr/i.test(result.error)) errorStep = "calibrate";
      else if (/solv|astrometr/i.test(result.error)) errorStep = "solve";
      else if (/upload/i.test(result.error)) errorStep = "upload";
      else if (/inspect/i.test(result.error)) errorStep = "inspect";
      else if (/photom/i.test(result.error)) errorStep = "photometry";
      else if (/subtract/i.test(result.error)) errorStep = "subtraction";
    }
    const isError = errorStep === step.key;
    // Steps before the errored step are done
    const stepOrderIdx = PIPELINE_STEPS.findIndex(s => s.key === step.key);
    const errorStepIdx = PIPELINE_STEPS.findIndex(s => s.key === errorStep);
    const isDone = errorStepIdx > 0 && stepOrderIdx < errorStepIdx;
    return { isDone, isActive: false, isError };
  }

  const curIdx  = STATUS_ORDER.indexOf(result.status);
  const stepIdx = Math.max(...step.statuses.map(s => STATUS_ORDER.indexOf(s)));
  const isActive = step.statuses.includes(result.status);
  // mpc_check has no status of its own — mark done when mpc_objects is set,
  // or when the pipeline has advanced past "solving" (upload or later)
  if (step.key === "mpc_check") {
    const pastSolve = curIdx >= STATUS_ORDER.indexOf("uploading");
    const ran = result.mpc_objects != null;
    return { isDone: ran || pastSolve, isActive: false, isError: false };
  }
  const isDone   = !isActive && curIdx > stepIdx;
  return { isDone, isActive, isError: false };
}

// Steps that can be retried via STDWeb (no local re-run needed)
const STDWEB_RETRY_STEPS = new Set(["inspect", "photometry", "subtraction"]);

function _renderStepsRow(cell, r) {
  const stepsHtml = PIPELINE_STEPS.map(step => {
    const { isDone, isActive, isError } = _getStepState(step, r);
    let icon, color, cursor = "default";
    // All errored steps get a retry button; the step key is sent directly to
    // the /resume endpoint which handles all cases (local and STDWeb steps).
    const retryIcon = isError ? ` <span class="btn-retry-step" data-result-id="${r.id}" data-step="${step.key}"
        style="cursor:pointer;color:#f59e0b;font-size:12px;margin-left:3px;" title="Retry from this step">↺</span>` : "";

    if (isDone) {
      icon = "✓"; color = "#34d399";
    } else if (isActive) {
      icon = "⟳"; color = "#60a5fa";
    } else if (isError) {
      icon = "✗"; color = "#f87171"; cursor = "pointer";
    } else {
      icon = "○"; color = "#4b5563";
    }

    return `<span style="display:inline-flex;align-items:center;margin-right:14px;font-size:11px;color:${color};cursor:${cursor};">
      <span style="margin-right:3px;font-size:12px;">${icon}</span>${step.label}${retryIcon}
    </span>`;
  }).join("");

  cell.innerHTML = stepsHtml;
}

function _renderPhotRow(photRow, d, result) {
  const cell = photRow.querySelector("td");
  if (d.error) {
    cell.innerHTML = `<span style="color:#f87171;">Error: ${d.error}</span>`;
    return;
  }

  const mjd    = d.mjd    != null ? d.mjd    : "?";
  const filter = d.filter || "?";
  const target = d.target || "";
  const stack  = (result?.n_frames && result?.exposure)
    ? `, ${result.n_frames}×${result.exposure}s` : "";

  const lines = [];
  if (d.direct) {
    lines.push(`${target ? target + ", " : ""}${mjd}, ${filter}${stack}, Aperture photometry, ${d.direct.mag.toFixed(2)}, ${d.direct.magerr.toFixed(2)}`);
  }
  if (d.sub) {
    lines.push(`${target ? target + ", " : ""}${mjd}, ${filter}${stack}, Template substraction aperture photometry, ${d.sub.mag.toFixed(2)}, ${d.sub.magerr.toFixed(2)}`);
  } else if (d.sub_ul) {
    lines.push(`${target ? target + ", " : ""}${mjd}, ${filter}${stack}, Template substraction aperture photometry, >${d.sub_ul.ul.toFixed(2)}, –`);
  }

  // Fallback: show transient candidates detected in difference image
  if (!lines.length && d.candidates && d.candidates.length) {
    const candLines = d.candidates.map((c, i) => {
      const filt = c.filter || filter;
      const magerr = c.magerr != null && !isNaN(c.magerr) ? ` ± ${c.magerr.toFixed(3)}` : "";
      const pos = (c.ra != null && c.dec != null) ? ` (RA ${c.ra.toFixed(5)}, Dec ${c.dec.toFixed(5)})` : "";
      return `Candidate ${i+1}${pos}: ${filt} = ${c.mag.toFixed(2)}${magerr}`;
    });
    const candEscaped = candLines.map(l => l.replace(/</g, "&lt;")).join("<br>");
    const candText    = candLines.join("\n");
    cell.innerHTML = `
      <div style="display:flex;align-items:flex-start;gap:10px;flex-direction:column;">
        <span style="color:#fbbf24;font-size:11px;">⚠ Target position unknown — showing ${d.candidates.length} transient candidate(s) from subtraction:</span>
        <div style="display:flex;align-items:flex-start;gap:10px;">
          <pre style="margin:0;font-size:11px;color:#fde68a;font-family:monospace;line-height:1.6;">${candEscaped}</pre>
          <button class="btn-phot-copy" title="Copy to clipboard"
                  style="flex-shrink:0;background:none;border:1px solid #78350f;color:#fbbf24;
                         border-radius:4px;padding:2px 7px;cursor:pointer;font-size:14px;line-height:1;"
                  data-text="${candText.replace(/"/g, "&quot;")}">⧉</button>
        </div>
      </div>`;
    cell.querySelector(".btn-phot-copy").addEventListener("click", (e) => {
      const copyBtn = e.currentTarget;
      navigator.clipboard.writeText(copyBtn.dataset.text).then(() => {
        copyBtn.textContent = "✓";
        setTimeout(() => { copyBtn.textContent = "⧉"; }, 1500);
      });
    });
    return;
  }

  if (!lines.length) {
    cell.innerHTML = `<span style="color:#9ca3af;font-size:11px;">No magnitude measurements found</span>`;
    return;
  }

  const text    = lines.join("\n");
  const escaped = lines.map(l => l.replace(/</g, "&lt;")).join("<br>");

  cell.innerHTML = `
    <div style="display:flex;align-items:flex-start;gap:10px;">
      <pre style="margin:0;font-size:11px;color:#a7f3d0;font-family:monospace;line-height:1.6;">${escaped}</pre>
      <button class="btn-phot-copy" title="Copy to clipboard"
              style="flex-shrink:0;background:none;border:1px solid #065f46;color:#34d399;
                     border-radius:4px;padding:2px 7px;cursor:pointer;font-size:14px;line-height:1;"
              data-text="${text.replace(/"/g, "&quot;")}">⧉</button>
    </div>`;

  cell.querySelector(".btn-phot-copy").addEventListener("click", (e) => {
    const copyBtn = e.currentTarget;
    const done = () => {
      copyBtn.textContent = "✓";
      setTimeout(() => { copyBtn.textContent = "⧉"; }, 1500);
    };
    const fail = () => {
      copyBtn.textContent = "✗";
      setTimeout(() => { copyBtn.textContent = "⧉"; }, 1500);
    };
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(copyBtn.dataset.text).then(done).catch(() => fallbackCopy(copyBtn.dataset.text, done));
    } else {
      fallbackCopy(copyBtn.dataset.text, done);
    }
  });
}

function renderPipelineJobs(jobs) {
  const tbody   = document.getElementById("pipeline-tbody");
  const table   = document.getElementById("pipeline-table");
  const empty   = document.getElementById("pipeline-empty");
  const statBar = document.getElementById("pipeline-status-bar");
  const statTxt = document.getElementById("pipeline-status-text");

  // Preserve expanded state and unchecked boxes across re-renders
  const openFrameRows = new Set(
    Array.from(tbody.querySelectorAll("tr[data-frames-for]"))
      .filter(r => r.style.display !== "none")
      .map(r => r.dataset.framesFor)
  );
  const openPhotRows = new Set(
    Array.from(tbody.querySelectorAll("tr[data-phot-for]"))
      .filter(r => r.style.display !== "none")
      .map(r => r.dataset.photFor)
  );
  const openStepsRows = new Set(
    Array.from(tbody.querySelectorAll("tr[data-steps-for]"))
      .filter(r => r.style.display !== "none")
      .map(r => r.dataset.stepsFor)
  );
  const openMpcRows = new Set(
    Array.from(tbody.querySelectorAll("tr[data-mpc-for]"))
      .filter(r => r.style.display !== "none")
      .map(r => r.dataset.mpcFor)
  );
  // Map: resultId → Set of source paths that are UNCHECKED
  const uncheckedSources = {};
  for (const previewRow of tbody.querySelectorAll("tr[data-frames-for]")) {
    const rid = previewRow.dataset.framesFor;
    uncheckedSources[rid] = new Set(
      Array.from(previewRow.querySelectorAll(".frame-check:not(:checked)"))
        .map(c => c.dataset.source)
    );
  }

  tbody.innerHTML = "";

  // Status bar: show active job + queue depth
  const activeJobs = jobs.filter((j) => !["done", "error", "queued"].includes(j.status));
  const queuedJobs = jobs.filter((j) => j.status === "queued");
  const running = activeJobs[0];
  if (running || queuedJobs.length) {
    statBar.style.display = "block";
    const isPending = running?.status === "selection_pending";
    const step = JOB_STATUS_STEPS[running?.status] || running?.status || "Queued…";
    const queueSuffix = queuedJobs.length
      ? ` — ${queuedJobs.length} job${queuedJobs.length > 1 ? "s" : ""} waiting`
      : "";
    statTxt.textContent = running
      ? `Job #${running.id} — ${step}${queueSuffix}`
      : `${queuedJobs.length} job${queuedJobs.length > 1 ? "s" : ""} queued…`;
    statBar.style.background = isPending ? "#1e3a2f" : "#1e3a5f";
    document.getElementById("pipeline-log-btn").dataset.jobId = running?.id || "";
    document.getElementById("pipeline-log-btn").style.display = running ? "" : "none";
  } else {
    statBar.style.display = "none";
  }

  if (!jobs.length) {
    table.style.display = "none";
    empty.style.display = "block";
    _btnCopyAllMeasures.style.display = "none";
    return;
  }
  table.style.display = "";
  empty.style.display = "none";

  // Show global 📊 button if any result has STDWeb measurements
  const SKIP_PHOT = new Set(["pending","running","uploading","error","uploaded"]);
  const anyPhot = jobs.some(j => (j.results||[]).some(r => r.stdweb_task_id && !SKIP_PHOT.has(r.status)));
  _btnCopyAllMeasures.style.display = anyPhot ? "" : "none";

  const RESULT_STATUS_CLASS = {
    processing:  "job-status-running",
    splitting:   "job-status-running",
    solving:     "job-status-running",
    uploading:   "job-status-running",
    uploaded:    "job-status-running",
    inspecting:  "job-status-running",
    photometry:  "job-status-running",
    subtraction: "job-status-running",
    done:        "job-status-done",
    error:       "job-status-error",
  };
  const RESULT_STATUS_LABEL = {
    processing:  "Calibrating…",
    splitting:   "Splitting…",
    solving:     "Solving…",
    uploading:   "Uploading…",
    uploaded:    "Uploaded",
    inspecting:  "Inspecting…",
    photometry:  "Photometry…",
    subtraction: "Subtraction…",
    done:        "Done",
    error:       "Error",
  };

  for (const j of jobs) {
    const cls        = JOB_STATUS_CLASS[j.status] || "job-status-queued";
    const date       = j.created_at ? j.created_at.slice(0, 16).replace("T", " ") : "–";
    const errorTitle = j.error ? ` title="${j.error.replace(/"/g, "&quot;")}"` : "";
    const folderShort = (j.fits_dir || "–").split("/").slice(-2).join("/");
    const results    = j.results || [];
    const hasResults = results.length > 0;

    // ── Job header row ────────────────────────────────────────────────────────

    const tr = document.createElement("tr");
    tr.style.borderTop = "1px solid #374151";
    tr.innerHTML = `
      <td style="font-size:12px;white-space:nowrap;">${date}</td>
      <td><strong>${j.target || "–"}</strong></td>
      <td style="font-size:11px;color:#6b7280;" title="${j.fits_dir || ""}">${folderShort}</td>
      <td><span class="job-status ${cls}"${errorTitle}>${JOB_STATUS_STEPS[j.status] || j.status}</span></td>
      <td style="color:#6b7280;font-size:12px;">${hasResults ? results.length + " filter" + (results.length > 1 ? "s" : "") : "–"}</td>
      <td style="white-space:nowrap;">
        <button class="btn-small btn-job-log" data-id="${j.id}" type="button">Log</button>
        ${j.status === "error" ? `<button class="btn-small btn-rerun-job" data-id="${j.id}" type="button"
                style="margin-left:4px;" title="Re-run with same parameters">↺ Re-run</button>` : ""}
        ${j.status === "selection_pending" ? `<button class="btn-small btn-select-frames" data-id="${j.id}" type="button"
                style="margin-left:4px;background:#1d4ed8;color:#fff;font-weight:600;"
                title="Review frames and choose which ones to stack">🖼 Review & select</button>` : ""}
        <button class="btn-small btn-del-job" data-id="${j.id}" type="button"
                style="color:#fca5a5;border-color:#7f1d1d;margin-left:4px;">✕</button>
      </td>
    `;
    tbody.appendChild(tr);

    // ── Per-filter result rows ─────────────────────────────────────────────────
    for (const r of results) {
      const rcls    = RESULT_STATUS_CLASS[r.status] || "job-status-queued";
      const rlabel  = RESULT_STATUS_LABEL[r.status] || r.status;
      const retitle = r.error ? ` title="${r.error.replace(/"/g, "&quot;")}"` : "";
      const link    = r.stdweb_url
        ? `<a href="${r.stdweb_url}" target="_blank" class="tns-link" style="font-size:12px;">
             #${r.stdweb_task_id} ↗</a>` : "–";

      // Parse frame previews (stored as JSON list of {preview, source} objects)
      let framePreviews = [];
      try { framePreviews = r.frame_previews ? JSON.parse(r.frame_previews) : []; } catch { /* ignore */ }
      const hasFrames = framePreviews.length > 0;
      const framesBtn = hasFrames
        ? `<button class="btn-small btn-frames-toggle" data-result-id="${r.id}"
                   style="margin-left:4px;color:#fbbf24;border-color:#78350f;" type="button">
             🖼 ${framePreviews.length} frame${framePreviews.length > 1 ? "s" : ""}
           </button>` : "";

      // Show phot button when STDWeb has processed past 'uploaded'
      const canFetchPhot = r.stdweb_task_id &&
        !["pending", "running", "uploading", "error", "uploaded"].includes(r.status);
      const photBtn = canFetchPhot
        ? `<button class="btn-small btn-phot-toggle" data-task-id="${r.stdweb_task_id}" data-result-id="${r.id}"
                   data-result-status="${r.status}"
                   style="margin-left:4px;color:#34d399;border-color:#065f46;" type="button"
                   title="Fetch photometry from STDWeb">📊</button>` : "";

      // FITS download + preview — shown whenever the stack file exists on disk (even if STDWeb failed)
      const stackDownload = r.stack_url
        ? `<a href="${r.stack_url}"
              download
              class="btn-small"
              style="margin-left:4px;color:#a78bfa;border:1px solid #5b21b6;border-radius:4px;
                     padding:2px 6px;font-size:12px;text-decoration:none;display:inline-block;"
              title="Download stacked FITS (${r.obs_date} · ${r.target} · ${r.filter} · ${r.exposure}s)">⭐ FITS</a>` : "";
      const previewThumb = r.preview_url
        ? `<a href="${r.preview_url}" target="_blank" style="margin-left:6px;display:inline-block;vertical-align:middle;"
              title="Preview stacked image">
             <img src="${r.preview_url}" style="height:28px;width:28px;object-fit:cover;border-radius:3px;
                        border:1px solid #374151;vertical-align:middle;" />
           </a>` : "";
      const fitsBtn = stackDownload + previewThumb;

      // MPC asteroid badge — shown when SkyBoT found objects in the field
      let mpcObjects = [];
      try { mpcObjects = r.mpc_objects ? JSON.parse(r.mpc_objects) : []; } catch { /* ignore */ }
      const mpcVisible = mpcObjects.filter(o => o.mag == null || o.mag <= _mpcMagLimit);
      const nNeo = mpcVisible.filter(o => {
        const cls = (o.class || "").toUpperCase().replace(/[>\s]/g, "");
        return ["AMOR","APOLLO","ATEN","IEO","AMO","APO","ATE","NEA","NEO"].some(c => cls.includes(c));
      }).length;
      const mpcBadge = mpcVisible.length > 0
        ? `<button class="btn-small btn-mpc-toggle" data-result-id="${r.id}"
                   style="margin-left:4px;color:${nNeo > 0 ? '#ff9500' : '#00e5ff'};
                          border-color:${nNeo > 0 ? '#7c3f00' : '#004d5c'};" type="button"
                   title="Solar system objects in field (SkyBoT) — V ≤ ${_mpcMagLimit}">
             &#x25CE; ${mpcVisible.length}${nNeo > 0 ? ' (' + nNeo + ' NEO)' : ''}
           </button>` : "";

      const rtr = document.createElement("tr");
      rtr.dataset.resultId = r.id;
      rtr.style.background = "#0d1929";
      rtr.innerHTML = `
        <td></td>
        <td style="font-size:12px;padding-left:12px;color:#93c5fd;">
          <button class="btn-steps-toggle" data-result-id="${r.id}"
                  style="background:none;border:none;color:#4b5563;cursor:pointer;
                         font-size:10px;padding:0 4px 0 0;line-height:1;vertical-align:middle;"
                  title="Show/hide steps">▶</button>${r.target ? `<span style="color:#a78bfa;">${r.target}</span> · ` : ""}${r.filter || "?"} · ${r.n_frames || "?"}×${r.exposure || "?"}s
        </td>
        <td></td>
        <td><span class="job-status ${rcls}"${retitle} style="font-size:11px;">${rlabel}</span></td>
        <td>${link}${framesBtn}${photBtn}${fitsBtn}${mpcBadge}</td>
        <td></td>
      `;
      tbody.appendChild(rtr);

      // ── Steps sub-row (hidden until triangle clicked) ─────────────────────
      const stepsRow = document.createElement("tr");
      stepsRow.dataset.stepsFor = r.id;
      stepsRow.style.cssText = "background:#060c18;display:none;";
      stepsRow.innerHTML = `<td colspan="6" style="padding:4px 16px 8px 36px;" class="steps-cell"></td>`;
      tbody.appendChild(stepsRow);
      _renderStepsRow(stepsRow.querySelector(".steps-cell"), r);

      // Restore open state
      if (openStepsRows.has(String(r.id))) {
        stepsRow.style.display = "";
        const tri = rtr.querySelector(".btn-steps-toggle");
        if (tri) { tri.textContent = "▼"; tri.style.color = "#60a5fa"; }
      }

      // Photometry row (hidden until fetched, restored if was open before re-render)
      const photRow = document.createElement("tr");
      photRow.dataset.photFor = r.id;
      photRow.style.cssText = "background:#060f1e;display:none;";
      photRow.innerHTML = `<td colspan="6" style="padding:6px 16px 10px 36px;"></td>`;
      tbody.appendChild(photRow);

      // Restore previously open phot row using cached data
      if (openPhotRows.has(String(r.id)) && r.stdweb_task_id && _photCache[r.stdweb_task_id]) {
        _renderPhotRow(photRow, _photCache[r.stdweb_task_id], r);
        photRow.style.display = "";
        const toggle = tbody.querySelector(`.btn-phot-toggle[data-result-id="${r.id}"]`);
        if (toggle) toggle.style.opacity = "1";
      }

      // MPC asteroid row (hidden until badge clicked)
      if (mpcObjects.length > 0) {
        const mpcRow = document.createElement("tr");
        mpcRow.dataset.mpcFor = r.id;
        mpcRow.style.cssText = "background:#040d18;display:none;";

        const mpcPreviewHtml = r.mpc_preview_url
          ? `<div style="display:flex;gap:12px;align-items:flex-start;">
               <div>
                 <div style="font-size:11px;color:#6b7280;margin-bottom:4px;">Annotated preview</div>
                 <a href="${r.mpc_preview_url}" target="_blank">
                   <img src="${r.mpc_preview_url}" alt="MPC preview"
                        style="max-height:260px;max-width:480px;border-radius:4px;
                               border:1px solid #1e3a5f;display:block;" />
                 </a>
               </div>
               <div style="min-width:220px;">{TABLE}</div>
             </div>`
          : `<div>{TABLE}</div>`;

        const tableRows = [...mpcVisible]
          .sort((a, b) => (a.mag ?? 99) - (b.mag ?? 99))
          .map(o => {
          const isNeo = ["AMOR","APOLLO","ATEN","IEO","AMO","APO","ATE","NEA","NEO"]
            .some(c => (o.class||"").toUpperCase().replace(/[>\s]/g,"").includes(c));
          const magStr = o.mag != null ? o.mag.toFixed(1) : "—";
          const uncStr = o.uncertainty_arcsec != null ? o.uncertainty_arcsec.toFixed(0) + '"' : "—";
          const nameColor = isNeo ? "#ff9500" : "#00e5ff";
          const tag = isNeo
            ? `<span style="color:#ff9500;font-size:10px;border:1px solid #7c3f00;
                            border-radius:3px;padding:0 3px;margin-left:4px;">NEO</span>` : "";
          return `<tr style="border-bottom:1px solid #0d2040;">
            <td style="padding:3px 6px;color:${nameColor};font-weight:${isNeo?'bold':'normal'};">${o.name}${tag}</td>
            <td style="padding:3px 6px;color:#6b7280;font-size:11px;">${o.class || "?"}</td>
            <td style="padding:3px 6px;color:#9ca3af;font-size:11px;">V=${magStr}</td>
            <td style="padding:3px 6px;color:#6b7280;font-size:11px;">±${uncStr}</td>
          </tr>`;
        }).join("");

        const tableHtml = `<table style="border-collapse:collapse;font-size:12px;width:100%;">
          <thead>
            <tr style="color:#4b5563;font-size:10px;border-bottom:1px solid #1e3a5f;">
              <th style="padding:2px 6px;text-align:left;">Object</th>
              <th style="padding:2px 6px;text-align:left;">Type</th>
              <th style="padding:2px 6px;text-align:left;">Mag</th>
              <th style="padding:2px 6px;text-align:left;">Uncert.</th>
            </tr>
          </thead>
          <tbody>${tableRows}</tbody>
        </table>`;

        mpcRow.innerHTML = `<td colspan="6" style="padding:6px 16px 12px 36px;">
          <div style="font-size:11px;color:#6b7280;margin-bottom:6px;">
            Solar system objects in field (SkyBoT) — ${mpcVisible.length} object(s) shown (V ≤ ${_mpcMagLimit}${mpcObjects.length > mpcVisible.length ? ', ' + (mpcObjects.length - mpcVisible.length) + ' filtered' : ''})
          </div>
          ${mpcPreviewHtml.replace("{TABLE}", tableHtml)}
        </td>`;
        tbody.appendChild(mpcRow);

        // Restore open state after re-render
        if (openMpcRows.has(String(r.id))) {
          mpcRow.style.display = "";
          const badge = rtr.querySelector(`.btn-mpc-toggle[data-result-id="${r.id}"]`);
          if (badge) badge.style.opacity = "1";
        }
      }

      // Hidden frame preview row — checkboxes + re-process button
      if (hasFrames) {
        const fitsDir = j.fits_dir || "";
        const jobTarget = j.target || "";

        const previewRow = document.createElement("tr");
        previewRow.dataset.framesFor = r.id;
        previewRow.style.cssText = "background:#060f1e;display:none;";

        const framesHtml = framePreviews.map((fp, i) => {
          // fp may be a {preview, source} object or a plain string (legacy)
          const previewPath = typeof fp === "object" ? fp.preview : fp;
          const sourcePath  = typeof fp === "object" ? fp.source : "";
          const previewUrl  = "/data/" + previewPath.split("/").map(encodeURIComponent).join("/");
          const fname       = sourcePath ? sourcePath.split("/").pop() : `Frame ${i + 1}`;
          return `
            <div style="text-align:center;">
              <label style="cursor:pointer;display:block;">
                <input type="checkbox" class="frame-check" checked
                       data-source="${sourcePath.replace(/"/g,'&quot;')}"
                       style="margin-bottom:4px;" />
                <div class="frame-thumb-wrap"
                     style="width:160px;height:160px;overflow:hidden;border-radius:4px;
                            border:1px solid #1e3a5f;display:block;cursor:crosshair;">
                  <img src="${previewUrl}" alt="Frame ${i+1}"
                       class="frame-zoom-img"
                       style="width:160px;height:160px;display:block;
                              transform-origin:center center;
                              will-change:transform;"
                       onclick="event.preventDefault();window.open(this.src,'_blank')" />
                </div>
                <div style="font-size:10px;color:#6b7280;margin-top:3px;word-break:break-all;">${fname}</div>
              </label>
            </div>`;
        }).join("");

        previewRow.innerHTML = `
          <td colspan="6" style="padding:10px 16px 16px;">
            <div style="font-size:11px;color:#f59e0b;margin-bottom:10px;">
              ⚠ Registration failed — no stars found. Check frames below for clouds or tracking issues.
              Uncheck bad frames, then click <strong>Re-process selected</strong>.
            </div>
            <div class="frame-grid" style="display:flex;flex-wrap:wrap;gap:10px;margin-bottom:12px;">
              ${framesHtml}
            </div>
            <div style="display:flex;align-items:center;gap:10px;">
              <button class="btn-small btn-reprocess"
                      data-fits-dir="${fitsDir.replace(/"/g,'&quot;')}"
                      data-target="${jobTarget.replace(/"/g,'&quot;')}"
                      data-row-id="${previewRow.dataset ? '' : ''}"
                      style="background:#1d4ed8;border-color:#3b82f6;color:#fff;" type="button">
                ▶ Re-process selected
              </button>
              <span class="reprocess-status" style="font-size:11px;color:#6b7280;"></span>
            </div>
          </td>`;

        tbody.appendChild(previewRow);

        // Restore expanded state after row is in the DOM
        if (openFrameRows.has(String(r.id))) {
          previewRow.style.display = "";
          const toggle = tbody.querySelector(`.btn-frames-toggle[data-result-id="${r.id}"]`);
          if (toggle) toggle.style.opacity = "1";
        }

        // Restore unchecked checkboxes
        const prevUnchecked = uncheckedSources[String(r.id)];
        if (prevUnchecked && prevUnchecked.size > 0) {
          for (const cb of previewRow.querySelectorAll(".frame-check")) {
            if (prevUnchecked.has(cb.dataset.source)) cb.checked = false;
          }
        }
      }
    }
  }

  for (const btn of tbody.querySelectorAll(".btn-del-job")) {
    btn.addEventListener("click", async () => {
      await fetch(`/api/pipeline/jobs/${btn.dataset.id}`, { method: "DELETE" });
      loadPipelineJobs();
    });
  }
  for (const btn of tbody.querySelectorAll(".btn-rerun-job")) {
    btn.addEventListener("click", async () => {
      btn.disabled = true;
      btn.textContent = "↺ …";
      try {
        const r = await fetch(`/api/pipeline/rerun/${btn.dataset.id}`, { method: "POST" });
        const data = await r.json();
        if (data.success) { loadPipelineJobs(); }
        else { btn.textContent = "✗ Failed"; btn.disabled = false; }
      } catch { btn.textContent = "✗ Error"; btn.disabled = false; }
    });
  }
  for (const btn of tbody.querySelectorAll(".btn-job-log")) {
    btn.addEventListener("click", () => openLogModal(btn.dataset.id));
  }
  for (const btn of tbody.querySelectorAll(".btn-select-frames")) {
    btn.addEventListener("click", () => openFrameSelectionModal(Number(btn.dataset.id)));
  }
  for (const btn of tbody.querySelectorAll(".btn-frames-toggle")) {
    btn.addEventListener("click", () => {
      const rid = btn.dataset.resultId;
      const row = tbody.querySelector(`tr[data-frames-for="${rid}"]`);
      if (!row) return;
      const hidden = row.style.display === "none";
      row.style.display = hidden ? "" : "none";
      btn.style.opacity = hidden ? "1" : "0.5";
    });
  }

  // ── Photometry fetch button ───────────────────────────────────────────────
  for (const btn of tbody.querySelectorAll(".btn-phot-toggle")) {
    btn.addEventListener("click", async () => {
      const rid     = btn.dataset.resultId;
      const taskId  = btn.dataset.taskId;
      const photRow = tbody.querySelector(`tr[data-phot-for="${rid}"]`);
      if (!photRow) return;

      // Toggle: hide if already visible
      if (photRow.style.display !== "none") {
        photRow.style.display = "none";
        btn.style.opacity = "0.5";
        return;
      }

      // Always force a live re-fetch from STDWeb when the user clicks the
      // measure button — they may have just changed photometry params on
      // STDWeb (color term, S/N threshold, catalog, …) and want the new
      // numbers, not whatever we cached earlier.
      delete _photCache[taskId];

      const cell = photRow.querySelector("td");
      btn.disabled = true;
      btn.textContent = "⏳";
      cell.textContent = "Fetching live from STDWeb…";
      photRow.style.display = "";

      try {
        const resp = await fetch(`/api/stdweb/task/${taskId}/photometry?refresh=1`);
        _photCache[taskId] = await resp.json();
      } catch (e) {
        cell.textContent = `Error: ${e.message}`;
        btn.disabled = false; btn.textContent = "📊";
        return;
      }
      btn.disabled = false; btn.textContent = "📊";

      const resultObj = jobs.flatMap(j => j.results || []).find(x => String(x.id) === String(rid));
      _renderPhotRow(photRow, _photCache[taskId], resultObj);
      photRow.style.display = "";
      btn.style.opacity = "1";
    });
  }

  // ── MPC asteroid badge toggle ─────────────────────────────────────────────
  for (const btn of tbody.querySelectorAll(".btn-mpc-toggle")) {
    btn.addEventListener("click", () => {
      const rid = btn.dataset.resultId;
      const row = tbody.querySelector(`tr[data-mpc-for="${rid}"]`);
      if (!row) return;
      const hidden = row.style.display === "none";
      row.style.display = hidden ? "" : "none";
      btn.style.opacity = hidden ? "1" : "0.5";
    });
  }

  // ── Steps triangle toggle ──────────────────────────────────────────────────
  for (const btn of tbody.querySelectorAll(".btn-steps-toggle")) {
    btn.addEventListener("click", (e) => {
      e.stopPropagation();
      const rid = btn.dataset.resultId;
      const row = tbody.querySelector(`tr[data-steps-for="${rid}"]`);
      if (!row) return;
      const hidden = row.style.display === "none";
      row.style.display = hidden ? "" : "none";
      btn.textContent = hidden ? "▼" : "▶";
      btn.style.color = hidden ? "#60a5fa" : "#4b5563";
    });
  }

  // ── Retry step button ──────────────────────────────────────────────────────
  for (const btn of tbody.querySelectorAll(".btn-retry-step")) {
    btn.addEventListener("click", async (e) => {
      e.stopPropagation();
      const rid  = btn.dataset.resultId;
      const step = btn.dataset.step;
      const orig = btn.textContent;
      btn.textContent = "⏳";
      btn.style.pointerEvents = "none";
      try {
        const resp = await fetch(`/api/pipeline/results/${rid}/retry-step`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ step }),
        });
        const data = await resp.json();
        if (data.success) {
          btn.textContent = "✓";
          btn.style.color = "#34d399";
          setTimeout(loadPipelineJobs, 1500);
        } else {
          btn.textContent = "✗";
          btn.style.color = "#f87171";
          btn.title = data.error || "Failed";
          btn.style.pointerEvents = "auto";
        }
      } catch (err) {
        btn.textContent = orig;
        btn.style.pointerEvents = "auto";
        console.error("retry-step error:", err);
      }
    });
  }

  for (const btn of tbody.querySelectorAll(".btn-reprocess")) {
    btn.addEventListener("click", async () => {
      const row = btn.closest("tr");
      const checks = row.querySelectorAll(".frame-check:checked");
      const selected = Array.from(checks)
        .map(c => c.dataset.source)
        .filter(s => s);

      if (!selected.length) {
        row.querySelector(".reprocess-status").textContent = "⚠ No frames selected";
        return;
      }

      const fits_dir = btn.dataset.fitsDir;
      const target   = btn.dataset.target;
      btn.disabled = true;
      btn.textContent = "Starting…";
      row.querySelector(".reprocess-status").textContent = `Queuing ${selected.length} frame(s)…`;

      try {
        const result = await postJson("/api/pipeline/trigger", {
          fits_dir,
          target,
          selected_files: selected,
        });
        if (result.success) {
          row.querySelector(".reprocess-status").textContent =
            `✓ Job #${result.id} queued with ${selected.length} frame(s)`;
          loadPipelineJobs();
        } else {
          row.querySelector(".reprocess-status").textContent = `✗ ${result.error}`;
          btn.disabled = false;
          btn.textContent = "▶ Re-process selected";
        }
      } catch (err) {
        row.querySelector(".reprocess-status").textContent = `✗ ${err.message}`;
        btn.disabled = false;
        btn.textContent = "▶ Re-process selected";
      }
    });
  }
}

// ── Frame in-place magnifier ───────────────────────────────────────────────
// Classic loupe: thumbnail stays still, a square magnifier window appears
// next to the cursor showing a zoomed crop of the image.

(function setupFrameZoom() {
  const ZOOM = 5;

  document.addEventListener("mouseenter", (e) => {
    const img = e.target;
    if (!img.classList || !img.classList.contains("frame-zoom-img")) return;
    const wrap = img.closest(".frame-thumb-wrap") || img.parentElement;
    wrap.style.position = "relative";

    const loupe = document.createElement("div");
    loupe.className = "_frame-loupe";
    Object.assign(loupe.style, {
      position: "fixed",
      width: "240px",
      height: "240px",
      border: "2px solid #4fc3f7",
      borderRadius: "4px",
      background: `url("${img.src}") no-repeat`,
      backgroundSize: `${img.clientWidth * ZOOM}px ${img.clientHeight * ZOOM}px`,
      pointerEvents: "none",
      zIndex: "9999",
      boxShadow: "0 4px 20px rgba(0,0,0,0.6)",
      display: "none",
    });
    document.body.appendChild(loupe);
    img._loupe = loupe;

    img.style.cursor = "crosshair";
  }, true);

  document.addEventListener("mousemove", (e) => {
    const img = e.target.closest ? e.target.closest(".frame-zoom-img") : null;
    if (!img || !img._loupe) return;
    const loupe = img._loupe;
    const rect = img.getBoundingClientRect();

    const cx = e.clientX - rect.left;
    const cy = e.clientY - rect.top;
    if (cx < 0 || cy < 0 || cx > rect.width || cy > rect.height) {
      loupe.style.display = "none";
      return;
    }

    const bgX = -(cx * ZOOM - 120);
    const bgY = -(cy * ZOOM - 120);
    loupe.style.backgroundPosition = `${bgX}px ${bgY}px`;
    loupe.style.backgroundSize = `${rect.width * ZOOM}px ${rect.height * ZOOM}px`;
    loupe.style.display = "block";

    // Position loupe to the right of cursor, or left if near edge
    let lx = e.clientX + 20;
    let ly = e.clientY - 120;
    if (lx + 250 > window.innerWidth) lx = e.clientX - 260;
    if (ly < 10) ly = 10;
    if (ly + 240 > window.innerHeight) ly = window.innerHeight - 250;
    loupe.style.left = lx + "px";
    loupe.style.top = ly + "px";
  });

  document.addEventListener("mouseleave", (e) => {
    const img = e.target;
    if (!img._loupe) return;
    img._loupe.remove();
    img._loupe = null;
  }, true);

  document.addEventListener("mouseout", (e) => {
    const img = e.target;
    if (img._loupe) {
      const rect = img.getBoundingClientRect();
      const cx = e.clientX - rect.left;
      const cy = e.clientY - rect.top;
      if (cx < 0 || cy < 0 || cx > rect.width || cy > rect.height) {
        img._loupe.remove();
        img._loupe = null;
      }
    }
  });
})();

let _allJobs = [];

async function loadPipelineJobs() {
  try {
    const res  = await fetch("/api/pipeline/jobs");
    const data = await res.json();
    _allJobs = data.jobs || [];
    renderPipelineJobs(_allJobs);
  } catch { /* ignore */ }
}


function startPipelinePolling() {
  if (pipelinePollInterval) return;
  loadNasDates();
  loadPipelineJobs();
  pipelinePollInterval = setInterval(() => { loadPipelineJobs(); }, 4000);
}

function stopPipelinePolling() {
  if (!pipelinePollInterval) return;
  clearInterval(pipelinePollInterval);
  pipelinePollInterval = null;
}

document.getElementById("pipeline-log-btn").addEventListener("click", (e) => {
  const id = e.currentTarget.dataset.jobId;
  if (id) openLogModal(id);
});

document.getElementById("pipeline-trigger-form").addEventListener("submit", async (e) => {
  e.preventDefault();
  const fits_dir       = document.getElementById("pipeline-fits-dir").value.trim();
  const target         = document.getElementById("pipeline-target").value.trim();
  const manualSelectCb = document.getElementById("pipeline-manual-select");
  const forceFreshCb   = document.getElementById("pipeline-force-fresh");
  const useColorCb     = document.getElementById("pipeline-use-color");
  const refineWcsCb    = document.getElementById("pipeline-refine-wcs");
  if (!fits_dir) return;
  const btn = document.getElementById("pipeline-trigger-btn");
  btn.disabled = true;
  btn.textContent = "Starting…";
  try {
    const targetSel = document.getElementById("pipeline-target-select");
    const selOpt = targetSel?.selectedOptions[0];
    const isSnapshot = selOpt?.dataset.snapshot === "1";
    const result = await postJson("/api/pipeline/trigger", {
      fits_dir,
      target: target || null,
      target_filter: isSnapshot ? (selOpt.dataset.name || target || null) : null,
      manual_selection: manualSelectCb?.checked ? true : undefined,
      force_fresh: forceFreshCb?.checked ? true : undefined,
      use_color: !!useColorCb?.checked,
      refine_wcs: !!refineWcsCb?.checked,
    });
    setLog(result);
    // Keep the date + target list in place so the next target can be queued
    // immediately without re-scanning — just reset the target selection fields.
    document.getElementById("pipeline-target").value = "";
    if (targetSel) targetSel.value = "";
    document.getElementById("pipeline-fits-dir").value = "";
    loadPipelineJobs();
  } catch (err) {
    setLog({ success: false, error: err.message });
  } finally {
    btn.disabled = false;
    btn.textContent = "▶ Trigger";
  }
});

// ── Frame selection modal ─────────────────────────────────────────────────────
(function initFrameSelectionModal() {
  const modal    = document.getElementById("frame-selection-modal");
  const grid     = document.getElementById("fsm-grid");
  const countEl  = document.getElementById("fsm-count");
  const btnAll   = document.getElementById("fsm-select-all");
  const btnNone  = document.getElementById("fsm-deselect-all");
  const btnConfirm = document.getElementById("fsm-confirm");
  const btnClose = document.getElementById("fsm-close");

  let _jobId  = null;
  let _frames = [];       // [{idx, file, name, preview, stars, fwhm, roundness, selected}]
  let _sortCol = "idx";
  let _sortDir = 1;       // 1 = asc, -1 = desc

  function updateCount() {
    const total   = _frames.length;
    const checked = _frames.filter(f => f._checked).length;
    countEl.textContent = `${checked} / ${total} selected`;
    btnConfirm.disabled = false;
    btnConfirm.textContent = checked === 0
      ? "⏭ Skip this filter"
      : "✓ Stack selected";
    btnConfirm.style.background = checked === 0 ? "#7f1d1d" : "#16a34a";
  }

  function qualityColor(val, goodThreshold, badThreshold, highIsGood) {
    if (val == null) return "";
    const isGood = highIsGood ? (val >= goodThreshold) : (val <= goodThreshold);
    const isBad  = highIsGood ? (val <= badThreshold)  : (val >= badThreshold);
    return isBad ? "fsm-bad" : isGood ? "fsm-good" : "";
  }

  function renderGrid() {
    const sorted = [..._frames].sort((a, b) => {
      let va = a[_sortCol], vb = b[_sortCol];
      if (va == null && vb == null) return 0;
      if (va == null) return 1;
      if (vb == null) return -1;
      return _sortDir * (va < vb ? -1 : va > vb ? 1 : 0);
    });
    grid.innerHTML = "";
    for (const f of sorted) {
      const fwhmCls  = qualityColor(f.fwhm, 3, 7, false);   // low FWHM = good
      const roundCls = qualityColor(f.roundness, 0.8, 0.6, true); // high roundness = good
      const card = document.createElement("div");
      card.className = "fsm-card" + (f._checked ? " selected" : " deselected");
      card.dataset.file = f.file;
      const previewUrl = f.preview ? `/data/${f.preview}` : "";
      card.innerHTML = `
        <label class="fsm-card-check">
          <input type="checkbox" class="fsm-cb" ${f._checked ? "checked" : ""} />
          <span title="${f.file}">#${f.idx} ${f.name.replace(/^.*[\\/]/, "")}</span>
        </label>
        <div class="fsm-thumb-wrap">
          ${previewUrl
            ? `<img class="fsm-thumb frame-zoom-img" src="${previewUrl}" alt="frame ${f.idx}" loading="lazy" />`
            : `<div style="display:flex;align-items:center;justify-content:center;height:100%;color:#4b5563;font-size:11px;">No preview</div>`
          }
        </div>
        <div class="fsm-stats">
          <div>Stars: <strong>${f.stars ?? "–"}</strong></div>
          <div>FWHM: <strong class="${fwhmCls}">${f.fwhm != null ? f.fwhm + " px" : "–"}</strong></div>
          <div>Roundness: <strong class="${roundCls}">${f.roundness != null ? f.roundness.toFixed(3) : "–"}</strong></div>
        </div>
      `;
      // Toggle selection
      card.querySelector(".fsm-cb").addEventListener("change", (ev) => {
        f._checked = ev.target.checked;
        card.className = "fsm-card" + (f._checked ? " selected" : " deselected");
        updateCount();
      });
      grid.appendChild(card);
    }
    updateCount();
  }

  // Attach zoom to newly added images
  function attachZoom() {
    for (const img of grid.querySelectorAll(".frame-zoom-img")) {
      if (img._zoomAttached) continue;
      img._zoomAttached = true;
      img.addEventListener("mousemove", (e) => {
        if (!img._loupe) {
          const loupe = document.createElement("div");
          Object.assign(loupe.style, {
            position: "fixed", width: "240px", height: "240px",
            border: "2px solid #4fc3f7", borderRadius: "4px",
            background: `url("${img.src}") no-repeat`,
            pointerEvents: "none", zIndex: "9999",
            boxShadow: "0 4px 20px rgba(0,0,0,0.6)",
          });
          document.body.appendChild(loupe);
          img._loupe = loupe;
        }
        const rect = img.getBoundingClientRect();
        const cx = e.clientX - rect.left;
        const cy = e.clientY - rect.top;
        const Z = 5;
        img._loupe.style.backgroundSize = `${rect.width * Z}px ${rect.height * Z}px`;
        img._loupe.style.backgroundPosition = `${-(cx * Z - 120)}px ${-(cy * Z - 120)}px`;
        img._loupe.style.display = "block";
        let lx = e.clientX + 20;
        let ly = e.clientY - 120;
        if (lx + 250 > window.innerWidth) lx = e.clientX - 260;
        if (ly < 10) ly = 10;
        img._loupe.style.left = lx + "px";
        img._loupe.style.top = ly + "px";
      });
      img.addEventListener("mouseleave", () => {
        if (img._loupe) { img._loupe.remove(); img._loupe = null; }
      });
    }
  }

  btnAll.addEventListener("click", () => {
    _frames.forEach(f => { f._checked = true; });
    renderGrid(); attachZoom();
  });
  btnNone.addEventListener("click", () => {
    _frames.forEach(f => { f._checked = false; });
    renderGrid(); attachZoom();
  });

  // Sort buttons
  document.querySelectorAll(".fsm-sort-btn").forEach(btn => {
    btn.addEventListener("click", () => {
      const col = btn.dataset.col;
      if (_sortCol === col) _sortDir *= -1;
      else { _sortCol = col; _sortDir = 1; }
      document.querySelectorAll(".fsm-sort-btn").forEach(b => b.style.fontWeight = "");
      btn.style.fontWeight = "700";
      renderGrid(); attachZoom();
    });
  });

  btnClose.addEventListener("click", () => { modal.style.display = "none"; });
  modal.addEventListener("click", (e) => { if (e.target === modal) modal.style.display = "none"; });

  btnConfirm.addEventListener("click", async () => {
    const selectedFiles = _frames.filter(f => f._checked).map(f => f.file);
    if (!selectedFiles.length) return;
    btnConfirm.disabled = true;
    btnConfirm.textContent = "Sending…";
    try {
      const r = await postJson(`/api/pipeline/jobs/${_jobId}/selection`, { files: selectedFiles });
      if (r.success) {
        modal.style.display = "none";
        loadPipelineJobs();
      } else {
        btnConfirm.textContent = `✗ ${r.error}`;
        btnConfirm.disabled = false;
      }
    } catch (err) {
      btnConfirm.textContent = `✗ ${err.message}`;
      btnConfirm.disabled = false;
    }
  });

  window.openFrameSelectionModal = async function(jobId) {
    _jobId = jobId;
    grid.innerHTML = '<p style="color:#9ca3af;padding:20px;">Loading frame data…</p>';
    modal.style.display = "block";
    countEl.textContent = "";
    btnConfirm.disabled = false;
    btnConfirm.textContent = "✓ Stack selected";
    btnConfirm.style.background = "#16a34a";
    try {
      const data = await (await fetch(`/api/pipeline/jobs/${jobId}/selection`)).json();
      if (!data.success) {
        grid.innerHTML = `<p style="color:#f87171;padding:20px;">Error: ${data.error}</p>`;
        return;
      }
      _frames = (data.frames || []).map(f => ({ ...f, _checked: f.selected !== false }));
      _sortCol = "idx"; _sortDir = 1;
      renderGrid();
      attachZoom();
    } catch (err) {
      grid.innerHTML = `<p style="color:#f87171;padding:20px;">Error: ${err.message}</p>`;
    }
  };
})();

setActiveTab("devices");
loadDefaults().then(() => {
  refreshStatus();
  startNinaStatusPolling();
});

let _ninaStatusPollTimer = null;
function startNinaStatusPolling() {
  if (_ninaStatusPollTimer) return;
  _ninaStatusPollTimer = setInterval(() => {
    if (document.visibilityState === "hidden") return;
    refreshStatus().catch(() => {});
  }, 12_000);
}
loadWatchdogConfig();
loadArbiterConfig();
loadHistory();

// ── Manual mode live toggle ───────────────────────────────────────────────────
const manualModeChk = document.getElementById("seq-manual-mode");
if (manualModeChk) {
  manualModeChk.addEventListener("change", async () => {
    try { await postJson("/api/sequence/manual-mode", { enabled: manualModeChk.checked }); }
    catch { /* ignore */ }
  });
}

// ── Measurement History ───────────────────────────────────────────────────────
let _histRows = [];   // last fetched rows for CSV copy

async function loadHistory() {
  const target   = document.getElementById("hist-target").value.trim();
  const filter   = document.getElementById("hist-filter").value;
  const dateFrom = document.getElementById("hist-date-from").value;
  const dateTo   = document.getElementById("hist-date-to").value;

  const params = new URLSearchParams();
  if (target)   params.set("target",    target);
  if (filter)   params.set("filter",    filter);
  if (dateFrom) params.set("date_from", dateFrom);
  if (dateTo)   params.set("date_to",   dateTo);

  const emptyEl = document.getElementById("hist-empty");
  const tableEl = document.getElementById("hist-table");
  const tbody   = document.getElementById("hist-tbody");
  emptyEl.textContent = "Loading…";
  emptyEl.style.display = "";
  tableEl.style.display = "none";

  try {
    const res  = await fetch(`/api/pipeline/history?${params}`);
    const data = await res.json();
    _histRows = data.results || [];

    if (!_histRows.length) {
      emptyEl.textContent = "No results found.";
      return;
    }

    tbody.innerHTML = _histRows.map(r => {
      const mag_ap  = r.mag_ap    != null ? `${r.mag_ap.toFixed(2)} ± ${(r.magerr_ap||0).toFixed(2)}` : "–";
      const mag_sub = r.mag_sub   != null ? `${r.mag_sub.toFixed(2)} ± ${(r.magerr_sub||0).toFixed(2)}`
                    : r.mag_sub_ul != null ? `&gt;${r.mag_sub_ul.toFixed(2)}` : "–";
      const mjd     = r.mjd != null ? r.mjd.toFixed(5) : "–";
      const link    = r.stdweb_url
        ? `<a href="${r.stdweb_url}" target="_blank" style="color:#60a5fa;font-size:11px;">#${r.stdweb_task_id} ↗</a>`
        : "–";
      const frames  = r.n_frames ? `${r.n_frames}×${r.exposure}s` : "–";
      const tgt     = (r.target || "").replace(/^(sn_|at_)/, m => m.slice(0,-1).toUpperCase() + " ");
      const hasMag  = r.mag_ap != null || r.mag_sub != null;
      return `<tr style="border-bottom:1px solid #0f2035;${hasMag ? '' : 'opacity:0.5;'}">
        <td style="padding:5px 8px;white-space:nowrap;">${r.obs_date || "–"}</td>
        <td style="padding:5px 8px;color:#a78bfa;">${tgt || r.target || "–"}</td>
        <td style="padding:5px 8px;text-align:center;">${r.filter || "–"}</td>
        <td style="padding:5px 8px;text-align:center;color:#6b7280;">${frames}</td>
        <td style="padding:5px 8px;text-align:right;font-family:monospace;font-size:11px;">${mjd}</td>
        <td style="padding:5px 8px;text-align:right;font-family:monospace;color:#6ee7b7;">${mag_ap}</td>
        <td style="padding:5px 8px;text-align:right;font-family:monospace;color:#93c5fd;">${mag_sub}</td>
        <td style="padding:5px 8px;text-align:center;">${link}</td>
      </tr>`;
    }).join("");

    emptyEl.style.display = "none";
    tableEl.style.display = "";
  } catch (e) {
    emptyEl.textContent = `Error: ${e.message}`;
  }
}

document.getElementById("hist-search-btn").addEventListener("click", loadHistory);
// Search on Enter in the target field
document.getElementById("hist-target").addEventListener("keydown", e => { if (e.key === "Enter") loadHistory(); });

document.getElementById("hist-copy-all-btn").addEventListener("click", async () => {
  if (!_histRows.length) return;
  const lines = _histRows
    .filter(r => r.mag_ap != null || r.mag_sub != null)
    .flatMap(r => {
      const rows = [];
      const tgt   = (r.target || "").replace(/_/g, " ").toUpperCase();
      const stack = (r.n_frames && r.exposure) ? `, ${r.n_frames}×${r.exposure}s` : "";
      if (r.mag_ap  != null) rows.push(`${tgt}, ${r.mjd?.toFixed(8) ?? ""}, ${r.filter}${stack}, Aperture photometry, ${r.mag_ap.toFixed(2)}, ${(r.magerr_ap||0).toFixed(2)}`);
      if (r.mag_sub != null) rows.push(`${tgt}, ${r.mjd?.toFixed(8) ?? ""}, ${r.filter}${stack}, Template substraction aperture photometry, ${r.mag_sub.toFixed(2)}, ${(r.magerr_sub||0).toFixed(2)}`);
      return rows;
    });
  const text = lines.join("\n");
  const btn = document.getElementById("hist-copy-all-btn");
  const done = () => {
    btn.textContent = "✓ Copied";
    setTimeout(() => { btn.textContent = "⧉ Copy CSV"; }, 2000);
  };
  const fail = () => {
    btn.textContent = "✗ Failed";
    setTimeout(() => { btn.textContent = "⧉ Copy CSV"; }, 2000);
  };
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(done).catch(() => fallbackCopy(text, done));
  } else {
    fallbackCopy(text, done);
  }
});

// ── OCS History Timeline Chart ────────────────────────────────────────────────

let _ocsChart = null;

/** Convert SQLite UTC string "2026-03-28 08:53:20" → epoch ms */
function ocsTs(isoStr) {
  return new Date(isoStr.replace(" ", "T") + "Z").getTime();
}

function formatTsMs(ms) {
  const d = new Date(ms);
  const mm  = String(d.getUTCMonth() + 1).padStart(2, "0");
  const dd  = String(d.getUTCDate()).padStart(2, "0");
  const hh  = String(d.getUTCHours()).padStart(2, "0");
  const min = String(d.getUTCMinutes()).padStart(2, "0");
  return `${mm}/${dd} ${hh}:${min}`;
}

function ocsRoofValue(roofStr) {
  if (!roofStr) return null;
  if (/open/i.test(roofStr)) return 1;
  if (/close/i.test(roofStr)) return 0;
  return null;
}

function ocsRainValue(rainStr) {
  if (!rainStr) return null;
  return /rain|yes|wet|true/i.test(rainStr) ? 1 : 0;
}

// Original X bounds for zoom reset
let _ocsXMin = null, _ocsXMax = null;
let _ocsDragStart = null; // canvas-relative X px where drag began

function _ocsDragZoomStart(e) {
  if (!_ocsChart) return;
  const area = _ocsChart.chartArea;
  const rect = _ocsChart.canvas.getBoundingClientRect();
  const x = e.clientX - rect.left;
  if (x < area.left || x > area.right) return; // only inside chart area
  _ocsDragStart = x;
  e.preventDefault();
}

function _ocsDragZoomMove(e) {
  if (!_ocsChart || _ocsDragStart === null) return;
  const canvas  = _ocsChart.canvas;
  const area    = _ocsChart.chartArea;
  const rect    = canvas.getBoundingClientRect();
  const xNow    = Math.max(area.left, Math.min(area.right, e.clientX - rect.left));

  // Draw selection overlay via the chart's selection plugin (redraws each move)
  _ocsChart._dragCurrent = xNow;
  _ocsChart.update("none");
}

function _ocsDragZoomEnd(e) {
  if (!_ocsChart || _ocsDragStart === null) return;
  const canvas  = _ocsChart.canvas;
  const area    = _ocsChart.chartArea;
  const rect    = canvas.getBoundingClientRect();
  const xNow    = Math.max(area.left, Math.min(area.right, e.clientX - rect.left));
  const x0      = Math.min(_ocsDragStart, xNow);
  const x1      = Math.max(_ocsDragStart, xNow);
  const start   = _ocsDragStart;
  _ocsDragStart = null;
  _ocsChart._dragCurrent = null;

  if (x1 - x0 < 8) { // too small = ignore (just a click)
    _ocsChart.update("none");
    return;
  }

  const xScale = _ocsChart.scales.x;
  const newMin = xScale.getValueForPixel(x0);
  const newMax = xScale.getValueForPixel(x1);
  if (newMax - newMin < 60000) return; // ignore sub-minute selections

  _ocsChart.options.scales.x.min = newMin;
  _ocsChart.options.scales.x.max = newMax;
  _ocsChart.update("none");
}

function renderOcsChart(history) {
  const canvas = document.getElementById("ocs-history-chart");
  if (!canvas) return;

  // Always remove listeners and destroy any chart instance on this canvas
  canvas.removeEventListener("mousedown",  _ocsDragZoomStart);
  canvas.removeEventListener("mousemove",  _ocsDragZoomMove);
  canvas.removeEventListener("mouseup",    _ocsDragZoomEnd);
  canvas.removeEventListener("mouseleave", _ocsDragZoomEnd);
  const existing = Chart.getChart(canvas);
  if (existing) existing.destroy();
  if (_ocsChart) { _ocsChart.destroy(); }
  _ocsChart = null;
  _ocsDragStart = null;

  if (!history || history.length === 0) {
    canvas.style.height = "80px";
    const ctx = canvas.getContext("2d");
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = "#555";
    ctx.font = "14px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("No OCS data yet — readings stored every 10 min", canvas.width / 2, 40);
    return;
  }
  canvas.style.height = "320px";

  // Convert each reading's timestamp to epoch ms for accurate X positioning
  const tsMs      = history.map(r => ocsTs(r.ts));
  const roofOpen  = history.map(r => ocsRoofValue(r.roof));
  const rainAct   = history.map(r => ocsRainValue(r.rain));
  const xMin      = tsMs[0];
  const xMax      = tsMs[tsMs.length - 1];
  _ocsXMin = xMin;
  _ocsXMax = xMax;
  const rangeMs   = Math.max(xMax - xMin, 60000); // at least 1 min to avoid division by zero

  // Half-interval for background band edges (use actual gap or 5 min fallback)
  function halfGap(i) {
    if (tsMs.length < 2) return 5 * 60000;
    const prev = i > 0 ? (tsMs[i] - tsMs[i - 1]) / 2 : (tsMs[1] - tsMs[0]) / 2;
    const next = i < tsMs.length - 1 ? (tsMs[i + 1] - tsMs[i]) / 2 : (tsMs[tsMs.length - 1] - tsMs[tsMs.length - 2]) / 2;
    return { left: prev, right: next };
  }

  const backgroundPlugin = {
    id: "ocsBands",
    beforeDraw(chart) {
      const { ctx, chartArea: { left: cLeft, right: cRight, top, bottom }, scales: { x } } = chart;
      const stripH   = (bottom - top) * 0.20;  // rain strip = bottom 20%
      const stripTop = bottom - stripH;
      const n = history.length;

      for (let i = 0; i < n; i++) {
        const g  = halfGap(i);
        const x1 = Math.max(cLeft,  x.getPixelForValue(tsMs[i] - g.left));
        const x2 = Math.min(cRight, x.getPixelForValue(tsMs[i] + g.right));
        const w  = x2 - x1;
        if (w <= 0) continue;

        // Roof open → subtle green tint full height
        if (roofOpen[i] === 1) {
          ctx.fillStyle = "rgba(34,197,94,0.10)";
          ctx.fillRect(x1, top, w, bottom - top);
        }
        // Rain → solid indigo bar in bottom 20%
        if (rainAct[i] === 1) {
          ctx.fillStyle = "rgba(99,102,241,0.75)";
          ctx.fillRect(x1, stripTop, w, stripH);
        }
      }

      // Rain strip separator line
      ctx.strokeStyle = "rgba(99,102,241,0.25)";
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(cLeft, stripTop);
      ctx.lineTo(cRight, stripTop);
      ctx.stroke();
    },
  };

  // X tick step: aim for ~8 labels across the range
  const nTicks   = 8;
  const stepMs   = Math.pow(10, Math.ceil(Math.log10(rangeMs / nTicks)));
  const niceStep = [1, 2, 3, 6, 12, 24].map(h => h * 3600000)
    .find(s => rangeMs / s <= nTicks) || stepMs;

  const toXY = (r) => ({ x: ocsTs(r.ts), y: null });  // placeholder

  const tempData     = history.map(r => ({ x: ocsTs(r.ts), y: r.temp     ?? null }));
  const humidityData = history.map(r => ({ x: ocsTs(r.ts), y: r.humidity ?? null }));
  const pressureData = history.map(r => ({ x: ocsTs(r.ts), y: r.pressure ?? null }));
  const skyData      = history.map(r => ({ x: ocsTs(r.ts), y: r.sky      ?? null }));
  const irSkyData    = history.map(r => ({ x: ocsTs(r.ts), y: r.ir_sky   ?? null }));

  const hasSky      = history.some(r => r.sky    != null);
  const hasPressure = history.some(r => r.pressure != null);
  const hasIrSky    = history.some(r => r.ir_sky  != null);

  const datasets = [
    {
      label: "Temp (\u00b0C)",
      data: tempData,
      borderColor: "#f59e0b",
      backgroundColor: "rgba(245,158,11,0.08)",
      pointRadius: 2,
      borderWidth: 2,
      tension: 0.3,
      yAxisID: "yTemp",
      spanGaps: true,
    },
    {
      label: "Humidity (%)",
      data: humidityData,
      borderColor: "#38bdf8",
      backgroundColor: "rgba(56,189,248,0.08)",
      pointRadius: 2,
      borderWidth: 2,
      tension: 0.3,
      yAxisID: "yHumid",
      spanGaps: true,
    },
  ];

  if (hasPressure) {
    datasets.push({
      label: "Pressure (hPa)",
      data: pressureData,
      borderColor: "#a78bfa",
      backgroundColor: "rgba(167,139,250,0.08)",
      pointRadius: 2,
      borderWidth: 2,
      tension: 0.3,
      yAxisID: "yPres",
      spanGaps: true,
    });
  }

  if (hasSky) {
    datasets.push({
      label: "Sky Quality (mag/\u2033\u00b2)",
      data: skyData,
      borderColor: "#22d3ee",
      backgroundColor: "rgba(34,211,238,0.08)",
      pointRadius: 2,
      borderWidth: 2,
      tension: 0.3,
      yAxisID: "ySky",
      spanGaps: true,
    });
  }

  if (hasIrSky) {
    datasets.push({
      label: "IR Sky (\u00b0C)",
      data: irSkyData,
      borderColor: "#f87171",
      backgroundColor: "rgba(248,113,113,0.08)",
      pointRadius: 2,
      borderWidth: 2,
      borderDash: [5, 3],
      tension: 0.3,
      yAxisID: "yTemp",  // shares °C axis with ambient temp
      spanGaps: true,
    });
  }

  const scales = {
    x: {
      type: "linear",
      min: xMin,
      max: xMax,
      ticks: {
        stepSize: niceStep,
        callback: (v) => formatTsMs(v),
        color: "#888",
        maxRotation: 35,
        maxTicksLimit: 10,
      },
      grid: { color: "rgba(255,255,255,0.05)" },
    },
    yTemp: {
      type: "linear",
      position: "left",
      title: { display: true, text: "\u00b0C / SQ", color: "#f59e0b" },
      ticks: { color: "#888" },
      grid: { color: "rgba(255,255,255,0.07)" },
    },
    yHumid: {
      type: "linear",
      position: "right",
      min: 0,
      max: 100,
      title: { display: true, text: "%", color: "#38bdf8" },
      ticks: { color: "#888" },
      grid: { display: false },
    },
  };

  if (hasPressure) {
    const pVals = history.map(r => r.pressure).filter(v => v != null);
    const pMin  = Math.floor(Math.min(...pVals) / 5) * 5 - 5;
    const pMax  = Math.ceil(Math.max(...pVals) / 5) * 5 + 5;
    scales.yPres = {
      type: "linear",
      position: "right",
      min: pMin,
      max: pMax,
      title: { display: true, text: "hPa", color: "#a78bfa" },
      ticks: { color: "#888" },
      grid: { display: false },
    };
  }

  if (hasSky) {
    const skyVals = history.map(r => r.sky).filter(v => v != null);
    const sMin = Math.floor(Math.min(...skyVals)) - 1;
    const sMax = Math.ceil(Math.max(...skyVals)) + 1;
    scales.ySky = {
      type: "linear",
      position: "left",
      min: sMin,
      max: sMax,
      title: { display: false },
      ticks: { display: false },
      grid: { display: false },
    };
  }

  const selectionPlugin = {
    id: "ocsDragSelect",
    afterDraw(chart) {
      if (_ocsDragStart === null || chart._dragCurrent === null) return;
      const { ctx, chartArea: { top, bottom } } = chart;
      const x0 = Math.min(_ocsDragStart, chart._dragCurrent);
      const x1 = Math.max(_ocsDragStart, chart._dragCurrent);
      ctx.save();
      ctx.fillStyle = "rgba(99,179,237,0.15)";
      ctx.fillRect(x0, top, x1 - x0, bottom - top);
      ctx.strokeStyle = "rgba(99,179,237,0.7)";
      ctx.lineWidth = 1;
      ctx.setLineDash([4, 3]);
      ctx.strokeRect(x0, top, x1 - x0, bottom - top);
      ctx.restore();
    },
  };

  _ocsChart = new Chart(canvas.getContext("2d"), {
    type: "line",
    data: { datasets },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      animation: false,
      interaction: { mode: "index", intersect: false },
      plugins: {
        legend: {
          labels: { color: "#aaa", boxWidth: 20, padding: 12, font: { size: 11 } },
        },
        tooltip: {
          backgroundColor: "#1a1a2e",
          titleColor: "#ccc",
          bodyColor: "#aaa",
          borderColor: "#444",
          borderWidth: 1,
          callbacks: {
            title: (items) => {
              const tsVal = items[0].parsed.x;
              // Find closest history row by timestamp
              let closest = history[0];
              let minDiff = Infinity;
              for (const r of history) {
                const d = Math.abs(ocsTs(r.ts) - tsVal);
                if (d < minDiff) { minDiff = d; closest = r; }
              }
              const parts = [formatTsMs(ocsTs(closest.ts))];
              if (closest.roof) parts.push(`Roof: ${closest.roof}`);
              if (closest.rain) parts.push(`Rain: ${closest.rain}`);
              if (closest.safe) parts.push(`Safe: ${closest.safe}`);
              return parts;
            },
          },
        },
      },
      scales,
    },
    plugins: [backgroundPlugin, selectionPlugin],
  });

  canvas.addEventListener("mousedown",  _ocsDragZoomStart);
  canvas.addEventListener("mousemove",  _ocsDragZoomMove);
  canvas.addEventListener("mouseup",    _ocsDragZoomEnd);
  canvas.addEventListener("mouseleave", _ocsDragZoomEnd);
  canvas.style.cursor = "crosshair";
}

async function loadOcsHistory() {
  const statusEl = document.getElementById("ocs-history-status");
  const hoursEl  = document.getElementById("ocs-history-hours");
  const hours = hoursEl ? hoursEl.value : 48;
  if (statusEl) statusEl.textContent = "Loading\u2026";
  try {
    const resp = await fetch(`/api/ocs/history?hours=${hours}`);
    const data = await resp.json();
    renderOcsChart(data.history || []);
    if (statusEl) {
      const n = data.history?.length ?? 0;
      statusEl.textContent = n
        ? `${n} readings over last ${hours}h`
        : "No data yet \u2014 check back after first polling cycle (~10 min after server start)";
    }
  } catch (e) {
    if (statusEl) statusEl.textContent = `Error: ${e.message}`;
  }
}

document.getElementById("ocs-history-refresh")?.addEventListener("click", loadOcsHistory);
document.getElementById("ocs-history-hours")?.addEventListener("change", loadOcsHistory);

// ── Rain radar forecast ───────────────────────────────────────────────────────
async function loadRainForecast() {
  const bar  = document.getElementById("rain-forecast-bar");
  const warn = document.getElementById("rain-forecast-warn");
  const upd  = document.getElementById("rain-forecast-updated");
  if (!bar) return;
  try {
    const r = await fetch("/api/weather/forecast");
    const d = await r.json();
    if (!d.success) { bar.innerHTML = `<span style="font-size:12px;color:#f87171;">${d.error}</span>`; return; }

    // Render mini bar chart
    const maxH = 40;
    bar.innerHTML = d.slots.map(s => {
      const p = s.prob ?? 0;
      const h = Math.max(4, Math.round(p / 100 * maxH));
      const col = p >= 70 ? "#ef4444" : p >= 40 ? "#f59e0b" : "#22c55e";
      const time = s.time.slice(11, 16);
      return `<div style="display:flex;flex-direction:column;align-items:center;gap:2px;min-width:28px;">
        <span style="font-size:10px;color:#94a3b8;">${p}%</span>
        <div style="width:22px;height:${h}px;background:${col};border-radius:2px;"></div>
        <span style="font-size:10px;color:#64748b;">${time}</span>
      </div>`;
    }).join("");

    // Warning banner
    const warnTexts = { imminent: "Rain imminent within 30 min", likely: "Rain likely within 1 h", possible: "Rain possible within 1 h", clear: null };
    const warnStyles = { imminent: "background:#7f1d1d;color:#fca5a5;border:1px solid #991b1b;", likely: "background:#78350f;color:#fcd34d;border:1px solid #92400e;", possible: "background:#1e3a5f;color:#93c5fd;border:1px solid #1e40af;", clear: "" };
    if (warn) {
      const txt = warnTexts[d.warnLevel];
      if (txt) { warn.style.display = ""; warn.style.cssText += warnStyles[d.warnLevel]; warn.textContent = txt + ` (${d.maxProb30}% in 30 min, ${d.maxProb60}% in 1 h)`; }
      else { warn.style.display = "none"; }
    }
    if (upd) upd.textContent = `Updated ${new Date(d.fetchedAt).toLocaleTimeString()}`;
  } catch (e) {
    if (bar) bar.innerHTML = `<span style="font-size:12px;color:#f87171;">Forecast unavailable</span>`;
  }
}

document.getElementById("ocs-backfill")?.addEventListener("click", async () => {
  const btn = document.getElementById("ocs-backfill");
  const statusEl = document.getElementById("ocs-history-status");
  btn.disabled = true;
  btn.textContent = "Backfilling…";
  if (statusEl) statusEl.textContent = "Scraping 60-min history from OCS…";
  try {
    const resp = await fetch("/api/ocs/backfill", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ ocsHost: currentOcsHost() }),
    });
    const data = await resp.json();
    if (data.success) {
      if (statusEl) statusEl.textContent = `Backfill done — ${data.inserted} new rows inserted (${data.total} points from OCS). Reloading…`;
      await loadOcsHistory();
    } else {
      if (statusEl) statusEl.textContent = `Backfill failed: ${data.error}`;
    }
  } catch (e) {
    if (statusEl) statusEl.textContent = `Error: ${e.message}`;
  } finally {
    btn.disabled = false;
    btn.textContent = "Backfill 60 min";
  }
});

document.getElementById("ocs-poll-now")?.addEventListener("click", async () => {
  const btn = document.getElementById("ocs-poll-now");
  const statusEl = document.getElementById("ocs-history-status");
  btn.disabled = true;
  btn.textContent = "Polling…";
  if (statusEl) statusEl.textContent = "Fetching OCS reading…";
  try {
    const resp = await fetch("/api/ocs/poll-now", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ ocsHost: currentOcsHost() }),
    });
    const data = await resp.json();
    if (data.success) {
      if (statusEl) statusEl.textContent = "Reading stored — reloading chart…";
      await loadOcsHistory();
    } else {
      if (statusEl) statusEl.textContent = `Poll failed: ${data.error}`;
    }
  } catch (e) {
    if (statusEl) statusEl.textContent = `Error: ${e.message}`;
  } finally {
    btn.disabled = false;
    btn.textContent = "Poll Now";
  }
});

// ── IR Thermal Camera heatmap ─────────────────────────────────────────────────

function renderIrHeatmap(pixels, avg, cols, rows) {
  const canvas = document.getElementById("ocs-ir-canvas");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  const W = canvas.width, H = canvas.height;
  const cw = W / cols, ch = H / rows;

  const temps = pixels.map(p => p.temp);
  const tMin = Math.min(...temps);
  const tMax = Math.max(...temps);
  const range = tMax - tMin || 1;

  // Inferno-like colormap: black → purple → red → orange → yellow
  function tempColor(t) {
    const frac = Math.max(0, Math.min(1, (t - tMin) / range));
    // 5-stop gradient
    const stops = [
      [0.0,  [10,  5,  30]],
      [0.25, [90,  0, 120]],
      [0.5,  [200, 30,  20]],
      [0.75, [240,120,   0]],
      [1.0,  [252,255, 160]],
    ];
    for (let i = 0; i < stops.length - 1; i++) {
      const [t0, c0] = stops[i];
      const [t1, c1] = stops[i + 1];
      if (frac >= t0 && frac <= t1) {
        const u = (frac - t0) / (t1 - t0);
        const r = Math.round(c0[0] + u * (c1[0] - c0[0]));
        const g = Math.round(c0[1] + u * (c1[1] - c0[1]));
        const b = Math.round(c0[2] + u * (c1[2] - c0[2]));
        return `rgb(${r},${g},${b})`;
      }
    }
    return "white";
  }

  ctx.clearRect(0, 0, W, H);
  for (const { col, row, temp } of pixels) {
    ctx.fillStyle = tempColor(temp);
    ctx.fillRect(col * cw, row * ch, cw, ch);
    // Temp label on each cell
    ctx.fillStyle = (temp - tMin) / range > 0.6 ? "#000" : "#fff";
    ctx.font = `${Math.min(cw * 0.38, 12)}px monospace`;
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText(temp.toFixed(1), col * cw + cw / 2, row * ch + ch / 2);
  }

  // Grid lines
  ctx.strokeStyle = "rgba(0,0,0,0.3)";
  ctx.lineWidth = 0.5;
  for (let c = 1; c < cols; c++) { ctx.beginPath(); ctx.moveTo(c * cw, 0); ctx.lineTo(c * cw, H); ctx.stroke(); }
  for (let r = 1; r < rows; r++) { ctx.beginPath(); ctx.moveTo(0, r * ch); ctx.lineTo(W, r * ch); ctx.stroke(); }

  const minEl = document.getElementById("ocs-ir-min");
  const maxEl = document.getElementById("ocs-ir-max");
  const avgEl = document.getElementById("ocs-ir-avg");
  if (minEl) minEl.textContent = `Min: ${tMin.toFixed(1)}°C`;
  if (maxEl) maxEl.textContent = `Max: ${tMax.toFixed(1)}°C`;
  if (avgEl) avgEl.textContent = avg != null ? `Avg: ${avg.toFixed(1)}°C` : "";
}

async function loadIrImage() {
  const btn     = document.getElementById("ocs-ir-refresh");
  const statusEl = document.getElementById("ocs-ir-status");
  if (btn) { btn.disabled = true; btn.textContent = "Loading…"; }
  if (statusEl) statusEl.textContent = "";
  try {
    const resp = await fetch(`/api/ocs/ir-image?ocsHost=${encodeURIComponent(currentOcsHost())}`);
    const data = await resp.json();
    if (data.success && data.pixels.length) {
      renderIrHeatmap(data.pixels, data.avg, data.cols, data.rows);
      if (statusEl) statusEl.textContent = `avg ${data.avg?.toFixed(1)}°C · ${data.pixels.length} pixels`;
    } else {
      if (statusEl) statusEl.textContent = data.error || "No IR data";
    }
  } catch (e) {
    if (statusEl) statusEl.textContent = `Error: ${e.message}`;
  } finally {
    if (btn) { btn.disabled = false; btn.textContent = "Refresh"; }
  }
}

document.getElementById("ocs-ir-refresh")?.addEventListener("click", loadIrImage);

// ── Alerts tab ───────────────────────────────────────────────────────────────

let _alertsCache = [];
let _alertsAutoTimer = null;

const BROKER_COLORS = {
  Swift: "#f59e0b", Fermi: "#f97316", IceCube: "#38bdf8", LVK: "#a78bfa",
  EP: "#34d399", SVOM: "#2dd4bf", INTEGRAL: "#e879f9", AGILE: "#fb923c",
  GECAM: "#facc15", MAXI: "#f472b6", IPN: "#94a3b8", "SN-nu": "#ef4444",
  GCN: "#6b7280",
};

const ACTION_LABELS = {
  queued:           { text: "QUEUED",   cls: "alert-action-queued"      },
  rejected:         { text: "REJECTED", cls: "alert-action-rejected"    },
  "no-position":    { text: "NO POS",   cls: "alert-action-no-position" },
  manual:           { text: "MANUAL",   cls: "alert-action-manual"      },
  ignored:          { text: "IGNORED",  cls: "alert-action-ignored"     },
  "too-observed":   { text: "ToO OBS",  cls: "alert-action-too"         },
};

function fmtDeg(v, prec = 3) { return v != null ? Number(v).toFixed(prec) : "—"; }
function fmtDeg1(v) { return v != null ? Number(v).toFixed(1) + "°" : "—"; }

async function loadAlertConfig() {
  try {
    const r = await fetch("/api/alerts/config");
    const d = await r.json();
    if (!d.success) return;
    const cb = document.getElementById("alert-ready-cb");
    const st = document.getElementById("alert-ready-status");
    if (cb) cb.checked = d.alertReady;
    if (st) {
      st.textContent = d.alertReady ? "ARMED" : "OFF";
      st.className = "alert-ready-status " + (d.alertReady ? "armed" : "");
    }
    // Global email notify checkbox
    const notifyCb = document.getElementById("alert-notify-email-cb");
    if (notifyCb) {
      const res = await fetch("/api/settings/alert_notify_email");
      if (res.ok) { const sd = await res.json(); notifyCb.checked = sd.value !== "false"; }
      else notifyCb.checked = true;
    }

    // Prep state
    const prepBtn    = document.getElementById("alert-prep-btn");
    const abortBtn   = document.getElementById("alert-prep-abort-btn");
    const prepStatus = document.getElementById("alert-prep-status");
    if (prepBtn && prepStatus) {
      if (d.alertPrepRunning) {
        prepBtn.disabled = true;
        prepBtn.classList.add("running");
        prepBtn.textContent = "Prep Running...";
        abortBtn.style.display = "";
        prepStatus.textContent = d.alertPrepStep || "Running...";
        prepStatus.className = "alert-prep-status";
        if (!_alertPrepPollTimer) {
          _alertPrepPollTimer = setInterval(loadAlertConfig, 3000);
        }
      } else {
        prepBtn.disabled = false;
        prepBtn.classList.remove("running");
        prepBtn.textContent = "Ready Obs for Alert";
        abortBtn.style.display = "none";
        if (_alertPrepPollTimer) { clearInterval(_alertPrepPollTimer); _alertPrepPollTimer = null; }
        if (d.alertPrepStep) {
          const isDone = d.alertPrepStep.startsWith("READY");
          const isErr  = d.alertPrepStep.startsWith("Error") || d.alertPrepStep.startsWith("ABORTED");
          prepStatus.textContent = d.alertPrepStep;
          prepStatus.className = "alert-prep-status" + (isDone ? " done" : isErr ? " error" : "");
        }
      }
    }
  } catch {}
}

document.getElementById("alert-ready-cb")?.addEventListener("change", async (e) => {
  const val = e.target.checked;
  try {
    await fetch("/api/alerts/config", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ alertReady: val }),
    });
    loadAlertConfig();
  } catch {}
});

document.getElementById("alert-notify-email-cb")?.addEventListener("change", async (e) => {
  await fetch("/api/settings", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ key: "alert_notify_email", value: e.target.checked ? "true" : "false" }),
  });
});

document.getElementById("alert-prep-btn")?.addEventListener("click", async () => {
  const hostEl = document.getElementById("nina-host");
  const portEl = document.getElementById("nina-port");
  const protocolEl = document.getElementById("nina-protocol");
  const host = hostEl?.value || "192.168.1.174";
  const port = portEl?.value || "1888";
  const protocol = protocolEl?.value || "http";
  try {
    const r = await fetch("/api/alerts/prep", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ host, port, protocol }),
    });
    const d = await r.json();
    if (!d.success) alert(d.error || "Failed to start alert prep");
    // Start polling the config to show progress
    if (!_alertPrepPollTimer) {
      _alertPrepPollTimer = setInterval(loadAlertConfig, 3000);
    }
    loadAlertConfig();
  } catch (e) { alert("Error: " + e.message); }
});

let _alertPrepPollTimer = null;

document.getElementById("alert-prep-abort-btn")?.addEventListener("click", async () => {
  try {
    await fetch("/api/alerts/prep/abort", { method: "POST" });
    loadAlertConfig();
  } catch {}
});

// ── Per-broker strategy table ─────────────────────────────────────────────────
let _strategiesCache = [];

async function loadStrategies() {
  try {
    const r = await fetch("/api/alert-strategies");
    const d = await r.json();
    if (d.success) { _strategiesCache = d.strategies || []; renderStrategies(); }
  } catch (e) { console.warn("[strategies-ui]", e); }
}

function renderStrategies() {
  const tbody = document.getElementById("strategy-tbody");
  if (!tbody) return;

  tbody.innerHTML = _strategiesCache.map(s => {
    const bc = BROKER_COLORS[s.broker] || "#6b7280";
    const dimCls = (!s.enabled || s.mode === "ignore") ? "strat-disabled-row" : "";
    return `<tr class="strat-row ${dimCls}" data-broker="${s.broker}">
      <td style="padding:5px 6px;white-space:nowrap;"><span class="alert-broker-dot" style="background:${bc};"></span>${s.broker}</td>
      <td style="padding:5px 6px;text-align:center;">
        <input type="checkbox" class="strat-cb" data-field="enabled" ${s.enabled ? "checked" : ""} />
      </td>
      <td style="padding:5px 6px;text-align:center;">
        <select class="strat-sel" data-field="mode" style="background:#0d1929;border:1px solid #1e3a5f;color:#e2e8f0;border-radius:3px;padding:2px 4px;font-size:11px;">
          <option value="too" ${s.mode==="too"?"selected":""}>too</option>
          <option value="queue" ${s.mode==="queue"?"selected":""}>queue</option>
          <option value="ignore" ${s.mode==="ignore"?"selected":""}>ignore</option>
        </select>
      </td>
      <td style="padding:5px 6px;"><input type="number" class="strat-inp" data-field="exposure" value="${s.exposure}" min="1" max="600" step="1" style="width:50px;" /></td>
      <td style="padding:5px 6px;"><input type="number" class="strat-inp" data-field="gain" value="${s.gain}" min="0" max="100" step="1" style="width:42px;" /></td>
      <td style="padding:5px 6px;"><input type="text" class="strat-inp" data-field="filter_cycle" value="${s.filter_cycle}" style="width:80px;" /></td>
      <td style="padding:5px 6px;"><input type="number" class="strat-inp" data-field="rapid_count" value="${s.rapid_count}" min="1" max="200" step="1" style="width:42px;" /></td>
      <td style="padding:5px 6px;text-align:center;"><input type="checkbox" class="strat-cb" data-field="do_af" ${s.do_af ? "checked" : ""} /></td>
      <td style="padding:5px 6px;text-align:center;"><input type="checkbox" class="strat-cb" data-field="do_guiding" ${s.do_guiding ? "checked" : ""} /></td>
      <td style="padding:5px 6px;text-align:center;"><input type="checkbox" class="strat-cb" data-field="do_center" ${s.do_center ? "checked" : ""} /></td>
      <td style="padding:5px 6px;"><input type="number" class="strat-inp" data-field="min_alt" value="${s.min_alt}" min="0" max="90" step="1" style="width:42px;" /></td>
      <td style="padding:5px 6px;"><input type="number" class="strat-inp" data-field="max_err_deg" value="${s.max_err_deg}" min="0.01" max="180" step="0.1" style="width:50px;" /></td>
      <td style="padding:5px 6px;text-align:center;" title="STDWeb target photometry (TNS-style). Off = transient detection.">
        <input type="checkbox" class="strat-cb" data-field="stdweb_use_target" ${s.stdweb_use_target ? "checked" : ""} />
      </td>
      <td style="padding:5px 6px;"><input type="text" class="strat-inp" data-field="notes" value="${(s.notes||"").replace(/"/g,"&quot;")}" style="width:160px;font-size:11px;" /></td>
      <td style="padding:5px 6px;"><button class="btn-small strat-reset-btn" data-broker="${s.broker}" title="Reset to defaults">↺</button></td>
    </tr>`;
  }).join("");

  tbody.querySelectorAll(".strat-cb, .strat-sel, .strat-inp").forEach(el => {
    const evt = el.tagName === "SELECT" || el.type === "checkbox" ? "change" : "change";
    el.addEventListener(evt, () => saveStrategyField(el));
  });
  tbody.querySelectorAll(".strat-reset-btn").forEach(btn => {
    btn.addEventListener("click", () => resetStrategy(btn.dataset.broker));
  });
}

async function saveStrategyField(el) {
  const row = el.closest("tr");
  const broker = row?.dataset.broker;
  if (!broker) return;
  const field = el.dataset.field;
  let value;
  if (el.type === "checkbox") value = el.checked ? 1 : 0;
  else value = el.value;

  const status = document.getElementById("strategy-save-status");
  try {
    const r = await fetch(`/api/alert-strategies/${encodeURIComponent(broker)}`, {
      method: "PUT",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ [field]: value }),
    });
    const d = await r.json();
    if (d.success) {
      const idx = _strategiesCache.findIndex(s => s.broker === broker);
      if (idx >= 0) _strategiesCache[idx] = d.strategy;
      if (field === "enabled" || field === "mode") renderStrategies();
      if (status) { status.textContent = `Saved ${broker}.${field}`; setTimeout(() => status.textContent = "", 2000); }
    } else {
      if (status) { status.textContent = `Error: ${d.error}`; status.style.color = "#f87171"; setTimeout(() => { status.textContent = ""; status.style.color = "#6b7280"; }, 3000); }
    }
  } catch (e) {
    if (status) { status.textContent = `Network error`; status.style.color = "#f87171"; setTimeout(() => { status.textContent = ""; status.style.color = "#6b7280"; }, 3000); }
  }
}

async function resetStrategy(broker) {
  try {
    const r = await fetch(`/api/alert-strategies/${encodeURIComponent(broker)}/reset`, { method: "POST" });
    const d = await r.json();
    if (d.success) {
      const idx = _strategiesCache.findIndex(s => s.broker === broker);
      if (idx >= 0) _strategiesCache[idx] = d.strategy;
      renderStrategies();
      const status = document.getElementById("strategy-save-status");
      if (status) { status.textContent = `${broker} reset to defaults`; setTimeout(() => status.textContent = "", 2000); }
    }
  } catch (e) { console.warn("[strategy-reset]", e); }
}

async function loadAlerts() {
  try {
    const r = await fetch("/api/alerts?limit=200");
    const d = await r.json();
    if (d.success) { _alertsCache = d.alerts || []; renderAlerts(); }
  } catch (e) { console.warn("[alerts-ui]", e); }
}

function renderAlerts() {
  const tbody = document.getElementById("alerts-tbody");
  const table = document.getElementById("alerts-table");
  const empty = document.getElementById("alerts-empty");
  const badge = document.getElementById("alerts-badge");
  if (!tbody) return;

  const brokerFilter = document.getElementById("alerts-filter-broker")?.value || "";
  const actionFilter = document.getElementById("alerts-filter-action")?.value || "";

  let rows = _alertsCache;
  if (brokerFilter) rows = rows.filter(a => a.broker === brokerFilter);
  if (actionFilter) rows = rows.filter(a => (a.action || "").startsWith(actionFilter));
  // Hide "ignored" alerts (disabled-strategy drops) unless explicitly requested
  else rows = rows.filter(a => a.action !== "ignored");

  const queuedCount = _alertsCache.filter(a => a.action === "queued" || a.action === "too-observed").length;
  if (badge) {
    if (queuedCount > 0) { badge.textContent = queuedCount; badge.style.display = "inline-block"; }
    else { badge.style.display = "none"; }
  }

  if (!rows.length) {
    table.style.display = "none";
    empty.style.display = "";
    empty.textContent = _alertsCache.length ? "No alerts match the current filter." : "No alerts received yet.";
    return;
  }
  table.style.display = "";
  empty.style.display = "none";

  tbody.innerHTML = rows.map(a => {
    const bc = BROKER_COLORS[a.broker] || "#6b7280";
    const al = ACTION_LABELS[a.action] || ACTION_LABELS.rejected;
    const actionText = a.action === "rejected" && a.action_reason
      ? `REJ: ${a.action_reason.split(":")[0]}`
      : a.action === "no-position" && a.action_reason
        ? `NO POS`
        : al.text;
    const age = a.received_at ? timeAgo(a.received_at) : "";
    return `<tr class="alert-row" data-id="${a.id}">
      <td style="padding:5px 6px;white-space:nowrap;" title="${a.received_at || ""}">${age}</td>
      <td style="padding:5px 6px;"><span class="alert-broker-dot" style="background:${bc};"></span>${a.broker || "?"}</td>
      <td style="padding:5px 6px;color:#93c5fd;">${a.trigger_id || "—"}</td>
      <td style="padding:5px 6px;text-align:right;font-family:monospace;">${fmtDeg(a.ra)}</td>
      <td style="padding:5px 6px;text-align:right;font-family:monospace;">${fmtDeg(a.dec)}</td>
      <td style="padding:5px 6px;text-align:right;">${a.err_deg != null ? (a.err_deg < 0.01667 ? (a.err_deg*3600).toFixed(0)+"″" : (a.err_deg*60).toFixed(1)+"′") : "—"}</td>
      <td style="padding:5px 6px;text-align:right;">${fmtDeg1(a.alt_now)}</td>
      <td style="padding:5px 6px;text-align:right;">${fmtDeg1(a.moon_sep)}</td>
      <td style="padding:5px 6px;text-align:center;"><span class="alert-action-pill ${al.cls}">${actionText}</span></td>
      <td style="padding:5px 6px;"><button class="btn-small alert-detail-btn" data-id="${a.id}">Detail</button></td>
    </tr>`;
  }).join("");

  tbody.querySelectorAll(".alert-detail-btn").forEach(btn => {
    btn.addEventListener("click", (e) => { e.stopPropagation(); openAlertDetail(Number(btn.dataset.id)); });
  });
}

function timeAgo(utcStr) {
  const d = new Date(utcStr + (utcStr.endsWith("Z") ? "" : "Z"));
  const sec = Math.floor((Date.now() - d.getTime()) / 1000);
  if (sec < 60)    return sec + "s ago";
  if (sec < 3600)  return Math.floor(sec / 60) + "m ago";
  if (sec < 86400) return Math.floor(sec / 3600) + "h ago";
  return Math.floor(sec / 86400) + "d ago";
}

async function openAlertDetail(id) {
  const modal = document.getElementById("alert-detail-modal");
  const title = document.getElementById("alert-detail-title");
  const tbody = document.getElementById("alert-detail-tbody");
  const rawEl = document.getElementById("alert-detail-raw");
  const qBtn  = document.getElementById("alert-detail-queue");
  const iBtn  = document.getElementById("alert-detail-ignore");
  if (!modal) return;
  try {
    const r = await fetch(`/api/alerts/${id}`);
    const d = await r.json();
    if (!d.success || !d.alert) return;
    const a = d.alert;
    title.textContent = `${a.broker || "?"} — ${a.trigger_id || "Alert #" + a.id}`;
    const kv = [
      ["Received", a.received_at], ["Event time", a.event_time], ["Broker", a.broker],
      ["Topic", a.topic], ["Trigger ID", a.trigger_id], ["Classification", a.classification],
      ["RA (deg)", fmtDeg(a.ra, 5)], ["Dec (deg)", fmtDeg(a.dec, 5)],
      ["Error radius", a.err_deg != null ? (a.err_deg*60).toFixed(2)+"′" : "—"],
      ["Altitude", fmtDeg1(a.alt_now)], ["Moon sep", fmtDeg1(a.moon_sep)],
      ["Action", a.action], ["Reason", a.action_reason || "—"],
      ["AstroColibri ID", a.colibri_id || "—"],
    ];
    tbody.innerHTML = kv.map(([k,v]) =>
      `<tr><td style="color:#6b7280;padding:4px 10px 4px 0;white-space:nowrap;">${k}</td><td style="color:#e2e8f0;padding:4px 0;">${v||"—"}</td></tr>`
    ).join("");
    try { rawEl.textContent = JSON.stringify(JSON.parse(a.raw), null, 2); }
    catch { rawEl.textContent = a.raw || "(empty)"; }
    const canQueue = (a.action !== "queued" && a.action !== "manual" && a.action !== "too-observed") && a.ra != null;
    qBtn.style.display = canQueue ? "" : "none";
    iBtn.style.display = a.action !== "ignored" ? "" : "none";
    qBtn.onclick = async () => { await fetch(`/api/alerts/${id}/queue`, { method: "POST" }); modal.style.display = "none"; loadAlerts(); };
    iBtn.onclick = async () => { await fetch(`/api/alerts/${id}/ignore`, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ reason: "operator-dismissed" }) }); modal.style.display = "none"; loadAlerts(); };
    modal.style.display = "";
  } catch (e) { console.warn("[alert-detail]", e); }
}

document.getElementById("alert-detail-close")?.addEventListener("click", () => { document.getElementById("alert-detail-modal").style.display = "none"; });
document.getElementById("alert-detail-modal")?.addEventListener("click", (e) => { if (e.target === e.currentTarget) e.currentTarget.style.display = "none"; });
document.getElementById("alerts-refresh")?.addEventListener("click", loadAlerts);
document.getElementById("alerts-filter-broker")?.addEventListener("change", renderAlerts);
document.getElementById("alerts-filter-action")?.addEventListener("change", renderAlerts);

function startAlertsAutoRefresh() { stopAlertsAutoRefresh(); const cb = document.getElementById("alerts-auto-refresh"); if (cb?.checked) _alertsAutoTimer = setInterval(loadAlerts, 30000); }
function stopAlertsAutoRefresh() { if (_alertsAutoTimer) { clearInterval(_alertsAutoTimer); _alertsAutoTimer = null; } }
document.getElementById("alerts-auto-refresh")?.addEventListener("change", (e) => { if (e.target.checked) startAlertsAutoRefresh(); else stopAlertsAutoRefresh(); });

loadAlerts();
loadAlertConfig();
loadStrategies();
startAlertsAutoRefresh();
