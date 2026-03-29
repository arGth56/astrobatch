const form = document.getElementById("config-form");
const connectionStatus = document.getElementById("connection-status");
const logEl = document.getElementById("log");
const refreshButton = document.getElementById("refresh-status");
const connectButtons = Array.from(document.querySelectorAll("button[data-device]"));

const tabDevicesBtn  = document.getElementById("tab-devices");
const tabActionsBtn  = document.getElementById("tab-actions");
const tabDomeBtn     = document.getElementById("tab-dome");
const tabTargetBtn   = document.getElementById("tab-target");
const tabTodoBtn      = document.getElementById("tab-todo");
const tabPipelineBtn  = document.getElementById("tab-pipeline");
const tabHistoryBtn   = document.getElementById("tab-history");
const panelDevices    = document.getElementById("panel-devices");
const panelActions    = document.getElementById("panel-actions");
const panelDome       = document.getElementById("panel-dome");
const panelTarget     = document.getElementById("panel-target");
const panelTodo       = document.getElementById("panel-todo");
const panelPipeline   = document.getElementById("panel-pipeline");
const panelHistory    = document.getElementById("panel-history");

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
let ocsConnected = false;
let observerLocation = { lat: null, lon: null, elevation: null };
let lastTargetCoords = { raDeg: null, decDeg: null };
let lastTnsTarget = null;

// Sequence polling
let seqPollInterval    = null;
let seqLogRenderedCount = 0;

// Drag-and-drop state (module-level so re-renders during polling don't lose it)
let _seqDragSrc        = null;
let _seqDragging       = false;
let _seqDragInsertAfter = false; // last hover half from dragover, used by drop

function currentConfig() {
  return {
    host: document.getElementById("host").value.trim(),
    port: document.getElementById("port").value.trim(),
    protocol: document.getElementById("protocol").value,
  };
}

function setLog(value) {
  logEl.textContent = typeof value === "string" ? value : JSON.stringify(value, null, 2);
}

function setActiveTab(tab) {
  const tabs = { devices: false, actions: false, dome: false, target: false, todo: false, pipeline: false, history: false };
  tabs[tab] = true;
  tabDevicesBtn.classList.toggle("active",  tabs.devices);
  tabActionsBtn.classList.toggle("active",  tabs.actions);
  tabDomeBtn.classList.toggle("active",     tabs.dome);
  tabTargetBtn.classList.toggle("active",   tabs.target);
  tabTodoBtn.classList.toggle("active",     tabs.todo);
  tabPipelineBtn.classList.toggle("active", tabs.pipeline);
  tabHistoryBtn.classList.toggle("active",  tabs.history);
  panelDevices.classList.toggle("active",   tabs.devices);
  panelActions.classList.toggle("active",   tabs.actions);
  panelDome.classList.toggle("active",      tabs.dome);
  panelTarget.classList.toggle("active",    tabs.target);
  panelTodo.classList.toggle("active",      tabs.todo);
  panelPipeline.classList.toggle("active",  tabs.pipeline);
  panelHistory.classList.toggle("active",   tabs.history);

  if (tab === "todo") {
    startSeqPolling();
    initNightPlan();
    refreshNightPlan();
  } else {
    stopSeqPolling();
  }
  if (tab === "pipeline") {
    startPipelinePolling();
  } else {
    stopPipelinePolling();
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
  if (ocsConnected) {
    el.textContent = `OCS ✓ ${currentOcsHost()}`;
    el.className = "ocs-status-indicator online";
  } else {
    el.textContent = `OCS ✗ ${currentOcsHost()} — not reachable`;
    el.className = "ocs-status-indicator offline";
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
    const payload = await postJson("/api/ocs/roof", { ocsHost: host, command });
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
  const filters = equipment?.FilterWheel?.AvailableFilters || [];
  const selectedId = equipment?.FilterWheel?.SelectedFilter?.Id;
  filterSelect.innerHTML = "";
  for (const filter of filters) {
    const option = document.createElement("option");
    option.value = String(filter.Id);
    const label = filter.Name && filter.Name.trim() ? filter.Name : `Filter ${filter.Id}`;
    option.textContent = `${label} (${filter.Id})`;
    if (selectedId === filter.Id) {
      option.selected = true;
    }
    filterSelect.appendChild(option);
  }

  if (!filters.length) {
    const option = document.createElement("option");
    option.value = "";
    option.textContent = "No filters available";
    filterSelect.appendChild(option);
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

  if (ninaResult.status === "fulfilled") {
    const payload = ninaResult.value;
    connectionStatus.textContent = payload.success
      ? `Connected to ${payload.target.protocol}://${payload.target.host}:${payload.target.port}`
      : `Connection failed (${payload.status || "unknown status"})`;
    setDeviceStatuses(payload.devices || {});
    lastEquipment = payload.equipment || null;
    updateFilterOptions(lastEquipment);

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

async function changeFilter() {
  const filterId = Number(filterSelect.value);
  setLog(`Changing filter to ${filterId}...`);
  try {
    const payload = await postJson("/api/nina/actions/filterwheel/change", {
      ...currentConfig(),
      filterId,
    });
    await refreshStatus();
    setLog(payload);
  } catch (error) {
    setLog({ success: false, action: "filter-change", error: error.message });
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
  await refreshStatus();
});

tabDevicesBtn.addEventListener("click", () => setActiveTab("devices"));
tabActionsBtn.addEventListener("click", () => setActiveTab("actions"));
tabDomeBtn.addEventListener("click", () => { setActiveTab("dome"); refreshOcsStatus(); loadOcsHistory(); });
tabTargetBtn.addEventListener("click", () => setActiveTab("target"));

document.getElementById("target-lookup-form").addEventListener("submit", async (event) => {
  event.preventDefault();
  await lookupTarget();
});

document.getElementById("dome-ocs-config-form").addEventListener("submit", async (event) => {
  event.preventDefault();
  syncOcsHostInputs(document.getElementById("dome-ocs-host").value.trim() || ocsHost);
  await refreshOcsStatus();
});

document.getElementById("ocs-connect-form").addEventListener("submit", async (event) => {
  event.preventDefault();
  syncOcsHostInputs(document.getElementById("ocs-host-top").value.trim() || ocsHost);
  await refreshOcsStatus();
});

document.getElementById("dome-open").addEventListener("click", async () => {
  if (confirm("Open the roof?")) await roofCommand("open");
});
document.getElementById("dome-close").addEventListener("click", async () => {
  // Check mount park state before allowing roof close
  const mount = lastEquipment?.Mount;
  if (mount?.Connected && mount?.AtPark !== true) {
    const detail = `Mount is ${mount.Slewing ? "slewing" : "not parked"} (AtPark=${mount.AtPark}).`;
    alert(`⛔ Cannot close roof — ${detail}\n\nPark the mount first.`);
    return;
  }
  if (!mount?.Connected) {
    const ok = confirm("⚠ Mount status unknown (NINA not connected).\nCannot verify park state.\n\nClose roof anyway?");
    if (!ok) return;
  } else {
    if (!confirm("Mount is parked ✓ — close the roof?")) return;
  }
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

focuserInBtn.addEventListener("click", async () => moveFocuser("in"));
focuserOutBtn.addEventListener("click", async () => moveFocuser("out"));

async function loadDefaults() {
  try {
    const response = await fetch("/api/config/defaults");
    const defaults = await response.json();
    document.getElementById("host").value = defaults.host || "192.168.1.174";
    document.getElementById("port").value = defaults.port || "1888";
    document.getElementById("protocol").value = defaults.protocol || "http";
    ocsHost = defaults.ocsHost || "192.168.1.220";
    syncOcsHostInputs(ocsHost);
    if (defaults.tns?.botId) document.getElementById("tns-bot-id").value = defaults.tns.botId;
    if (defaults.tns?.botName) document.getElementById("tns-bot-name").value = defaults.tns.botName;
    if (defaults.tns?.hasApiKey) document.getElementById("tns-api-key").placeholder = "Loaded from server env";
  } catch {
    document.getElementById("host").value = "192.168.1.174";
    document.getElementById("port").value = "1888";
    document.getElementById("protocol").value = "http";
    syncOcsHostInputs("192.168.1.220");
  }
}

// ── Sequence / ToDo tab ──────────────────────────────────────────────────────

function seqConfig() {
  const filterRaw = document.getElementById("seq-filters").value;
  return {
    duration:       Number(document.getElementById("seq-duration").value)       || 120,
    gain:           Number(document.getElementById("seq-gain").value)           || 10,
    count:          Number(document.getElementById("seq-count").value)          || 10,
    filters:        filterRaw.split(",").map(s => s.trim()).filter(Boolean),
    solveEnabled:              document.getElementById("seq-solve-enable").checked,
    solveExp:                  Number(document.getElementById("seq-solve-exp").value)                  || 5,
    solveThreshold:            Number(document.getElementById("seq-solve-threshold").value)            || 60,
    frameCheckEnabled:         document.getElementById("seq-frame-check-enable").checked,
    frameCheckThresholdArcmin: Number(document.getElementById("seq-frame-check-threshold").value)     || 5,
    manualMode:                document.getElementById("seq-manual-mode").checked,
  };
}

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
      const raH     = (t.raDeg || 0) / 15;
      const raHH    = Math.floor(raH);
      const raM     = String(Math.floor((raH - raHH) * 60)).padStart(2, "0");
      const decSign = (t.decDeg || 0) >= 0 ? "+" : "";
      li.innerHTML = `
        <span class="seq-queue-num">${i + 1}</span>
        <span class="seq-queue-name">${t.name}</span>
        <span class="seq-queue-coords">α&nbsp;${raHH}h${raM}m &nbsp;δ&nbsp;${decSign}${(t.decDeg || 0).toFixed(1)}°</span>
        ${t.done
          ? `<span class="seq-queue-done-badge">✓ Done</span>`
          : `<button class="seq-queue-remove" data-idx="${i}" type="button">×</button>`}
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
      _seqDragInsertAfter = e.offsetY >= li.offsetHeight / 2;
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
    });
  }

  // Wait-time edits: update the wait item's time on change
  for (const inp of list.querySelectorAll(".seq-queue-wait-input")) {
    inp.addEventListener("change", async () => {
      const idx = parseInt(inp.dataset.idx);
      const newTime = inp.value; // HH:MM
      await fetch(`/api/sequence/queue/wait/${idx}`, {
        method: "PATCH",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ waitTime: newTime }),
      });
      await refreshSeqState();
    });
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
    renderSeqQueue(state.queue);
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
  queued:      "Queued…",
  scanning:    "Scanning FITS files…",
  splitting:   "Copying to work dir…",
  calibrating: "Calibrating + stacking (Siril)…",
  solving:     "Plate solving…",
  uploading:   "Uploading to STDWeb…",
  done:        "Complete",
  error:       "Error",
};

// ── NAS date / target picker ──────────────────────────────────────────────────

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

  // ── Step 2: when a night is chosen, scan that folder for targets ──────────
  dateSel.addEventListener("change", async () => {
    const date = dateSel.value;
    targetSel.innerHTML = '<option value="">⏳ Scanning night…</option>';
    targetSel.style.display = "";
    targetSel.disabled = true;
    fitsInput.value = "";

    if (!date) {
      targetSel.style.display = "none";
      targetSel.disabled = false;
      return;
    }

    let targets = [];
    try {
      const res  = await fetch(`/api/pipeline/nas-dates/${date}`);
      const data = await res.json();
      targets = data.targets || [];
    } catch {
      targetSel.innerHTML = '<option value="">Scan failed</option>';
      targetSel.disabled = false;
      return;
    }

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
    const step = JOB_STATUS_STEPS[running?.status] || running?.status || "Queued…";
    const queueSuffix = queuedJobs.length
      ? ` — ${queuedJobs.length} job${queuedJobs.length > 1 ? "s" : ""} waiting`
      : "";
    statTxt.textContent = running
      ? `Job #${running.id} — ${step}${queueSuffix}`
      : `${queuedJobs.length} job${queuedJobs.length > 1 ? "s" : ""} queued…`;
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

      // FITS download button (only when result is done)
      const fitsBtn = r.status === "done"
        ? `<a href="/api/pipeline/result/${r.id}/download"
              class="btn-small"
              style="margin-left:4px;color:#a78bfa;border:1px solid #5b21b6;border-radius:4px;
                     padding:2px 6px;font-size:12px;text-decoration:none;display:inline-block;"
              title="Download stacked FITS (${r.obs_date} · ${r.target} · ${r.filter} · ${r.exposure}s)">⭐</a>` : "";

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
        <td>${link}${framesBtn}${photBtn}${fitsBtn}</td>
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

      // Bust cache if it has no subtraction data — re-fetch to get the full result
      if (_photCache[taskId] && !_photCache[taskId].sub && !_photCache[taskId].sub_ul) {
        delete _photCache[taskId];
      }

      // Fetch if not yet cached
      const cell = photRow.querySelector("td");
      if (!_photCache[taskId]) {
        btn.disabled = true;
        btn.textContent = "⏳";
        cell.textContent = "Fetching…";
        photRow.style.display = "";

        try {
          const resp = await fetch(`/api/stdweb/task/${taskId}/photometry`);
          _photCache[taskId] = await resp.json();
        } catch (e) {
          cell.textContent = `Error: ${e.message}`;
          btn.disabled = false; btn.textContent = "📊";
          return;
        }
        btn.disabled = false; btn.textContent = "📊";
      }

      const resultObj = jobs.flatMap(j => j.results || []).find(x => String(x.id) === String(rid));
      _renderPhotRow(photRow, _photCache[taskId], resultObj);
      photRow.style.display = "";
      btn.style.opacity = "1";
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
// Zooms into the image under the cursor, clipped to the thumbnail rectangle.
// Scale factor brings it to 1 image-pixel = 1 screen pixel (naturalWidth/displayWidth).

(function setupFrameZoom() {
  let activeImg = null;

  function applyZoom(img, e) {
    const rect = img.getBoundingClientRect();
    // Cursor position relative to the displayed image
    const cx = Math.max(0, Math.min(e.clientX - rect.left,  rect.width));
    const cy = Math.max(0, Math.min(e.clientY - rect.top,   rect.height));
    // Scale so that 1 image pixel = 1 screen pixel
    const scale = (img.naturalWidth || 800) / rect.width;
    img.style.transformOrigin = `${cx}px ${cy}px`;
    img.style.transform = `scale(${scale})`;
  }

  function resetZoom(img) {
    img.style.transform = "";
    img.style.transformOrigin = "center center";
  }

  document.addEventListener("mousemove", (e) => {
    const img = e.target.closest(".frame-zoom-img");
    if (!img) {
      if (activeImg) { resetZoom(activeImg); activeImg = null; }
      return;
    }
    activeImg = img;
    applyZoom(img, e);
  });

  document.addEventListener("mouseleave", (e) => {
    if (activeImg) { resetZoom(activeImg); activeImg = null; }
  }, true);

  // Reset when cursor leaves the wrapper div
  document.addEventListener("mouseout", (e) => {
    const wrap = e.target.closest(".frame-thumb-wrap");
    if (!wrap) return;
    if (e.relatedTarget && wrap.contains(e.relatedTarget)) return;
    const img = wrap.querySelector(".frame-zoom-img");
    if (img) { resetZoom(img); activeImg = null; }
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
  const fits_dir = document.getElementById("pipeline-fits-dir").value.trim();
  const target   = document.getElementById("pipeline-target").value.trim();
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
    });
    setLog(result);
    document.getElementById("pipeline-target").value = "";
    document.getElementById("pipeline-date-select").value = "";
    loadPipelineJobs();
  } catch (err) {
    setLog({ success: false, error: err.message });
  } finally {
    btn.disabled = false;
    btn.textContent = "▶ Trigger";
  }
});

setActiveTab("devices");
loadDefaults().then(refreshStatus);
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

  const hasSky      = history.some(r => r.sky      != null);
  const hasPressure = history.some(r => r.pressure != null);

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
