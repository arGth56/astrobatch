const form = document.getElementById("config-form");
const connectionStatus = document.getElementById("connection-status");
const logEl = document.getElementById("log");
const refreshButton = document.getElementById("refresh-status");
const connectButtons = Array.from(document.querySelectorAll("button[data-device]"));

const tabDevicesBtn  = document.getElementById("tab-devices");
const tabActionsBtn  = document.getElementById("tab-actions");
const tabDomeBtn     = document.getElementById("tab-dome");
const tabTargetBtn   = document.getElementById("tab-target");
const tabTodoBtn     = document.getElementById("tab-todo");
const tabPipelineBtn  = document.getElementById("tab-pipeline");
const tabNightPlanBtn = document.getElementById("tab-nightplan");
const panelDevices    = document.getElementById("panel-devices");
const panelActions    = document.getElementById("panel-actions");
const panelDome       = document.getElementById("panel-dome");
const panelTarget     = document.getElementById("panel-target");
const panelTodo       = document.getElementById("panel-todo");
const panelPipeline   = document.getElementById("panel-pipeline");
const panelNightPlan  = document.getElementById("panel-nightplan");

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
let ocsHost = "ocs.local";
let observerLocation = { lat: null, lon: null, elevation: null };
let lastTargetCoords = { raDeg: null, decDeg: null };
let lastTnsTarget = null;

// Sequence polling
let seqPollInterval = null;
let seqLogRenderedCount = 0;

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
  const tabs = { devices: false, actions: false, dome: false, target: false, todo: false, pipeline: false, nightplan: false };
  tabs[tab] = true;
  tabDevicesBtn.classList.toggle("active",   tabs.devices);
  tabActionsBtn.classList.toggle("active",   tabs.actions);
  tabDomeBtn.classList.toggle("active",      tabs.dome);
  tabTargetBtn.classList.toggle("active",    tabs.target);
  tabTodoBtn.classList.toggle("active",      tabs.todo);
  tabPipelineBtn.classList.toggle("active",  tabs.pipeline);
  tabNightPlanBtn.classList.toggle("active", tabs.nightplan);
  panelDevices.classList.toggle("active",    tabs.devices);
  panelActions.classList.toggle("active",    tabs.actions);
  panelDome.classList.toggle("active",       tabs.dome);
  panelTarget.classList.toggle("active",     tabs.target);
  panelTodo.classList.toggle("active",       tabs.todo);
  panelPipeline.classList.toggle("active",   tabs.pipeline);
  panelNightPlan.classList.toggle("active",  tabs.nightplan);

  if (tab === "nightplan") {
    initNightPlan();
    refreshNightPlan();
  }

  if (tab === "todo") {
    startSeqPolling();
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
  return document.getElementById("dome-ocs-host").value.trim() || ocsHost;
}

function setOcsStatus(roof = {}) {
  // Dome tab detail view
  const roofEl = document.getElementById("dome-roof-status");
  const safeEl = document.getElementById("dome-safe");

  const roofStatus = (roof.status || "Unknown").replace(/<[^>]*>/g, "");
  roofEl.textContent = roofStatus;
  const isOpen = /open/i.test(roofStatus);
  const isClosed = /close/i.test(roofStatus);
  roofEl.className = isOpen ? "online" : isClosed ? "offline" : "";

  const safe = roof.safe || "Unknown";
  safeEl.textContent = safe;
  safeEl.className = /safe/i.test(safe) ? "online" : "offline";

  document.getElementById("dome-rain").textContent = roof.rain || "–";
  document.getElementById("dome-temp").textContent = roof.temp || "–";
  document.getElementById("dome-humidity").textContent = roof.humidity || "–";
  document.getElementById("dome-pressure").textContent = roof.pressure || "–";

  // Devices tab — Dome row (roof open/closed)
  ocsDomeEl.textContent = roofStatus;
  ocsDomeEl.className = isOpen ? "online" : isClosed ? "offline" : "";
  if (ocsDomeDetailEl) {
    ocsDomeDetailEl.textContent = roof.rain ? `Rain: ${roof.rain}` : "";
  }

  // Devices tab — Safety Monitor row
  ocsSafetyEl.textContent = safe;
  ocsSafetyEl.className = /safe/i.test(safe) ? "online" : "offline";
  if (ocsSafetyDetailEl) {
    const details = [];
    if (roof.temp && roof.temp !== "Invalid") details.push(`${roof.temp}°C`);
    if (roof.humidity && roof.humidity !== "Invalid") details.push(`${roof.humidity}% RH`);
    ocsSafetyDetailEl.textContent = details.join("  ");
  }
}

async function refreshOcsStatus() {
  const host = currentOcsHost();
  setLog(`Fetching OCS status from ${host}...`);
  try {
    const payload = await postJson("/api/ocs/status", { ocsHost: host });
    setOcsStatus(payload.roof || {});
    setLog(payload);
  } catch (error) {
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
    setOcsStatus(ocsResult.value.roof || {});
  } else {
    ocsDomeEl.textContent = "Unreachable";
    ocsDomeEl.className = "offline";
    ocsSafetyEl.textContent = "Unreachable";
    ocsSafetyEl.className = "offline";
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
tabDomeBtn.addEventListener("click", () => { setActiveTab("dome"); refreshOcsStatus(); });
tabTargetBtn.addEventListener("click", () => setActiveTab("target"));

document.getElementById("target-lookup-form").addEventListener("submit", async (event) => {
  event.preventDefault();
  await lookupTarget();
});

document.getElementById("dome-ocs-config-form").addEventListener("submit", async (event) => {
  event.preventDefault();
  await refreshOcsStatus();
});

document.getElementById("dome-open").addEventListener("click", async () => {
  if (confirm("Open the roof?")) await roofCommand("open");
});
document.getElementById("dome-close").addEventListener("click", async () => {
  if (confirm("Close the roof?")) await roofCommand("close");
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
    ocsHost = defaults.ocsHost || "ocs.local";
    document.getElementById("dome-ocs-host").value = ocsHost;
    if (defaults.tns?.botId) document.getElementById("tns-bot-id").value = defaults.tns.botId;
    if (defaults.tns?.botName) document.getElementById("tns-bot-name").value = defaults.tns.botName;
    if (defaults.tns?.hasApiKey) document.getElementById("tns-api-key").placeholder = "Loaded from server env";
  } catch {
    document.getElementById("host").value = "192.168.1.174";
    document.getElementById("port").value = "1888";
    document.getElementById("protocol").value = "http";
    document.getElementById("dome-ocs-host").value = "ocs.local";
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
    solveEnabled:   document.getElementById("seq-solve-enable").checked,
    solveExp:       Number(document.getElementById("seq-solve-exp").value)       || 5,
    solveThreshold: Number(document.getElementById("seq-solve-threshold").value) || 60,
    manualMode:     document.getElementById("seq-manual-mode").checked,
  };
}

function renderSeqQueue(queue) {
  const list  = document.getElementById("seq-queue-list");
  const empty = document.getElementById("seq-queue-empty");
  list.innerHTML = "";

  if (!queue.length) {
    empty.style.display = "block";
    return;
  }
  empty.style.display = "none";

  queue.forEach((t, i) => {
    const raH  = t.raDeg / 15;
    const raHH = Math.floor(raH);
    const raM  = String(Math.floor((raH - raHH) * 60)).padStart(2, "0");
    const decSign = t.decDeg >= 0 ? "+" : "";

    const li = document.createElement("li");
    li.className = "seq-queue-item" + (t.done ? " done" : "");

    li.innerHTML = `
      <span class="seq-queue-num">${i + 1}</span>
      <span class="seq-queue-name">${t.name}</span>
      <span class="seq-queue-coords">α&nbsp;${raHH}h${raM}m &nbsp;δ&nbsp;${decSign}${t.decDeg.toFixed(1)}°</span>
      ${t.done
        ? `<span class="seq-queue-done-badge">✓ Done</span>`
        : `<button class="seq-queue-remove" data-idx="${i}" type="button">×</button>`
      }
    `;
    list.appendChild(li);
  });

  for (const btn of list.querySelectorAll(".seq-queue-remove")) {
    btn.addEventListener("click", async () => {
      await fetch(`/api/sequence/queue/${btn.dataset.idx}`, { method: "DELETE" });
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
    runBtn.disabled = false;
    abortBtn.disabled = true;
    nextBtn.disabled = true;
  } else {
    const allDone = state.queue.length > 0 && state.queue.every(t => t.done);
    dot.className = "seq-dot " + (allDone ? "seq-dot-done" : "seq-dot-idle");
    label.textContent = allDone ? "All targets complete ✓" : "Idle";
    runBtn.disabled = false;
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
tabNightPlanBtn.addEventListener("click", () => setActiveTab("nightplan"));

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

let nasTree = [];

async function loadNasDates() {
  const dateSel   = document.getElementById("pipeline-date-select");
  const targetSel = document.getElementById("pipeline-target-select");
  const fitsInput = document.getElementById("pipeline-fits-dir");
  const targetInput = document.getElementById("pipeline-target");

  dateSel.innerHTML = '<option value="">Loading…</option>';

  try {
    const res  = await fetch("/api/pipeline/nas-dates");
    const data = await res.json();
    nasTree = data.tree || [];

    dateSel.innerHTML = '<option value="">— pick a date —</option>';
    for (const entry of nasTree) {
      const opt = document.createElement("option");
      opt.value = entry.date;
      opt.textContent = entry.date + (entry.targets.length ? ` (${entry.targets.length} target${entry.targets.length > 1 ? "s" : ""})` : "");
      dateSel.appendChild(opt);
    }
    if (!nasTree.length) {
      dateSel.innerHTML = '<option value="">No data on NAS</option>';
    }
  } catch {
    dateSel.innerHTML = '<option value="">NAS unavailable</option>';
  }

  // When date changes, populate target selector
  dateSel.addEventListener("change", () => {
    const entry = nasTree.find((e) => e.date === dateSel.value);
    targetSel.innerHTML = '<option value="">— pick a target —</option>';
    if (!entry) { targetSel.style.display = "none"; return; }

    const hasRealTargets = entry.targets.some((t) => t.name !== "SNAPSHOT");
    for (const t of entry.targets) {
      // Hide the SNAPSHOT fallback folder when real per-target folders exist
      if (t.name === "SNAPSHOT" && hasRealTargets) continue;
      const opt = document.createElement("option");
      opt.value = t.path;
      opt.dataset.name = t.name;
      opt.textContent = t.name === "SNAPSHOT" ? "SNAPSHOT (unsorted)" : t.name;
      targetSel.appendChild(opt);
    }
    targetSel.style.display = "";

    // Auto-select if only one target
    if (entry.targets.length === 1) {
      targetSel.value = entry.targets[0].path;
      fitsInput.value = entry.targets[0].path;
      if (entry.targets[0].name !== "SNAPSHOT") {
        targetInput.value = entry.targets[0].name;
      }
    }
  });

  targetSel.addEventListener("change", () => {
    if (!targetSel.value) return;
    fitsInput.value = targetSel.value;
    const name = targetSel.selectedOptions[0]?.dataset.name || "";
    if (name && name !== "SNAPSHOT") {
      targetInput.value = name;
    }
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
  navigator.clipboard.writeText(text).then(() => {
    const btn = document.getElementById("pipeline-log-copy");
    btn.textContent = "✓ Copied";
    setTimeout(() => { btn.textContent = "⎘ Copy"; }, 1500);
  });
});

// ── Render jobs ───────────────────────────────────────────────────────────────

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

  // Status bar: show latest running job
  const running = jobs.find((j) => !["done", "error"].includes(j.status));
  if (running) {
    statBar.style.display = "block";
    const step = JOB_STATUS_STEPS[running.status] || running.status;
    statTxt.textContent = `Job #${running.id} — ${step}`;
    document.getElementById("pipeline-log-btn").dataset.jobId = running.id;
    document.getElementById("pipeline-log-btn").style.display = "";
  } else {
    statBar.style.display = "none";
  }

  if (!jobs.length) {
    table.style.display = "none";
    empty.style.display = "block";
    return;
  }
  table.style.display = "";
  empty.style.display = "none";

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

      const rtr = document.createElement("tr");
      rtr.dataset.resultId = r.id;
      rtr.style.background = "#0d1929";
      rtr.innerHTML = `
        <td></td>
        <td style="font-size:12px;padding-left:20px;color:#93c5fd;">
          ${r.target ? `<span style="color:#a78bfa;">${r.target}</span> · ` : ""}${r.filter || "?"} · ${r.n_frames || "?"}×${r.exposure || "?"}s
        </td>
        <td></td>
        <td><span class="job-status ${rcls}"${retitle} style="font-size:11px;">${rlabel}</span></td>
        <td>${link}${framesBtn}</td>
        <td></td>
      `;
      tbody.appendChild(rtr);

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

async function loadPipelineJobs() {
  try {
    const res  = await fetch("/api/pipeline/jobs");
    const data = await res.json();
    renderPipelineJobs(data.jobs || []);
  } catch { /* ignore */ }
}

async function loadStackedResults() {
  try {
    const res   = await fetch("/api/data/results");
    const items = await res.json();
    const tbody = document.getElementById("results-tbody");
    const empty = document.getElementById("results-empty");
    const table = document.getElementById("results-table");
    if (!items.length) { table.style.display = "none"; empty.style.display = ""; return; }
    table.style.display = "";
    empty.style.display = "none";
    tbody.innerHTML = items.map(r => `
      <tr>
        <td>${r.date}</td>
        <td>${r.target}</td>
        <td>${r.filter}</td>
        <td>${r.exp}</td>
        <td>${r.size_mb} MB</td>
        <td><a href="${r.url}" download style="color:#60a5fa;">⬇ res.fit</a></td>
      </tr>
      ${r.preview_url ? `<tr>
        <td colspan="6" style="padding:8px 12px 16px;">
          <img src="${r.preview_url}" alt="${r.target} ${r.filter}"
               style="max-width:100%;border-radius:6px;cursor:zoom-in;display:block;"
               onclick="this.style.maxWidth=this.style.maxWidth==='100%'?'none':'100%'" />
        </td>
      </tr>` : ""}`).join("");
  } catch { /* ignore */ }
}

function startPipelinePolling() {
  if (pipelinePollInterval) return;
  loadNasDates();
  loadPipelineJobs();
  loadStackedResults();
  pipelinePollInterval = setInterval(() => { loadPipelineJobs(); loadStackedResults(); }, 4000);
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
    const result = await postJson("/api/pipeline/trigger", {
      fits_dir,
      target: target || null,
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
