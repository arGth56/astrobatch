/**
 * Night Plan — multi-target altitude chart, scheduling, and optimize-order.
 * Depends on astronomy.js being loaded first (window.Astronomy).
 */

const NP_COLORS = [
  "#60a5fa", "#34d399", "#f59e0b", "#f87171",
  "#a78bfa", "#fb923c", "#22d3ee", "#e879f9",
  "#a3e635", "#2dd4bf",
];

let _npInitialized      = false;
let _npQueue            = [];
let _npSchedule         = [];
let _wtData             = null; // latest 7Timer ASTRO data

// ── Drag state ────────────────────────────────────────────────────────────────
let _npDrag             = null;  // { slotIdx, mouseStartX, deltaH }
let _npChartMeta        = null;  // { padLeft, pW, T, ROW_H, nRows, ganttTop }
let _npScheduleForChart = [];    // last enriched schedule (for hit-testing)
let _npWindowStart      = null;  // last chartStartDate
let _npLastDrawArgs     = null;  // cached args for re-draw during drag

// ── Formatting helpers ────────────────────────────────────────────────────────

function npFmt(date) {
  return date.toLocaleTimeString("fr-FR", { hour: "2-digit", minute: "2-digit" });
}

function npDurFmt(min) {
  const h = Math.floor(min / 60);
  const m = Math.round(min % 60);
  return h > 0 ? `${h}h${String(m).padStart(2, "0")}m` : `${m}m`;
}

// ── Sun altitude at an instant (self-contained, no Astronomy dependency) ─────

function _sunAltAt(lat, lon, date) {
  const D2R = Math.PI / 180;
  const jd  = date.getTime() / 86400000 + 2440587.5;
  const n   = jd - 2451545.0;
  const L   = ((280.46  + 0.9856474 * n) % 360 + 360) % 360;
  const g   = ((357.528 + 0.9856003 * n) % 360 + 360) % 360;
  const lam = L + 1.915 * Math.sin(g * D2R) + 0.02 * Math.sin(2 * g * D2R);
  const eps = 23.439 - 0.0000004 * n;
  const lamR = lam * D2R, epsR = eps * D2R;
  const raDeg  = ((Math.atan2(Math.cos(epsR) * Math.sin(lamR), Math.cos(lamR)) * 180 / Math.PI) + 360) % 360;
  const decDeg = Math.asin(Math.sin(epsR) * Math.sin(lamR)) * 180 / Math.PI;
  return Astronomy.computeAltAz(raDeg, decDeg, lat, lon, date).alt;
}

// ── Find the first sun crossing of targetAlt within `hours` of startDate ─────
// wantRising = true → find when sun rises above targetAlt, false → falls below

function npFindSunCrossing(lat, lon, startDate, hours, targetAlt, wantRising) {
  const totalMs = hours * 3600000;
  const STEPS   = 288; // every 5 min
  let prevAlt   = _sunAltAt(lat, lon, startDate);

  for (let i = 1; i <= STEPS; i++) {
    const ms  = (i / STEPS) * totalMs;
    const alt = _sunAltAt(lat, lon, new Date(startDate.getTime() + ms));
    const cross = wantRising
      ? (prevAlt < targetAlt && alt >= targetAlt)
      : (prevAlt > targetAlt && alt <= targetAlt);

    if (cross) {
      let lo = (i - 1) / STEPS * totalMs;
      let hi = ms;
      for (let j = 0; j < 48; j++) {
        const mid = (lo + hi) / 2;
        const a   = _sunAltAt(lat, lon, new Date(startDate.getTime() + mid));
        (wantRising ? a < targetAlt : a > targetAlt) ? (lo = mid) : (hi = mid);
      }
      return new Date(startDate.getTime() + (lo + hi) / 2);
    }
    prevAlt = alt;
  }
  return null;
}

// ── Peak altitude time for a target within a window ──────────────────────────

function npFindPeakTime(raDeg, decDeg, lat, lon, startDate, hours) {
  const steps = Math.round(hours * 6); // every 10 min
  let bestAlt  = -Infinity;
  let bestTime = startDate;
  for (let i = 0; i <= steps; i++) {
    const t   = new Date(startDate.getTime() + i * 600000);
    const alt = Astronomy.computeAltAz(raDeg, decDeg, lat, lon, t).alt;
    if (alt > bestAlt) { bestAlt = alt; bestTime = t; }
  }
  return { time: bestTime, alt: bestAlt };
}

// ── Schedule computation ──────────────────────────────────────────────────────

function _parseWaitTime(hhMM, afterDate) {
  const [h, m] = (hhMM || "00:00").split(":").map(Number);
  const d = new Date(afterDate);
  d.setHours(h, m, 0, 0);
  if (d.getTime() <= afterDate.getTime()) d.setDate(d.getDate() + 1);
  return d;
}

function npComputeSchedule(queue, params, nightStart, nightEnd) {
  let current = new Date(nightStart);
  const limitMs = nightEnd
    ? nightEnd.getTime()
    : nightStart.getTime() + 12 * 3600000;

  const slots = [];

  // ── Init block (open roof, unpark, etc.) ────────────────────────────────
  if (params.initMin > 0) {
    const end = new Date(current.getTime() + params.initMin * 60000);
    slots.push({ type: "init", target: { name: "Init", color: "#6b7280" }, start: new Date(current), end, durationMin: params.initMin, overflows: false });
    current = end;
  }

  // ── Queue items ──────────────────────────────────────────────────────────
  for (const t of queue) {
    if (t.done) continue;

    if (t.itemType === "wait") {
      const waitUntil = _parseWaitTime(t.waitTime, current);
      if (waitUntil.getTime() > current.getTime()) {
        const waitMin = (waitUntil.getTime() - current.getTime()) / 60000;
        slots.push({ type: "wait", target: { name: `⏳ ${t.waitTime}`, color: "#374151" }, start: new Date(current), end: waitUntil, durationMin: waitMin, overflows: false });
        current = waitUntil;
      }
      continue;
    }

    if (!t.raDeg) continue;

    const imagingMin  = (params.filters * params.count * (params.exposure + params.save)) / 60;
    const durationMin = params.af + params.slew + imagingMin;
    const start       = new Date(current);
    const end         = new Date(current.getTime() + durationMin * 60000);
    current           = end;
    slots.push({ type: "target", target: t, start, end, durationMin, overflows: end.getTime() > limitMs });
  }

  // ── End block (park mount, close roof, etc.) ─────────────────────────────
  if (params.endMin > 0) {
    const end = new Date(current.getTime() + params.endMin * 60000);
    slots.push({ type: "end", target: { name: "End", color: "#6b7280" }, start: new Date(current), end, durationMin: params.endMin, overflows: false });
  }

  return slots;
}

// ── Arbiter-based schedule simulation ────────────────────────────────────────
//
// Mirrors the server-side pickBestTarget + runSequence arbiter logic.
// At each decision point we pick the highest-altitude target that is inside
// all safety limits, with an urgency tie-break (closest to leaving the window).
// When nothing is available we advance time in 5-min steps (shown as a "gap"
// slot) and retry.  Max 45 min of total gap before giving up.

/**
 * Pick the best target from candidates at a given date/time.
 * Returns { target, alt, ha, minutesUntilInvalid } or null.
 */
function npPickBest(candidates, lat, lon, date, { minAlt = 20, zenithLimit = 70, meridianGap = 10 } = {}) {
  const scored = candidates
    .filter(t => t.raDeg != null && t.decDeg != null)
    .map(t => {
      const { alt, ha } = Astronomy.computeAltAz(t.raDeg, t.decDeg, lat, lon, date);
      const absHA  = Math.abs(ha);
      const valid  = alt >= minAlt && alt <= zenithLimit && absHA >= meridianGap;

      // Sample forward to find when the target leaves the valid window
      let minutesUntilInvalid = Infinity;
      if (valid) {
        for (let m = 5; m <= 240; m += 5) {
          const future = new Date(date.getTime() + m * 60000);
          const { alt: fAlt, ha: fHA } = Astronomy.computeAltAz(t.raDeg, t.decDeg, lat, lon, future);
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

  // Primary: highest altitude. Tie-break (< 2°): most urgent (sets soonest)
  scored.sort((a, b) => {
    if (Math.abs(b.alt - a.alt) > 2) return b.alt - a.alt;
    return a.minutesUntilInvalid - b.minutesUntilInvalid;
  });

  return scored[0];
}

/**
 * Arbiter simulation: dynamically builds the schedule by choosing the best
 * available target at each decision point, exactly as runSequence does.
 */
function npComputeScheduleArbiter(queue, params, nightStart, nightEnd, lat, lon) {
  const limitMs = nightEnd
    ? nightEnd.getTime()
    : nightStart.getTime() + 12 * 3600000;

  const imagingMin  = (params.filters * params.count * (params.exposure + params.save)) / 60;
  const targetDur   = params.af + params.slew + imagingMin; // minutes per target
  const limits      = { minAlt: params.minAlt, zenithLimit: params.zenithLimit, meridianGap: params.meridianGap };

  const MAX_GAP_MS  = 45 * 60 * 1000;   // give up if nothing available for 45 min
  const POLL_MS     = 5  * 60 * 1000;   // poll interval while waiting

  let current = new Date(nightStart);
  const slots = [];

  // ── Init block ────────────────────────────────────────────────────────────
  if (params.initMin > 0) {
    const end = new Date(current.getTime() + params.initMin * 60000);
    slots.push({ type: "init", target: { name: "Init", color: "#6b7280" }, start: new Date(current), end, durationMin: params.initMin, overflows: false });
    current = end;
  }

  // Handle user-inserted "wait" markers: collect them as sorted checkpoints
  const waitMarkers = queue
    .filter(t => t.itemType === "wait" && t.waitTime)
    .map(t => _parseWaitTime(t.waitTime, nightStart))
    .filter(d => d > current)
    .sort((a, b) => a - b);

  // Remaining observation targets (no wait markers, no done)
  const remaining = new Set(
    queue.filter(t => !t.done && t.itemType !== "wait" && t.raDeg != null)
  );

  let gapMs = 0; // accumulated wait time without picking a target

  while (remaining.size > 0 && current.getTime() < limitMs - params.endMin * 60000) {
    // Respect user wait markers: if we've reached one, insert a wait slot
    if (waitMarkers.length && waitMarkers[0] <= current) {
      waitMarkers.shift(); // already passed — discard
    }
    if (waitMarkers.length && waitMarkers[0] > current && waitMarkers[0] < new Date(current.getTime() + targetDur * 60000)) {
      const wUntil  = waitMarkers.shift();
      const waitMin = (wUntil.getTime() - current.getTime()) / 60000;
      slots.push({ type: "wait", target: { name: `⏳ wait`, color: "#374151" }, start: new Date(current), end: wUntil, durationMin: waitMin, overflows: false });
      current = wUntil;
      continue;
    }

    // Ask the arbiter
    const pick = npPickBest([...remaining], lat, lon, current, limits);

    if (!pick) {
      // No target available — open a gap slot (merge with running gap if possible)
      if (gapMs >= MAX_GAP_MS) break; // give up

      const gapEnd = new Date(current.getTime() + POLL_MS);
      const lastSlot = slots[slots.length - 1];
      if (lastSlot && lastSlot.type === "gap") {
        // Extend the existing gap bar instead of creating a new one
        lastSlot.end = gapEnd;
        lastSlot.durationMin += POLL_MS / 60000;
        lastSlot.endH = undefined; // will be recalculated in enrichment
      } else {
        slots.push({
          type: "gap",
          target: { name: "⏳ waiting…", color: "#374151" },
          start: new Date(current), end: gapEnd,
          durationMin: POLL_MS / 60000,
          overflows: false,
          gapReason: "No target in valid sky zone",
        });
      }
      current  = gapEnd;
      gapMs   += POLL_MS;
      continue;
    }

    gapMs = 0; // reset gap timer when a target is picked
    const t   = pick.target;
    const end = new Date(current.getTime() + targetDur * 60000);
    slots.push({
      type: "target", target: t, start: new Date(current), end,
      durationMin: targetDur,
      overflows: end.getTime() > limitMs,
      // Extra info for tooltip / table
      altAtStart:           pick.alt,
      haAtStart:            pick.ha,
      minutesUntilInvalid:  pick.minutesUntilInvalid,
    });
    remaining.delete(t);
    current = end;
  }

  // ── End block ──────────────────────────────────────────────────────────────
  if (params.endMin > 0 && current.getTime() < limitMs) {
    const end = new Date(current.getTime() + params.endMin * 60000);
    slots.push({ type: "end", target: { name: "End", color: "#6b7280" }, start: new Date(current), end, durationMin: params.endMin, overflows: false });
  }

  return slots;
}

// ── Canvas: multi-target altitude chart + Gantt bars ─────────────────────────

function npDrawChart(canvas, targets, sunCurve, schedule, chartStartDate, totalHours, minAlt, weatherData = null) {
  // Cache for drag interaction
  _npLastDrawArgs     = [canvas, targets, sunCurve, schedule, chartStartDate, totalHours, minAlt, weatherData];
  _npWindowStart      = chartStartDate;
  _npScheduleForChart = schedule;

  const dpr  = window.devicePixelRatio || 1;
  const cssW = canvas.clientWidth  || 800;
  const cssH = canvas.clientHeight || 480;
  canvas.width  = cssW * dpr;
  canvas.height = cssH * dpr;

  const ctx = canvas.getContext("2d");
  ctx.scale(dpr, dpr);
  const W = cssW, H = cssH;

  const nRows     = Math.min(schedule.length, 10);
  const ROW_H     = nRows > 0 ? Math.max(20, Math.min(34, 140 / Math.max(1, nRows))) : 28;
  const GANTT_H   = nRows * ROW_H + 36 + 8;
  const WEATHER_H = weatherData ? 88 : 0; // two-row panel: cloud + rain
  const ganttTopY = H - GANTT_H - WEATHER_H + 6;

  const pad = { left: 46, right: 92, top: 18, bottom: GANTT_H + WEATHER_H + 24 };
  const pW  = W - pad.left - pad.right;
  const pH  = H - pad.top  - pad.bottom;
  const T   = totalHours;

  const altToY = alt => pad.top + pH * (1 - (alt + 5) / 95);
  const hToX   = h   => pad.left + pW * (h / T);

  // Save layout meta for mouse hit-testing
  _npChartMeta = { padLeft: pad.left, pW, T, ROW_H, nRows, ganttTop: ganttTopY };

  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = "#090d18";
  ctx.fillRect(0, 0, W, H);

  // ── Cloud cover strip (top 10 px) ────────────────────────────────────────
  if (weatherData && weatherData.points) {
    const SH = 10;
    for (const p of weatherData.points) {
      const h = (new Date(p.time).getTime() - chartStartDate.getTime()) / 3_600_000;
      if (h < -0.5 || h > T + 0.5) continue;
      const x1 = hToX(Math.max(0, h - 0.5));
      const x2 = hToX(Math.min(T,  h + 0.5));
      const t  = p.cloudCover / 100;
      let r, g, b, a;
      if (p.precipType !== "none") {
        r = 40; g = 80; b = 200; a = 0.55 + t * 0.3;
      } else {
        r = Math.round(34 + t * 100); g = Math.round(197 - t * 130);
        b = Math.round(94 + t * 100); a = 0.2 + t * 0.6;
      }
      ctx.fillStyle = `rgba(${r},${g},${b},${a})`;
      ctx.fillRect(x1, 0, x2 - x1, SH);
    }
    ctx.fillStyle = "#3d5070"; ctx.font = "8px Arial";
    ctx.textAlign = "right"; ctx.textBaseline = "middle";
    ctx.fillText("☁", pad.left - 2, SH / 2);
  }

  // ── Twilight shading ──────────────────────────────────────────────────────
  for (let i = 0; i < sunCurve.length - 1; i++) {
    const mid = (sunCurve[i].alt + sunCurve[i + 1].alt) / 2;
    const x1  = hToX(sunCurve[i].h);
    const x2  = hToX(sunCurve[i + 1].h);
    let color;
    if      (mid > 0)   color = "rgba(140,90,10,0.30)";
    else if (mid > -6)  color = "rgba(100,60,10,0.18)";
    else if (mid > -12) color = "rgba(60,40,80,0.10)";
    else if (mid > -18) color = "rgba(30,20,55,0.06)";
    if (color) {
      ctx.fillStyle = color;
      ctx.fillRect(x1, pad.top, x2 - x1, pH);
    }
  }

  // ── Altitude grid ─────────────────────────────────────────────────────────
  for (const alt of [0, 20, 30, 45, 60, 90]) {
    const y = altToY(alt);
    if (y < pad.top - 2 || y > pad.top + pH + 2) continue;
    ctx.strokeStyle = alt === 0 ? "#374151" : "#1a2535";
    ctx.lineWidth   = alt === 0 ? 1 : 0.6;
    ctx.beginPath(); ctx.moveTo(pad.left, y); ctx.lineTo(W - pad.right, y); ctx.stroke();
    ctx.fillStyle   = "#4b6080";
    ctx.font        = "10px Arial";
    ctx.textAlign   = "right";
    ctx.textBaseline = "middle";
    ctx.fillText(`${alt}°`, pad.left - 4, y);
  }

  // ── Hour grid ─────────────────────────────────────────────────────────────
  for (let h = 0; h <= T; h += 2) {
    const x = hToX(h);
    ctx.strokeStyle = "#1a2535"; ctx.lineWidth = 0.5;
    ctx.beginPath(); ctx.moveTo(x, pad.top); ctx.lineTo(x, pad.top + pH); ctx.stroke();
    const t = new Date(chartStartDate.getTime() + h * 3600000);
    ctx.fillStyle = "#4b6080"; ctx.font = "10px Arial";
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    ctx.fillText(`${String(t.getHours()).padStart(2, "0")}h`, x, pad.top + pH + 4);
  }

  // ── Min altitude dashed line ──────────────────────────────────────────────
  const yMin = altToY(minAlt);
  if (yMin >= pad.top && yMin <= pad.top + pH) {
    ctx.strokeStyle = "#ef444488"; ctx.lineWidth = 1.2;
    ctx.setLineDash([4, 4]);
    ctx.beginPath(); ctx.moveTo(pad.left, yMin); ctx.lineTo(W - pad.right, yMin); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = "#ef4444"; ctx.font = "10px Arial";
    ctx.textAlign = "left"; ctx.textBaseline = "middle";
    ctx.fillText(`${minAlt}° min`, W - pad.right + 4, yMin);
  }

  // ── Sun curve ─────────────────────────────────────────────────────────────
  if (sunCurve.length > 1) {
    ctx.beginPath();
    sunCurve.forEach((p, i) => {
      const x = hToX(p.h), y = altToY(p.alt);
      i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    });
    ctx.strokeStyle = "rgba(251,191,36,0.40)"; ctx.lineWidth = 1.5; ctx.stroke();
  }

  // ── Target curves ─────────────────────────────────────────────────────────
  for (const tgt of targets) {
    const pts = tgt.points;
    if (pts.length < 2) continue;

    // Soft fill above minAlt
    ctx.beginPath();
    ctx.moveTo(hToX(pts[0].h), Math.min(altToY(pts[0].alt), yMin));
    pts.forEach(p => ctx.lineTo(hToX(p.h), Math.min(altToY(p.alt), yMin)));
    ctx.lineTo(hToX(pts[pts.length - 1].h), yMin);
    ctx.closePath();
    ctx.fillStyle = tgt.color + "18"; ctx.fill();

    // Curve line
    ctx.beginPath();
    pts.forEach((p, i) => {
      const x = hToX(p.h), y = altToY(p.alt);
      i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    });
    ctx.strokeStyle = tgt.color; ctx.lineWidth = 1.8; ctx.stroke();

    // Label at peak
    let peak = pts[0];
    pts.forEach(p => { if (p.alt > peak.alt) peak = p; });
    const px = hToX(peak.h);
    const py = altToY(peak.alt) - 4;
    ctx.font = "bold 11px Arial";
    ctx.textAlign = "center"; ctx.textBaseline = "bottom";
    const tw = ctx.measureText(tgt.name).width + 8;
    // Background pill
    ctx.fillStyle = "#090d18cc";
    ctx.beginPath();
    if (ctx.roundRect) {
      ctx.roundRect(px - tw / 2, py - 14, tw, 14, 4);
    } else {
      ctx.rect(px - tw / 2, py - 14, tw, 14);
    }
    ctx.fill();
    ctx.fillStyle = tgt.color;
    ctx.fillText(tgt.name, px, py);
  }

  // ── "Now" vertical line ───────────────────────────────────────────────────
  const nowH = (Date.now() - chartStartDate.getTime()) / 3600000;
  if (nowH >= 0 && nowH <= T) {
    const nowX = hToX(nowH);
    ctx.strokeStyle = "#fbbf24"; ctx.lineWidth = 1.5;
    ctx.setLineDash([4, 3]);
    ctx.beginPath(); ctx.moveTo(nowX, pad.top); ctx.lineTo(nowX, pad.top + pH); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = "#fbbf2490"; ctx.font = "bold 10px Arial";
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    ctx.fillText("now", nowX, pad.top + pH + 4);
  }

  // ── Gantt section ─────────────────────────────────────────────────────────
  const ganttTop = ganttTopY;

  ctx.fillStyle = "#3a5070"; ctx.font = "10px Arial";
  ctx.textAlign = "left"; ctx.textBaseline = "top";
  ctx.fillText("Observation schedule", pad.left, ganttTop);

  // Gantt hour ticks
  for (let h = 0; h <= T; h += 2) {
    const x = hToX(h);
    ctx.strokeStyle = "#1a2535"; ctx.lineWidth = 0.4;
    ctx.beginPath(); ctx.moveTo(x, ganttTop + 14); ctx.lineTo(x, ganttTop + 14 + nRows * ROW_H); ctx.stroke();
  }

  schedule.slice(0, 10).forEach((slot, ri) => {
    const y0 = ganttTop + 14 + ri * ROW_H + 2;
    const bH = ROW_H - 4;
    const x1 = hToX(slot.startH);
    const x2 = hToX(slot.endH);

    let barColor;
    if (slot.type === "init" || slot.type === "end") {
      barColor = "#6b7280";
    } else if (slot.type === "wait" || slot.type === "gap") {
      barColor = "#374151";
    } else {
      barColor = slot.status === "red" ? "#ef4444" : slot.status === "orange" ? "#f59e0b" : (slot.target.color || "#60a5fa");
    }

    // Dashed border for wait / gap items
    if (slot.type === "wait" || slot.type === "gap") ctx.setLineDash([4, 3]);
    ctx.fillStyle   = barColor + "44";
    ctx.strokeStyle = barColor;
    ctx.lineWidth   = slot.type === "init" || slot.type === "end" ? 1.5 : 1;
    ctx.beginPath();
    if (ctx.roundRect) {
      ctx.roundRect(x1, y0, Math.max(3, x2 - x1), bH, 3);
    } else {
      ctx.rect(x1, y0, Math.max(3, x2 - x1), bH);
    }
    ctx.fill(); ctx.stroke();
    if (slot.type === "wait" || slot.type === "gap") ctx.setLineDash([]);

    // Label inside bar
    if (x2 - x1 > 28) {
      const label = slot.target.name.length > 12 ? slot.target.name.slice(0, 11) + "…" : slot.target.name;
      ctx.fillStyle = barColor;
      ctx.font      = `${Math.min(11, bH - 4)}px Arial`;
      ctx.textAlign = "left"; ctx.textBaseline = "middle";
      ctx.fillText(label, x1 + 4, y0 + bH / 2);
    }

    // Right-side label
    ctx.fillStyle = "#6b8090"; ctx.font = "10px Arial";
    ctx.textAlign = "left"; ctx.textBaseline = "middle";
    ctx.fillText(slot.target.name, W - pad.right + 6, y0 + bH / 2);
  });

  // ── Drag overlay ─────────────────────────────────────────────────────────────
  if (_npDrag !== null && _npDrag.slotIdx < schedule.length) {
    const ri   = _npDrag.slotIdx;
    const slot = schedule[ri];
    // Only render overlay for target slots (no `type` property)
    if (!slot.type) {
      const dH        = _npDrag.deltaH;
      const origW     = slot.endH - slot.startH;
      const newStartH = Math.max(0, slot.startH + dH);
      const newEndH   = newStartH + origW;
      const y0        = ganttTop + 14 + ri * ROW_H + 2;
      const bH        = ROW_H - 4;
      const ox1       = hToX(slot.startH);
      const ox2       = hToX(slot.endH);
      const nx1       = hToX(newStartH);
      const nx2       = hToX(newEndH);
      const col       = slot.target?.color || "#60a5fa";

      // Dim original position
      ctx.fillStyle = "rgba(9,13,24,0.72)";
      ctx.fillRect(ox1, y0, Math.max(3, ox2 - ox1), bH);

      // Ghost bar at new position
      ctx.setLineDash([]);
      ctx.fillStyle   = col + "55";
      ctx.strokeStyle = col;
      ctx.lineWidth   = 2;
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(nx1, y0, Math.max(3, nx2 - nx1), bH, 3);
      else ctx.rect(nx1, y0, Math.max(3, nx2 - nx1), bH);
      ctx.fill(); ctx.stroke();

      // Label in ghost bar
      if (nx2 - nx1 > 28) {
        const lbl = slot.target.name.length > 12 ? slot.target.name.slice(0, 11) + "…" : slot.target.name;
        ctx.fillStyle = col; ctx.font = `${Math.min(11, bH - 4)}px Arial`;
        ctx.textAlign = "left"; ctx.textBaseline = "middle";
        ctx.fillText(lbl, nx1 + 4, y0 + bH / 2);
      }

      // Tooltip bubble above ghost bar
      const newStartDate = new Date(chartStartDate.getTime() + newStartH * 3600000);
      const hh      = String(newStartDate.getHours()).padStart(2, "0");
      const mm      = String(newStartDate.getMinutes()).padStart(2, "0");
      const deltaMin = Math.round(dH * 60);
      const sign    = deltaMin >= 0 ? "+" : "";
      const tip     = `${hh}:${mm}  (${sign}${deltaMin} min)`;
      ctx.font = "bold 11px Arial";
      const tw  = ctx.measureText(tip).width + 14;
      const tx  = Math.min(W - pad.right - tw / 2, Math.max(pad.left + tw / 2, (nx1 + nx2) / 2));
      const ty  = Math.max(4, y0 - 24);
      ctx.fillStyle   = "#0f172aee";
      ctx.strokeStyle = col;
      ctx.lineWidth   = 1.2;
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(tx - tw / 2, ty, tw, 18, 5);
      else ctx.rect(tx - tw / 2, ty, tw, 18);
      ctx.fill(); ctx.stroke();
      ctx.fillStyle    = col;
      ctx.textAlign    = "center";
      ctx.textBaseline = "middle";
      ctx.fillText(tip, tx, ty + 9);

      // Vertical dashed line at new start
      ctx.strokeStyle = col + "88"; ctx.lineWidth = 1;
      ctx.setLineDash([3, 3]);
      ctx.beginPath(); ctx.moveTo(nx1, pad.top); ctx.lineTo(nx1, ganttTop + 14 + nRows * ROW_H); ctx.stroke();
      ctx.setLineDash([]);
    }
  }

  // ── Weather panel: two separate rows inside the canvas ──────────────────────
  //   Row 1 (top half):    ☁ CLOUD   – filled bar, 0-100 %
  //   Row 2 (bottom half): 🌧 RAIN    – solid accent only when > 0
  if (WEATHER_H > 0 && weatherData && weatherData.points) {
    const wTop  = H - WEATHER_H;
    const ROW_H = (WEATHER_H - 4) / 2;     // height of each sub-row
    const LBL_W = pad.left;                 // reuse left padding for labels

    // Map hourly data to chart-time offsets
    const wpts = weatherData.points.map(p => ({
      h:     (new Date(p.time).getTime() - chartStartDate.getTime()) / 3_600_000,
      cc:    p.cloudCover / 100,           // 0–1
      mm:    p.precip    ?? 0,             // mm/h
      snow:  p.precipType === "snow",
    })).filter(p => p.h > -2 && p.h < T + 2);

    // Panel background
    ctx.fillStyle = "#07090f";
    ctx.fillRect(0, wTop, W, WEATHER_H);

    // Outer separator
    ctx.strokeStyle = "#1a2535"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(0, wTop + 0.5); ctx.lineTo(W, wTop + 0.5); ctx.stroke();

    // Mid separator between the two rows
    const midY = wTop + ROW_H + 2;
    ctx.strokeStyle = "#12192a"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(pad.left, midY); ctx.lineTo(W - pad.right, midY); ctx.stroke();

    // Row labels (right-aligned into left padding)
    ctx.font = "bold 8px Arial"; ctx.textAlign = "right"; ctx.textBaseline = "middle";
    ctx.fillStyle = "#3a5070";
    ctx.fillText("☁", pad.left - 4, wTop + ROW_H / 2 + 1);
    ctx.fillStyle = "#1e4a7a";
    ctx.fillText("🌧", pad.left - 4, midY + ROW_H / 2 + 1);

    if (wpts.length >= 2) {
      for (let i = 0; i < wpts.length; i++) {
        const h1 = i === 0             ? wpts[0].h - 0.5              : (wpts[i-1].h + wpts[i].h) / 2;
        const h2 = i === wpts.length-1 ? wpts[wpts.length-1].h + 0.5  : (wpts[i].h + wpts[i+1].h) / 2;
        const x1 = hToX(Math.max(0, h1)) + 1;
        const x2 = hToX(Math.min(T, h2)) - 1;
        if (x2 - x1 < 1) continue;
        const bW = Math.max(1, x2 - x1);

        // ── Cloud row ──
        const cc   = wpts[i].cc;
        const cTop = wTop + 2;
        const cH   = ROW_H - 4;
        // background always dark
        ctx.fillStyle = "#0b0f1a";
        ctx.fillRect(x1, cTop, bW, cH);
        // filled bar proportional to cloud cover
        if (cc > 0.02) {
          const shade = Math.round(30 + cc * 80);   // 30..110
          ctx.fillStyle = `rgb(${shade},${Math.round(shade*1.3)},${Math.round(shade*2)})`;
          ctx.fillRect(x1, cTop + cH * (1 - cc), bW, cH * cc);
        }
        // percentage label
        ctx.font = "8px Arial"; ctx.textAlign = "center"; ctx.textBaseline = "middle";
        ctx.fillStyle = cc > 0.5 ? "#e2e8f0" : "#475569";
        ctx.fillText(`${Math.round(cc * 100)}%`, (x1 + x2) / 2, cTop + cH / 2);

        // ── Rain row ──
        const mm   = wpts[i].mm;
        const rTop = midY + 2;
        const rH   = ROW_H - 4;
        if (mm > 0) {
          const intensity = Math.min(1, mm / 5);     // saturate at 5 mm/h
          const bg = wpts[i].snow
            ? `rgba(120,180,255,${0.4 + intensity * 0.5})`
            : `rgba(30,100,220,${0.4 + intensity * 0.5})`;
          ctx.fillStyle = bg;
          ctx.fillRect(x1, rTop, bW, rH);
          ctx.font = "8px Arial"; ctx.textAlign = "center"; ctx.textBaseline = "middle";
          ctx.fillStyle = "#e2e8f0";
          ctx.fillText(wpts[i].snow ? `❄${mm.toFixed(1)}` : `${mm.toFixed(1)}`, (x1 + x2) / 2, rTop + rH / 2);
        } else {
          ctx.fillStyle = "#080c14";
          ctx.fillRect(x1, rTop, bW, rH);
        }
      }
    } else {
      ctx.fillStyle = "#2a3a50"; ctx.font = "10px Arial";
      ctx.textAlign = "center"; ctx.textBaseline = "middle";
      ctx.fillText("No weather data", pad.left + pW / 2, wTop + WEATHER_H / 2);
    }
  }
}

// ── Main: refresh Night Plan ──────────────────────────────────────────────────

async function refreshNightPlan() {
  // observerLocation comes from app.js (already in scope)
  const lat = (typeof observerLocation !== "undefined") ? observerLocation.lat : null;
  const lon = (typeof observerLocation !== "undefined") ? observerLocation.lon : null;

  // Refresh weather in parallel (non-blocking)
  fetchWeather();

  const bar = document.getElementById("nightplan-info-bar");

  if (!lat || !lon) {
    bar.textContent = "⚠ Set observer location first (Target tab → gear icon ⚙).";
    bar.className   = "nightplan-info-bar nightplan-warn";
    return;
  }

  bar.textContent = "Computing…";
  bar.className   = "nightplan-info-bar";

  // Fetch current queue
  let queue = [];
  try {
    const state = await fetch("/api/sequence/state").then(r => r.json());
    queue = state.queue || [];
  } catch (e) {
    console.warn("Night Plan: could not fetch queue", e);
  }
  _npQueue = queue;

  // Read overhead params
  const params = {
    af:           parseFloat(document.getElementById("np-af").value)           || 3,
    slew:         parseFloat(document.getElementById("np-slew").value)         || 3,
    save:         parseFloat(document.getElementById("np-save").value)         || 5,
    minAlt:       parseFloat(document.getElementById("set-minalt")?.value || document.getElementById("np-minalt")?.value) || 20,
    initMin:      parseFloat(document.getElementById("np-init-min").value)     || 0,
    endMin:       parseFloat(document.getElementById("np-end-min").value)      || 0,
    meridianGap:  parseFloat(document.getElementById("set-meridian-gap")?.value) || 10,
    zenithLimit:  parseFloat(document.getElementById("set-zenith-limit")?.value) || 70,
    // read from ToDo tab
    exposure: parseFloat(document.getElementById("seq-duration").value) || 120,
    count:    parseInt(document.getElementById("seq-count").value)    || 10,
    filters:  (document.getElementById("seq-filters").value || "G,BP,RP").split(",").filter(Boolean).length,
  };

  // Chart window: start at 18:00 local (or yesterday if morning), show 14 h
  const now = new Date();
  const windowStart = new Date(now);
  windowStart.setHours(18, 0, 0, 0);
  if (now.getHours() < 6) windowStart.setDate(windowStart.getDate() - 1);
  // If it's well past 18:00 but still dark, anchor window 2h before now
  if (now > windowStart && (now.getTime() - windowStart.getTime()) > 2 * 3600000) {
    windowStart.setTime(now.getTime() - 2 * 3600000);
  }
  const totalHours = 14;

  // Sun events
  const astroStart = npFindSunCrossing(lat, lon, windowStart, totalHours, -18, false);
  const astroEnd   = npFindSunCrossing(lat, lon, windowStart, totalHours, -18, true);
  const civilDusk  = npFindSunCrossing(lat, lon, windowStart, totalHours, -6,  false);
  const nightStart = astroStart || civilDusk || windowStart;
  const nightEnd   = astroEnd;

  // Night info bar
  const nowMs = Date.now();
  if (astroStart && astroEnd) {
    const durMin = (astroEnd - astroStart) / 60000;
    if (nowMs < astroStart.getTime()) {
      const inMin = Math.round((astroStart.getTime() - nowMs) / 60000);
      const inH   = Math.floor(inMin / 60);
      const inM   = inMin % 60;
      bar.innerHTML = `🌆 Astronomical night in ${inH > 0 ? inH + "h " : ""}${inM}m &nbsp;·&nbsp; Duration: ${npDurFmt(durMin)} &nbsp;·&nbsp; ${npFmt(astroStart)} → ${npFmt(astroEnd)}`;
      bar.className = "nightplan-info-bar nightplan-day";
    } else if (nowMs < astroEnd.getTime()) {
      bar.innerHTML = `🌙 Astronomical night active &nbsp;·&nbsp; Ends: ${npFmt(astroEnd)} &nbsp;·&nbsp; Remaining: ${npDurFmt((astroEnd.getTime() - nowMs) / 60000)}`;
      bar.className = "nightplan-info-bar nightplan-night";
    } else {
      bar.innerHTML = `☀️ Night over &nbsp;·&nbsp; Was ${npFmt(astroStart)} → ${npFmt(astroEnd)}`;
      bar.className = "nightplan-info-bar nightplan-day";
    }
  } else {
    bar.textContent = "⚠ Could not determine astronomical night (check observer location).";
    bar.className   = "nightplan-info-bar nightplan-warn";
  }

  // Target curves
  const targetData = queue.filter(t => t.raDeg && !t.done).map((t, i) => ({
    ...t,
    color:  NP_COLORS[i % NP_COLORS.length],
    points: Astronomy.computeCurve(t.raDeg, t.decDeg, lat, lon, windowStart, totalHours, 10),
  }));

  // Sun curve
  const sunCurve = Astronomy.computeSunCurve(lat, lon, windowStart, totalHours, 10);

  // Schedule start: beginning of astronomical night, or now if already in night
  let schedStart = new Date(nightStart);
  if (nightEnd && nowMs > schedStart.getTime() && nowMs < nightEnd.getTime()) {
    schedStart = new Date(nowMs);
  }

  // Use arbiter simulation when observer location is known, otherwise fall back
  // to the simple fixed-order schedule.
  const canArbiter = lat != null && lon != null;
  const rawSchedule = canArbiter
    ? npComputeScheduleArbiter(queue, params, schedStart, nightEnd, lat, lon)
    : npComputeSchedule(queue, params, schedStart, nightEnd);
  _npSchedule = rawSchedule;

  // Enrich schedule with startH/endH and altitude status
  const scheduleForChart = rawSchedule.map(slot => {
    const startH = (slot.start.getTime() - windowStart.getTime()) / 3600000;
    const endH   = (slot.end.getTime()   - windowStart.getTime()) / 3600000;

    if (slot.type !== "target") {
      return { ...slot, startH, endH, minAlt: 0, maxAlt: 0, status: slot.type };
    }

    const steps  = Math.max(2, Math.round(slot.durationMin / 5));
    let minA = Infinity, maxA = -Infinity;
    let minAbsHA = Infinity;
    let zenithViolation = false, meridianViolation = false;
    for (let s = 0; s <= steps; s++) {
      const t   = new Date(slot.start.getTime() + s * 5 * 60000);
      if (t > slot.end) break;
      const { alt, ha } = Astronomy.computeAltAz(slot.target.raDeg, slot.target.decDeg, lat, lon, t);
      if (alt < minA) minA = alt;
      if (alt > maxA) maxA = alt;
      if (Math.abs(ha) < minAbsHA) minAbsHA = Math.abs(ha);
      if (alt  > params.zenithLimit)  zenithViolation   = true;
      if (Math.abs(ha) < params.meridianGap) meridianViolation = true;
    }
    if (!isFinite(minA)) { minA = 0; maxA = 0; }
    // Status: red = definitely bad, orange = marginal, green = good
    let status;
    if (maxA < params.minAlt || (zenithViolation && meridianViolation))
      status = "red";
    else if (minA < params.minAlt || zenithViolation || meridianViolation)
      status = "orange";
    else
      status = "green";
    const tgt = targetData.find(t => t.name === slot.target.name) || { ...slot.target, color: NP_COLORS[0] };
    return {
      target: tgt, startH, endH, status, minAlt: minA, maxAlt: maxA,
      overflows: slot.overflows, zenithViolation, meridianViolation, minAbsHA,
      // Carry through arbiter metadata if present
      altAtStart: slot.altAtStart, haAtStart: slot.haAtStart,
      minutesUntilInvalid: slot.minutesUntilInvalid,
    };
  });

  // Draw chart
  const canvas = document.getElementById("nightplan-canvas");
  npDrawChart(canvas, targetData, sunCurve, scheduleForChart, windowStart, totalHours, params.minAlt, _wtData);

  // Show/hide arbiter badge
  const arbiterBadge = document.getElementById("np-arbiter-badge");
  if (arbiterBadge) {
    arbiterBadge.style.display = canArbiter ? "inline-flex" : "none";
    arbiterBadge.title = canArbiter
      ? `Order optimised by arbiter (alt ≥${params.minAlt}°, alt ≤${params.zenithLimit}°, |HA| ≥${params.meridianGap}°)`
      : "Set observer location to enable arbiter ordering";
  }

  // Render table
  const tbody   = document.getElementById("nightplan-tbody");
  const emptyEl = document.getElementById("nightplan-empty");
  const tableEl = document.getElementById("nightplan-table");

  if (rawSchedule.length === 0) {
    emptyEl.style.display = "block";
    tableEl.style.display = "none";
    return;
  }
  emptyEl.style.display = "none";
  tableEl.style.display = "";
  tbody.innerHTML = "";

  rawSchedule.forEach((slot, i) => {
    const sc  = scheduleForChart[i];
    const c   = slot.target.color || "#60a5fa";
    const tr  = document.createElement("tr");

    if (slot.type === "init" || slot.type === "end") {
      tr.style.opacity = "0.6";
      tr.innerHTML = `
        <td>${i + 1}</td>
        <td><span style="color:#9ca3af;">🔧 ${slot.target.name}</span></td>
        <td>${npFmt(slot.start)}</td>
        <td>${npFmt(slot.end)}</td>
        <td>${npDurFmt(slot.durationMin)}</td>
        <td>—</td><td>—</td>
      `;
    } else if (slot.type === "wait") {
      tr.style.opacity = "0.7";
      const waiting = slot.durationMin > 0
        ? `${npDurFmt(slot.durationMin)} idle`
        : "already past";
      tr.innerHTML = `
        <td>${i + 1}</td>
        <td><span style="color:#60a5fa;font-style:italic;">${slot.target.name}</span></td>
        <td>${npFmt(slot.start)}</td>
        <td>${npFmt(slot.end)}</td>
        <td>${waiting}</td>
        <td>—</td><td>—</td>
      `;
    } else if (slot.type === "gap") {
      tr.style.opacity = "0.55";
      tr.innerHTML = `
        <td>${i + 1}</td>
        <td><span style="color:#6b7280;font-style:italic;">⏳ Waiting</span></td>
        <td>${npFmt(slot.start)}</td>
        <td>${npFmt(slot.end)}</td>
        <td>${npDurFmt(slot.durationMin)}</td>
        <td>—</td>
        <td style="color:#6b7280;">${slot.gapReason || "No valid target"}</td>
      `;
    } else {
      const icon = sc.status === "green" ? "✅" : sc.status === "orange" ? "⚠️" : "❌";
      const altTxt = `${sc.minAlt.toFixed(0)}°–${sc.maxAlt.toFixed(0)}°`;
      const warnings = [];
      if (sc.zenithViolation)   warnings.push(`alt > ${params.zenithLimit}° (zenith limit)`);
      if (sc.meridianViolation) warnings.push(`|HA| < ${params.meridianGap}° (meridian gap)`);
      if (sc.minAlt < params.minAlt) warnings.push(`alt < ${params.minAlt}° (too low)`);

      // Arbiter extra info: altitude and time-to-leave at start of slot
      let arbiterTip = "";
      if (sc.altAtStart != null) {
        const haSide = sc.haAtStart >= 0 ? "W" : "E";
        const timeLeft = sc.minutesUntilInvalid === Infinity ? ">4h" : `${sc.minutesUntilInvalid}min`;
        arbiterTip = ` · picked at alt ${sc.altAtStart.toFixed(0)}° HA${Math.abs(sc.haAtStart).toFixed(0)}°${haSide} valid for ${timeLeft}`;
      }

      const statusTxt = sc.status === "green"
        ? `OK (alt ${sc.minAlt.toFixed(0)}–${sc.maxAlt.toFixed(0)}°, |HA|min ${isFinite(sc.minAbsHA) ? sc.minAbsHA.toFixed(0) : "?"}°)${arbiterTip}`
        : warnings.join("; ") + arbiterTip;
      tr.innerHTML = `
        <td>${i + 1}</td>
        <td><span style="color:${c};font-weight:600;">${slot.target.name}</span></td>
        <td>${npFmt(slot.start)}</td>
        <td>${npFmt(slot.end)}${slot.overflows ? " <span title='Ends after night' style='color:#f59e0b;'>⚠</span>" : ""}</td>
        <td>${npDurFmt(slot.durationMin)}</td>
        <td>${altTxt}</td>
        <td>${icon} ${statusTxt}</td>
      `;
    }
    tbody.appendChild(tr);
  });
}

// ── Optimize Order ────────────────────────────────────────────────────────────

async function npOptimizeOrder() {
  const lat = (typeof observerLocation !== "undefined") ? observerLocation.lat : null;
  const lon = (typeof observerLocation !== "undefined") ? observerLocation.lon : null;
  if (!lat || !lon) { alert("Set observer location first."); return; }

  const btn = document.getElementById("np-optimize-btn");
  btn.disabled    = true;
  btn.textContent = "Computing…";

  try {
    const now         = new Date();
    const windowStart = new Date(now);
    windowStart.setHours(18, 0, 0, 0);
    if (now.getHours() < 6) windowStart.setDate(windowStart.getDate() - 1);

    const nightStart = npFindSunCrossing(lat, lon, windowStart, 14, -18, false) || windowStart;

    const targets = _npQueue.filter(t => t.raDeg && !t.done);
    const withPeak = targets.map(t => {
      const { time } = npFindPeakTime(t.raDeg, t.decDeg, lat, lon, nightStart, 14);
      return { ...t, peakMs: time ? time.getTime() : Infinity };
    });
    withPeak.sort((a, b) => a.peakMs - b.peakMs);

    const resp = await fetch("/api/sequence/reorder", {
      method:  "POST",
      headers: { "Content-Type": "application/json" },
      body:    JSON.stringify({ order: withPeak.map(t => t.name) }),
    });
    const data = await resp.json();
    if (!data.success) throw new Error(data.error || "reorder failed");

    await refreshNightPlan();
    // Also refresh the ToDo tab queue display if visible
    if (typeof refreshSeqState === "function") refreshSeqState();
  } catch (e) {
    alert("Failed to optimize order: " + e.message);
  } finally {
    btn.disabled    = false;
    btn.textContent = "⚡ Optimize Order";
  }
}

// ── 7Timer weather ────────────────────────────────────────────────────────────

const _7T_CC_PCT   = [0, 3, 12, 25, 37, 50, 62, 75, 87, 97];
const _7T_SEE_LBL  = ["", "<0.5\"", "0.6\"", "0.9\"", "1.1\"", "1.4\"", "1.8\"", "2.3\"", ">2.5\""];
const _7T_TRP_LBL  = ["", "Exc", "VGood", "Good", "Good", "Avg", "Poor", "Poor", "Bad"];
const _7T_WIND_MPS = [0, 0.2, 1.8, 5.7, 9.4, 14.0, 20.8, 28.5, 36];
const _7T_WIND_ARR = { N:"↓", NNE:"↙", NE:"↙", ENE:"←", E:"←", ESE:"↖", SE:"↖", SSE:"↑", S:"↑", SSW:"↑", SW:"↗", WSW:"↗", W:"→", WNW:"↘", NW:"↘", NNW:"↓" };
const _7T_PREC_ICO = { none:"", rain:"🌧", snow:"❄", frzr:"🌨", icep:"🌨", tsra:"⛈" };

function _wtBg(v, maxBad = 9) {
  if (!v || v === -9999) return "#1a2535";
  const f = v / maxBad;
  if (f <= 0.22) return "#14532d";
  if (f <= 0.44) return "#166534";
  if (f <= 0.60) return "#713f12";
  if (f <= 0.78) return "#7f1d1d";
  return "#450a0a";
}

function _wtWindBg(s) {
  if (!s || s === -9999) return "#1a2535";
  if (s <= 2) return "#14532d";
  if (s <= 3) return "#166534";
  if (s <= 4) return "#713f12";
  return "#7f1d1d";
}

async function fetchWeather() {
  const lat = (typeof observerLocation !== "undefined") ? observerLocation.lat : null;
  const lon = (typeof observerLocation !== "undefined") ? observerLocation.lon : null;
  const status    = document.getElementById("np-weather-status");
  const container = document.getElementById("np-weather-container");

  if (!lat || !lon) {
    if (status) status.textContent = "Set observer location to fetch forecast.";
    return;
  }

  if (status) { status.textContent = "Loading…"; status.style.display = ""; }

  try {
    const resp = await fetch(`/api/weather/openmeteo?lat=${lat}&lon=${lon}`);
    if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
    const data = await resp.json();
    if (data.error) throw new Error(data.error);
    _wtData = data;
    if (status) status.style.display = "none";
    _renderWeather(container, data);
    // Redraw canvas with weather data now that it's available
    if (_npLastDrawArgs) {
      const args = [..._npLastDrawArgs];
      args[7] = data;
      npDrawChart(...args);
    }
  } catch (e) {
    if (status) { status.textContent = `⚠ Weather unavailable: ${e.message}`; status.style.display = ""; }
    const old = container.querySelector(".weather-grid");
    if (old) old.remove();
  }
}

function _omWindArrow(deg) {
  return ["↓","↙","←","↖","↑","↗","→","↘"][Math.round(deg / 45) % 8];
}
function _omCloudBg(pct) {
  if (pct < 10)  return "#14532d";
  if (pct < 30)  return "#166534";
  if (pct < 55)  return "#713f12";
  if (pct < 75)  return "#7f1d1d";
  return "#450a0a";
}
function _omHumBg(pct)  { return pct > 90 ? "#7f1d1d" : pct > 75 ? "#713f12" : "#1a2535"; }
function _omWindBg(mps) { return mps > 10 ? "#7f1d1d" : mps > 6 ? "#713f12" : mps > 3 ? "#166534" : "#1a2535"; }

function _renderWeather(container, data) {
  // Weather is rendered entirely inside the canvas panel — remove any old table
  const old = container.querySelector(".weather-grid");
  if (old) old.remove();
}

// ── Init (called from app.js on tab switch) ───────────────────────────────────

function initNightPlan() {
  if (_npInitialized) return;
  _npInitialized = true;

  document.getElementById("np-optimize-btn")?.addEventListener("click", npOptimizeOrder);
  document.getElementById("np-refresh-btn").addEventListener("click", refreshNightPlan);

  const simBtn = document.getElementById("np-simulate-btn");
  if (simBtn) {
    simBtn.addEventListener("click", async () => {
      const orig = simBtn.textContent;
      simBtn.disabled    = true;
      simBtn.textContent = "⏳ Simulating…";
      try {
        await refreshNightPlan();
      } finally {
        simBtn.disabled    = false;
        simBtn.textContent = orig;
      }
    });
  }
  document.getElementById("np-weather-btn").addEventListener("click", fetchWeather);

  let _npResizeTimer;
  window.addEventListener("resize", () => {
    clearTimeout(_npResizeTimer);
    _npResizeTimer = setTimeout(() => {
      const panel = document.getElementById("panel-todo");
      if (panel && panel.classList.contains("active")) refreshNightPlan();
    }, 200);
  });

  // Auto-refresh every 60 s when the plan tab is active
  setInterval(() => {
    const panel = document.getElementById("panel-todo");
    if (panel && panel.classList.contains("active")) refreshNightPlan();
  }, 60_000);

  // ── Gantt drag-to-reschedule ──────────────────────────────────────────────
  const canvas = document.getElementById("nightplan-canvas");

  function _ganttHitSlot(x, y) {
    if (!_npChartMeta || !_npScheduleForChart.length) return -1;
    const { padLeft, pW, T, ROW_H, nRows, ganttTop } = _npChartMeta;
    const hToX = h => padLeft + pW * (h / T);
    for (let ri = 0; ri < Math.min(_npScheduleForChart.length, 10); ri++) {
      const slot = _npScheduleForChart[ri];
      if (slot.type) continue; // skip init/end/wait; only targets have no type
      const y0 = ganttTop + 14 + ri * ROW_H + 2;
      const bH = ROW_H - 4;
      const x1 = hToX(slot.startH);
      const x2 = hToX(slot.endH);
      if (y >= y0 && y <= y0 + bH && x >= x1 && x <= x2) return ri;
    }
    return -1;
  }

  canvas.addEventListener("mousedown", e => {
    const rect = canvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    const ri = _ganttHitSlot(x, y);
    if (ri >= 0) {
      _npDrag = { slotIdx: ri, mouseStartX: x, deltaH: 0 };
      canvas.style.cursor = "grabbing";
      e.preventDefault();
    }
  });

  canvas.addEventListener("mousemove", e => {
    const rect = canvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;

    if (_npDrag) {
      const { pW, T } = _npChartMeta;
      let dH = (x - _npDrag.mouseStartX) / pW * T;
      // Snap to 5-minute increments
      dH = Math.round(dH * 12) / 12;
      // Clamp: can't go before chart start
      const slot = _npScheduleForChart[_npDrag.slotIdx];
      if (slot && slot.startH + dH < 0) dH = -slot.startH;
      _npDrag.deltaH = dH;
      if (_npLastDrawArgs) npDrawChart(..._npLastDrawArgs);
      return;
    }

    // Hover cursor: show grab over draggable bars
    const over = _ganttHitSlot(x, y) >= 0;
    canvas.style.cursor = over ? "grab" : "default";
  });

  canvas.addEventListener("mouseup", async e => {
    if (!_npDrag) return;
    const drag = { ..._npDrag };
    _npDrag = null;
    canvas.style.cursor = "default";

    // Ignore tiny drags (< 1 min)
    if (Math.abs(drag.deltaH * 60) < 1) {
      if (_npLastDrawArgs) npDrawChart(..._npLastDrawArgs);
      return;
    }

    const slot = _npScheduleForChart[drag.slotIdx];
    if (!slot || slot.type) { if (_npLastDrawArgs) npDrawChart(..._npLastDrawArgs); return; }

    const newStartH    = Math.max(0, slot.startH + drag.deltaH);
    const newStartDate = new Date(_npWindowStart.getTime() + newStartH * 3600000);
    const hh       = String(newStartDate.getHours()).padStart(2, "0");
    const mm       = String(newStartDate.getMinutes()).padStart(2, "0");
    const waitTime = `${hh}:${mm}`;

    // Find this target's index in the queue
    const targetName = slot.target?.name;
    const queueIdx   = _npQueue.findIndex(t => t.name === targetName && !t.done && t.itemType !== "wait");
    if (queueIdx < 0) { await refreshNightPlan(); return; }

    if (queueIdx > 0 && _npQueue[queueIdx - 1].itemType === "wait") {
      // Update existing wait item immediately before this target
      await fetch(`/api/sequence/queue/wait/${queueIdx - 1}`, {
        method: "PATCH",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ waitTime }),
      });
    } else {
      // Insert a new wait item before this target
      await fetch("/api/sequence/queue/wait-at", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ waitTime, index: queueIdx }),
      });
    }

    await refreshNightPlan();
    if (typeof refreshSeqState === "function") refreshSeqState();
  });

  canvas.addEventListener("mouseleave", () => {
    if (_npDrag) {
      _npDrag = null;
      if (_npLastDrawArgs) npDrawChart(..._npLastDrawArgs);
    }
    canvas.style.cursor = "default";
  });
}
