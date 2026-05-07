/**
 * Browser-side spherical astronomy utilities.
 * Computes alt/az from RA/Dec + observer location, and draws sky charts on Canvas.
 */
(function (global) {
  const D2R = Math.PI / 180;
  const R2D = 180 / Math.PI;

  // ── Math helpers ────────────────────────────────────────────────────────────

  function toJD(date) {
    return date.getTime() / 86400000 + 2440587.5;
  }

  /** Greenwich Mean Sidereal Time in degrees */
  function gmstDeg(jd) {
    const T = (jd - 2451545.0) / 36525.0;
    const g =
      280.46061837 +
      360.98564736629 * (jd - 2451545.0) +
      0.000387933 * T * T -
      (T * T * T) / 38710000;
    return ((g % 360) + 360) % 360;
  }

  /** Altitude and Azimuth (degrees) for a target at RA/Dec from an observer */
  function computeAltAz(raDeg, decDeg, latDeg, lonDeg, date) {
    const jd = toJD(date);
    const lst = (gmstDeg(jd) + lonDeg + 360) % 360;
    const haRaw    = (lst - raDeg + 360) % 360;          // 0–360°
    const haSigned = haRaw > 180 ? haRaw - 360 : haRaw; // −180…+180° (+ = West)
    const haR = haRaw * D2R;
    const decR = decDeg * D2R;
    const latR = latDeg * D2R;

    const sinAlt =
      Math.sin(decR) * Math.sin(latR) +
      Math.cos(decR) * Math.cos(latR) * Math.cos(haR);
    const altR = Math.asin(Math.max(-1, Math.min(1, sinAlt)));

    const cosAz =
      (Math.sin(decR) - Math.sin(altR) * Math.sin(latR)) /
      (Math.cos(altR) * Math.cos(latR));
    let azR = Math.acos(Math.max(-1, Math.min(1, cosAz)));
    if (Math.sin(haR) > 0) azR = 2 * Math.PI - azR;

    return { alt: altR * R2D, az: azR * R2D, ha: haSigned };
  }

  /** Approximate Sun RA/Dec (degrees) — good to ~0.01° */
  function sunRaDec(jd) {
    const n = jd - 2451545.0;
    const L = ((280.46 + 0.9856474 * n) % 360 + 360) % 360;
    const g = ((357.528 + 0.9856003 * n) % 360 + 360) % 360;
    const lam = L + 1.915 * Math.sin(g * D2R) + 0.02 * Math.sin(2 * g * D2R);
    const eps = 23.439 - 0.0000004 * n;
    const lamR = lam * D2R;
    const epsR = eps * D2R;
    return {
      raDeg:
        ((Math.atan2(Math.cos(epsR) * Math.sin(lamR), Math.cos(lamR)) * R2D) +
          360) %
        360,
      decDeg: Math.asin(Math.sin(epsR) * Math.sin(lamR)) * R2D,
    };
  }

  /** Compute alt/az curve for a target over `hours` hours, every `stepMin` minutes */
  function computeCurve(raDeg, decDeg, latDeg, lonDeg, startDate, hours, stepMin) {
    const steps = Math.round((hours * 60) / stepMin);
    return Array.from({ length: steps + 1 }, (_, i) => {
      const t = new Date(startDate.getTime() + i * stepMin * 60000);
      const { alt, az } = computeAltAz(raDeg, decDeg, latDeg, lonDeg, t);
      return { h: i * stepMin / 60, alt, az };
    });
  }

  /** Compute sun altitude curve for twilight shading */
  function computeSunCurve(latDeg, lonDeg, startDate, hours, stepMin) {
    const steps = Math.round((hours * 60) / stepMin);
    return Array.from({ length: steps + 1 }, (_, i) => {
      const t = new Date(startDate.getTime() + i * stepMin * 60000);
      const sun = sunRaDec(toJD(t));
      const { alt } = computeAltAz(sun.raDeg, sun.decDeg, latDeg, lonDeg, t);
      return { h: i * stepMin / 60, alt };
    });
  }

  // ── Canvas drawing ───────────────────────────────────────────────────────────

  /**
   * Draw an azimuth-altitude dome chart.
   * N at top, E at LEFT (standard sky atlas orientation, looking up).
   */
  function drawSkyDome(canvas, altDeg, azDeg) {
    const ctx = canvas.getContext("2d");
    const W = canvas.width;
    const H = canvas.height;
    const cx = W / 2;
    const cy = H / 2;
    const R = Math.min(W, H) / 2 - 24;

    ctx.clearRect(0, 0, W, H);

    // Sky gradient background
    const skyGrad = ctx.createRadialGradient(cx, cy, 0, cx, cy, R);
    skyGrad.addColorStop(0, "#0b1030");
    skyGrad.addColorStop(1, "#1a2248");
    ctx.beginPath();
    ctx.arc(cx, cy, R, 0, 2 * Math.PI);
    ctx.fillStyle = skyGrad;
    ctx.fill();
    ctx.strokeStyle = "#3a5070";
    ctx.lineWidth = 1.5;
    ctx.stroke();

    // Altitude circles at 30° and 60°
    for (const a of [30, 60]) {
      const r = R * ((90 - a) / 90);
      ctx.beginPath();
      ctx.arc(cx, cy, r, 0, 2 * Math.PI);
      ctx.strokeStyle = "#1e3050";
      ctx.lineWidth = 0.8;
      ctx.stroke();
      ctx.fillStyle = "#2a4065";
      ctx.font = "9px Arial";
      ctx.textAlign = "left";
      ctx.textBaseline = "top";
      ctx.fillText(`${a}°`, cx + 4, cy - r + 2);
    }

    // Cardinal dashed lines
    ctx.setLineDash([3, 5]);
    ctx.strokeStyle = "#1a3050";
    ctx.lineWidth = 0.8;
    ctx.beginPath();
    ctx.moveTo(cx, cy - R);
    ctx.lineTo(cx, cy + R);
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(cx - R, cy);
    ctx.lineTo(cx + R, cy);
    ctx.stroke();
    ctx.setLineDash([]);

    // Cardinal labels (N up, E left = sky-atlas convention)
    const cardinals = [
      ["N", 0],
      ["E", 90],
      ["S", 180],
      ["W", 270],
    ];
    ctx.font = "bold 11px Arial";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    for (const [label, az] of cardinals) {
      const azR = az * D2R;
      // Mirror E/W by using negative sin: x = cx - sin(az) * offset
      const x = cx - Math.sin(azR) * (R + 14);
      const y = cy - Math.cos(azR) * (R + 14);
      ctx.fillStyle = label === "N" ? "#8ab4d4" : "#5a7a9a";
      ctx.fillText(label, x, y);
    }

    // Target marker
    if (altDeg !== null && azDeg !== null) {
      const clampedAlt = Math.max(-20, Math.min(90, altDeg));
      const r = R * ((90 - clampedAlt) / 90);
      const azR = azDeg * D2R;
      const tx = cx - Math.sin(azR) * r;
      const ty = cy - Math.cos(azR) * r;

      const color =
        altDeg >= 30 ? "#4ade80" : altDeg >= 0 ? "#fbbf24" : "#f87171";

      // Glow
      if (altDeg > -10) {
        const glow = ctx.createRadialGradient(tx, ty, 0, tx, ty, 18);
        glow.addColorStop(0, color + "44");
        glow.addColorStop(1, color + "00");
        ctx.beginPath();
        ctx.arc(tx, ty, 18, 0, 2 * Math.PI);
        ctx.fillStyle = glow;
        ctx.fill();
      }

      // Crosshair
      ctx.strokeStyle = color;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(tx - 10, ty);
      ctx.lineTo(tx + 10, ty);
      ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(tx, ty - 10);
      ctx.lineTo(tx, ty + 10);
      ctx.stroke();

      // Center dot
      ctx.beginPath();
      ctx.arc(tx, ty, 4, 0, 2 * Math.PI);
      ctx.fillStyle = color;
      ctx.fill();

      // Alt / Az readout
      ctx.fillStyle = color;
      ctx.font = "11px monospace";
      ctx.textAlign = "center";
      ctx.textBaseline = "bottom";
      ctx.fillText(
        `Alt ${altDeg >= 0 ? "+" : ""}${altDeg.toFixed(1)}°  Az ${azDeg.toFixed(1)}°`,
        cx,
        H - 5
      );
    }
  }

  /**
   * Draw altitude vs time chart with twilight shading.
   * targetPoints and sunPoints are arrays of { h, alt } where h = hours from now.
   */
  function drawAltChart(canvas, targetPoints, sunPoints) {
    const ctx = canvas.getContext("2d");
    const W = canvas.width;
    const H = canvas.height;
    const pad = { left: 36, right: 10, top: 14, bottom: 28 };
    const pW = W - pad.left - pad.right;
    const pH = H - pad.top - pad.bottom;
    const T =
      targetPoints.length > 0
        ? targetPoints[targetPoints.length - 1].h
        : 12;

    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = "#090d18";
    ctx.fillRect(0, 0, W, H);

    // -20° to 90° range mapped to canvas height
    const altToY = (alt) => pad.top + pH * (1 - (alt + 20) / 110);
    const hToX = (h) => pad.left + pW * (h / T);

    // Twilight shading: bright overlay for day/twilight, dark background for night
    if (sunPoints.length > 1) {
      for (let i = 0; i < sunPoints.length - 1; i++) {
        const mid = (sunPoints[i].alt + sunPoints[i + 1].alt) / 2;
        const x1 = hToX(sunPoints[i].h);
        const x2 = hToX(sunPoints[i + 1].h);
        let color;
        if (mid > 0) color = "rgba(140,90,10,0.30)";
        else if (mid > -6) color = "rgba(100,60,10,0.18)";
        else if (mid > -12) color = "rgba(60,40,80,0.10)";
        else if (mid > -18) color = "rgba(30,20,55,0.05)";
        if (color) {
          ctx.fillStyle = color;
          ctx.fillRect(x1, pad.top, x2 - x1, pH);
        }
      }
    }

    // Altitude grid lines
    for (const alt of [0, 30, 60, 90]) {
      const y = altToY(alt);
      if (y < pad.top - 1 || y > pad.top + pH + 1) continue;
      ctx.strokeStyle = alt === 0 ? "#374151" : "#1a2535";
      ctx.lineWidth = alt === 0 ? 1 : 0.6;
      ctx.beginPath();
      ctx.moveTo(pad.left, y);
      ctx.lineTo(W - pad.right, y);
      ctx.stroke();
      ctx.fillStyle = "#4b6080";
      ctx.font = "10px Arial";
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      ctx.fillText(`${alt}°`, pad.left - 3, y);
    }

    // Hour grid & time labels
    const hourStep = T <= 12 ? 2 : 4;
    for (let h = 0; h <= T; h += hourStep) {
      const x = hToX(h);
      ctx.strokeStyle = "#1a2535";
      ctx.lineWidth = 0.6;
      ctx.beginPath();
      ctx.moveTo(x, pad.top);
      ctx.lineTo(x, pad.top + pH);
      ctx.stroke();
      if (h > 0 && h < T) {
        const t = new Date(Date.now() + h * 3600000);
        ctx.fillStyle = "#4b6080";
        ctx.font = "10px Arial";
        ctx.textAlign = "center";
        ctx.textBaseline = "top";
        ctx.fillText(
          `${String(t.getHours()).padStart(2, "0")}h`,
          x,
          pad.top + pH + 4
        );
      }
    }
    // "now" label at x=0
    ctx.fillStyle = "#6b8090";
    ctx.font = "10px Arial";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.fillText("now", hToX(0), pad.top + pH + 4);

    if (targetPoints.length < 2) return;

    const hY = altToY(0);

    // Blue fill above horizon
    ctx.beginPath();
    ctx.moveTo(hToX(targetPoints[0].h), hY);
    for (const p of targetPoints) {
      ctx.lineTo(hToX(p.h), Math.min(altToY(p.alt), hY));
    }
    ctx.lineTo(hToX(targetPoints[targetPoints.length - 1].h), hY);
    ctx.closePath();
    ctx.fillStyle = "rgba(96,165,250,0.12)";
    ctx.fill();

    // Altitude curve
    ctx.beginPath();
    for (let i = 0; i < targetPoints.length; i++) {
      const x = hToX(targetPoints[i].h);
      const y = altToY(targetPoints[i].alt);
      i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    }
    ctx.strokeStyle = "#60a5fa";
    ctx.lineWidth = 2;
    ctx.stroke();

    // 30° "good visibility" dashed line
    const y30 = altToY(30);
    if (y30 >= pad.top && y30 <= pad.top + pH) {
      ctx.strokeStyle = "#22c55e40";
      ctx.lineWidth = 0.8;
      ctx.setLineDash([2, 6]);
      ctx.beginPath();
      ctx.moveTo(pad.left, y30);
      ctx.lineTo(W - pad.right, y30);
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // "Now" vertical line
    const nowX = hToX(0);
    ctx.strokeStyle = "#fbbf24";
    ctx.lineWidth = 1.5;
    ctx.setLineDash([4, 3]);
    ctx.beginPath();
    ctx.moveTo(nowX, pad.top);
    ctx.lineTo(nowX, pad.top + pH);
    ctx.stroke();
    ctx.setLineDash([]);
  }

  global.Astronomy = {
    computeAltAz,
    computeCurve,
    computeSunCurve,
    drawSkyDome,
    drawAltChart,
    toJD,
    gmstDeg,
    sunRaDec,
  };
})(window);
