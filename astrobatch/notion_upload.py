"""Upload daily light-curve PNGs to a Notion page."""

from __future__ import annotations

import logging
import os
import time
from pathlib import Path

import requests

log = logging.getLogger(__name__)

NOTION_API = "https://api.notion.com/v1"
NOTION_VERSION = os.environ.get("NOTION_VERSION", "2022-06-28")
DEFAULT_PARENT_PAGE_ID = "36c8e68e17a880f6a89be6ab3731d873"


def notion_page_id(page_id: str) -> str:
    """Normalize 32-char hex to UUID with dashes."""
    s = page_id.replace("-", "").strip()
    if len(s) != 32:
        return page_id
    return f"{s[:8]}-{s[8:12]}-{s[12:16]}-{s[16:20]}-{s[20:]}"


def _token() -> str | None:
    return os.environ.get("NOTION_TOKEN") or os.environ.get("NOTION_API_KEY")


def _parent_page_id() -> str:
    raw = os.environ.get(
        "NOTION_LIGHTCURVES_PAGE_ID",
        DEFAULT_PARENT_PAGE_ID,
    )
    return notion_page_id(raw)


def _headers(*, json_body: bool = True) -> dict[str, str]:
    token = _token()
    if not token:
        raise RuntimeError("NOTION_TOKEN is not set")
    h = {
        "Authorization": f"Bearer {token}",
        "Notion-Version": NOTION_VERSION,
    }
    if json_body:
        h["Content-Type"] = "application/json"
    return h


def _request(method: str, path: str, **kwargs) -> dict:
    url = f"{NOTION_API}{path}"
    resp = requests.request(method, url, headers=_headers(), timeout=120, **kwargs)
    if not resp.ok:
        log.error("Notion %s %s → %s: %s", method, path, resp.status_code, resp.text[:500])
        if resp.status_code == 404 and "shared with your integration" in resp.text:
            raise RuntimeError(
                "Notion page not accessible. Open the Lightcurves page in Notion → "
                "⋯ → Connections → add your Astrobatch integration."
            ) from None
        resp.raise_for_status()
    return resp.json()


def upload_png(file_path: Path) -> str:
    """Upload a PNG; return file_upload id for use in image blocks."""
    created = _request(
        "POST",
        "/file_uploads",
        json={
            "filename": file_path.name,
            "content_type": "image/png",
        },
    )
    upload_id = created["id"]
    with file_path.open("rb") as fh:
        send = requests.post(
            f"{NOTION_API}/file_uploads/{upload_id}/send",
            headers={
                "Authorization": _headers(json_body=False)["Authorization"],
                "Notion-Version": NOTION_VERSION,
            },
            files={"file": (file_path.name, fh, "image/png")},
            timeout=180,
        )
    if not send.ok:
        log.error("Notion file send → %s: %s", send.status_code, send.text[:500])
        send.raise_for_status()

    for _ in range(30):
        status = _request("GET", f"/file_uploads/{upload_id}")
        state = status.get("status")
        if state == "uploaded":
            return upload_id
        if state in ("failed", "expired"):
            raise RuntimeError(f"Notion file upload {upload_id} ended with status={state}")
        time.sleep(0.5)
    raise RuntimeError(f"Notion file upload {upload_id} did not reach uploaded state")


def _find_child_page(parent_id: str, title: str) -> str | None:
    cursor = None
    while True:
        params: dict = {"page_size": 100}
        if cursor:
            params["start_cursor"] = cursor
        data = _request("GET", f"/blocks/{parent_id}/children", params=params)
        for block in data.get("results", []):
            if block.get("type") != "child_page":
                continue
            cp = block.get("child_page", {})
            title_field = cp.get("title", "")
            if isinstance(title_field, str):
                plain = title_field
            else:
                plain = "".join(t.get("plain_text", "") for t in title_field)
            if plain.strip() == title:
                return block["id"]
        if not data.get("has_more"):
            break
        cursor = data.get("next_cursor")
    return None


def _create_day_page(parent_id: str, title: str) -> str:
    existing = _find_child_page(parent_id, title)
    if existing:
        return existing
    page = _request(
        "POST",
        "/pages",
        json={
            "parent": {"page_id": parent_id},
            "properties": {
                "title": [{"type": "text", "text": {"content": title}}],
            },
        },
    )
    return page["id"]


def _append_blocks(page_id: str, children: list[dict]) -> None:
    chunk = 100
    for i in range(0, len(children), chunk):
        _request(
            "PATCH",
            f"/blocks/{page_id}/children",
            json={"children": children[i : i + chunk]},
        )


def _heading(level: int, text: str) -> dict:
    key = f"heading_{level}"
    return {
        "object": "block",
        "type": key,
        key: {
            "rich_text": [{"type": "text", "text": {"content": text}}],
        },
    }


def _image_block(file_upload_id: str) -> dict:
    return {
        "object": "block",
        "type": "image",
        "image": {
            "type": "file_upload",
            "file_upload": {"id": file_upload_id},
        },
    }


def upload_daily_summary(summary: dict) -> dict | None:
    """Create/update a dated Notion sub-page and attach today's light curves."""
    token = _token()
    if not token:
        log.info("NOTION_TOKEN not set — skipping Notion upload")
        return None

    targets = summary.get("targets") or []
    if not targets:
        log.info("No light curves to upload to Notion")
        return None

    parent_id = _parent_page_id()
    day_dir = Path(summary["output_dir"])
    day_title = day_dir.name

    try:
        day_page_id = _create_day_page(parent_id, day_title)
        blocks: list[dict] = [
            _heading(2, f"Run {summary.get('generated_at', '')[:16]}"),
        ]

        for entry in targets:
            display = entry.get("display") or entry.get("target", "")
            blocks.append(_heading(2, display))

            ap = entry.get("aperture")
            if ap and Path(ap).is_file():
                blocks.append(_heading(3, "Aperture photometry"))
                blocks.append(_image_block(upload_png(Path(ap))))

            sub = entry.get("subtraction")
            if sub and Path(sub).is_file():
                blocks.append(_heading(3, "Template subtraction"))
                blocks.append(_image_block(upload_png(Path(sub))))

        _append_blocks(day_page_id, blocks)
        url = f"https://www.notion.so/{day_page_id.replace('-', '')}"
        result = {
            "notion_day_page_id": day_page_id,
            "notion_url": url,
            "parent_page_id": parent_id,
        }
        log.info("Uploaded %d target(s) to Notion → %s", len(targets), url)
        return result
    except Exception as exc:
        log.exception("Notion upload failed: %s", exc)
        raise


def is_configured() -> bool:
    return bool(_token())
