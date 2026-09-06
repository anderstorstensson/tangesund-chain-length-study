"""Asterionellopsis glacialis validation set, stage 2 of 3: the Label Studio counting project.

Creates a blind project over the PNGs from R/23f_asterionellopsis_select.R: one box per
cell body (rostrum excluded) (manual count = number of boxes) plus three exclusion flags. `--preannotate`
adds the detector's boxes (conf 0.25, IoU 0.30); the project is then not blind.

The API helpers come from ls_preannotate_api.py in the IFCBChainCounts repository's
tools/ folder; set IFCBCHAINCOUNTS_TOOLS (default ../IFCBChainCounts/tools).

Start Label Studio with local file serving and a document root above the images:

    LABEL_STUDIO_LOCAL_FILES_SERVING_ENABLED=true \
    LABEL_STUDIO_LOCAL_FILES_DOCUMENT_ROOT=/path/to/parent/of/this/repo \
    label-studio start

Reads: outputs/tangesund/asterionellopsis_validation/images/*.png; with --preannotate,
    data/chains/chains_asterionellopsis_yolo11s/weights/best.pt
Writes: a Label Studio project and its local-files storage. Annotate, export JSON,
    then run Rscript R/23h_asterionellopsis_analysis.R <export.json>.
Deps: requests; ultralytics with --preannotate
Run (from the repo root): python3 R/23g_asterionellopsis_ls.py --token-file ~/.ls_token
"""
import argparse
import os
import sys
import time
from pathlib import Path

import requests

REPO = Path(__file__).resolve().parents[1]
TOOLS = Path(os.environ.get("IFCBCHAINCOUNTS_TOOLS", REPO.parent / "IFCBChainCounts" / "tools"))
if not (TOOLS / "ls_preannotate_api.py").exists():
    raise SystemExit(f"ls_preannotate_api.py not found in {TOOLS}; clone "
                     "https://github.com/EuropeanIFCBGroup/IFCBChainCounts next to this repository "
                     "or point IFCBCHAINCOUNTS_TOOLS at its tools/ folder")
sys.path.insert(0, str(TOOLS))
from ls_preannotate_api import authenticate, get_token, request, detect_config, image_filename, iter_tasks  # noqa: E402

IMAGES = REPO / "outputs/tangesund/asterionellopsis_validation/images"
WEIGHTS = REPO / "data/chains/chains_asterionellopsis_yolo11s/weights/best.pt"  # the V6.1 production counter

LABEL_CONFIG = """<View>
  <Header value="Count the Asterionellopsis glacialis cells: one box per cell body (the triangular base; exclude the long rostrum). Box partly hidden cells too. Flag the frame instead of boxing if it is not A. glacialis or cannot be counted."/>
  <Image name="image" value="$image" zoom="true" zoomControl="true" brightnessControl="true" contrastControl="true"/>
  <RectangleLabels name="label" toName="image">
    <Label value="cell" background="#39CC7D" hotkey="1"/>
  </RectangleLabels>
  <Choices name="flag" toName="image" choice="multiple" showInline="true">
    <Choice value="not Asterionellopsis glacialis"/>
    <Choice value="colony cut off by frame"/>
    <Choice value="cells not resolvable"/>
  </Choices>
</View>"""


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--url", default="http://localhost:8080")
    p.add_argument("--token")
    p.add_argument("--token-file", dest="token_file")
    p.add_argument("--title", default="A. glacialis validation (blind)")
    p.add_argument("--images", default=str(IMAGES))
    p.add_argument("--preannotate", action="store_true", help="attach production detector boxes (NOT blind)")
    p.add_argument("--sync-timeout", type=int, default=300)
    args = p.parse_args()
    if len(args.title) > 50:
        raise SystemExit(f"Label Studio caps project titles at 50 characters ({len(args.title)} given)")

    images = Path(args.images)
    pngs = sorted(images.glob("*.png"))
    if not pngs:
        raise SystemExit(f"no PNGs in {images}; run Rscript R/23f_asterionellopsis_select.R first")
    print(f"{len(pngs)} images in {images}")

    session = requests.Session()
    token = get_token(args)
    authenticate(session, args.url, token)
    reauth = lambda: authenticate(session, args.url, token)  # noqa: E731

    existing = request(session, "GET", f"{args.url}/api/projects", reauth, params={"page_size": 200}).json()
    existing = existing.get("results", existing) if isinstance(existing, dict) else existing
    for pr in existing:
        if pr.get("title") == args.title:
            raise SystemExit(f"a project titled {args.title!r} already exists (id {pr['id']}); "
                             "delete it in Label Studio or pass another --title")

    r = request(session, "POST", f"{args.url}/api/projects", reauth,
                json={"title": args.title, "label_config": LABEL_CONFIG,
                      "show_skip_button": False, "sampling": "Uniform sampling"})
    r.raise_for_status()
    pid = r.json()["id"]
    print(f"created project {pid}: {args.title!r}")

    r = request(session, "POST", f"{args.url}/api/storages/localfiles", reauth, json={
        "project": pid, "title": "ag validation images", "path": str(images.resolve()),
        "regex_filter": r".*\.png$", "use_blob_urls": True})
    if not r.ok:
        raise SystemExit(f"storage create failed ({r.status_code}): {r.text}\n"
                         "Is LABEL_STUDIO_LOCAL_FILES_SERVING_ENABLED=true and the image folder under "
                         "LABEL_STUDIO_LOCAL_FILES_DOCUMENT_ROOT?")
    sid = r.json()["id"]
    request(session, "POST", f"{args.url}/api/storages/localfiles/{sid}/sync", reauth).raise_for_status()
    deadline = time.time() + args.sync_timeout
    n_tasks = 0
    while time.time() < deadline:
        n_tasks = request(session, "GET", f"{args.url}/api/projects/{pid}", reauth).json().get("task_number") or 0
        if n_tasks >= len(pngs):
            break
        time.sleep(3)
    print(f"synced {n_tasks} tasks")

    if args.preannotate:
        from ultralytics import YOLO
        model = YOLO(str(WEIGHTS))
        project = request(session, "GET", f"{args.url}/api/projects/{pid}", reauth).json()
        from_name, to_name, image_key, label = detect_config(project)
        local = {f.name: f for f in pngs}
        n = 0
        for task in iter_tasks(session, args.url, pid, reauth):
            fname = image_filename(task.get("data", {}), image_key)
            img = local.get(fname) if fname else None
            if img is None:
                continue
            res = model(str(img), conf=0.25, iou=0.30, verbose=False)[0]
            h, w = res.orig_shape
            boxes = [{
                "type": "rectanglelabels", "from_name": from_name, "to_name": to_name,
                "original_width": int(w), "original_height": int(h), "image_rotation": 0,
                "value": {"x": float(x1) * 100, "y": float(y1) * 100, "width": float(x2 - x1) * 100,
                          "height": float(y2 - y1) * 100, "rotation": 0, "rectanglelabels": [label]},
            } for (x1, y1, x2, y2) in res.boxes.xyxyn.cpu().numpy()]
            request(session, "POST", f"{args.url}/api/predictions", reauth,
                    json={"task": task["id"], "result": boxes, "model_version": "asterionellopsis-production",
                          "score": 1.0}).raise_for_status()
            n += 1
        print(f"pre-annotated {n} tasks (this project is NOT blind)")

    print(f"\nopen {args.url}/projects/{pid}/ and annotate; then Export -> JSON and run\n"
          f"  Rscript R/23h_asterionellopsis_analysis.R <export.json>")


if __name__ == "__main__":
    main()
