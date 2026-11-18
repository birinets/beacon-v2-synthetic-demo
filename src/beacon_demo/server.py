"""Serve the Beacon over HTTP with the Python standard library.

Paths work with or without the ``/api`` prefix: ``/api/g_variants`` and ``/g_variants`` are the same endpoint.
Responses carry ``Access-Control-Allow-Origin: *`` so that a web page on the same machine can query the node.
"""

from __future__ import annotations

import json
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlsplit

MAX_BODY = 1_000_000


def make_handler(beacon, prefix="/api", pretty=False, quiet=True):
    class Handler(BaseHTTPRequestHandler):
        server_version = "beacon-v2-synthetic-demo"

        def _send(self, status, payload):
            body = json.dumps(payload, indent=2 if pretty else None).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Access-Control-Allow-Origin", "*")
            self.send_header("Access-Control-Allow-Headers", "Authorization, Content-Type")
            self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _relative(self, path):
            if path == prefix or path.startswith(prefix + "/"):
                return path[len(prefix):] or "/"
            return path

        def do_OPTIONS(self):  # CORS preflight
            self.send_response(204)
            self.send_header("Access-Control-Allow-Origin", "*")
            self.send_header("Access-Control-Allow-Headers", "Authorization, Content-Type")
            self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
            self.end_headers()

        def do_GET(self):
            url = urlsplit(self.path)
            self._send(*beacon.handle("GET", self._relative(url.path), url.query, None, dict(self.headers)))

        def do_POST(self):
            url = urlsplit(self.path)
            size = int(self.headers.get("Content-Length") or 0)
            if size > MAX_BODY:
                self._send(*beacon.error(413, "request body too large"))
                return
            raw = self.rfile.read(size) if size else b""
            try:
                body = json.loads(raw.decode("utf-8")) if raw.strip() else None
            except (UnicodeDecodeError, json.JSONDecodeError) as exc:
                self._send(*beacon.error(400, f"the request body is not valid JSON: {exc}"))
                return
            self._send(*beacon.handle("POST", self._relative(url.path), url.query, body, dict(self.headers)))

        def log_message(self, fmt, *args):
            if not quiet:
                sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))

    return Handler


def make_server(beacon, host="127.0.0.1", port=5050, prefix="/api", pretty=False, quiet=True):
    return ThreadingHTTPServer((host, port), make_handler(beacon, prefix, pretty, quiet))
