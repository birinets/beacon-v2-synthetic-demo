"""Command line: beacon-demo synth | serve | query | examples."""

from __future__ import annotations

import argparse
import json
import sys

from . import __version__, data, synth
from .beacon import Beacon


def _beacon(args):
    indexed = data.load(args.data) if args.data else data.from_cohort(synth.make_cohort(args.seed))
    for note in indexed.skipped:
        print(f"skipped: {note}", file=sys.stderr)
    return Beacon(indexed, base_url=args.base_url)


def main(argv=None):
    parser = argparse.ArgumentParser(prog="beacon-demo", description="A Beacon v2 node serving synthetic genomes.")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("synth", help="write the synthetic cohort: datasets.tsv, samples.tsv, terms.tsv and VCFs")
    p.add_argument("outdir")
    p.add_argument("--seed", type=int, default=synth.DEFAULT_SEED)

    for name, text in (("serve", "serve the Beacon over HTTP"),
                       ("query", "answer one request without a server and print the JSON"),
                       ("examples", "answer the example requests and write them, with their responses, as JSON")):
        p = sub.add_parser(name, help=text)
        p.add_argument("--data", help="a data directory; without it, the synthetic cohort is made in memory")
        p.add_argument("--seed", type=int, default=synth.DEFAULT_SEED, help="seed for the in-memory cohort")
        p.add_argument("--base-url", default="http://127.0.0.1:5050/api", help="the URL the node is reached at")
        if name == "serve":
            p.add_argument("--host", default="127.0.0.1")
            p.add_argument("--port", type=int, default=5050)
            p.add_argument("--pretty", action="store_true", help="indent the JSON")
            p.add_argument("--log", action="store_true", help="log each request to standard error")
        elif name == "query":
            p.add_argument("path", help='for example "/g_variants?referenceName=10&start=94781858&referenceBases=G'
                                        '&alternateBases=A"')
            p.add_argument("--post", help="a JSON request body; the request is then a POST")
            p.add_argument("--token", help="registered-demo or controlled-demo")
        else:
            p.add_argument("out", help="output JSON file")

    args = parser.parse_args(argv)
    if args.command == "synth":
        for path in synth.write(synth.make_cohort(args.seed), args.outdir):
            print(path)
        return 0
    beacon = _beacon(args)
    if args.command == "serve":
        from .server import make_server
        server = make_server(beacon, args.host, args.port, pretty=args.pretty, quiet=not args.log)
        print(f"Beacon v2 demo node on http://{args.host}:{args.port}/api (synthetic data). Ctrl-C stops it.",
              file=sys.stderr)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass
        finally:
            server.server_close()
        return 0
    if args.command == "query":
        path, _, query = args.path.partition("?")
        headers = {"Authorization": f"Bearer {args.token}"} if args.token else {}
        method, body = ("POST", json.loads(args.post)) if args.post else ("GET", None)
        status, payload = beacon.handle(method, path, query, body, headers)
        print(json.dumps(payload, indent=2))
        return 0 if status == 200 else 1
    from .examples import run
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(run(beacon), fh, indent=1)
        fh.write("\n")
    print(args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
