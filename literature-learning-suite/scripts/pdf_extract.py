#!/usr/bin/env python3
"""Extract text or Markdown from a local PDF using optional PyMuPDF."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pdf", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument("--markdown", action="store_true")
    parser.add_argument("--metadata", action="store_true")
    args = parser.parse_args()
    try:
        import fitz
    except ImportError as exc:
        raise SystemExit("Install optional dependency: pip install pymupdf") from exc

    document = fitz.open(args.pdf)
    if args.metadata:
        value = {"pages": document.page_count, **document.metadata}
        content = json.dumps(value, ensure_ascii=False, indent=2)
    elif args.markdown:
        try:
            import pymupdf4llm
        except ImportError as exc:
            raise SystemExit(
                "Install optional dependency: pip install pymupdf4llm"
            ) from exc
        content = pymupdf4llm.to_markdown(document)
    else:
        content = "\n\n".join(page.get_text("text") for page in document)
    if args.output:
        args.output.write_text(content, encoding="utf-8")
        print(args.output)
    else:
        print(content)


if __name__ == "__main__":
    main()
