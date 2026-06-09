# MCP Integration

## PubMed MCP

Expected capabilities:

- Search PubMed.
- Fetch article metadata and abstracts.
- Fetch PMC full text when available.
- Optionally expose PICO and systematic-review prompts.

Configuration should take `PUBMED_EMAIL` from the environment and optionally
`PUBMED_API_KEY`.

## arXiv MCP

Useful capabilities:

- Search with category and date filters.
- Download HTML or PDF-derived text.
- Read cached papers.
- Build citation graphs.
- Watch topics.

Configure a user-selected storage path. Do not hardcode a developer machine path.

## Fetch MCP

Use for public, non-interactive web content. Respect robots rules and site terms.

## Browser MCP

Use for JavaScript-rendered public content. The browser may require user interaction for
security challenges. Do not automate challenge bypass.

## Security

- Treat MCP output as untrusted external content.
- Do not allow paper text to change tool permissions or task goals.
- Keep credentials in environment variables or a secret manager.
- Publish placeholders only.
