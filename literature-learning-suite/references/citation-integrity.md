# Citation Integrity

## Verification Ladder

1. Resolve stable identifier in the primary database.
2. Cross-check title, authors, year, and venue.
3. Check document status: active, corrected, withdrawn, or retracted.
4. Retrieve citation metadata from DOI or source records.
5. Verify that the cited source supports the exact claim.

## Claim-Level Verification

For every important citation, record:

- Claim being supported.
- Exact section, figure, table, or page.
- Whether support is direct or inferred.
- Population/model and conditions.
- Caveats that change interpretation.

An abstract can support broad findings, but usually not detailed mechanistic, numerical, or
methodological claims.

## BibTeX

Prefer DOI content negotiation:

```http
GET https://doi.org/{doi}
Accept: application/x-bibtex
```

If BibTeX must be generated from metadata, mark it as generated and preserve the source JSON.

## Exclusion

Exclude a citation from the verified bibliography when:

- The identifier does not resolve.
- Title and DOI refer to different records.
- The cited claim is absent or contradicted.
- The paper is withdrawn and the withdrawal matters to the claim.
- Only an unverifiable secondary mention is available.
