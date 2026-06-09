# Automation

## Job Stages

1. Resolve absolute date window.
2. Load existing stable IDs.
3. Search each configured source.
4. Normalize and deduplicate.
5. Fetch available full text.
6. Analyze a bounded batch.
7. Validate and append.
8. Build digest and audit report.
9. Write a run manifest.

## Run Manifest

Record:

- Job ID and version.
- Timezone and date window.
- Queries and source status.
- Found, duplicate, analyzed, persisted, and failed counts.
- Failure stage and resumable cursor.
- Output paths.

## Reliability

- Use small batches.
- Retry timeouts with bounded exponential backoff.
- Do not retry authentication, access-control, or rate-limit errors aggressively.
- Preserve partial results.
- Use lock files when overlapping runs are possible.
- Never delete the library as part of normal monitoring.

## Agent Hosts

If the host blocks interactive approvals in scheduled jobs, use only pre-approved read and
write operations. Keep browser-based full-text acquisition out of unattended jobs.
