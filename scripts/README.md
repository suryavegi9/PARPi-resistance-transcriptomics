# Analysis scripts

The staged copies of 53 numbered R scripts and one alternative D1 script are imported. The staging manifest records original and staged SHA-256 digests; portability changes are limited to the working-directory setup. See ../provenance/DEPENDENCY_AUDIT.md for workflow mapping, missing historical sources and input requirements.

Scripts 62 and 63 have their frozen input tables available. They have not been executed during this audit. Both produce output files; testing must occur in an isolated copy. Do not execute all scripts indiscriminately in numerical order or overwrite the frozen results.

Upstream scripts require source expression files and intermediate results. The missing original D4 primary regional model and historical D1–D6 master-generation code must not be represented as recovered. Historical environment records are incomplete.
