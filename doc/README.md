# AstroBatch Documentation

AstroBatch is an observatory automation and photometric pipeline system for transient astronomy. It combines:

- **Night Manager** — browser-based dashboard to control NINA, the mount, camera, dome (OCS), and run automated multi-target imaging sequences
- **Processing Service** — Python pipeline that takes raw FITS from NAS and produces calibrated, stacked, plate-solved images uploaded to STDWeb for photometry and template subtraction
- **CLI tools** — standalone batch processing and candidate labelling

---

## Documentation Index

| Document | Contents |
|----------|---------|
| [architecture.md](./architecture.md) | System overview, components, data flow |
| [ui.md](./ui.md) | Web dashboard — all tabs and controls |
| [pipeline.md](./pipeline.md) | Processing pipeline — from raw FITS to photometry results |
| [api.md](./api.md) | Full REST API reference (Node + Python service) |
| [configuration.md](./configuration.md) | Environment variables, paths, systemd services |
| [database.md](./database.md) | SQLite schema — all tables and columns |
| [integrations.md](./integrations.md) | NINA, OCS, STDWeb, Astrometry.net, Siril, TNS |
| [sequence.md](./sequence.md) | Automated imaging sequence — step by step |
| [cli.md](./cli.md) | Command-line tools |
